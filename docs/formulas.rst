Formulas
========

This page gives the exact mathematical formula behind every parameter the
package can produce, together with the source code that implements it and
the reasoning behind each design choice. Formulas are transcribed directly
from the implementation (not textbook restatements) -- the code blocks below
are pulled live from the source via Sphinx's ``literalinclude``, so they
always match the installed version of the package. Notation follows the
source code: :math:`\varepsilon` denotes a small constant (typically
:math:`10^{-12}`, sometimes :math:`10^{-30}`) added to denominators to avoid
division by zero.

Position / contact-force pipeline (:mod:`footprint_attributes.position`)
--------------------------------------------------------------------------

Pipeline overview
~~~~~~~~~~~~~~~~~~

For every touching pair of footprints, the shared-boundary segments are
extracted as individual 2-point edges. For each edge segment :math:`e`
(belonging to building :math:`i`, touching neighbour :math:`j`), the edge
midpoint and force vector are

.. math::

   \mathbf{f}_e = h_i \, \ell_e \, \hat{\mathbf{n}}_e

where :math:`h_i` is building :math:`i`'s height (1.0 if no height column is
given), :math:`\ell_e` the edge length, and :math:`\hat{\mathbf{n}}_e` the
unit outward normal of the segment (perpendicular to the shared wall). Each
building's raw force vectors are summed to a resultant.

``contact_force``
~~~~~~~~~~~~~~~~~~

.. math::

   \mathbf{F}_i = \sum_{e \in E_i} \mathbf{f}_e
   = \sum_{e \in E_i} h_i \, \ell_e \, \hat{\mathbf{n}}_e

.. math::

   \text{contact\_force}_i = \frac{\lVert \mathbf{F}_i \rVert}{\sqrt{A_i + \varepsilon}}

where :math:`E_i` is the set of touching boundary edge segments of building
:math:`i` and :math:`A_i` its footprint area. The per-edge force magnitude is
proportional to *wall area* (height x length, a "virtual unit pressure"
analogy); the resultant is normalised by :math:`\sqrt{A_i}` (confinement/
pressure per unit boundary, dimensionally like force/length). Buildings with
no touching neighbours get ``contact_force = 0``. Classification divides
``force`` by ``height`` before comparing to thresholds, since raw force
scales linearly with height but thresholds are calibrated for ``height = 1``.

.. literalinclude:: ../src/footprint_attributes/geometry.py
   :pyobject: edge_normal
   :language: python

.. literalinclude:: ../src/footprint_attributes/position.py
   :pyobject: contact_forces_df
   :language: python

The full ``contact_forces_df`` function above computes ``force``,
``confinementRatio``, ``angularAcc`` and ``angle`` together in one pass (they
share the same per-edge force vectors), so it is quoted once here rather than
split by column; the sections below explain each of its outputs in turn.

``contact_confinementRatio``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \text{contact\_confinementRatio}_i =
   \frac{\left(\sum_{e \in E_i} \lVert \mathbf{f}_e \rVert\right) - \lVert \mathbf{F}_i \rVert}
        {\left(\sum_{e \in E_i} \lVert \mathbf{f}_e \rVert\right) + \varepsilon}

Ranges :math:`[0, 1]`: 0 when all forces point the same way (nothing cancels,
fully lateral); tends to 1 when forces cancel almost completely in the vector
sum (fully enclosed / surrounded on opposing sides). This is the metric that
captures "opposing-wall cancellation" -- see the note on ``contact_angle``
below for why that matters. It is height-invariant (ratio of two force sums).

``contact_angularAcc``
~~~~~~~~~~~~~~~~~~~~~~~

For each edge, the 2-D torque of that edge's force about the building
centroid :math:`\mathbf{c}_i` is

.. math::

   M_e = r_{e,x} f_{e,y} - r_{e,y} f_{e,x}, \qquad \mathbf{r}_e = \text{midpoint}_e - \mathbf{c}_i

If the edge midpoint lies closer to the centroid than
:math:`\text{min\_dist} = \text{minRadius} \cdot \sqrt{A_i/\pi}` (a fraction
of the equivalent-circle radius), only the sign-consistent ("conservative")
component of that edge's torque is counted, i.e. a close-in edge's moment is
included only if it reduces net torque magnitude rather than amplifying it.
The building-level moment magnitude is the minimum in absolute value across
four such filtered accumulations:

.. math::

   \text{momentum\_mag}_i = \min_k \left| \sum_{e \in E_i} M_e^{(k)} \right|, \quad k \in \{1,2,3,4\}

.. math::

   \text{contact\_angularAcc}_i = \frac{\text{momentum\_mag}_i}{I_{z,i} + \varepsilon} \cdot A_i

where :math:`I_{z,i}` is the polygon's polar second moment of area about its
own centroid. This behaves like a moment/area-normalised angular
acceleration. It scales linearly with height (through the force magnitudes);
classification divides by ``height`` before comparing to its threshold.

.. literalinclude:: ../src/footprint_attributes/geometry.py
   :pyobject: edge_momentum
   :language: python

``contact_angle``
~~~~~~~~~~~~~~~~~~

``contact_angle`` is the force-weighted mean deviation of each touching
edge's own force from the building's resultant force.

**Step 1 -- per-edge unfolded angle to the resultant:**

.. math::

   \cos\theta_e = \operatorname{clip}\!\left(
     \frac{\mathbf{f}_e \cdot \mathbf{F}_i}{\lVert \mathbf{f}_e \rVert \, \lVert \mathbf{F}_i \rVert + \varepsilon},
     -1, 1
   \right), \qquad
   \theta_e^{\text{raw}} = \arccos(\cos\theta_e) \in [0, \pi]

**Step 2 -- fold to** :math:`[0, \pi/2]`:

.. math::

   \theta_e = \min\!\left(\theta_e^{\text{raw}}, \; \pi - \theta_e^{\text{raw}}\right) \in \left[0, \tfrac{\pi}{2}\right]

If either :math:`\mathbf{f}_e` or :math:`\mathbf{F}_i` has near-zero norm
(:math:`< 10^{-12}`), :math:`\theta_e = 0` by convention.

**Step 3 -- force-weighted aggregation:**

.. math::

   \text{contact\_angle}_i =
   \frac{2\sum_{e \in E_i} \theta_e \, \lVert \mathbf{f}_e \rVert}{\sum_{e \in E_i} \lVert \mathbf{f}_e \rVert + \varepsilon}

i.e. twice the plain force-weighted mean of :math:`\theta_e` (the factor of
2 is applied at aggregation time and is a fixed scaling of the module, not
part of the fold).

.. note::

   **Why the angle is folded to** :math:`[0, \pi/2]` **(v1.0.0 correction).**
   Previously, ``contact_angle`` used the unfolded :math:`\theta_e^{\text{raw}}
   \in [0, \pi]` directly. This produced invalid averages whenever a building
   was touched on two *opposite* sides by neighbours of very unequal wall
   length -- e.g. a full 10 m wall on the North side and only a 3 m wall on
   the South side. The resultant :math:`\mathbf{F}_i` points mostly toward
   the larger (North) force, so the *minority* South edge's own force vector
   is nearly anti-parallel to :math:`\mathbf{F}_i`
   (:math:`\theta_e^{\text{raw}} \approx \pi`), even though geometrically the
   South neighbour sits directly *opposite* the North one, not
   perpendicular/corner-like at all. Averaging an angle near 0 deg with one
   near 180 deg does not sensibly represent "how spread out are these
   forces" -- the mean is not meaningful for angles measured in that way.

   Without the fold, that near-:math:`\pi` angle inflated the force-weighted
   mean by exactly as much as a genuinely perpendicular (corner-like)
   neighbour would, making "opposite, unequal walls" and "true corner"
   indistinguishable via ``contact_angle`` alone -- and buildings like this
   were misclassified as ``"corner"`` instead of ``"lateral"``/``"confined"``.

   Opposition/cancellation between unequal opposing forces is already
   captured by ``contact_confinementRatio`` above. Folding
   :math:`\theta_e^{\text{raw}} \mapsto \min(\theta_e^{\text{raw}}, \pi -
   \theta_e^{\text{raw}})` treats a near-180 deg deviation the same as a
   near-0 deg one, removing the double-counting of the same cancellation
   effect across two different statistics. Only genuinely lateral spread
   (forces pointing in different, non-opposite directions -- the true
   "corner" signature) now drives ``contact_angle`` up. This is the same
   fold used by ``angle_between_0_90`` in :mod:`footprint_attributes.geometry`
   for undirected-axis comparisons (e.g. ``ASCE7_parallelityAngle``,
   ``bearing``).

.. literalinclude:: ../src/footprint_attributes/position.py
   :pyobject: resultant_angle
   :language: python

``contact_height``
~~~~~~~~~~~~~~~~~~~~

.. math::

   \text{contact\_height}_i = h_i

The height value used for the computation (1.0 if no height column is
given); carried through so downstream reuse of the prefixed columns can
undo the height scaling of ``force``/``angularAcc``.

``relativePosition``
~~~~~~~~~~~~~~~~~~~~~~

Classification, evaluated on height-normalised force/angularAcc when a
height column is available:

.. math::

   \tilde{F}_i = \frac{\text{force}_i}{h_i}, \qquad
   \tilde{a}_i = \frac{\text{angularAcc}_i}{h_i}

Priority order (later rules override earlier ones):

1. Default: ``"isolated"``.
2. ``"lateral"`` if :math:`\tilde{F}_i > \text{minForce}`.
3. ``"corner"`` if currently ``"lateral"`` **and** ``contact_angle``
   :math:`> \text{minAngle}`.
4. ``"confined"`` if ``confinementRatio`` :math:`> \text{minConfinement}`
   (overrides lateral/corner).
5. ``"torque"`` if currently ``"corner"`` or ``"confined"`` **and**
   :math:`\tilde{a}_i > \text{minAngularAcc}`.

``"isolated"`` buildings (no touching edges) never satisfy rule 2 since
:math:`\tilde F_i = 0`.

.. literalinclude:: ../src/footprint_attributes/position.py
   :pyobject: _Position._classify
   :language: python

Shape irregularity indices (:mod:`footprint_attributes.shape`)
------------------------------------------------------------------

Code-independent indices
~~~~~~~~~~~~~~~~~~~~~~~~~~

**polsby_popper** (compactness):

.. math::

   \text{polsby\_popper} = \frac{4\pi A}{P^2 + \varepsilon}

:math:`A` is the area and :math:`P` the perimeter of the hole-filled
footprint. Range :math:`(0, 1]`; 1 for a perfect circle.

.. literalinclude:: ../src/footprint_attributes/shape.py
   :pyobject: polsby_popper
   :language: python

**convex_hull_irregularity**:

.. math::

   \text{convex\_hull\_irregularity} = \frac{A_{\text{hull}} - A}{A + \varepsilon}

:math:`A` is the hole-filled footprint's area, :math:`A_{\text{hull}}` its
convex hull's area. 0 for a convex shape; larger for deeper/larger setbacks
relative to the footprint's own area.

.. literalinclude:: ../src/footprint_attributes/shape.py
   :pyobject: convex_hull_irregularity
   :language: python

**inertia_circle_ratio**:

.. math::

   \text{inertia\_circle\_ratio} = \frac{I_{z,\text{circle}}}{|I_z| + \varepsilon}, \qquad
   I_{z,\text{circle}} = \frac{A^2}{2\pi}

:math:`I_z` is the polar second moment of area of the real (not hole-filled)
footprint about its own centroid; :math:`I_{z,\text{circle}}` is the polar
moment of a circle of the same area. Range :math:`(0, 1]`, 1 for a circle
(the circle maximises :math:`I_z` for a given area).

.. literalinclude:: ../src/footprint_attributes/shape.py
   :pyobject: inertia_circle_ratio
   :language: python

EC8 (Eurocode 8)
~~~~~~~~~~~~~~~~~~

All EC8 quantities use the centre-of-mass/centre-of-stiffness model and
principal-inertia computation below, plus the Mohr's-circle worst-case
optimisation described in *Eccentricity optimisation*.

**EC8_eccentricityRatio** (:math:`e/r_t`), limit :math:`\le 0.30`:

.. math::

   \text{EC8\_eccentricityRatio} = \frac{e \, |\cos(x_{\text{opt}} - b)|}{r_t(x_{\text{opt}}) + \varepsilon}

.. literalinclude:: ../src/footprint_attributes/shape.py
   :pyobject: EC8EccentricityRatio
   :language: python

**EC8_radiusRatio** (:math:`r_t/r_g`), limit :math:`\ge 1.0`:

.. math::

   \text{EC8\_radiusRatio} = \frac{r_t(x_{\text{opt}})}{r_g + \varepsilon}, \qquad
   r_g = \sqrt{\frac{I_1+I_2}{A + \varepsilon}}

.. literalinclude:: ../src/footprint_attributes/shape.py
   :pyobject: EC8RadiusRatio
   :language: python

**EC8_compactness**, limit :math:`\ge 0.95`:

.. math::

   \text{EC8\_compactness} =
   \begin{cases}
     1 - \dfrac{A_{\text{gap}}}{A_{\text{filled}} + \varepsilon} & A_{\text{filled}} > \varepsilon \\
     1 & \text{otherwise}
   \end{cases}

:math:`A_{\text{gap}}` is the area of the single largest connected component
of ``hull(filled).difference(filled)`` -- not the sum of all setback pieces.

.. literalinclude:: ../src/footprint_attributes/shape.py
   :pyobject: EC8Compactness
   :language: python

ASCE 7
~~~~~~~~

**ASCE7_setbackRatio**, limit :math:`\le 0.20`:

.. math::

   \text{ASCE7\_setbackRatio} = \min\!\left(\frac{b_1}{L_1 + \varepsilon}, \frac{b_2}{L_2 + \varepsilon}\right)

using the dominant setback piece and bounding-box directions (see
*GNDT setback construction* below).

.. literalinclude:: ../src/footprint_attributes/shape.py
   :pyobject: ASCE7SetbackRatio
   :language: python

**ASCE7_holeRatio**, limit :math:`\le 0.25`:

.. math::

   \text{ASCE7\_holeRatio} = \max_{\text{holes}} \frac{A_{\text{filled}} - A_{\text{poly}}}{A_{\text{filled}} + \varepsilon}

The largest single hole's area fraction (holes below 0.1% of the filled
area are ignored).

.. literalinclude:: ../src/footprint_attributes/shape.py
   :pyobject: ASCE7HoleRatio
   :language: python

.. literalinclude:: ../src/footprint_attributes/geometry.py
   :pyobject: max_hole_area_ratio
   :language: python

**ASCE7_parallelityAngle**, limit :math:`\le 5°`:

.. math::

   \text{ASCE7\_parallelityAngle} = \frac{180}{\pi}\,\arccos\!\bigl(|\hat{\mathbf{x}} \cdot \text{dir}_1|\bigr)

folded to :math:`[0°, 90°]`; :math:`\text{dir}_1` is the building's longer
bounding-box axis and :math:`\hat{\mathbf{x}} = (1, 0)`.

.. literalinclude:: ../src/footprint_attributes/shape.py
   :pyobject: ASCE7ParalelityAngle
   :language: python

GNDTII (Italian GNDT Level II)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All use the dominant :math:`(a, L)` configuration from the inscribed-circle
construction, picked per building as whichever of :math:`(L_1, a_1)` /
:math:`(L_2, a_2)` maximises :math:`L \cdot a`:

.. math::

   (a_{\text{dom}}, L_{\text{dom}}) =
   \begin{cases}
     (a_1, L_1) & L_1 a_1 \ge L_2 a_2 \\
     (a_2, L_2) & \text{otherwise}
   \end{cases}

.. literalinclude:: ../src/footprint_attributes/shape.py
   :pyobject: _gndt_dominant_a
   :language: python

**GNDTII_beta1_mainShapeSlenderness**:

.. math::

   \beta_1 = \frac{a_{\text{dom}}}{L_{\text{dom}} + \varepsilon}

.. literalinclude:: ../src/footprint_attributes/shape.py
   :pyobject: GNDTIIBeta1MainShapeSlenderness
   :language: python

**GNDTII_beta2_setbackRatio**: same formula as ``ASCE7_setbackRatio`` above,
under a different compliance table.

.. literalinclude:: ../src/footprint_attributes/shape.py
   :pyobject: GNDTIIBeta2SetbackRatio
   :language: python

**GNDTII_beta4_eccentricityRatio**:

.. math::

   \beta_4 = \frac{\lVert \mathbf{cm} - \mathbf{cs} \rVert}{a_{\text{dom}} + \varepsilon}

.. literalinclude:: ../src/footprint_attributes/shape.py
   :pyobject: GNDTIIBeta4EccentricityRatio
   :language: python

**GNDTII_beta6_setbackSlenderness**:

.. math::

   \beta_6 = \frac{c}{b + \varepsilon}

where :math:`b` is the winning setback configuration's own width and
:math:`c` the protrusion width of the solid footprint measured perpendicular
to :math:`b`, through the setback piece's centroid.

.. literalinclude:: ../src/footprint_attributes/shape.py
   :pyobject: GNDTIIBeta6SetbackSlenderness
   :language: python

CSCR 2010 (Costa Rica)
~~~~~~~~~~~~~~~~~~~~~~~~

**CSCR2010_eccentricityRatio** (:math:`e/l`):

.. math::

   \text{CSCR2010\_eccentricityRatio} =
   \frac{|e \cos(x_{\text{opt}} - b)|}{l(x_{\text{opt}}) + \varepsilon}

Full derivation of :math:`e`, :math:`b`, :math:`x_{\text{opt}}`, :math:`l`
in *Eccentricity optimisation* below.

.. literalinclude:: ../src/footprint_attributes/shape.py
   :pyobject: CSCR2010EccentricityRatio
   :language: python

NTC-23 (Mexico)
~~~~~~~~~~~~~~~~~

**NTC23_setbackRatio**: identical formula to ``ASCE7_setbackRatio`` /
``GNDTII_beta2_setbackRatio`` above, against a more lenient limit
(:math:`\le 0.40`).

.. literalinclude:: ../src/footprint_attributes/shape.py
   :pyobject: NTC23SetbackRatio
   :language: python

**NTC23_holeRatio**, limit :math:`\le 0.40`:

.. math::

   \text{NTC23\_holeRatio} = \max_{\text{holes}} \frac{h}{L + \varepsilon}

:math:`h` is the shorter side of the hole's own minimum rotated bounding box
(independent of the building's own axes); :math:`L` is the length of the
intersection segment of the building's filled footprint with a line through
the hole's centroid, drawn along the hole's own longer axis direction.

.. literalinclude:: ../src/footprint_attributes/shape.py
   :pyobject: NTC23HoleRatio
   :language: python

.. literalinclude:: ../src/footprint_attributes/geometry.py
   :pyobject: hole_h_over_l
   :language: python

Slenderness
~~~~~~~~~~~~~

Plan slenderness (``bbox`` or ``inertia`` method):

.. math::

   \text{slenderness\_bbox} = \frac{L_1^{\text{bbox}}}{L_2^{\text{bbox}} + \varepsilon}, \qquad
   \text{slenderness\_inertia} = \frac{L_1^{\text{inertia}}}{L_2^{\text{inertia}} + \varepsilon}

where :math:`L_1, L_2` come from the bounding-box / inertia direction
methods below. Checked against the EC8 slenderness limit table for both
variants.

Vertical slenderness (``vertical=True``, either method):

.. math::

   \text{slenderness}_{\text{vertical}} = \frac{h}{L_2 + \varepsilon}

:math:`h` is from the given height column, :math:`L_2` the shorter plan
dimension for the selected method.

.. literalinclude:: ../src/footprint_attributes/shape.py
   :pyobject: _SlendernessMethod.compute
   :language: python

Direction / orientation (:mod:`footprint_attributes.direction`)
---------------------------------------------------------------------

Two methods, both returning :math:`(L_1, \text{dir}_1, L_2, \text{dir}_2, \text{bearing})`.

**direction.inertia** (used for the ``bearing`` column and EC8-style
eccentricity): the inertia tensor's larger eigenvalue corresponds to the
axis of greater bending stiffness, which is physically the *shorter* plan
dimension (a moment-of-inertia tensor is like a covariance matrix with x/y
roles swapped), so the eigenvector roles are swapped:

.. math::

   \text{dir}_1 = \text{eig\_dir}_2 \quad (\text{physically the longer side}), \qquad
   \text{dir}_2 = \text{eig\_dir}_1 \quad (\text{the weak/shorter axis})

Side lengths via a closed-form approximation:

.. math::

   s = \sqrt{\frac{I_1}{I_2 + \varepsilon}} \quad (\text{slenderness}), \qquad
   L_1 = \sqrt{A \cdot s}, \qquad L_2 = \frac{A}{L_1 + \varepsilon}

exact for a true rectangle (reduces to the real side lengths), an average
side length otherwise.

.. literalinclude:: ../src/footprint_attributes/direction.py
   :pyobject: inertia
   :language: python

.. literalinclude:: ../src/footprint_attributes/geometry.py
   :pyobject: inertia_side_lengths
   :language: python

**direction.bbox**: :math:`(L_1, \text{dir}_1, L_2, \text{dir}_2)` directly
from the exact minimum rotated bounding rectangle (rotating calipers), with
:math:`L_1 \ge L_2` enforced by swapping if needed.

.. literalinclude:: ../src/footprint_attributes/direction.py
   :pyobject: bbox
   :language: python

.. literalinclude:: ../src/footprint_attributes/geometry.py
   :pyobject: min_bounding_box
   :language: python

**bearing** column:

.. math::

   \phi = \arctan_2(\text{dir}_{2,x}, \, \text{dir}_{2,y}), \qquad
   \phi' = ((\phi + \pi) \bmod 2\pi) - \pi

.. math::

   \text{bearing} =
   \frac{180}{\pi}
   \begin{cases}
     \phi' - \pi & \phi' > \pi/2 \\
     \phi' + \pi & \phi' < -\pi/2 \\
     \phi' & \text{otherwise}
   \end{cases}
   \; \in [-90°, 90°]

Applied to :math:`\text{dir}_2` -- the weak (shorter-dimension) axis.
Clockwise from geographic North (+y in UTM); folded to :math:`\pm 90°` since
an axis is undirected (same fold family as ``contact_angle`` and
``ASCE7_parallelityAngle`` above).

.. literalinclude:: ../src/footprint_attributes/geometry.py
   :pyobject: bearing_from_dir
   :language: python

Eccentricity optimisation (:mod:`footprint_attributes.eccentricity`)
--------------------------------------------------------------------------

Shared setup for both EC8 and CSCR2010, based on a Mohr's-circle
parameterisation of the principal moments of inertia :math:`I_1 \ge I_2`:

.. math::

   c = \frac{I_1 + I_2}{2} \quad (\text{Mohr's-circle centre}), \qquad
   r = \frac{I_1 - I_2}{2} \quad (\text{Mohr's-circle radius})

.. math::

   I(x) = c - r\cos(2x) \quad (\text{moment of inertia at analysis angle } x)

.. math::

   b = \operatorname{atan2}\bigl(\text{dir}_1 \times \mathbf{e}, \; \text{dir}_1 \cdot \mathbf{e}\bigr)
   \quad (\text{signed angle from dir}_1 \text{ to eccentricity vector } \mathbf{e} = \mathbf{cm}-\mathbf{cs})

:math:`b = 0` if :math:`\lVert \mathbf{e} \rVert < 10^{-10}` (negligible
eccentricity). :math:`x_{\text{opt}}` is found by a coarse grid search over
:math:`x \in [0, \pi)` followed by golden-section refinement, maximising the
relevant objective below.

.. literalinclude:: ../src/footprint_attributes/eccentricity.py
   :pyobject: mohr_params
   :language: python

.. literalinclude:: ../src/footprint_attributes/eccentricity.py
   :pyobject: _signed_angle
   :language: python

**EC8 objective and outputs**:

.. math::

   \text{objective}_{\text{EC8}}(x) = \cos^2(x-b)\,\bigl(c - r\cos(2x)\bigr), \qquad
   x_{\text{opt}} = \arg\max_x \; \text{objective}_{\text{EC8}}(x)

.. math::

   I_t = I_1 + I_2 + A\,e^2, \qquad e = \lVert \mathbf{e} \rVert, \qquad
   I_j(x_{\text{opt}}) = c - r\cos(2 x_{\text{opt}})

.. math::

   r_t = \sqrt{\frac{I_t}{I_j(x_{\text{opt}}) + \varepsilon}}, \qquad
   r_g = \sqrt{\frac{I_1+I_2}{A + \varepsilon}}

:math:`x_{\text{opt}} = 0` trivially when :math:`e < 10^{-10}`.

.. literalinclude:: ../src/footprint_attributes/eccentricity.py
   :pyobject: optimise_ec8
   :language: python

**CSCR 2010 objective and outputs**:

.. math::

   I_{\text{max}}(x) = c + r\cos(2x), \qquad I_{\text{min}}(x) = c - r\cos(2x)

.. math::

   \text{objective}_{\text{CSCR}}(x) =
   \begin{cases}
     0 & |I_{\text{min}}(x)| < 10^{-30} \\
     \cos^4(x-b) \cdot \dfrac{I_{\text{max}}(x)}{I_{\text{min}}(x)} & \text{otherwise}
   \end{cases}

.. math::

   x_{\text{opt}} = \arg\max_x \; \text{objective}_{\text{CSCR}}(x), \qquad
   e_{\text{proj}} = |e \cos(x_{\text{opt}} - b)|

.. math::

   l = \sqrt{A + \varepsilon} \left( \frac{I_{\text{max}}(x_{\text{opt}})}{I_{\text{min}}(x_{\text{opt}}) + \varepsilon} \right)^{0.25}

.. literalinclude:: ../src/footprint_attributes/eccentricity.py
   :pyobject: optimise_cscr
   :language: python

Feeder geometry (:mod:`footprint_attributes.geometry`)
------------------------------------------------------------

These quantities are not reported directly but feed into the parameters
above.

**Centre of mass / centre of stiffness** (hollow-box assumption):

.. math::

   \mathbf{cm}_i = \frac{A_i\,\mathbf{c}_{\text{slab},i} + w_i\,\mathbf{c}_{\text{wall},i}}{A_i + w_i + \varepsilon},
   \qquad w_i = P_i \cdot h_{\text{wall}}, \qquad \mathbf{cs}_i = \mathbf{c}_{\text{wall},i}

:math:`\mathbf{c}_{\text{slab},i}` is the polygon (area) centroid,
:math:`\mathbf{c}_{\text{wall},i}` the boundary (perimeter) centroid,
:math:`P_i` the perimeter, and :math:`h_{\text{wall}} = 3.0` m the default
wall height (cancels out of CM for identical additional storeys, sets the
relative weight of wall vs. slab mass for one storey). Degenerate (zero
area & perimeter) geometries fall back to the slab centroid.

.. literalinclude:: ../src/footprint_attributes/geometry.py
   :pyobject: centre_of_mass_and_stiffness
   :language: python

**Principal second moment of area** -- exact closed form per polygon about
its own centroid, via the shoelace/Green's-theorem identity over polygon
edges:

.. math::

   I_{xx} = \frac{1}{12}\sum_k \text{cross}_k\,(y_k^2 + y_k y_{k+1} + y_{k+1}^2), \qquad
   I_{yy} = \frac{1}{12}\sum_k \text{cross}_k\,(x_k^2 + x_k x_{k+1} + x_{k+1}^2)

.. math::

   I_{xy} = -\frac{1}{24}\sum_k \text{cross}_k\,(x_k y_{k+1} + 2x_k y_k + 2x_{k+1}y_{k+1} + x_{k+1}y_k)

where :math:`\text{cross}_k = x_k y_{k+1} - x_{k+1} y_k` (sign-canonicalised
to positive orientation). :math:`I_1 \ge I_2` and eigenvectors
:math:`\text{dir}_1, \text{dir}_2` come from eigendecomposition of
:math:`\begin{pmatrix}I_{xx} & I_{xy} \\ I_{xy} & I_{yy}\end{pmatrix}`, each
eigenvector's sign canonicalised so its largest-magnitude component is
positive.

.. literalinclude:: ../src/footprint_attributes/geometry.py
   :pyobject: calc_principal_inertia
   :language: python

**GNDT setback construction**: for every disconnected piece of
``hull(filled).difference(filled)``:

.. math::

   \text{ratio}_1 = \frac{b_1}{L_1 + \varepsilon}, \qquad \text{ratio}_2 = \frac{b_2}{L_2 + \varepsilon}

where :math:`b_1, b_2` are that piece's own circumscribed-rectangle extents
along :math:`\text{dir}_1, \text{dir}_2`. The dominant piece/configuration
per building is the one maximising :math:`\min(\text{ratio}_1,
\text{ratio}_2)` across pieces; :math:`b` is whichever of :math:`b_1/b_2` is
the binding (smaller-ratio) one, and :math:`c` is the length of the
intersection of the solid footprint with a line through the setback piece's
centroid, cast perpendicular to :math:`b`'s own direction. Pieces with
compactness :math:`1 - A_{\text{piece}}/A_{\text{hull}} \le 0.001` are
discarded (treated as convex / no setback).

.. literalinclude:: ../src/footprint_attributes/geometry.py
   :pyobject: setback_gndt_metrics
   :language: python

**GNDT inscribed-circle "a" construction**: the largest circle
:math:`(c_x, c_y, r)` fitting inside the polygon is found, then boundary
points within a tolerance of the circle are projected onto
:math:`\text{dir}_1, \text{dir}_2`:

.. math::

   a_2 = \max(P \cdot \text{dir}_1) - \min(P \cdot \text{dir}_1), \qquad
   a_1 = \max(P \cdot \text{dir}_2) - \min(P \cdot \text{dir}_2)

(:math:`a_1` pairs with :math:`L_1`, measured along :math:`\text{dir}_2`;
:math:`a_2` pairs with :math:`L_2`, measured along :math:`\text{dir}_1`). If
there are at most two tangent points (plain rectangle case),
:math:`a_1 = a_2 = 2r`.

.. literalinclude:: ../src/footprint_attributes/geometry.py
   :pyobject: max_inscribed_circle
   :language: python

.. literalinclude:: ../src/footprint_attributes/geometry.py
   :pyobject: main_element_a_lengths
   :language: python

**Polar second moment of area** (shoelace, ring-by-ring, holes subtracted
via the parallel-axis theorem):

.. math::

   I_z^{\text{ring}} = \left| \frac{1}{12}\sum_k \text{cross}_k \,
   (x_k^2 + x_kx_{k+1} + x_{k+1}^2 + y_k^2 + y_ky_{k+1} + y_{k+1}^2) \right|

.. literalinclude:: ../src/footprint_attributes/geometry.py
   :pyobject: ring_inertia_z
   :language: python

.. literalinclude:: ../src/footprint_attributes/geometry.py
   :pyobject: calc_inertia_z
   :language: python

**Undirected-axis angle fold** (used by ``ASCE7_parallelityAngle``):

.. math::

   \theta = \arccos\bigl(|\hat v_0 \cdot \hat v_1|\bigr) \in [0, \pi/2]

Same :math:`[0, \pi/2]` fold family as ``contact_angle`` and ``bearing``
above, applied here to compare undirected axes (a line has no inherent
"positive" direction) rather than to remove double-counted cancellation --
different motivation, same mechanism.

.. literalinclude:: ../src/footprint_attributes/geometry.py
   :pyobject: angle_between_0_90
   :language: python
