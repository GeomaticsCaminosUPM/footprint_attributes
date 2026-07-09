run
===

.. automodule:: footprint_attributes.runner
   :members:
   :undoc-members:
   :show-inheritance:

Complete column reference
==========================

Every column :func:`~footprint_attributes.runner.run` (or :func:`~footprint_attributes.shape.shape`) can produce, grouped by seismic code. Each parameter also gets a matching ``compliance_{column}`` column (0-100 score against the code's own limit) unless noted otherwise.

ASCE7
-----

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - Column
     - Compliance
     - Description
   * - ``ASCE7_setbackRatio``
     - yes
     - min(b1/L1, b2/L2) dual-configuration setback ratio (ASCE 7)
   * - ``ASCE7_holeRatio``
     - yes
     - Max interior hole area / filled area (ASCE 7)
   * - ``ASCE7_parallelityAngle``
     - yes
     - Angle between bounding-box sides and cardinal axes (ASCE 7)

Request all of ASCE7's columns at once with ``"ASCE7"`` in ``config["columns"]``, or ``shape.ASCE7(gdf)``.

CSCR2010
--------

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - Column
     - Compliance
     - Description
   * - ``CSCR2010_eccentricityRatio``
     - yes
     - Worst-case eccentricity / building dimension (CSCR 2010)

Request all of CSCR2010's columns at once with ``"CSCR2010"`` in ``config["columns"]``, or ``shape.CSCR2010(gdf)``.

EC8
---

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - Column
     - Compliance
     - Description
   * - ``EC8_eccentricityRatio``
     - yes
     - Ratio of worst-case eccentricity to torsional radius (EC8)
   * - ``EC8_radiusRatio``
     - yes
     - Ratio of torsional radius to radius of gyration (EC8)
   * - ``EC8_compactness``
     - yes
     - 1 – (largest convex-hull setback area / footprint area) (EC8)

Request all of EC8's columns at once with ``"EC8"`` in ``config["columns"]``, or ``shape.EC8(gdf)``.

GNDTII
------

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - Column
     - Compliance
     - Description
   * - ``GNDTII_beta1_mainShapeSlenderness``
     - yes
     - a / L, dominant (L, a) configuration by max(L*a) (GNDTII)
   * - ``GNDTII_beta2_setbackRatio``
     - yes
     - min(b1/L1, b2/L2) dual-configuration setback ratio (GNDTII)
   * - ``GNDTII_beta4_eccentricityRatio``
     - yes
     - Eccentricity / dominant-configuration a (GNDTII)
   * - ``GNDTII_beta6_setbackSlenderness``
     - yes
     - Protrusion depth c / winning-configuration setback width b (GNDTII)

Request all of GNDTII's columns at once with ``"GNDTII"`` in ``config["columns"]``, or ``shape.GNDTII(gdf)``.

NTC23
-----

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - Column
     - Compliance
     - Description
   * - ``NTC23_setbackRatio``
     - yes
     - min(b1/L1, b2/L2) dual-configuration setback ratio (NTC-23)
   * - ``NTC23_holeRatio``
     - yes
     - Worst hole's own-MBB width / through-centroid footprint length (NTC-23)

Request all of NTC23's columns at once with ``"NTC23"`` in ``config["columns"]``, or ``shape.NTC23(gdf)``.

Plan/vertical slenderness
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - Column
     - Compliance
     - Description
   * - ``slenderness_inertia``
     - yes
     - Plan slenderness via the inertia axis convention (EC8 limit)
   * - ``slenderness_bbox``
     - yes
     - Plan slenderness via the bbox axis convention (EC8 limit)

Request both with ``"slenderness"``, or ``shape.slenderness(gdf)``.

Code-independent shape indices
-------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Column
     - Description
   * - ``convex_hull_irregularity``
     - (hull_area - footprint_area) / footprint_area -- 0 = convex
   * - ``inertia_circle_ratio``
     - I_z(equal-area circle) / I_z(footprint)
   * - ``polsby_popper``
     - 4π A / P² -- compactness, 1.0 = perfect circle

No compliance columns (these are not tied to a specific code's limit).

Position pipeline
------------------

Requesting ``"position"`` (or any one of the five columns below, which always triggers the same atomic pipeline run) adds:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Column
     - Description
   * - ``contact_force``
     - Magnitude of the net contact-force resultant from touching neighbours
   * - ``contact_confinementRatio``
     - How evenly touching neighbours surround the building (0-1)
   * - ``contact_angularAcc``
     - Angular acceleration from an uneven (off-centre) contact-force distribution
   * - ``contact_angle``
     - Direction of the net contact-force resultant
   * - ``contact_height``
     - Height used for the contact-force computation
   * - ``relativePosition``
     - Categorical class: isolated / lateral / corner / confined / torque

Direction
---------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Column
     - Description
   * - ``bearing``
     - Building orientation, degrees clockwise from North, range [-90, 90] (``direction.inertia``)

