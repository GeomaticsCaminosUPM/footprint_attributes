Examples
========

Each notebook loads the real sample datasets under ``examples/data/``,
checks every computation on hand-built shapes with a known answer, then runs
it on the three pilot regions. Every notebook page below includes both of
its interactive maps: the 2-D Folium map with one layer per column, and the
3-D MapLibre + deck.gl viewer, opened on that notebook's topic.

The same 3-D viewer can be opened on any attribute or overlay from its URL
(``?dataset=``, ``?attribute=``, ``?overlays=``, ``?zoom=``, ``?tour=``);
see ``docs/docs_maps_and_plots.py``. For example, `block position in
Guatemala City with the contact-force arrows
<_static/maps/index.html?dataset=guatemala&attribute=blockPosition&overlays=position_arrows&zoom=17>`_
or `the shape index in San José with convex hulls
<_static/maps/index.html?dataset=san_jose&attribute=shape_index&overlays=convex_hull>`_
opens full-screen.

Interactive map
----------------

Runs entirely client-side (MapLibre + deck.gl, no backend) and computes
nothing itself -- everything it shows is precomputed from the pilot-region
footprints by :func:`footprint_attributes.visualization.build_map` (see
``docs/docs_maps_and_plots.py``), using this package's own
:func:`~footprint_attributes.shape.shape` and
:func:`~footprint_attributes.position.position`. Requires the ``visualization``
extra: ``pip install "footprint-attributes[visualization]"``.

Switch dataset with the top-left dropdown (Guatemala Zona 10, San José Mata
Redonda, Santo Domingo Ensanche Quisquella); color by block position,
the 4-category shape index, or any of the 15 EC8/ASCE7/GNDT-II/CSCR2010/
NTC-23 shape-irregularity metrics described in :doc:`formulas`; click a
building for its full value breakdown. The "Overlays" checkboxes draw
per-building geometry on top of the map: convex hull, minimum bounding box,
inertia axis, the L1/L2/a1/a2/b/c basic-length dimensions, and the net
contact-force resultant arrow. The auto-playing tour (▶ button, top right)
cycles block position, then the shape index, then each shape metric in
turn, one every 10s.

.. raw:: html

   <iframe src="_static/maps/index.html" width="100%" height="640" style="border:1px solid #cbd5e0;border-radius:6px;" loading="lazy"></iframe>

.. image:: figures/interactive_map_block_position.jpg
   :width: 49%
   :alt: Interactive map colored by block position

.. image:: figures/interactive_map_shape_index.jpg
   :width: 49%
   :alt: Interactive map colored by shape index

.. image:: figures/interactive_map_ec8_eccentricity.jpg
   :width: 49%
   :alt: Interactive map colored by EC8 eccentricity ratio, with the norm-exceedance chart

.. image:: figures/interactive_map_overlays.jpg
   :width: 49%
   :alt: Interactive map with all five geometry overlays enabled

Build it yourself against this repo's own pilot-region datasets with
``examples/generate_interactive_map.py``, or against your own footprints
with :func:`footprint_attributes.visualization.build_map` directly.

Static maps
-----------

Every static map and plot in these docs is generated from the same
pilot-region footprints by ``docs/docs_maps_and_plots.py`` (see
:doc:`concepts` and :doc:`formulas` for all of them). Two examples -- the
contact forces behind ``blockPosition`` and the convex-hull setback
pieces behind the shape indices, drawn on real buildings:

.. image:: maps/detail_contact_forces.jpg
   :width: 100%
   :alt: Close-up map with per-wall contact forces

.. image:: maps/detail_convex_hull.jpg
   :width: 100%
   :alt: Close-up map with convex hulls and setback pieces

Notebooks
---------

.. toctree::
   :maxdepth: 1

   examples/direction.ipynb
   examples/position.ipynb
   examples/shape.ipynb
   examples/building_sizes.ipynb
   examples/run.ipynb
