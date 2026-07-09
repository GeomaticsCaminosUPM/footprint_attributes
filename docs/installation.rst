Installation
============

Pinned to the `v1.0.0 release
<https://github.com/GeomaticsCaminosUPM/footprint_attributes/releases/tag/v1.0.0>`_:

.. code-block:: bash

    pip install "footprint-attributes @ git+https://github.com/GeomaticsCaminosUPM/footprint_attributes.git@v1.0.0"

Dependencies: ``geopandas``, ``shapely>=2.0``, ``numpy``, ``pandas``,
``scipy``, ``scikit-learn``, ``statsmodels``, ``tabulate``.

Optional, for the example notebooks: ``matplotlib``, ``folium``,
``ipykernel``, ``mapclassify``.

Quick start
-----------

.. code-block:: python

    import geopandas as gpd
    from footprint_attributes import direction, shape, position

    footprints = gpd.read_file("footprints.gpkg")

    bearing = direction.inertia(footprints)
    ec8 = shape.EC8(footprints)
    result = position(footprints, buffer=0.1)

See :doc:`examples` for full walkthroughs of every module.
