Installation
============

.. code-block:: bash

    pip install "footprint-attributes @ git+https://github.com/GeomaticsCaminosUPM/SeismicBuildingExposure.git"

Dependencies: ``geopandas``, ``shapely>=2.0``, ``numpy``, ``pandas``,
``scipy``, ``packaging``.

Optional, for the example notebooks: ``matplotlib``, ``ipykernel``.

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
