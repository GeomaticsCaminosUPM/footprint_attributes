"""
footprint_attributes.visualization
===================================

Interactive 3D map builder (MapLibre + deck.gl, client-side, no backend) for
this package's outputs: the geometric relative-position classification, the
4-category shape index, and the 15 EC8/ASCE7/GNDT-II/CSCR2010/NTC-23
shape-irregularity metrics, plus optional convex-hull / bounding-box /
inertia-axis / basic-length / contact-force-arrow geometry overlays.

Requires the ``vis`` extra::

    pip install "footprint-attributes[vis]"

Example
-------
>>> import geopandas as gpd
>>> from footprint_attributes.visualization import build_map
>>>
>>> datasets = {
...     "guatemala": gpd.read_file("guatemala_pilot_region.gpkg"),
...     "san_jose": gpd.read_file("san_jose_pilot_region.gpkg"),
... }
>>> build_map(datasets, "output/map", title="Pilot regions")
>>> # then: python -m http.server --directory output/map
"""

from __future__ import annotations

from .maps import build_map
from .overlays import build_overlays

__all__ = ["build_map", "build_overlays"]
