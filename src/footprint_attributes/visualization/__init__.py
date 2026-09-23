"""
footprint_attributes.visualization
===================================

Two complementary map builders for this package's outputs (the geometric
relative-position classification, the 4-category shape index, and the 15
EC8/ASCE7/GNDT-II/CSCR2010/NTC-23 shape-irregularity metrics):

- :func:`build_map` -- an interactive 3D map (MapLibre + deck.gl,
  client-side, no backend), buildings extruded by height, with an
  auto-playing showcase tour and optional convex-hull / bounding-box /
  inertia-axis / basic-length / contact-force-arrow geometry overlays.
- :func:`build_comparison_map` -- a 2D FancyFolium map (folium/Leaflet)
  with a dataset switcher and one flat choropleth layer per column, using
  FancyFolium's own layer-control panel and popups.

Requires the ``visualization`` extra::

    pip install "footprint-attributes[visualization]"

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

from .comparison_map import build_comparison_map
from .maps import build_map
from .overlays import build_overlays

__all__ = ["build_map", "build_comparison_map", "build_overlays"]
