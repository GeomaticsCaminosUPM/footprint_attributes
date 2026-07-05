"""Idealised building shapes with known-by-construction properties.

Plain functions (not pytest fixtures), so they can be shared between the
pytest suite (via thin wrappers in ``conftest.py``) and the example
notebooks under ``examples/`` -- both exercise the exact same synthetic
geometries, and the notebooks import these rather than re-defining their
own copies.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
from shapely.affinity import rotate
from shapely.geometry import Polygon, box

CRS = "EPSG:32630"  # arbitrary projected CRS (metres)


def gdf_of(*geoms, crs: str = CRS) -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(geometry=list(geoms), crs=crs)


# ─────────────────────────────────────────────────────────────────────────────
# Direction / shape polygons
# ─────────────────────────────────────────────────────────────────────────────


def square_polygon(side: float = 10.0) -> Polygon:
    return box(0, 0, side, side)


def rect_polygon(angle_deg: float = 0.0, w: float = 20.0, h: float = 10.0) -> Polygon:
    """A w x h rectangle, optionally rotated angle_deg (CCW) about its centroid."""
    poly = box(0, 0, w, h)
    if angle_deg:
        poly = rotate(poly, angle_deg, origin="centroid")
    return poly


def l_shape_polygon() -> Polygon:
    """20x20 square minus a 10x10 corner notch -> setback ratio = 10/20 = 0.5."""
    return box(0, 0, 20, 20).difference(box(10, 10, 20, 20))


def asymmetric_l_shape_polygon() -> Polygon:
    """30x20 rectangle minus a 20x8 notch (L1 != L2, notch off-centre and
    not symmetric under any axis swap) -- catches bugs a symmetric L-shape
    (L1==L2) would hide, e.g. confusing c/L1 with min(b1/L1, b2/L2)."""
    return box(0, 0, 30, 20).difference(box(10, 12, 30, 20))


def t_shape_polygon() -> Polygon:
    """A T-shape: 4m-wide stem + 20m-wide top bar -> 2 symmetric setbacks."""
    return box(8, 0, 12, 15).union(box(0, 15, 20, 20))


def x_shape_polygon() -> Polygon:
    """An asymmetric cross: a 30m vertical bar (one part of the X longer
    than the other) crossing a 20m horizontal bar -> 4 setbacks."""
    return box(8, 0, 12, 30).union(box(0, 10, 20, 20))


def square_with_hole_polygon() -> Polygon:
    """20x20 square with a 5x5 central hole -> hole ratio = 25/400 = 0.0625."""
    return box(0, 0, 20, 20).difference(box(7.5, 7.5, 12.5, 12.5))


def rotated_hole_building_polygon() -> Polygon:
    """30x20 rectangle with an 8x4 rectangular hole rotated 25deg -- the
    hole is deliberately NOT parallel to the building's own sides."""
    building = box(0, 0, 30, 20)
    hole = rotate(box(11, 7, 19, 11), 25, origin="centroid")
    return building.difference(hole)


def circle_polygon(r: float = 10.0, n: int = 64) -> Polygon:
    """~regular n-gon approximating a circle of radius r."""
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    coords = np.column_stack([r * np.cos(theta), r * np.sin(theta)])
    return Polygon(coords)


def thin_cross_polygon() -> Polygon:
    """A plus-sign / cross shape: highly non-convex, low polsby-popper."""
    return box(4, 0, 6, 10).union(box(0, 4, 10, 6))


# ─────────────────────────────────────────────────────────────────────────────
# Position (contact-force) scenarios -- row 0 is the "building of interest"
# ─────────────────────────────────────────────────────────────────────────────


def isolated_building_gdf() -> gpd.GeoDataFrame:
    return gdf_of(box(0, 0, 10, 10))


def two_isolated_buildings_gdf() -> gpd.GeoDataFrame:
    """Two buildings far apart -> both isolated."""
    return gdf_of(box(0, 0, 10, 10), box(1000, 1000, 1010, 1010))


def lateral_pair_gdf() -> gpd.GeoDataFrame:
    """Two 10x10 squares sharing one full 10m wall -> both 'lateral'."""
    return gdf_of(box(0, 0, 10, 10), box(10, 0, 20, 10))


def corner_triplet_gdf() -> gpd.GeoDataFrame:
    """Centre square touched on two perpendicular sides -> 'corner'."""
    center = box(10, 10, 20, 20)
    right = box(20, 10, 30, 20)
    top = box(10, 20, 20, 30)
    return gdf_of(center, right, top)


def confined_quartet_gdf() -> gpd.GeoDataFrame:
    """Centre square touched on 3 of its 4 sides (N, S, E; W left open) ->
    'confined'. Both members of one opposite pair (N/S) touch -- which is
    what confinement actually requires -- so the 4th side is deliberately
    left open to confirm the classifier doesn't need every side occupied,
    just one fully cancelling opposite pair."""
    c = box(10, 10, 20, 20)
    n = box(10, 20, 20, 30)
    s = box(10, 0, 20, 10)
    e = box(20, 10, 30, 20)
    return gdf_of(c, n, s, e)


def torque_triplet_gdf(frac: float = 0.34) -> gpd.GeoDataFrame:
    """Centre 20x10 rectangle (slenderness 2, matching POSITION_DEFAULTS'
    minAngularAcc derivation) touched on opposite long sides by structures
    covering ~34% of each side, positioned at OPPOSITE ends (top-left and
    bottom-right) to maximise net torque -- the two forces nearly cancel
    (reads as 'confined') but their offset produces enough angular
    acceleration to be upgraded to 'torque'."""
    center = box(0, 0, 20, 10)
    top = box(0, 10, 20 * frac, 15)
    bottom = box(20 - 20 * frac, -5, 20, 0)
    return gdf_of(center, top, bottom)
