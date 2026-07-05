"""Pytest fixtures: thin wrappers around the plain factory functions in
``shapes.py`` (shared with the example notebooks)."""

from __future__ import annotations

import pytest

from shapes import (
    CRS,
    gdf_of,
    asymmetric_l_shape_polygon,
    circle_polygon,
    confined_quartet_gdf,
    corner_triplet_gdf,
    isolated_building_gdf,
    l_shape_polygon,
    lateral_pair_gdf,
    rect_polygon,
    rotated_hole_building_polygon,
    square_polygon,
    square_with_hole_polygon,
    thin_cross_polygon,
    torque_triplet_gdf,
    two_isolated_buildings_gdf,
    t_shape_polygon,
    x_shape_polygon,
)

__all__ = ["CRS", "gdf_of"]  # re-exported for tests that build custom shapes


# ─────────────────────────────────────────────────────────────────────────────
# Direction fixtures
# ─────────────────────────────────────────────────────────────────────────────


@pytest.fixture
def square_10():
    return gdf_of(square_polygon(10.0))


@pytest.fixture
def rect_20x10():
    """Axis-aligned rectangle: L1=20 (x), L2=10 (y), bearing = 0."""
    return gdf_of(rect_polygon(0.0))


@pytest.fixture
def rect_20x10_rot30():
    """Same rectangle rotated 30 deg CCW about its centroid."""
    return gdf_of(rect_polygon(30.0))


@pytest.fixture
def rect_20x10_rot_neg40():
    return gdf_of(rect_polygon(-40.0))


# ─────────────────────────────────────────────────────────────────────────────
# Shape fixtures
# ─────────────────────────────────────────────────────────────────────────────


@pytest.fixture
def l_shape_notch_half():
    """20x20 square minus a 10x10 corner notch -> setback ratio = 10/20 = 0.5."""
    return gdf_of(l_shape_polygon())


@pytest.fixture
def asymmetric_l_shape():
    """30x20 rectangle minus a 20x8 notch (L1 != L2, notch off-centre and
    not symmetric under any axis swap) -- catches bugs that a symmetric
    L-shape (L1==L2) would hide, e.g. confusing c/L1 with min(b1/L1,b2/L2)."""
    return gdf_of(asymmetric_l_shape_polygon())


@pytest.fixture
def t_shape():
    """A T-shape: 4m-wide stem + 20m-wide top bar -> 2 symmetric setbacks."""
    return gdf_of(t_shape_polygon())


@pytest.fixture
def x_shape():
    """An asymmetric cross: one part of the X longer than the other -> 4 setbacks."""
    return gdf_of(x_shape_polygon())


@pytest.fixture
def square_with_hole():
    """20x20 square with a 5x5 central hole -> hole ratio = 25/400 = 0.0625."""
    return gdf_of(square_with_hole_polygon())


@pytest.fixture
def rotated_hole_building():
    """30x20 rectangle with an 8x4 rectangular hole rotated 25deg -- the
    hole is deliberately NOT parallel to the building's own sides."""
    return gdf_of(rotated_hole_building_polygon())


@pytest.fixture
def circle_approx():
    """~regular 64-gon approximating a circle -> polsby_popper close to 1."""
    return gdf_of(circle_polygon(10.0, 64))


@pytest.fixture
def thin_cross():
    """A plus-sign / cross shape: highly non-convex, low polsby-popper."""
    return gdf_of(thin_cross_polygon())


# ─────────────────────────────────────────────────────────────────────────────
# Position fixtures
# ─────────────────────────────────────────────────────────────────────────────


@pytest.fixture
def isolated_building():
    return isolated_building_gdf()


@pytest.fixture
def two_isolated_buildings():
    """Two buildings far apart -> both isolated. Regression test for the
    'no contacts at all in the batch' crash."""
    return two_isolated_buildings_gdf()


@pytest.fixture
def lateral_pair():
    """Two 10x10 squares sharing one full 10 m wall -> both 'lateral'."""
    return lateral_pair_gdf()


@pytest.fixture
def corner_triplet():
    """Centre square touched on two perpendicular sides -> 'corner'."""
    return corner_triplet_gdf()


@pytest.fixture
def confined_quartet():
    """Centre square touched on 3 of its 4 sides (one opposite pair fully
    cancels) -> 'confined'."""
    return confined_quartet_gdf()


@pytest.fixture
def torque_triplet():
    """Centre 20x10 rectangle touched on opposite long sides at opposite
    ends -> upgraded from 'confined' to 'torque'."""
    return torque_triplet_gdf()
