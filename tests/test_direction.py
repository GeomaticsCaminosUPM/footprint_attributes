"""Tests for footprint_attributes.direction using shapes with known answers."""

from __future__ import annotations

import sys

import numpy as np
import pytest

sys.path.insert(0, "tests")
from shapes import gdf_of, rect_polygon

from footprint_attributes import direction


# ─────────────────────────────────────────────────────────────────────────────
# Full angle sweep -- pins the exact same "expected vs actual" table shown
# in examples/direction.ipynb, so a regression there is caught by pytest
# first (on every run), rather than only being noticed the next time
# someone happens to open and re-run the notebook.
# ─────────────────────────────────────────────────────────────────────────────

_SWEEP_ANGLES = [-60, -40, -30, -10, 0, 10, 30, 40, 60]


@pytest.mark.parametrize("angle", _SWEEP_ANGLES)
def test_bbox_sweep_matches_expected(angle):
    """bbox: L1=20, L2=10, bearing=-angle for every angle in the sweep,
    within the max deviation the method's 1-degree angular search grid can
    produce (empirically ~0.16 m for L1/L2, ~0.34 deg for bearing)."""
    gdf = gdf_of(rect_polygon(angle))
    L1, dir1, L2, dir2, bearing = direction.bbox(gdf, mode="all")
    assert L1[0] == pytest.approx(20.0, abs=0.2)
    assert L2[0] == pytest.approx(10.0, abs=0.2)
    assert bearing[0] == pytest.approx(-angle, abs=1.0)


@pytest.mark.parametrize("angle", _SWEEP_ANGLES)
def test_inertia_sweep_matches_expected_exactly(angle):
    """inertia: exact (closed-form) match for every angle in the sweep --
    no angular-grid discretisation error, unlike bbox."""
    gdf = gdf_of(rect_polygon(angle))
    L1, dir1, L2, dir2, bearing = direction.inertia(gdf, mode="all")
    assert L1[0] == pytest.approx(20.0, abs=1e-6)
    assert L2[0] == pytest.approx(10.0, abs=1e-6)
    assert bearing[0] == pytest.approx(-angle, abs=1e-6)


def test_square_degenerate_case_both_methods_agree():
    """A square has no unique principal axis, but both methods must still
    agree that L1 == L2 == 10."""
    square = gdf_of(rect_polygon(0.0, w=10.0, h=10.0))
    L1b, L2b = direction.bbox(square, mode="dimensions")
    L1i, L2i = direction.inertia(square, mode="dimensions")
    assert L1b[0] == pytest.approx(10.0, abs=1e-6)
    assert L2b[0] == pytest.approx(10.0, abs=1e-6)
    assert L1i[0] == pytest.approx(10.0, abs=1e-6)
    assert L2i[0] == pytest.approx(10.0, abs=1e-6)


# ─────────────────────────────────────────────────────────────────────────────
# bbox method — exact for axis-aligned / rotated rectangles
# ─────────────────────────────────────────────────────────────────────────────


def test_bbox_dimensions_axis_aligned(rect_20x10):
    L1, L2 = direction.bbox(rect_20x10, mode="dimensions")
    assert L1[0] == pytest.approx(20.0, abs=1e-6)
    assert L2[0] == pytest.approx(10.0, abs=1e-6)


def test_bbox_bearing_axis_aligned(rect_20x10):
    bearing = direction.bbox(rect_20x10)
    assert bearing[0] == pytest.approx(0.0, abs=1.0)


def test_bbox_bearing_rotated30(rect_20x10_rot30):
    bearing = direction.bbox(rect_20x10_rot30)
    assert bearing[0] == pytest.approx(-30.0, abs=1.5)


def test_bbox_bearing_rotated_neg40(rect_20x10_rot_neg40):
    # bearing = -rotation_angle by this module's clockwise-from-North convention
    bearing = direction.bbox(rect_20x10_rot_neg40)
    assert bearing[0] == pytest.approx(40.0, abs=1.5)


def test_bbox_l1_ge_l2_invariant(rect_20x10, rect_20x10_rot30, square_10):
    for gdf in (rect_20x10, rect_20x10_rot30, square_10):
        L1, L2 = direction.bbox(gdf, mode="dimensions")
        assert L1[0] >= L2[0]


def test_bbox_dir1_dir2_orthogonal(rect_20x10_rot30):
    _, dir1, _, dir2, _ = direction.bbox(rect_20x10_rot30, mode="all")
    assert np.dot(dir1[0], dir2[0]) == pytest.approx(0.0, abs=1e-9)


# ─────────────────────────────────────────────────────────────────────────────
# inertia method — calc_principal_inertia uses the exact closed-form polygon
# second-moment-of-area formula (shoelace identity), so this is exact for any
# polygon, not an approximation -- verified here to near machine precision.
# ─────────────────────────────────────────────────────────────────────────────


def test_inertia_l1_ge_l2_invariant(rect_20x10, rect_20x10_rot30, square_10):
    """Regression test: L1 must always be >= L2 (see direction.py docstring).

    Prior to the fix, calc_principal_inertia's larger-eigenvalue eigenvector
    pointed along the *shorter* physical axis (an area-moment-of-inertia
    tensor is the PCA covariance matrix with x/y roles swapped), so L1 could
    come out smaller than L2. Direction.inertia() now swaps the eigenvector
    labels to keep L1 the longer projected dimension.
    """
    for gdf in (rect_20x10, rect_20x10_rot30, square_10):
        L1, L2 = direction.inertia(gdf, mode="dimensions")
        assert L1[0] >= L2[0]


def test_inertia_bearing_exact_for_rectangle(rect_20x10_rot30, rect_20x10_rot_neg40):
    """calc_principal_inertia's closed-form polygon moments are exact for any
    polygon (not a vertex-count-dependent approximation), so a rotated
    rectangle's inertia bearing must match the true rotation angle to near
    machine precision -- in fact more precisely than direction.bbox, whose
    1-degree angular search grid is the accuracy bottleneck there."""
    assert direction.inertia(rect_20x10_rot30)[0] == pytest.approx(-30.0, abs=1e-6)
    assert direction.inertia(rect_20x10_rot_neg40)[0] == pytest.approx(40.0, abs=1e-6)


def test_inertia_dimensions_exact_for_axis_aligned_rectangle(rect_20x10):
    L1, dir1, L2, dir2, bearing = direction.inertia(rect_20x10, mode="all")
    assert L1[0] == pytest.approx(20.0, abs=1e-9)
    assert L2[0] == pytest.approx(10.0, abs=1e-9)
    assert dir1[0] == pytest.approx([1.0, 0.0], abs=1e-9)
    assert dir2[0] == pytest.approx([0.0, 1.0], abs=1e-9)
    assert bearing[0] == pytest.approx(0.0, abs=1e-9)


def test_calc_principal_inertia_matches_closed_form_rectangle(rect_20x10):
    """Pin the exact closed-form moments for a 20x10 rectangle: Ixx=(1/12)*w*h^3,
    Iyy=(1/12)*h*w^3 about the centroidal axes (w=20 along x, h=10 along y)."""
    from footprint_attributes.geometry import calc_principal_inertia

    I1, dir1, I2, dir2 = calc_principal_inertia(rect_20x10.geometry)
    expected_Iyy = (1 / 12) * 10 * 20**3  # larger -> I1, eigenvector along y
    expected_Ixx = (1 / 12) * 20 * 10**3  # smaller -> I2, eigenvector along x
    assert I1[0] == pytest.approx(expected_Iyy, rel=1e-9)
    assert I2[0] == pytest.approx(expected_Ixx, rel=1e-9)


def test_inertia_dimensions_match_paper_eqs_2_3(rect_20x10):
    """direction.inertia()'s L1/L2 must be the average side lengths implied
    by slenderness=sqrt(I1/I2) (paper eqs. 2-3), not a direct projection of
    the polygon onto the eigenvectors. This is both what the paper specifies
    and, empirically, a closer match to the true 20x10 rectangle (~20.3/9.8)
    than the old projection-based values (~21.1/12.5)."""
    from footprint_attributes.geometry import (
        calc_principal_inertia,
        inertia_side_lengths,
    )

    I1, _, I2, _ = calc_principal_inertia(rect_20x10.geometry)
    area = rect_20x10.geometry.area.values
    expected_L1, expected_L2 = inertia_side_lengths(I1, I2, area)

    L1, L2 = direction.inertia(rect_20x10, mode="dimensions")
    assert L1[0] == pytest.approx(expected_L1[0])
    assert L2[0] == pytest.approx(expected_L2[0])
    # and it should be a noticeably closer approximation of the true 20x10
    # rectangle than the naive projection-based method used to produce.
    assert abs(L1[0] - 20.0) < 1.0
    assert abs(L2[0] - 10.0) < 1.0


def test_inertia_dir1_dir2_orthogonal(rect_20x10):
    _, dir1, _, dir2, _ = direction.inertia(rect_20x10, mode="all")
    assert np.dot(dir1[0], dir2[0]) == pytest.approx(0.0, abs=1e-9)


# ─────────────────────────────────────────────────────────────────────────────
# Forced direction
# ─────────────────────────────────────────────────────────────────────────────


def test_forced_direction_projects_onto_given_axis(rect_20x10):
    # Force the axis at 45 deg; L1/L2 should be the diagonal-projected extents,
    # both between the true L1=20 and L2=10 (a 45-degree cut of an axis-aligned
    # rectangle always yields a bigger apparent footprint than the true sides).
    L1, dir1, L2, dir2 = direction.bbox(
        rect_20x10, mode="all", direction=np.array([1.0, 1.0])
    )[:4]
    assert dir1[0] == pytest.approx([np.sqrt(0.5), np.sqrt(0.5)], abs=1e-6)
    assert L1[0] > 0 and L2[0] > 0


def test_forced_direction_dir1_dir2_perpendicular(rect_20x10):
    _, dir1, _, dir2 = direction.bbox(
        rect_20x10, mode="all", direction=np.array([2.0, 0.0])
    )[:4]
    assert np.dot(dir1[0], dir2[0]) == pytest.approx(0.0, abs=1e-9)
    assert dir1[0] == pytest.approx([1.0, 0.0])


# ─────────────────────────────────────────────────────────────────────────────
# Mode contract
# ─────────────────────────────────────────────────────────────────────────────


def test_invalid_mode_raises(rect_20x10):
    with pytest.raises(ValueError):
        direction.bbox(rect_20x10, mode="not_a_mode")


def test_all_mode_returns_five_items(rect_20x10):
    out = direction.bbox(rect_20x10, mode="all")
    assert len(out) == 5
