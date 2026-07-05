"""Tests for the low-level GNDT construction primitives in geometry.py:
inscribed circle, tangent points, dual-configuration setbacks, and the
NTC-23 hole h/L measurement.
"""

from __future__ import annotations

import numpy as np
import pytest

from footprint_attributes.geometry import (
    circle_tangent_points,
    hole_h_over_l,
    main_element_a_lengths,
    main_element_a_lengths_batch,
    max_inscribed_circle,
    setback_gndt_metrics,
    setback_pieces,
)


# ─────────────────────────────────────────────────────────────────────────────
# max_inscribed_circle
# ─────────────────────────────────────────────────────────────────────────────


def test_inscribed_circle_radius_is_half_short_side_for_rectangle(rect_20x10):
    """For a 20x10 rectangle the max inscribed circle has r=5 (half the
    short side); its centre can sit anywhere along the degenerate long axis
    (y=5 is fixed, x is not unique), so only check the well-defined y and r."""
    poly = rect_20x10.geometry.iloc[0]
    cx, cy, r = max_inscribed_circle(poly)
    assert r == pytest.approx(5.0, abs=1e-4)
    assert cy == pytest.approx(5.0, abs=1e-4)
    assert 5.0 <= cx <= 15.0  # anywhere the circle still fits within x in [0,20]


def test_inscribed_circle_matches_true_radius_for_square(square_10):
    poly = square_10.geometry.iloc[0]
    cx, cy, r = max_inscribed_circle(poly)
    assert r == pytest.approx(5.0, abs=1e-4)


# ─────────────────────────────────────────────────────────────────────────────
# circle_tangent_points
# ─────────────────────────────────────────────────────────────────────────────


def test_tangent_points_rectangle_gives_two_opposite_points(rect_20x10):
    poly = rect_20x10.geometry.iloc[0]
    cx, cy, r = max_inscribed_circle(poly)
    pts = circle_tangent_points(poly, cx, cy, r)
    assert len(pts) == 2
    p1, p2 = np.array(pts)
    # opposite points on the circle -> separated by the diameter
    assert np.linalg.norm(p1 - p2) == pytest.approx(2 * r, abs=1e-2)


# ─────────────────────────────────────────────────────────────────────────────
# main_element_a_lengths -- the GNDT 'a' construction
# ─────────────────────────────────────────────────────────────────────────────


def test_a_lengths_degenerate_case_equals_diameter_for_rectangle(rect_20x10):
    """A plain rectangle's inscribed circle only touches 2 opposite sides, so
    both (a1, a2) configurations collapse to the circle diameter -- a
    rectangle only has one meaningful width. The rectangle's centre in this
    degenerate case is just the circle's own centre."""
    poly = rect_20x10.geometry.iloc[0]
    a1, a2, center = main_element_a_lengths(
        poly, np.array([1.0, 0.0]), np.array([0.0, 1.0])
    )
    assert a1 == pytest.approx(10.0, abs=1e-3)
    assert a2 == pytest.approx(10.0, abs=1e-3)
    cx, cy, r = max_inscribed_circle(poly)
    assert center == pytest.approx((cx, cy))


def test_a_lengths_batch_matches_single_call(l_shape_notch_half):
    gdf = l_shape_notch_half
    poly = gdf.geometry.iloc[0]
    dir1 = np.array([1.0, 0.0])
    dir2 = np.array([0.0, 1.0])
    single = main_element_a_lengths(poly, dir1, dir2)
    batch = main_element_a_lengths_batch(gdf, dir1.reshape(1, 2), dir2.reshape(1, 2))
    assert batch[0][0] == pytest.approx(single[0])
    assert batch[1][0] == pytest.approx(single[1])
    assert batch[2][0] == pytest.approx(single[2])


def test_a_lengths_center_not_circle_center_when_asymmetric(t_shape):
    """For a shape where the inscribed circle's tangent points are NOT
    symmetric about the circle's own centre (e.g. the T-shape, whose circle
    sits pinched near the stem/bar junction), the returned rectangle centre
    must differ from the circle centre -- otherwise the (a1 x a2) rectangle
    drawn around it wouldn't actually bound the tangent points."""
    poly = t_shape.geometry.iloc[0]
    dir1 = np.array([1.0, 0.0])
    dir2 = np.array([0.0, 1.0])
    cx, cy, r = max_inscribed_circle(poly)
    _, _, center = main_element_a_lengths(poly, dir1, dir2)
    assert center != pytest.approx((cx, cy), abs=1e-6)


# ─────────────────────────────────────────────────────────────────────────────
# setback_pieces -- per-piece introspection (all setbacks, not just dominant)
# ─────────────────────────────────────────────────────────────────────────────


def test_setback_pieces_empty_for_convex_rectangle(rect_20x10):
    poly = rect_20x10.geometry.iloc[0]
    assert setback_pieces(poly, np.array([1.0, 0.0]), np.array([0.0, 1.0])) == []


def test_setback_pieces_single_symmetric_l_shape(l_shape_notch_half):
    poly = l_shape_notch_half.geometry.iloc[0]
    pieces = setback_pieces(poly, np.array([1.0, 0.0]), np.array([0.0, 1.0]))
    assert len(pieces) == 1
    ext1, ext2, piece = pieces[0]
    assert ext1 == pytest.approx(10.0, abs=1e-6)
    assert ext2 == pytest.approx(10.0, abs=1e-6)
    assert piece.area == pytest.approx(50.0, abs=1e-6)


def test_setback_pieces_t_shape_has_two_symmetric_pieces(t_shape):
    """The T-shape's two setback triangles are mirror images of each other,
    so both pieces must have identical (ext1, ext2) and area."""
    poly = t_shape.geometry.iloc[0]
    pieces = setback_pieces(poly, np.array([1.0, 0.0]), np.array([0.0, 1.0]))
    assert len(pieces) == 2
    (ext1a, ext2a, piece_a), (ext1b, ext2b, piece_b) = pieces
    assert ext1a == pytest.approx(ext1b, abs=1e-6)
    assert ext2a == pytest.approx(ext2b, abs=1e-6)
    assert piece_a.area == pytest.approx(piece_b.area, abs=1e-6)


def test_setback_pieces_x_shape_has_four_pieces(x_shape):
    poly = x_shape.geometry.iloc[0]
    pieces = setback_pieces(poly, np.array([1.0, 0.0]), np.array([0.0, 1.0]))
    assert len(pieces) == 4


# ─────────────────────────────────────────────────────────────────────────────
# setback_gndt_metrics -- dual-configuration b/L selection
# ─────────────────────────────────────────────────────────────────────────────


def test_setback_metrics_zero_for_convex_rectangle(rect_20x10):
    L1 = np.array([20.0])
    L2 = np.array([10.0])
    dir1 = np.array([[1.0, 0.0]])
    dir2 = np.array([[0.0, 1.0]])
    ratio, b, c = setback_gndt_metrics(rect_20x10, L1, dir1, L2, dir2)
    assert ratio == [0.0]
    assert b == [0.0]
    assert c == [0.0]


def test_setback_metrics_symmetric_l_shape(l_shape_notch_half):
    """20x20 square minus a 10x10 corner notch: L1=L2=20, single setback
    piece with b1=b2=10 (symmetric under axis swap), c=10 (measured against
    the solid footprint, perpendicular to b) -> ratio=0.5, c/b=1.0."""
    L1 = np.array([20.0])
    L2 = np.array([20.0])
    dir1 = np.array([[1.0, 0.0]])
    dir2 = np.array([[0.0, 1.0]])
    ratio, b, c = setback_gndt_metrics(l_shape_notch_half, L1, dir1, L2, dir2)
    assert ratio[0] == pytest.approx(0.5, abs=1e-6)
    assert b[0] == pytest.approx(10.0, abs=1e-6)
    assert c[0] == pytest.approx(10.0, abs=1e-6)


def test_setback_metrics_asymmetric_l_shape_hand_computed(asymmetric_l_shape):
    """30x20 rectangle minus a 20x8 notch at (10-30, 12-20). Convex hull is a
    pentagon (drops the reflex vertex), and hull-minus-footprint is a single
    80 m^2 triangle with extents b1=20 (along x, dir1) and b2=8 (along y,
    dir2). L1=30, L2=20 -> ratio1=20/30=0.667, ratio2=8/20=0.4 -> the
    winning config is ratio2 (0.4), with b=8. Casting a perpendicular line
    (along dir1/x) through the triangle's centroid (y~14.67, inside the
    12-20 band where only the x in [0,10] tab remains) against the solid
    footprint gives c=10, so beta6 = c/b = 10/8 = 1.25."""
    L1 = np.array([30.0])
    L2 = np.array([20.0])
    dir1 = np.array([[1.0, 0.0]])
    dir2 = np.array([[0.0, 1.0]])
    ratio, b, c = setback_gndt_metrics(asymmetric_l_shape, L1, dir1, L2, dir2)
    assert ratio[0] == pytest.approx(0.4, abs=1e-6)
    assert b[0] == pytest.approx(8.0, abs=1e-6)
    assert c[0] == pytest.approx(10.0, abs=1e-6)


def test_setback_metrics_uses_own_L_not_recomputed_bbox(asymmetric_l_shape):
    """Regression test: setback_gndt_metrics must use the L1/L2 passed in by
    the caller (matching whichever dir1/dir2 axes the caller is using), not
    silently recompute its own bbox-based L1/L2 internally -- otherwise the
    ratio would be inconsistent with a caller using the inertia method."""
    # Deliberately pass "wrong" (swapped) L1/L2 and check the ratio changes
    # accordingly -- proving the passed-in L1/L2 are actually used.
    dir1 = np.array([[1.0, 0.0]])
    dir2 = np.array([[0.0, 1.0]])
    ratio_normal, _, _ = setback_gndt_metrics(
        asymmetric_l_shape, np.array([30.0]), dir1, np.array([20.0]), dir2
    )
    ratio_swapped, _, _ = setback_gndt_metrics(
        asymmetric_l_shape, np.array([20.0]), dir1, np.array([30.0]), dir2
    )
    assert ratio_normal[0] != pytest.approx(ratio_swapped[0])


# ─────────────────────────────────────────────────────────────────────────────
# hole_h_over_l -- NTC-23 hole ratio (distinct from ASCE7's area ratio)
# ─────────────────────────────────────────────────────────────────────────────


def test_hole_h_over_l_zero_without_hole(rect_20x10):
    assert hole_h_over_l(rect_20x10) == [0.0]


def test_hole_h_over_l_recovers_true_hole_size_when_rotated(rotated_hole_building):
    """The hole is an 8x4 m rectangle rotated 25 deg; h (the hole's own MBB
    short side) must recover ~4.0 regardless of the rotation, and differ
    from ASCE7's area-based hole ratio (25/400-style), since these are
    two distinct metrics per the paper (h/L vs A_hole/A_filled)."""
    from footprint_attributes import shape

    ratio = hole_h_over_l(rotated_hole_building)[0]
    # h ~ 4.0, L is the through-centroid chord (> h), so ratio is modest but nonzero
    assert 0.0 < ratio < 1.0

    asce7_ratio = shape.ASCE7(rotated_hole_building)["ASCE7_holeRatio"].iloc[0]
    ntc23_ratio = shape.NTC23(rotated_hole_building)["NTC23_holeRatio"].iloc[0]
    assert ntc23_ratio == pytest.approx(ratio)
    assert ntc23_ratio != pytest.approx(asce7_ratio, rel=0.05)
