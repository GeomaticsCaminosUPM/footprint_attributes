"""Tests for footprint_attributes.shape using shapes with known answers."""

from __future__ import annotations

import sys

import numpy as np
import pytest

sys.path.insert(0, "tests")
from shapes import (
    circle_polygon,
    gdf_of,
    l_shape_polygon,
    rect_polygon,
    square_polygon,
    square_with_hole_polygon,
    thin_cross_polygon,
)

from footprint_attributes import shape


# ─────────────────────────────────────────────────────────────────────────────
# Ideal-shape expected-value table -- pins the exact same "expected vs
# actual" checks shown in examples/shape.ipynb, so a regression there is
# caught by pytest on every run, not only when someone happens to re-run
# the notebook by hand.
# ─────────────────────────────────────────────────────────────────────────────

_IDEAL_SHAPES = {
    "square": square_polygon(),
    "rect": rect_polygon(),
    "l_shape": l_shape_polygon(),
    "with_hole": square_with_hole_polygon(),
    "circle": circle_polygon(),
    "cross": thin_cross_polygon(),
}

_EXPECTED_INDICES = {
    "square": dict(polsby_popper=np.pi / 4, convex_hull_irregularity=0.0),
    "rect": dict(convex_hull_irregularity=0.0),
    "l_shape": dict(convex_hull_irregularity=50 / 300),
    "with_hole": dict(convex_hull_irregularity=0.0),  # hull ignores the hole
    "circle": dict(
        polsby_popper=1.0, inertia_circle_ratio=1.0, convex_hull_irregularity=0.0
    ),
}

_EXPECTED_BBOX_SLENDERNESS = {
    "square": 1.0,
    "rect": 2.0,
    "l_shape": 1.0,
    "with_hole": 1.0,
    "circle": 1.0,
    "cross": 1.0,
}


@pytest.mark.parametrize(
    "name, metric, expected",
    [
        (name, metric, expected)
        for name, metrics in _EXPECTED_INDICES.items()
        for metric, expected in metrics.items()
    ],
)
def test_ideal_shape_code_independent_indices(name, metric, expected):
    gdf = gdf_of(_IDEAL_SHAPES[name])
    fn = getattr(shape, metric)
    got = fn(gdf)[0]
    assert got == pytest.approx(expected, abs=1e-3)


@pytest.mark.parametrize("name, expected", list(_EXPECTED_BBOX_SLENDERNESS.items()))
def test_ideal_shape_bbox_slenderness(name, expected):
    gdf = gdf_of(_IDEAL_SHAPES[name])
    got = shape.slenderness.bbox(gdf)[0]
    assert got == pytest.approx(expected, abs=1e-3)


def test_ideal_shape_asce7_vs_ntc23_hole_ratio_are_distinct_metrics():
    """20x20 square with a centred 5x5 hole: ASCE7's area ratio (25/400 =
    0.0625) and NTC-23's length ratio (h=5 / L=20 = 0.25) must differ --
    they are documented as two distinct metrics (area vs. length), not
    variants of the same formula."""
    gdf = gdf_of(square_with_hole_polygon())
    asce7 = shape.ASCE7(gdf)["ASCE7_holeRatio"].iloc[0]
    ntc23 = shape.NTC23(gdf)["NTC23_holeRatio"].iloc[0]
    assert asce7 == pytest.approx(25 / 400, abs=1e-3)
    assert ntc23 == pytest.approx(5 / 20, abs=1e-3)
    assert asce7 != pytest.approx(ntc23, rel=0.05)


# ─────────────────────────────────────────────────────────────────────────────
# Code-independent indices
# ─────────────────────────────────────────────────────────────────────────────


def test_polsby_popper_square(square_10):
    # 4*pi*A / P^2 = 4*pi*100 / 1600 = pi/4
    assert shape.polsby_popper(square_10)[0] == pytest.approx(np.pi / 4, abs=1e-9)


def test_polsby_popper_circle_near_one(circle_approx):
    assert shape.polsby_popper(circle_approx)[0] == pytest.approx(1.0, abs=1e-3)


def test_polsby_popper_cross_lower_than_square(thin_cross, square_10):
    assert shape.polsby_popper(thin_cross)[0] < shape.polsby_popper(square_10)[0]


def test_inertia_circle_ratio_circle_near_one(circle_approx):
    assert shape.inertia_circle_ratio(circle_approx)[0] == pytest.approx(1.0, abs=1e-2)


def test_inertia_circle_ratio_bounded(square_10, thin_cross):
    for gdf in (square_10, thin_cross):
        val = shape.inertia_circle_ratio(gdf)[0]
        assert 0.0 < val <= 1.0 + 1e-9


def test_convex_hull_irregularity_zero_for_convex_square(square_10):
    assert shape.convex_hull_irregularity(square_10)[0] == pytest.approx(0.0, abs=1e-9)


def test_convex_hull_irregularity_l_shape_matches_hand_computation(l_shape_notch_half):
    """Fixed: convex_hull_irregularity = (hull_area - footprint_area) / footprint_area.
    For the L-shape (area=300, hull=350, see the EC8 compactness hand
    computation elsewhere), this is (350-300)/300 = 50/300."""
    assert shape.convex_hull_irregularity(l_shape_notch_half)[0] == pytest.approx(
        50 / 300, abs=1e-6
    )


def test_convex_hull_irregularity_positive_for_notched_shapes(
    l_shape_notch_half, thin_cross
):
    for gdf in (l_shape_notch_half, thin_cross):
        assert shape.convex_hull_irregularity(gdf)[0] > 0.0


def test_convex_hull_irregularity_ignores_holes_not_setbacks(square_with_hole):
    """A hole makes the footprint area smaller without changing its (already
    convex) outer boundary or hull -- this index should stay 0, since holes
    are a separate concern (ASCE7/NTC23 hole ratio), not a plan-shape setback."""
    assert shape.convex_hull_irregularity(square_with_hole)[0] == pytest.approx(
        0.0, abs=1e-9
    )


# ─────────────────────────────────────────────────────────────────────────────
# Slenderness
# ─────────────────────────────────────────────────────────────────────────────


def test_slenderness_bbox_square_is_one(square_10):
    assert shape.slenderness.bbox(square_10)[0] == pytest.approx(1.0, abs=1e-6)


def test_slenderness_bbox_rect_is_two(rect_20x10):
    assert shape.slenderness.bbox(rect_20x10)[0] == pytest.approx(2.0, abs=1e-6)


def test_slenderness_ge_one(rect_20x10, square_10, l_shape_notch_half):
    for gdf in (rect_20x10, square_10, l_shape_notch_half):
        assert shape.slenderness.bbox(gdf)[0] >= 1.0 - 1e-9
        assert shape.slenderness.inertia(gdf)[0] >= 1.0 - 1e-9


def test_slenderness_inertia_exact_for_rectangle(rect_20x10):
    """Per paper eq. 2, slenderness = sqrt(I1/I2), which is exactly 2.0 for
    a true 20x10 rectangle: calc_principal_inertia now uses the exact
    closed-form polygon second-moment formula, not a vertex-count-dependent
    approximation, so this must match to near machine precision."""
    assert shape.slenderness.inertia(rect_20x10)[0] == pytest.approx(2.0, abs=1e-9)


def test_slenderness_vertical_requires_height_column(rect_20x10):
    with pytest.raises(ValueError):
        shape.slenderness.bbox(rect_20x10, vertical=True)


def test_slenderness_vertical_height_over_l2(rect_20x10):
    gdf = rect_20x10.copy()
    gdf["height"] = 30.0  # L2 = 10 -> vertical slenderness = 3
    val = shape.slenderness.bbox(gdf, vertical=True, height_column="height")[0]
    assert val == pytest.approx(3.0, abs=1e-6)


# ─────────────────────────────────────────────────────────────────────────────
# ASCE 7 / GNDTII setback & hole ratios
# ─────────────────────────────────────────────────────────────────────────────


def test_asce7_setback_ratio_zero_for_rectangle(rect_20x10):
    out = shape.ASCE7(rect_20x10)
    assert out["ASCE7_setbackRatio"].iloc[0] == pytest.approx(0.0, abs=1e-6)


def test_asce7_setback_ratio_l_shape(l_shape_notch_half):
    out = shape.ASCE7(l_shape_notch_half)
    assert out["ASCE7_setbackRatio"].iloc[0] == pytest.approx(0.5, abs=1e-3)


def test_asce7_hole_ratio_known_value(square_with_hole):
    out = shape.ASCE7(square_with_hole)
    # 5x5 hole in a 20x20 square -> 25 / 400
    assert out["ASCE7_holeRatio"].iloc[0] == pytest.approx(25 / 400, abs=1e-3)


def test_asce7_hole_ratio_zero_without_hole(rect_20x10):
    out = shape.ASCE7(rect_20x10)
    assert out["ASCE7_holeRatio"].iloc[0] == pytest.approx(0.0, abs=1e-9)


def test_gndtii_beta2_setback_matches_l_shape(l_shape_notch_half):
    out = shape.GNDTII(l_shape_notch_half)
    assert out["GNDTII_beta2_setbackRatio"].iloc[0] == pytest.approx(0.5, abs=1e-3)


def test_gndtii_beta1_slenderness_square_is_one(square_10):
    out = shape.GNDTII(square_10)
    assert out["GNDTII_beta1_mainShapeSlenderness"].iloc[0] == pytest.approx(
        1.0, abs=1e-6
    )


# ─────────────────────────────────────────────────────────────────────────────
# EC8 / CSCR2010 / GNDTII eccentricity — CM/CS from the hollow-box model
# ─────────────────────────────────────────────────────────────────────────────
#
# CM (centre of mass) is the area/perimeter-weighted average of the slab
# centroid (footprint centroid) and the wall centroid (boundary centroid),
# per the paper's eq. (1); CS (centre of stiffness) is the wall centroid
# alone. For a symmetric shape (square, axis-aligned rectangle) both slab
# and wall centroids coincide with the geometric centre, so CM == CS and
# eccentricity is genuinely zero. For an asymmetric shape (an L-notch) the
# wall mass is concentrated away from the slab centroid, so CM != CS and a
# non-zero eccentricity ratio must come out.


def test_eccentricity_zero_for_symmetric_shapes(rect_20x10, square_10):
    for gdf in (rect_20x10, square_10):
        out = shape.EC8(gdf)
        assert out["EC8_eccentricityRatio"].iloc[0] == pytest.approx(0.0, abs=1e-9)


def test_eccentricity_nonzero_for_asymmetric_l_shape(l_shape_notch_half):
    out = shape.EC8(l_shape_notch_half)
    assert out["EC8_eccentricityRatio"].iloc[0] > 0.0

    cscr = shape.CSCR2010(l_shape_notch_half)
    assert cscr["CSCR2010_eccentricityRatio"].iloc[0] > 0.0

    gndt = shape.GNDTII(l_shape_notch_half)
    assert gndt["GNDTII_beta4_eccentricityRatio"].iloc[0] > 0.0


def test_centre_of_mass_and_stiffness_matches_hand_computation(l_shape_notch_half):
    """L-shape = 20x20 square minus a 10x10 corner notch (area=300, perimeter=100).
    Slab centroid (of the notched polygon) and wall centroid (of its boundary,
    perimeter-weighted) are independently verifiable by hand; CM must be their
    area/(perimeter*3)-weighted average, and CS must equal the wall centroid."""
    from footprint_attributes.geometry import centre_of_mass_and_stiffness

    cm, cs = centre_of_mass_and_stiffness(l_shape_notch_half.geometry)
    slab_centroid = l_shape_notch_half.geometry.iloc[0].centroid
    wall_centroid = l_shape_notch_half.geometry.iloc[0].boundary.centroid

    assert cs[0] == pytest.approx([wall_centroid.x, wall_centroid.y], abs=1e-9)

    area = l_shape_notch_half.geometry.iloc[0].area
    perimeter = l_shape_notch_half.geometry.iloc[0].boundary.length
    wall_weight = perimeter * 3.0
    expected_cm = (
        (area * slab_centroid.x + wall_weight * wall_centroid.x) / (area + wall_weight),
        (area * slab_centroid.y + wall_weight * wall_centroid.y) / (area + wall_weight),
    )
    assert cm[0] == pytest.approx(expected_cm, abs=1e-6)


def test_ec8_radius_ratio_positive(rect_20x10):
    out = shape.EC8(rect_20x10)
    assert out["EC8_radiusRatio"].iloc[0] > 0.0


def test_ec8_compactness_full_for_rectangle(rect_20x10):
    out = shape.EC8(rect_20x10)
    assert out["EC8_compactness"].iloc[0] == pytest.approx(1.0, abs=1e-6)


def test_ec8_compactness_reduced_for_l_shape(l_shape_notch_half):
    """Per the paper: compactness = 1 - (largest hull-minus-footprint piece's
    area) / footprint area. For our L-shape (20x20 square minus a 10x10
    corner notch, area=300), the convex hull is a pentagon of area 350 (it
    already excludes the 50 m^2 triangle beyond the (10,20)-(20,10) diagonal,
    since that triangle is outside the footprint AND outside its own convex
    hull). hull.difference(footprint) is therefore a single 50 m^2 triangle
    -- not the full 100 m^2 notch -- so compactness = 1 - 50/300 = 0.8333."""
    out = shape.EC8(l_shape_notch_half)
    assert out["EC8_compactness"].iloc[0] == pytest.approx(1 - 50 / 300, abs=1e-6)


def test_largest_convex_hull_gap_area_matches_hand_computation(
    l_shape_notch_half, square_10
):
    from footprint_attributes.geometry import largest_convex_hull_gap_area

    assert largest_convex_hull_gap_area(l_shape_notch_half.geometry)[
        0
    ] == pytest.approx(50.0, abs=1e-6)
    assert largest_convex_hull_gap_area(square_10.geometry)[0] == pytest.approx(
        0.0, abs=1e-9
    )


# ─────────────────────────────────────────────────────────────────────────────
# Batch shape() call — shared-work caching must not change results
# ─────────────────────────────────────────────────────────────────────────────


def test_batch_call_matches_individual_calls(rect_20x10):
    batch = shape(rect_20x10, ["EC8_compactness", "slenderness_bbox", "polsby_popper"])
    assert batch["EC8_compactness"].iloc[0] == pytest.approx(
        shape.EC8(rect_20x10)["EC8_compactness"].iloc[0]
    )
    assert batch["slenderness_bbox"].iloc[0] == pytest.approx(
        shape.slenderness.bbox(rect_20x10)[0]
    )
    assert batch["polsby_popper"].iloc[0] == pytest.approx(
        shape.polsby_popper(rect_20x10)[0]
    )


# ─────────────────────────────────────────────────────────────────────────────
# ASCE7 / NTC23 / GNDTII setback ratios — fixed dual-configuration formula
# ─────────────────────────────────────────────────────────────────────────────
#
# Regression tests: ASCE7SetbackRatio, NTC23SetbackRatio, GNDTIIBeta2SetbackRatio
# and GNDTIIBeta6SetbackSlenderness used to divide the *perpendicular
# protrusion* `c` by L1 (or L2) directly, instead of the documented
# min(b1/L1, b2/L2) ratio. On a symmetric shape (L1==L2, b1==b2==c) this
# happens to give the same number, which is exactly why the earlier tests
# (all on a symmetric L-shape) never caught it. `asymmetric_l_shape` (L1=30,
# L2=20, notch b1=20/b2=8) makes the two formulas diverge (0.4 vs 20/30).


def test_asce7_setback_ratio_asymmetric_shape(asymmetric_l_shape):
    out = shape.ASCE7(asymmetric_l_shape)
    assert out["ASCE7_setbackRatio"].iloc[0] == pytest.approx(0.4, abs=1e-6)


def test_ntc23_setback_ratio_matches_asce7_formula(asymmetric_l_shape):
    """NTC-23 uses the identical formula to ASCE7, just a looser limit."""
    asce7 = shape.ASCE7(asymmetric_l_shape)["ASCE7_setbackRatio"].iloc[0]
    ntc23 = shape.NTC23(asymmetric_l_shape)["NTC23_setbackRatio"].iloc[0]
    assert ntc23 == pytest.approx(asce7)


def test_gndtii_beta2_and_beta6_asymmetric_shape(asymmetric_l_shape):
    out = shape.GNDTII(asymmetric_l_shape)
    assert out["GNDTII_beta2_setbackRatio"].iloc[0] == pytest.approx(0.4, abs=1e-6)
    # winning b=8, c=10 (see test_geometry.py hand computation) -> beta6=c/b=1.25
    assert out["GNDTII_beta6_setbackSlenderness"].iloc[0] == pytest.approx(
        1.25, abs=1e-6
    )


def test_setback_ratio_zero_for_convex_shapes(rect_20x10, square_10):
    for gdf in (rect_20x10, square_10):
        assert shape.ASCE7(gdf)["ASCE7_setbackRatio"].iloc[0] == pytest.approx(
            0.0, abs=1e-9
        )
        assert shape.NTC23(gdf)["NTC23_setbackRatio"].iloc[0] == pytest.approx(
            0.0, abs=1e-9
        )
        assert shape.GNDTII(gdf)["GNDTII_beta2_setbackRatio"].iloc[0] == pytest.approx(
            0.0, abs=1e-9
        )


# ─────────────────────────────────────────────────────────────────────────────
# GNDTII beta1 / beta4 — real inscribed-circle 'a', not the bbox L2/L1 proxy
# ─────────────────────────────────────────────────────────────────────────────


def test_gndtii_beta1_rectangle_uses_true_a(rect_20x10):
    """For a rectangle, the inscribed-circle 'a' degenerates to the true
    short side (10), so beta1 = a/L = 10/20 = 0.5 exactly, same as the old
    L2/L1 proxy -- this shape can't distinguish the two implementations,
    it's a basic sanity check."""
    out = shape.GNDTII(rect_20x10)
    assert out["GNDTII_beta1_mainShapeSlenderness"].iloc[0] == pytest.approx(
        0.5, abs=1e-3
    )


def test_gndtii_beta4_uses_a_not_L2(l_shape_notch_half):
    """beta4 = e/a must use the dominant-configuration inscribed-circle 'a'
    (~10 for this L-shape), not simply L2 -- pinned via cross-check against
    the independently computed eccentricity magnitude."""
    from footprint_attributes.geometry import centre_of_mass_and_stiffness

    out = shape.GNDTII(l_shape_notch_half)
    cm, cs = centre_of_mass_and_stiffness(l_shape_notch_half.geometry)
    e_mag = float(((cm - cs) ** 2).sum() ** 0.5)
    beta4 = out["GNDTII_beta4_eccentricityRatio"].iloc[0]
    assert beta4 == pytest.approx(e_mag / 10.0, abs=0.05)


# ─────────────────────────────────────────────────────────────────────────────
# method="bbox" vs method="inertia" — explicit switch, sensible defaults
# ─────────────────────────────────────────────────────────────────────────────


def test_setback_parameters_default_to_bbox_method(rect_20x10_rot30):
    """Per the paper, the GNDT/ASCE7/NTC dual-configuration construction is
    always described as MBB-aligned, so the default method must be "bbox"."""
    default = shape.ASCE7(rect_20x10_rot30)["ASCE7_setbackRatio"].iloc[0]
    explicit_bbox = shape.ASCE7.setbackRatio(rect_20x10_rot30, method="bbox")
    assert default == pytest.approx(explicit_bbox[0])


def test_method_kwarg_is_accepted_and_changes_axes(asymmetric_l_shape):
    """method="inertia" must be a valid override that actually changes which
    axes (and hence which numeric values) are used -- not silently ignored."""
    bbox_ratio = shape.ASCE7.setbackRatio(asymmetric_l_shape, method="bbox")
    inertia_ratio = shape.ASCE7.setbackRatio(asymmetric_l_shape, method="inertia")
    assert bbox_ratio[0] != pytest.approx(inertia_ratio[0])


def test_invalid_method_raises(rect_20x10):
    with pytest.raises(ValueError):
        shape.ASCE7.setbackRatio(rect_20x10, method="not_a_method")


# ─────────────────────────────────────────────────────────────────────────────
# GNDTII cross-parameter caching — must not change results
# ─────────────────────────────────────────────────────────────────────────────
#
# shape.GNDTII(gdf) used to run the expensive inscribed-circle construction
# twice (once each for beta1 and beta4) and the setback construction twice
# (once each for beta2 and beta6). NormAggregate now caches both per method
# and shares them across the matching params within one GNDTII(gdf) call --
# this must be numerically identical to computing each parameter alone.


def test_gndtii_batch_matches_individual_params(asymmetric_l_shape):
    batch = shape.GNDTII(asymmetric_l_shape)
    individual = {
        "GNDTII_beta1_mainShapeSlenderness": shape.GNDTII.beta1_mainShapeSlenderness(
            asymmetric_l_shape
        ),
        "GNDTII_beta2_setbackRatio": shape.GNDTII.beta2_setbackRatio(
            asymmetric_l_shape
        ),
        "GNDTII_beta4_eccentricityRatio": shape.GNDTII.beta4_eccentricityRatio(
            asymmetric_l_shape
        ),
        "GNDTII_beta6_setbackSlenderness": shape.GNDTII.beta6_setbackSlenderness(
            asymmetric_l_shape
        ),
    }
    for col, vals in individual.items():
        assert batch[col].iloc[0] == pytest.approx(vals[0])


def test_gndtii_batch_respects_method_override(asymmetric_l_shape):
    """A caller-supplied method= override must reach every GNDTII parameter
    uniformly, including through the shared cache."""
    bbox_batch = shape.GNDTII(asymmetric_l_shape, method="bbox")
    inertia_batch = shape.GNDTII(asymmetric_l_shape, method="inertia")
    assert bbox_batch["GNDTII_beta2_setbackRatio"].iloc[0] != pytest.approx(
        inertia_batch["GNDTII_beta2_setbackRatio"].iloc[0]
    )
