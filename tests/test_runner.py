"""Tests for footprint_attributes.run() -- the single entry-point orchestrator."""

from __future__ import annotations

import pytest

import footprint_attributes
from footprint_attributes import position as position_module
from footprint_attributes import shape as shape_module


def test_run_norm_shorthand_expands_to_all_columns(l_shape_notch_half):
    result = footprint_attributes.run(l_shape_notch_half, config={"columns": ["EC8"]})
    for col in ("EC8_eccentricityRatio", "EC8_radiusRatio", "EC8_compactness"):
        assert col in result.columns
        assert f"compliance_{col}" in result.columns


def test_run_exact_column_name(l_shape_notch_half):
    result = footprint_attributes.run(
        l_shape_notch_half, config={"columns": ["polsby_popper"]}
    )
    assert result["polsby_popper"].iloc[0] == pytest.approx(
        shape_module.polsby_popper(l_shape_notch_half)[0]
    )
    # only the requested column (+ geometry) should be added, not the whole norm set
    assert "EC8_eccentricityRatio" not in result.columns


def test_run_slenderness_shorthand(rect_20x10):
    result = footprint_attributes.run(rect_20x10, config={"columns": ["slenderness"]})
    assert "slenderness_bbox" in result.columns
    assert "slenderness_inertia" in result.columns
    assert result["slenderness_bbox"].iloc[0] == pytest.approx(2.0, abs=1e-6)


def test_run_position_shorthand(lateral_pair):
    result = footprint_attributes.run(lateral_pair, config={"columns": ["position"]})
    for col in (
        "contact_force",
        "contact_confinementRatio",
        "contact_angularAcc",
        "contact_angle",
        "relativePosition",
    ):
        assert col in result.columns
    assert result["relativePosition"].tolist() == ["lateral", "lateral"]


def test_run_position_exact_column_name_triggers_full_pipeline(lateral_pair):
    """Requesting just 'relativePosition' (not the "position" shorthand) must
    still run the full (atomic) position pipeline, since the classification
    can't be computed without the contact-force columns."""
    result = footprint_attributes.run(
        lateral_pair, config={"columns": ["relativePosition"]}
    )
    assert result["relativePosition"].tolist() == ["lateral", "lateral"]


def test_run_bearing_shorthand(rect_20x10):
    result = footprint_attributes.run(rect_20x10, config={"columns": ["bearing"]})
    assert "bearing" in result.columns
    assert result["bearing"].iloc[0] == pytest.approx(0.0, abs=1e-6)


def test_run_returns_original_gdf_with_columns_added(l_shape_notch_half):
    result = footprint_attributes.run(
        l_shape_notch_half, config={"columns": ["polsby_popper"]}
    )
    assert len(result) == len(l_shape_notch_half)
    assert result.geometry.iloc[0].equals(l_shape_notch_half.geometry.iloc[0])


def test_run_no_columns_is_a_no_op(l_shape_notch_half):
    result = footprint_attributes.run(l_shape_notch_half, config={"columns": []})
    assert list(result.columns) == list(l_shape_notch_half.columns)
    result_none_config = footprint_attributes.run(l_shape_notch_half, config=None)
    assert list(result_none_config.columns) == list(l_shape_notch_half.columns)


# ─────────────────────────────────────────────────────────────────────────────
# overwrite semantics
# ─────────────────────────────────────────────────────────────────────────────


def test_run_does_not_recompute_existing_column_by_default(l_shape_notch_half):
    gdf = l_shape_notch_half.copy()
    gdf["polsby_popper"] = 999.0  # obviously-fake pre-existing value
    result = footprint_attributes.run(gdf, config={"columns": ["polsby_popper"]})
    assert result["polsby_popper"].iloc[0] == 999.0


def test_run_overwrite_true_recomputes_existing_column(l_shape_notch_half):
    gdf = l_shape_notch_half.copy()
    gdf["polsby_popper"] = 999.0
    result = footprint_attributes.run(
        gdf, config={"columns": ["polsby_popper"]}, overwrite=True
    )
    assert result["polsby_popper"].iloc[0] == pytest.approx(
        shape_module.polsby_popper(l_shape_notch_half)[0]
    )


def test_run_position_skipped_when_all_columns_present_and_not_overwriting(
    lateral_pair,
):
    full = position_module(lateral_pair, buffer=0.1)
    result = footprint_attributes.run(
        full, config={"columns": ["position"], "position": {"buffer": 0.9}}
    )
    # buffer=0.9 would change the force values if recomputed; unchanged means it was skipped
    assert (result["contact_force"] == full["contact_force"]).all()


def test_run_position_overwrite_true_uses_new_config(lateral_pair):
    cfg = {"columns": ["position"], "position": {"buffer": 0.1}}
    first = footprint_attributes.run(lateral_pair, config=cfg)
    second = footprint_attributes.run(first, config=cfg, overwrite=True)
    assert (second["contact_force"] == first["contact_force"]).all()


# ─────────────────────────────────────────────────────────────────────────────
# method is never configurable through run()
# ─────────────────────────────────────────────────────────────────────────────


def test_run_ignores_method_override_in_shape_config(rect_20x10_rot30):
    """config["shape"]["method"] must be silently stripped -- run() always
    uses each shape parameter's own package default (bbox for setbacks)."""
    default = shape_module.ASCE7.setbackRatio(rect_20x10_rot30)  # method="bbox" default
    via_run = footprint_attributes.run(
        rect_20x10_rot30,
        config={"columns": ["ASCE7_setbackRatio"], "shape": {"method": "inertia"}},
    )
    assert via_run["ASCE7_setbackRatio"].iloc[0] == pytest.approx(default[0])


# ─────────────────────────────────────────────────────────────────────────────
# forwarding kwargs to position()/shape()
# ─────────────────────────────────────────────────────────────────────────────


def test_run_forwards_position_kwargs(torque_triplet):
    """A stricter minAngularAcc threshold should prevent the torque
    upgrade, proving config["position"] kwargs actually reach position()."""
    default = footprint_attributes.run(torque_triplet, config={"columns": ["position"]})
    assert default["relativePosition"].iloc[0] == "torque"

    stricter = footprint_attributes.run(
        torque_triplet,
        config={"columns": ["position"], "position": {"minAngularAcc": 100.0}},
    )
    assert stricter["relativePosition"].iloc[0] != "torque"
