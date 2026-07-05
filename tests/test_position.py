"""Tests for footprint_attributes.position using synthetic contact scenarios."""

from __future__ import annotations

import pytest

from footprint_attributes import position


def test_isolated_building_classified_isolated(isolated_building):
    out = position(isolated_building)
    assert out["relativePosition"].tolist() == ["isolated"]
    assert out["contact_force"].iloc[0] == pytest.approx(0.0)


def test_no_contacts_in_whole_batch_does_not_crash(two_isolated_buildings):
    """Regression test: contact_forces_df used to raise a KeyError whenever
    *no* building in the batch touched any other (explode_edges returns an
    empty frame, and geoms['edge_center'] = normal_results[0] fails because
    an empty apply(result_type='expand') has no column 0)."""
    out = position(two_isolated_buildings)
    assert out["relativePosition"].tolist() == ["isolated", "isolated"]


def test_lateral_pair_classified_lateral(lateral_pair):
    out = position(lateral_pair)
    assert out["relativePosition"].tolist() == ["lateral", "lateral"]
    assert (out["contact_force"] > 0).all()


def test_corner_building_classified_corner(corner_triplet):
    out = position(corner_triplet)
    # index 0 = centre building, touched on two perpendicular sides
    assert out["relativePosition"].iloc[0] == "corner"
    assert out["relativePosition"].iloc[1] == "lateral"
    assert out["relativePosition"].iloc[2] == "lateral"


def test_confined_building_classified_confined(confined_quartet):
    out = position(confined_quartet)
    # index 0 = centre building, touched on 3 of 4 sides (N, S, E)
    assert out["relativePosition"].iloc[0] == "confined"
    assert all(p == "lateral" for p in out["relativePosition"].iloc[1:])


def test_relative_position_reuses_existing_columns(lateral_pair):
    """position.relative_position() must reuse prefixed force columns rather
    than recompute them (a fast path documented in the module)."""
    full = position(lateral_pair)
    labels = position.relative_position(full)
    assert labels == full["relativePosition"].tolist()


def test_relative_position_computes_when_columns_absent(lateral_pair):
    labels = position.relative_position(lateral_pair)
    assert labels == ["lateral", "lateral"]


def test_buffer_zero_vs_positive_buffer_isolated_stays_isolated(isolated_building):
    out_default = position(isolated_building, buffer=0.0)
    out_buffered = position(isolated_building, buffer=0.15)
    assert out_default["relativePosition"].tolist() == ["isolated"]
    assert out_buffered["relativePosition"].tolist() == ["isolated"]


def test_torque_building_classified_torque(torque_triplet):
    """A confined building with high angular acceleration (opposing forces
    positioned to maximise net torque, per the exact worst-case scenario
    POSITION_DEFAULTS['minAngularAcc'] is derived from) must be upgraded
    from 'confined' to 'torque'."""
    out = position(torque_triplet)
    assert out["relativePosition"].iloc[0] == "torque"
    assert out["contact_angularAcc"].iloc[0] > 2.133


def test_edge_normal_handles_looped_segment_without_nan():
    """Regression test: edge_normal used to compute the tangent as
    coords[-1] - coords[0], which is the zero vector (-> NaN after
    normalising) whenever a multi-point contact segment loops back to its
    starting point despite having non-zero length -- a real occurrence in
    messy digitised footprint boundaries."""
    import numpy as np
    from shapely.geometry import LineString
    from footprint_attributes.geometry import edge_normal

    looped = LineString([(0, 0), (1, 0), (1, 1), (0, 0)])
    _, force = edge_normal(looped, scale=1.0)
    assert not np.isnan(force).any()
    assert np.linalg.norm(force) > 0.0


def test_all_five_position_classes_are_reachable(
    isolated_building, lateral_pair, corner_triplet, confined_quartet, torque_triplet
):
    """Sanity check that every documented relativePosition category
    (isolated, lateral, corner, confined, torque) is actually reachable with
    a concrete, idealised geometric configuration."""
    seen = set()
    for gdf in (
        isolated_building,
        lateral_pair,
        corner_triplet,
        confined_quartet,
        torque_triplet,
    ):
        seen.update(position(gdf)["relativePosition"].tolist())
    assert seen == {"isolated", "lateral", "corner", "confined", "torque"}
