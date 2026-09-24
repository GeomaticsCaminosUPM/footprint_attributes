"""Geometry-overlay GeoDataFrames for :func:`.maps.build_map`'s interactive
map: convex hull, minimum-rotated bounding box, inertia axis, the
L1/L2/a1/a2/b/c "basic lengths" dimension set, and the net contact-force
resultant arrow per building
(:func:`~footprint_attributes.position.contact_force_vectors`).

Each overlay's geometry is built at "arrow_gdf" scale (i.e. real physical
units in the footprints' own projected CRS), matching what the example
notebooks already draw with ``FancyFolium``/matplotlib -- this module only
adapts that same geometry into plain GeoDataFrames for the deck.gl map.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd

from .. import direction
from ..geometry import ensure_projected, to_gdf
from ..notebook_utils import ARROW_STYLE, basic_length_axis
from ..position import contact_force_vectors
from ..viz import arrow_gdf


def _convex_hull(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    out = gdf[["id"]].copy()
    out["geometry"] = gdf.geometry.convex_hull
    return gpd.GeoDataFrame(out, geometry="geometry", crs=gdf.crs)


def _bounding_box(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    out = gdf[["id"]].copy()
    out["geometry"] = gdf.geometry.apply(lambda g: g.minimum_rotated_rectangle)
    return gpd.GeoDataFrame(out, geometry="geometry", crs=gdf.crs)


def _inertia_axis(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    # Reimplements notebook_utils.direction_arrow_gdfs's two specs directly
    # (rather than calling it) so `id` can ride along as an
    # arrow_gdf extra_column, correctly sliced to whichever rows actually
    # get an arrow -- see the _basic_lengths comment below for why that
    # matters.
    L1, dir1, L2, dir2, _ = direction.inertia(gdf, mode="all")
    centroids = np.column_stack([gdf.geometry.centroid.x, gdf.geometry.centroid.y])
    uid = gdf["id"].values
    L1, L2 = np.asarray(L1, dtype=float), np.asarray(L2, dtype=float)
    dir1, dir2 = np.asarray(dir1, dtype=float), np.asarray(dir2, dtype=float)
    arrows = [
        arrow_gdf(
            centroids,
            dir1 * L1[:, None],
            gdf.crs,
            caps="arrow",
            kind="L1",
            value=L1,
            id=uid,
        ),
        arrow_gdf(
            centroids,
            dir2 * L2[:, None],
            gdf.crs,
            caps="arrow",
            kind="L2",
            value=L2,
            id=uid,
        ),
    ]
    return gpd.GeoDataFrame(pd.concat(arrows, ignore_index=True), crs=gdf.crs)


def _bbox_axis(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    # Same as _inertia_axis, but the L1/L2 axes come from the minimum
    # rotated bounding box (direction.bbox) instead of the principal
    # moments of inertia -- the other of the two "direction" conventions,
    # for side-by-side comparison against the inertia-based axis.
    L1, dir1, L2, dir2, _ = direction.bbox(gdf, mode="all")
    centroids = np.column_stack([gdf.geometry.centroid.x, gdf.geometry.centroid.y])
    uid = gdf["id"].values
    L1, L2 = np.asarray(L1, dtype=float), np.asarray(L2, dtype=float)
    dir1, dir2 = np.asarray(dir1, dtype=float), np.asarray(dir2, dtype=float)
    arrows = [
        arrow_gdf(
            centroids,
            dir1 * L1[:, None],
            gdf.crs,
            caps="arrow",
            kind="L1",
            value=L1,
            id=uid,
        ),
        arrow_gdf(
            centroids,
            dir2 * L2[:, None],
            gdf.crs,
            caps="arrow",
            kind="L2",
            value=L2,
            id=uid,
        ),
    ]
    return gpd.GeoDataFrame(pd.concat(arrows, ignore_index=True), crs=gdf.crs)


def _basic_lengths(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    # Reimplements notebook_utils.dimension_arrow_gdfs's specs directly
    # (rather than calling it) so `id` can ride along as an
    # arrow_gdf extra_column: a1/a2/b/c aren't defined for every building
    # (e.g. no setback -> no b/c), and arrow_gdf drops those rows via its
    # own min_length filter -- assigning `id` *after* the fact,
    # positionally, would then misalign with whichever rows survived.
    L1, dir1, L2, dir2, _ = direction.bbox(gdf, mode="all")
    axis = basic_length_axis(gdf, np.asarray(L1), dir1, np.asarray(L2), dir2)
    uid = gdf["id"].values
    centered_specs = {
        "L1": (axis["centroids"], axis["dir1"] * axis["L1"][:, None]),
        "L2": (axis["centroids"], axis["dir2"] * axis["L2"][:, None]),
        "a1": (axis["a_centers"], axis["dir2"] * axis["a1"][:, None]),
        "a2": (axis["a_centers"], axis["dir1"] * axis["a2"][:, None]),
    }
    arrows = [
        arrow_gdf(
            anchor,
            vec,
            gdf.crs,
            caps=ARROW_STYLE[kind]["caps"],
            head_frac=ARROW_STYLE[kind]["head_frac"],
            kind=kind,
            value=np.linalg.norm(vec, axis=1),
            id=uid,
        )
        for kind, (anchor, vec) in centered_specs.items()
    ]
    for kind in ("b", "c"):
        arrows.append(
            arrow_gdf(
                axis[f"{kind}_start"],
                axis[f"{kind}_vector"],
                gdf.crs,
                centered=False,
                caps=ARROW_STYLE[kind]["caps"],
                head_frac=ARROW_STYLE[kind]["head_frac"],
                kind=kind,
                value=np.linalg.norm(axis[f"{kind}_vector"], axis=1),
                id=uid,
            )
        )
    return gpd.GeoDataFrame(pd.concat(arrows, ignore_index=True), crs=gdf.crs)


def _basic_lengths_inertia(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    # Same as _basic_lengths, but axes come from the principal-inertia
    # directions (direction.inertia) instead of the minimum-rotated bounding
    # box (direction.bbox) -- an alternate "which axes count as L1/L2" basis
    # for the same a1/a2/b/c setback construction, for side-by-side
    # comparison against the bbox-based version.
    L1, dir1, L2, dir2, _ = direction.inertia(gdf, mode="all")
    axis = basic_length_axis(gdf, np.asarray(L1), dir1, np.asarray(L2), dir2)
    uid = gdf["id"].values
    centered_specs = {
        "L1": (axis["centroids"], axis["dir1"] * axis["L1"][:, None]),
        "L2": (axis["centroids"], axis["dir2"] * axis["L2"][:, None]),
        "a1": (axis["a_centers"], axis["dir2"] * axis["a1"][:, None]),
        "a2": (axis["a_centers"], axis["dir1"] * axis["a2"][:, None]),
    }
    arrows = [
        arrow_gdf(
            anchor,
            vec,
            gdf.crs,
            caps=ARROW_STYLE[kind]["caps"],
            head_frac=ARROW_STYLE[kind]["head_frac"],
            kind=kind,
            value=np.linalg.norm(vec, axis=1),
            id=uid,
        )
        for kind, (anchor, vec) in centered_specs.items()
    ]
    for kind in ("b", "c"):
        arrows.append(
            arrow_gdf(
                axis[f"{kind}_start"],
                axis[f"{kind}_vector"],
                gdf.crs,
                centered=False,
                caps=ARROW_STYLE[kind]["caps"],
                head_frac=ARROW_STYLE[kind]["head_frac"],
                kind=kind,
                value=np.linalg.norm(axis[f"{kind}_vector"], axis=1),
                id=uid,
            )
        )
    return gpd.GeoDataFrame(pd.concat(arrows, ignore_index=True), crs=gdf.crs)


def _position_arrows(
    gdf: gpd.GeoDataFrame, height_column: str | None
) -> gpd.GeoDataFrame:
    # Mirrors examples/position.ipynb's own force-arrow drawing: edges and
    # resultants come back scaled by building height (force = height x wall
    # length -- see position.contact_force_vectors), so divide back out by
    # height to recover a vector in real, physical metres -- otherwise a
    # tall building's arrows blow up to many times the footprint's own size.
    # Both layers are included: "contact_force_edge" (one per touching wall
    # segment) and "contact_force_resultant" (net force per building).
    edges, resultants = contact_force_vectors(gdf, height_column=height_column)
    if len(resultants) == 0:
        return gpd.GeoDataFrame({"id": [], "kind": []}, geometry=[], crs=gdf.crs)

    height_by_id = (
        gdf[height_column].astype(float)
        if height_column
        else pd.Series(1.0, index=gdf.index)
    )

    edge_height = height_by_id.loc[edges["geom_id"]].to_numpy()
    edge_anchors = np.stack(edges["anchor"].to_numpy())
    edge_vectors = np.stack(edges["vector"].to_numpy()) / edge_height[:, None]
    edge_uid = gdf.loc[edges["geom_id"], "id"].to_numpy()
    edge_arrows = arrow_gdf(
        edge_anchors,
        edge_vectors,
        gdf.crs,
        centered=False,
        caps="arrow",
        kind="contact_force_edge",
        id=edge_uid,
    )

    res_height = height_by_id.loc[resultants["geom_id"]].to_numpy()
    res_anchors = np.stack(resultants["anchor"].to_numpy())
    res_vectors = np.stack(resultants["vector"].to_numpy()) / res_height[:, None]
    res_uid = gdf.loc[resultants["geom_id"], "id"].to_numpy()
    resultant_arrows = arrow_gdf(
        res_anchors,
        res_vectors,
        gdf.crs,
        centered=False,
        caps="arrow",
        kind="contact_force_resultant",
        id=res_uid,
    )

    return gpd.GeoDataFrame(
        pd.concat([edge_arrows, resultant_arrows], ignore_index=True), crs=gdf.crs
    )


def build_overlays(
    gdf: gpd.GeoDataFrame,
    *,
    height_column: str | None = "height",
) -> dict[str, gpd.GeoDataFrame]:
    """Build all five map overlays for one dataset.

    Args:
        gdf: Footprints with an ``id`` column (see
            :func:`.maps.build_map`), any CRS.
        height_column: Column with building heights in metres, forwarded to
            :func:`~footprint_attributes.position.contact_force_vectors`
            (``None`` treats every building as height 1, i.e. force
            proportional to touching-edge length only).

    Returns:
        ``{"convex_hull": gdf, "bounding_box": gdf, "inertia_axis": gdf,
        "basic_lengths": gdf, "basic_lengths_inertia": gdf,
        "position_arrows": gdf}``, each in *gdf*'s CRS.
    """
    gdf = ensure_projected(to_gdf(gdf))
    return {
        "convex_hull": _convex_hull(gdf),
        "bounding_box": _bounding_box(gdf),
        "inertia_axis": _inertia_axis(gdf),
        "bbox_axis": _bbox_axis(gdf),
        "basic_lengths": _basic_lengths(gdf),
        "basic_lengths_inertia": _basic_lengths_inertia(gdf),
        "position_arrows": _position_arrows(gdf, height_column),
    }
