#!/usr/bin/env python3
"""One-off preprocessing: examples/data/*.gpkg -> docs/_static/interactive/data/.

Writes, per pilot region (guatemala, san_jose, santo_domingo):
  data/<name>/buildings.geojson              -- blockPosition + shape_index + 15 shape metrics
  data/<name>/overlays/position_arrows.geojson    -- contact-force arrows (edges + resultant), split
                                                      into shaft/head parts with filled-triangle tips
  data/<name>/overlays/bbox_axis.geojson          -- L1/L2 direction arrows (bounding-box method)
  data/<name>/overlays/inertia_axis.geojson       -- L1/L2 direction arrows (inertia method)
  data/<name>/overlays/basic_lengths.geojson      -- L1/L2/a1/a2/b/c dimension arrows (bbox method)
  data/<name>/overlays/basic_lengths_inertia.geojson -- same, inertia method

Mirrors SevillaConference/data/prepare_data.py -- see that script's comments
for why the contact-force arrows are shortened+flipped and split into a
dashed shaft / filled-triangle head.

Run with the package's own venv:
  .venv/bin/python docs/prepare_interactive_data.py
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

import geopandas as gpd
import pandas as pd
import shapely
from shapely.affinity import scale as shapely_scale
from shapely.geometry import Polygon
from footprint_attributes import direction, position, shape
from footprint_attributes.geometry import ensure_projected, to_gdf
from footprint_attributes.visualization.overlays import build_overlays

DATA_ROOT = os.path.join(ROOT, "examples", "data")
OUT_ROOT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "_static", "interactive", "data"
)

# dataset id -> (source gpkg, display label)
DATASETS = {
    "guatemala": ("guatemala_pilot_region.gpkg", "Guatemala City (Zona 10)"),
    "san_jose": ("san_jose_pilot_region.gpkg", "San José (Mata Redonda)"),
    "santo_domingo": (
        "santo_domingo_pilot_region.gpkg",
        "Santo Domingo (Ensanche Quisquella)",
    ),
}

SHAPE_COLUMNS = [
    "EC8_eccentricityRatio",
    "EC8_radiusRatio",
    "EC8_compactness",
    "ASCE7_setbackRatio",
    "ASCE7_holeRatio",
    "ASCE7_parallelityAngle",
    "GNDTII_beta1_mainShapeSlenderness",
    "GNDTII_beta2_setbackRatio",
    "GNDTII_beta4_eccentricityRatio",
    "GNDTII_beta6_setbackSlenderness",
    "CSCR2010_eccentricityRatio",
    "NTC23_setbackRatio",
    "NTC23_holeRatio",
]
SLENDERNESS_COLUMNS = ["slenderness_bbox", "slenderness_inertia"]

# Both the resultant and per-edge contact-force vectors can be many times a
# building's own size -- shrink them toward their own start point so they
# read as arrows on the building. Negative factors both shrink AND flip the
# vector 180 degrees (a point reflection through the anchor): raw
# contact_force_vectors point toward the neighbour pressing on a building,
# which reads backwards on a map -- people expect the arrow to point away
# from what's pushing them (see SevillaConference/data/prepare_data.py for
# the full derivation/verification of this).
RESULTANT_LENGTH_SCALE = -0.3
EDGE_LENGTH_SCALE = -0.35


def shorten_geom(geom, scale):
    origin = shapely.get_parts(geom)[0].coords[0]
    return shapely_scale(geom, xfact=scale, yfact=scale, origin=origin)


def split_shaft_head(geom):
    """Split an arrow_gdf caps="arrow" MultiLineString into its plain shaft
    (part 0) and a FILLED triangle Polygon built from the tip + the two
    arrowhead barb endpoints (parts 1-2) -- a solid triangle degrades
    gracefully at low zoom, where two thin diverging strokes blur into a
    formless smudge."""
    parts = shapely.get_parts(geom)
    shaft = parts[0]
    tip = shapely.get_coordinates(parts[0])[-1]
    back1 = shapely.get_coordinates(parts[1])[-1]
    back2 = shapely.get_coordinates(parts[2])[-1]
    head = Polygon([tip, back1, back2])
    return shaft, head


def shaft_head_rows(rows, scale):
    rows = rows.copy()
    rows["geometry"] = rows["geometry"].apply(shorten_geom, scale=scale)
    shaft_and_head = rows["geometry"].apply(split_shaft_head)
    shaft_rows = rows.copy()
    shaft_rows["geometry"] = [s for s, _ in shaft_and_head]
    shaft_rows["part"] = "shaft"
    head_rows = rows.copy()
    head_rows["geometry"] = [h for _, h in shaft_and_head]
    head_rows["part"] = "head"
    return shaft_rows, head_rows


def to_largest_polygon(geom):
    """make_valid() on a self-intersecting footprint can return a
    GeometryCollection (polygon + a degenerate line/point sliver);
    footprint_attributes requires pure Polygon geometry, so keep only the
    largest polygonal part."""
    if geom.geom_type == "Polygon":
        return geom
    parts = [
        g for g in shapely.get_parts(geom) if g.geom_type in ("Polygon", "MultiPolygon")
    ]
    if not parts:
        return geom
    polys = []
    for p in parts:
        polys.extend(list(p.geoms) if p.geom_type == "MultiPolygon" else [p])
    return max(polys, key=lambda p: p.area)


def main():
    for name, (src, label) in DATASETS.items():
        print(f"=== {name} ({label}) ===")
        raw = gpd.read_file(os.path.join(DATA_ROOT, src))
        n_before = len(raw)
        raw = raw[raw.geometry.notna() & ~raw.geometry.is_empty].reset_index(drop=True)
        if len(raw) != n_before:
            print(f"  dropped {n_before - len(raw)} rows with missing/empty geometry")

        gdf = raw[["id", "geometry"]].copy()
        gdf.geometry = gdf.geometry.make_valid()
        gdf = ensure_projected(to_gdf(gdf))
        gdf.geometry = gdf.geometry.make_valid()
        gdf.geometry = gdf.geometry.apply(to_largest_polygon)

        result = shape(gdf, columns=[*SHAPE_COLUMNS, *SLENDERNESS_COLUMNS])
        result["id"] = gdf["id"].values
        result["blockPosition"] = position(gdf)["blockPosition"].values
        # Bearing under both direction conventions, for the "direction" map's
        # bbox-vs-inertia coloring toggle.
        result["bearing_bbox"] = direction.bbox(gdf, mode="bearing")
        result["bearing_inertia"] = direction.inertia(gdf, mode="bearing")

        for c in [*SHAPE_COLUMNS, *SLENDERNESS_COLUMNS]:
            print(f"  {c}: {result[c].notna().sum()}/{len(result)}")
        print(
            "  blockPosition counts:", result["blockPosition"].value_counts().to_dict()
        )

        result = result.to_crs(epsg=4326)

        extra = raw[["id", "height"]].copy()
        merged = result.merge(extra, on="id", how="left")
        merged = merged.rename(columns={"id": "building_uid"})

        out_dir = f"{OUT_ROOT}/{name}"
        os.makedirs(out_dir, exist_ok=True)
        merged.to_file(f"{out_dir}/buildings.geojson", driver="GeoJSON")
        print(f"  wrote {out_dir}/buildings.geojson")

        overlays = build_overlays(
            gdf.merge(extra, on="id", how="left"), height_column="height"
        )

        overlay_dir = f"{out_dir}/overlays"
        os.makedirs(overlay_dir, exist_ok=True)

        # Contact-force arrows: shorten+flip, split into dashed shaft + filled head.
        position_arrows = overlays["position_arrows"].rename(
            columns={"id": "building_uid"}
        )
        is_resultant = position_arrows["kind"] == "contact_force_resultant"
        is_edge = position_arrows["kind"] == "contact_force_edge"

        resultant_shaft, resultant_head = shaft_head_rows(
            position_arrows.loc[is_resultant], RESULTANT_LENGTH_SCALE
        )
        edge_shaft, edge_head = shaft_head_rows(
            position_arrows.loc[is_edge], EDGE_LENGTH_SCALE
        )
        position_arrows = gpd.GeoDataFrame(
            pd.concat(
                [resultant_shaft, resultant_head, edge_shaft, edge_head],
                ignore_index=True,
            ),
            crs=position_arrows.crs,
        )
        position_arrows.to_crs(epsg=4326).to_file(
            f"{overlay_dir}/position_arrows.geojson", driver="GeoJSON"
        )
        print(f"  wrote {overlay_dir}/position_arrows.geojson")

        for overlay_id in ["inertia_axis", "basic_lengths", "basic_lengths_inertia"]:
            layer = overlays[overlay_id].rename(columns={"id": "building_uid"})
            layer.to_crs(epsg=4326).to_file(
                f"{overlay_dir}/{overlay_id}.geojson", driver="GeoJSON"
            )
            print(f"  wrote {overlay_dir}/{overlay_id}.geojson")

        # bbox_axis is split into a dashed shaft + solid filled-triangle head
        # (scale=1, no shrink) for the "direction" map, which draws the bbox
        # axis dotted (except its tip) against the inertia axis's solid line.
        bbox_shaft, bbox_head = shaft_head_rows(
            overlays["bbox_axis"].rename(columns={"id": "building_uid"}), 1.0
        )
        bbox_axis = gpd.GeoDataFrame(
            pd.concat([bbox_shaft, bbox_head], ignore_index=True),
            crs=overlays["bbox_axis"].crs,
        )
        bbox_axis.to_crs(epsg=4326).to_file(
            f"{overlay_dir}/bbox_axis.geojson", driver="GeoJSON"
        )
        print(f"  wrote {overlay_dir}/bbox_axis.geojson")


if __name__ == "__main__":
    main()
