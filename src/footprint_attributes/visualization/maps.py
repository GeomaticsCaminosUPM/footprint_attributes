"""Build the interactive 3D pilot-region map (MapLibre + deck.gl, no
backend) -- one dataset per pilot region, colored by relative position,
shape index, or any of the 15 raw shape/plan-irregularity metrics, with
optional convex-hull / bounding-box / inertia-axis / basic-length /
contact-force-arrow overlays.

Everything here only *writes files* (GeoJSON data + a static HTML/JS/CSS
page rendered from the ``_assets/map`` Jinja2 templates); the map itself
runs entirely in the browser. Needs the ``vis`` extra
(``pip install "footprint-attributes[vis]"``) for ``jinja2``.

See ``examples/generate_interactive_map.py`` for a full worked example
against this package's own pilot-region datasets.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import geopandas as gpd

from . import overlays as _overlays
from ..position import position
from ..shape import shape

_ASSETS_DIR = Path(__file__).parent / "_assets" / "map"

#: The 15 raw, physical-unit shape/plan-irregularity metrics shown by the
#: map (norm limits mirrored into ``main.js``'s own ``NORM_INFO``).
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

#: Overlay ids the map's "Overlays" checkboxes expect a matching
#: ``overlays/<id>.geojson`` file for -- see :mod:`.overlays`.
OVERLAY_IDS = [
    "convex_hull",
    "bounding_box",
    "inertia_axis",
    "basic_lengths",
    "position_arrows",
]


def _default_label(dataset_id: str) -> str:
    return dataset_id.replace("_", " ").title()


def _prepare_gdf(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    gdf = gdf.copy()
    if "id" in gdf.columns:
        gdf = gdf.rename(columns={"id": "building_uid"})
    elif "building_uid" not in gdf.columns:
        gdf["building_uid"] = gdf.index.astype(str)
    return gdf


def _to_geojson_dict(gdf: gpd.GeoDataFrame) -> dict:
    return json.loads(gdf.to_crs(epsg=4326).to_json())


def _write_geojson(gdf: gpd.GeoDataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(_to_geojson_dict(gdf)))


def _datasets_js(labels: dict[str, str]) -> str:
    entries = ",\n  ".join(
        f'{key}: {{ label: {json.dumps(label)}, dir: "data/{key}" }}'
        for key, label in labels.items()
    )
    return "{\n  " + entries + ",\n}"


def _render_page(output_dir: Path, context: dict) -> None:
    try:
        from jinja2 import Environment, FileSystemLoader
    except ImportError as exc:
        raise ImportError(
            'footprint_attributes.visualization requires the "vis" extra: '
            'pip install "footprint-attributes[vis]"'
        ) from exc

    output_dir.mkdir(parents=True, exist_ok=True)
    env = Environment(loader=FileSystemLoader(_ASSETS_DIR), keep_trailing_newline=True)
    for name in ("index.html.j2", "main.js.j2"):
        rendered = env.get_template(name).render(**context)
        (output_dir / name.removesuffix(".j2")).write_text(rendered)
    shutil.copyfile(_ASSETS_DIR / "style.css", output_dir / "style.css")


def build_map(
    datasets: dict[str, gpd.GeoDataFrame],
    output_dir: str | Path,
    *,
    labels: dict[str, str] | None = None,
    default_dataset: str | None = None,
    title: str = "Pilot regions",
    with_overlays: bool = True,
) -> None:
    """Build the interactive map for one or more footprint datasets.

    Computes every column the map can color by (the geometric relative
    position via :func:`~footprint_attributes.position.position`, and the
    15 shape/plan-irregularity metrics + slenderness via
    :func:`~footprint_attributes.shape.shape`) plus, unless
    *with_overlays* is ``False``, the five geometry overlays (convex hull,
    bounding box, inertia axis, basic lengths, contact-force resultant --
    see :mod:`.overlays`), and writes everything ``output_dir`` needs to
    serve the map standalone:

    - ``output_dir/index.html`` / ``main.js`` / ``style.css`` -- the page.
    - ``output_dir/data/<dataset_id>/buildings.geojson`` -- per-dataset
      building attributes.
    - ``output_dir/data/<dataset_id>/overlays/<overlay_id>.geojson`` --
      per-dataset overlay geometry (if *with_overlays*).

    Serve ``output_dir`` over HTTP (``python -m http.server``) -- the page
    loads its data via ``fetch()``, which ``file://`` does not allow.

    Args:
        datasets: Mapping of dataset id -> footprint GeoDataFrame. Each
            must have a ``height`` column (metres) for the 3D extrusion
            and, ideally, an ``id`` column (used as the building's stable
            identifier; falls back to the row index).
        output_dir: Directory to write the map into (created if missing).
        labels: Optional ``{dataset_id: display label}`` override (default:
            title-cased dataset id, e.g. ``"san_jose"`` -> ``"San Jose"``).
        default_dataset: Dataset id shown on first load (default: the
            first key of *datasets*).
        title: Browser tab title.
        with_overlays: Whether to compute and write the five geometry
            overlays. Set ``False`` to skip the extra computation when you
            only need the colored buildings.
    """
    output_dir = Path(output_dir)
    labels = {k: (labels or {}).get(k, _default_label(k)) for k in datasets}
    default_dataset = default_dataset or next(iter(datasets))

    for dataset_id, raw_gdf in datasets.items():
        gdf = _prepare_gdf(raw_gdf)

        result = shape(gdf, columns=[*SHAPE_COLUMNS, *SLENDERNESS_COLUMNS])
        result["building_uid"] = gdf["building_uid"].values
        if "height" in gdf.columns:
            result["height"] = gdf["height"].values
        result["relativePosition"] = position(gdf)["relativePosition"].values
        result = result.set_geometry(gdf.geometry.values, crs=gdf.crs)

        _write_geojson(result, output_dir / "data" / dataset_id / "buildings.geojson")

        if with_overlays:
            height_col = "height" if "height" in gdf.columns else None
            for overlay_id, overlay_gdf in _overlays.build_overlays(
                gdf, height_column=height_col
            ).items():
                _write_geojson(
                    overlay_gdf,
                    output_dir
                    / "data"
                    / dataset_id
                    / "overlays"
                    / f"{overlay_id}.geojson",
                )

    _render_page(
        output_dir,
        {
            "datasets_js": _datasets_js(labels),
            "default_dataset": default_dataset,
            "title": title,
            "subtitle": labels[default_dataset],
        },
    )
