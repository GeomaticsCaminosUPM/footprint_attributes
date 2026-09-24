"""Build the interactive 3D pilot-region map -- one dataset per pilot
region, colored by block position, shape index, or any of the 15 raw
shape/plan-irregularity metrics, with optional convex-hull / bounding-box /
inertia-axis / basic-length / contact-force-arrow overlays.

This module only supplies footprint-specific data (which columns, their
norm limits/colors, the overlay geometry) and delegates the actual 3D
MapLibre + deck.gl page to :mod:`FancyFolium.deck3d`, a reusable engine
that knows nothing about seismic codes or footprints. Needs the ``visualization``
extra (``pip install "footprint-attributes[visualization]"``).

See ``examples/generate_interactive_map.py`` for a full worked example
against this package's own pilot-region datasets.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np

from . import overlays as _overlays
from ..position import position
from ..shape import shape

#: The 15 raw, physical-unit shape/plan-irregularity metrics shown by the
#: map, plus each one's code limit/direction/unit -- mirrors
#: footprint_attributes/config.py's own limits (kept separate since the
#: map's coloring needs a flat, JS-friendly table, not config.py's
#: multi-limit-list structure used for compliance scoring).
_NORM_SPECS = [
    (
        "EC8_eccentricityRatio",
        "EC8 eccentricity ratio",
        0.30,
        True,
        "",
        "EC8: eccentricity <= 0.30",
        "EC8",
    ),
    (
        "EC8_radiusRatio",
        "EC8 radius ratio",
        1.0,
        False,
        "",
        "EC8: radius ratio >= 1.00",
        "EC8",
    ),
    (
        "EC8_compactness",
        "EC8 compactness",
        0.95,
        False,
        "",
        "EC8: compactness >= 0.95",
        "EC8",
    ),
    (
        "ASCE7_setbackRatio",
        "ASCE 7 setback ratio",
        0.2,
        True,
        "",
        "ASCE 7: setback ratio <= 0.20",
        "ASCE7",
    ),
    (
        "ASCE7_holeRatio",
        "ASCE 7 hole ratio",
        0.25,
        True,
        "",
        "ASCE 7: hole ratio <= 0.25",
        "ASCE7",
    ),
    (
        "ASCE7_parallelityAngle",
        "ASCE 7 parallelity angle",
        10,
        True,
        "°",
        "ASCE 7: angle <= 5°",
        "ASCE7",
    ),
    (
        "GNDTII_beta1_mainShapeSlenderness",
        "GNDT-II β1 (slenderness)",
        0.4,
        False,
        "",
        "GNDT-II β1 >= 0.8",
        "GNDTII",
    ),
    (
        "GNDTII_beta2_setbackRatio",
        "GNDT-II β2 (setback ratio)",
        0.3,
        True,
        "",
        "GNDT-II β2 <= 0.1",
        "GNDTII",
    ),
    (
        "GNDTII_beta4_eccentricityRatio",
        "GNDT-II β4 (eccentricity ratio)",
        0.4,
        True,
        "",
        "GNDT-II β4 <= 0.2",
        "GNDTII",
    ),
    (
        "GNDTII_beta6_setbackSlenderness",
        "GNDT-II β6 (setback slenderness)",
        0.25,
        False,
        "",
        "GNDT-II β6 >= 0.5",
        "GNDTII",
    ),
    (
        "CSCR2010_eccentricityRatio",
        "CSCR-2010 eccentricity ratio",
        0.25,
        True,
        "",
        "CSCR-2010: eccentricity <= 0.05",
        "CSCR2010",
    ),
    (
        "NTC23_setbackRatio",
        "NTC-23 setback ratio",
        0.4,
        True,
        "",
        "NTC-23: setback ratio <= 0.40",
        "NTC23",
    ),
    (
        "NTC23_holeRatio",
        "NTC-23 hole ratio",
        0.4,
        True,
        "",
        "NTC-23: hole ratio <= 0.40",
        "NTC23",
    ),
    (
        "slenderness_bbox",
        "Slenderness (bbox)",
        4.0,
        True,
        "",
        "EC8: slenderness (bbox) <= 4.0",
        "slenderness",
    ),
    (
        "slenderness_inertia",
        "Slenderness (inertia)",
        4.0,
        True,
        "",
        "EC8: slenderness (inertia) <= 4.0",
        "slenderness",
    ),
]
SHAPE_COLUMNS = [name for name, *_ in _NORM_SPECS if not name.startswith("slenderness")]
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

_POSITION_COLORS = {
    "isolated": "#4299e1",
    "lateral": "#e2a33f",
    "corner": "#4fbf8f",
    "confined": "#8744ad",
    "torque": "#c0392b",
    "unlabeled": "#6b7280",
}
_POSITION_ORDER = ["isolated", "lateral", "corner", "confined", "torque", "unlabeled"]
_POSITION_LABELS = {
    "isolated": "Isolated",
    "lateral": "Lateral",
    "corner": "Corner",
    "confined": "Confined",
    "torque": "Torque",
    "unlabeled": "Unlabeled",
}

_FSI_COLORS = {
    "regular": "#4fbf8f",
    "shape": "#e2793f",
    "eccentricity": "#e2c23f",
    "slenderness": "#c0392b",
}
_FSI_ORDER = ["regular", "shape", "eccentricity", "slenderness"]
_FSI_LABELS = {
    "regular": "Regular",
    "shape": "Irregular: shape",
    "eccentricity": "Irregular: eccentricity",
    "slenderness": "Irregular: slenderness",
}

_OVERLAY_STYLE = {
    "convex_hull": ("Convex hull", "#4299e1", "polygon"),
    "bounding_box": ("Bounding box", "#e2a33f", "polygon"),
    "inertia_axis": ("Inertia axis", "#e2938a", "line"),
    "basic_lengths": ("Building lengths (L1/L2/a1/a2/b/c)", "#111111", "line"),
    "position_arrows": ("Position force arrows", "#8744ad", "line"),
}


def _default_label(dataset_id: str) -> str:
    return dataset_id.replace("_", " ").title()


def _prepare_gdf(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    gdf = gdf.copy()
    if "id" not in gdf.columns:
        gdf["id"] = gdf.index.astype(str)
    return gdf


def _shape_index(result: gpd.GeoDataFrame) -> np.ndarray:
    """The 4-category shape index shown on the map: start "regular",
    eccentricity ratio > 0.3 -> "eccentricity", setback ratio > 0.2 ->
    "shape" (overrides eccentricity), slenderness > 4.0 -> "slenderness"
    (overrides everything) -- so slenderness always wins if triggered
    (worst), eccentricity is overridden by either of the other two
    (mildest). Matches the paper's own analysis notebooks' exact
    thresholds/override order.
    """
    key = np.full(len(result), "regular", dtype=object)
    key[result["EC8_eccentricityRatio"] > 0.3] = "eccentricity"
    key[result["ASCE7_setbackRatio"] > 0.2] = "shape"
    key[result["slenderness_inertia"] > 4.0] = "slenderness"
    return key


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
    position via :func:`~footprint_attributes.position.position`, the 15
    shape/plan-irregularity metrics + slenderness via
    :func:`~footprint_attributes.shape.shape`, and the shape index derived
    from 3 of those) plus, unless *with_overlays* is ``False``, the five
    geometry overlays (convex hull, bounding box, inertia axis, basic
    lengths, contact-force resultant -- see :mod:`.overlays`), then calls
    :func:`FancyFolium.deck3d.build_interactive_map` to render the page.

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
    try:
        from FancyFolium.deck3d import (
            CategoricalAttribute,
            NormAttribute,
            OverlayDef,
            build_interactive_map,
        )
    except ImportError as exc:
        raise ImportError(
            'footprint_attributes.visualization requires the "vis" extra: '
            'pip install "footprint-attributes[visualization]"'
        ) from exc

    prepared = {name: _prepare_gdf(gdf) for name, gdf in datasets.items()}
    map_data: dict[str, gpd.GeoDataFrame] = {}
    overlay_data: dict[str, dict[str, gpd.GeoDataFrame]] = {}

    for dataset_id, gdf in prepared.items():
        result = shape(gdf, columns=[*SHAPE_COLUMNS, *SLENDERNESS_COLUMNS])
        result["id"] = gdf["id"].values
        if "height" in gdf.columns:
            result["height"] = gdf["height"].values
        result["blockPosition"] = position(gdf)["blockPosition"].values
        result["shape_index"] = _shape_index(result)
        map_data[dataset_id] = result.set_geometry(gdf.geometry.values, crs=gdf.crs)

        if with_overlays:
            height_col = "height" if "height" in gdf.columns else None
            overlay_data[dataset_id] = _overlays.build_overlays(
                gdf, height_column=height_col
            )

    attributes = [
        CategoricalAttribute(
            "blockPosition",
            "Block position",
            _POSITION_COLORS,
            _POSITION_ORDER,
            _POSITION_LABELS,
        ),
        CategoricalAttribute(
            "shape_index", "Shape index", _FSI_COLORS, _FSI_ORDER, _FSI_LABELS
        ),
        *(
            NormAttribute(
                name,
                label,
                limit,
                worse_is_high=worse_is_high,
                unit=unit,
                criteria=criteria,
                family=family,
            )
            for name, label, limit, worse_is_high, unit, criteria, family in _NORM_SPECS
        ),
    ]
    overlay_defs = (
        [OverlayDef(oid, *_OVERLAY_STYLE[oid]) for oid in OVERLAY_IDS]
        if with_overlays
        else []
    )

    build_interactive_map(
        map_data,
        attributes,
        output_dir,
        overlays=overlay_data,
        overlay_defs=overlay_defs,
        labels=labels,
        default_dataset=default_dataset,
        title=title,
    )
