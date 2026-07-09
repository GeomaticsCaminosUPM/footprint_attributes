"""Helpers shared by the example notebooks under ``examples/`` (not used by
the core attribute-computation pipeline).

Every notebook that loads one of the pilot-region files, or draws L1/L2/
a1/a2/b/c as arrows on a plot or an interactive map, does it through the
functions here instead of redefining its own copy -- keeps the notebooks
short and guarantees all of them treat "loading a city" and "drawing a
dimension" the same way.

Drawing convention (deliberately black-only, distinguished by line style,
not by colour -- easy to read on a busy satellite basemap and in print):

- ``L1``/``L2`` (and a building's own axis ``dir1``/``dir2``) are *directed*
  quantities -> drawn as a real arrow (one arrowhead at the tip).
- ``a1``/``a2``/``b``/``c`` are plain *measurements* (no direction of their
  own) -> drawn as a dimension line: a shaft with a small perpendicular
  tick ("|") at each end, no arrowhead.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import geopandas as gpd
import numpy as np
import shapely
from shapely.geometry import LineString, Point

from .geometry import ensure_projected
from .viz import arrow_gdf

if TYPE_CHECKING:
    # Only needed for type hints -- both are optional, notebook-only
    # dependencies (see the `visualization` extra in pyproject.toml) that
    # this module otherwise imports lazily, inside the functions that use
    # them, so importing them here unconditionally would turn them into a
    # hard runtime dependency of the whole package.
    import folium
    import matplotlib.axes

# Per-"kind" drawing style: whether it gets an arrowhead ("arrow") or two
# end ticks ("bar"), the line-dash pattern used to tell kinds apart now that
# every line is black, and how long the end ticks are (`head_frac`, as a
# fraction of the line's own length -- see `viz.arrow_gdf`). `mpl` is a
# matplotlib linestyle (or dash tuple); `dashArray` is the equivalent for
# FancyFolium/Leaflet (SVG dash pattern, in pixels).
#: Each kind's pattern is chosen to be unmistakably different from every
#: other one, not just "technically distinct" -- e.g. a1's sparse dots and
#: c's old pattern used to both render as plain fine dots, so a1 (main-
#: element width, drawn near the footprint's own body) and c (a setback's
#: protrusion depth, drawn near the notch) could look like the same line
#: continuing in the wrong place, even though the underlying numbers were
#: correct all along. a2/b had the same near-collision (both dash-dot-ish).
ARROW_STYLE: dict[str, dict[str, Any]] = {
    "L1": {"caps": "arrow", "mpl": "-", "dashArray": None, "head_frac": 0.2},
    "L2": {"caps": "arrow", "mpl": "--", "dashArray": "10,6", "head_frac": 0.2},
    "a1": {"caps": "bar", "mpl": ":", "dashArray": "1,5", "head_frac": 0.6},
    "a2": {
        "caps": "bar",
        "mpl": (0, (3, 5, 1, 5)),
        "dashArray": "6,4,1,4",
        "head_frac": 0.6,
    },
    "b": {"caps": "bar", "mpl": (0, (5, 1)), "dashArray": "8,2", "head_frac": 0.2},
    "c": {
        "caps": "bar",
        "mpl": (0, (3, 5, 1, 5, 1, 5)),
        "dashArray": "6,2,1,2,1,2",
        "head_frac": 0.2,
    },
}

# Matplotlib-only emphasis: `c` (a setback's protrusion depth) sits right on
# top of the building outline and, drawn black like everything else, is easy
# to mistake for a piece of `L1`/`L2` crossing behind it. Give it its own
# colour and a heavier line just in the static (non-map) plots below --
# ``add_dimension_layers`` (the interactive map) intentionally stays
# black-only per the module-level convention, since it isn't printed and the
# layer control already lets you isolate one kind at a time.
MPL_COLOR: dict[str, str] = {
    "L1": "black",
    "L2": "black",
    "a1": "black",
    "a2": "black",
    "b": "black",
    "c": "crimson",
}
MPL_LINEWIDTH: dict[str, float] = {
    "L1": 1.5,
    "L2": 1.5,
    "a1": 1.5,
    "a2": 1.5,
    "b": 1.5,
    "c": 2.5,
}

# Axis-method line style (bbox vs. inertia), used where the *same* kind is
# drawn twice, once per method, and needs its own way to tell the two apart.
METHOD_STYLE: dict[str, dict[str, Any]] = {
    "bbox": {"mpl": "-", "dashArray": None},
    "inertia": {"mpl": "--", "dashArray": "8,5"},
}


def load_city(path: str, fill_height_column: str | None = None) -> gpd.GeoDataFrame:
    """Load a pilot-region file into a clean, single-Polygon-per-row GeoDataFrame.

    Always explodes multi-part geometries, then drops (with a warning) any
    row that is null/empty or fails a validity check, instead of silently
    repairing it -- ``footprint_attributes`` functions reject invalid/
    multi-part geometries outright, so bad rows must be removed here, not
    patched. Finally reprojects to a projected (metric) CRS via
    :func:`footprint_attributes.geometry.ensure_projected`, since every
    length/area/direction computation in this package requires one.

    Args:
        path: Path to a vector file readable by ``geopandas.read_file``
            (e.g. a ``.gpkg``).
        fill_height_column: If given, missing values in this column are
            filled with the column's own median (e.g. ``"height"`` --
            :func:`footprint_attributes.position` needs a real height for
            every building).

    Returns:
        A clean GeoDataFrame in a projected CRS.
    """
    gdf = gpd.read_file(path)

    null_mask = gdf.geometry.isna() | gdf.geometry.is_empty
    if null_mask.any():
        print(f"WARNING: {path}: dropping {null_mask.sum()} null/empty geometries")
    gdf = gdf[~null_mask].reset_index(drop=True)

    gdf = gdf.explode(index_parts=False).reset_index(drop=True)

    non_polygon = gdf.geometry.type != "Polygon"
    if non_polygon.any():
        print(f"WARNING: {path}: dropping {non_polygon.sum()} non-Polygon geometries")
    gdf = gdf[~non_polygon].reset_index(drop=True)

    invalid = ~gdf.geometry.is_valid
    if invalid.any():
        print(f"WARNING: {path}: dropping {invalid.sum()} invalid geometries")
    gdf = gdf[~invalid].reset_index(drop=True)

    gdf = ensure_projected(gdf)

    if fill_height_column is not None and fill_height_column in gdf.columns:
        gdf[fill_height_column] = gdf[fill_height_column].fillna(
            gdf[fill_height_column].median()
        )

    return gdf


def dimension_arrow_gdfs(
    gdf: gpd.GeoDataFrame, axis: dict[str, np.ndarray]
) -> dict[str, gpd.GeoDataFrame]:
    """Build one arrow/dimension-line GeoDataFrame per kind in :data:`ARROW_STYLE`.

    Args:
        gdf: The footprints these dimensions were computed for (only its
            ``crs`` is used).
        axis: Dict with the arrays needed to place each kind, as produced
            by :func:`basic_length_axis`: ``dir1``, ``dir2`` (unit
            vectors), ``L1``, ``L2``, ``a1``, ``a2`` (lengths),
            ``centroids`` (anchor for L1/L2), ``a_centers`` (anchor for
            a1/a2), and ``b_start``/``b_vector``, ``c_start``/``c_vector``
            (b's and c's own exact, non-centered spans -- see
            :func:`_setback_b_extent`/:func:`_setback_c_extent`).

    Returns:
        ``{"L1": gdf, "L2": gdf, "a1": gdf, "a2": gdf, "b": gdf, "c": gdf}``.
    """
    centered_specs = {
        "L1": (axis["centroids"], axis["dir1"] * axis["L1"][:, None]),
        "L2": (axis["centroids"], axis["dir2"] * axis["L2"][:, None]),
        "a1": (axis["a_centers"], axis["dir2"] * axis["a1"][:, None]),
        "a2": (axis["a_centers"], axis["dir1"] * axis["a2"][:, None]),
    }
    result = {
        kind: arrow_gdf(
            anchor,
            vec,
            gdf.crs,
            caps=ARROW_STYLE[kind]["caps"],
            head_frac=ARROW_STYLE[kind]["head_frac"],
            kind=kind,
            value=np.linalg.norm(vec, axis=1),
        )
        for kind, (anchor, vec) in centered_specs.items()
    }
    # b and c are NOT centered: b_start/b_vector and c_start/c_vector already
    # span their own true, generally-asymmetric extents (see
    # _setback_b_extent/_setback_c_extent) -- centering them again here
    # would shift them off their real, measured location.
    for kind in ("b", "c"):
        result[kind] = arrow_gdf(
            axis[f"{kind}_start"],
            axis[f"{kind}_vector"],
            gdf.crs,
            centered=False,
            caps=ARROW_STYLE[kind]["caps"],
            head_frac=ARROW_STYLE[kind]["head_frac"],
            kind=kind,
            value=np.linalg.norm(axis[f"{kind}_vector"], axis=1),
        )
    return result


def add_dimension_layers(
    m: "folium.Map | None",
    gdf: gpd.GeoDataFrame,
    axis: dict[str, np.ndarray],
    *,
    method: str,
    weight: int = 2,
) -> "folium.Map":
    """Add one black FancyFolium layer per L1/L2/a1/a2/b/c kind to map *m*.

    Each layer is named ``"<kind> (<method>)"`` (e.g. ``"L1 (bbox)"``) and
    starts hidden (``active=False``) so the map isn't cluttered by default;
    toggle individual kinds on from the layer control.

    Args:
        m: Existing FancyFolium/folium map to add layers to, or ``None``
            to start a new one.
        gdf: The footprints these dimensions were computed for (its
            ``crs`` is used to build the arrow layers).
        axis: ``basic_length_axis()`` output for this axis convention.
        method: Axis convention label used in each layer's name (e.g.
            ``"bbox"`` or ``"inertia"``).
        weight: Stroke width (px) for every dimension line.

    Returns:
        The map *m*, with the six new layers added.
    """
    import FancyFolium

    for kind, arrows in dimension_arrow_gdfs(gdf, axis).items():
        m = FancyFolium.vector_layer(
            gdf=arrows,
            layer_name=f"{kind} ({method})",
            column="value",
            color="black",
            color_by_column=False,
            overlay=True,
            active=False,
            legend=False,
            style={
                "stroke_color": "black",
                "weight": weight,
                "dashArray": ARROW_STYLE[kind]["dashArray"],
            },
            m=m,
        )
    return m


def plot_dimension_arrows(
    ax: "matplotlib.axes.Axes", gdf: gpd.GeoDataFrame, axis: dict[str, np.ndarray]
) -> None:
    """Draw L1/L2/a1/a2/b/c on an existing matplotlib ``ax`` -- the Part 1
    "idealized shapes" plots -- one linestyle per kind, all black except
    ``c`` (see :data:`MPL_COLOR`).

    Args:
        ax: Matplotlib axes to draw on.
        gdf: The footprints these dimensions were computed for.
        axis: ``basic_length_axis()`` output for this axis convention.
    """
    for kind, arrows in dimension_arrow_gdfs(gdf, axis).items():
        arrows.plot(
            ax=ax,
            color=MPL_COLOR[kind],
            linestyle=ARROW_STYLE[kind]["mpl"],
            linewidth=MPL_LINEWIDTH[kind],
        )


def dimension_legend_handles() -> list:
    """``Line2D`` handles for L1/L2/a1/a2/b/c, matching :func:`plot_dimension_arrows`'s
    styling -- pass to ``ax.legend(handles=...)`` on any plot built with it.
    """
    from matplotlib.lines import Line2D

    return [
        Line2D(
            [0],
            [0],
            color=MPL_COLOR[kind],
            linestyle=ARROW_STYLE[kind]["mpl"],
            linewidth=MPL_LINEWIDTH[kind],
            label=kind,
        )
        for kind in ARROW_STYLE
    ]


def direction_arrow_gdfs(
    gdf: gpd.GeoDataFrame,
    dir1: np.ndarray,
    L1: np.ndarray,
    dir2: np.ndarray,
    L2: np.ndarray,
) -> dict[str, gpd.GeoDataFrame]:
    """Build the ``{"L1": gdf, "L2": gdf}`` pair of direction-arrow GeoDataFrames
    used throughout ``direction.ipynb`` -- one real arrow per axis, anchored
    at each footprint's centroid.

    Args:
        gdf: The footprints these directions were computed for (only its
            ``crs`` and centroids are used).
        dir1: (N, 2) unit vectors for the L1 axis.
        L1: (N,) lengths along ``dir1``.
        dir2: (N, 2) unit vectors for the L2 axis.
        L2: (N,) lengths along ``dir2``.

    Returns:
        ``{"L1": gdf, "L2": gdf}`` of arrow GeoDataFrames.
    """
    centroids = np.column_stack([gdf.geometry.centroid.x, gdf.geometry.centroid.y])
    dir1 = np.asarray(dir1, dtype=float)
    dir2 = np.asarray(dir2, dtype=float)
    L1 = np.asarray(L1, dtype=float)
    L2 = np.asarray(L2, dtype=float)
    return {
        "L1": arrow_gdf(
            centroids,
            dir1 * L1[:, None],
            gdf.crs,
            caps="arrow",
            kind="L1",
            value=L1,
        ),
        "L2": arrow_gdf(
            centroids,
            dir2 * L2[:, None],
            gdf.crs,
            caps="arrow",
            kind="L2",
            value=L2,
        ),
    }


def add_direction_layers(
    m: "folium.Map | None",
    gdf: gpd.GeoDataFrame,
    dir1: np.ndarray,
    L1: np.ndarray,
    dir2: np.ndarray,
    L2: np.ndarray,
    *,
    method: str,
    weight: int = 2,
) -> "folium.Map":
    """Add one black ``"L1 (<method>)"``/``"L2 (<method>)"`` FancyFolium layer.

    Args:
        m: Existing FancyFolium/folium map to add layers to, or ``None``
            to start a new one.
        gdf: The footprints these directions were computed for.
        dir1: (N, 2) unit vectors for the L1 axis.
        L1: (N,) lengths along ``dir1``.
        dir2: (N, 2) unit vectors for the L2 axis.
        L2: (N,) lengths along ``dir2``.
        method: Axis convention label used in each layer's name (e.g.
            ``"bbox"`` or ``"inertia"``).
        weight: Stroke width (px) for both arrows.

    Returns:
        The map *m*, with the two new layers added.
    """
    import FancyFolium

    for kind, arrows in direction_arrow_gdfs(gdf, dir1, L1, dir2, L2).items():
        m = FancyFolium.vector_layer(
            gdf=arrows,
            layer_name=f"{kind} ({method})",
            column="value",
            color="black",
            color_by_column=False,
            overlay=True,
            active=False,
            legend=False,
            style={
                "stroke_color": "black",
                "weight": weight,
                "dashArray": METHOD_STYLE[method]["dashArray"],
            },
            m=m,
        )
    return m


def plot_direction_arrows(
    ax: "matplotlib.axes.Axes",
    gdf: gpd.GeoDataFrame,
    dir1: np.ndarray,
    L1: np.ndarray,
    dir2: np.ndarray,
    L2: np.ndarray,
    *,
    method: str = "bbox",
) -> None:
    """Draw dir1/dir2 as black arrows on an existing matplotlib ``ax``, one
    linestyle per *method* (so bbox vs. inertia stay visually distinct when
    plotted on the same axes).

    Args:
        ax: Matplotlib axes to draw on.
        gdf: The footprints these directions were computed for.
        dir1: (N, 2) unit vectors for the L1 axis.
        L1: (N,) lengths along ``dir1``.
        dir2: (N, 2) unit vectors for the L2 axis.
        L2: (N,) lengths along ``dir2``.
        method: Which :data:`METHOD_STYLE` linestyle to draw both arrows
            with (e.g. ``"bbox"`` or ``"inertia"``).
    """
    for arrows in direction_arrow_gdfs(gdf, dir1, L1, dir2, L2).values():
        arrows.plot(
            ax=ax, color="black", linestyle=METHOD_STYLE[method]["mpl"], linewidth=1.5
        )


def _setback_b_extent(
    gdf: gpd.GeoDataFrame,
    dir1: np.ndarray,
    dir2: np.ndarray,
    b1: np.ndarray,
    b2: np.ndarray,
    b_dir: np.ndarray,
    centroids: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """For each row, find the dominant setback piece (matched to *b1*/*b2*
    from :func:`footprint_attributes.geometry.setback_gndt_metrics`) and
    return ``(start, vector)`` so the ``b`` bar spans *exactly* that piece's
    true projected extent along *b_dir* -- no more, no less -- while passing
    through *centroids* (the same point :func:`setback_gndt_metrics` casts
    its ``c`` ray through), at whatever position along *b_dir* that point
    projects to.

    This -- not a symmetric "anchor +/- b/2" bar -- is required for
    correctness: the true setback centroid is generally *not* at the
    midpoint of the piece's own bounding box (e.g. a triangular hull-gap
    piece), so centering a length-``b`` bar there would make the bar
    overhang past the piece on one side and fall short on the other, while
    anchoring a bbox-centered bar there (an earlier version of this
    function) would draw ``b``/``c`` at a *different* point than the one
    the real ``c`` measurement was taken from. Using the true centroid as
    the (fixed) cross-axis reference and the piece's real projected
    min/max along `b_dir` as the span avoids both problems at once.
    """
    from .geometry import setback_pieces

    n = len(gdf)
    start = centroids.copy()
    vector = np.zeros((n, 2))
    for i, poly in enumerate(gdf.geometry):
        if b1[i] <= 0 and b2[i] <= 0:
            continue
        pieces = setback_pieces(poly, dir1[i], dir2[i])
        if not pieces:
            continue
        match = next(
            (
                piece
                for ext1, ext2, piece in pieces
                if np.isclose(ext1, b1[i], rtol=1e-3, atol=1e-6)
                and np.isclose(ext2, b2[i], rtol=1e-3, atol=1e-6)
            ),
            pieces[0][2],  # fallback: dominant-by-area piece
        )
        coords = shapely.get_coordinates(match)
        proj_b = coords @ b_dir[i]
        proj_min, proj_max = proj_b.min(), proj_b.max()
        anchor_proj = centroids[i] @ b_dir[i]
        start_point = centroids[i] + (proj_min - anchor_proj) * b_dir[i]
        end_point = centroids[i] + (proj_max - anchor_proj) * b_dir[i]
        start[i] = start_point
        vector[i] = end_point - start_point
    return start, vector


def _setback_c_extent(
    gdf: gpd.GeoDataFrame,
    b1: np.ndarray,
    b2: np.ndarray,
    centroids: np.ndarray,
    c_dir: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """For each row, replicate :func:`footprint_attributes.geometry.setback_gndt_metrics`'s
    own ray-cast for ``c`` and return ``(start, vector)`` for the *actual*
    solid segment it measures -- not a segment centered on the setback
    centroid.

    ``c`` is the length of whichever solid-footprint segment a ray through
    the setback centroid, cast along *c_dir*, happens to hit closest to
    that centroid. That segment is generally **not** centered on (and may
    not even contain) the centroid itself -- the centroid sits in the empty
    setback piece, not on the solid material being measured, so a ray
    through it can land anywhere along the solid segment it eventually
    crosses. Drawing ``c`` as a length-``c`` segment centered on the
    centroid (an earlier version of this function) draws the right
    *length* at the wrong *place* -- this recovers the true segment
    directly instead.
    """
    from .geometry import fill_holes

    n = len(gdf)
    start = centroids.copy()
    vector = np.zeros((n, 2))
    filled = fill_holes(gdf.geometry)
    for i, poly in enumerate(filled):
        if b1[i] <= 0 and b2[i] <= 0:
            continue
        cx, cy = centroids[i]
        perp = c_dir[i]
        minx, miny, maxx, maxy = poly.bounds
        diag = float(np.hypot(maxx - minx, maxy - miny))
        line = LineString(
            [
                (cx - perp[0] * (diag + 1), cy - perp[1] * (diag + 1)),
                (cx + perp[0] * (diag + 1), cy + perp[1] * (diag + 1)),
            ]
        )
        intersec = poly.intersection(line)
        parts = list(intersec.geoms) if hasattr(intersec, "geoms") else [intersec]
        pt_c = Point(cx, cy)
        closest = min(
            (p for p in parts if not p.is_empty),
            key=lambda g: g.centroid.distance(pt_c),
            default=None,
        )
        if closest is None:
            continue
        coords = list(closest.coords)
        start[i] = np.asarray(coords[0], dtype=float)
        vector[i] = np.asarray(coords[-1], dtype=float) - start[i]
    return start, vector


def basic_length_axis(
    gdf: gpd.GeoDataFrame,
    L1: np.ndarray,
    dir1: np.ndarray,
    L2: np.ndarray,
    dir2: np.ndarray,
) -> dict[str, np.ndarray]:
    """Compute the full L1/L2/a1/a2/b/c bundle (anchors + directions) needed
    by :func:`dimension_arrow_gdfs`/:func:`add_dimension_layers`/
    :func:`plot_dimension_arrows`, for one axis convention (``bbox`` or
    ``inertia``).

    ``b``/``c`` are drawn from *only* the single winning setback
    configuration (whichever of ``b1``/``b2`` actually set ``beta2`` --
    see :func:`footprint_attributes.geometry.setback_gndt_metrics`), each
    anchored at the setback's own centroid and pointing along its own
    axis -- drawing both ``b1`` and ``b2`` unconditionally, as earlier
    versions of this notebook did, mislabels the non-winning one as if it
    were also a real, used dimension.

    Args:
        gdf: The footprints these dimensions are computed for.
        L1: (N,) longer plan dimension per footprint.
        dir1: (N, 2) unit vectors for the L1 axis.
        L2: (N,) shorter plan dimension per footprint.
        dir2: (N, 2) unit vectors for the L2 axis.

    Returns:
        Dict with ``dir1``, ``dir2``, ``L1``, ``L2``, ``a1``, ``a2``,
        ``b``, ``c``, ``ratio``, ``centroids``, ``a_centers``,
        ``setback_centroids``, ``b_start``, ``b_vector``, ``c_start``,
        ``c_vector``, ``b_dir``, ``c_dir`` -- see :func:`dimension_arrow_gdfs`
        for how these are turned into drawable geometries.
    """
    from .geometry import main_element_a_lengths_batch, setback_gndt_metrics

    L1 = np.asarray(L1, dtype=float)
    L2 = np.asarray(L2, dtype=float)
    dir1 = np.asarray(dir1, dtype=float)
    dir2 = np.asarray(dir2, dtype=float)

    centroids = np.column_stack([gdf.geometry.centroid.x, gdf.geometry.centroid.y])

    a1, a2, a_centers = main_element_a_lengths_batch(gdf, dir1, dir2)
    a1 = np.asarray(a1, dtype=float)
    a2 = np.asarray(a2, dtype=float)
    a_centers = np.asarray(a_centers, dtype=float)

    ratio, b, c, b1, b2, cx, cy = setback_gndt_metrics(
        gdf, L1, dir1, L2, dir2, full_output=True
    )
    b = np.asarray(b, dtype=float)
    c = np.asarray(c, dtype=float)
    b1 = np.asarray(b1, dtype=float)
    b2 = np.asarray(b2, dtype=float)
    # The exact point setback_gndt_metrics casts its "c" ray through --
    # the true setback-piece centroid, NOT the centre of its bounding box.
    setback_centroids = np.column_stack([cx, cy]).astype(float)

    # b is exactly a copy of whichever of b1/b2 won (see setback_gndt_metrics);
    # recover which one, so b/c can be drawn along their own real axes
    # instead of both candidates being drawn regardless of which was used.
    b1_won = np.isclose(b, b1)
    b_dir = np.where(b1_won[:, None], dir1, dir2)
    c_dir = np.where(b1_won[:, None], dir2, dir1)  # c is measured perpendicular to b

    # b spans exactly the winning piece's true extent along b_dir (so it
    # can't overhang past the setback), passing through the same true
    # centroid the "c" ray was cast through -- see _setback_b_extent.
    b_start, b_vector = _setback_b_extent(
        gdf, dir1, dir2, b1, b2, b_dir, setback_centroids
    )
    # c likewise spans the *actual* solid segment the ray hits, not a
    # length-c segment centered on the centroid -- see _setback_c_extent.
    c_start, c_vector = _setback_c_extent(gdf, b1, b2, setback_centroids, c_dir)

    return dict(
        dir1=dir1,
        dir2=dir2,
        L1=L1,
        L2=L2,
        a1=a1,
        a2=a2,
        b=b,
        c=c,
        ratio=np.asarray(ratio, dtype=float),
        centroids=centroids,
        a_centers=a_centers,
        setback_centroids=setback_centroids,
        b_start=b_start,
        b_vector=b_vector,
        c_start=c_start,
        c_vector=c_vector,
        b_dir=b_dir,
        c_dir=c_dir,
    )
