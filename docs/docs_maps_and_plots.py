#!/usr/bin/env python3
"""Generate every map and plot the documentation (and README) embeds.

Outputs, all regenerated from scratch on every run:

- ``docs/figures/*.png`` -- plots: the example notebooks' charts, annotated
  diagrams on the hand-built test shapes (``footprint_attributes.testing_shapes``),
  and distribution / code-exceedance / sensitivity charts for the three
  pilot regions.
- ``docs/maps/*.jpg`` -- static maps of the three pilot regions (one panel
  per region, muted OpenStreetMap basemap, north arrow, scale bar) for every
  attribute the package computes, plus close-up "detail" maps with the
  direction / basic-length / contact-force / convex-hull constructions
  drawn on the buildings.
- ``docs/_static/maps/`` -- the interactive MapLibre + deck.gl 3D map
  embedded in the docs (via :func:`footprint_attributes.visualization.build_map`).
  It lives under ``_static`` because Sphinx only serves raw HTML/JS from
  its static path. Skip it with ``--skip-interactive``.

Hand-made paper figures (graphical abstract, the inscribed-circle/setback
construction steps, ...) are *not* generated here -- they live alongside
the generated plots in ``docs/figures/`` and are simply left untouched.

Every value drawn is computed live by this package from the raw
geometry-only pilot-region footprints under ``examples/data/``. Basemap
tiles are downloaded once from ``tile.openstreetmap.org`` (two connections
at most, per the OSM tile usage policy), cached under
``~/.cache/footprint_attributes/tiles``, and desaturated so the buildings
stay the focus -- the first run needs network access.

Run with this package's own venv (needs the ``visualization`` extra)::

    uv sync --all-groups
    uv run python docs/docs_maps_and_plots.py              # everything
    uv run python docs/docs_maps_and_plots.py --only maps  # figures | maps | interactive
"""

from __future__ import annotations

import argparse
import io
import math
import os
import time
import urllib.request
import warnings
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import colors as mcolors  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Circle, Patch, Rectangle  # noqa: E402
from PIL import Image  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

from footprint_attributes import config, direction, position, shape  # noqa: E402
from footprint_attributes.geometry import (  # noqa: E402
    centre_of_mass_and_stiffness,
    ensure_projected,
    max_inscribed_circle,
    min_bounding_box,
    setback_pieces,
)
from footprint_attributes.notebook_utils import (  # noqa: E402
    ARROW_STYLE,
    MPL_COLOR,
    MPL_LINEWIDTH,
    basic_length_axis,
    dimension_arrow_gdfs,
    dimension_legend_handles,
    direction_arrow_gdfs,
    load_city,
)
from footprint_attributes.position import contact_force_vectors  # noqa: E402
from footprint_attributes.testing_shapes import (  # noqa: E402
    asymmetric_l_shape_polygon,
    circle_polygon,
    confined_quartet_gdf,
    corner_triplet_gdf,
    gdf_of,
    isolated_building_gdf,
    l_shape_polygon,
    lateral_pair_gdf,
    rect_polygon,
    square_polygon,
    square_with_hole_polygon,
    t_shape_polygon,
    thin_cross_polygon,
    torque_triplet_gdf,
    x_shape_polygon,
)
from footprint_attributes.visualization.maps import (  # noqa: E402
    _FSI_COLORS,
    _FSI_LABELS,
    _FSI_ORDER,
    _NORM_SPECS,
    _POSITION_COLORS,
    _POSITION_LABELS,
    SHAPE_COLUMNS,
    SLENDERNESS_COLUMNS,
    _shape_index,
)
from footprint_attributes.viz import arrow_gdf  # noqa: E402

warnings.filterwarnings("ignore", category=UserWarning)

DOCS = Path(__file__).resolve().parent
ROOT = DOCS.parent
DATA = ROOT / "examples" / "data"
FIGURES = DOCS / "figures"
MAPS = DOCS / "maps"
INTERACTIVE = DOCS / "_static" / "maps"
TILE_CACHE = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
TILE_CACHE = TILE_CACHE / "footprint_attributes" / "tiles"

CITIES = {
    "guatemala": ("Guatemala City — Zona 10", "guatemala_pilot_region.gpkg"),
    "san_jose": ("San José — Mata Redonda", "san_jose_pilot_region.gpkg"),
    "santo_domingo": (
        "Santo Domingo — Ensanche Quisquella",
        "santo_domingo_pilot_region.gpkg",
    ),
}
CITY_COLORS = {
    "guatemala": "#2b6cb0",
    "san_jose": "#dd6b20",
    "santo_domingo": "#2f855a",
}
POSITION_ORDER = ["isolated", "lateral", "corner", "confined", "torque"]
EXTRA_SHAPE_COLUMNS = [
    "polsby_popper",
    "convex_hull_irregularity",
    "inertia_circle_ratio",
]

INK = "#1a202c"
MUTED = "#718096"
FILL = "#e2e8f0"
HULL = "#3182ce"
BBOX = "#dd6b20"
ATTRIBUTION = "© OpenStreetMap contributors"

plt.rcParams.update(
    {
        "figure.dpi": 110,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.15,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.labelcolor": INK,
        "axes.edgecolor": "#a0aec0",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": "#edf2f7",
        "grid.linewidth": 0.8,
        "xtick.color": "#4a5568",
        "ytick.color": "#4a5568",
        "legend.frameon": False,
        "figure.titlesize": 13,
        "figure.titleweight": "bold",
    }
)


# ─────────────────────────────────────────────────────────────────────────────
# Data
# ─────────────────────────────────────────────────────────────────────────────


def load_all() -> dict[str, gpd.GeoDataFrame]:
    """Every pilot region with every column the docs draw, in its own UTM CRS."""
    out = {}
    for city, (_, fname) in CITIES.items():
        gdf = load_city(str(DATA / fname))
        gdf = gdf.drop(columns="height", errors="ignore")
        cols = [*SHAPE_COLUMNS, *SLENDERNESS_COLUMNS, *EXTRA_SHAPE_COLUMNS]
        res = shape(gdf, columns=cols)
        for c in cols:
            gdf[c] = res[c].values
        pos = position(gdf, buffer=0.1)
        for c in pos.columns:
            if c != "geometry":
                gdf[c] = pos[c].values
        L1, dir1, L2, dir2, bearing = direction.inertia(gdf, mode="all")
        gdf["bearing_inertia"] = bearing
        L1b, _, L2b, _, bearing_b = direction.bbox(gdf, mode="all")
        gdf["bearing_bbox"] = bearing_b
        gdf["L1_bbox"], gdf["L2_bbox"] = L1b, L2b
        gdf["shape_index"] = _shape_index(gdf)
        out[city] = gdf
        print(f"  {city:14s} {len(gdf):5d} footprints")
    return out


def grade_pivot(column: str) -> tuple[float, bool]:
    """The code limit of *column* (boundary of its score-100 grade in
    :mod:`footprint_attributes.config`) and whether higher values are worse.
    """
    if column.startswith("slenderness"):
        grades = config.SLENDERNESS_LIMITS["EC8"]
    else:
        norm, param = column.split("_", 1)
        grades = getattr(config, f"{norm}_LIMITS")[param]
    if grades[0]["score"] == 100:
        return grades[0]["max"], True
    first_ok = next(i for i, g in enumerate(grades) if g["score"] == 100)
    return grades[first_ok - 1]["max"], False


def is_compliant(values: pd.Series, column: str) -> pd.Series:
    limit, worse_high = grade_pivot(column)
    ok = values <= limit if worse_high else values > limit
    return ok.where(values.notna())


NORM_LABELS = {
    name: (label, criteria) for name, label, _, _, _, criteria, _ in _NORM_SPECS
}
NORM_CRITERIA = {
    # config.py's own score-100 boundary, spelled out for each column
    "ASCE7_parallelityAngle": "ASCE 7: angle ≤ 5°",
    "CSCR2010_eccentricityRatio": "CSCR-2010: eccentricity ≤ 0.05",
    "GNDTII_beta1_mainShapeSlenderness": "GNDT-II β1 ≥ 0.8 (grade A)",
    "GNDTII_beta2_setbackRatio": "GNDT-II β2 ≤ 0.1 (grade A)",
    "GNDTII_beta4_eccentricityRatio": "GNDT-II β4 ≤ 0.2 (grade A)",
    "GNDTII_beta6_setbackSlenderness": "GNDT-II β6 ≥ 0.5 (grade A)",
}


def criteria(column: str) -> str:
    return (
        NORM_CRITERIA.get(column, NORM_LABELS[column][1])
        .replace("<=", "≤")
        .replace(">=", "≥")
    )


# ─────────────────────────────────────────────────────────────────────────────
# Basemap tiles
# ─────────────────────────────────────────────────────────────────────────────

WEB_MERCATOR_HALF = 20037508.342789244
TILE_URL = "https://tile.openstreetmap.org/{z}/{x}/{y}.png"


def _tile(z: int, x: int, y: int) -> np.ndarray | None:
    path = TILE_CACHE / str(z) / str(x) / f"{y}.png"
    if not path.exists():
        url = TILE_URL.format(z=z, x=x, y=y)
        req = urllib.request.Request(
            url, headers={"User-Agent": "footprint_attributes-docs/1.0"}
        )
        for attempt in range(3):
            try:
                with urllib.request.urlopen(req, timeout=20) as r:
                    data = r.read()
                break
            except OSError:
                time.sleep(1 + attempt)
        else:
            return None
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    return np.asarray(Image.open(path).convert("RGB"))


def add_basemap(ax, xmin, ymin, xmax, ymax, panel_px: int) -> None:
    """Draw muted OpenStreetMap tiles under *ax* for an EPSG:3857 extent."""
    span = max(xmax - xmin, ymax - ymin)
    z = int(
        np.clip(
            math.ceil(math.log2(2 * WEB_MERCATOR_HALF * panel_px / 256 / span)), 1, 19
        )
    )
    size = 2 * WEB_MERCATOR_HALF / 2**z

    def col(x):
        return int((x + WEB_MERCATOR_HALF) // size)

    def row(y):
        return int((WEB_MERCATOR_HALF - y) // size)

    xs = range(col(xmin), col(xmax) + 1)
    ys = range(row(ymax), row(ymin) + 1)
    keys = [(x, y) for y in ys for x in xs]
    with ThreadPoolExecutor(2) as pool:
        tiles = dict(zip(keys, pool.map(lambda k: _tile(z, *k), keys)))
    px = next((t.shape[0] for t in tiles.values() if t is not None), 256)
    mosaic = np.full((len(ys) * px, len(xs) * px, 3), 245, dtype=np.uint8)
    for (x, y), t in tiles.items():
        if t is not None:
            i, j = (y - ys.start) * px, (x - xs.start) * px
            mosaic[i : i + px, j : j + px] = t
    # desaturate + lighten: a quiet, light-grey street map under coloured data
    rgb = mosaic.astype(float)
    grey = rgb @ np.array([0.299, 0.587, 0.114])
    rgb = 0.25 * rgb + 0.75 * grey[..., None]
    mosaic = (255 - 0.45 * (255 - rgb)).clip(0, 255).astype(np.uint8)
    extent = (
        xs.start * size - WEB_MERCATOR_HALF,
        (xs.stop) * size - WEB_MERCATOR_HALF,
        WEB_MERCATOR_HALF - ys.stop * size,
        WEB_MERCATOR_HALF - ys.start * size,
    )
    ax.imshow(mosaic, extent=extent, interpolation="lanczos", zorder=0)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)


# ─────────────────────────────────────────────────────────────────────────────
# Map furniture
# ─────────────────────────────────────────────────────────────────────────────


def _nice(x: float) -> float:
    for n in [5, 10, 20, 25, 50, 100, 150, 200, 250, 500, 1000]:
        if n >= x:
            return n
    return 1000


def decorate_map(ax, lat: float, *, scale: bool = True) -> None:
    """North arrow, scale bar, attribution, thin frame."""
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_color("#a0aec0")
        s.set_linewidth(0.8)
    ax.annotate(
        "N",
        xy=(0.07, 0.95),
        xytext=(0.07, 0.83),
        xycoords="axes fraction",
        ha="center",
        va="center",
        fontsize=11,
        fontweight="bold",
        color=INK,
        arrowprops=dict(
            facecolor=INK, edgecolor=INK, width=3, headwidth=10, headlength=9
        ),
        zorder=10,
    )
    ax.text(
        0.01,
        0.01,
        ATTRIBUTION,
        transform=ax.transAxes,
        fontsize=5.5,
        color="#4a5568",
        ha="left",
        va="bottom",
        zorder=10,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=1),
    )
    if not scale:
        return
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    k = math.cos(math.radians(lat))  # web-mercator units -> ground metres
    ground = _nice((xmax - xmin) * k / 5)
    length = ground / k
    x1 = xmax - 0.05 * (xmax - xmin)
    x0 = x1 - length
    y0 = ymin + 0.06 * (ymax - ymin)
    h = 0.012 * (ymax - ymin)
    ax.add_patch(
        Rectangle(
            (x0 - 0.02 * (xmax - xmin), y0 - 0.035 * (ymax - ymin)),
            length + 0.04 * (xmax - xmin),
            0.09 * (ymax - ymin),
            facecolor="white",
            edgecolor="none",
            alpha=0.8,
            zorder=9,
        )
    )
    ax.add_patch(
        Rectangle((x0, y0), length, h, facecolor=INK, edgecolor=INK, zorder=10)
    )
    ax.text(
        (x0 + x1) / 2,
        y0 + 2.2 * h,
        f"{ground:g} m",
        ha="center",
        va="bottom",
        fontsize=8,
        color=INK,
        zorder=10,
    )


def square_extent(bounds, pad: float = 0.04) -> tuple[float, float, float, float]:
    xmin, ymin, xmax, ymax = bounds
    cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
    half = max(xmax - xmin, ymax - ymin) / 2 * (1 + pad)
    return cx - half, cy - half, cx + half, cy + half


def region_panels(n: int = 3, legend_height: float = 0.1, width: float = 15.0):
    fig = plt.figure(figsize=(width, width / n * (1 + legend_height) + 1.0))
    gs = fig.add_gridspec(
        2,
        n,
        height_ratios=[1, legend_height],
        hspace=0.06,
        wspace=0.03,
        left=0.01,
        right=0.99,
        top=0.86,
        bottom=0.01,
    )
    axes = [fig.add_subplot(gs[0, i]) for i in range(n)]
    lax = fig.add_subplot(gs[1, :])
    lax.set_axis_off()
    return fig, axes, lax


def save(fig, path: Path, **kw) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".jpg":
        buf = io.BytesIO()
        fig.savefig(buf, format="png", **kw)
        buf.seek(0)
        Image.open(buf).convert("RGB").save(
            path, quality=84, optimize=True, progressive=True
        )
    else:
        fig.savefig(path, **kw)
        Image.open(path).save(path, optimize=True)
    plt.close(fig)
    print(f"  wrote {path.relative_to(ROOT)}  ({path.stat().st_size / 1024:.0f} KB)")


# ─────────────────────────────────────────────────────────────────────────────
# Region maps (docs/maps)
# ─────────────────────────────────────────────────────────────────────────────


def _draw_region(ax, gdf3857, colors, lat, extent=None, panel_px=640, lw=0.25):
    extent = extent or square_extent(gdf3857.total_bounds)
    add_basemap(ax, *extent, panel_px=panel_px)
    gdf3857.plot(ax=ax, color=colors, edgecolor="#2d3748", linewidth=lw, zorder=2)
    decorate_map(ax, lat)


def region_map_categorical(data, column, palette, labels, order, title, path):
    fig, axes, lax = region_panels()
    for ax, (city, gdf) in zip(axes, data.items()):
        g = gdf.to_crs(3857)
        colors = g[column].map(palette).fillna("#a0aec0")
        _draw_region(ax, g, colors.values, gdf.to_crs(4326).geometry.centroid.y.mean())
        ax.set_title(f"{CITIES[city][0]}\n{len(g)} buildings", fontsize=10.5)
    handles = [
        Patch(facecolor=palette[k], edgecolor="#2d3748", linewidth=0.5, label=labels[k])
        for k in order
    ]
    lax.legend(
        handles=handles, loc="center", ncol=len(handles), fontsize=10, handlelength=1.6
    )
    fig.suptitle(title, y=0.98)
    save(fig, path)


def _continuous_norm(values: np.ndarray, column: str | None, vmin=None, vmax=None):
    v = values[np.isfinite(values)]
    if column in NORM_LABELS:
        limit, worse_high = grade_pivot(column)
        lo = min(np.percentile(v, 2) if len(v) else limit, limit * 0.8)
        hi = max(np.percentile(v, 98) if len(v) else limit, limit * 1.2)
        if column == "EC8_compactness":
            lo, hi = min(lo, 0.7), 1.0
        norm = mcolors.TwoSlopeNorm(vcenter=limit, vmin=lo, vmax=hi)
        cmap = plt.get_cmap("RdYlGn_r" if worse_high else "RdYlGn")
        return norm, cmap, limit
    lo = vmin if vmin is not None else np.percentile(v, 2)
    hi = vmax if vmax is not None else np.percentile(v, 98)
    return mcolors.Normalize(lo, hi), None, None


def region_map_continuous(
    data, column, title, path, *, cmap=None, vmin=None, vmax=None, unit="", ticks=None
):
    allv = np.concatenate([d[column].to_numpy(dtype=float) for d in data.values()])
    norm, auto_cmap, limit = _continuous_norm(allv, column, vmin, vmax)
    cmap = plt.get_cmap(cmap) if cmap else auto_cmap
    fig, axes, lax = region_panels(legend_height=0.13)
    for ax, (city, gdf) in zip(axes, data.items()):
        g = gdf.to_crs(3857)
        v = g[column].to_numpy(dtype=float)
        colors = [cmap(norm(x)) if np.isfinite(x) else (0.63, 0.68, 0.75, 1) for x in v]
        _draw_region(ax, g, colors, gdf.to_crs(4326).geometry.centroid.y.mean())
        sub = f"{CITIES[city][0]}"
        if limit is not None:
            ok = is_compliant(g[column], column)
            sub += f"\n{100 * (1 - ok.mean()):.0f}% of buildings exceed the code limit"
        ax.set_title(sub, fontsize=10.5)
    cax = lax.inset_axes([0.3, 0.35, 0.4, 0.3])
    cb = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal"
    )
    cb.outline.set_visible(False)
    if ticks is not None:
        cb.set_ticks(ticks)
    label = unit
    if limit is not None:
        cax.axvline(limit, color=INK, lw=2)
        label = (
            f"code limit ({criteria(column)}) marked ▮ — green complies, red exceeds"
        )
    cb.set_label(label, fontsize=9, color="#4a5568")
    fig.suptitle(title, y=0.98)
    save(fig, path)


def dense_window(
    gdf3857: gpd.GeoDataFrame, half: float
) -> tuple[float, float, float, float]:
    """A square window centred on the densest cluster of footprints."""
    c = np.column_stack([gdf3857.centroid.x, gdf3857.centroid.y])
    counts = np.array([len(x) for x in cKDTree(c).query_ball_point(c, r=half)])
    cx, cy = c[np.argmax(counts)]
    return cx - half, cy - half, cx + half, cy + half


def _clip_to(gdf, extent):
    xmin, ymin, xmax, ymax = extent
    return gdf.cx[xmin:xmax, ymin:ymax]


def detail_maps(data) -> None:
    """Close-up maps with each construction drawn on real buildings."""
    windows = {}
    for city, gdf in data.items():
        g = gdf.to_crs(3857)
        k = 1 / math.cos(math.radians(gdf.to_crs(4326).geometry.centroid.y.mean()))
        windows[city] = dense_window(g, 65 * k)

    def panels(title, draw, legend_handles, path, ncol=None, zoom=1.0):
        fig, axes, lax = region_panels(legend_height=0.09)
        for ax, (city, gdf) in zip(axes, data.items()):
            x0, y0, x1, y1 = windows[city]
            cx, cy, h = (x0 + x1) / 2, (y0 + y1) / 2, (x1 - x0) / 2 * zoom
            ext = (cx - h, cy - h, cx + h, cy + h)
            # work in UTM (metres), select the window, then draw in 3857
            utm_ext = gpd.GeoSeries.from_xy(
                [ext[0], ext[2]], [ext[1], ext[3]], crs=3857
            ).to_crs(gdf.crs)
            sel = gdf.cx[
                utm_ext.x.min() - 30 : utm_ext.x.max() + 30,
                utm_ext.y.min() - 30 : utm_ext.y.max() + 30,
            ].copy()
            add_basemap(ax, *ext, panel_px=900)
            draw(ax, sel)
            ax.set_xlim(ext[0], ext[2])
            ax.set_ylim(ext[1], ext[3])
            decorate_map(ax, gdf.to_crs(4326).geometry.centroid.y.mean())
            ax.set_title(CITIES[city][0], fontsize=10.5)
        lax.legend(
            handles=legend_handles,
            loc="center",
            ncol=ncol or len(legend_handles),
            fontsize=9.5,
            handlelength=2.6,
        )
        fig.suptitle(title, y=0.98)
        save(fig, path)

    # -- direction: bounding boxes + inertia axes, buildings coloured by bearing
    twilight = plt.get_cmap("twilight_shifted")

    def draw_direction(ax, sel):
        g = sel.to_crs(3857)
        colors = [twilight((b + 90) / 180) for b in sel["bearing_inertia"]]
        g.plot(
            ax=ax,
            color=colors,
            alpha=0.55,
            edgecolor="#2d3748",
            linewidth=0.4,
            zorder=2,
        )
        boxes = gpd.GeoSeries(
            sel.geometry.apply(lambda p: p.minimum_rotated_rectangle), crs=sel.crs
        )
        boxes.to_crs(3857).plot(
            ax=ax, facecolor="none", edgecolor=BBOX, linewidth=0.9, zorder=3
        )
        L1, d1, L2, d2, _ = direction.inertia(sel, mode="all")
        for kind, arrows in direction_arrow_gdfs(sel, d1, L1, d2, L2).items():
            arrows.to_crs(3857).plot(
                ax=ax,
                color=INK,
                linewidth=1.3 if kind == "L1" else 1.0,
                linestyle="-" if kind == "L1" else "--",
                zorder=4,
            )

    panels(
        "Building direction — minimum bounding box and principal axes of inertia",
        draw_direction,
        [
            Patch(
                facecolor="none", edgecolor=BBOX, label="minimum bounding box (bbox)"
            ),
            Line2D([0], [0], color=INK, lw=1.3, label="inertia L1 / dir1"),
            Line2D([0], [0], color=INK, lw=1.0, ls="--", label="inertia L2 / dir2"),
            Patch(facecolor=twilight(0.25), alpha=0.6, label="fill: bearing (cyclic)"),
        ],
        MAPS / "detail_direction.jpg",
    )

    # -- basic lengths: L1/L2/a1/a2/b/c dimension lines
    def draw_lengths(ax, sel):
        g = sel.to_crs(3857)
        g.plot(
            ax=ax,
            color="#edf2f7",
            alpha=0.85,
            edgecolor="#2d3748",
            linewidth=0.5,
            zorder=2,
        )
        # only the ~10 largest irregular buildings near the centre, or the
        # dimension lines of neighbouring buildings drown each other out
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        c = g.centroid
        inner = (
            c.x.between(x0 + 0.15 * (x1 - x0), x1 - 0.15 * (x1 - x0))
            & c.y.between(y0 + 0.15 * (y1 - y0), y1 - 0.15 * (y1 - y0))
            & (sel["ASCE7_setbackRatio"] > 0.05)
        )
        sel = sel[inner.values]
        sel = sel.loc[sel.area.sort_values(ascending=False).index[:10]]
        sel.to_crs(3857).plot(
            ax=ax, color="#fefcbf", edgecolor=INK, linewidth=0.9, zorder=3
        )
        L1, d1, L2, d2 = min_bounding_box(sel)
        axis = basic_length_axis(sel, L1, d1, L2, d2)
        for kind, arrows in dimension_arrow_gdfs(sel, axis).items():
            if len(arrows):
                arrows.to_crs(3857).plot(
                    ax=ax,
                    color=MPL_COLOR[kind],
                    linestyle=ARROW_STYLE[kind]["mpl"],
                    linewidth=MPL_LINEWIDTH[kind],
                    zorder=4,
                )

    panels(
        "Basic plan dimensions — L1, L2 (bbox), main element a1/a2, setback b/c",
        draw_lengths,
        [
            *dimension_legend_handles(),
            Patch(facecolor="#fefcbf", edgecolor=INK, label="measured building"),
        ],
        MAPS / "detail_basic_lengths.jpg",
        zoom=0.6,
    )

    # -- contact forces: per-wall arrows + resultant, coloured by class
    def draw_position(ax, sel):
        g = sel.to_crs(3857)
        colors = sel["relativePosition"].map(_POSITION_COLORS).fillna("#a0aec0")
        g.plot(
            ax=ax,
            color=colors.values,
            alpha=0.75,
            edgecolor="#2d3748",
            linewidth=0.4,
            zorder=2,
        )
        sel = sel.reset_index(drop=True)
        edges, res = contact_force_vectors(sel, buffer=0.1)
        if len(edges):
            a = np.stack(edges["anchor"].to_numpy())
            v = np.stack(edges["vector"].to_numpy()) * 0.3
            arrow_gdf(a, v, sel.crs, centered=False, head_frac=0.3).to_crs(3857).plot(
                ax=ax, color="#4a5568", linewidth=0.6, zorder=3
            )
        if len(res):
            a = np.stack(res["anchor"].to_numpy())
            v = np.stack(res["vector"].to_numpy()) * 0.3
            arrow_gdf(a, v, sel.crs, centered=False, head_frac=0.3).to_crs(3857).plot(
                ax=ax, color=INK, linewidth=1.6, zorder=4
            )

    panels(
        "Relative position — contact force on every shared wall and its resultant",
        draw_position,
        [
            *(
                Patch(
                    facecolor=_POSITION_COLORS[k], alpha=0.8, label=_POSITION_LABELS[k]
                )
                for k in POSITION_ORDER
            ),
            Line2D([0], [0], color="#4a5568", lw=0.8, label="wall contact force"),
            Line2D([0], [0], color=INK, lw=1.8, label="resultant"),
        ],
        MAPS / "detail_contact_forces.jpg",
    )

    # -- convex hull and setback pieces
    def draw_hull(ax, sel):
        g = sel.to_crs(3857)
        norm, cmap, _ = _continuous_norm(
            sel["ASCE7_setbackRatio"].to_numpy(float), "ASCE7_setbackRatio"
        )
        colors = [
            cmap(norm(x)) if np.isfinite(x) else "#cbd5e0"
            for x in sel["ASCE7_setbackRatio"]
        ]
        g.plot(
            ax=ax, color=colors, alpha=0.7, edgecolor="#2d3748", linewidth=0.4, zorder=2
        )
        hulls = sel.geometry.convex_hull
        gaps = gpd.GeoSeries(hulls.difference(sel.geometry), crs=sel.crs)
        gaps.to_crs(3857).plot(
            ax=ax, facecolor="none", edgecolor=HULL, hatch="////", linewidth=0, zorder=3
        )
        gpd.GeoSeries(hulls, crs=sel.crs).to_crs(3857).plot(
            ax=ax,
            facecolor="none",
            edgecolor=HULL,
            linewidth=0.9,
            linestyle="--",
            zorder=4,
        )

    panels(
        "Shape — convex hull and setback pieces (hull − footprint), fill = ASCE 7 setback ratio",
        draw_hull,
        [
            Line2D([0], [0], color=HULL, ls="--", lw=1, label="convex hull"),
            Patch(
                facecolor="none", edgecolor=HULL, hatch="////", label="setback piece"
            ),
            Patch(
                facecolor=plt.get_cmap("RdYlGn_r")(0.1), label="setback ratio ≤ 0.20"
            ),
            Patch(
                facecolor=plt.get_cmap("RdYlGn_r")(0.9), label="setback ratio > 0.20"
            ),
        ],
        MAPS / "detail_convex_hull.jpg",
    )


def all_region_maps(data) -> None:
    region_map_categorical(
        data,
        "relativePosition",
        _POSITION_COLORS,
        _POSITION_LABELS,
        POSITION_ORDER,
        "Relative position within the block",
        MAPS / "relative_position.jpg",
    )
    region_map_categorical(
        data,
        "shape_index",
        _FSI_COLORS,
        _FSI_LABELS,
        _FSI_ORDER,
        "Shape index (EC8 eccentricity · ASCE 7 setback · slenderness)",
        MAPS / "shape_index.jpg",
    )
    region_map_continuous(
        data,
        "bearing_inertia",
        "Building bearing (principal axis of inertia)",
        MAPS / "bearing.jpg",
        cmap="twilight_shifted",
        vmin=-90,
        vmax=90,
        unit="bearing of the short axis, degrees clockwise from North (cyclic: −90° ≡ +90°)",
        ticks=[-90, -45, 0, 45, 90],
    )
    l1max = float(pd.concat([d["L1_bbox"] for d in data.values()]).quantile(0.95))
    region_map_continuous(
        data,
        "L1_bbox",
        "Longer plan dimension L1 (bounding box)",
        MAPS / "L1.jpg",
        cmap="viridis",
        vmin=0,
        vmax=l1max,
        unit="L1 (m), capped at the 95th percentile",
    )
    region_map_continuous(
        data,
        "contact_confinementRatio",
        "Contact confinement ratio",
        MAPS / "contact_confinement_ratio.jpg",
        cmap="magma_r",
        vmin=0,
        vmax=1,
        unit="0 = one-sided push · 1 = fully balanced by opposite walls",
    )
    fmax = float(pd.concat([d["contact_force"] for d in data.values()]).quantile(0.95))
    region_map_continuous(
        data,
        "contact_force",
        "Contact force (resultant / √area)",
        MAPS / "contact_force.jpg",
        cmap="magma_r",
        vmin=0,
        vmax=fmax,
        unit="contact_force (height = 1), capped at the 95th percentile",
    )
    region_map_continuous(
        data,
        "polsby_popper",
        "Polsby–Popper compactness",
        MAPS / "polsby_popper.jpg",
        cmap="RdYlGn",
        vmin=0.3,
        vmax=0.8,
        unit="4πA / P²  (1 = circle, square ≈ 0.785)",
    )
    region_map_continuous(
        data,
        "convex_hull_irregularity",
        "Convex-hull irregularity",
        MAPS / "convex_hull_irregularity.jpg",
        cmap="RdYlGn_r",
        vmin=0,
        vmax=0.3,
        unit="(hull area − area) / area  (0 = convex)",
    )
    for column in [*SHAPE_COLUMNS, *SLENDERNESS_COLUMNS]:
        label = NORM_LABELS[column][0]
        region_map_continuous(data, column, label, MAPS / f"{column}.jpg")


# ─────────────────────────────────────────────────────────────────────────────
# Diagrams on test shapes (docs/figures)
# ─────────────────────────────────────────────────────────────────────────────


def _shape_ax(ax, g, title=None, face=FILL):
    g.plot(ax=ax, facecolor=face, edgecolor=INK, linewidth=1.2, zorder=2)
    ax.set_aspect("equal")
    ax.set_axis_off()
    if title:
        ax.set_title(title, fontsize=10)


def fig_direction(data) -> None:
    # 1. rotated rectangle: both methods recover the known bearing
    angles = np.arange(-80, 81, 10)
    rects = gpd.GeoDataFrame(
        {"angle": angles},
        geometry=[rect_polygon(a) for a in angles],
        crs=gdf_of(square_polygon()).crs,
    )
    bb = direction.bbox(rects)
    ine = direction.inertia(rects)
    fig, (ax0, ax) = plt.subplots(
        1, 2, figsize=(11, 4), gridspec_kw={"width_ratios": [1, 1.3]}
    )
    for a in [-60, -30, 0, 30, 60]:
        g = gdf_of(rect_polygon(a))
        g = g.set_geometry(g.translate(xoff=a * 0.9))
        g.plot(ax=ax0, facecolor=FILL, edgecolor=INK, lw=1)
        L1, d1, L2, d2, _ = direction.inertia(g, mode="all")
        for kind, arr in direction_arrow_gdfs(g, d1, L1, d2, L2).items():
            arr.plot(ax=ax0, color=INK, lw=1.2, ls="-" if kind == "L1" else "--")
        ax0.text(g.centroid.x[0], -17, f"{a}°", ha="center", fontsize=9, color=MUTED)
    ax0.set_aspect("equal")
    ax0.set_axis_off()
    ax0.set_title("20 × 10 m rectangle at known angles")
    ax.plot(
        angles, -angles, color="#e53e3e", lw=6, alpha=0.25, label="expected (−angle)"
    )
    ax.plot(angles, bb, "o-", color=BBOX, label="bbox", ms=5)
    ax.plot(angles, ine, "s--", color=INK, label="inertia", ms=4, mfc="white")
    ax.set_xlabel("true rotation angle (°)")
    ax.set_ylabel("reported bearing (°)")
    ax.set_title("Both methods recover the known bearing")
    ax.legend()
    save(fig, FIGURES / "direction_rectangle_bearing.png")

    # 2. both methods on the hand-built test shapes
    shapes = {
        "L-shape": l_shape_polygon(),
        "T-shape": t_shape_polygon(),
        "X-shape": x_shape_polygon(),
        "asymmetric L": asymmetric_l_shape_polygon(),
    }
    fig, axes = plt.subplots(1, 4, figsize=(14, 4))
    for ax, (name, poly) in zip(axes, shapes.items()):
        g = gdf_of(poly)
        _shape_ax(ax, g)
        gpd.GeoSeries([poly.minimum_rotated_rectangle]).plot(
            ax=ax, facecolor="none", edgecolor=BBOX, lw=1.2, zorder=3
        )
        L1, d1, L2, d2, bb = direction.bbox(g, mode="all")
        for kind, arr in direction_arrow_gdfs(g, d1, L1, d2, L2).items():
            arr.plot(ax=ax, color=BBOX, lw=1.6, zorder=4)
        Li1, di1, Li2, di2, bi = direction.inertia(g, mode="all")
        for kind, arr in direction_arrow_gdfs(g, di1, Li1, di2, Li2).items():
            arr.plot(ax=ax, color=INK, lw=1.4, ls="--", zorder=5)
        ax.set_title(
            f"{name}\nbearing: bbox {bb[0]:.0f}° · inertia {bi[0]:.0f}°", fontsize=10
        )
    fig.legend(
        handles=[
            Patch(facecolor="none", edgecolor=BBOX, label="minimum bounding box"),
            Line2D([0], [0], color=BBOX, lw=1.6, label="bbox axes (L1, L2)"),
            Line2D([0], [0], color=INK, lw=1.4, ls="--", label="inertia axes (L1, L2)"),
        ],
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, -0.06),
    )
    save(fig, FIGURES / "direction_methods_test_shapes.png")

    # 3. bbox vs inertia on real data
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    for ax, (city, gdf) in zip(axes, data.items()):
        diff = np.abs(((gdf["bearing_bbox"] - gdf["bearing_inertia"]) + 90) % 180 - 90)
        ax.plot([-90, 90], [-90, 90], color="#e53e3e", lw=1, ls="--", zorder=1)
        ax.scatter(
            gdf["bearing_bbox"],
            gdf["bearing_inertia"],
            s=7,
            alpha=0.45,
            color=CITY_COLORS[city],
            edgecolor="none",
            zorder=2,
        )
        ax.set_title(
            f"{CITIES[city][0]}\n{(diff < 5).mean() * 100:.0f}% agree within 5°",
            fontsize=10,
        )
        ax.set_xlabel("bearing — bbox (°)")
        ax.set_xlim(-92, 92)
        ax.set_ylim(-92, 92)
        ax.set_aspect("equal")
    axes[0].set_ylabel("bearing — inertia (°)")
    save(fig, FIGURES / "direction_bbox_vs_inertia.png")

    # 4. orientation roses
    fig, axes = plt.subplots(
        1, 3, figsize=(12, 4.4), subplot_kw={"projection": "polar"}
    )
    bins = np.radians(np.arange(0, 361, 10))
    for ax, (city, gdf) in zip(axes, data.items()):
        b = gdf["bearing_inertia"].to_numpy()
        both = np.radians(
            np.concatenate([b % 360, (b + 180) % 360])
        )  # axes are undirected
        h, _ = np.histogram(both, bins=bins)
        ax.bar(
            bins[:-1],
            h,
            width=np.radians(10),
            align="edge",
            color=CITY_COLORS[city],
            edgecolor="white",
            linewidth=0.6,
            alpha=0.9,
        )
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_yticklabels([])
        ax.set_title(CITIES[city][0], fontsize=10, pad=14)
        ax.grid(color="#e2e8f0")
    fig.suptitle(
        "Orientation of the short axis (inertia bearing) — the street grid shows through"
    )
    save(fig, FIGURES / "direction_bearing_rose.png")


def fig_position(data) -> None:
    scenarios = {
        "isolated": isolated_building_gdf(),
        "lateral": lateral_pair_gdf(),
        "corner": corner_triplet_gdf(),
        "confined": confined_quartet_gdf(),
        "torque": torque_triplet_gdf(),
    }
    fig, axes = plt.subplots(1, 5, figsize=(16, 4.2))
    for ax, (name, gdf) in zip(axes, scenarios.items()):
        gdf = gdf.reset_index(drop=True)
        res = position(gdf)
        cls = res["relativePosition"].iloc[0]
        gdf.iloc[1:].plot(ax=ax, facecolor=FILL, edgecolor="#a0aec0", lw=1)
        gdf.iloc[[0]].plot(
            ax=ax, facecolor=_POSITION_COLORS[cls], edgecolor=INK, lw=1.3, alpha=0.85
        )
        edges, resultants = contact_force_vectors(gdf)
        scale = 0.35 * math.sqrt(gdf.geometry.iloc[0].area)
        e0 = edges[edges["geom_id"] == 0]
        r0 = resultants[resultants["geom_id"] == 0]
        vmax = max([np.linalg.norm(v) for v in e0["vector"]] + [1e-9])
        if len(e0):
            arrow_gdf(
                np.stack(e0["anchor"]),
                np.stack(e0["vector"]) / vmax * scale,
                gdf.crs,
                centered=False,
                head_frac=0.3,
            ).plot(ax=ax, color="#4a5568", lw=1.2)
        if len(r0) and np.linalg.norm(r0["vector"].iloc[0]) > 1e-6:
            arrow_gdf(
                np.stack(r0["anchor"]),
                np.stack(r0["vector"]) / vmax * scale,
                gdf.crs,
                centered=False,
                head_frac=0.3,
            ).plot(ax=ax, color=INK, lw=2.4)
        ax.set_aspect("equal")
        ax.set_axis_off()
        r = res.iloc[0]
        ax.set_title(f"{name}  →  {cls}", color=_POSITION_COLORS[cls], fontsize=11)
        ax.text(
            0.5,
            -0.04,
            f"force {r['contact_force']:.2f} · confinement {r['contact_confinementRatio']:.2f}\n"
            f"angular acc. {r['contact_angularAcc']:.2f}",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=8.5,
            color="#4a5568",
        )
    fig.legend(
        handles=[
            Line2D(
                [0],
                [0],
                color="#4a5568",
                lw=1.2,
                label="unit pressure on each shared wall",
            ),
            Line2D([0], [0], color=INK, lw=2.4, label="resultant contact force"),
            Patch(facecolor=FILL, edgecolor="#a0aec0", label="neighbours"),
        ],
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, -0.1),
    )
    save(fig, FIGURES / "position_scenarios.png")

    # class shares per region
    fig, ax = plt.subplots(figsize=(10, 2.8))
    ax.grid(False)
    for i, (city, gdf) in enumerate(data.items()):
        share = gdf["relativePosition"].value_counts(normalize=True)
        left = 0
        for k in POSITION_ORDER:
            w = share.get(k, 0) * 100
            ax.barh(i, w, left=left, color=_POSITION_COLORS[k], edgecolor="white")
            if w > 4:
                ax.text(
                    left + w / 2,
                    i,
                    f"{w:.0f}%",
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=9,
                    fontweight="bold",
                )
            left += w
    ax.set_yticks(range(len(data)), [CITIES[c][0] for c in data])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("share of buildings (%)")
    ax.legend(
        handles=[
            Patch(color=_POSITION_COLORS[k], label=_POSITION_LABELS[k])
            for k in POSITION_ORDER
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, 1.28),
        ncol=5,
    )
    save(fig, FIGURES / "position_class_shares.png")

    # metric space: where the classes sit
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    allg = pd.concat([pd.DataFrame(g.drop(columns="geometry")) for g in data.values()])
    for k in POSITION_ORDER:
        s = allg[allg["relativePosition"] == k]
        axes[0].scatter(
            s["contact_confinementRatio"],
            s["contact_force"],
            s=7,
            alpha=0.5,
            color=_POSITION_COLORS[k],
            label=_POSITION_LABELS[k],
            edgecolor="none",
        )
        axes[1].scatter(
            s["contact_angle"],
            s["contact_angularAcc"],
            s=7,
            alpha=0.5,
            color=_POSITION_COLORS[k],
            edgecolor="none",
        )
    axes[0].set_xlabel("contact_confinementRatio")
    axes[0].set_ylabel("contact_force")
    axes[0].set_ylim(0, allg["contact_force"].quantile(0.99))
    axes[0].set_title("Force vs. confinement")
    axes[1].set_xlabel("contact_angle (°)")
    axes[1].set_ylabel("contact_angularAcc")
    axes[1].set_ylim(0, allg["contact_angularAcc"].quantile(0.99))
    axes[1].set_title("Angular acceleration vs. contact angle")
    axes[0].legend(markerscale=2.5, loc="upper right")
    fig.suptitle("The contact metrics behind each class — all three pilot regions")
    save(fig, FIGURES / "position_metric_space.png")


def _position_sweep(buf: float) -> dict[str, float]:
    parts = []
    for _, fname in CITIES.values():
        gdf = load_city(str(DATA / fname)).drop(columns="height", errors="ignore")
        parts.append(pd.DataFrame(position(gdf, buffer=buf).drop(columns="geometry")))
    df = pd.concat(parts, ignore_index=True)
    vc = df["relativePosition"].value_counts(normalize=True)
    out = {k: vc.get(k, 0.0) for k in POSITION_ORDER}
    out["force"] = df["contact_force"].mean()
    out["confinement"] = df["contact_confinementRatio"].mean()
    return out


def _shape_sweep(buf: float) -> dict[str, float]:
    parts = []
    for _, fname in CITIES.values():
        g = ensure_projected(
            load_city(str(DATA / fname)).drop(columns="height", errors="ignore")
        )
        if buf > 0:
            g.geometry = g.geometry.buffer(buf, join_style="mitre").buffer(
                -buf, join_style="mitre"
            )
            g = g.explode(index_parts=False).reset_index(drop=True)
            bad = (
                (g.geometry.type != "Polygon")
                | ~g.geometry.is_valid
                | g.geometry.is_empty
            )
            g = g[~bad].reset_index(drop=True)
        parts.append(
            shape(
                g,
                columns=[
                    "polsby_popper",
                    "convex_hull_irregularity",
                    "EC8_compactness",
                ],
            )
        )
    df = pd.concat([pd.DataFrame(p.drop(columns="geometry")) for p in parts])
    return df.mean().to_dict()


def fig_sensitivity() -> None:
    buffers = np.round(np.arange(0.0, 2.0 + 1e-9, 0.1), 2)
    with ProcessPoolExecutor(min(8, os.cpu_count() or 1)) as pool:
        pos = list(pool.map(_position_sweep, buffers))
        shp = list(pool.map(_shape_sweep, buffers))

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    for k in POSITION_ORDER:
        axes[0].plot(
            buffers,
            [p[k] * 100 for p in pos],
            marker=".",
            color=_POSITION_COLORS[k],
            label=_POSITION_LABELS[k],
            lw=2,
        )
    axes[0].axvline(0.1, color=MUTED, ls=":", lw=1)
    axes[0].text(
        0.12,
        axes[0].get_ylim()[1] * 0.95,
        "default 0.1 m",
        color=MUTED,
        fontsize=8.5,
        va="top",
    )
    axes[0].set_xlabel("contact buffer (m)")
    axes[0].set_ylabel("share of buildings (%)")
    axes[0].set_title("relativePosition vs. buffer")
    axes[0].legend(ncol=2)
    axes[1].plot(
        buffers,
        [p["force"] for p in pos],
        marker=".",
        color="#805ad5",
        lw=2,
        label="mean contact_force",
    )
    axes[1].plot(
        buffers,
        [p["confinement"] for p in pos],
        marker=".",
        color="#38a169",
        lw=2,
        label="mean contact_confinementRatio",
    )
    axes[1].set_xlabel("contact buffer (m)")
    axes[1].set_title("Mean contact metrics vs. buffer")
    axes[1].legend()
    save(fig, FIGURES / "position_buffer_sensitivity.png")

    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for key, color, label in [
        ("polsby_popper", "#3182ce", "polsby_popper"),
        ("convex_hull_irregularity", "#dd6b20", "convex_hull_irregularity"),
        ("EC8_compactness", "#38a169", "EC8_compactness"),
    ]:
        ax.plot(
            buffers, [s[key] for s in shp], marker=".", lw=2, color=color, label=label
        )
    ax.set_xlabel("smoothing buffer b (m): buffer(+b).buffer(−b)")
    ax.set_ylabel("mean over all footprints")
    ax.set_title("Shape indices vs. digitisation smoothing")
    ax.legend()
    save(fig, FIGURES / "shape_buffer_sensitivity.png")


def fig_shape(data) -> None:
    shapes = {
        "square": square_polygon(),
        "rectangle": rect_polygon(),
        "L-shape": l_shape_polygon(),
        "with hole": square_with_hole_polygon(),
        "circle": circle_polygon(),
        "cross": thin_cross_polygon(),
    }
    # 1. gallery with convex hull + inscribed circle + indices
    rows = []
    fig, axes = plt.subplots(1, 6, figsize=(16, 3.9))
    for ax, (name, poly) in zip(axes, shapes.items()):
        g = gdf_of(poly)
        _shape_ax(ax, g)
        gpd.GeoSeries([poly.convex_hull]).plot(
            ax=ax, facecolor="none", edgecolor=HULL, lw=1.2, ls="--", zorder=3
        )
        gap = gpd.GeoSeries([poly.convex_hull.difference(poly)])
        if not gap.iloc[0].is_empty:
            gap.plot(
                ax=ax, facecolor="none", edgecolor=HULL, hatch="////", lw=0, zorder=3
            )
        x, y, r = max_inscribed_circle(poly)
        ax.add_patch(
            Circle((x, y), r, facecolor="none", edgecolor="#d53f8c", lw=1.1, zorder=4)
        )
        pp = shape.polsby_popper(g)[0]
        ch = shape.convex_hull_irregularity(g)[0]
        ic = shape.inertia_circle_ratio(g)[0]
        rows.append(
            {
                "shape": name,
                "polsby_popper": pp,
                "inertia_circle_ratio": ic,
                "convex_hull_irregularity": ch,
            }
        )
        ax.set_title(f"{name}\nPP {pp:.2f} · hull irr. {ch:.2f}", fontsize=10)
    fig.legend(
        handles=[
            Line2D([0], [0], color=HULL, ls="--", label="convex hull"),
            Patch(
                facecolor="none", edgecolor=HULL, hatch="////", label="hull − footprint"
            ),
            Line2D([0], [0], color="#d53f8c", label="largest inscribed circle"),
        ],
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, -0.08),
    )
    save(fig, FIGURES / "shape_idealized_shapes.png")

    df = pd.DataFrame(rows).set_index("shape")
    fig, ax = plt.subplots(figsize=(10, 3.8))
    x = np.arange(len(df))
    for i, (col, color) in enumerate(
        [
            ("polsby_popper", "#3182ce"),
            ("inertia_circle_ratio", "#805ad5"),
            ("convex_hull_irregularity", "#dd6b20"),
        ]
    ):
        bars = ax.bar(x + (i - 1) * 0.26, df[col], width=0.26, color=color, label=col)
        ax.bar_label(bars, fmt="%.2f", fontsize=7.5, padding=2, color="#4a5568")
    ax.set_xticks(x, df.index)
    ax.set_ylim(0, 1.15)
    ax.set_title("Code-independent indices on the idealised shapes")
    ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.02))
    save(fig, FIGURES / "shape_idealized_indices.png")

    # 2. hollow-box model: centre of mass vs centre of stiffness
    ecc_shapes = {
        "rectangle": rect_polygon(),
        "L-shape": l_shape_polygon(),
        "asymmetric L": asymmetric_l_shape_polygon(),
        "T-shape": t_shape_polygon(),
    }
    fig, axes = plt.subplots(1, 4, figsize=(14, 4))
    for ax, (name, poly) in zip(axes, ecc_shapes.items()):
        g = gdf_of(poly)
        _shape_ax(ax, g)
        cm, cs = centre_of_mass_and_stiffness(g)
        ax.plot(*cm[0], "o", color="#e53e3e", ms=8, zorder=5)
        ax.plot(*cs[0], "X", color="#3182ce", ms=9, zorder=5)
        ax.plot([cm[0][0], cs[0][0]], [cm[0][1], cs[0][1]], color=INK, lw=1, zorder=4)
        ec8 = shape.EC8(g)["EC8_eccentricityRatio"].iloc[0]
        cscr = shape.CSCR2010(g)["CSCR2010_eccentricityRatio"].iloc[0]
        ax.set_title(f"{name}\nEC8 e/r {ec8:.3f} · CSCR e/l {cscr:.3f}", fontsize=10)
    fig.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                color="#e53e3e",
                lw=0,
                ms=8,
                label="centre of mass (slab + walls)",
            ),
            Line2D(
                [0],
                [0],
                marker="X",
                color="#3182ce",
                lw=0,
                ms=9,
                label="centre of stiffness (walls)",
            ),
            Line2D([0], [0], color=INK, lw=1, label="eccentricity e"),
        ],
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, -0.07),
    )
    save(fig, FIGURES / "shape_centre_of_mass_stiffness.png")

    # 3. distributions of every code metric, all regions, with the code limit
    cols = [*SHAPE_COLUMNS, *SLENDERNESS_COLUMNS]
    ncol = 5
    nrow = math.ceil(len(cols) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(17, 3.1 * nrow))
    for ax, col in zip(axes.flat, cols):
        allv = pd.concat([d[col] for d in data.values()]).dropna()
        limit, worse_high = grade_pivot(col)
        hi = max(allv.quantile(0.98), limit * 1.3)
        lo = min(allv.quantile(0.0), limit * 0.7) if not worse_high else 0
        if col == "GNDTII_beta6_setbackSlenderness":  # c/b has a long tail
            lo, hi = 0, 3
        bins = np.linspace(lo, hi, 30)
        for city, gdf in data.items():
            ax.hist(
                gdf[col].dropna().clip(lo, hi),
                bins=bins,
                histtype="step",
                lw=1.6,
                color=CITY_COLORS[city],
                density=True,
            )
        ax.axvspan(
            limit,
            hi,
            color="#fed7d7" if worse_high else "#c6f6d5",
            alpha=0.5,
            zorder=0,
            lw=0,
        )
        ax.axvspan(
            lo,
            limit,
            color="#c6f6d5" if worse_high else "#fed7d7",
            alpha=0.5,
            zorder=0,
            lw=0,
        )
        ax.axvline(limit, color=INK, lw=1.2, ls="--")
        ax.set_xlim(lo, hi)
        ax.set_yticks([])
        ax.set_title(NORM_LABELS[col][0], fontsize=9.5)
        ax.text(
            0.98,
            0.95,
            criteria(col),
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=7.5,
            color="#4a5568",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1),
        )
    for ax in list(axes.flat)[len(cols) :]:
        ax.set_visible(False)
    fig.legend(
        handles=[
            *(
                Line2D([0], [0], color=CITY_COLORS[c], lw=2, label=CITIES[c][0])
                for c in data
            ),
            Patch(facecolor="#c6f6d5", label="complies"),
            Patch(facecolor="#fed7d7", label="exceeds the code limit"),
        ],
        loc="lower center",
        ncol=5,
        bbox_to_anchor=(0.5, -0.03),
    )
    fig.suptitle("Distribution of every code shape metric in the pilot regions")
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    save(fig, FIGURES / "shape_metric_distributions.png")

    # 4. share exceeding each code limit
    fig, ax = plt.subplots(figsize=(10, 6.4))
    y = np.arange(len(cols))
    for i, (city, gdf) in enumerate(data.items()):
        share = [100 * (1 - is_compliant(gdf[c], c).mean()) for c in cols]
        ax.barh(
            y + (i - 1) * 0.27,
            share,
            height=0.27,
            color=CITY_COLORS[city],
            label=CITIES[city][0],
        )
    ax.set_yticks(
        y,
        [f"{NORM_LABELS[c][0]}   ({criteria(c).split(': ')[-1]})" for c in cols],
        fontsize=9,
    )
    ax.invert_yaxis()
    ax.set_xlabel("buildings exceeding the code limit (%)")
    ax.set_xlim(0, 100)
    ax.legend(loc="lower right")
    ax.set_title("How often each code limit is exceeded")
    save(fig, FIGURES / "shape_code_exceedance.png")

    # 5. shape index shares
    fig, ax = plt.subplots(figsize=(10, 2.8))
    ax.grid(False)
    for i, (city, gdf) in enumerate(data.items()):
        share = gdf["shape_index"].value_counts(normalize=True)
        left = 0
        for k in _FSI_ORDER:
            w = share.get(k, 0) * 100
            ax.barh(i, w, left=left, color=_FSI_COLORS[k], edgecolor="white")
            if w > 4:
                ax.text(
                    left + w / 2,
                    i,
                    f"{w:.0f}%",
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=9,
                    fontweight="bold",
                )
            left += w
    ax.set_yticks(range(len(data)), [CITIES[c][0] for c in data])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("share of buildings (%)")
    ax.legend(
        handles=[Patch(color=_FSI_COLORS[k], label=_FSI_LABELS[k]) for k in _FSI_ORDER],
        loc="upper center",
        bbox_to_anchor=(0.5, 1.28),
        ncol=4,
    )
    save(fig, FIGURES / "shape_index_shares.png")


def fig_basic_lengths(data) -> None:
    shapes = {
        "L-shape": l_shape_polygon(),
        "T-shape": t_shape_polygon(),
        "X-shape": x_shape_polygon(),
        "asymmetric L": asymmetric_l_shape_polygon(),
    }
    for method in ["bbox", "inertia"]:
        fig, axes = plt.subplots(1, 4, figsize=(15, 4.3))
        rows = []
        for ax, (name, poly) in zip(axes, shapes.items()):
            g = gdf_of(poly)
            _shape_ax(ax, g)
            if method == "bbox":
                L1, d1, L2, d2 = min_bounding_box(g)
                gpd.GeoSeries([poly.minimum_rotated_rectangle]).plot(
                    ax=ax, facecolor="none", edgecolor=BBOX, lw=1, zorder=3
                )
            else:
                L1, d1, L2, d2, _ = direction.inertia(g, mode="all")
            for _, _, piece in setback_pieces(
                poly, np.asarray(d1)[0], np.asarray(d2)[0]
            ):
                gpd.GeoSeries([piece]).plot(
                    ax=ax,
                    facecolor="#bee3f8",
                    edgecolor=HULL,
                    lw=0.6,
                    hatch="///",
                    zorder=1,
                )
            x, y, r = max_inscribed_circle(poly)
            ax.add_patch(
                Circle(
                    (x, y), r, facecolor="none", edgecolor="#d53f8c", lw=0.9, zorder=3
                )
            )
            axis = basic_length_axis(g, L1, d1, L2, d2)
            for kind, arrows in dimension_arrow_gdfs(g, axis).items():
                if len(arrows):
                    arrows.plot(
                        ax=ax,
                        color=MPL_COLOR[kind],
                        linestyle=ARROW_STYLE[kind]["mpl"],
                        linewidth=MPL_LINEWIDTH[kind],
                        zorder=5,
                    )
            rows.append(axis)
            ax.set_title(
                f"{name}\nL1 {axis['L1'][0]:.1f} · L2 {axis['L2'][0]:.1f} · a {axis['a1'][0]:.1f}/"
                f"{axis['a2'][0]:.1f}\nb {axis['b'][0]:.1f} · c {axis['c'][0]:.1f} · b/L {axis['ratio'][0]:.2f}",
                fontsize=9.5,
            )
        handles = [
            *dimension_legend_handles(),
            Patch(
                facecolor="#bee3f8", edgecolor=HULL, hatch="///", label="setback piece"
            ),
            Line2D([0], [0], color="#d53f8c", label="inscribed circle"),
        ]
        fig.legend(
            handles=handles,
            loc="lower center",
            ncol=8,
            bbox_to_anchor=(0.5, -0.07),
            handlelength=2.6,
        )
        conv = (
            "minimum bounding box" if method == "bbox" else "principal axes of inertia"
        )
        fig.suptitle(f"Basic lengths — axis convention: {conv}", y=1.04)
        save(fig, FIGURES / f"basic_lengths_test_shapes_{method}.png")

    # GNDT ratios on real data
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.6))
    for ax, col in zip(
        axes,
        [
            "GNDTII_beta1_mainShapeSlenderness",
            "GNDTII_beta2_setbackRatio",
            "GNDTII_beta6_setbackSlenderness",
        ],
    ):
        limit, worse_high = grade_pivot(col)
        bins = np.linspace(0, 1, 26)
        for city, gdf in data.items():
            ax.hist(
                gdf[col].dropna().clip(0, 1),
                bins=bins,
                histtype="step",
                lw=1.8,
                color=CITY_COLORS[city],
                density=True,
                label=CITIES[city][0],
            )
        ax.axvline(limit, color=INK, ls="--", lw=1.2)
        ax.set_title(NORM_LABELS[col][0])
        ax.text(
            0.98,
            0.95,
            criteria(col),
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8,
            color="#4a5568",
        )
        ax.set_yticks([])
    axes[0].legend(loc="upper left", fontsize=8)
    save(fig, FIGURES / "basic_lengths_gndt_ratios.png")


# ─────────────────────────────────────────────────────────────────────────────
# Interactive map (docs/_static/maps)
# ─────────────────────────────────────────────────────────────────────────────


def interactive_map() -> None:
    from footprint_attributes.visualization import build_map

    datasets = {
        city: gpd.read_file(DATA / fname) for city, (_, fname) in CITIES.items()
    }
    build_map(
        datasets,
        INTERACTIVE,
        labels={
            city: label.replace("City ", "") for city, (label, _) in CITIES.items()
        },
        default_dataset="guatemala",
        title="Pilot regions — footprint_attributes",
    )
    add_embed_params(INTERACTIVE / "main.js")
    print(f"  wrote {INTERACTIVE.relative_to(ROOT)}/")


# The docs embed the *same* viewer (one copy of its ~15 MB of data) in
# several places, each opened on the attribute/overlays that section is
# about, via URL parameters:
#
#   index.html?dataset=san_jose&attribute=relativePosition&overlays=position_arrows&zoom=17
#
#   dataset    guatemala | san_jose | santo_domingo
#   attribute  any "Color by" attribute name (e.g. shape_index, EC8_compactness)
#   overlays   comma-separated overlay ids (convex_hull, bounding_box,
#              inertia_axis, basic_lengths, position_arrows)
#   zoom       MapLibre zoom level (default: fit the whole dataset)
#   view       2d to start flat (ground-level overlays are hidden by 3D
#              extrusions at street zoom)
#   tour       1 to auto-play the tour; default 1 without parameters, 0 with
#
# FancyFolium.deck3d's viewer has no such hook, so it is appended to the
# generated main.js here. The anchors are asserted, so a FancyFolium change
# that moves them fails loudly instead of silently dropping the feature.
_EMBED_JS = """
// ---------------------------------------------------------------------------
// Docs embed parameters (appended by footprint_attributes'
// docs/docs_maps_and_plots.py): ?dataset=&attribute=&overlays=&zoom=&view=&tour=
const EMBED = (() => {
  const keys = ["dataset", "attribute", "overlays", "zoom", "view", "tour"];
  const any = keys.some((k) => URL_PARAMS.has(k));
  return { any, tour: URL_PARAMS.has("tour") ? URL_PARAMS.get("tour") === "1" : !any };
})();

async function applyEmbedParams() {
  const datasetId = URL_PARAMS.get("dataset");
  if (datasetId && DATASETS[datasetId] && datasetId !== state.datasetId) {
    await setDataset(datasetId, { fromShowcase: true });
  }
  const overlays = (URL_PARAMS.get("overlays") || "").split(",").filter((id) => id in state.overlaysActive);
  if (overlays.length) {
    overlays.forEach((id) => (state.overlaysActive[id] = true));
    renderOverlayCheckboxes();
    document.getElementById("controls-fields").classList.remove("hidden");
    document.getElementById("controls-toggle").classList.remove("collapsed");
  }
  const attribute = ATTRIBUTES.find((a) => a.name === URL_PARAMS.get("attribute"));
  if (attribute) {
    showcaseAttributeIndex = ATTRIBUTES.indexOf(attribute);
    setAttribute(attribute, { fromShowcase: true });
  }
  const zoom = parseFloat(URL_PARAMS.get("zoom"));
  if (Number.isFinite(zoom)) map.jumpTo({ zoom });
  if (URL_PARAMS.get("view") === "2d") toggle3D(false);
  if (EMBED.tour) startShowcase();
  else renderLayer();
}
"""
_EMBED_PATCHES = [
    # bootstrap(): honour the parameters instead of always starting the tour
    (
        '  startShowcase();\n}\n\nmap.on("load"',
        '  await applyEmbedParams();\n}\n\nmap.on("load"',
    ),
    # a pinned view must not drift back into the tour after 30 s idle
    (
        "if (resumeAfterIdle) showcaseIdleTimer = setTimeout(",
        "if (resumeAfterIdle && EMBED.tour) showcaseIdleTimer = setTimeout(",
    ),
]


def add_embed_params(main_js: Path) -> None:
    src = main_js.read_text(encoding="utf-8")
    for old, new in _EMBED_PATCHES:
        if src.count(old) != 1:
            raise RuntimeError(
                f"{main_js}: expected exactly one {old!r} -- FancyFolium.deck3d's "
                "template changed; update _EMBED_PATCHES in docs_maps_and_plots.py"
            )
        src = src.replace(old, new)
    main_js.write_text(src + _EMBED_JS, encoding="utf-8")


# ─────────────────────────────────────────────────────────────────────────────


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--only", choices=["figures", "maps", "interactive"], action="append"
    )
    parser.add_argument("--skip-interactive", action="store_true")
    args = parser.parse_args()
    todo = set(args.only or ["figures", "maps", "interactive"])
    if args.skip_interactive:
        todo.discard("interactive")

    t0 = time.time()
    if todo & {"figures", "maps"}:
        print("computing attributes for the pilot regions ...")
        data = load_all()
    if "figures" in todo:
        print("plots -> docs/figures")
        fig_direction(data)
        fig_position(data)
        fig_shape(data)
        fig_basic_lengths(data)
        fig_sensitivity()
    if "maps" in todo:
        print("maps -> docs/maps")
        all_region_maps(data)
        detail_maps(data)
    if "interactive" in todo:
        print("interactive map -> docs/_static/maps")
        interactive_map()
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
