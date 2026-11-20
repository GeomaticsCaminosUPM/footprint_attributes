"""Static-figure style shared by this package's maps and plots.

The look follows the figures of the structural-system paper built on this
package: one sans-serif font (Liberation Sans, falling back to Arial,
Helvetica and DejaVu Sans), dark-gray ink instead of pure black, no top or
right spines, no grid, frameless legends with square swatches, thin light
gray map frames, and a plain line scale bar with a white halo. The block
position classes keep the paper's colors.

Needs ``matplotlib`` (the ``visualization`` extra). Typical use::

    import matplotlib.pyplot as plt
    from footprint_attributes.visualization import style

    style.apply_style()
    fig, ax = plt.subplots()
    gdf.plot(ax=ax, color=gdf["blockPosition"].map(style.POSITION_COLORS))
    style.map_frame(ax)
    style.scale_bar(ax)
    ax.legend(handles=style.category_handles(style.POSITION_COLORS, style.POSITION_LABELS),
              **style.LEGEND_KW)
"""

from __future__ import annotations

import math

#: Ink of text, outlines and arrows (dark gray, as in the paper).
INK = "#262626"
#: Secondary text (captions, attribution).
MUTED = "#6B6B6B"
#: Grid lines, when a plot needs one.
GRID = "#E6E6E6"
#: Light fill for context footprints.
LIGHT = "#EDEDED"
#: Thin frame around map panels.
FRAME = "#BDBDBD"
#: Blue ramp, light to dark: the paper's single accent hue.
BLUES = ["#F2F6FA", "#D6E4F0", "#A9C6DF", "#6F9CC3", "#3E6F9E", "#1F4E79"]
ACCENT = BLUES[5]

#: Block position classes, in the paper's colors (distinct muted hues).
POSITION_COLORS = {
    "isolated": "#E4E4E4",
    "lateral": "#A8CBE6",
    "corner": "#BCA9D8",
    "confined": "#7DBBAE",
    "torque": "#DE95A3",
    "unlabeled": "#FFFFFF",
}
POSITION_ORDER = ["isolated", "lateral", "corner", "confined", "torque", "unlabeled"]
POSITION_LABELS = {k: k.capitalize() for k in POSITION_ORDER}

FONT_FAMILY = ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"]

#: Keyword arguments for a frameless legend with square swatches.
LEGEND_KW = dict(
    frameon=False,
    handlelength=1.2,
    handleheight=1.0,
    handletextpad=0.5,
    columnspacing=1.4,
    labelspacing=0.5,
    borderaxespad=0.0,
)


def apply_style(font_size: float = 10.0) -> None:
    """Set matplotlib's rcParams to the paper style.

    Args:
        font_size: Base font size in points. The paper uses 7.5 pt on a
            5.3 in wide page; wider figures (web images) need a larger base.
    """
    import matplotlib

    matplotlib.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": FONT_FAMILY,
            "font.size": font_size,
            "axes.titlesize": font_size + 0.5,
            "axes.titleweight": "bold",
            "axes.titlelocation": "left",
            "axes.labelsize": font_size,
            "xtick.labelsize": font_size - 0.5,
            "ytick.labelsize": font_size - 0.5,
            "legend.fontsize": font_size - 0.5,
            "legend.frameon": False,
            "figure.titlesize": font_size + 2,
            "figure.titleweight": "bold",
            "text.color": INK,
            "axes.edgecolor": INK,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
            "grid.color": GRID,
            "grid.linewidth": 0.5,
            "axes.axisbelow": True,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def halo(width: float = 2.2):
    """White outline that keeps text and lines readable over a basemap."""
    from matplotlib import patheffects

    return [patheffects.withStroke(linewidth=width, foreground="white")]


def map_frame(ax) -> None:
    """No ticks, a thin light gray frame on all four sides."""
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.6)
        spine.set_color(FRAME)


def nice_length(target_m: float) -> float:
    """The round distance (1, 2, 2.5, 5 x 10^k metres) closest below *target_m*."""
    if target_m <= 0:
        return 1.0
    k = 10 ** math.floor(math.log10(target_m))
    for step in (5, 2.5, 2, 1):
        if step * k <= target_m:
            return step * k
    return k


def scale_bar(
    ax,
    length_m: float | None = None,
    *,
    units_per_metre: float = 1.0,
    fraction: float = 0.25,
    fontsize: float | None = None,
) -> float:
    """Plain line scale bar in the lower left corner, labeled above.

    Args:
        ax: Map axes (its limits must already be set).
        length_m: Ground length of the bar. Default: a round length close
            to *fraction* of the panel width.
        units_per_metre: Axis units per ground metre (1 for a projected
            metric CRS, ``1 / cos(latitude)`` for Web Mercator).
        fraction: Target share of the panel width when *length_m* is None.
        fontsize: Label size (default: the legend font size).

    Returns:
        The ground length drawn, in metres.
    """
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    if length_m is None:
        length_m = nice_length((x1 - x0) * fraction / units_per_metre)
    length = length_m * units_per_metre
    xs, ys = x0 + 0.05 * (x1 - x0), y0 + 0.05 * (y1 - y0)
    ax.plot(
        [xs, xs + length],
        [ys, ys],
        color=INK,
        lw=1.6,
        solid_capstyle="butt",
        zorder=20,
        path_effects=halo(3.5),
    )
    ax.text(
        xs + length / 2,
        ys + 0.02 * (y1 - y0),
        f"{length_m:g} m",
        ha="center",
        va="bottom",
        fontsize=fontsize,
        color=INK,
        zorder=20,
        path_effects=halo(),
    )
    return length_m


def north_arrow(
    ax, *, x: float = 0.95, y: float = 0.95, fontsize: float | None = None
) -> None:
    """Small "N" with an arrow, in axes fraction (default: upper right)."""
    ax.annotate(
        "N",
        xy=(x, y),
        xytext=(x, y - 0.09),
        xycoords="axes fraction",
        ha="center",
        va="center",
        fontsize=fontsize,
        fontweight="bold",
        color=INK,
        arrowprops=dict(
            arrowstyle="-|>", color=INK, lw=1.2, mutation_scale=10, shrinkA=7, shrinkB=0
        ),
        zorder=20,
        path_effects=halo(),
    )


def category_handles(
    colors: dict,
    labels: dict,
    order=None,
    *,
    edgecolor: str = INK,
    linewidth: float = 0.4,
):
    """Square legend swatches, one per category in *order* (default: *labels*)."""
    from matplotlib.lines import Line2D

    return [
        Line2D(
            [0],
            [0],
            marker="s",
            linestyle="none",
            markersize=9,
            markerfacecolor=colors[k],
            markeredgecolor=edgecolor,
            markeredgewidth=linewidth,
            label=labels[k],
        )
        for k in (order or labels)
    ]
