"""Plotting helpers shared by the example notebooks (not used by the core
attribute-computation pipeline).

The one function here, :func:`arrow_gdf`, turns a set of 2-D vectors (e.g.
a building's principal axis, a contact force, a resultant) into a
GeoDataFrame of arrow-shaped ``MultiLineString`` geometries (a shaft plus a
two-line arrowhead) that can be handed straight to ``FancyFolium.vector_layer``
so the vectors show up as literal arrows on a map.
"""

from __future__ import annotations

from typing import Any

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, MultiLineString


def arrow_gdf(
    anchors: np.ndarray,
    vectors: np.ndarray,
    crs: Any,
    *,
    centered: bool = True,
    caps: str = "arrow",
    head_frac: float = 0.2,
    head_angle_deg: float = 25.0,
    min_length: float = 1e-9,
    **extra_columns: Any,
) -> gpd.GeoDataFrame:
    """Build a GeoDataFrame of arrow/dimension-line ``MultiLineString`` geometries.

    Each line is a shaft plus, depending on *caps*, either a single
    arrowhead at its tip (``caps="arrow"``) or a short perpendicular tick
    mark ("|") at *both* ends (``caps="bar"``, like a drafting dimension
    line). Either renders as more than a plain line once styled with
    ``FancyFolium.vector_layer``.

    Args:
        anchors: ``(N, 2)`` array-like of (x, y) anchor points.
        vectors: ``(N, 2)`` array-like of (vx, vy) vectors -- direction
            *and* length are both encoded here (unlike a unit direction
            vector, magnitude matters: it sets the arrow's length).
        crs: CRS to assign to the output GeoDataFrame (should match the
            CRS *anchors*/*vectors* are expressed in).
        centered: If ``True`` (default), the anchor is the arrow's
            midpoint (spans ``anchor - vector/2`` to ``anchor + vector/2``)
            -- natural for a building's own axis. If ``False``, the anchor
            is the arrow's tail (spans ``anchor`` to ``anchor + vector``)
            -- natural for a force/vector "pushing" from a point.
        caps: ``"arrow"`` (default) draws one arrowhead at the tip --
            use this for a true *direction* (e.g. ``L1``/``L2``/``dir1``).
            ``"bar"`` draws a plain shaft with a perpendicular tick at each
            end, like a dimension line -- use this for a *measurement*
            that has no direction of its own (e.g. ``a1``, ``a2``, ``b``,
            ``c``).
        head_frac: Arrowhead barb length (``caps="arrow"``) or end-tick
            length (``caps="bar"``), as a fraction of the line's own
            length.
        head_angle_deg: Half-angle (degrees) between the two arrowhead
            barbs and the shaft. Only used for ``caps="arrow"``.
        min_length: Rows whose vector magnitude is below this are dropped
            (e.g. isolated buildings with zero contact force) rather than
            drawn as a degenerate zero-length arrow.
        **extra_columns: Extra columns to attach to the output (e.g.
            ``kind=["L1", "L1", ...]``), sliced to match the rows kept
            after the ``min_length`` filter.

    Returns:
        GeoDataFrame of ``MultiLineString`` lines, plus any *extra_columns*.
    """
    if caps not in ("arrow", "bar"):
        raise ValueError(f"caps must be 'arrow' or 'bar', got {caps!r}")

    anchors = np.asarray(anchors, dtype=float)
    vectors = np.asarray(vectors, dtype=float)
    lengths = np.linalg.norm(vectors, axis=1)
    keep = lengths >= min_length

    ang = np.radians(head_angle_deg)
    cos_a, sin_a = np.cos(ang), np.sin(ang)

    def _rotate(v: np.ndarray, cos_t: float, sin_t: float) -> np.ndarray:
        """Rotate 2-D vector *v* by the angle whose cosine/sine are given."""
        return np.array([cos_t * v[0] - sin_t * v[1], sin_t * v[0] + cos_t * v[1]])

    lines = []
    for anchor, vec in zip(anchors[keep], vectors[keep]):
        if centered:
            p0 = anchor - vec / 2
            p1 = anchor + vec / 2
        else:
            p0 = anchor
            p1 = anchor + vec

        direction = vec / np.linalg.norm(vec)

        if caps == "arrow":
            head_len = np.linalg.norm(vec) * head_frac
            back1 = p1 - _rotate(direction, cos_a, sin_a) * head_len
            back2 = p1 - _rotate(direction, cos_a, -sin_a) * head_len
            lines.append(
                MultiLineString(
                    [
                        LineString([tuple(p0), tuple(p1)]),
                        LineString([tuple(p1), tuple(back1)]),
                        LineString([tuple(p1), tuple(back2)]),
                    ]
                )
            )
        else:  # caps == "bar": a tick mark perpendicular to the shaft at each end
            perp = np.array([-direction[1], direction[0]])
            tick_len = np.linalg.norm(vec) * head_frac / 2
            tick0a, tick0b = p0 - perp * tick_len, p0 + perp * tick_len
            tick1a, tick1b = p1 - perp * tick_len, p1 + perp * tick_len
            lines.append(
                MultiLineString(
                    [
                        LineString([tuple(p0), tuple(p1)]),
                        LineString([tuple(tick0a), tuple(tick0b)]),
                        LineString([tuple(tick1a), tuple(tick1b)]),
                    ]
                )
            )

    data = {}
    for name, values in extra_columns.items():
        values = np.asarray(values, dtype=object)
        if values.ndim == 0:
            values = np.full(len(anchors), values.item(), dtype=object)
        data[name] = values[keep]

    return gpd.GeoDataFrame(data, geometry=lines, crs=crs)
