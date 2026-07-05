"""
Low-level geometry primitives.

All functions here operate on Shapely objects or plain NumPy arrays.
They have no dependency on the rest of the package.

Naming conventions used throughout the package
-----------------------------------------------
- dir1, dir2  : unit vectors (length 1), dir1 always along L1 (longer side)
- L1          : longer plan dimension (m)
- L2          : shorter plan dimension (m)
- I1          : principal moment of inertia associated with dir1 (larger I)
- I2          : principal moment of inertia associated with dir2 (smaller I)
- bearing     : angle from geographic North (UTM +y) to dir2 (weak axis),
                clockwise positive, in degrees, range [-90, 90]
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from packaging.version import Version
from shapely.geometry import LineString, Point, Polygon


# ─────────────────────────────────────────────────────────────────────────────
# CRS helpers
# ─────────────────────────────────────────────────────────────────────────────


def ensure_projected(
    gdf: gpd.GeoDataFrame | gpd.GeoSeries,
) -> gpd.GeoDataFrame | gpd.GeoSeries:
    """Return *gdf* in a projected CRS (UTM estimate if needed).

    Args:
        gdf: Input GeoDataFrame or GeoSeries.

    Returns:
        The same object guaranteed to have a projected CRS.
    """
    if not gdf.crs.is_projected:
        return gdf.to_crs(gdf.geometry.estimate_utm_crs())
    return gdf


def to_gdf(geoms: gpd.GeoDataFrame | gpd.GeoSeries) -> gpd.GeoDataFrame:
    """Ensure *geoms* is a GeoDataFrame with a clean integer index.

    Args:
        geoms: GeoDataFrame or GeoSeries.

    Returns:
        GeoDataFrame with integer index starting at 0.
    """
    geoms = geoms.copy().reset_index(drop=True)
    if isinstance(geoms, gpd.GeoSeries):
        geoms = gpd.GeoDataFrame(geometry=geoms, crs=geoms.crs)
    return geoms


def fill_holes(gs: gpd.GeoSeries) -> gpd.GeoSeries:
    """Return a GeoSeries with interior rings (holes / courtyards) removed.

    Args:
        gs: GeoSeries of Polygon geometries.

    Returns:
        GeoSeries where every geometry is the exterior ring only.
    """
    return gs.apply(lambda g: Polygon(g.exterior))


# ─────────────────────────────────────────────────────────────────────────────
# Angle helpers (scalar, operate on numpy arrays)
# ─────────────────────────────────────────────────────────────────────────────


def angle_between_0_90(v0: np.ndarray, v1: np.ndarray) -> float:
    """Unsigned angle between two 2-D vectors, folded to [0, π/2].

    Returns 0 when either vector is zero.

    Args:
        v0: First 2-D vector (need not be unit).
        v1: Second 2-D vector (need not be unit).

    Returns:
        Angle in radians, range [0, π/2].
    """
    n0, n1 = np.linalg.norm(v0), np.linalg.norm(v1)
    if n0 == 0 or n1 == 0:
        return 0.0
    dot = float(np.clip(np.dot(v0 / n0, v1 / n1), -1.0, 1.0))
    return float(np.arccos(abs(dot)))


def angle_signed(v0: np.ndarray, v1: np.ndarray) -> float:
    """Signed angle from *v0* to *v1*, range (-π, π].

    Args:
        v0: Reference vector.
        v1: Target vector.

    Returns:
        Signed angle in radians.
    """
    _v0 = v0 / np.linalg.norm(v0)
    _v1 = v1 / np.linalg.norm(v1)
    dot = float(np.clip(np.dot(_v0, _v1), -1.0, 1.0))
    cross = float(_v0[0] * _v1[1] - _v0[1] * _v1[0])
    return float(np.arctan2(cross, dot))


def bearing_from_dir(dir_vec: np.ndarray) -> float:
    """Convert a UTM 2-D direction vector to a geographic bearing (degrees).

    Bearing is measured clockwise from geographic North (+y in UTM) to the
    vector, restricted to [-90, 90] by folding (axes are undirected).

    Args:
        dir_vec: 2-D unit vector (x, y) in UTM coordinates.

    Returns:
        Bearing in degrees, range [-90, 90].
    """
    # arctan2(x, y) gives clockwise angle from +y (North)
    angle_rad = np.arctan2(dir_vec[0], dir_vec[1])
    # Normalise to [-π, π]
    angle_rad = (angle_rad + np.pi) % (2 * np.pi) - np.pi
    # Fold to [-π/2, π/2]  (axis symmetry)
    if angle_rad > np.pi / 2:
        angle_rad -= np.pi
    if angle_rad < -np.pi / 2:
        angle_rad += np.pi
    return float(np.degrees(angle_rad))


def normalize_vec(v: np.ndarray) -> np.ndarray:
    """Return a unit vector in the direction of *v*.

    Args:
        v: 2-D vector (may have any magnitude ≠ 0).

    Returns:
        Unit vector.

    Raises:
        ValueError: If *v* is the zero vector.
    """
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError("Cannot normalise a zero vector.")
    return v / n


# ─────────────────────────────────────────────────────────────────────────────
# Inertia calculations (single polygon → scalar)
# ─────────────────────────────────────────────────────────────────────────────


def validate_geodataframe(gdf: gpd.GeoDataFrame, *, context: str = "") -> None:
    """Raise a descriptive error if *gdf* contains invalid or multi-part geometries.

    All public-facing functions in this package call this before doing any
    computation.  Catching problems early produces a clear error message rather
    than a confusing shape mismatch or silent wrong result later.

    Args:
        gdf:     GeoDataFrame to check.
        context: Optional caller name shown in the error message.

    Raises:
        ValueError: If any geometry is None, empty, invalid (self-intersecting
            etc.), or a MultiPolygon / GeometryCollection.
    """
    prefix = f"[{context}] " if context else ""

    geoms = gdf.geometry

    # --- null / missing ---
    null_mask = geoms.isna()
    if null_mask.any():
        bad = list(gdf.index[null_mask])
        raise ValueError(
            f"{prefix}Found {null_mask.sum()} null geometry/geometries at "
            f"index positions {bad[:10]}{'…' if len(bad) > 10 else ''}. "
            "Remove or repair before calling footprint_attributes functions."
        )

    # --- empty ---
    empty_mask = geoms.is_empty
    if empty_mask.any():
        bad = list(gdf.index[empty_mask])
        raise ValueError(
            f"{prefix}Found {empty_mask.sum()} empty geometry/geometries at "
            f"index positions {bad[:10]}{'…' if len(bad) > 10 else ''}. "
            "Remove or repair before calling footprint_attributes functions."
        )

    # --- multi-part / non-polygon ---
    geom_types = geoms.geom_type
    bad_type_mask = ~geom_types.isin(["Polygon"])
    if bad_type_mask.any():
        bad = list(gdf.index[bad_type_mask])
        types_found = geom_types[bad_type_mask].unique().tolist()
        raise ValueError(
            f"{prefix}Found non-Polygon geometry types {types_found} at "
            f"index positions {bad[:10]}{'…' if len(bad) > 10 else ''}. "
            "Explode MultiPolygons with gdf.explode(index_parts=False) and "
            "keep only Polygon rows before calling footprint_attributes functions."
        )

    # --- invalid (self-intersecting etc.) ---
    invalid_mask = ~geoms.is_valid
    if invalid_mask.any():
        bad = list(gdf.index[invalid_mask])
        raise ValueError(
            f"{prefix}Found {invalid_mask.sum()} invalid (e.g. self-intersecting) "
            f"geometry/geometries at index positions {bad[:10]}{'…' if len(bad) > 10 else ''}. "
            "Repair with gdf.geometry = gdf.geometry.make_valid() before calling "
            "footprint_attributes functions."
        )


def cast(collection) -> np.ndarray:
    """Cast a geometry collection to a NumPy array of Shapely objects."""
    if Version(shapely.__version__) < Version("2"):
        raise ImportError("Shapely >= 2.0 is required.")
    if isinstance(collection, (gpd.GeoSeries, gpd.GeoDataFrame)):
        return np.asarray(collection.geometry.array)
    if isinstance(collection, (np.ndarray, list)):
        return np.asarray(collection)
    return np.array([collection])


def ring_inertia_z(polygon: shapely.Geometry) -> float:
    """Polar second moment of area (I_z) of a single ring about its centroid.

    Uses the shoelace formula for arbitrary polygons.

    Args:
        polygon: A simple (one-ring) Shapely Polygon.

    Returns:
        Unsigned I_z value (m⁴ when coordinates are in metres).
    """
    coords = shapely.get_coordinates(polygon)
    cx, cy = shapely.get_coordinates(shapely.centroid(polygon))[0]
    pts = coords - np.array([cx, cy])
    cross = pts[:-1, 0] * pts[1:, 1] - pts[1:, 0] * pts[:-1, 1]
    quad = (
        pts[1:, 0] ** 2
        + pts[1:, 0] * pts[:-1, 0]
        + pts[:-1, 0] ** 2
        + pts[1:, 1] ** 2
        + pts[1:, 1] * pts[:-1, 1]
        + pts[:-1, 1] ** 2
    )
    return float(np.abs(np.sum(cross * quad) / 12))


def calc_inertia_z(collection) -> np.ndarray:
    """Polar second moment of area (I_z) for each geometry in *collection*.

    Handles multi-part geometries and interior rings (holes) via the
    parallel-axis theorem, subtracting hole contributions.

    Args:
        collection: GeoDataFrame, GeoSeries, list, or array of Shapely Polygons.

    Returns:
        1-D NumPy array with one I_z per geometry (m⁴).
    """
    ga = cast(collection)
    parts, coll_ix = shapely.get_parts(ga, return_index=True)
    rings, ring_ix = shapely.get_rings(parts, return_index=True)
    coll_ix = np.repeat(coll_ix, shapely.get_num_interior_rings(parts) + 1)

    poly_rings = shapely.polygons(rings)
    is_ext = np.zeros_like(coll_ix, dtype=bool)
    is_ext[0] = True
    is_ext[1:] = ring_ix[1:] != ring_ix[:-1]

    df = gpd.GeoDataFrame(
        dict(coll_ix=coll_ix, ring_ix=ring_ix, is_ext=is_ext),
        geometry=poly_rings,
    )
    df["moa"] = df.geometry.apply(ring_inertia_z)
    df["sign"] = (1 - df["is_ext"].astype(int) * 2) * -1
    orig_centroids = shapely.centroid(ga)
    df["coll_centroid"] = orig_centroids[coll_ix]
    df["radius"] = shapely.distance(
        shapely.centroid(df.geometry.values), df["coll_centroid"].values
    )
    # groupby only produces a row for indices that appear in coll_ix.  Empty
    # geometries are skipped by get_parts, so their index never appears and
    # the raw .values would be shorter than the input.  Reindexing to the full
    # range guarantees one value per input geometry (0.0 for any empties).
    per_geom = (
        df.groupby("coll_ix")
        .apply(
            lambda g: np.sum(g["moa"] + g["sign"] * (g["radius"] ** 2)),
            include_groups=False,
        )
        .reindex(range(len(ga)), fill_value=0.0)
    )
    return per_geom.values


# ─────────────────────────────────────────────────────────────────────────────
# Principal inertia (2×2 tensor)
# ─────────────────────────────────────────────────────────────────────────────


def calc_principal_inertia(
    geoms: gpd.GeoSeries | gpd.GeoDataFrame,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute principal moments and eigenvectors of each footprint.

    Uses the exact closed-form polygon second-moment-of-area formulas
    (``Ixx = ∫y²dA``, ``Iyy = ∫x²dA``, ``Pxy = ∫xy dA``, via the shoelace
    identity -- the same family as :func:`ring_inertia_z`'s polar moment),
    evaluated about the polygon's own centroid. This is *exact* for any
    simple polygon, not an approximation: a previous version summed
    ``x²``/``y²``/``xy`` over the polygon's *vertices* (a discrete
    point-mass tensor of the vertex coordinates), which is not the
    continuous area integral the paper defines and gave visibly wrong
    principal axes even for a plain rectangle.

    Args:
        geoms: GeoSeries or GeoDataFrame of Polygon geometries.

    Returns:
        ``(I1, dir1, I2, dir2)`` where:
        - I1, I2: (N,) arrays of principal moments (I1 ≥ I2).
        - dir1, dir2: (N, 2) arrays of unit eigenvectors.
    """
    validate_geodataframe(geoms, context="calc_principal_inertia")
    if isinstance(geoms, gpd.GeoDataFrame):
        geoms = geoms.geometry
    n = len(geoms)
    I1_arr = np.zeros(n)
    I2_arr = np.zeros(n)
    dir1 = np.zeros((n, 2))
    dir2 = np.zeros((n, 2))

    for i, poly in enumerate(geoms):
        coords = shapely.get_coordinates(poly)
        cx, cy = shapely.get_coordinates(shapely.centroid(poly))[0]
        pts = coords - np.array([cx, cy])

        x0, y0 = pts[:-1, 0], pts[:-1, 1]
        x1, y1 = pts[1:, 0], pts[1:, 1]
        cross = x0 * y1 - x1 * y0
        if cross.sum() < 0:
            # Enforce a consistent (CCW) orientation so Ixx/Iyy come out
            # positive and Pxy has the correct sign, regardless of the
            # input ring's winding order.
            cross = -cross

        Ixx = np.sum(cross * (y0**2 + y0 * y1 + y1**2)) / 12.0
        Iyy = np.sum(cross * (x0**2 + x0 * x1 + x1**2)) / 12.0
        Pxy = np.sum(cross * (x0 * y1 + 2 * x0 * y0 + 2 * x1 * y1 + x1 * y0)) / 24.0
        Ixy = -Pxy  # matches the [[Ixx, Ixy], [Ixy, Iyy]] tensor convention below

        # Eigenvalues and eigenvectors
        evals, evecs = np.linalg.eigh(np.array([[Ixx, Ixy], [Ixy, Iyy]]))
        idx = np.argsort(evals)[::-1]  # Descending order
        I1_arr[i] = evals[idx[0]]
        I2_arr[i] = evals[idx[1]]
        v1 = evecs[:, idx[0]]
        v2 = evecs[:, idx[1]]
        # eigh's sign for each eigenvector is arbitrary (an axis has no
        # inherent direction), so without a fixed convention the same
        # physical axis can come back pointing to either corner from one
        # call to the next. Canonicalise by flipping so the largest-
        # magnitude component is positive (the same convention used by
        # e.g. scikit-learn's PCA sign-flip) -- this makes dir1/dir2 stable
        # and lets them be compared/plotted directly against direction.bbox.
        if v1[np.argmax(np.abs(v1))] < 0:
            v1 = -v1
        if v2[np.argmax(np.abs(v2))] < 0:
            v2 = -v2
        dir1[i] = v1
        dir2[i] = v2

    return I1_arr, dir1, I2_arr, dir2


# ─────────────────────────────────────────────────────────────────────────────
# Bounding box and rectangle projection
# ─────────────────────────────────────────────────────────────────────────────


def min_bounding_box(
    gdf: gpd.GeoDataFrame | gpd.GeoSeries,
) -> tuple[list, np.ndarray, list, np.ndarray]:
    """Minimum rotated bounding box for each footprint.

    Args:
        gdf: GeoDataFrame or GeoSeries of Polygon geometries (projected CRS).

    Returns:
        ``(L1, dir1, L2, dir2)`` where L1 ≥ L2 are side lengths and
        dir1, dir2 are unit vectors.
    """
    validate_geodataframe(gdf, context="min_bounding_box")
    if isinstance(gdf, gpd.GeoDataFrame):
        geoms = gdf.geometry
    else:
        geoms = gdf
    n = len(geoms)
    L1_list = []
    L2_list = []
    dir1 = np.zeros((n, 2))
    dir2 = np.zeros((n, 2))

    for i, poly in enumerate(geoms):
        # Compute rotated bounding box at different angles
        best_area = np.inf
        best_l1, best_l2 = 0, 0
        best_d1 = np.array([1, 0])
        best_d2 = np.array([0, 1])

        for angle in np.linspace(0, np.pi / 2, 90):
            cos_a, sin_a = np.cos(angle), np.sin(angle)
            rot_mat = np.array([[cos_a, sin_a], [-sin_a, cos_a]])
            coords = shapely.get_coordinates(poly)
            rot_coords = coords @ rot_mat.T
            min_x, min_y = rot_coords.min(axis=0)
            max_x, max_y = rot_coords.max(axis=0)
            l1, l2 = max_x - min_x, max_y - min_y
            area = l1 * l2
            if area < best_area:
                best_area = area
                best_l1, best_l2 = max(l1, l2), min(l1, l2)
                # Determine which axis is longer
                if l1 >= l2:
                    best_d1 = np.array([cos_a, sin_a])
                    best_d2 = np.array([-sin_a, cos_a])
                else:
                    best_d1 = np.array([-sin_a, cos_a])
                    best_d2 = np.array([cos_a, sin_a])

        L1_list.append(best_l1)
        L2_list.append(best_l2)
        dir1[i] = best_d1
        dir2[i] = best_d2

    return L1_list, dir1, L2_list, dir2


def circumscribed_rectangle_lengths(
    gdf: gpd.GeoDataFrame | gpd.GeoSeries,
    dir1: np.ndarray,
    dir2: np.ndarray,
) -> tuple[list, list]:
    """Project footprints onto given axes to get L1, L2.

    Args:
        gdf: GeoDataFrame or GeoSeries (projected CRS).
        dir1: (N, 2) unit vectors for L1 axis.
        dir2: (N, 2) unit vectors for L2 axis.

    Returns:
        ``(L1_list, L2_list)`` lists of lengths per building.
    """
    validate_geodataframe(gdf, context="circumscribed_rectangle_lengths")
    if isinstance(gdf, gpd.GeoDataFrame):
        geoms = gdf.geometry.values
    else:
        geoms = gdf.values
    L1_list = []
    L2_list = []

    for i, poly in enumerate(geoms):
        coords = shapely.get_coordinates(poly)
        proj1 = np.dot(coords, dir1[i])
        proj2 = np.dot(coords, dir2[i])
        L1 = proj1.max() - proj1.min()
        L2 = proj2.max() - proj2.min()
        L1_list.append(L1)
        L2_list.append(L2)

    return L1_list, L2_list


def inertia_side_lengths(
    I1: np.ndarray,
    I2: np.ndarray,
    area: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Average plan side lengths from principal moments (paper eqs. 2-3).

    ``slenderness = sqrt(I1 / I2)``; the side length perpendicular to the I1
    axis (i.e. running along the I2 eigenvector -- the geometrically longer
    side) is ``sqrt(area * slenderness)``; the other side is ``area`` divided
    by that length. For a true rectangle with sides ``a > b``,
    ``I1 = (1/12) b a^3``, ``I2 = (1/12) a b^3``, so this reduces exactly to
    ``L1 = a``, ``L2 = b``. For non-rectangular footprints this is an
    *average* side length, not a measured projection.

    Args:
        I1: (N,) larger principal moments.
        I2: (N,) smaller principal moments.
        area: (N,) footprint areas.

    Returns:
        ``(L1, L2)`` arrays, the longer and shorter average side lengths.
    """
    slenderness = np.sqrt(I1 / (I2 + 1e-30))
    L1 = np.sqrt(area * slenderness)
    L2 = area / (L1 + 1e-30)
    return L1, L2


def largest_convex_hull_gap_area(geoms: gpd.GeoSeries | gpd.GeoDataFrame) -> np.ndarray:
    """Area of the single largest setback for each footprint (paper compactness).

    Per the paper: compactness is found by taking the set difference between
    the (hole-filled) footprint and its convex hull, measuring the area of
    each resulting disconnected piece, and keeping the *largest* one -- not
    the sum of all of them.

    Args:
        geoms: GeoSeries or GeoDataFrame of Polygon geometries.

    Returns:
        (N,) array with the area of the largest hull-minus-footprint
        component per building (0.0 if the footprint is already convex).
    """
    validate_geodataframe(geoms, context="largest_convex_hull_gap_area")
    if isinstance(geoms, gpd.GeoDataFrame):
        geoms = geoms.geometry

    filled = fill_holes(geoms)
    hull = filled.convex_hull
    gaps = hull.difference(filled)

    def _largest_part_area(g):
        if g is None or g.is_empty:
            return 0.0
        if hasattr(g, "geoms"):
            return max((p.area for p in g.geoms), default=0.0)
        return g.area

    return np.array([_largest_part_area(g) for g in gaps])


# ─────────────────────────────────────────────────────────────────────────────
# Setback and hole metrics
# ─────────────────────────────────────────────────────────────────────────────


def max_hole_area_ratio(
    gdf: gpd.GeoDataFrame,
    min_area_fraction: float = 0.001,
) -> list[float]:
    """Ratio of largest hole to filled (solid) area for each footprint.

    Args:
        gdf: GeoDataFrame of building footprints.
        min_area_fraction: Minimum fractional area to consider a hole.

    Returns:
        List of hole-area ratios (0–1).
    """
    validate_geodataframe(gdf, context="max_hole_area_ratio")
    exterior = fill_holes(gdf.geometry)
    ratios = []
    for i, poly in enumerate(gdf.geometry):
        filled_area = exterior.iloc[i].area
        if filled_area < 1e-12:
            ratios.append(0.0)
            continue
        hole_area = filled_area - poly.area
        if filled_area * min_area_fraction > hole_area:
            ratios.append(0.0)
        else:
            ratios.append(hole_area / filled_area)
    return ratios


def setback_gndt_metrics(
    gdf: gpd.GeoDataFrame,
    L1: np.ndarray,
    dir1: np.ndarray,
    L2: np.ndarray,
    dir2: np.ndarray,
    min_length: float = 0.0,
    min_area_fraction: float = 0.001,
) -> tuple[list[float], list[float], list[float]]:
    """GNDT-style dual-configuration setback metrics per building.

    For every disconnected setback piece (``hull.difference(footprint)``),
    the piece's own circumscribed-rectangle extents ``b1`` (along *dir1*)
    and ``b2`` (along *dir2*) are computed. Per the paper: *"a and b must be
    measured along the same direction [as L]... this results in two
    possible configurations: either (L1, a1, b1) or (L2, a2, b2). For
    beta2 = b/L, use the minimum between b1/L1 and b2/L2."* This function
    picks, across all setback pieces, whichever single piece/configuration
    minimises ``b/L``, and additionally measures ``c`` -- the true
    protrusion width of the *solid* footprint (not just the setback piece)
    measured perpendicular to that winning ``b``, via a ray cast through
    the piece's centroid.

    Args:
        gdf: GeoDataFrame of building footprints (projected CRS).
        L1: (N,) footprint extents along *dir1* (the L1/longer axis).
        dir1: (N, 2) unit vectors for the L1 axis.
        L2: (N,) footprint extents along *dir2* (the L2/shorter axis).
        dir2: (N, 2) unit vectors for the L2 axis.
        min_length: Minimum setback side length (m) to be considered.
        min_area_fraction: Minimum fractional setback area to be considered.

    Returns:
        ``(ratio, b, c)`` lists, one entry per building:
        - ``ratio``: ``min(b1/L1, b2/L2)`` for the dominant setback (the
          value to use directly for beta2/setback_ratio); 0 for convex
          footprints.
        - ``b``: the winning configuration's own b1 or b2 (m); 0 if convex.
        - ``c``: protrusion width of the solid footprint perpendicular to
          ``b`` (m), for beta6 = c/b; 0 if convex.
    """
    validate_geodataframe(gdf, context="setback_gndt_metrics")
    gdf = ensure_projected(gdf).reset_index(drop=True)
    exterior = fill_holes(gdf.geometry)

    sb_geom = exterior.convex_hull.difference(exterior)
    sb = gpd.GeoDataFrame(
        {
            "orig_id": np.arange(len(gdf)),
            "footprint_no_holes": exterior.values,
            "L1": np.asarray(L1, dtype=float),
            "L2": np.asarray(L2, dtype=float),
            "d1x": dir1[:, 0],
            "d1y": dir1[:, 1],
            "d2x": dir2[:, 0],
            "d2y": dir2[:, 1],
        },
        geometry=sb_geom.values,
        crs=gdf.crs,
    ).explode("geometry", ignore_index=True)

    sb["compact"] = (
        1 - sb["footprint_no_holes"].area / sb["footprint_no_holes"].convex_hull.area
    )
    sb.loc[sb["compact"] <= min_area_fraction, "geometry"] = Polygon()

    non_empty = sb[~sb.geometry.is_empty]
    if len(non_empty) == 0:
        zeros = [0.0] * len(gdf)
        return zeros, list(zeros), list(zeros)

    _d1 = np.column_stack([non_empty["d1x"].values, non_empty["d1y"].values])
    _d2 = np.column_stack([non_empty["d2x"].values, non_empty["d2y"].values])
    b1_list, b2_list = circumscribed_rectangle_lengths(non_empty, _d1, _d2)

    sb["b1"] = 0.0
    sb["b2"] = 0.0
    sb.loc[~sb.geometry.is_empty, "b1"] = b1_list
    sb.loc[~sb.geometry.is_empty, "b2"] = b2_list
    sb.loc[(sb["b1"] < min_length) | (sb["b2"] < min_length), ["b1", "b2"]] = 0.0

    sb["ratio1"] = sb["b1"] / (sb["L1"] + 1e-12)
    sb["ratio2"] = sb["b2"] / (sb["L2"] + 1e-12)
    sb["ratio"] = np.minimum(sb["ratio1"], sb["ratio2"])
    idx = sb.groupby("orig_id")["ratio"].idxmax()
    dominant = sb.loc[idx].copy().reset_index(drop=True)

    # For each dominant setback, measure c perpendicular to the winning b.
    # The winning config is whichever of ratio1/ratio2 is the (smaller)
    # binding one, NOT whichever of b1/b2 is numerically larger (b1/b2 alone
    # ignores that L1 != L2 in general).
    b_values = []
    c_values = []
    for _, row in dominant.iterrows():
        if row.geometry.is_empty:
            b_values.append(0.0)
            c_values.append(0.0)
            continue
        if row["ratio1"] <= row["ratio2"]:
            b = row["b1"]
            perp = np.array([row["d2x"], row["d2y"]])  # perpendicular to b1 -> dir2
        else:
            b = row["b2"]
            perp = np.array([row["d1x"], row["d1y"]])  # perpendicular to b2 -> dir1
        b_values.append(float(b))

        # Cast a line through the setback centroid in perp direction
        cx, cy = row.geometry.centroid.x, row.geometry.centroid.y
        diag = np.hypot(
            row["footprint_no_holes"].bounds[2] - row["footprint_no_holes"].bounds[0],
            row["footprint_no_holes"].bounds[3] - row["footprint_no_holes"].bounds[1],
        )
        line = LineString(
            [
                (cx - perp[0] * (diag + 1), cy - perp[1] * (diag + 1)),
                (cx + perp[0] * (diag + 1), cy + perp[1] * (diag + 1)),
            ]
        )
        intersec = row["footprint_no_holes"].intersection(line)
        # Find the segment closest to the setback centroid
        parts = list(intersec.geoms) if hasattr(intersec, "geoms") else [intersec]
        pt_c = Point(cx, cy)
        closest = min(parts, key=lambda g: g.centroid.distance(pt_c), default=None)
        c_values.append(float(closest.length) if closest is not None else 0.0)

    dominant["b"] = b_values
    dominant["c"] = c_values
    out = (
        pd.DataFrame({"orig_id": np.arange(len(gdf))})
        .merge(dominant[["orig_id", "ratio", "b", "c"]], on="orig_id", how="left")
        .fillna(0.0)
    )
    return list(out["ratio"]), list(out["b"]), list(out["c"])


# ─────────────────────────────────────────────────────────────────────────────
# GNDT 'a' — inscribed-circle main-element construction
# ─────────────────────────────────────────────────────────────────────────────


def max_inscribed_circle(polygon, grid_n: int = 15) -> tuple[float, float, float]:
    """Largest circle fitting inside a single polygon (coarse grid + refine).

    Args:
        polygon: A single Shapely Polygon.
        grid_n: Grid resolution per axis for the initial coarse search.

    Returns:
        ``(cx, cy, r)`` -- centre and radius of the largest inscribed circle.
        ``r=0.0`` (centred on the centroid) for a degenerate (zero-area)
        polygon.
    """
    from scipy.optimize import minimize

    minx, miny, maxx, maxy = polygon.bounds
    boundary = polygon.boundary
    xs = np.linspace(minx, maxx, grid_n)
    ys = np.linspace(miny, maxy, grid_n)
    best, best_d = None, -1.0
    for x in xs:
        for y in ys:
            p = Point(x, y)
            if polygon.contains(p):
                d = p.distance(boundary)
                if d > best_d:
                    best_d, best = d, (x, y)

    if best is None:
        c = polygon.centroid
        return float(c.x), float(c.y), 0.0

    def _neg_dist(pt):
        p = Point(pt)
        if not polygon.contains(p):
            return 1e6
        return -p.distance(boundary)

    res = minimize(
        _neg_dist, best, method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-6}
    )
    cx, cy = res.x
    r = -res.fun
    return float(cx), float(cy), float(r)


def circle_tangent_points(
    polygon,
    cx: float,
    cy: float,
    r: float,
    n: int = 1440,
    tol_frac: float = 0.01,
) -> list[tuple[float, float]]:
    """Points where the inscribed circle touches the polygon boundary.

    Samples the circle at *n* angles, keeps samples within ``tol_frac * r``
    of the boundary, and collapses each contiguous run of close samples to
    one representative point (a true tangency is a single point; a flush
    contact along a straight edge segment collapses to its midpoint).

    Args:
        polygon: The Shapely Polygon the circle is inscribed in.
        cx, cy, r: Centre and radius of the inscribed circle.
        n: Number of angular samples.
        tol_frac: Distance-to-boundary tolerance, as a fraction of *r*.

    Returns:
        List of ``(x, y)`` tangent points.
    """
    if r < 1e-9:
        return []
    boundary = polygon.boundary
    tol = max(r * tol_frac, 1e-9)
    thetas = np.linspace(0, 2 * np.pi, n, endpoint=False)
    pts = [(cx + r * np.cos(t), cy + r * np.sin(t)) for t in thetas]
    dists = np.array([Point(p).distance(boundary) for p in pts])
    close = dists < tol

    reps: list[tuple[float, float]] = []
    visited = np.zeros(len(close), dtype=bool)
    n_pts = len(close)
    for i in range(n_pts):
        if close[i] and not visited[i]:
            j, group = i, []
            while close[j % n_pts] and not visited[j % n_pts]:
                visited[j % n_pts] = True
                group.append(j % n_pts)
                j += 1
                if j % n_pts == i:
                    break
            reps.append(pts[group[len(group) // 2]])
    return reps


def main_element_a_lengths(
    polygon,
    dir1: np.ndarray,
    dir2: np.ndarray,
    grid_n: int = 15,
) -> tuple[float, float, tuple[float, float]]:
    """GNDT 'a' lengths for both (L, a) configurations of a single footprint.

    Implements the paper's 3-step process (fig. 7): (1) inscribe the
    largest possible circle; (2) find its tangent points with the footprint
    boundary; (3) measure the extent of those tangent points along the
    footprint's own principal axes -- this is the circumscribed rectangle
    "with sides parallel to the directions of the minimum bounding box of
    the footprint" the paper describes.

    A plain rectangle's inscribed circle only touches 2 opposite sides
    (collinear tangent points, degenerate along the other axis), so with
    <=2 tangent points both configurations collapse to the circle diameter
    directly rather than projecting (a numerically safer choice than
    projecting near-collinear points, which -- especially for the inertia
    method, whose axes are not always exactly axis-aligned -- can pick up a
    small spurious nonzero extent instead of a clean zero).

    Args:
        polygon: A single Shapely Polygon.
        dir1: Unit vector for the L1 axis.
        dir2: Unit vector for the L2 axis.
        grid_n: Grid resolution passed to :func:`max_inscribed_circle`.

    Returns:
        ``(a1, a2, center)``: a1 pairs with L1 (measured along dir2,
        perpendicular to dir1); a2 pairs with L2 (measured along dir1);
        ``center`` is the world-coordinate centre of the (a1 x a2)
        rectangle -- the midpoint of the tangent points' projected bounds,
        which is generally NOT the inscribed circle's own centre (that's
        only true when the tangent points happen to be symmetric about it,
        e.g. for a plain rectangle).
    """
    cx, cy, r = max_inscribed_circle(polygon, grid_n=grid_n)
    touch_pts = circle_tangent_points(polygon, cx, cy, r)
    if len(touch_pts) <= 2:
        return 2 * r, 2 * r, (cx, cy)
    P = np.array(touch_pts)
    proj1 = P @ dir1
    proj2 = P @ dir2
    a2 = float(proj1.max() - proj1.min())
    a1 = float(proj2.max() - proj2.min())
    mid1 = (proj1.max() + proj1.min()) / 2
    mid2 = (proj2.max() + proj2.min()) / 2
    center = (mid1 * dir1[0] + mid2 * dir2[0], mid1 * dir1[1] + mid2 * dir2[1])
    return a1, a2, center


def main_element_a_lengths_batch(
    gdf: gpd.GeoDataFrame,
    dir1: np.ndarray,
    dir2: np.ndarray,
    grid_n: int = 15,
) -> tuple[list[float], list[float], list[tuple[float, float]]]:
    """Batched :func:`main_element_a_lengths` over a GeoDataFrame.

    Args:
        gdf: GeoDataFrame of building footprints (projected CRS).
        dir1: (N, 2) unit vectors for the L1 axis.
        dir2: (N, 2) unit vectors for the L2 axis.
        grid_n: Grid resolution passed to :func:`max_inscribed_circle`.

    Returns:
        ``(a1_list, a2_list, center_list)``.
    """
    validate_geodataframe(gdf, context="main_element_a_lengths_batch")
    a1_list, a2_list, center_list = [], [], []
    filled = fill_holes(gdf.geometry)
    for i, poly in enumerate(filled):
        a1, a2, center = main_element_a_lengths(poly, dir1[i], dir2[i], grid_n=grid_n)
        a1_list.append(a1)
        a2_list.append(a2)
        center_list.append(center)
    return a1_list, a2_list, center_list


def setback_pieces(
    polygon,
    dir1: np.ndarray,
    dir2: np.ndarray,
    min_area_fraction: float = 0.001,
) -> list[tuple[float, float, object]]:
    """All disconnected setback pieces of a single footprint.

    A setback is a connected component of ``hull.difference(footprint)``
    (paper fig. 8a). Each piece's own circumscribed-rectangle extents along
    *dir1* and *dir2* are computed -- these are the GNDT-style "sides of the
    circumscribed rectangle" (fig. 8b) for that individual piece.

    This is the per-piece building block: :func:`setback_gndt_metrics` uses
    it internally (per building, vectorised) to pick the dominant
    configuration for beta2/beta6; this function exposes every piece, e.g.
    for visualising all of a T- or X-shape's setbacks individually.

    Args:
        polygon: A single Shapely Polygon.
        dir1: Unit vector for the L1 axis.
        dir2: Unit vector for the L2 axis.
        min_area_fraction: Minimum fractional area (of the convex hull) for
            a piece to be considered.

    Returns:
        List of ``(ext1, ext2, piece)`` tuples, sorted by piece area
        descending: ``ext1``/``ext2`` are the piece's own extents along
        *dir1*/*dir2*, and ``piece`` is its Shapely (Multi)Polygon geometry.
        Empty list for a convex footprint.
    """
    filled = Polygon(polygon.exterior)
    hull = filled.convex_hull
    if hull.area < 1e-12:
        return []
    gap = hull.difference(filled)
    parts = (
        list(gap.geoms)
        if hasattr(gap, "geoms")
        else ([gap] if not gap.is_empty else [])
    )

    results = []
    for part in parts:
        if part.area < hull.area * min_area_fraction:
            continue
        coords = shapely.get_coordinates(part)
        proj1 = coords @ dir1
        proj2 = coords @ dir2
        ext1 = float(proj1.max() - proj1.min())
        ext2 = float(proj2.max() - proj2.min())
        results.append((ext1, ext2, part))
    results.sort(key=lambda t: -t[2].area)
    return results


# ─────────────────────────────────────────────────────────────────────────────
# NTC-23 hole ratio — hole's own minimum bounding box
# ─────────────────────────────────────────────────────────────────────────────


def hole_h_over_l(
    gdf: gpd.GeoDataFrame, min_area_fraction: float = 0.001
) -> list[float]:
    """NTC-23 hole ratio ``h / L`` per footprint (worst hole).

    Per the paper: for each hole, a bounding rectangle is drawn aligned with
    the hole's OWN minimum bounding box (independent of, and generally not
    parallel to, the building's own principal axes); ``h`` is the shorter
    side of that box. ``L`` is the length of the building's cross-section
    through the hole's centroid, measured along the hole's longer axis.
    This is a distinct construction from ASCE7's area-based hole_ratio
    (``A_hole / A_filled``) -- the two must not be confused.

    Args:
        gdf: GeoDataFrame of building footprints (projected CRS).
        min_area_fraction: Minimum hole area (relative to the filled
            footprint) to be considered.

    Returns:
        List of h/L ratios (0.0 for footprints with no significant holes).
    """
    validate_geodataframe(gdf, context="hole_h_over_l")
    ratios = []
    for poly in gdf.geometry:
        filled = Polygon(poly.exterior)
        filled_area = filled.area
        best_ratio = 0.0
        for ring in poly.interiors:
            hole_poly = Polygon(ring)
            if filled_area < 1e-12 or hole_poly.area < filled_area * min_area_fraction:
                continue
            hole_gs = gpd.GeoSeries([hole_poly], crs=gdf.crs)
            L1h, dir1h, L2h, _ = min_bounding_box(hole_gs)
            h = L2h[0]
            long_dir = dir1h[0]

            cx, cy = hole_poly.centroid.coords[0]
            diag = np.hypot(
                filled.bounds[2] - filled.bounds[0], filled.bounds[3] - filled.bounds[1]
            )
            line = LineString(
                [
                    (cx - long_dir[0] * (diag + 1), cy - long_dir[1] * (diag + 1)),
                    (cx + long_dir[0] * (diag + 1), cy + long_dir[1] * (diag + 1)),
                ]
            )
            intersec = filled.intersection(line)
            parts = list(intersec.geoms) if hasattr(intersec, "geoms") else [intersec]
            pt_c = Point(cx, cy)
            closest = min(parts, key=lambda g: g.centroid.distance(pt_c), default=None)
            L = closest.length if closest is not None else 0.0

            ratio = h / (L + 1e-12) if L > 1e-12 else 0.0
            best_ratio = max(best_ratio, ratio)
        ratios.append(best_ratio)
    return ratios


# ─────────────────────────────────────────────────────────────────────────────
# Contact-edge helpers (used by position module)
# ─────────────────────────────────────────────────────────────────────────────


def split_linestring_to_segments(ls) -> shapely.MultiLineString:
    """Split any linear geometry into individual two-point segments.

    Accepts LineString, LinearRing, or MultiLineString.  Each output segment
    has exactly 2 coordinate pairs so that :func:`edge_normal` can always
    compute a well-defined normal.

    Args:
        ls: A Shapely linear geometry.

    Returns:
        MultiLineString whose parts are all two-point LineStrings.
    """
    segs = []
    geom_type = shapely.get_type_id(ls)
    # MultiLineString (type id 5) → recurse on parts
    if geom_type == 5:
        for part in shapely.get_parts(ls):
            segs.extend(_two_point_segments(part))
    else:
        segs.extend(_two_point_segments(ls))
    if not segs:
        return shapely.MultiLineString()
    return shapely.MultiLineString(segs)


def _two_point_segments(ls) -> list:
    """Return a list of two-point (coord[i], coord[i+1]) LineStrings."""
    coords = shapely.get_coordinates(ls)
    return [
        LineString([coords[i], coords[i + 1]])
        for i in range(len(coords) - 1)
        if np.linalg.norm(coords[i] - coords[i + 1]) > 1e-9  # skip zero-length
    ]


def select_touching_edges(gdf: gpd.GeoDataFrame, buffer: float = 0) -> gpd.GeoDataFrame:
    """Keep only boundary segments of each footprint that touch neighbours.

    Args:
        gdf: GeoDataFrame of building footprints (projected CRS).
        buffer: Contact detection buffer (metres).

    Returns:
        Copy of *gdf* with geometry replaced by touching boundary segments.
    """
    out = gdf.copy()
    buf = max(buffer, 0.0001)
    union = shapely.buffer(
        out.geometry.union_all(), buf, cap_style="square", join_style="mitre"
    )
    union = shapely.buffer(
        union, -buffer - 0.001, cap_style="square", join_style="mitre"
    )
    out.geometry = (
        out.geometry.buffer(max(buffer, 0.001), cap_style="square", join_style="mitre")
        .buffer(min(-buffer, -0.001), cap_style="square", join_style="mitre")
        .boundary.intersection(union)
    )
    return out


def explode_edges(gdf: gpd.GeoDataFrame, min_length: float = 0.0) -> gpd.GeoDataFrame:
    """Explode boundary geometries into individual two-point segments.

    Args:
        gdf: GeoDataFrame with LineString/MultiLineString geometry.
        min_length: Minimum segment length to retain (metres).

    Returns:
        Exploded GeoDataFrame with an ``edges`` geometry column. The
        original (pre-split) geometry column is dropped: leaving it in
        place is a footgun for any later ``.apply(..., axis=1)`` call,
        since a plain row ``Series`` resolves ``r.geometry`` by column
        label, not by the GeoDataFrame's active-geometry column -- it
        would silently return the stale, un-split geometry instead of the
        two-point ``edges`` segment.
    """
    out = gdf.copy()
    crs = gdf.crs
    geom_col = out.geometry.name
    out = out[~out.geometry.is_empty].explode(index_parts=False).reset_index(drop=True)
    out["edges"] = out.geometry.apply(split_linestring_to_segments)
    out = (
        out.drop(columns=geom_col)
        .set_geometry("edges", crs=crs)
        .explode()
        .reset_index(drop=True)
    )
    out = out[out["edges"].length > max(min_length, 0.001)]
    return out


def centre_of_mass_and_stiffness(
    geoms: gpd.GeoSeries | gpd.GeoDataFrame,
    wall_height: float = 3.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Centre of mass (CM) and centre of stiffness (CS) under the hollow-box
    building assumption (uniform walls + a slab, see paper eq. 1).

    The building is modelled as a hollow box: a horizontal slab (ceiling),
    whose mass is distributed like the footprint polygon's own area, plus
    perimeter walls of height *wall_height*, whose mass is distributed along
    the polygon boundary (exterior + interior/hole rings). CM is the
    area/perimeter-weighted average of the two centroids; CS -- the centre of
    the lateral-force-resisting elements (the walls) -- is the boundary
    centroid alone.

    Args:
        geoms: GeoSeries or GeoDataFrame of Polygon geometries (projected CRS).
        wall_height: Assumed storey height (m). Cancels out of the CM formula
            for additional identical storeys, but is kept as a parameter
            since it is the physical weight of the wall mass relative to the
            slab mass for a single storey.

    Returns:
        ``(cm, cs)`` each an (N, 2) array of (x, y) coordinates.
    """
    validate_geodataframe(geoms, context="centre_of_mass_and_stiffness")
    if isinstance(geoms, gpd.GeoDataFrame):
        geoms = geoms.geometry

    area = geoms.area.values
    boundary = geoms.boundary
    perimeter = boundary.length.values

    slab_centroid = np.column_stack([geoms.centroid.x.values, geoms.centroid.y.values])
    wall_centroid = np.column_stack(
        [boundary.centroid.x.values, boundary.centroid.y.values]
    )

    wall_weight = perimeter * wall_height
    total_weight = area + wall_weight
    # Degenerate (zero-area / zero-perimeter) geometries fall back to the slab centroid.
    safe_total = np.where(total_weight > 1e-12, total_weight, 1.0)

    cm = (
        slab_centroid * area[:, None] + wall_centroid * wall_weight[:, None]
    ) / safe_total[:, None]
    cs = wall_centroid

    return cm, cs


def eq_circle_inertia(area: np.ndarray | float) -> np.ndarray:
    """Polar moment of inertia of a circle with the given area.

    For a circle of radius r (area = π r²), I_z = π r⁴ / 2 = area² / (2π).

    Args:
        area: Footprint area (m²).  Scalar or array.

    Returns:
        I_z of the equivalent circle (m⁴).
    """
    area = np.asarray(area, dtype=float)
    return area**2 / (2 * np.pi)


def explode_exterior_rings(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Explode each polygon into its exterior ring LineString only (no holes).

    Used by :func:`convex_hull_irregularity` to iterate over boundary edges
    while ignoring interior courtyards.

    Args:
        gdf: GeoDataFrame of Polygon geometries.

    Returns:
        GeoDataFrame with geometry replaced by the exterior LineString,
        one row per input polygon (no exploding of multi-parts).
    """
    out = gdf.copy()
    crs = gdf.crs
    out.geometry = out.geometry.apply(
        lambda g: LineString(g.exterior.coords) if hasattr(g, "exterior") else g
    )
    out.crs = crs
    return out


def explode_rings(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Explode each polygon into its exterior + interior ring LineStrings.

    Args:
        gdf: GeoDataFrame of Polygon geometries.

    Returns:
        GeoDataFrame with one row per ring.
    """
    out = gdf.copy()
    crs = gdf.crs
    out["_ext"] = out.exterior
    out["_int"] = out.interiors
    out.geometry = out.apply(
        lambda r: shapely.MultiLineString(list(r["_int"]) + [r["_ext"]]), axis=1
    )
    out.crs = crs
    return out.drop(columns=["_ext", "_int"]).explode().reset_index(drop=True)


def edge_normal(segment: LineString, scale: float = 1.0) -> tuple[Point, np.ndarray]:
    """Unit normal of a segment, optionally scaled by contact area.

    Args:
        segment: Shapely LineString.  Only the first and last coordinate are
            used, so multi-point LineStrings are handled correctly.
        scale: If > 0, scale the normal by ``scale × length`` (i.e., scale is
            the building height, giving force ∝ wall area).
            If 0, return a pure unit normal.

    Returns:
        ``(midpoint, force_vector)`` where ``force_vector`` has magnitude
        ``scale × length`` (or 1 when scale=0).
    """
    midpoint = segment.interpolate(0.5, normalized=True)
    coords = shapely.get_coordinates(segment)
    p1 = coords[0]
    # Use whichever vertex is farthest from p1 to define the tangent: for a
    # simple 2-point segment this is just coords[-1], but for a multi-point
    # segment whose first and last points happen to coincide (e.g. an
    # out-and-back or looped contact boundary in messy real-world data),
    # coords[-1] - p1 would be the zero vector even though the segment has
    # non-zero length, giving a NaN normal after normalisation.
    dists = np.linalg.norm(coords - p1, axis=1)
    far_idx = int(np.argmax(dists))
    if dists[far_idx] < 1e-12:
        # Genuinely a single point: no well-defined direction, no force.
        return midpoint, np.zeros(2)
    tangent = coords[far_idx] - p1
    normal = np.array([-tangent[1], tangent[0]])
    normal /= np.linalg.norm(normal)
    if scale == 0:
        return midpoint, normal
    return midpoint, normal * (scale * segment.length)


def edge_momentum(
    midpoint: Point,
    force: np.ndarray,
    centroid: Point,
    min_dist: float = 0.0,
) -> np.ndarray:
    """Moment of a contact force about the footprint centroid.

    Returns a 4-element array encoding whether the moment helps or hinders
    torsion.  Elements: [total, only_if_pos, only_if_neg_from_close, only_if_pos_from_close].

    Args:
        midpoint: Application point of the force.
        force: 2-D force vector.
        centroid: Reference point (footprint centroid).
        min_dist: If the edge is closer than this to the centroid, the
            momentum is only counted when it reduces the net torque.

    Returns:
        4-element NumPy float array.
    """
    r = np.array([midpoint.x - centroid.x, midpoint.y - centroid.y])
    M = float(r[0] * force[1] - r[1] * force[0])  # 2-D cross product = torque
    dist = float(np.linalg.norm(r))
    if dist < min_dist:
        if M > 0:
            return np.array([M, 0.0, M, 0.0])
        else:
            return np.array([M, 0.0, 0.0, M])
    return np.array([M, M, M, M])
