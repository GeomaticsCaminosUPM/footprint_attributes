"""
Building direction computation.

A building's *direction* is the orientation of its plan relative to geographic
North.  Two principal methods are provided, both returning the same four
quantities (L1, dir1, L2, dir2):

- **inertia** : principal axes of the second moment of area of the footprint.
  L1 corresponds to the *larger* principal moment I1 (the axis with more mass
  spread out), which is the geometrically *longer* dimension.
  The *weak* axis (L2 / dir2) is the direction with less resistance to lateral
  forces and is used as the building's bearing.

- **bbox** : axes of the minimum rotated bounding box (MBB).
  L1 is always the *longer* side; L2 is the shorter side.

Convention
----------
- **L1** (float, m) : longer plan dimension.
- **L2** (float, m) : shorter plan dimension.  Always L2 ≤ L1.
- **dir1** : unit vector parallel to L1.
- **dir2** : unit vector parallel to L2 (perpendicular to dir1).
- **bearing** : clockwise angle from geographic North (+y in UTM) to dir2,
  restricted to [-90°, 90°].  This is the direction of the *weak* axis.

When a ``direction`` keyword is supplied (a 2-element vector), it is used as
dir1 and the computation of L1/L2 is forced along that axis.  This is needed
for some shape functions that must project dimensions onto a prescribed axis.

Usage
-----
>>> import geopandas as gpd
>>> from footprint_attributes import direction
>>>
>>> footprints = gpd.read_file("footprints.gpkg")
>>>
>>> # Default: inertia method, return bearing only
>>> bearing = direction.inertia(footprints)
>>>
>>> # Bounding-box method, all outputs
>>> L1, dir1, L2, dir2, bearing = direction.bbox(footprints, mode="all")
>>>
>>> # Force a specific dir1
>>> L1, dir1, L2, dir2, bearing = direction.inertia(
...     footprints, mode="all", direction=np.array([1, 0])
... )
"""

from __future__ import annotations

import numpy as np
import geopandas as gpd

from .geometry import (
    ensure_projected,
    to_gdf,
    validate_geodataframe,
    calc_principal_inertia,
    min_bounding_box,
    circumscribed_rectangle_lengths,
    inertia_side_lengths,
    bearing_from_dir,
    normalize_vec,
)


# Valid mode strings
_MODES = {"bearing", "dimensions", "directions", "all"}


# ─────────────────────────────────────────────────────────────────────────────
# Internal: forced-direction projection
# ─────────────────────────────────────────────────────────────────────────────


def _apply_forced_direction(
    gdf: gpd.GeoDataFrame,
    forced_dir1: np.ndarray,
) -> tuple[list, np.ndarray, list, np.ndarray]:
    """Project footprints onto a prescribed dir1 axis to get L1, L2.

    Args:
        gdf: Projected GeoDataFrame of footprints.
        forced_dir1: 2-D direction vector (will be normalised to unit vector).

    Returns:
        ``(L1, dir1, L2, dir2)`` where dir1 == forced_dir1_unit and
        dir2 is perpendicular.
    """
    n = len(gdf)
    dir1_unit = normalize_vec(np.asarray(forced_dir1, dtype=float))
    dir2_unit = np.array([-dir1_unit[1], dir1_unit[0]])  # perpendicular

    dir1_arr = np.tile(dir1_unit, (n, 1))
    dir2_arr = np.tile(dir2_unit, (n, 1))

    L1_list, L2_list = circumscribed_rectangle_lengths(gdf, dir1_arr, dir2_arr)
    return L1_list, dir1_arr, L2_list, dir2_arr


# ─────────────────────────────────────────────────────────────────────────────
# Internal: pack outputs according to mode
# ─────────────────────────────────────────────────────────────────────────────


def _pack(
    L1: list | np.ndarray,
    dir1: np.ndarray,
    L2: list | np.ndarray,
    dir2: np.ndarray,
    mode: str,
) -> (
    list[float]
    | tuple[list, list]
    | tuple[np.ndarray, np.ndarray]
    | tuple[list, np.ndarray, list, np.ndarray, list[float]]
):
    """Select and return only the quantities requested by *mode*.

    Args:
        L1:    List/array of longer dimensions.
        dir1:  (N, 2) unit vectors along L1.
        L2:    List/array of shorter dimensions.
        dir2:  (N, 2) unit vectors along L2.
        mode:  One of "bearing", "dimensions", "directions", "all".

    Returns:
        - "bearing"    → list of floats (degrees)
        - "dimensions" → (L1_list, L2_list)
        - "directions" → (dir1_array, dir2_array)
        - "all"        → (L1_list, dir1_array, L2_list, dir2_array, bearing_list)
    """
    if mode not in _MODES:
        raise ValueError(f"mode must be one of {_MODES}, got {mode!r}")

    bearings = [bearing_from_dir(dir2[i]) for i in range(len(dir2))]
    L1 = list(L1)
    L2 = list(L2)

    if mode == "bearing":
        return bearings
    if mode == "dimensions":
        return L1, L2
    if mode == "directions":
        return dir1, dir2
    # mode == "all"
    return L1, dir1, L2, dir2, bearings


# ─────────────────────────────────────────────────────────────────────────────
# Public methods
# ─────────────────────────────────────────────────────────────────────────────


def inertia(
    footprints_gdf: gpd.GeoDataFrame | gpd.GeoSeries,
    mode: str = "bearing",
    direction: np.ndarray | None = None,
) -> (
    list[float]
    | tuple[list, list]
    | tuple[np.ndarray, np.ndarray]
    | tuple[list, np.ndarray, list, np.ndarray, list[float]]
):
    """Building direction from the principal axes of inertia.

    The *larger* principal moment I1 is associated with dir1 and L1 (the
    direction in which mass is most spread out = geometrically longer).
    The *smaller* principal moment I2 is associated with dir2 and L2, which
    is the *weak* lateral direction and defines the bearing.

    Args:
        footprints_gdf: GeoDataFrame or GeoSeries of Polygon footprints.
        mode: What to return. One of:
            - ``"bearing"`` (default): list of bearings (degrees from N).
            - ``"dimensions"``: ``(L1_list, L2_list)``.
            - ``"directions"``: ``(dir1_array, dir2_array)``.
            - ``"all"``: ``(L1_list, dir1_array, L2_list, dir2_array, bearing_list)``.
        direction: Optional 2-D vector forcing dir1.  When given, L1 and L2
            are projected onto that axis and its perpendicular.

    Returns:
        As described under *mode*.
    """
    gdf = to_gdf(footprints_gdf)
    gdf = ensure_projected(gdf)
    validate_geodataframe(gdf, context="direction.inertia")

    if direction is not None:
        L1, dir1, L2, dir2 = _apply_forced_direction(gdf, direction)
    else:
        # calc_principal_inertia's eigenvector for the *larger* eigenvalue
        # (I1) is the axis of greater bending stiffness, which points along
        # the *shorter* physical plan dimension (an area-moment-of-inertia
        # tensor is the covariance matrix with x/y roles swapped). So the
        # smaller-eigenvalue eigenvector (I2, dir2) is the one that runs
        # along the geometrically longer side -> maps to dir1 here.
        I1_arr, eig_dir1, I2_arr, eig_dir2 = calc_principal_inertia(gdf.geometry)
        area = gdf.geometry.area.values
        # Paper eqs. 2-3: L1/L2 are the average side lengths implied by the
        # slenderness sqrt(I1/I2), not the polygon's measured projection onto
        # the eigenvectors -- this matches a true rectangle exactly and is
        # more stable than the raw vertex-based projection for other shapes.
        L1, L2 = inertia_side_lengths(I1_arr, I2_arr, area)
        dir1, dir2 = eig_dir2, eig_dir1

    return _pack(L1, dir1, L2, dir2, mode)


def bbox(
    footprints_gdf: gpd.GeoDataFrame | gpd.GeoSeries,
    mode: str = "bearing",
    direction: np.ndarray | None = None,
) -> (
    list[float]
    | tuple[list, list]
    | tuple[np.ndarray, np.ndarray]
    | tuple[list, np.ndarray, list, np.ndarray, list[float]]
):
    """Building direction from the minimum rotated bounding box.

    The longer bounding-box side is L1 / dir1; the shorter is L2 / dir2.

    Args:
        footprints_gdf: GeoDataFrame or GeoSeries of Polygon footprints.
        mode: As for :func:`inertia`.
        direction: As for :func:`inertia`.

    Returns:
        As described under *mode*.
    """
    gdf = to_gdf(footprints_gdf)
    gdf = ensure_projected(gdf)
    validate_geodataframe(gdf, context="direction.bbox")

    if direction is not None:
        L1, dir1, L2, dir2 = _apply_forced_direction(gdf, direction)
    else:
        L1, dir1, L2, dir2 = min_bounding_box(gdf)

    return _pack(L1, dir1, L2, dir2, mode)


# ─────────────────────────────────────────────────────────────────────────────
# Eccentricity direction sub-namespace
# ─────────────────────────────────────────────────────────────────────────────


class _EccentricityDirectionMethod:
    """Worst-case eccentricity *direction* for a specific seismic code.

    Returns the optimum analysis angle *x_opt* (radians from principal axis)
    that maximises the eccentricity ratio, as a list aligned with the input gdf.

    Keyword arguments ``I1``, ``dir1``, ``I2`` allow passing precomputed
    principal inertia values to avoid redundant work.

    Args:
        norm: ``"EC8"`` or ``"CSCR2010"``.
    """

    def __init__(self, norm: str):
        """Bind this callable to one seismic code's eccentricity formula.

        Args:
            norm: ``"EC8"`` or ``"CSCR2010"``.
        """
        if norm not in ("EC8", "CSCR2010"):
            raise ValueError(f"Unsupported eccentricity norm: {norm!r}")
        self.norm = norm

    def __call__(
        self,
        footprints_gdf: gpd.GeoDataFrame | gpd.GeoSeries,
        *,
        I1: np.ndarray | None = None,
        dir1: np.ndarray | None = None,
        I2: np.ndarray | None = None,
        **kwargs,
    ) -> list[float]:
        """Return worst-case eccentricity analysis angle (radians) per building.

        Args:
            footprints_gdf: GeoDataFrame or GeoSeries of Polygon footprints.
            I1:   Precomputed larger principal moments (N,). Optional.
            dir1: Precomputed eigenvectors for I1 (N, 2). Optional.
            I2:   Precomputed smaller principal moments (N,). Optional.

        Returns:
            List of optimum analysis angles in radians, same order as input.
        """
        from .eccentricity import optimise_ec8, optimise_cscr
        from .geometry import (
            ensure_projected,
            to_gdf,
            calc_principal_inertia,
            centre_of_mass_and_stiffness,
        )

        gdf = ensure_projected(to_gdf(footprints_gdf))

        if I1 is None or dir1 is None or I2 is None:
            I1, dir1, I2, _ = calc_principal_inertia(gdf.geometry)

        cm, cs = centre_of_mass_and_stiffness(gdf.geometry)
        e_vec = cm - cs
        area = gdf.geometry.area.values

        if self.norm == "EC8":
            _, _, x_opt, _ = optimise_ec8(I1, dir1, I2, e_vec, area)
        else:  # CSCR2010
            _, x_opt = optimise_cscr(I1, dir1, I2, e_vec, area)

        return list(x_opt)


class _EccentricityDirectionNamespace:
    """``direction.eccentricity`` namespace.

    Usage::

        direction.eccentricity.EC8(gdf)
        direction.eccentricity.CSCR2010(gdf)
        direction.eccentricity.EC8(gdf, I1=I1, dir1=dir1, I2=I2)
    """

    EC8 = _EccentricityDirectionMethod("EC8")
    CSCR2010 = _EccentricityDirectionMethod("CSCR2010")


#: Sub-namespace for eccentricity-direction queries.
eccentricity = _EccentricityDirectionNamespace()
