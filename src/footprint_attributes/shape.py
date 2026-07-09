"""
Shape irregularity indices — code-independent and seismic-code-specific.

Access pattern
--------------
``shape(footprints_gdf, columns)``   → batch call, returns gdf with requested columns.
``shape.EC8(footprints_gdf)``        → all EC8 parameters returned in gdf.
``shape.EC8.eccentricityRatio(gdf)`` → single parameter as list.
``shape.EC8.eccentricityRatio.limits`` → code limit table.
``shape.polsby_popper(gdf)``         → code-independent index as list.
``shape.slenderness(gdf)``           → all slenderness columns in gdf.
``shape.slenderness.inertia(gdf)``   → single-method slenderness as list.

Column naming: ``{NORM}_{code}_{paramName}``
Compliance:    ``compliance_{NORM}_{code}_{paramName}``

Seismic codes supported
-----------------------
- EC8      Eurocode 8
- ASCE7    ASCE 7
- GNDTII   Italian GNDT Level II
- CSCR2010 Costa Rica Seismic Code
- NTC23    Mexican NTC-23

Code-independent indices
------------------------
- polsby_popper          — 4π A / P²  (compactness, 0–1)
- convex_hull_irregularity — momentum-weighted deviation from convex hull
- inertia_circle_ratio   — I_z(circle) / I_z(footprint)
"""

from __future__ import annotations

import numpy as np
import geopandas as gpd

from .geometry import (
    ensure_projected,
    to_gdf,
    validate_geodataframe,
    fill_holes,
    calc_inertia_z,
    calc_principal_inertia,
    centre_of_mass_and_stiffness,
    largest_convex_hull_gap_area,
    main_element_a_lengths_batch,
    hole_h_over_l,
    max_hole_area_ratio,
    setback_gndt_metrics,
    angle_between_0_90,
    eq_circle_inertia,
    inertia_side_lengths,
    bearing_from_dir,
)
from .direction import inertia as compute_inertia_direction
from .direction import bbox as compute_bbox_direction
from .config import (
    EC8_LIMITS,
    ASCE7_LIMITS,
    GNDTII_LIMITS,
    CSCR2010_LIMITS,
    NTC23_LIMITS,
    SLENDERNESS_LIMITS,
    compliance_score,
)
from .eccentricity import optimise_ec8, optimise_cscr


# ─────────────────────────────────────────────────────────────────────────────
# Direction-method dispatch
# ─────────────────────────────────────────────────────────────────────────────

#: Valid values for the ``method`` keyword accepted by size-dependent
#: parameters below (setbacks, GNDT a/b/c, slenderness-derived ratios).
_DIRECTION_METHODS = {"bbox", "inertia"}


def _basic_lengths(
    gdf: gpd.GeoDataFrame, method: str
) -> tuple[list, np.ndarray, list, np.ndarray]:
    """Dispatch to bbox or inertia direction, returning ``(L1, dir1, L2, dir2)``.

    Args:
        gdf: Projected GeoDataFrame of footprints.
        method: ``"bbox"`` (minimum rotated bounding box) or ``"inertia"``
            (exact closed-form second moment of area, paper eqs. 2-3).

    Returns:
        ``(L1, dir1, L2, dir2)`` as produced by :func:`direction.bbox` /
        :func:`direction.inertia` with ``mode="all"`` (bearing dropped).
    """
    if method not in _DIRECTION_METHODS:
        raise ValueError(f"method must be one of {_DIRECTION_METHODS}, got {method!r}")
    fn = compute_bbox_direction if method == "bbox" else compute_inertia_direction
    L1, dir1, L2, dir2, _ = fn(gdf, mode="all")
    return L1, dir1, L2, dir2


# ─────────────────────────────────────────────────────────────────────────────
# Parameter base class
# ─────────────────────────────────────────────────────────────────────────────


class Parameter:
    """Base class for a single computed parameter.

    Subclasses override :meth:`compute`.  Calling the instance is a shorthand
    for ``compute()``.

    Attributes:
        name:        Short camelCase name of the parameter.
        column_name: Full column name following ``{NORM}_{code}_{paramName}``.
        limits:      Compliance table (list of dicts) or ``None``.
        description: Human-readable description.
    """

    def __init__(
        self,
        name: str,
        column_name: str,
        limits: list[dict] | None = None,
        description: str = "",
    ):
        """Register this parameter's name, column, limits, and description."""
        self.name = name
        self.column_name = column_name
        self.limits = limits
        self.description = description

    def compute(self, gdf: gpd.GeoDataFrame, **kwargs) -> list:
        """Compute parameter values for each building.  Override in subclass."""
        raise NotImplementedError

    def __call__(self, gdf: gpd.GeoDataFrame, **kwargs) -> list:
        """Shorthand for :meth:`compute`."""
        return self.compute(gdf, **kwargs)


# ─────────────────────────────────────────────────────────────────────────────
# EC8  (Eurocode 8)
# ─────────────────────────────────────────────────────────────────────────────


class EC8EccentricityRatio(Parameter):
    """EC8 §4.2.3.2 – eccentricity ratio  e / r_t  (≤ 0.30 for regular).

    The worst-case analysis direction is found via Mohr's-circle optimisation
    (see :mod:`eccentricity`).
    """

    def __init__(self):
        """Register this parameter's name, column, limits, and description."""
        super().__init__(
            "eccentricityRatio",
            "EC8_eccentricityRatio",
            EC8_LIMITS["eccentricityRatio"],
            "Ratio of worst-case eccentricity to torsional radius (EC8)",
        )

    def compute(
        self,
        gdf: gpd.GeoDataFrame,
        *,
        I1: np.ndarray | None = None,
        dir1: np.ndarray | None = None,
        I2: np.ndarray | None = None,
        **kwargs,
    ) -> list:
        """Compute this parameter's values for each building; see the class docstring for the formula."""
        gdf = ensure_projected(to_gdf(gdf))
        if I1 is None or dir1 is None or I2 is None:
            I1, dir1, I2, _ = calc_principal_inertia(gdf.geometry)
        cm, cs = centre_of_mass_and_stiffness(gdf.geometry)
        e_vec = cm - cs
        area = gdf.geometry.area.values
        ecc_ratio, _, _, _ = optimise_ec8(I1, dir1, I2, e_vec, area)
        return list(ecc_ratio)


class EC8RadiusRatio(Parameter):
    """EC8 §4.2.3.2 – radius ratio  r_t / r_g  (≥ 1.0 for regular)."""

    def __init__(self):
        """Register this parameter's name, column, limits, and description."""
        super().__init__(
            "radiusRatio",
            "EC8_radiusRatio",
            EC8_LIMITS["radiusRatio"],
            "Ratio of torsional radius to radius of gyration (EC8)",
        )

    def compute(
        self,
        gdf: gpd.GeoDataFrame,
        *,
        I1: np.ndarray | None = None,
        dir1: np.ndarray | None = None,
        I2: np.ndarray | None = None,
        **kwargs,
    ) -> list:
        """Compute this parameter's values for each building; see the class docstring for the formula."""
        gdf = ensure_projected(to_gdf(gdf))
        if I1 is None or dir1 is None or I2 is None:
            I1, dir1, I2, _ = calc_principal_inertia(gdf.geometry)
        cm, cs = centre_of_mass_and_stiffness(gdf.geometry)
        e_vec = cm - cs
        area = gdf.geometry.area.values
        _, rad_ratio, _, _ = optimise_ec8(I1, dir1, I2, e_vec, area)
        return list(rad_ratio)


class EC8Compactness(Parameter):
    """EC8 §4.2.3.2 – compactness  1 – (A_setback / A_total)  (≥ 0.95)."""

    def __init__(self):
        """Register this parameter's name, column, limits, and description."""
        super().__init__(
            "compactness",
            "EC8_compactness",
            EC8_LIMITS["compactness"],
            "1 – (largest convex-hull setback area / footprint area) (EC8)",
        )

    def compute(self, gdf: gpd.GeoDataFrame, **kwargs) -> list:
        """Compute this parameter's values for each building; see the class docstring for the formula."""
        gdf = ensure_projected(to_gdf(gdf))
        # Per the paper: take the set difference between the (hole-filled)
        # footprint and its convex hull, and use the area of the single
        # LARGEST resulting piece (not the sum of all setback areas).
        gap_areas = largest_convex_hull_gap_area(gdf.geometry)
        filled_areas = fill_holes(gdf.geometry).area.values
        return [
            1.0 - (gap / (a + 1e-12)) if a > 1e-12 else 1.0
            for gap, a in zip(gap_areas, filled_areas)
        ]


# ─────────────────────────────────────────────────────────────────────────────
# ASCE 7
# ─────────────────────────────────────────────────────────────────────────────


class ASCE7SetbackRatio(Parameter):
    """ASCE 7 Table 12.3-1 – setback ratio  min(b1/L1, b2/L2)  (≤ 0.20).

    ``method`` selects which principal-axis convention (``"bbox"`` default,
    or ``"inertia"``) defines dir1/dir2/L1/L2 for the dual-configuration
    setback construction (paper §3.4.6 / fig. 10a).
    """

    def __init__(self):
        """Register this parameter's name, column, limits, and description."""
        super().__init__(
            "setbackRatio",
            "ASCE7_setbackRatio",
            ASCE7_LIMITS["setbackRatio"],
            "min(b1/L1, b2/L2) dual-configuration setback ratio (ASCE 7)",
        )

    def compute(
        self,
        gdf: gpd.GeoDataFrame,
        *,
        method: str = "bbox",
        _gndt_setback: tuple | None = None,
        **kwargs,
    ) -> list:
        """Compute this parameter's values for each building; see the class docstring for the formula."""
        gdf = ensure_projected(to_gdf(gdf))
        if _gndt_setback is not None:
            ratio, _, _ = _gndt_setback
        else:
            L1, dir1, L2, dir2 = _basic_lengths(gdf, method)
            ratio, _, _ = setback_gndt_metrics(gdf, L1, dir1, L2, dir2)
        return ratio


class ASCE7HoleRatio(Parameter):
    """ASCE 7 Table 12.3-1 – hole ratio  A_hole / A_filled  (≤ 0.25)."""

    def __init__(self):
        """Register this parameter's name, column, limits, and description."""
        super().__init__(
            "holeRatio",
            "ASCE7_holeRatio",
            ASCE7_LIMITS["holeRatio"],
            "Max interior hole area / filled area (ASCE 7)",
        )

    def compute(self, gdf: gpd.GeoDataFrame, **kwargs) -> list:
        """Compute this parameter's values for each building; see the class docstring for the formula."""
        return max_hole_area_ratio(ensure_projected(to_gdf(gdf)))


class ASCE7ParalelityAngle(Parameter):
    """ASCE 7 – parallelity angle (degrees)  (≤ 5° for regular)."""

    def __init__(self):
        """Register this parameter's name, column, limits, and description."""
        super().__init__(
            "parallelityAngle",
            "ASCE7_parallelityAngle",
            ASCE7_LIMITS["parallelityAngle"],
            "Angle between bounding-box sides and cardinal axes (ASCE 7)",
        )

    def compute(
        self,
        gdf: gpd.GeoDataFrame,
        *,
        dir1: np.ndarray | None = None,
        **kwargs,
    ) -> list:
        """Compute this parameter's values for each building; see the class docstring for the formula."""
        gdf = ensure_projected(to_gdf(gdf))
        if dir1 is None:
            _, dir1, _, _, _ = compute_bbox_direction(gdf, mode="all")
        return [np.degrees(angle_between_0_90(np.array([1.0, 0.0]), d)) for d in dir1]


# ─────────────────────────────────────────────────────────────────────────────
# GNDTII  (Italian GNDT Level II)
# ─────────────────────────────────────────────────────────────────────────────


def _gndt_dominant_a(
    gdf: gpd.GeoDataFrame, method: str, _basic: tuple | None = None
) -> tuple[list, list, np.ndarray, np.ndarray]:
    """Shared (a, L, dir1, dir2) dominant-configuration pick for β1/β4.

    Computes a1 (paired with L1) and a2 (paired with L2) via the inscribed-
    circle construction (paper fig. 7), then selects whichever configuration
    maximises ``L * a`` (paper §3.4.5), per building.

    Args:
        gdf: Projected GeoDataFrame of footprints.
        method: ``"bbox"`` or ``"inertia"``.
        _basic: Optional precomputed ``(L1, dir1, L2, dir2)`` for *method*
            (see :class:`_SharedGeometryCache`), to avoid re-deriving the
            bbox/inertia axes when a caller already has them cached.

    Returns:
        ``(a_dom, L_dom, dir1, dir2)`` -- the winning a and L per building,
        plus the L1/L2 axis vectors (kept for callers that also need them).
    """
    L1, dir1, L2, dir2 = _basic if _basic is not None else _basic_lengths(gdf, method)
    a1, a2, _centers = main_element_a_lengths_batch(gdf, dir1, dir2)
    a_dom, L_dom = [], []
    for L1v, a1v, L2v, a2v in zip(L1, a1, L2, a2):
        if L1v * a1v >= L2v * a2v:
            a_dom.append(a1v)
            L_dom.append(L1v)
        else:
            a_dom.append(a2v)
            L_dom.append(L2v)
    return a_dom, L_dom, dir1, dir2


# ─────────────────────────────────────────────────────────────────────────────
# Shared-computation cache -- lets a single shape()/run() call compute each
# expensive per-building construction (bbox/inertia axes, GNDT 'a',
# GNDT setback) at most once, no matter how many requested columns/norms
# need it. Norm-scoped caching alone (the previous approach) still redoes
# the same construction once per norm: ASCE7/GNDTII/NTC-23's setback ratios
# all default to the same bbox-method setback_gndt_metrics() call, so a
# config requesting all three (a realistic combination -- see
# runner.run()'s docstring example) used to run it 3 times over.
# ─────────────────────────────────────────────────────────────────────────────


class _SharedGeometryCache:
    """Lazily computes and memoises bbox/inertia axes and GNDT constructions
    for one GeoDataFrame, so every :class:`Parameter` (across every norm)
    asked to compute against it shares the same underlying work.

    A fresh instance is cheap to create (it does nothing until first used),
    so :meth:`NormAggregate.__call__` creates its own when called directly
    (e.g. ``shape.EC8(gdf)``) and :func:`_ShapeModule.__call__` (``shape()``,
    and hence ``run()``) creates exactly one and shares it across every norm
    /slenderness/bearing computation in that call.
    """

    def __init__(self, gdf: gpd.GeoDataFrame):
        """Register this parameter's name, column, limits, and description."""
        self._gdf = gdf
        self._bbox: dict | None = None
        self._inertia: dict | None = None
        self._gndt_a: dict[str, tuple] = {}
        self._gndt_setback: dict[str, tuple] = {}

    def bbox(self) -> dict:
        """Cached ``{L1, dir1, L2, dir2}`` from the minimum bounding box."""
        if self._bbox is None:
            L1, dir1, L2, dir2, _ = compute_bbox_direction(self._gdf, mode="all")
            self._bbox = dict(L1=L1, dir1=dir1, L2=L2, dir2=dir2)
        return self._bbox

    def inertia(self) -> dict:
        """Raw ``calc_principal_inertia`` output (I1/dir1 <-> larger
        eigenvalue), the convention EC8's eccentricity/radius-ratio
        parameters expect -- NOT the physical-length ``dir1`` swap that
        :func:`direction.inertia` applies (see :meth:`bearing`, which
        applies that swap itself when deriving the bearing from this same
        cached computation).
        """
        if self._inertia is None:
            I1, dir1, I2, dir2 = calc_principal_inertia(self._gdf.geometry)
            self._inertia = dict(I1=I1, dir1=dir1, I2=I2, dir2=dir2)
        return self._inertia

    def basic_lengths(self, method: str) -> tuple:
        """``(L1, dir1, L2, dir2)`` for *method*, from the cached bbox/inertia
        computation -- matches :func:`_basic_lengths` but never triggers a
        second ``calc_principal_inertia`` call for the inertia method.
        """
        if method == "bbox":
            b = self.bbox()
            return b["L1"], b["dir1"], b["L2"], b["dir2"]
        i = self.inertia()
        area = self._gdf.geometry.area.values
        L1, L2 = inertia_side_lengths(i["I1"], i["I2"], area)
        # direction.inertia()'s physical-length swap: the larger-eigenvalue
        # eigenvector points along the geometrically *shorter* side.
        return L1, i["dir2"], L2, i["dir1"]

    def gndt_a(self, method: str) -> tuple:
        """Cached ``_gndt_dominant_a(gdf, method)`` result for *method*.

        Args:
            method: ``"bbox"`` or ``"inertia"``.

        Returns:
            ``(a, L)`` as returned by :func:`_gndt_dominant_a`.
        """
        if method not in self._gndt_a:
            self._gndt_a[method] = _gndt_dominant_a(
                self._gdf, method, _basic=self.basic_lengths(method)
            )
        return self._gndt_a[method]

    def gndt_setback(self, method: str) -> tuple:
        """Cached ``setback_gndt_metrics(gdf, ...)`` result for *method*.

        Args:
            method: ``"bbox"`` or ``"inertia"``.

        Returns:
            The tuple returned by
            :func:`footprint_attributes.geometry.setback_gndt_metrics`.
        """
        if method not in self._gndt_setback:
            L1, dir1, L2, dir2 = self.basic_lengths(method)
            self._gndt_setback[method] = setback_gndt_metrics(
                self._gdf, L1, dir1, L2, dir2
            )
        return self._gndt_setback[method]

    def bearing(self) -> list:
        """Building bearing (degrees from North), exactly as
        :func:`direction.inertia` computes it, but reusing this cache's
        already-computed ``calc_principal_inertia`` result instead of
        recomputing it from scratch.
        """
        dir1_raw = self.inertia()["dir1"]
        return [bearing_from_dir(d) for d in dir1_raw]


class GNDTIIBeta1MainShapeSlenderness(Parameter):
    """GNDTII β₁ – a / L : masonry main-shape slenderness.

    ``a`` comes from the inscribed-circle construction (paper fig. 7);
    ``method`` (``"bbox"`` default, or ``"inertia"``) selects which
    principal-axis convention defines dir1/dir2/L1/L2, and hence which of
    the two ``(L, a)`` configurations is available. The dominant
    configuration is the one maximising ``L * a``.
    """

    def __init__(self):
        """Register this parameter's name, column, limits, and description."""
        super().__init__(
            "beta1_mainShapeSlenderness",
            "GNDTII_beta1_mainShapeSlenderness",
            GNDTII_LIMITS["beta1_mainShapeSlenderness"],
            "a / L, dominant (L, a) configuration by max(L*a) (GNDTII)",
        )

    def compute(
        self,
        gdf: gpd.GeoDataFrame,
        *,
        method: str = "bbox",
        _gndt_a: tuple | None = None,
        **kwargs,
    ) -> list:
        """Compute this parameter's values for each building; see the class docstring for the formula."""
        gdf = ensure_projected(to_gdf(gdf))
        a_dom, L_dom, _, _ = (
            _gndt_a if _gndt_a is not None else _gndt_dominant_a(gdf, method)
        )
        return [a / (L + 1e-12) if L > 1e-12 else 0.0 for a, L in zip(a_dom, L_dom)]


class GNDTIIBeta2SetbackRatio(Parameter):
    """GNDTII β₂ – min(b1/L1, b2/L2) : masonry setback ratio.

    ``method`` (``"bbox"`` default, or ``"inertia"``) selects which
    principal-axis convention defines dir1/dir2/L1/L2 for the dual-
    configuration setback construction.
    """

    def __init__(self):
        """Register this parameter's name, column, limits, and description."""
        super().__init__(
            "beta2_setbackRatio",
            "GNDTII_beta2_setbackRatio",
            GNDTII_LIMITS["beta2_setbackRatio"],
            "min(b1/L1, b2/L2) dual-configuration setback ratio (GNDTII)",
        )

    def compute(
        self,
        gdf: gpd.GeoDataFrame,
        *,
        method: str = "bbox",
        _gndt_setback: tuple | None = None,
        **kwargs,
    ) -> list:
        """Compute this parameter's values for each building; see the class docstring for the formula."""
        gdf = ensure_projected(to_gdf(gdf))
        if _gndt_setback is not None:
            ratio, _, _ = _gndt_setback
        else:
            L1, dir1, L2, dir2 = _basic_lengths(gdf, method)
            ratio, _, _ = setback_gndt_metrics(gdf, L1, dir1, L2, dir2)
        return ratio


class GNDTIIBeta4EccentricityRatio(Parameter):
    """GNDTII β₄ – e / a : concrete eccentricity ratio.

    ``a`` is the same dominant-configuration inscribed-circle length used
    by β₁ (paper §3.4.5); ``method`` (``"bbox"`` default, or ``"inertia"``)
    selects the underlying principal-axis convention.
    """

    def __init__(self):
        """Register this parameter's name, column, limits, and description."""
        super().__init__(
            "beta4_eccentricityRatio",
            "GNDTII_beta4_eccentricityRatio",
            GNDTII_LIMITS["beta4_eccentricityRatio"],
            "Eccentricity / dominant-configuration a (GNDTII)",
        )

    def compute(
        self,
        gdf: gpd.GeoDataFrame,
        *,
        method: str = "bbox",
        _gndt_a: tuple | None = None,
        **kwargs,
    ) -> list:
        """Compute this parameter's values for each building; see the class docstring for the formula."""
        gdf = ensure_projected(to_gdf(gdf))
        a_dom, _, _, _ = (
            _gndt_a if _gndt_a is not None else _gndt_dominant_a(gdf, method)
        )
        cm, cs = centre_of_mass_and_stiffness(gdf.geometry)
        e_mag = np.linalg.norm(cm - cs, axis=1)
        return [e / (a + 1e-12) if a > 1e-12 else 0.0 for e, a in zip(e_mag, a_dom)]


class GNDTIIBeta6SetbackSlenderness(Parameter):
    """GNDTII β₆ – c / b : setback slenderness (protrusion depth / setback width).

    ``method`` (``"bbox"`` default, or ``"inertia"``) selects which
    principal-axis convention defines dir1/dir2/L1/L2 for the dual-
    configuration setback construction; b/c come from the same winning
    configuration as β₂.
    """

    def __init__(self):
        """Register this parameter's name, column, limits, and description."""
        super().__init__(
            "beta6_setbackSlenderness",
            "GNDTII_beta6_setbackSlenderness",
            GNDTII_LIMITS["beta6_setbackSlenderness"],
            "Protrusion depth c / winning-configuration setback width b (GNDTII)",
        )

    def compute(
        self,
        gdf: gpd.GeoDataFrame,
        *,
        method: str = "bbox",
        _gndt_setback: tuple | None = None,
        **kwargs,
    ) -> list:
        """Compute this parameter's values for each building; see the class docstring for the formula."""
        gdf = ensure_projected(to_gdf(gdf))
        if _gndt_setback is not None:
            _, b, c = _gndt_setback
        else:
            L1, dir1, L2, dir2 = _basic_lengths(gdf, method)
            _, b, c = setback_gndt_metrics(gdf, L1, dir1, L2, dir2)
        return [cv / (bv + 1e-12) if bv > 1e-12 else 0.0 for bv, cv in zip(b, c)]


# ─────────────────────────────────────────────────────────────────────────────
# CSCR 2010  (Costa Rica)
# ─────────────────────────────────────────────────────────────────────────────


class CSCR2010EccentricityRatio(Parameter):
    """CSCR 2010 – worst-case eccentricity ratio  e / l."""

    def __init__(self):
        """Register this parameter's name, column, limits, and description."""
        super().__init__(
            "eccentricityRatio",
            "CSCR2010_eccentricityRatio",
            CSCR2010_LIMITS["eccentricityRatio"],
            "Worst-case eccentricity / building dimension (CSCR 2010)",
        )

    def compute(
        self,
        gdf: gpd.GeoDataFrame,
        *,
        I1: np.ndarray | None = None,
        dir1: np.ndarray | None = None,
        I2: np.ndarray | None = None,
        **kwargs,
    ) -> list:
        """Compute this parameter's values for each building; see the class docstring for the formula."""
        gdf = ensure_projected(to_gdf(gdf))
        if I1 is None or dir1 is None or I2 is None:
            I1, dir1, I2, _ = calc_principal_inertia(gdf.geometry)
        cm, cs = centre_of_mass_and_stiffness(gdf.geometry)
        e_vec = cm - cs
        area = gdf.geometry.area.values
        ecc_ratio, _ = optimise_cscr(I1, dir1, I2, e_vec, area)
        return list(ecc_ratio)


# ─────────────────────────────────────────────────────────────────────────────
# NTC-23  (Mexico)
# ─────────────────────────────────────────────────────────────────────────────


class NTC23SetbackRatio(Parameter):
    """NTC-23 – setback ratio  min(b1/L1, b2/L2)  (≤ 0.40).

    Same dual-configuration construction as ASCE7's setback ratio, just
    against the more lenient NTC-23 limit. ``method`` (``"bbox"`` default,
    or ``"inertia"``) selects the underlying principal-axis convention.
    """

    def __init__(self):
        """Register this parameter's name, column, limits, and description."""
        super().__init__(
            "setbackRatio",
            "NTC23_setbackRatio",
            NTC23_LIMITS["setbackRatio"],
            "min(b1/L1, b2/L2) dual-configuration setback ratio (NTC-23)",
        )

    def compute(
        self,
        gdf: gpd.GeoDataFrame,
        *,
        method: str = "bbox",
        _gndt_setback: tuple | None = None,
        **kwargs,
    ) -> list:
        """Compute this parameter's values for each building; see the class docstring for the formula."""
        gdf = ensure_projected(to_gdf(gdf))
        if _gndt_setback is not None:
            ratio, _, _ = _gndt_setback
        else:
            L1, dir1, L2, dir2 = _basic_lengths(gdf, method)
            ratio, _, _ = setback_gndt_metrics(gdf, L1, dir1, L2, dir2)
        return ratio


class NTC23HoleRatio(Parameter):
    """NTC-23 – hole ratio  h / L  (≤ 0.40).

    Distinct from ASCE7's area-based hole ratio: here ``h`` is the shorter
    side of each hole's OWN minimum bounding box (generally not parallel to
    the building's own axes), and ``L`` is the length of the building's
    cross-section through the hole's centroid along the hole's longer axis
    (paper §3.4.7, fig. 10b).
    """

    def __init__(self):
        """Register this parameter's name, column, limits, and description."""
        super().__init__(
            "holeRatio",
            "NTC23_holeRatio",
            NTC23_LIMITS["holeRatio"],
            "Worst hole's own-MBB width / through-centroid footprint length (NTC-23)",
        )

    def compute(self, gdf: gpd.GeoDataFrame, **kwargs) -> list:
        """Compute this parameter's values for each building; see the class docstring for the formula."""
        return hole_h_over_l(ensure_projected(to_gdf(gdf)))


# ─────────────────────────────────────────────────────────────────────────────
# Code-independent shape indices
# ─────────────────────────────────────────────────────────────────────────────


def polsby_popper(
    geoms: gpd.GeoDataFrame | gpd.GeoSeries,
    fill_holes_flag: bool = True,
) -> list[float]:
    """Polsby-Popper compactness index  4π A / P².

    Values near 1 indicate compact (near-circular) shapes; lower values
    indicate elongated or irregular footprints.

    Args:
        geoms: GeoDataFrame or GeoSeries of footprint polygons.
        fill_holes_flag: If ``True`` (default), interior courtyards are ignored.

    Returns:
        List of Polsby-Popper values in (0, 1], same order as input.
    """
    gdf = ensure_projected(to_gdf(geoms))
    validate_geodataframe(gdf, context="shape.polsby_popper")
    gs = fill_holes(gdf.geometry) if fill_holes_flag else gdf.geometry
    return list((4.0 * np.pi * gs.area) / (gs.boundary.length**2 + 1e-30))


def convex_hull_irregularity(
    geoms: gpd.GeoDataFrame | gpd.GeoSeries,
) -> list[float]:
    """Convex-hull area-excess ratio: how much bigger the hull is than the footprint.

    ``(hull_area - footprint_area) / footprint_area``. 0 for a convex shape
    (hull == footprint); higher values indicate deeper and/or larger
    setbacks relative to the footprint's own area. Interior holes are
    filled before computing both areas, since this index is about plan
    *shape* (setbacks), not holes -- see ``ASCE7_holeRatio`` /
    ``NTC23_holeRatio`` for hole-specific metrics.

    Args:
        geoms: GeoDataFrame or GeoSeries of footprint polygons.

    Returns:
        List of irregularity values (dimensionless, >= 0), same order as input.
    """
    gdf = ensure_projected(to_gdf(geoms))
    validate_geodataframe(gdf, context="shape.convex_hull_irregularity")
    filled = fill_holes(gdf.geometry)
    hull_area = filled.convex_hull.area
    area = filled.area
    return list((hull_area - area) / (area + 1e-30))


def inertia_circle_ratio(
    geoms: gpd.GeoDataFrame | gpd.GeoSeries,
) -> list[float]:
    """Ratio of the equivalent-circle polar inertia to the footprint's I_z.

    A circle with the same area as the footprint has the maximum possible I_z
    for that area.  The ratio I_circle / I_footprint falls in (0, 1] and
    equals 1 for a perfect circle; irregular shapes have lower values.

    Args:
        geoms: GeoDataFrame or GeoSeries of footprint polygons.

    Returns:
        List of inertia-circle ratios in (0, 1], same order as input.
    """
    gdf = ensure_projected(to_gdf(geoms))
    validate_geodataframe(gdf, context="shape.inertia_circle_ratio")
    areas = gdf.geometry.area.values
    iz = np.abs(calc_inertia_z(gdf.geometry))
    circle_iz = eq_circle_inertia(areas)
    return list(circle_iz / (iz + 1e-30))


# ─────────────────────────────────────────────────────────────────────────────
# Slenderness (separate namespace, applies to multiple codes)
# ─────────────────────────────────────────────────────────────────────────────


class _SlendernessMethod(Parameter):
    """A single slenderness computation using one direction method."""

    def __init__(self, direction_method: str):
        """Register this parameter's name, column, limits, and description."""
        super().__init__(
            "planSlenderness",
            f"slenderness_{direction_method}",
            None,
            f"Plan slenderness L1/L2 via {direction_method} method",
        )
        self.direction_method = direction_method
        # Compliance limits come from EC8 slenderness table
        self._compliance_limits = SLENDERNESS_LIMITS.get("EC8")

    def compute(
        self,
        gdf: gpd.GeoDataFrame,
        *,
        L1: list | None = None,
        L2: list | None = None,
        vertical: bool = False,
        height_column: str | None = None,
        **kwargs,
    ) -> list:
        """Compute this parameter's values for each building; see the class docstring for the formula."""
        gdf = ensure_projected(to_gdf(gdf))

        if vertical:
            # Vertical slenderness = height / L2
            if L2 is None:
                if self.direction_method == "inertia":
                    _, _, L2, _, _ = compute_inertia_direction(gdf, mode="all")
                else:
                    _, _, L2, _, _ = compute_bbox_direction(gdf, mode="all")
            col = height_column or "height"
            if col not in gdf.columns:
                raise ValueError(
                    f"height_column '{col}' not found in GeoDataFrame. "
                    "Provide a height column or pass height_column= keyword."
                )
            heights = gdf[col].astype(float).values
            return [h / (L + 1e-12) if L > 1e-12 else 0.0 for h, L in zip(heights, L2)]

        # Plan slenderness = L1 / L2
        if L1 is None or L2 is None:
            if self.direction_method == "inertia":
                L1, _, L2, _, _ = compute_inertia_direction(gdf, mode="all")
            else:
                L1, _, L2, _, _ = compute_bbox_direction(gdf, mode="all")
        return [L1v / (L2v + 1e-12) if L2v > 1e-12 else 0.0 for L1v, L2v in zip(L1, L2)]


class SlendernessAccessor:
    """Namespace for plan and vertical slenderness computation.

    Usage::

        shape.slenderness(gdf)               # all methods, returns gdf
        shape.slenderness.inertia(gdf)       # list, plan slenderness via inertia
        shape.slenderness.bbox(gdf)          # list, plan slenderness via bbox
        shape.slenderness.inertia(gdf, vertical=True, height_column="height")
    """

    def __init__(self):
        """Register this parameter's name, column, limits, and description."""
        self._methods: dict[str, _SlendernessMethod] = {}
        # Register built-in methods
        self._register("inertia", _SlendernessMethod("inertia"))
        self._register("bbox", _SlendernessMethod("bbox"))

    def _register(self, name: str, method: _SlendernessMethod) -> None:
        """Add *method* to the method registry and expose it as ``self.<name>``.

        Args:
            name: Attribute name the method is exposed under (e.g. ``"bbox"``).
            method: The :class:`_SlendernessMethod` instance to register.
        """
        self._methods[name] = method
        setattr(self, name, method)

    def __call__(
        self,
        gdf: gpd.GeoDataFrame,
        **kwargs,
    ) -> gpd.GeoDataFrame:
        """Compute all registered slenderness methods and add columns to gdf."""
        gdf = ensure_projected(to_gdf(gdf))
        validate_geodataframe(gdf, context="shape.slenderness")
        result = gdf.copy()

        for method_name, method_obj in self._methods.items():
            col = method_obj.column_name
            if col not in result.columns:
                result[col] = method_obj.compute(result, **kwargs)

            # Compliance against EC8 slenderness limit
            comp_col = f"compliance_EC8_{col}"
            if method_obj._compliance_limits is not None:
                result[comp_col] = [
                    compliance_score(v, method_obj._compliance_limits)[0]
                    for v in result[col]
                ]

        return result


slenderness = SlendernessAccessor()


# ─────────────────────────────────────────────────────────────────────────────
# NormAggregate — container for all parameters of a single seismic code
# ─────────────────────────────────────────────────────────────────────────────


class NormAggregate:
    """Container for all parameters of a single seismic code.

    Calling the instance computes all parameters and adds them (plus compliance
    columns) to the input GeoDataFrame.  Individual parameters are accessible
    as attributes.

    Args:
        name:       Norm identifier (e.g. ``"EC8"``).
        parameters: Mapping of attribute-name → :class:`Parameter` instance.
    """

    def __init__(self, name: str, parameters: dict[str, Parameter]):
        """Register this parameter's name, column, limits, and description."""
        self.name = name
        self.parameters = parameters
        for attr, param_obj in parameters.items():
            setattr(self, attr, param_obj)

    def __call__(
        self,
        gdf: gpd.GeoDataFrame,
        _cache: "_SharedGeometryCache | None" = None,
        _columns: set[str] | None = None,
        **kwargs,
    ) -> gpd.GeoDataFrame:
        """Compute all parameters; skip any column already present in *gdf*.

        Args:
            _cache: Internal. A :class:`_SharedGeometryCache` to reuse
                (passed by :func:`_ShapeModule.__call__` when several norms
                are being computed in the same ``shape()``/``run()`` call,
                so they all share one bbox/inertia/GNDT computation instead
                of each norm redoing it). Callers using a norm directly
                (e.g. ``shape.EC8(gdf)``) leave this as ``None`` and get a
                cache scoped to just this call, same as before.
            _columns: Internal. If given, only parameters whose own column
                (or compliance column) is in this set are computed --
                lets :func:`_ShapeModule.__call__` request e.g. just
                ``GNDTII_beta2_setbackRatio`` without also paying for
                beta1/beta4's expensive inscribed-circle construction just
                because they belong to the same norm. ``None`` (the
                default, and always the case for a direct ``shape.EC8(gdf)``
                call) computes every parameter in the norm, as before.
        """
        gdf = ensure_projected(to_gdf(gdf))
        validate_geodataframe(gdf, context=f"shape.{self.name}")
        result = gdf.copy()

        cache = _cache if _cache is not None else _SharedGeometryCache(result)

        for param_name, param_obj in self.parameters.items():
            col_name = param_obj.column_name
            comp_col = f"compliance_{col_name}"

            if (
                _columns is not None
                and col_name not in _columns
                and comp_col not in _columns
            ):
                continue

            if col_name not in result.columns:
                # Pass precomputed values where the parameter supports it
                import inspect

                sig = inspect.signature(param_obj.compute)
                extra = {}
                if sig.parameters.keys() & {"L1", "dir1", "L2", "dir2"}:
                    extra.update(cache.bbox())
                if sig.parameters.keys() & {"I1", "I2"}:
                    extra.update(cache.inertia())
                if "_gndt_a" in sig.parameters or "_gndt_setback" in sig.parameters:
                    # Resolve the method this param will actually run with:
                    # the caller's override (applies uniformly to every
                    # param in this norm call) if given, else this param's
                    # own default.
                    method = kwargs.get("method", sig.parameters["method"].default)
                    if "_gndt_a" in sig.parameters:
                        extra["_gndt_a"] = cache.gndt_a(method)
                    if "_gndt_setback" in sig.parameters:
                        extra["_gndt_setback"] = cache.gndt_setback(method)
                extra.update(kwargs)  # caller overrides take priority
                values = param_obj.compute(result, **extra)
                result[col_name] = values
            else:
                values = list(result[col_name])

            if param_obj.limits is not None and comp_col not in result.columns:
                result[comp_col] = [
                    compliance_score(v, param_obj.limits)[0] for v in values
                ]

        return result


# ─────────────────────────────────────────────────────────────────────────────
# Norm instances (these are the public API objects)
# ─────────────────────────────────────────────────────────────────────────────

EC8 = NormAggregate(
    "EC8",
    {
        "eccentricityRatio": EC8EccentricityRatio(),
        "radiusRatio": EC8RadiusRatio(),
        "compactness": EC8Compactness(),
    },
)

ASCE7 = NormAggregate(
    "ASCE7",
    {
        "setbackRatio": ASCE7SetbackRatio(),
        "holeRatio": ASCE7HoleRatio(),
        "parallelityAngle": ASCE7ParalelityAngle(),
    },
)

GNDTII = NormAggregate(
    "GNDTII",
    {
        "beta1_mainShapeSlenderness": GNDTIIBeta1MainShapeSlenderness(),
        "beta2_setbackRatio": GNDTIIBeta2SetbackRatio(),
        "beta4_eccentricityRatio": GNDTIIBeta4EccentricityRatio(),
        "beta6_setbackSlenderness": GNDTIIBeta6SetbackSlenderness(),
    },
)

CSCR2010 = NormAggregate(
    "CSCR2010",
    {
        "eccentricityRatio": CSCR2010EccentricityRatio(),
    },
)

NTC23 = NormAggregate(
    "NTC23",
    {
        "setbackRatio": NTC23SetbackRatio(),
        "holeRatio": NTC23HoleRatio(),
    },
)


# ─────────────────────────────────────────────────────────────────────────────
# Top-level shape() callable
# ─────────────────────────────────────────────────────────────────────────────

#: All norm aggregators in one place, for the batch shape() function.
_ALL_NORMS = [EC8, ASCE7, GNDTII, CSCR2010, NTC23]

#: Code-independent functions available as ``shape.polsby_popper`` etc.
#  (they are module-level functions, exposed here for autodoc convenience)

_CODE_INDEPENDENT = {
    "polsby_popper": polsby_popper,
    "convex_hull_irregularity": convex_hull_irregularity,
    "inertia_circle_ratio": inertia_circle_ratio,
}


class _ShapeModule:
    """Top-level ``shape`` object.

    Provides:
    - ``shape(gdf, columns)`` – batch call for any mix of columns.
    - ``shape.EC8(gdf)``      – all EC8 parameters.
    - ``shape.polsby_popper(gdf)`` – code-independent index.
    - ``shape.slenderness(gdf)``   – all slenderness variants.
    """

    # Norm aggregators
    EC8 = EC8
    ASCE7 = ASCE7
    GNDTII = GNDTII
    CSCR2010 = CSCR2010
    NTC23 = NTC23
    slenderness = slenderness

    # Code-independent indices (callable directly on the module object)
    polsby_popper = staticmethod(polsby_popper)
    convex_hull_irregularity = staticmethod(convex_hull_irregularity)
    inertia_circle_ratio = staticmethod(inertia_circle_ratio)

    def __call__(
        self,
        gdf: gpd.GeoDataFrame,
        columns: list[str] | None = None,
        **kwargs,
    ) -> gpd.GeoDataFrame:
        """Compute any mix of shape columns for *gdf*.

        Shared quantities (bounding-box dimensions, principal inertia, the
        GNDT inscribed-circle 'a' construction, the GNDT dual-configuration
        setback construction) are each computed at most once for the whole
        call via one shared :class:`_SharedGeometryCache`, no matter how
        many norms/columns need them -- e.g. requesting ``"ASCE7"``,
        ``"GNDTII"``, and ``"NTC23_setbackRatio"`` together (a realistic
        combination) used to run the setback construction 3 times over;
        now it runs once.

        Args:
            gdf:     GeoDataFrame of building footprints.
            columns: List of column names to compute.  Passing ``None``
                     computes all parameters from all norms, all
                     slenderness variants, the bearing, and all
                     code-independent indices. Columns already present in
                     *gdf* are not recomputed.
            **kwargs: Forwarded to individual parameter compute() calls
                     (e.g. ``height_column="h"``).

        Returns:
            GeoDataFrame with requested columns added.
        """
        gdf = ensure_projected(to_gdf(gdf))
        validate_geodataframe(gdf, context="shape")
        result = gdf.copy()
        cache = _SharedGeometryCache(result)

        # Build a map: column_name → (param_obj or callable, kind)
        column_map: dict[str, tuple] = {}

        for norm in _ALL_NORMS:
            for param_obj in norm.parameters.values():
                column_map[param_obj.column_name] = ("param", norm, param_obj)
                comp = f"compliance_{param_obj.column_name}"
                column_map[comp] = ("compliance", norm, param_obj)

        for method_name, method_obj in slenderness._methods.items():
            col = method_obj.column_name
            column_map[col] = ("slenderness_param", method_obj)
            comp = f"compliance_EC8_{col}"
            column_map[comp] = ("slenderness_compliance", method_obj)

        for name, fn in _CODE_INDEPENDENT.items():
            column_map[name] = ("independent", fn)

        column_map["bearing"] = ("bearing",)

        if columns is None:
            columns = list(column_map.keys())

        # Determine which norms / methods are needed -- and, per norm, the
        # *exact* columns wanted from it, so e.g. requesting only
        # "GNDTII_beta2_setbackRatio" doesn't also force beta1/beta4's
        # unrelated (and much more expensive) inscribed-circle construction
        # just because they happen to belong to the same norm.
        needed_norm_columns: dict[str, set[str]] = {}
        needed_slenderness: set[str] = set()
        needed_independent: set[str] = set()
        needed_bearing = False

        for col in columns:
            if col not in result.columns and col in column_map:
                entry = column_map[col]
                if entry[0] in ("param", "compliance"):
                    needed_norm_columns.setdefault(entry[1].name, set()).add(col)
                elif entry[0] in ("slenderness_param", "slenderness_compliance"):
                    needed_slenderness.add(entry[1].direction_method)
                elif entry[0] == "independent":
                    needed_independent.add(col)
                elif entry[0] == "bearing":
                    needed_bearing = True

        # Compute norms, all sharing one cache
        for norm in _ALL_NORMS:
            if norm.name in needed_norm_columns:
                result = norm(
                    result,
                    _cache=cache,
                    _columns=needed_norm_columns[norm.name],
                    **kwargs,
                )

        # Compute slenderness methods -- L1/L2 come from the same shared
        # bbox/inertia cache the norms above may have already populated.
        for method_name, method_obj in slenderness._methods.items():
            if method_name in needed_slenderness:
                col = method_obj.column_name
                if col not in result.columns:
                    extra = dict(kwargs)
                    if "L1" not in extra and "L2" not in extra:
                        L1, _, L2, _ = cache.basic_lengths(method_obj.direction_method)
                        extra["L1"], extra["L2"] = L1, L2
                    result[col] = method_obj.compute(result, **extra)
                comp = f"compliance_EC8_{col}"
                if comp not in result.columns and method_obj._compliance_limits:
                    result[comp] = [
                        compliance_score(v, method_obj._compliance_limits)[0]
                        for v in result[col]
                    ]

        # Compute code-independent indices
        for name in needed_independent:
            if name not in result.columns:
                fn = _CODE_INDEPENDENT[name]
                result[name] = fn(result, **kwargs)

        # Bearing -- reuses the cache's calc_principal_inertia() result
        # rather than a separate direction.inertia() call, when it was (or
        # is about to be) computed anyway for an inertia-based norm/column.
        if needed_bearing and "bearing" not in result.columns:
            result["bearing"] = cache.bearing()

        # Return only requested columns (plus geometry)
        final_cols = ["geometry"] + [c for c in columns if c in result.columns]
        return result[[c for c in final_cols if c in result.columns]]


# The module-level ``shape`` object that users import
shape = _ShapeModule()
