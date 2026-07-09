"""
Eccentricity optimisation helpers.

Both Eurocode 8 and CSCR 2010 require finding the *worst-case* direction of
eccentricity across all possible horizontal analysis axes.  This is a
one-variable optimisation over the rotation angle *x* that maximises the
eccentricity ratio.

The seismic inertia tensor (about the centre of stiffness) in the rotated
frame is parameterised by:

    I(x) = c - r·cos(2x)

where:
    c = (I1 + I2) / 2   (centre of Mohr's circle)
    r = (I1 - I2) / 2   (radius of Mohr's circle)
    I1, I2               principal moments of inertia (I1 ≥ I2)

The angle *b* is the angle between the eigenvector dir1 and the eccentricity
vector.

EC8 optimises:
    max  cos²(x−b) · (c − r·cos(2x))          [torsional radius²]

CSCR optimises:
    max  cos⁴(x−b) · (c − r·cos(2x)) / (c + r·cos(2x))   [e/l]
"""

from __future__ import annotations

from collections.abc import Callable, Iterator

import numpy as np


def mohr_params(I1: np.ndarray, I2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Mohr's circle centre *c* and radius *r* from principal moments.

    Args:
        I1: Larger principal moments (N × 1 array).
        I2: Smaller principal moments (N × 1 array).

    Returns:
        ``(c, r)`` arrays.
    """
    c = 0.5 * (I1 + I2)
    r = 0.5 * (I1 - I2)
    return c, r


def _signed_angle(dir1: np.ndarray, e_vec: np.ndarray) -> np.ndarray:
    """Signed angle (radians) from *dir1* to *e_vec*, one building per row.

    Plain vectorised ``arctan2`` over the whole array -- the earlier version
    of this built the same cross/dot products but fed them into
    ``arctan2`` one building at a time inside a Python list comprehension,
    which is pure overhead: cross/dot/arctan2 are already elementwise
    NumPy ufuncs, so there's nothing a per-row loop buys here.
    """
    cross = dir1[:, 0] * e_vec[:, 1] - dir1[:, 1] * e_vec[:, 0]
    dot = np.einsum("ij,ij->i", dir1, e_vec)
    return np.arctan2(cross, dot)


def _golden_section_max(
    f: Callable[[np.ndarray], np.ndarray],
    lo: np.ndarray,
    hi: np.ndarray,
    iters: int = 60,
) -> np.ndarray:
    """Vectorised golden-section search for the maximiser of *f* on ``[lo, hi]``.

    ``f`` is called with an ``(N,)`` array of candidate angles and must
    return an ``(N,)`` array of objective values (elementwise NumPy ops
    only -- no Python-level per-row branching). Every building is refined
    in lockstep: each of the ``iters`` rounds is exactly 1-2 vectorised
    function evaluations over *all* buildings at once, instead of the
    previous ``scipy.optimize.fmin`` (Nelder-Mead) called once per
    building inside a Python loop, which paid full per-call Python/SciPy
    dispatch overhead N times over for what is, per building, a single
    smooth 1-D optimisation. Assumes ``f`` is unimodal on ``[lo, hi]``,
    which holds here because ``lo``/``hi`` bracket one cell of the coarse
    grid search that seeds them (see :func:`optimise_ec8`/
    :func:`optimise_cscr`) -- narrow enough that the objective's known
    smooth, low-frequency shape can't have a second local max inside it.

    Returns:
        ``(N,)`` array of angles maximising ``f`` within ``[lo, hi]``.
    """
    invphi = (np.sqrt(5.0) - 1.0) / 2.0  # ~0.618
    a, b = lo.astype(float).copy(), hi.astype(float).copy()
    for _ in range(iters):
        c = b - invphi * (b - a)
        d = a + invphi * (b - a)
        take_left = f(c) > f(d)  # maximiser lies in [a, d]
        a = np.where(take_left, a, c)
        b = np.where(take_left, d, b)
    return 0.5 * (a + b)


def _maximize_periodic(
    f: Callable[[np.ndarray], np.ndarray],
    n: int,
    n_grid: int = 180,
    refine_iters: int = 40,
) -> np.ndarray:
    """Maximise a period-π objective ``f`` for every one of *n* buildings at once.

    Coarse vectorised grid search (all buildings, all grid angles, in one
    ``(n_grid,) x (n,)`` broadcast -- no per-building work) picks which
    grid cell each building's maximum falls in, then
    :func:`_golden_section_max` polishes every building's estimate inside
    its own cell in lockstep. Both stages are pure array ops, so the whole
    search across an entire GeoDataFrame costs a fixed, small number of
    vectorised objective evaluations regardless of *n* -- unlike the
    previous implementation's ``scipy.optimize.fmin`` call issued
    separately for every building inside a Python ``for`` loop, where
    per-call SciPy/Python dispatch overhead (not the trivial trig math
    itself) dominated runtime for any dataset of realistic size.

    Args:
        f: Objective; called with an ``(n_grid, n)``-broadcastable array
            of angles and must return values of the same shape.
        n: Number of buildings (rows).
        n_grid: Number of coarse grid angles spanning one period (π).
        refine_iters: Golden-section iterations to polish the grid winner.

    Returns:
        ``(n,)`` array of angles (radians) maximising ``f``.
    """
    grid = np.linspace(0.0, np.pi, n_grid, endpoint=False)
    values = f(grid[:, None])  # (n_grid, n)
    best_idx = np.argmax(values, axis=0)
    x0 = grid[best_idx]
    step = np.pi / n_grid
    lo, hi = x0 - step, x0 + step
    return _golden_section_max(f, lo, hi, iters=refine_iters)


#: Row-chunk size the grid search in :func:`optimise_ec8`/:func:`optimise_cscr`
#: is capped to. The grid step builds one ``(180, chunk)`` array (plus a
#: handful of same-shaped temporaries) per chunk instead of per whole
#: dataset -- at the default this caps that step to roughly the size of the
#: 911-building San Jose pilot dataset's arrays x ~20, i.e. a few hundred
#: MB regardless of how many buildings are passed in, instead of scaling
#: linearly with N (~440MB measured at N=100,000 unchunked). Large enough
#: that ordinary-sized datasets (well under this) still run as a single
#: chunk -- one vectorised pass, same speed as before -- so this only
#: kicks in for datasets actually big enough to need it.
_DEFAULT_CHUNK_SIZE = 10_000


def _chunk_slices(n: int, chunk_size: int) -> Iterator[slice]:
    """Yield consecutive ``slice(start, end)`` row-chunks covering ``range(n)``.

    Args:
        n: Total number of rows to cover.
        chunk_size: Rows per chunk (the final chunk may be shorter).

    Yields:
        One ``slice`` per chunk, in order, together covering ``[0, n)``.
    """
    for start in range(0, n, chunk_size):
        yield slice(start, min(start + chunk_size, n))


def optimise_ec8(
    I1: np.ndarray,
    dir1: np.ndarray,
    I2: np.ndarray,
    e_vec: np.ndarray,
    area: np.ndarray,
    chunk_size: int = _DEFAULT_CHUNK_SIZE,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Find the worst-case EC8 eccentricity ratio for each building.

    EC8 defines::

        r_t^2(x) = I_t / I_j(x)
        eccentricity_ratio = e * abs(cos(x - b)) / r_t(x)

    where I_t = I_0 + A · e² is the torsional inertia about the centre of
    stiffness, and I_j(x) = c − r·cos(2x) is the moment in direction x.

    Args:
        I1:    (N,) array of larger principal moments.
        dir1:  (N,2) array of unit eigenvectors for I1.
        I2:    (N,) array of smaller principal moments.
        e_vec: (N,2) eccentricity vectors (CM − CS).
        area:  (N,) footprint areas.
        chunk_size: Buildings processed per vectorised grid-search batch
            -- bounds peak memory (see :data:`_DEFAULT_CHUNK_SIZE`)
            independent of the total number of buildings.

    Returns:
        ``(ecc_ratio, radius_ratio, x_opt, b)`` all 1-D arrays of length N.
    """
    n = len(I1)
    if n > chunk_size:
        parts = [
            optimise_ec8(
                I1[s], dir1[s], I2[s], e_vec[s], area[s], chunk_size=chunk_size
            )
            for s in _chunk_slices(n, chunk_size)
        ]
        return tuple(np.concatenate(arrs) for arrs in zip(*parts))

    c, r = mohr_params(I1, I2)
    I0 = I1 + I2
    e_mag = np.linalg.norm(e_vec, axis=1)
    has_ecc = e_mag >= 1e-10

    b = np.where(has_ecc, _signed_angle(dir1, e_vec), 0.0)

    # I_t: torsional inertia = I0 + A * e²
    I_t = I0 + area * e_mag**2

    # EC8 wants the worst ratio, i.e. the x that *maximises*
    # cos²(x-b) * (c - r*cos(2x)) (this minimises the torsional radius that
    # ratio is divided by). Solved for every building at once; skip
    # buildings with no eccentricity (x_opt=0 trivially, same as fmin's
    # unconverged x0=0.0 there previously).
    def _objective(x: np.ndarray) -> np.ndarray:
        """EC8's torsional-radius objective, maximised over analysis angle ``x``."""
        return (np.cos(x - b) ** 2) * (c - r * np.cos(-2.0 * x))

    x_opt = np.where(has_ecc, _maximize_periodic(_objective, len(I1)), 0.0)

    I_j = c - r * np.cos(-2.0 * x_opt)  # moment in worst direction
    r_t = np.sqrt(I_t / (I_j + 1e-30))  # torsional radius
    r_g = np.sqrt(I0 / (area + 1e-30))  # radius of gyration

    ecc_ratio = e_mag * np.abs(np.cos(x_opt - b)) / (r_t + 1e-30)
    rad_ratio = r_t / (r_g + 1e-30)

    return ecc_ratio, rad_ratio, x_opt, b


def optimise_cscr(
    I1: np.ndarray,
    dir1: np.ndarray,
    I2: np.ndarray,
    e_vec: np.ndarray,
    area: np.ndarray,
    chunk_size: int = _DEFAULT_CHUNK_SIZE,
) -> tuple[np.ndarray, np.ndarray]:
    """Find the worst-case CSCR 2010 eccentricity ratio for each building.

    CSCR uses e / l where l is proportional to the building dimension in the
    analysis direction:

        l(x) = sqrt(area) · ((c + r·cos(2x)) / (c − r·cos(2x)))^0.25

    Args:
        I1, dir1, I2: As for :func:`optimise_ec8`.
        e_vec: (N,2) eccentricity vectors (CM − CS_polygon_centroid).
        area:  (N,) footprint areas.
        chunk_size: See :func:`optimise_ec8`.

    Returns:
        ``(ecc_ratio, x_opt)`` both 1-D arrays of length N.
    """
    n = len(I1)
    if n > chunk_size:
        parts = [
            optimise_cscr(
                I1[s], dir1[s], I2[s], e_vec[s], area[s], chunk_size=chunk_size
            )
            for s in _chunk_slices(n, chunk_size)
        ]
        return tuple(np.concatenate(arrs) for arrs in zip(*parts))

    c, r = mohr_params(I1, I2)
    e_mag = np.linalg.norm(e_vec, axis=1)
    has_ecc = e_mag >= 1e-10

    b = np.where(has_ecc, _signed_angle(dir1, e_vec), 0.0)

    def _objective(x: np.ndarray) -> np.ndarray:
        """CSCR 2010's e/l objective, maximised over analysis angle ``x``."""
        Ij_max = c + r * np.cos(-2.0 * x)
        Ij_min = c - r * np.cos(-2.0 * x)
        # Guard the pole at Ij_min == 0 the same way the old per-row
        # optimiser did (treat it as a non-improving, zero-valued point
        # rather than propagating inf/nan into the grid/golden-section
        # search); Ij_min is c ∓ r, always >= 0 for I1 >= I2, so this only
        # triggers exactly at the (measure-zero) degenerate angle.
        # `np.where` evaluates both branches eagerly, so the division is
        # still computed (and would warn) at those points even though its
        # result is discarded -- silence just that expected warning rather
        # than the whole numpy error-state, so an unrelated real
        # divide-by-zero elsewhere still surfaces normally.
        safe_min = np.where(np.abs(Ij_min) < 1e-30, 1.0, Ij_min)
        return np.where(
            np.abs(Ij_min) < 1e-30, 0.0, np.cos(x - b) ** 4 * Ij_max / safe_min
        )

    x_opt = np.where(has_ecc, _maximize_periodic(_objective, len(I1)), 0.0)

    e_proj = np.abs(e_mag * np.cos(x_opt - b))
    I_long = c + r * np.cos(-2.0 * x_opt)  # moment along long side
    I_short = c - r * np.cos(-2.0 * x_opt)  # moment along short side
    l_dim = np.sqrt(area + 1e-30) * (I_long / (I_short + 1e-30)) ** 0.25
    ecc_ratio = e_proj / (l_dim + 1e-30)

    return ecc_ratio, x_opt
