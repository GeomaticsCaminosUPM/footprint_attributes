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
import numpy as np
from scipy.optimize import fmin


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


def optimise_ec8(
    I1: np.ndarray,
    dir1: np.ndarray,
    I2: np.ndarray,
    e_vec: np.ndarray,
    area: np.ndarray,
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

    Returns:
        ``(ecc_ratio, radius_ratio, x_opt, b)`` all 1-D arrays of length N.
    """
    c, r = mohr_params(I1, I2)
    I0 = I1 + I2
    e_mag = np.linalg.norm(e_vec, axis=1)

    # Angle b: signed angle from dir1 to eccentricity vector
    b = np.array(
        [
            0.0
            if e_mag[i] < 1e-10
            else float(
                np.arctan2(
                    dir1[i, 0] * e_vec[i, 1] - dir1[i, 1] * e_vec[i, 0],  # cross
                    np.dot(dir1[i], e_vec[i]),  # dot
                )
            )
            for i in range(len(I1))
        ]
    )

    # I_t: torsional inertia = I0 + A * e²
    I_t = I0 + area * e_mag**2

    x_opt = np.zeros(len(I1))
    for i in range(len(I1)):
        if e_mag[i] < 1e-10:
            continue

        def _neg_objective(x, _c=c[i], _r=r[i], _b=b[i]):
            # We maximise (= minimise negative) the torsional radius
            # which minimises the eccentricity ratio.
            # EC8 wants the worst ratio, so we minimise -(cos²(x-b) * I_j(x))
            return -(np.cos(x - _b) ** 2) * (_c - _r * np.cos(-2.0 * x))

        x_opt[i] = fmin(_neg_objective, x0=0.0, xtol=1e-5, ftol=1e-5, disp=False)[0]

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
) -> tuple[np.ndarray, np.ndarray]:
    """Find the worst-case CSCR 2010 eccentricity ratio for each building.

    CSCR uses e / l where l is proportional to the building dimension in the
    analysis direction:

        l(x) = sqrt(area) · ((c + r·cos(2x)) / (c − r·cos(2x)))^0.25

    Args:
        I1, dir1, I2: As for :func:`optimise_ec8`.
        e_vec: (N,2) eccentricity vectors (CM − CS_polygon_centroid).
        area:  (N,) footprint areas.

    Returns:
        ``(ecc_ratio, x_opt)`` both 1-D arrays of length N.
    """
    c, r = mohr_params(I1, I2)
    e_mag = np.linalg.norm(e_vec, axis=1)

    b = np.array(
        [
            0.0
            if e_mag[i] < 1e-10
            else float(
                np.arctan2(
                    dir1[i, 0] * e_vec[i, 1] - dir1[i, 1] * e_vec[i, 0],
                    np.dot(dir1[i], e_vec[i]),
                )
            )
            for i in range(len(I1))
        ]
    )

    x_opt = np.zeros(len(I1))
    for i in range(len(I1)):
        if e_mag[i] < 1e-10:
            continue

        def _neg_obj(x, _c=c[i], _r=r[i], _b=b[i]):
            Ij_max = _c + _r * np.cos(-2.0 * x)
            Ij_min = _c - _r * np.cos(-2.0 * x)
            # Avoid division by zero
            if abs(Ij_min) < 1e-30:
                return 0.0
            return -(np.cos(x - _b) ** 4 * Ij_max / Ij_min)

        x_opt[i] = fmin(_neg_obj, x0=0.0, xtol=1e-5, ftol=1e-5, disp=False)[0]

    e_proj = np.abs(e_mag * np.cos(x_opt - b))
    I_long = c + r * np.cos(-2.0 * x_opt)  # moment along long side
    I_short = c - r * np.cos(-2.0 * x_opt)  # moment along short side
    l_dim = np.sqrt(area + 1e-30) * (I_long / (I_short + 1e-30)) ** 0.25
    ecc_ratio = e_proj / (l_dim + 1e-30)

    return ecc_ratio, x_opt
