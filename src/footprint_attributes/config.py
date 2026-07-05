"""
Default configuration: code limits, compliance grades, and thresholds.

All default values used throughout the package are defined here so they
can be inspected and overridden from a single place.

Compliance grade structure
--------------------------
Each norm parameter is described by an ordered list of threshold entries::

    [
        {"max": <value>,  "label": "<text>",  "score": <0–100>},
        ...
        {"max": inf,      "label": "<worst>", "score": 0},
    ]

To evaluate a parameter *x*:
  - Walk the list until ``x <= entry["max"]``.
  - The matching ``score`` (0–100) and ``label`` are assigned.
  - score=100 means the code limit is satisfied;
    score=0 means the worst violation.

The column name pattern is  ``{NORM}_{code}_{name}``  where ``code`` is the
parameter identifier within the norm (e.g. "beta2") and ``name`` is the
parameter's camelCase name.  The compliance column is
``compliance_{NORM}_{code}_{name}``.
"""

from __future__ import annotations
import math


# ─────────────────────────────────────────────────────────────────────────────
# Eurocode 8 (EN 1998-1)
# ─────────────────────────────────────────────────────────────────────────────

EC8_LIMITS: dict[str, list[dict]] = {
    # eccentricity_ratio  e / r_t  ≤ 0.30
    "eccentricityRatio": [
        {"max": 0.30, "label": "regular", "score": 100},
        {"max": math.inf, "label": "irregular", "score": 0},
    ],
    # radius_ratio  r_t / r_g  ≥ 1.0  (inverted: non-compliance if < 1)
    "radiusRatio": [
        {"max": 1.00, "label": "irregular", "score": 0},  # < 1 → fail
        {"max": math.inf, "label": "regular", "score": 100},
    ],
    # compactness  ≥ 0.95  (stored as 1 − max_setback/area)
    "compactness": [
        {"max": 0.95, "label": "irregular", "score": 0},
        {"max": math.inf, "label": "regular", "score": 100},
    ],
}

EC8_WALL_HEIGHT: float = 3.0  # metres assumed per storey


# ─────────────────────────────────────────────────────────────────────────────
# Slenderness (used by multiple codes, treated separately)
# ─────────────────────────────────────────────────────────────────────────────

SLENDERNESS_LIMITS: dict[str, list[dict]] = {
    # EC8 plan slenderness ≤ 4.0
    "EC8": [
        {"max": 4.0, "label": "regular", "score": 100},
        {"max": math.inf, "label": "irregular", "score": 0},
    ],
    # GNDTII β₃ (concrete) ≥ 0.4  i.e. L2/L1 ≥ 0.4  → slenderness = L1/L2 ≤ 2.5
    "GNDTII_beta3": [
        {"max": 2.5, "label": "A", "score": 100},
        {"max": 5.0, "label": "B", "score": 50},
        {"max": math.inf, "label": "C", "score": 0},
    ],
}


# ─────────────────────────────────────────────────────────────────────────────
# Costa Rica Seismic Code (CSCR 2010)
# ─────────────────────────────────────────────────────────────────────────────

CSCR2010_LIMITS: dict[str, list[dict]] = {
    # eccentricity_ratio  e / l
    "eccentricityRatio": [
        {"max": 0.05, "label": "regular", "score": 100},
        {"max": 0.25, "label": "moderate", "score": 50},
        {"max": math.inf, "label": "high", "score": 0},
    ],
}


# ─────────────────────────────────────────────────────────────────────────────
# Italian GNDT Level II
# ─────────────────────────────────────────────────────────────────────────────

GNDTII_LIMITS: dict[str, list[dict]] = {
    # β₁  a / L  (masonry)  — higher is better
    "beta1_mainShapeSlenderness": [
        {"max": 0.4, "label": "D", "score": 0},
        {"max": 0.6, "label": "C", "score": 33},
        {"max": 0.8, "label": "B", "score": 66},
        {"max": math.inf, "label": "A", "score": 100},
    ],
    # β₂  b / L  (masonry setback ratio)  — lower is better
    "beta2_setbackRatio": [
        {"max": 0.1, "label": "A", "score": 100},
        {"max": 0.2, "label": "B", "score": 66},
        {"max": 0.3, "label": "C", "score": 33},
        {"max": math.inf, "label": "D", "score": 0},
    ],
    # β₄  e / a  (concrete eccentricity)  — lower is better
    "beta4_eccentricityRatio": [
        {"max": 0.2, "label": "A", "score": 100},
        {"max": 0.4, "label": "B", "score": 50},
        {"max": math.inf, "label": "C", "score": 0},
    ],
    # β₆  c / b  (setback slenderness)  — higher is better
    "beta6_setbackSlenderness": [
        {"max": 0.25, "label": "C", "score": 0},
        {"max": 0.5, "label": "B", "score": 50},
        {"max": math.inf, "label": "A", "score": 100},
    ],
}


# ─────────────────────────────────────────────────────────────────────────────
# ASCE 7
# ─────────────────────────────────────────────────────────────────────────────

ASCE7_LIMITS: dict[str, list[dict]] = {
    # setback_ratio  b / L  ≤ 0.20
    "setbackRatio": [
        {"max": 0.20, "label": "regular", "score": 100},
        {"max": math.inf, "label": "irregular", "score": 0},
    ],
    # hole_ratio  A_hole / A_filled  ≤ 0.25
    "holeRatio": [
        {"max": 0.25, "label": "regular", "score": 100},
        {"max": math.inf, "label": "irregular", "score": 0},
    ],
    # parallelity_angle  — no quantitative limit; flag > 5° as non-rectangular
    "parallelityAngle": [
        {"max": 5.0, "label": "regular", "score": 100},
        {"max": 10.0, "label": "skewed", "score": 50},
        {"max": math.inf, "label": "triangular", "score": 0},
    ],
}


# ─────────────────────────────────────────────────────────────────────────────
# Mexican NTC-23
# ─────────────────────────────────────────────────────────────────────────────

NTC23_LIMITS: dict[str, list[dict]] = {
    # setback_ratio  b / L  ≤ 0.40
    "setbackRatio": [
        {"max": 0.40, "label": "regular", "score": 100},
        {"max": math.inf, "label": "irregular", "score": 0},
    ],
    # hole_ratio  h / L  ≤ 0.40
    "holeRatio": [
        {"max": 0.40, "label": "regular", "score": 100},
        {"max": math.inf, "label": "irregular", "score": 0},
    ],
}


# ─────────────────────────────────────────────────────────────────────────────
# Position thresholds
# ─────────────────────────────────────────────────────────────────────────────

POSITION_DEFAULTS: dict = {
    # Minimum resultant force / sqrt(area) to count as "lateral".
    # Default: a square building with 1/6 of one side touching → force = 1/6.
    "minForce": 0.166,
    # Minimum weighted angle (rad) between individual forces and resultant
    # to upgrade lateral → corner. Default: π/4 ≈ 45°.
    "minAngle": 0.78,
    # Minimum confinement ratio to classify as "confined". Default: 0.5.
    "minConfinement": 0.5,
    # Minimum angular-acceleration proxy to upgrade corner/confined → torque.
    # Default: a 1×0.5 rectangle with 1/3 of two sides touching in worst case.
    "minAngularAcc": 2.133,
    # Contact detection buffer (metres).
    "buffer": 0.0,
    # Minimum radius fraction for conservative momentum counting.
    "minRadius": 0.0,
}


# ─────────────────────────────────────────────────────────────────────────────
# Direction method defaults
# ─────────────────────────────────────────────────────────────────────────────

DIRECTION_DEFAULTS: dict = {
    "method": "inertia",  # or "bbox"
    "mode": "bearing",  # or "dimensions", "directions", "all"
}


# ─────────────────────────────────────────────────────────────────────────────
# Helper: evaluate a parameter value against a compliance table
# ─────────────────────────────────────────────────────────────────────────────


def compliance_score(value: float, table: list[dict]) -> tuple[int, str]:
    """Look up a parameter value in a compliance table.

    Args:
        value: The parameter value to evaluate.
        table: Ordered list of ``{"max": ..., "label": ..., "score": ...}``
            dicts.  Must be sorted by ``max`` ascending.

    Returns:
        Tuple ``(score, label)`` for the first entry whose ``max >= value``.
    """
    for entry in table:
        if value <= entry["max"]:
            return entry["score"], entry["label"]
    # Should not reach here if table ends with inf
    last = table[-1]
    return last["score"], last["label"]
