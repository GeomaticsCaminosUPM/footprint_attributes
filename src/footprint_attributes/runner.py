"""Single entry-point orchestration: ``footprint_attributes.run(...)``.

Wires together :mod:`direction`, :mod:`shape`, and :mod:`position` behind
one call so a caller can request columns by name (or by norm) without
having to know which submodule produces them.
"""

from __future__ import annotations

import geopandas as gpd

from .geometry import ensure_projected, to_gdf
from .shape import shape as _shape, EC8, ASCE7, GNDTII, CSCR2010, NTC23, slenderness
from .position import position as _position
from .direction import inertia as _direction_inertia

#: Norm-name shorthand -> NormAggregate instance (for expanding e.g. "EC8"
#: into all of that norm's own columns, including its compliance columns).
_NORMS_BY_NAME = {
    "EC8": EC8,
    "ASCE7": ASCE7,
    "GNDTII": GNDTII,
    "CSCR2010": CSCR2010,
    "NTC23": NTC23,
}

#: Code-independent shape indices, requestable by their plain column name.
_CODE_INDEPENDENT_SHAPE_COLUMNS = {
    "polsby_popper",
    "convex_hull_irregularity",
    "inertia_circle_ratio",
}

#: Columns produced by the position pipeline, as a single atomic group
#: (position.__call__ always computes all of them together).
_POSITION_COLUMNS = [
    "contact_force",
    "contact_confinementRatio",
    "contact_angularAcc",
    "contact_angle",
    "contact_height",
    "relativePosition",
]

#: Keyword arguments accepted by position()/shape() that `run()` will NOT
#: forward from config, because they select the bbox-vs-inertia axis
#: convention -- `run()` always uses each parameter's own package default
#: for that choice (see the "method" note in the docstring below).
_METHOD_KEYS = {"method"}

#: Group/shorthand names accepted in ``columns`` besides exact shape.py
#: column names (norm names, "slenderness", "position", its five exact
#: column names, and "bearing").
_GROUP_NAMES = (
    set(_NORMS_BY_NAME)
    | {"slenderness", "position", "bearing"}
    | set(_POSITION_COLUMNS)
)

#: Every exact `shape.py` column name this package can produce, used to
#: validate columns that don't match a group name above.
_ALL_SHAPE_COLUMNS: set[str] = set(_CODE_INDEPENDENT_SHAPE_COLUMNS)
for _norm in _NORMS_BY_NAME.values():
    for _param in _norm.parameters.values():
        _ALL_SHAPE_COLUMNS.add(_param.column_name)
        if _param.limits is not None:
            _ALL_SHAPE_COLUMNS.add(f"compliance_{_param.column_name}")
for _method_obj in slenderness._methods.values():
    _ALL_SHAPE_COLUMNS.add(_method_obj.column_name)
    if _method_obj._compliance_limits is not None:
        _ALL_SHAPE_COLUMNS.add(f"compliance_EC8_{_method_obj.column_name}")
del _norm, _param, _method_obj

#: Every valid ``columns`` entry: group names plus exact shape columns.
_ALL_VALID_COLUMNS = _GROUP_NAMES | _ALL_SHAPE_COLUMNS


def _expand_columns(columns: list[str]) -> tuple[set[str], bool, bool]:
    """Resolve the requested ``columns`` into (shape_columns, want_position, want_bearing).

    Args:
        columns: Raw column/group names as passed in ``config["columns"]``.

    Returns:
        ``(shape_columns, want_position, want_bearing)``: the resolved set
        of exact `shape.py` column names to compute, whether the position
        pipeline was requested, and whether the direction bearing was
        requested.
    """
    shape_columns: set[str] = set()
    want_position = False
    want_bearing = False

    for col in columns:
        if col in _NORMS_BY_NAME:
            norm = _NORMS_BY_NAME[col]
            for param in norm.parameters.values():
                shape_columns.add(param.column_name)
                if param.limits is not None:
                    shape_columns.add(f"compliance_{param.column_name}")
        elif col == "slenderness":
            for method_obj in slenderness._methods.values():
                shape_columns.add(method_obj.column_name)
                if method_obj._compliance_limits is not None:
                    shape_columns.add(f"compliance_EC8_{method_obj.column_name}")
        elif col == "position":
            want_position = True
        elif col in _POSITION_COLUMNS:
            want_position = True
        elif col == "bearing":
            want_bearing = True
        elif col in _ALL_SHAPE_COLUMNS:
            # Exact shape.py column name (e.g. "EC8_eccentricityRatio",
            # "compliance_ASCE7_setbackRatio", "slenderness_bbox",
            # "polsby_popper", ...).
            shape_columns.add(col)
        else:
            print(
                f"WARNING: run(): unrecognised column/code name {col!r} -- "
                f"skipping. Valid names are norm names ({sorted(_NORMS_BY_NAME)}), "
                "'slenderness', 'position', 'bearing', or an exact "
                "shape.py column name."
            )

    return shape_columns, want_position, want_bearing


def run(
    footprints_gdf: gpd.GeoDataFrame,
    config: dict | None = None,
    overwrite: bool = False,
) -> gpd.GeoDataFrame:
    """Compute a requested set of footprint attributes in one call.

    This is the single entry point that wires together
    :mod:`footprint_attributes.direction`, :mod:`footprint_attributes.shape`,
    and :mod:`footprint_attributes.position`: tell it which columns you want
    (exact column names, or a norm name as shorthand for "all of that
    norm's columns"), optionally override the underlying functions'
    parameters, and get back the original GeoDataFrame with those columns
    added.

    Args:
        footprints_gdf: GeoDataFrame of building footprint polygons (any
            CRS; reprojected to a projected CRS internally if needed).
        config: Dictionary describing what to compute. See "Config dict
            shape" below. If ``None`` or ``{"columns": []}``, no columns
            are computed and the input is returned unchanged (aside from
            CRS/index normalisation).
        overwrite: If ``False`` (default), any requested column that
            already exists in ``footprints_gdf`` is left untouched (not
            recomputed) -- this makes ``run()`` safe to call repeatedly on
            a GeoDataFrame you're incrementally building up. If ``True``,
            every requested column is recomputed and replaced, even if
            already present.

    Config dict shape
    ------------------
    ``config`` has one required key and up to three optional
    per-module keyword-argument sub-dictionaries::

        config = {
            "columns": [...],       # required: what to compute (see below)
            "position": {...},      # optional: kwargs forwarded to position()
            "shape": {...},         # optional: kwargs forwarded to shape()
            "direction": {...},     # optional: kwargs forwarded to direction.inertia()
        }

    ``columns`` (required, list of str)
        Each entry is one of:

        - An **exact `shape.py` column name**, e.g. ``"EC8_eccentricityRatio"``,
          ``"compliance_ASCE7_setbackRatio"``, ``"slenderness_bbox"``,
          ``"polsby_popper"``, ``"convex_hull_irregularity"``,
          ``"inertia_circle_ratio"``.
        - A **norm name** as shorthand for all of that norm's own columns
          (including its compliance columns): ``"EC8"``, ``"ASCE7"``,
          ``"GNDTII"``, ``"CSCR2010"``, or ``"NTC23"``.
        - ``"slenderness"`` -- shorthand for both `slenderness_bbox` and
          `slenderness_inertia` (plus their EC8 compliance columns).
        - ``"position"`` -- shorthand for the full position pipeline's
          five columns: ``contact_force``, ``contact_confinementRatio``,
          ``contact_angularAcc``, ``contact_angle``, ``relativePosition``.
          Any of these five exact names also works on their own and
          triggers the same (atomic -- they're always computed together)
          pipeline run.
        - ``"bearing"`` -- the building direction bearing (degrees from
          geographic North), via ``direction.inertia()``.

    ``config["position"]`` (optional dict)
        Forwarded as keyword arguments to ``position()``. Accepted keys
        (see ``config.POSITION_DEFAULTS`` for the package defaults):
        ``buffer`` (contact-detection buffer, metres), ``height_column``
        (name of a column with building heights), ``minRadius``,
        ``minForce``, ``minAngle``, ``minConfinement``, ``minAngularAcc``
        (classification thresholds).

    ``config["shape"]`` (optional dict)
        Forwarded as keyword arguments to ``shape()``. The main use is
        ``height_column`` (for `slenderness`'s ``vertical=True`` mode --
        pass ``{"height_column": "height", "vertical": True}`` if you also
        request vertical slenderness explicitly via
        ``shape.slenderness.bbox(gdf, vertical=True, ...)`` outside of
        `run()`, since vertical slenderness is not part of the
        `"slenderness"` column shorthand above).

        **Note:** ``method`` (``"bbox"`` vs. ``"inertia"``, i.e. which
        principal-axis convention underlies setback/GNDT constructions) is
        **not** configurable through `run()` -- each parameter always uses
        its own package default (`"bbox"` for setback/GNDT parameters,
        `"inertia"`-only for eccentricity parameters). Any `"method"` key
        found in `config["shape"]` is silently ignored. If you need a
        specific method, call `shape.<Norm>.<param>(gdf, method=...)`
        directly instead of `run()`.

    ``config["direction"]`` (optional dict)
        Forwarded as keyword arguments to ``direction.inertia()`` when
        ``"bearing"`` is requested (e.g. nothing meaningful to override
        today besides ``direction=`` to force a specific axis; kept for
        forward-compatibility).

    Returns:
        A copy of ``footprints_gdf`` with the requested columns added (or
        left as-is where ``overwrite=False`` and the column already
        existed).

    Example:
        >>> import footprint_attributes
        >>> result = footprint_attributes.run(
        ...     footprints,
        ...     config={
        ...         "columns": ["EC8", "position", "bearing"],
        ...         "position": {"buffer": 0.1},
        ...     },
        ... )
    """
    config = config or {}
    columns = config.get("columns") or []
    position_kwargs = dict(config.get("position", {}))
    shape_kwargs = dict(config.get("shape", {}))
    direction_kwargs = dict(config.get("direction", {}))
    for key in _METHOD_KEYS:
        shape_kwargs.pop(key, None)

    gdf = ensure_projected(to_gdf(footprints_gdf))
    result = gdf.copy()

    shape_columns, want_position, want_bearing = _expand_columns(columns)

    if shape_columns:
        working = result.copy()
        if overwrite:
            working = working.drop(
                columns=[c for c in shape_columns if c in working.columns]
            )
            cols_to_compute = shape_columns
        else:
            cols_to_compute = {c for c in shape_columns if c not in result.columns}

        if cols_to_compute:
            computed = _shape(working, columns=list(cols_to_compute), **shape_kwargs)
            for col in cols_to_compute:
                if col in computed.columns:
                    result[col] = computed[col].values

    if want_position:
        already_present = all(c in result.columns for c in _POSITION_COLUMNS)
        if overwrite or not already_present:
            pos_result = _position(result, **position_kwargs)
            for col in _POSITION_COLUMNS:
                result[col] = pos_result[col].values

    if want_bearing:
        if overwrite or "bearing" not in result.columns:
            result["bearing"] = _direction_inertia(result, **direction_kwargs)

    return result
