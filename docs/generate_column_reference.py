"""Regenerate the "Complete column reference" section of ``docs/api/run.rst``.

Introspects the package's own norm/parameter registries (the same ones
``footprint_attributes.run()`` validates ``columns`` against) so the table
can never drift out of sync with what the code actually produces. Run this
and paste its output into ``docs/api/run.rst`` (replacing everything from
the ``Complete column reference`` heading onward) whenever a column,
parameter, or norm is added, renamed, or removed.

Usage::

    uv run python docs/generate_column_reference.py > /tmp/run_columns.rst
    # then copy /tmp/run_columns.rst into docs/api/run.rst
"""

from __future__ import annotations

from footprint_attributes import runner
from footprint_attributes.shape import slenderness

_CODE_INDEPENDENT_DESCRIPTIONS = {
    "polsby_popper": "4π A / P² -- compactness, 1.0 = perfect circle",
    "convex_hull_irregularity": (
        "(hull_area - footprint_area) / footprint_area -- 0 = convex"
    ),
    "inertia_circle_ratio": "I_z(equal-area circle) / I_z(footprint)",
}

_POSITION_DESCRIPTIONS = {
    "contact_force": "Magnitude of the net contact-force resultant from touching neighbours",
    "contact_confinementRatio": "How evenly touching neighbours surround the building (0-1)",
    "contact_angularAcc": (
        "Angular acceleration from an uneven (off-centre) contact-force distribution"
    ),
    "contact_angle": "Direction of the net contact-force resultant",
    "contact_height": "Height used for the contact-force computation",
    "relativePosition": "Categorical class: isolated / lateral / corner / confined / torque",
}


def _param_table(rows: list[tuple[str, str, str]]) -> list[str]:
    """Render ``(column, compliance, description)`` rows as an RST list-table."""
    lines = [
        ".. list-table::",
        "   :header-rows: 1",
        "   :widths: 30 10 60",
        "",
        "   * - Column",
        "     - Compliance",
        "     - Description",
    ]
    for column, compliance, description in rows:
        lines.append(f"   * - ``{column}``")
        lines.append(f"     - {compliance}")
        lines.append(f"     - {description}")
    return lines


def _plain_table(rows: list[tuple[str, str]]) -> list[str]:
    """Render ``(column, description)`` rows as an RST list-table (no compliance column)."""
    lines = [
        ".. list-table::",
        "   :header-rows: 1",
        "   :widths: 30 70",
        "",
        "   * - Column",
        "     - Description",
    ]
    for column, description in rows:
        lines.append(f"   * - ``{column}``")
        lines.append(f"     - {description}")
    return lines


def build_rst() -> str:
    """Render the full "Complete column reference" RST section."""
    lines = [
        "Complete column reference",
        "==========================",
        "",
        "Every column :func:`~footprint_attributes.runner.run` (or"
        " :func:`~footprint_attributes.shape.shape`) can produce, grouped by"
        " seismic code. Each parameter also gets a matching"
        " ``compliance_{column}`` column (0-100 score against the code's own"
        " limit) unless noted otherwise.",
        "",
    ]

    for name, norm in sorted(runner._NORMS_BY_NAME.items()):
        rows = [
            (
                param.column_name,
                "yes" if param.limits is not None else "no",
                param.description,
            )
            for param in norm.parameters.values()
        ]
        lines += [name, "-" * len(name), ""]
        lines += _param_table(rows)
        lines += [
            "",
            f'Request all of {name}\'s columns at once with ``"{name}"`` in'
            f' ``config["columns"]``, or ``shape.{name}(gdf)``.',
            "",
        ]

    slenderness_rows = [
        (
            m.column_name,
            "yes" if m._compliance_limits is not None else "no",
            f"Plan slenderness via the {mname} axis convention (EC8 limit)",
        )
        for mname, m in slenderness._methods.items()
    ]
    lines += ["Plan/vertical slenderness", "--------------------------", ""]
    lines += _param_table(slenderness_rows)
    lines += [
        "",
        'Request both with ``"slenderness"``, or ``shape.slenderness(gdf)``.',
        "",
    ]

    code_independent_rows = [
        (c, _CODE_INDEPENDENT_DESCRIPTIONS[c])
        for c in sorted(runner._CODE_INDEPENDENT_SHAPE_COLUMNS)
    ]
    lines += ["Code-independent shape indices", "-------------------------------", ""]
    lines += _plain_table(code_independent_rows)
    lines += [
        "",
        "No compliance columns (these are not tied to a specific code's limit).",
        "",
    ]

    position_rows = [(c, _POSITION_DESCRIPTIONS[c]) for c in runner._POSITION_COLUMNS]
    lines += ["Position pipeline", "------------------", ""]
    lines += [
        'Requesting ``"position"`` (or any one of the five columns below,'
        " which always triggers the same atomic pipeline run) adds:",
        "",
    ]
    lines += _plain_table(position_rows)
    lines += [""]

    lines += ["Direction", "---------", ""]
    lines += _plain_table(
        [
            (
                "bearing",
                "Building orientation, degrees clockwise from North, range"
                " [-90, 90] (``direction.inertia``)",
            )
        ]
    )
    lines += [""]

    return "\n".join(lines)


if __name__ == "__main__":
    print(build_rst())
