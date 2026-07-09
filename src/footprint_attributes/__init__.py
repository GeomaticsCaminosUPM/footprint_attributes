"""
footprint_attributes
====================

Automated computation of seismic behaviour modifiers from 2-D building
footprint polygons.

Main modules
-----------
- **direction**: Building orientation (principal axes, bounding box, eccentricity axes)
- **shape**:     Plan-irregularity indices (EC8, ASCE7, GNDTII, CSCR2010, NTC23,
                 slenderness, polsby_popper, convex_hull_irregularity,
                 inertia_circle_ratio)
- **position**:  Relative position within urban block (isolated / lateral /
                 corner / confined / torque)
- **run**:       Single entry point -- request columns by name (or by norm)
                 across direction/shape/position in one call.

Intended import style
---------------------
All functionality is accessed through the three top-level module objects::

    from footprint_attributes import direction, shape, position

Or, for a one-shot "give me these columns" call::

    import footprint_attributes
    result = footprint_attributes.run(footprints, config={"columns": ["EC8", "position"]})

Do **not** import norm aggregators (``EC8``, ``ASCE7`` …) directly from the
package; use ``shape.EC8``, ``shape.ASCE7``, etc. instead.

Example usage
-------------
>>> import geopandas as gpd
>>> from footprint_attributes import direction, shape, position
>>>
>>> footprints = gpd.read_file("footprints.gpkg")
>>>
>>> # Building direction — several methods
>>> bearings = direction.inertia(footprints)
>>> L1, dir1, L2, dir2, bearings = direction.bbox(footprints, mode="all")
>>>
>>> # Worst-case eccentricity analysis direction (EC8)
>>> x_opt = direction.eccentricity.EC8(footprints)
>>>
>>> # All EC8 shape parameters
>>> result = shape.EC8(footprints)
>>>
>>> # Single parameter (returns list)
>>> ecc_ratio = shape.EC8.eccentricityRatio(footprints)
>>>
>>> # Code-independent compactness index
>>> pp = shape.polsby_popper(footprints)
>>>
>>> # Plan slenderness
>>> result = shape.slenderness(footprints)          # all methods → gdf
>>> vals   = shape.slenderness.inertia(footprints)  # single method → list
>>>
>>> # Batch call — any mix of columns, shared work computed once
>>> result = shape(footprints, ["EC8_eccentricityRatio", "slenderness_bbox"])
>>>
>>> # Position classification with contact forces
>>> result = position(footprints)
"""

__version__ = "1.0.0"

# ── Public module objects ───────────────────────────────────────────────────
from . import direction  # noqa: F401  direction.inertia / direction.bbox / direction.eccentricity
from . import shape as _shape_module  # noqa: F401
from . import position as _position_module  # noqa: F401

# Re-export the *callable* shape and position objects under their natural names
from .shape import shape  # noqa: F401  shape(gdf) / shape.EC8 / shape.slenderness …
from .position import position  # noqa: F401  position(gdf)
from .runner import run  # noqa: F401  run(gdf, config={"columns": [...]})

# ── __all__ ─────────────────────────────────────────────────────────────────
# Only expose the three namespace objects (plus run); internal modules
# (geometry, config, eccentricity, runner) are not part of the public API.
__all__ = [
    "direction",
    "shape",
    "position",
    "run",
]
