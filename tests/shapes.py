"""Re-export shim: the actual factory functions now live in
``footprint_attributes.testing_shapes`` so example notebooks can import them
without a ``sys.path``/``tests`` hack. Kept here so existing test imports
(``from shapes import ...``) keep working unchanged.
"""

from __future__ import annotations

from footprint_attributes.testing_shapes import (  # noqa: F401
    CRS,
    gdf_of,
    asymmetric_l_shape_polygon,
    circle_polygon,
    confined_quartet_gdf,
    corner_triplet_gdf,
    isolated_building_gdf,
    l_shape_polygon,
    lateral_pair_gdf,
    rect_polygon,
    rotated_hole_building_polygon,
    square_polygon,
    square_with_hole_polygon,
    thin_cross_polygon,
    torque_triplet_gdf,
    two_isolated_buildings_gdf,
    t_shape_polygon,
    x_shape_polygon,
)
