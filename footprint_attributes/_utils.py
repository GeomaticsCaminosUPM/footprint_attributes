import  geopandas as gpd 
import pandas as pd
import shapely 
import numpy as np
from shapely.geometry import LineString, Point
from packaging.version import Version
import shapely


"""
def get_angle_sharp(normal_0,normal_1,geom_id_0,geom_id_1):
    if geom_id_0 != geom_id_1:
        return 0 
     
    _normal_0 = normal_0 / np.sqrt(np.sum(normal_0 ** 2))
    _normal_1 = normal_1 / np.sqrt(np.sum(normal_1 ** 2))

    # Compute the dot and cross products
    dot = np.dot(_normal_0, _normal_1)  # Dot product

    if np.abs(dot) > 1:
        if np.abs(dot) > 1.0000001:
            warnings.warn(f"invalid value encountered in np.arccos({dot})", RuntimeWarning)

        dot = 1 * dot/np.abs(dot)

    cross = np.cross(_normal_0, _normal_1)  # Cross product (2D equivalent)
    
    if cross < 0:
        return 0 

    if np.abs(cross) > 1:
        if np.abs(cross) > 1.0000001:
            warnings.warn(f"invalid value encountered in np.arcsin({cross})", RuntimeWarning)

        cross = 1 * cross/np.abs(cross)

    angle = np.arccos(np.abs(dot))
    return angle
"""


def get_angle_90(normal_0,normal_1,geom_id_0=0,geom_id_1=0):
    if geom_id_0 != geom_id_1:
        return 0 

    l_normal_0 = np.sqrt(np.sum(normal_0 ** 2))
    l_normal_1 = np.sqrt(np.sum(normal_1 ** 2))
     
    if (l_normal_0 == 0) or (l_normal_1 == 0):
        return 0

    _normal_0 = normal_0 / l_normal_0
    _normal_1 = normal_1 / l_normal_1
    dot = _normal_0[0] * _normal_1[0] + _normal_0[1] * _normal_1[1]
    if np.abs(dot) > 1:
        if np.abs(dot) > 1.0000001:
            warnings.warn(f"invalid value encountered in np.arccos({np.abs(dot)})", RuntimeWarning)

        dot = 1

    angle = np.arccos(np.abs(dot))

    return angle


def get_angle(normal_0, normal_1, geom_id_0=0, geom_id_1=0):
    # If geometry IDs are different, the angle is undefined (return 0)
    if geom_id_0 != geom_id_1:
        return 0
    
    # Normalize the input vectors
    _normal_0 = normal_0 / np.linalg.norm(normal_0)
    _normal_1 = normal_1 / np.linalg.norm(normal_1)
    
    # Compute the dot and cross products
    dot = np.dot(_normal_0, _normal_1)  # Dot product
    if np.abs(dot) > 1:
        if np.abs(dot) > 1.0000001:
            warnings.warn(f"invalid value encountered in np.arccos({dot})", RuntimeWarning)

        dot = 1 * dot/np.abs(dot)

    cross = np.cross(_normal_0, _normal_1)  # Cross product (2D equivalent)
    
    if np.abs(cross) > 1:
        if np.abs(cross) > 1.0000001:
            warnings.warn(f"invalid value encountered in np.arcsin({cross})", RuntimeWarning)

        cross = 1 * cross/np.abs(cross)

    # Calculate the angle using arctan2 for accurate quadrant handling
    angle = np.arctan2(cross, dot)
    
    return angle


def linestring_to_multilinestring(linestring):
    return shapely.MultiLineString(
        [
            LineString(
                (Point(linestring.coords[i]),
                Point(linestring.coords[i+1]))
            ) for i in range(len(linestring.coords)-1)
        ]
    )

def select_touching_edges(gdf,buffer=0):
    gdf_copy = gdf.copy()
    crs = gdf_copy.geometry.crs
    gdf_union = shapely.buffer(gdf_copy.geometry.union_all(),max(buffer,0.0001),cap_style='square',join_style='mitre')
    gdf_union = shapely.buffer(gdf_union,-buffer-0.001,cap_style='square',join_style='mitre')

    gdf_copy.geometry = gdf_copy.geometry.buffer(max(buffer,0.001),cap_style='square',join_style='mitre')
    gdf_copy.geometry = gdf_copy.geometry.buffer(min(-buffer,-0.001),cap_style='square',join_style='mitre')
    gdf_copy.geometry = gdf_copy.geometry.boundary.intersection(gdf_union)

    return gdf_copy


def explode_edges(gdf,min_length=0):
    gdf_copy = gdf.copy()
    crs = gdf.crs
    gdf_copy = gdf_copy.loc[gdf_copy.geometry.is_empty == False]
    gdf_copy = gdf_copy.explode(index_parts=False).reset_index(drop=True)
    gdf_copy['edges'] = gdf_copy.apply(lambda x: linestring_to_multilinestring(x['geometry']),axis=1)
    gdf_copy = gdf_copy.set_geometry('edges',crs=crs).explode().reset_index(drop=True)
    gdf_copy = gdf_copy.loc[gdf_copy['edges'].length > max(min_length,0.001),:]
    return gdf_copy


# Function to get normal vector at the center of a LineString with scaling factor x
def get_normal(linestring, scale=1):
    # 1. Get center point (by length interpolation)
    center_point = linestring.interpolate(0.5, normalized=True)
    
    # 2. Find the two closest points around the center for the tangent
    p1 = Point(linestring.coords[0])
    p2 = Point(linestring.coords[1])
    if len(linestring.coords) > 2:
        raise Exception("Linestring has more than 2 points.")
    # Tangent vector between p1 and p2
    tangent_vector = np.array([p2.x - p1.x, p2.y - p1.y])
    
    # 3. Compute normal vector (perpendicular to tangent)
    normal_vector = np.array([-tangent_vector[1], tangent_vector[0]])
    normal_vector /= np.sqrt(np.sum(normal_vector**2))
    
    # 4. Scale the normal vector
    length_of_linestring = linestring.length
    if scale == 0:
        scaled_normal_vector = normal_vector
    else: 
        scaled_normal_vector = normal_vector * (scale * length_of_linestring)
    
    return center_point, scaled_normal_vector


# Function to calculate momentum relative to a reference point
def calculate_momentum(center_point, normal_vector, reference_point,min_dist=0):
    r = np.array([center_point.x - reference_point.x, center_point.y - reference_point.y])
    cross_product = float(r[0] * normal_vector[1] - r[1] * normal_vector[0])
    if np.sqrt(np.sum(r**2)) < min_dist:
        if cross_product > 0:
            return np.array([cross_product,0,cross_product, 0])
        else:
            return np.array([cross_product,0,0, cross_product])
    else:
        return np.array([cross_product,cross_product,cross_product, cross_product])


def _cast(collection):
    """
    Cast a collection to a shapely geometry array.
    """
    try:
        import  geopandas as gpd
        import shapely
    except (ImportError, ModuleNotFoundError) as exception:
        raise type(exception)(
            "shapely and gpd are required for shape statistics."
        ) from None

    if Version(shapely.__version__) < Version("2"):
        raise ImportError("Shapely 2.0 or newer is required.")

    if isinstance(collection, gpd.GeoSeries | gpd.GeoDataFrame):
        return np.asarray(collection.geometry.array)
    else:
        if isinstance(collection, np.ndarray | list):
            return np.asarray(collection)
        else:
            return np.array([collection])

def _second_moa_ring_xplusy(points):
    """
    implementation of the moment of area for a single ring
    """
    return sum(
        (points[:-1,0] * points[1:,1] - points[1:,0] * points[:-1,1]) * (
            points[1:,0]**2
            + points[1:,0] * points[:-1,0]
            + points[:-1,0]**2
            + points[1:,1]**2
            + points[1:,1] * points[:-1,1]
            + points[:-1,1]**2
        )
    ) / 12 


def _second_moment_of_area_polygon(polygon):
    """
    Compute the absolute value of the moment of area (i.e. ignoring winding direction)
    for an input polygon.
    """
    coordinates = shapely.get_coordinates(polygon)
    centroid = shapely.centroid(polygon)
    centroid_coords = shapely.get_coordinates(centroid)
    moi = _second_moa_ring_xplusy(coordinates - centroid_coords)
    return abs(moi)

def calc_inertia(collection):
    """
    Using equation listed on en.wikipedia.org/wiki/Second_moment_of_area#Any_polygon, the second
    moment of area is the sum of the inertia across the x and y axes:

    The :math:`x` axis is given by:

    .. math::

        I_x = (1/12)\\sum^{N}_{i=1} (x_i y_{i+1} - x_{i+1}y_i) (x_i^2 + x_ix_{i+1} + x_{i+1}^2)

    While the :math:`y` axis is in a similar form:

    .. math::

        I_y = (1/12)\\sum^{N}_{i=1} (x_i y_{i+1} - x_{i+1}y_i) (y_i^2 + y_iy_{i+1} + y_{i+1}^2)

    where :math:`x_i`, :math:`y_i` is the current point and :math:`x_{i+1}`, :math:`y_{i+1}` is the next point,
    and where :math:`x_{n+1} = x_1, y_{n+1} = y_1`. For multipart polygons with holes,
    all parts are treated as separate contributions to the overall centroid, which
    provides the same result as if all parts with holes are separately computed, and then
    merged together using the parallel axis theorem.

    References
    ----------
    Hally, D. 1987. The calculations of the moments of polygons. Canadian National
    Defense Research and Development Technical Memorandum 87/209.
    https://apps.dtic.mil/dtic/tr/fulltext/u2/a183444.pdf

    """  # noqa: E501
    ga = _cast(collection)
    import  geopandas as gpd  # function level, to follow module design

    # construct a dataframe of the fundamental parts of all input polygons
    parts, collection_ix = shapely.get_parts(ga, return_index=True)
    rings, ring_ix = shapely.get_rings(parts, return_index=True)
    # get_rings() always returns the exterior first, then the interiors
    collection_ix = np.repeat(
        collection_ix, shapely.get_num_interior_rings(parts) + 1
    )
    # we need to work in polygon-space for the algorithms
    # (centroid, shoelace calculation) to work
    polygon_rings = shapely.polygons(rings)
    is_external = np.zeros_like(collection_ix).astype(bool)
    # the first element is always external
    is_external[0] = True
    # and each subsequent element is external iff
    # it is different from the preceeding index
    is_external[1:] = ring_ix[1:] != ring_ix[:-1]
    # now, our analysis frame contains a bunch of (guaranteed-to-be-simple) polygons
    # that represent either exterior rings or holes
    polygon_rings = gpd.GeoDataFrame(
        dict(
            collection_ix=collection_ix,
            ring_within_geom_ix=ring_ix,
            is_external_ring=is_external,
        ),
        geometry=polygon_rings,
    )
    # the polygonal moi can be calculated using the same ring-based strategy,
    # and this could be parallelized if necessary over the elemental shapes with:

    # from joblib import Parallel, parallel_backend, delayed
    # with parallel_backend('loky', n_jobs=-1):
    #     engine = Parallel()
    #     promise = delayed(_second_moment_of_area_polygon)
    #     result = engine(promise(geom) for geom in polygon_rings.geometry.values)

    # but we will keep simple for now
    polygon_rings["moa"] = polygon_rings.geometry.apply(_second_moment_of_area_polygon)
    # the above algorithm computes an unsigned moa
    # to be insensitive to winding direction.
    # however, we need to subtract the moa of holes. Hence, the sign of the moa is
    # -1 when the polygon is an internal ring and 1 otherwise:
    polygon_rings["sign"] = (1 - polygon_rings.is_external_ring * 2) * -1
    # shapely already uses the correct formulation for centroids
    polygon_rings["centroids"] = shapely.centroid(polygon_rings.geometry)
    # the inertia of parts applies to the overall center of mass:
    original_centroids = shapely.centroid(ga)
    polygon_rings["collection_centroid"] = original_centroids[collection_ix]
    # proportional to the squared distance between the original and part centroids:
    polygon_rings["radius"] = shapely.distance(
        polygon_rings.centroid.values, polygon_rings.collection_centroid.values
    )
    # now, we take the sum of (I+Ar^2) for each ring, treating the
    # contribution of holes as negative.
    # Then, we take the sum of all of the contributions
    return (
        polygon_rings.groupby(["collection_ix", "ring_within_geom_ix"])
        .apply(
            lambda ring_in_part: (
                (ring_in_part.moa + ring_in_part.radius**2 * ring_in_part.area)
                * ring_in_part.sign
            ).sum()
        ,include_groups=False)
        .groupby(level="collection_ix")
        .sum()
        .values
    )

def resultant_angle(gdf,vector_column='normal_vector',id_column='geom_id'):
    gdf_copy = gdf.copy()
    resultant_vector = gdf_copy[[vector_column,id_column,'geometry']].copy()
    resultant_vector = resultant_vector.groupby(id_column).agg({vector_column:'sum'}).reset_index()
    resultant_vector = resultant_vector.rename(columns={vector_column:'resultant_vector'})
    gdf_copy = gdf_copy.merge(resultant_vector,on=id_column,how='left')

    gdf_copy['angle'] = gdf_copy.apply(lambda x: pd.Series(get_angle_90(x[vector_column],x['resultant_vector'],x[id_column],x[id_column])),axis=1)
    return gdf_copy


def explode_exterior_and_interior_rings(gdf):
    gdf_copy = gdf.copy()
    crs = gdf.crs
    gdf_copy['exterior'] = gdf_copy.exterior
    gdf_copy['interiors'] = gdf_copy.interiors
    gdf_copy.geometry = gdf_copy.apply(lambda x: shapely.MultiLineString(x['interiors'] + [x['exterior']]),axis=1)
    gdf_copy.crs = crs
    gdf_copy = gdf_copy.drop(columns=['exterior','interiors'])
    gdf_copy = gdf_copy.explode().reset_index(drop=True)
    return gdf_copy

def eq_circle_intertia(area):
    r = np.sqrt(area / np.pi)
    return 0.25*np.pi*r**4