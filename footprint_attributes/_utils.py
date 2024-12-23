import geopandas as gpd 
import pandas as pd
import shapely 
import numpy as np
from shapely.geometry import LineString, Point


"""
def get_angle_sharp(normal_0,normal_1,geom_id_0,geom_id_1):
    if geom_id_0 != geom_id_1:
        return 0 
     
    _normal_0 = normal_0 / np.sqrt(np.sum(normal_0 ** 2))
    _normal_1 = normal_1 / np.sqrt(np.sum(normal_1 ** 2))
    dot = _normal_0[0] * _normal_1[0] + _normal_0[1] * _normal_1[1]
    cross = _normal_0[0] * _normal_1[1] - _normal_0[1] * _normal_1[0]
    if cross < 0:
        return 0 

    angle = np.arccos(np.abs(dot))
    return angle
"""


def get_angle_90(normal_0,normal_1,geom_id_0=0,geom_id_1=0):
    if geom_id_0 != geom_id_1:
        return 0 
     
    _normal_0 = normal_0 / np.sqrt(np.sum(normal_0 ** 2))
    _normal_1 = normal_1 / np.sqrt(np.sum(normal_1 ** 2))
    dot = _normal_0[0] * _normal_1[0] + _normal_0[1] * _normal_1[1]
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
    cross = np.cross(_normal_0, _normal_1)  # Cross product (2D equivalent)
    
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
    gdf_union = shapely.buffer(gdf_copy.geometry.union_all(),max(buffer,0.001),cap_style='square',join_style='mitre')
    gdf_union = shapely.buffer(gdf_union,min(-buffer,-0.001),cap_style='square',join_style='mitre')

    gdf_copy.geometry = gdf_copy.geometry.buffer(buffer,cap_style='square',join_style='mitre')
    gdf_copy.geometry = gdf_copy.geometry.buffer(-buffer,cap_style='square',join_style='mitre')
    gdf_copy.geometry = gdf_copy.geometry.boundary.intersection(gdf_union)

    return gdf_copy


def explode_edges(gdf,min_length=0):
    gdf_copy = gdf.copy()
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
def calculate_momentum(center_point, normal_vector, reference_point):
    r = np.array([center_point.x - reference_point.x, center_point.y - reference_point.y])
    cross_product = r[0] * normal_vector[1] - r[1] * normal_vector[0]
    return cross_product


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

