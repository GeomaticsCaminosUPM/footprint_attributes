import geopandas as gpd 
import pandas as pd
import shapely 
from shapely.geometry import Polygon
import numpy as np
import warnings
from ._utils import get_normal, explode_edges, explode_exterior_and_interior_rings, calc_inertia_z, eq_circle_intertia, calc_inertia_all, calc_inertia_principal, get_angle


def shape_irregularity(geoms:gpd.GeoDataFrame) -> list:
    """TODO: Explore normalize by the boundary and not by the hull and some type of normalization 0-1. Explore polsby + shape"""
    """
    Calculates an index to quantify the irregularity of building footprints.

    The irregularity index is computed using the formula:
        Irregularity = (l * d) / L
    where:
        - `l`: Length of the segments outside the convex hull of the shape.
        - `d`: Distance from the center of gravity of the segments outside the hull to the hull.
        - `L`: Total perimeter length of the convex hull.

    Parameters:
        geoms (gpd.GeoDataFrame): A GeoDataFrame containing building footprint geometries as polygons.

    Returns:
        list: A list with the same order as geoms which contains the computed irregularity index for each geometry.
    """

    #if "shape_irregularity" in geoms.columns:
    #    warnings.warn("The 'shape_irregularity' column already exists and will be overwritten.", UserWarning)
 
    geoms = geoms.drop(columns=['shape_irregularity'],errors='ignore')
    geoms_copy = geoms.copy()
    geoms_copy['geom_id'] = geoms_copy.index.copy()

    if not geoms_copy.crs.is_projected:
        geoms_copy = geoms_copy.to_crs(geoms_copy.geometry.estimate_utm_crs())


    geoms_copy = explode_exterior_and_interior_rings(geoms_copy)
    boundary = geoms_copy.geometry.copy()
    convex_hull = geoms_copy.geometry.convex_hull.boundary.copy()
    geoms_copy['hull_length'] = convex_hull.length
    geoms_copy['boundary_length'] = boundary.length

    geoms_copy.geometry = shapely.difference(boundary, convex_hull.buffer(0.005)) #Path through irregularity
    geoms_copy['hull_geom'] = shapely.difference(convex_hull, boundary.buffer(0.005)) #Path through hull
    geoms_copy = geoms_copy.loc[geoms_copy.geometry.is_empty == False]

    geoms_copy = geoms_copy.explode(index_parts=False).reset_index(drop=True)

    geoms_copy = explode_edges(geoms_copy,min_length=0)

    geoms_copy[['edge_center','normal_vector']] = geoms_copy.apply(lambda x: pd.Series(get_normal(x['edges'],1)),axis=1)

    geoms_copy['distance_to_hull'] = geoms_copy.apply(lambda x: shapely.distance(x['edge_center'],x['hull_geom']),axis=1)
    geoms_copy['edge_length'] = geoms_copy['edges'].length

    geoms_copy['shape_irregularity'] = (geoms_copy['distance_to_hull'] * geoms_copy['edge_length']) / geoms_copy['hull_length']

    #geoms_copy.loc[0:len(geoms_copy)-1,'normal_1'] = geoms_copy.loc[1:len(geoms_copy),'normal'].reset_index(drop=True)
    #geoms_copy.loc[0:len(geoms_copy)-1,'geom_id_1'] = geoms_copy.loc[1:len(geoms_copy),'geom_id'].reset_index(drop=True)
    #geoms_copy['angle_irregularity'] = geoms_copy.apply(lambda x: pd.Series(get_angle_sharp(x['normal'],x['normal_1'],x['geom_id'],x['geom_id_1'])),axis=1)

    geoms_copy = geoms_copy.groupby('geom_id').agg({'shape_irregularity':'sum'})

    result = geoms.merge(geoms_copy[['shape_irregularity']],right_index=True,left_index=True,how='left')
    result.loc[result['shape_irregularity'].isna(),'shape_irregularity'] = 0 
    result['shape_irregularity'] = result['shape_irregularity'].astype(float)

    return list(result['shape_irregularity'])

def polsby_popper(geoms:gpd.GeoDataFrame, fill_holes:bool=True) -> list:
    """TODO: Polsby popper donut shape. boundary.length takes both inner and outer circle into accout so that the perimeter is very large"""
    """
    Calculates the Polsby-Popper index for building footprint polygons.

    The Polsby-Popper index is a measure of shape compactness, defined by the formula:
        Polsby-Popper Index = (4 * π * A) / P²
    where:
        - `A`: The area of the polygon.
        - `P`: The perimeter of the polygon.

    Parameters:
        geoms (gpd.GeoDataFrame): A GeoDataFrame containing building footprint geometries as polygons.
        fill_holes (bool, optional): Fill all holes of the geometries (interior patios, etc.) to compute the index. Defaults to True.

    Returns:
        list: A list in the same order as geoms which contains the calculated Polsby-Popper index for each geometry.
    """
    geoms = geoms.copy()
    # Ensure the geometries are in a projected CRS for accurate area and length calculations
    if not geoms.crs.is_projected:
        geoms = geoms.to_crs(geoms.geometry.estimate_utm_crs())
    
    # Calculate the Polsby-Popper compactness score
    if fill_holes:
        geoms_holes_filled = geoms.geometry.apply(
            lambda x: Polygon(x.exterior)
        )
        geoms['polsby_popper'] = (4 * np.pi * geoms_holes_filled.geometry.area) / (geoms_holes_filled.geometry.boundary.length ** 2)
    else:
        geoms['polsby_popper'] = (4 * np.pi * geoms.geometry.area) / (geoms.geometry.boundary.length ** 2)
    
    return list(geoms['polsby_popper']) 
    
def inertia_slenderness(geoms:gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Calculates the inertia irregularity index for building footprint polygons comparing the principal components of the inertia tensor (max and min).

    Parameters:
        geoms (gpd.GeoDataFrame): A GeoDataFrame containing building footprint geometries as polygons.

    Returns:
        list: A list with the same order as geoms which contains the calculated Polsby-Popper index for each geometry.
    """
    geoms = geoms.copy() 
    # Ensure the geometries are in a projected CRS for accurate area and length calculations
    if not geoms.crs.is_projected:
        geoms = geoms.to_crs(geoms.geometry.estimate_utm_crs())

    I_max, I_min = calc_inertia_principal(geoms.geometry,principal_dirs=False)
    return np.sqrt(I_min / I_max)
    
def inertia_circle(geoms:gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Calculates the inertia irregularity index for building footprint polygons comparing inertia with the inertia of a circle with the same area.

    The inertia irregularity index is a measure of shape compactness, defined by the formula:
        inertia irregularity index = eq circle inertia / polygon inertia

    Parameters:
        geoms (gpd.GeoDataFrame): A GeoDataFrame containing building footprint geometries as polygons.

    Returns:
        list: A list with the same order as geoms which contains the calculated Polsby-Popper index for each geometry.
    """
    geoms = geoms.copy()
    
    # Ensure the geometries are in a projected CRS for accurate area and length calculations
    if not geoms.crs.is_projected:
        geoms = geoms.to_crs(geoms.geometry.estimate_utm_crs())
    
    # Calculate the inertia irregularity score comparing to a circle with the same area
    geoms['inertia_circle'] = eq_circle_intertia(geoms.geometry.area) / np.abs(calc_inertia_z(geoms.geometry)) 
    
    return list(geoms['inertia_circle'])

def compactness(geoms:gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    convex_hull = geoms.geometry.convex_hull
    geoms_holes_filled = geoms.geometry.apply(
        lambda x: Polygon(x.exterior)
    )
    return (convex_hull.area - geoms_holes_filled.area) / convex_hull.area

def excentricity_ratio(geoms):
    """
    Computes the excentricity ratio for a set of geometries following EC8 and choosing the worst axis.

    Parameters:
        geoms (GeoDataFrame): A GeoDataFrame containing geometry objects.

    Returns:
        list: A list with the same order as geoms which contains the the computed excentricity ratios.
    """
    import scipy
    
    # Compute principal moments of inertia and their corresponding eigenvectors
    inertia_df = calc_inertia_principal(geoms, principal_dirs=True)

    # Compute eccentricity vectors (difference between centroid and boundary centroid)
    e_vect = geoms.geometry.apply(lambda geom: np.array([
        geom.centroid.x - geom.boundary.centroid.x,
        geom.centroid.y - geom.boundary.centroid.y
    ]))

    # Compute magnitude of eccentricity vectors
    e_magnitude = np.sqrt(np.sum(np.array([*(e_vect * e_vect)]), axis=1))

    # Create DataFrame with necessary parameters
    df = pd.DataFrame({
        'e_vect': e_vect,
        'e_magnitude': e_magnitude,
        'I_1': inertia_df[0],  # First principal moment of inertia
        'I_2': inertia_df[2],  # Second principal moment of inertia
        'vect_1': inertia_df[1],  # First principal axis
        'vect_2': inertia_df[3],  # Second principal axis
        'area': geoms.geometry.to_crs(geoms.geometry.estimate_utm_crs()).area  # Projected area
    })

    # Compute angle `b` based on eccentricity direction and principal axis
    df['b'] = df.apply(
        lambda row: 0 if row['e_magnitude'] <= 10**-10 else get_angle(row['vect_1'], row['e_vect']),
        axis=1
    )

    # Optimize `x` to minimize the function, setting `x_opt` to 0 when `e_magnitude` is small
    df['x_opt'] = df.apply(
        lambda row: 0 if row['e_magnitude'] <= 10**-10 else scipy.optimize.fmin(
            lambda x: - np.cos(x - row['b']) ** 2 * (
                0.5 * (row['I_1'] - row['I_2']) * np.cos(2 * x) + 
                0.5 * (row['I_1'] + row['I_2'])
            ),
            x0=0,
            xtol=1e-5,
            ftol=1e-5,
            disp=False
        )[0],
        axis=1
    )

    # Compute excentricity ratio
    result = np.sqrt(
        df['e_magnitude'] ** 2 * np.cos(df['x_opt'] - df['b']) ** 2 * (
            0.5 * (df['I_1'] - df['I_2']) * np.cos(2 * df['x_opt']) + 
            0.5 * (df['I_1'] + df['I_2'])
        ) / (df['I_1'] + df['I_2'] + df['area'] * df['e_magnitude'] ** 2)
    )

    return list(result)
    
    
