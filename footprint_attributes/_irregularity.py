import geopandas as gpd 
import pandas as pd
import shapely 
import numpy as np
import warnings
from ._utils import get_normal, explode_edges, explode_exterior_and_interior_rings


def calc_shape_irregularity(geoms:gpd.GeoDataFrame) -> gpd.GeoDataFrame:
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
        gpd.GeoDataFrame: The input GeoDataFrame with an added column `"shape_irregularity"`,
                          which contains the computed irregularity index for each geometry.
    """

    if "shape_irregularity" in geoms.columns:
        warnings.warn("The 'shape_irregularity' column already exists and will be overwritten.", UserWarning)
        
    if "geom_id" in geoms.columns:
        warnings.warn("The 'geom_id' column already exists and will be overwritten.", UserWarning)
        
    geoms['geom_id'] = geoms.index.copy()
    geoms_copy = geoms.copy()
    if not geoms_copy.crs.is_projected:
        geoms_copy = geoms_copy.to_crs(geoms_copy.geometry.estimate_utm_crs())


    geoms_copy = explode_exterior_and_interior_rings(irregularity)
    boundary = geoms_copy.geometry.copy()
    convex_hull = geoms_copy.geometry.convex_hull.boundary.copy()

    irregularity.geometry = shapely.difference(boundary, convex_hull.buffer(0.005)) #Path through irregularity
    irregularity['hull_geom'] = shapely.difference(convex_hull, boundary.buffer(0.005)) #Path through hull
    irregularity = irregularity.loc[irregularity.geometry.is_empty == False]

    geoms_copy['irregularity_geom'] = shapely.difference(boundary,convex_hull) #Path through irregularity
    geoms_copy = geoms_copy.loc[geoms_copy['irregularity_geom'].is_empty == False,:]
    geoms_copy['hull'] = shapely.difference(convex_hull,boundary) #Path through hull
    geoms_copy['irregularity_length'] = geoms_copy['irregularity_geom'].length - irregularity['hull'].length
    
    geoms_copy = geoms_copy.explode(index_parts=False).reset_index(drop=True)

    geoms_copy = explode_edges(geoms_copy,min_length=0)

    geoms_copy[['edge_center','normal_vector']] = geoms_copy.apply(lambda x: pd.Series(get_normal(x['edges'],1)),axis=1)

    geoms_copy['distance_to_hull'] = geoms_copy.apply(lambda x: shapely.distance(x['edge_center'],x['hull_geom']),axis=1)
    geoms_copy['edge_length'] = geoms_copy['edges'].length

    geoms_copy['shape_irregularity'] = geoms_copy['distance_to_hull'] * geoms_copy['edge_length']

    #geoms_copy.loc[0:len(geoms_copy)-1,'normal_1'] = geoms_copy.loc[1:len(geoms_copy),'normal'].reset_index(drop=True)
    #geoms_copy.loc[0:len(geoms_copy)-1,'geom_id_1'] = geoms_copy.loc[1:len(geoms_copy),'geom_id'].reset_index(drop=True)
    #geoms_copy['angle_irregularity'] = geoms_copy.apply(lambda x: pd.Series(get_angle_sharp(x['normal'],x['normal_1'],x['geom_id'],x['geom_id_1'])),axis=1)

    geoms_copy = geoms_copy.groupby('geom_id').agg({
        'shape_irregularity':'sum',
        'irregularity_length':'first',
        'dist_to_hull':'sum',
        #'angle_irregularity':'sum',
        'hull':'first',
        'irregularity_geom':'first'
    })

    #geoms_copy['angle_irregularity'] = geoms_copy['angle_irregularity'] / (np.pi / 2)
    #geoms_copy['angle_irregularity_length'] = geoms_copy['angle_irregularity'] * geoms_copy['irregularity_length']
    geoms_copy = geoms_copy.reset_index(drop=True)

    result = geoms.merge(geoms_copy[['shape_irregularity','geom_id']],on='geom_id',how='left')
    result.loc[result['shape_irregularity'].isna(),'shape_irregularity'] = 0 
    result['shape_irregularity'] = result['shape_irregularity'].astype(float)
    result.index = result['geom_id'].copy()
    result = result.drop(columns='geom_id')

    return result

def calc_polsby_popper(geoms:gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Calculates the Polsby-Popper index for building footprint polygons.

    The Polsby-Popper index is a measure of shape compactness, defined by the formula:
        Polsby-Popper Index = (4 * π * A) / P²
    where:
        - `A`: The area of the polygon.
        - `P`: The perimeter of the polygon.

    Parameters:
        geoms (gpd.GeoDataFrame): A GeoDataFrame containing building footprint geometries as polygons.

    Returns:
        gpd.GeoDataFrame: The input GeoDataFrame with an added column `"polsby_popper"`,
                          which contains the calculated Polsby-Popper index for each geometry.
    """

    # Warn if Polsby-Popper column already exists
    if "polsby_popper" in geoms.columns:
        warnings.warn("The 'polsby_popper' column already exists and will be overwritten.", UserWarning)
        
    # Preserve original CRS for re-projection if needed
    orig_crs = geoms.crs
    
    # Ensure the geometries are in a projected CRS for accurate area and length calculations
    if not geoms.crs.is_projected:
        geoms = geoms.to_crs(geoms.geometry.estimate_utm_crs())
    
    # Calculate the Polsby-Popper compactness score
    geoms['polsby_popper'] = (4 * np.pi * geoms.geometry.area) / (geoms.geometry.boundary.length ** 2)
    
    # Return to the original CRS
    if orig_crs is not None and orig_crs != geoms.crs:
        geoms = geoms.to_crs(orig_crs)
    
    return geoms
    