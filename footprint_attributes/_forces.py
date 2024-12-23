import geopandas as gpd 
import pandas as pd
import shapely 
import numpy as np
import warnings
from ._utils import get_normal, explode_edges, calculate_momentum, select_touching_edges, resultant_angle

def calc_forces(geoms:gpd.GeoDataFrame,buffer:float=0,height_column:str=None) -> gpd.GeoDataFrame:
    """
    Calculates force-based metrics for building footprints based on their geometry and proximity.

    Parameters:
        geoms (gpd.GeoDataFrame): A GeoDataFrame containing building footprints as polygon geometries.
        buffer (float): Buffer distance in meters to determine if two buildings are considered touching.
        height_column (str, optional): Column name in `geoms` specifying the building height in meters.
                                       If `None`, all buildings are assumed to have a uniform height of 1.

    Returns:
        gpd.GeoDataFrame: The input GeoDataFrame with the following new columns:
            - "momentum": Momentum of the resultant force with respect to the footprint centroid.
              Formula: sum(distance * |force_i|)
            - "force": Magnitude of the resultant force acting on the footprint.
              Formula: |sum(force_i)|
            - "confinement": Proportion of total forces that are confined (counterbalanced by opposing forces).
              Formula: (sum(|force_i|) - |sum(force_i)|) / |sum(force_i)|
            - "angle": Normalized sum of the angles between individual forces and the resultant force.
              Formula: sum(|force_i| * angle(force_i, sum(force_j))) / |sum(force_i)|

    Notes:
        - This method assumes that forces are calculated based on proximity and building height.
        - The buffer determines which buildings are close enough to interact.
    """

    if "force" in geoms.columns:
        warnings.warn("The 'force' column already exists and will be overwritten.", UserWarning)
        
    if "confinement" in geoms.columns:
        warnings.warn("The 'confinement' column already exists and will be overwritten.", UserWarning)
        
    if "momentum" in geoms.columns:
        warnings.warn("The 'momentum' column already exists and will be overwritten.", UserWarning)
        
    if "angle" in geoms.columns:
        warnings.warn("The 'angle' column already exists and will be overwritten.", UserWarning)
        
    if "geom_id" in geoms.columns:
        warnings.warn("The 'geom_id' column already exists and will be overwritten.", UserWarning)
        
    orig_crs = geoms.crs
    geoms['geom_id'] = geoms.index.copy()
    geoms_copy = geoms.copy() 
    geoms_copy.geometry = geoms_copy.geometry.force_2d()

    if type(height_column) == type(None):
        geoms_copy['height'] = 1
    else:
        geoms_copy['height'] = geoms_copy[height_column].astype(float)

    if not geoms_copy.crs.is_projected:
        geoms_copy = geoms_copy.to_crs(geoms_copy.geometry.estimate_utm_crs())

    geoms_copy = select_touching_edges(geoms_copy,buffer=buffer)
    geoms_copy = explode_edges(geoms_copy,min_length=buffer)

    geoms_copy[['edge_center','normal_vector']] = geoms_copy.apply(lambda x: pd.Series(get_normal(x['edges'],scale=x['height'])),axis=1)

    geoms_copy['momentum'] = geoms_copy.apply(lambda x:calculate_momentum(x['edge_center'],x['normal_vector'],x['centroid']),axis=1)

    geoms_copy['abs_foce'] = geoms_copy.apply(lambda x: np.sqrt(x['normal_vector'][0]**2+x['normal_vector'][1]**2),axis=1)
    geoms_copy['sum_of_squares'] = geoms_copy['abs_foce'] ** 2

    geoms_copy = resultant_angle(geoms_copy,vector_column='normal_vector',id_column='geom_id')

    geoms_copy['angle'] = geoms_copy['angle'] * 2
    geoms_copy['angle'] = geoms_copy['angle'] * geoms_copy['abs_foce']

    geoms_copy = geoms_copy.groupby('geom_id').agg({
        'height':'first',
        'normal':'sum',
        'abs_foce':'sum',
        'angle':'sum',
        'angle':'sum',
        'sum_of_squares':'sum',
        'momentum':'sum',
        'polygon':'first'
    })

    geoms_copy = geoms_copy.rename(columns={'normal':'res_normal'})
    geoms_copy['sum_of_squares'] = np.sqrt(geoms_copy['sum_of_squares'])
    geoms_copy = geoms_copy.set_geometry('polygon',crs=crs)

    geoms_copy['force'] = geoms_copy.apply(lambda x: np.sqrt(x['res_normal'][0]**2+x['res_normal'][1]**2),axis=1)
    geoms_copy['confinement'] = (geoms_copy['abs_force'] - geoms_copy['force']) / geoms_copy['force'] #(geoms_copy['sum_of_squares'] - geoms_copy['force'])
    geoms_copy['momentum'] = np.abs(geoms_copy['momentum'])
    geoms_copy['angle'] = geoms_copy['angle'] / geoms_copy['abs_foce']

    #geoms_copy['area'] = geoms_copy['polygon'].area 
    #geoms_copy['area_sqrt'] = np.sqrt(geoms_copy['polygon'].area)
    #geoms_copy['envelope'] = geoms_copy['polygon'].envelope.length
    geoms_copy = geoms_copy.reset_index()

    result = geoms.merge(geoms_copy[['force','confinement','momentum','angle']],on='geom_id',how='left')
    result.loc[result['force'].isna(),'force'] = 0 
    result.loc[result['confinement'].isna(),'confinement'] = 0 
    result.loc[result['momentum'].isna(),'momentum'] = 0 
    result.loc[result['angle'].isna(),'angle'] = 0 
    result[['force','confinement','momentum','angle']] = result[['force','confinement','momentum','angle']].astype(float)
    result.index = result['geom_id'].copy()
    result = result.drop(columns='geom_id')

    return result

def relative_position(footprints: gpd.GeoDataFrame, min_momentum: float = 0.0825, min_confinement: float = 1, min_angle: float = 0.78, min_force: float = 0.166) -> gpd.GeoDataFrame:
    """
    Classifies building footprints based on their relative position and interaction with surrounding structures.

    Parameters:
        footprints (gpd.GeoDataFrame): GeoDataFrame outputted by `calc_forces()` containing `force`, `confinement`, and `angle` columns.
        min_force (float, optional): Threshold for the significance of the resultant force.
                                     Default: 0.166 (e.g., for a square building with height 1 and side length 1, 
                                     a touching structure covering 1/6 of one side would have a resultant force of 1/6).
        min_angle (float, optional): Angle threshold in radians. Default: pi / 4 (45 degrees).
        min_confinement (float, optional): Threshold for confinement. Default: 1 (equal amounts of confined and resultant forces).
        min_momentum (float, optional): Threshold for momentum. Default: 0.0825 
                                        (e.g., for a square building with height 1 and side length 1, 
                                        a touching structure covering 1/6 of two sides in the worst case 
                                        would have a momentum of 0.0825).

    Returns:
        gpd.GeoDataFrame: The input GeoDataFrame with a new column `relative_position`, classifying buildings into the following categories (priority order):
            1. "torque": Buildings classified as *confined* or *corner* with momentum exceeding the minimum momentum.
            2. "confined": Structures touching on both the left and right lateral sides.
            3. "corner": Structures touching at a corner, based on force and angle thresholds.
            4. "partial": Structures touching on either the left or right side.
            5. "isolated": Structures with no touching neighbors.

    Notes:
        - The classification prioritizes categories in the order listed above.
    """

    #footprints['relative_position'] = 'isolated'
    #footprints.loc[(footprints['force'] / footprints['area_sqrt']) > 0.05,'relative_position'] = 'lateral'
    #footprints.loc[(footprints['angle_normalized'] > 0.6) & ((footprints['force'] / footprints['area_sqrt']) > 0.35),'relative_position'] = 'corner'
    #footprints.loc[(footprints['confinement'] / footprints['area_sqrt']) > 0.07,'relative_position'] = 'confined'
    # Warn if relative_position column already exists

    if "relative_position" in footprints.columns:
        warnings.warn(
            "The 'relative_position' column already exists and will be overwritten.", UserWarning
        )
    
    # Preserve original CRS for re-projection if needed
    orig_crs = footprints.crs
    
    # Ensure the geometries are in a projected CRS for accurate calculations
    if not footprints.crs.is_projected:
        footprints = footprints.to_crs(footprints.geometry.estimate_utm_crs())
    
    # Initialize 'relative_position' column
    footprints['relative_position'] = 'isolated'
    
    # Precompute common values for efficiency
    normalized_force = footprints['force'] / np.sqrt(footprints.geometry.area)
    normalized_momentum = footprints['momentum'] / np.sqrt(footprints.geometry.area)
    confinement_ratio = footprints['confinement']# / np.sqrt(footprints.geometry.area)
    
    # Update 'relative_position' based on criteria
    footprints.loc[normalized_force > min_force, 'relative_position'] = 'partial'
    
    footprints.loc[
        (
            footprints['angle'] > angle_significance
        ) & (
            normalized_force > min_force
        ),
        'relative_position'
    ] = 'corner'
    
    footprints.loc[confinement_ratio > confinement_significance, 'relative_position'] = 'confined'

    footprints.loc[
        ((footprints['relative_position'] == 'corner') | (footprints['relative_position'] == 'confined')
        ) & (
            normalized_momentum > min_momentum
        ), 
        'relative_position'
    ] = 'torque'
    
    # Return to the original CRS if needed
    if orig_crs is not None and orig_crs != footprints.crs:
        footprints = footprints.to_crs(orig_crs)
    
    return footprints
