import geopandas as gpd 
import pandas as pd
import shapely 
import numpy as np
import warnings
from .utils import get_scaled_normal_vector_at_center, get_angle_90, explode_edges, calculate_momentum

def calc_forces(geoms:gpd.GeoDataFrame,buffer:float=0,height_column:str=None) -> gpd.GeoDataFrame:
    """
    - **Parameters**:
        - `geoms`: GeoDataFrame containing building footprints as polygon geometries.
        - `buffer`: Buffer in meters to consider when determining if two buildings are touching each other.
        - `height_column`: Column in the `geoms` DataFrame specifying the height of the building in meters. If `None`, all buildings are assumed to have height 1.

    - **Output**:
        Returns a GeoDataFrame with the following new columns:
            - **"moment"**: Moment of the resultant force with respect to the centroid of the footprint.
            - **"force"**: Magnitude of the sum of all forces on each footprint (resultant force). \( |\sum(\text{force}_i)| \)
            - **"confinement"**: A measure of the amount of force from the total forces that is confined (has an opposing force).  
                Formula:  
                \[
                \frac{\sum(|\text{force}_i|) - |\sum(\text{force}_i)|}{|\sum(\text{force}_i)|}
                \]
            - **"angle"**: Sum of angles of forces with respect to the resultant force, normalized by force magnitude.  
                Formula:  
                \[
                \frac{\sum(|\text{force}_i| \cdot \text{angle}(\text{force}_i, \sum(\text{force}_j)))}{|\sum(\text{force}_i)|}
                \]
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

    geoms_union = geoms_copy.geometry.union_all().buffer(buffer,cap_style='square',join_style='mitre')
    geoms_union = geoms_union.buffer(-buffer,cap_style='square',join_style='mitre')
    
    geoms_copy.geometry = geoms_copy.geometry.buffer(buffer,cap_style='square',join_style='mitre')
    geoms_copy.geometry = geoms_copy.geometry.buffer(-buffer,cap_style='square',join_style='mitre')
    geoms_copy.geometry = geoms_copy.geometry.boundary.intersection(geoms_union)

    geoms_copy = geoms_copy.loc[geoms_copy.geometry.is_empty == False]

    geoms_copy = explode_edges(geoms_copy,geometry_column='geometry')
    geoms_copy = geoms_copy.loc[geoms_copy['edges'].length > buffer,:]

    geoms_copy[['edge_center','normal']] = geoms_copy.apply(lambda x: pd.Series(get_scaled_normal_vector_at_center(x['edges'],x['height'])),axis=1)

    geoms_copy['momentum'] = geoms_copy.apply(lambda x:calculate_momentum(x['edge_center'],x['normal'],x['centroid']),axis=1)

    geoms_copy['abs_foce'] = geoms_copy.apply(lambda x: np.sqrt(x['normal'][0]**2+x['normal'][1]**2),axis=1)
    geoms_copy['sum_of_squares'] = geoms_copy['abs_foce'] ** 2

    res_normal = geoms_copy.copy()
    res_normal = res_normal.groupby('geom_id').agg({'normal':'sum'}).reset_index()
    res_normal = res_normal.rename(columns={'normal':'res_normal'})
    geoms_copy = geoms_copy.merge(res_normal,on='geom_id',how='left')

    geoms_copy['angle'] = geoms_copy.apply(lambda x: pd.Series(get_angle_90(x['normal'],x['res_normal'],x['geom_id'],x['geom_id'])),axis=1)
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
    - **Parameters**:
        - `footprints`: GeoDataFrame outputted by `calc_forces()` with `force`, `confinement`, and `angle` columns.
        - `min_force`: Significance threshold for the resultant force. Default: `0.166`. (E.g., for a square building with height 1 and side length 1, if a touching structure covers only 1/6 of one side, the resultant force would be 1/6.)
        - `angle_significance`: Angle threshold (in radians). Default: \( \pi / 4 \) (45 degrees).
        - `confinement_significance`: Threshold for confinement. Default: `1` (indicating equal amounts of confined and resultant forces).
        - `min_momentum`: Threshold for momentum. Default: `0.0825` (E.g., for a square building with height 1 and side length 1, if a touching structure covers only 1/6 two sides. In the worst case the momentum would be 0.0825.) 

    - **Output**:
        Returns the same GeoDataFrame with a new column `"relative_position"`. Classifies buildings into the following categories (priority order):
        1. **"confined"**: Structures touching on both the left and right lateral sides.
        2. **"corner"**: Structures touching at a corner (determined by force and angle thresholds).
        3. **"partial"**: Structures touching on either the left or right side.
        4. **"isolated"**: No touching structures.
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
