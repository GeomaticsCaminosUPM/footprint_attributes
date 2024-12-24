# Footprint Attributes

**Footprint Attributes** is a Python package built using GeoPandas to calculate various geometric indices related to building footprint geometry and seismic risk.

## Installation 

`pip install git+https://github.com/GeomaticsCaminosUPM/footprint_attributes.git`

## Features

The package provides functions for:

### 1. **Shape Irregularity**
- **Polsby-Popper Index**: A measure of shape compactness (similarity to a circle).
  - #### Function: `calc_polsby_popper(geoms: gpd.GeoDataFrame) -> gpd.GeoDataFrame`
   - Formula: $\text{Polsby-Popper Index} = \frac{4 \pi A}{P^2}$
    where:
    - $A$: Area of the polygon.
    - $P$: Perimeter of the polygon.
    - **Parameters**:
      - `geoms`: GeoDataFrame with building footprint polygon geometries.
    - **Output**: Returns the same GeoDataFrame with a new column `"polsby_popper"`.

- **Custom Irregularity Index**: Our own index to quantify the irregularity of building footprints.
  - Formula: $\frac{l \cdot d}{L}$, where:
    - $l$: Length of the shapes outside the convex hull.
    - $d$: Distance of the center of gravity of the shapes outside the hull to the hull.
    - $L$: Total length of the convex hull.
  - #### Function: `calc_shape_irregularity(geoms: gpd.GeoDataFrame) -> gpd.GeoDataFrame`
    - **Parameters**:
      - `geoms`: GeoDataFrame with building footprint polygon geometries.
    - **Output**: Returns the same GeoDataFrame with a new column `"shape_irregularity"`.
   
- **Inertia Irregularity**: Inertia of a circle with the same area as the polygon geometry divided by the inertia of the polygon.
    - #### Function: `calc_inertia_irregularity(geoms:gpd.GeoDataFrame) -> gpd.GeoDataFrame`
        - Formula: $\text{intertia irregularity} = \frac{\text{intertia eq circle}}{\text{itertia}}$
        - **Parameters**:
          - `geoms`: GeoDataFrame with building footprint polygon geometries.
        - **Output**: Returns the same GeoDataFrame with a new column `"inertia_irregularity"`.

---

### 2. **Relative Position of Buildings**
Determines if the building touches other structures (relative position in the city block). This is done by calculating "forces" that neighboring structures exert on the building. The force is proportional to the contact area (length of touching footprints multiplied by building height) in the normal direction of the touching plane.

#### Function: `calc_forces(geoms: gpd.GeoDataFrame, buffer: float = 0, height_column: str = None) -> gpd.GeoDataFrame`
##### Parameters

- **`geoms` (`gpd.GeoDataFrame`)**  
  A GeoDataFrame containing building footprints as polygon geometries.

- **`buffer` (`float`)**  
  Buffer distance in meters to determine if two buildings are considered touching.

- **`height_column` (`str`, optional)**  
  Column name in `geoms` specifying the building height in meters.  
  If `None`, all buildings are assumed to have a uniform height of 1.

- **`min_radius` (`float`, optional)**  
  Minimum distance multiplier used to exclude forces that would otherwise increase momentum. Forces with a distance below a threshold  
  (`min_radius * equivalent radius`) will contribute to the momentum calculation only if they decrease the momentum. The equivalent radius of a building is defined as the radius of a circle with the same area as the building's footprint.

##### Output

The function returns the input `gpd.GeoDataFrame`, which includes the following new columns:

1. **`angular_acc`**  
   Angular acceleration calculated as: $\text{angular acc} = \frac{\text{momentum} \cdot \text{area}}{\text{inertia}}$
   - **Momentum** is calculated as: $\text{momentum} = \sum (\text{distance} \cdot |\text{force}_i|)$

2. **`force`**  
   Magnitude of the resultant force acting on the footprint, normalized by the square root of the area: $\text{force} = \left| \sum \text{force}_i \right|$

3. **`confinement_ratio`**  
   Proportion of total forces that are confined (counterbalanced by opposing forces): $\text{confinement ratio} = \frac{\sum |\text{force}_i| - \left| \sum \text{force}_i \right|}{\left| \sum \text{force}_i \right|}$

4. **`angle`**  
   Normalized sum of the angles between individual forces and the resultant force: $\text{angle} = \frac{\sum \left( |\text{force}_i| \cdot \text{angle}(\text{force}_i, \sum \text{force}_j) \right)}{\left| \sum \text{force}_i \right|}$

#### Function: `relative_position(footprints: gpd.GeoDataFrame, min_angular_acc: float = 0.0825, min_confinement: float = 1, min_angle: float = 0.78, min_force: float = 0.166) -> gpd.GeoDataFrame`

##### Parameters
  - **`footprints` (`gpd.GeoDataFrame`)**: GeoDataFrame outputted by `calc_forces()` with `force`, `confinement`, and `angle` columns.
  - **`min_force` (`float`, optional)**: Significance threshold for the resultant force. Default: `0.166`. (E.g., for a square building with height 1 and side length 1, if a touching structure covers only 1/6 of one side, the resultant force would be 1/6.)
  - **`min_angle` (`float`, optional)**: Angle threshold (in radians). Default: $\pi / 4$ (45 degrees).
  - **`min_confinement` (`float`, optional)**: Threshold for confinement. Default: `1` (indicating equal amounts of confined and resultant forces).
  - **`min_angular_acc` (`float`, optional)**: Threshold for angular acceleration $\frac{momentum * area}{inertia}$. Default: 2.133 (e.g., for a rectangular building with height 1 and sides of length 1 and 0.5, 
                            a touching structure covering 1/3 of two sides in the worst case would have an anuglar acceleration of 2.133)

##### Output
  Returns the input `gpd.GeoDataFrame` with a new column `"relative_position"`. Classifies buildings into the following categories (priority order):
  1. **"torque"**: Buildings of class **confined** or **corner** with an angular acceleration exceeding the minimum.
  2. **"confined"**: Structures touching on both the left and right lateral sides.
  3. **"corner"**: Structures touching at a corner (determined by force and angle thresholds).
  4. **"partial"**: Structures touching on either the left or right side.
  5. **"isolated"**: No touching structures.




