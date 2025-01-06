```markdown
# Footprint Attributes

**Footprint Attributes** is a Python package built using GeoPandas to calculate various geometric indices related to building footprint geometry and seismic risk.

---

## Installation

To install the package, use the following command:

```bash
pip install git+https://github.com/GeomaticsCaminosUPM/footprint_attributes.git
```

---

## Features

### 1. **Relative Position of Buildings**
This feature determines if a building touches other structures (relative position within the city block). It calculates "forces" that neighboring structures exert on the building, proportional to the contact area (length of touching footprints multiplied by building height) in the normal direction of the touching plane.

#### **Function: `get_forces_gdf`**
```python
get_forces_gdf(geoms: gpd.GeoDataFrame, buffer: float = 0, height_column: str = None) -> gpd.GeoDataFrame
```

##### Parameters:
- **`geoms`** (`gpd.GeoDataFrame`): A GeoDataFrame containing building footprints as polygon geometries.
- **`buffer`** (`float`): Buffer distance in meters to determine if two buildings are considered touching.
- **`height_column`** (`str`, optional): Column name specifying building height. Defaults to `1` if `None`.

##### Output:
Returns a `gpd.GeoDataFrame` with additional columns:
- **`height`**: Building height.
- **`angular_acc`**: Angular acceleration.
- **`force`**: Magnitude of the resultant force.
- **`confinement_ratio`**: Proportion of confined forces.
- **`angle`**: Normalized sum of force angles.
- **`geometry`**: Original building footprint geometries.

---

#### **Function: `relative_position`**
```python
relative_position(
    forces: gpd.GeoDataFrame,
    min_angular_acc: float = 0.0825,
    min_confinement: float = 1,
    min_angle: float = 0.78,
    min_force: float = 0.166
) -> list
```

##### Parameters:
- **`forces`** (`gpd.GeoDataFrame`): Output GeoDataFrame from `get_forces_gdf`.
- **`min_force`** (`float`): Minimum force threshold (default: `0.166`).
- **`min_angle`** (`float`): Minimum angle threshold (default: `π/4` radians or 45 degrees).
- **`min_confinement`** (`float`): Minimum confinement threshold (default: `1`).
- **`min_angular_acc`** (`float`): Minimum angular acceleration threshold.

##### Output:
Returns a list of relative positions for buildings, classified as:
1. **"torque"**: High angular acceleration.
2. **"confined"**: Touches on both lateral sides.
3. **"corner"**: Touches at a corner.
4. **"partial"**: Touches on one side.
5. **"isolated"**: No touching structures.

---

### 2. **Shape Irregularity**
Measures geometric irregularity of building footprints using various indices.

#### **Polsby-Popper Index**
Measures shape compactness (similarity to a circle).
- **Formula**:  
  \[
  \text{Polsby-Popper Index} = \frac{4 \pi A}{P^2}
  \]
  where:
  - \( A \): Area of the polygon.
  - \( P \): Perimeter of the polygon.

##### Function: `calc_polsby_popper`
```python
calc_polsby_popper(geoms: gpd.GeoDataFrame) -> list
```

- **Parameters**:  
  - `geoms`: GeoDataFrame with building footprint geometries.
- **Output**: List of `polsby_popper` indices corresponding to `geoms` rows.

---

#### **Custom Irregularity Index**
Quantifies the irregularity of footprints using convex hull analysis.
- **Formula**:  
  \[
  \text{Custom Irregularity Index} = \frac{l \cdot d}{L}
  \]
  where:
  - \( l \): Length outside the convex hull.
  - \( d \): Distance of the center of gravity outside the hull.
  - \( L \): Total convex hull length.

##### Function: `calc_shape_irregularity`
```python
calc_shape_irregularity(geoms: gpd.GeoDataFrame) -> list
```

- **Parameters**:  
  - `geoms`: GeoDataFrame with building footprint geometries.
- **Output**: List of `shape_irregularity` indices.

---

#### **Inertia Irregularity**
Compares the inertia of a polygon to a circle with the same area.
- **Formula**:  
  \[
  \text{Inertia Irregularity} = \frac{\text{Inertia of Equivalent Circle}}{\text{Inertia of Polygon}}
  \]

##### Function: `calc_inertia_irregularity`
```python
calc_inertia_irregularity(geoms: gpd.GeoDataFrame) -> list
```

- **Parameters**:  
  - `geoms`: GeoDataFrame with building footprint geometries.
- **Output**: List of `inertia_irregularity` indices.

---



# Footprint Attributes

**Footprint Attributes** is a Python package built using GeoPandas to calculate various geometric indices related to building footprint geometry and seismic risk.

## Installation 

`pip install git+https://github.com/GeomaticsCaminosUPM/footprint_attributes.git`

## Features

The package provides functions for:


### 1. **Relative Position of Buildings**
Determines if the building touches other structures (relative position in the city block). This is done by calculating "forces" that neighboring structures exert on the building. The force is proportional to the contact area (length of touching footprints multiplied by building height) in the normal direction of the touching plane.

#### Function: `get_forces_gdf(geoms: gpd.GeoDataFrame, buffer: float = 0, height_column: str = None) -> gpd.GeoDataFrame`
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

The function returns a `gpd.GeoDataFrame`, which includes the following columns:

- **`height`**
    Building height. If `height_column` is set to `None` the height will be `1`.
  
- **`angular_acc`**  
   Angular acceleration calculated as: $\text{angular acc} = \frac{\text{momentum} \cdot \text{area}}{\text{inertia}}$
   - **Momentum** is calculated as: $\text{momentum} = \sum (\text{distance} \cdot |\text{force}_i|)$

- **`force`**  
   Magnitude of the resultant force acting on the footprint, normalized by the square root of the area: $\text{force} = \left| \sum \text{force}_i \right|$

- **`confinement_ratio`**  
   Proportion of total forces that are confined (counterbalanced by opposing forces): $\text{confinement ratio} = \frac{\sum |\text{force}_i| - \left| \sum \text{force}_i \right|}{\left| \sum \text{force}_i \right|}$

- **`angle`**  
   Normalized sum of the angles between individual forces and the resultant force: $\text{angle} = \frac{\sum \left( |\text{force}_i| \cdot \text{angle}(\text{force}_i, \sum \text{force}_j) \right)}{\left| \sum \text{force}_i \right|}$

- **`geometry`**
  Geopandas geometry column with the building footprint shapely polygons from the input.

The row indices are the same as the input `geoms` GeoDataframe.

#### Function: `relative_position(forces: gpd.GeoDataFrame, min_angular_acc: float = 0.0825, min_confinement: float = 1, min_angle: float = 0.78, min_force: float = 0.166) -> list`

##### Parameters
  - **`forces` (`gpd.GeoDataFrame`)**: GeoDataFrame outputted by `calc_forces()` with `force`, `confinement`, and `angle` columns.
  - **`min_force` (`float`, optional)**: Significance threshold for the resultant force. Default: `0.166`. (E.g., for a square building with height 1 and side length 1, if a touching structure covers only 1/6 of one side, the resultant force would be 1/6.)
  - **`min_angle` (`float`, optional)**: Angle threshold (in radians). Default: $\pi / 4$ (45 degrees).
  - **`min_confinement` (`float`, optional)**: Threshold for confinement. Default: `1` (indicating equal amounts of confined and resultant forces).
  - **`min_angular_acc` (`float`, optional)**: Threshold for angular acceleration $\frac{momentum * area}{inertia}$. Default: 2.133 (e.g., for a rectangular building with height 1 and sides of length 1 and 0.5, 
                            a touching structure covering 1/3 of two sides in the worst case would have an anuglar acceleration of 2.133)

##### Output
  Returns a list in the same ordes as forces rows with the `relative_position`. Classifies buildings into the following categories (priority order):
  1. **"torque"**: Buildings of class **confined** or **corner** with an angular acceleration exceeding the minimum.
  2. **"confined"**: Structures touching on both the left and right lateral sides.
  3. **"corner"**: Structures touching at a corner (determined by force and angle thresholds).
  4. **"partial"**: Structures touching on either the left or right side.
  5. **"isolated"**: No touching structures.

---

### 2. **Shape Irregularity**
- **Polsby-Popper Index**: A measure of shape compactness (similarity to a circle).
  - **Formula**: $\text{Polsby-Popper Index} = \frac{4 \pi A}{P^2}$
    where:
    - $A$: Area of the polygon.
    - $P$: Perimeter of the polygon.
  - #### Function: `calc_polsby_popper(geoms: gpd.GeoDataFrame) -> list`
    - **Parameters**:
      - `geoms`: GeoDataFrame with building footprint polygon geometries.
    - **Output**: A list in the same order as `geoms` rows with the `polsby_popper` index.

- **Custom Irregularity Index**: Our own index to quantify the irregularity of building footprints.
  - **Formula**: $\frac{l \cdot d}{L}$, where:
    - $l$: Length of the shapes outside the convex hull.
    - $d$: Distance of the center of gravity of the shapes outside the hull to the hull.
    - $L$: Total length of the convex hull.
  - #### Function: `calc_shape_irregularity(geoms: gpd.GeoDataFrame) -> list`
    - **Parameters**:
      - `geoms`: GeoDataFrame with building footprint polygon geometries.
    - **Output**: A list in the same order as `geoms` rows with the `shape_irregularity` index. 
   
- **Inertia Irregularity**: Inertia of a circle with the same area as the polygon geometry divided by the inertia of the polygon.
    - **Formula**: $\text{intertia irregularity} = \frac{\text{intertia eq circle}}{\text{itertia}}$
    - #### Function: `calc_inertia_irregularity(geoms:gpd.GeoDataFrame) -> list`
        - **Parameters**:
          - `geoms`: GeoDataFrame with building footprint polygon geometries.
        - **Output**: A list in the same order as `geoms` rows with the `inertia_irregularity` index.



