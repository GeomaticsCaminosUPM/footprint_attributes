# Footprint Attributes

**Footprint Attributes** is a Python package built using GeoPandas to calculate various geometric indices related to building footprint geometry and seismic risk.

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

---

### 2. **Relative Position of Buildings**
Determines if the building touches other structures (relative position in the city block). This is done by calculating "forces" that neighboring structures exert on the building. The force is proportional to the contact area (length of touching footprints multiplied by building height) in the normal direction of the touching plane.

#### Function: `calc_forces(geoms: gpd.GeoDataFrame, buffer: float = 0, height_column: str = None) -> gpd.GeoDataFrame`
- **Parameters**:
  - `geoms`: GeoDataFrame containing building footprints as polygon geometries.
  - `buffer`: Buffer in meters to consider when determining if two buildings are touching each other.
  - `height_column`: Column in the `geoms` DataFrame specifying the height of the building in meters. If `None`, all buildings are assumed to have height 1.

- **Output**:
  Returns a GeoDataFrame with the following new columns:
  - **"momentum"**: Momentum of the resultant force with respect to the centroid of the footprint. $\sum(d*|\text{force}_i|)$
  - **"force"**: Magnitude of the sum of all forces on each footprint (resultant force). $|\sum(\text{force}_i)|$
  - **"confinement"**: A measure of the amount of force from the total forces that is confined (has an opposing force).  
    Formula: $\frac{\sum(|\text{force}_i|) - |\sum(\text{force}_i)|}{|\sum(\text{force}_i)|}$
  - **"angle"**: Sum of angles of forces with respect to the resultant force, normalized by force magnitude.  
    Formula: $\frac{\sum(|\text{force}_i| \cdot \text{angle}(\text{force}_i, \sum(\text{force}_j)))}{|\sum(\text{force}_i)|}$

#### Function: `relative_position(footprints: gpd.GeoDataFrame, min_momentum: float = 0.0825, min_confinement: float = 1, min_angle: float = 0.78, min_force: float = 0.166) -> gpd.GeoDataFrame`

- **Parameters**:
    - `footprints`: GeoDataFrame outputted by `calc_forces()` with `force`, `confinement`, and `angle` columns.
    - `min_force`: Significance threshold for the resultant force. Default: `0.166`. (E.g., for a square building with height 1 and side length 1, if a touching structure covers only 1/6 of one side, the resultant force would be 1/6.)
    - `min_angle`: Angle threshold (in radians). Default: $\pi / 4$ (45 degrees).
    - `min_confinement`: Threshold for confinement. Default: `1` (indicating equal amounts of confined and resultant forces).
    - `min_momentum`: Threshold for momentum. Default: `0.0825`. (E.g., for a square building with height 1 and side length 1, if a touching structure covers only 1/6 of two sides, in the worst case the momentum would be 0.0825.)

- **Output**:
  Returns the same GeoDataFrame with a new column `"relative_position"`. Classifies buildings into the following categories (priority order):
  1. **"torque"**: Buildings of class **confined** or **corner** with a momentum exceeding the minimum momentum.
  2. **"confined"**: Structures touching on both the left and right lateral sides.
  3. **"corner"**: Structures touching at a corner (determined by force and angle thresholds).
  4. **"partial"**: Structures touching on either the left or right side.
  5. **"isolated"**: No touching structures.




