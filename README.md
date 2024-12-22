# Footprint Attributes

**Footprint Attributes** is a Python package built using GeoPandas to calculate various geometric indices related to building footprint geometry and seismic risk. 

## Features

The package provides functions for:

1. **Shape Irregularity**:
   - **Polsby-Popper Index**: A measure of shape compactness.
   - **Custom Irregularity Index**: Our own index to quantify the irregularity of building footprints.

2. **Relative Position of Buildings**:
   - **Confined**: Buildings with structures touching on both the left and right lateral sides.
   - **Lateral**: Buildings with a structure on either the left or right side.
   - **Corner**: Buildings positioned at a corner.
