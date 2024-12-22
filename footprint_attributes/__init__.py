from .forces import calc_forces, relative_position
from .irregularity import calc_shape_irregularity, calc_polsby_popper

__all__ = ['calc_forces', 'relative_position', 'calc_shape_irregularity', 'calc_polsby_popper']

# Explicitly hide the 'forces' module to prevent it from being directly imported
import sys
import warnings

# Delete 'forces' from the globals to ensure it cannot be accessed
del sys.modules['.forces']
del sys.modules['.irregularity']
del sys.modules['.utils']