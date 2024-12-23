from ._forces import calc_forces, relative_position
from ._irregularity import calc_shape_irregularity, calc_polsby_popper
import sys
__all__ = ['calc_forces', 'relative_position', 'calc_shape_irregularity', 'calc_polsby_popper']

# Remove the private submodules from the package namespace
del sys.modules[__name__ + '._forces']
del sys.modules[__name__ + '._irregularity']
del sys.modules[__name__ + '._utils']
