from ._forces import calc_forces, relative_position
from ._irregularity import shape_irregularity, polsby_popper, inertia_irregularity
import sys
__all__ = ['calc_forces', 'relative_position', 'shape_irregularity', 'polsby_popper', 'inertia_irregularity']

# Remove the private submodules from the package namespace
del sys.modules[__name__ + '._forces']
del sys.modules[__name__ + '._irregularity']
del sys.modules[__name__ + '._utils']
