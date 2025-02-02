from ._forces import get_forces_gdf, relative_position
from ._irregularity import shape_irregularity, polsby_popper, inertia_irregularity
import sys

__version__ = '0.2.0'
__all__ = ['get_forces_gdf', 'relative_position', 'shape_irregularity', 'polsby_popper', 'inertia_circle','inertia_slenderness','excentricity_EC_8','__version__']


# Remove the private submodules from the package namespace
del sys.modules[__name__ + '._forces']
del sys.modules[__name__ + '._irregularity']
del sys.modules[__name__ + '._utils']
