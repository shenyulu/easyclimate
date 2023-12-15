# from importlib.metadata import version, PackageNotFoundError

from .core import *
from .interp import *
from .filter import *

from . import windspharm
from . import interp
from . import plot
from . import wavelet
from . import field

# Version number
__version__ = "2023.12.01"

# try:
#     __version__ = version("package-name")
# except PackageNotFoundError:
#     # package is not installed
#     pass