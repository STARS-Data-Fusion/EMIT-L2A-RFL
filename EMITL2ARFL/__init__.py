"""
EMIT L2A Reflectance Tools

This package provides tools for working with NASA EMIT L2A reflectance data.

HPC Compatibility:
    The package automatically sets HDF5_USE_FILE_LOCKING=FALSE during import
    to prevent file locking issues on network filesystems (NFS/Lustre/GPFS)
    commonly found on HPC systems. This happens automatically - no user
    configuration needed. If you've already set this environment variable,
    your setting will be respected.
"""

import os
import sys

# CRITICAL: Set HDF5_USE_FILE_LOCKING before any netCDF4/HDF5 imports
# This prevents file locking issues on network filesystems (NFS/Lustre/GPFS)
# commonly found on HPC systems.
if 'HDF5_USE_FILE_LOCKING' not in os.environ:
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

# Check if netCDF4 was already imported (it shouldn't be at this point)
if 'netCDF4' in sys.modules:
    # If it was already imported, try to reload it with the new setting
    import importlib
    import warnings
    warnings.warn(
        "netCDF4 was already imported before EMITL2ARFL initialization. "
        "HDF5_USE_FILE_LOCKING setting may not take effect. "
        "For best results, import EMITL2ARFL before any packages that use netCDF4.",
        RuntimeWarning
    )
    try:
        importlib.reload(sys.modules['netCDF4'])
    except Exception:
        pass

# Now safe to import all submodules
from .EMITL2ARFL import *
from .version import __version__

