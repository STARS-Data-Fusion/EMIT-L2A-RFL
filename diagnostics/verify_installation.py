#!/usr/bin/env python
"""
Verify that EMITL2ARFL is using the correct HDF5/NetCDF libraries.
This ensures no conflicts with system libraries on HPC systems.
"""

import sys
import os
from pathlib import Path

print("=" * 70)
print("EMITL2ARFL Installation Verification")
print("=" * 70)
print()

# Check Python location
print("Python:")
print(f"  Executable: {sys.executable}")
print(f"  Version: {sys.version.split()[0]}")

# Check if in conda environment
if 'CONDA_PREFIX' in os.environ:
    print(f"  Conda env: {os.environ['CONDA_PREFIX']}")
elif 'VIRTUAL_ENV' in os.environ:
    print(f"  Venv: {os.environ['VIRTUAL_ENV']}")
else:
    print("  ⚠️  WARNING: Not in conda env or virtualenv")
print()

# Check critical environment variable
print("Environment:")
hdf5_locking = os.environ.get('HDF5_USE_FILE_LOCKING', '<not set>')
if hdf5_locking == 'FALSE':
    print(f"  ✓ HDF5_USE_FILE_LOCKING: {hdf5_locking}")
else:
    print(f"  ⚠️  HDF5_USE_FILE_LOCKING: {hdf5_locking}")
    print(f"     Should be 'FALSE' for HPC systems")
print()

# Check h5py
print("h5py:")
try:
    import h5py
    h5py_path = Path(h5py.__file__)
    python_path = Path(sys.executable).parent.parent
    
    is_in_env = str(python_path) in str(h5py_path)
    status = "✓" if is_in_env else "⚠️ "
    
    print(f"  {status} Version: {h5py.__version__}")
    print(f"  {status} HDF5 lib: {h5py.version.hdf5_version}")
    print(f"  {status} Location: {h5py_path}")
    
    if not is_in_env:
        print(f"     WARNING: h5py not from your Python environment!")
        print(f"     Expected path containing: {python_path}")
except ImportError as e:
    print(f"  ✗ Not installed: {e}")
print()

# Check netCDF4
print("netCDF4:")
try:
    import netCDF4
    nc4_path = Path(netCDF4.__file__)
    python_path = Path(sys.executable).parent.parent
    
    is_in_env = str(python_path) in str(nc4_path)
    status = "✓" if is_in_env else "⚠️ "
    
    print(f"  {status} Version: {netCDF4.__version__}")
    print(f"  {status} Location: {nc4_path}")
    
    # Try to get HDF5 version from netCDF4
    try:
        import netCDF4._netCDF4 as nc4
        if hasattr(nc4, '__hdf5libversion__'):
            print(f"  {status} HDF5 lib: {nc4.__hdf5libversion__}")
        if hasattr(nc4, '__netcdf4libversion__'):
            print(f"  {status} NetCDF lib: {nc4.__netcdf4libversion__}")
    except Exception:
        pass
    
    if not is_in_env:
        print(f"     WARNING: netCDF4 not from your Python environment!")
        print(f"     Expected path containing: {python_path}")
except ImportError as e:
    print(f"  ✗ Not installed: {e}")
print()

# Check EMITL2ARFL
print("EMITL2ARFL:")
try:
    import EMITL2ARFL
    emit_path = Path(EMITL2ARFL.__file__)
    
    print(f"  ✓ Version: {EMITL2ARFL.__version__}")
    print(f"  ✓ Location: {emit_path}")
    
    # Check if it's editable install
    if 'site-packages' not in str(emit_path):
        print(f"     (editable/development install)")
except ImportError as e:
    print(f"  ✗ Not installed: {e}")
print()

# Summary
print("=" * 70)
print("Summary:")
print("=" * 70)

all_good = True

# Check 1: In environment
if 'CONDA_PREFIX' not in os.environ and 'VIRTUAL_ENV' not in os.environ:
    print("✗ Not in conda environment or virtualenv")
    print("  Recommendation: conda create -n EMITL2ARFL python=3.10")
    all_good = False

# Check 2: HDF5 locking
if hdf5_locking != 'FALSE':
    print("✗ HDF5_USE_FILE_LOCKING not set to FALSE")
    print("  Recommendation: set -Ux HDF5_USE_FILE_LOCKING FALSE  (fish)")
    print("  Or: export HDF5_USE_FILE_LOCKING=FALSE  (bash, in ~/.bashrc)")
    all_good = False

# Check 3: Libraries installed
try:
    import h5py
    import netCDF4
    import EMITL2ARFL
    
    # Check if from environment
    h5py_path = Path(h5py.__file__)
    nc4_path = Path(netCDF4.__file__)
    python_path = Path(sys.executable).parent.parent
    
    if str(python_path) not in str(h5py_path):
        print("✗ h5py not from your Python environment")
        print("  Recommendation: conda install -c conda-forge h5py")
        all_good = False
    
    if str(python_path) not in str(nc4_path):
        print("✗ netCDF4 not from your Python environment")
        print("  Recommendation: conda install -c conda-forge netcdf4")
        all_good = False
        
except ImportError as e:
    print(f"✗ Missing required package: {e}")
    print("  Recommendation: Follow INSTALL_HPC.md instructions")
    all_good = False

if all_good:
    print("✓ Installation looks good!")
    print()
    print("Next steps:")
    print("  1. Test with: python check_hdf5_version.py <netcdf_file>")
    print("  2. Or run: python test_single_granule.py")
else:
    print()
    print("⚠️  Installation has issues - see recommendations above")
    print()
    print("For complete setup instructions, see: INSTALL_HPC.md")

print()
