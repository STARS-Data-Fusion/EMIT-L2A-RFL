#!/usr/bin/env python
"""
Check HDF5 and NetCDF library versions and compatibility.
"""

import sys
import os

print("=" * 60)
print("HDF5/NetCDF Library Diagnostic")
print("=" * 60)
print()

# Check environment
print("Environment Variables:")
print(f"  HDF5_USE_FILE_LOCKING: {os.environ.get('HDF5_USE_FILE_LOCKING', '<not set>')}")
print()

# Check h5py
try:
    import h5py
    print("h5py:")
    print(f"  Version: {h5py.__version__}")
    print(f"  HDF5 Version: {h5py.version.hdf5_version}")
    print(f"  API Version: {h5py.version.api_version}")
    print(f"  Location: {h5py.__file__}")
    print()
except ImportError as e:
    print(f"h5py: Not installed ({e})")
    print()

# Check netCDF4
try:
    import netCDF4
    print("netCDF4:")
    print(f"  Version: {netCDF4.__version__}")
    print(f"  Location: {netCDF4.__file__}")
    
    # Get underlying HDF5 version from netCDF4
    try:
        import netCDF4._netCDF4 as nc4
        if hasattr(nc4, '__hdf5libversion__'):
            print(f"  HDF5 lib version: {nc4.__hdf5libversion__}")
        if hasattr(nc4, '__netcdf4libversion__'):
            print(f"  NetCDF lib version: {nc4.__netcdf4libversion__}")
    except Exception:
        pass
    print()
except ImportError as e:
    print(f"netCDF4: Not installed ({e})")
    print()

# Check Python version
print("Python:")
print(f"  Version: {sys.version}")
print(f"  Executable: {sys.executable}")
print()

# Test file opening with h5py directly
if len(sys.argv) > 1:
    test_file = sys.argv[1]
    print("=" * 60)
    print(f"Testing file: {test_file}")
    print("=" * 60)
    print()
    
    # Check file exists and size
    if os.path.exists(test_file):
        size = os.path.getsize(test_file)
        print(f"File exists: {test_file}")
        print(f"File size: {size:,} bytes ({size/1024/1024:.1f} MB)")
        print()
        
        # Try h5py first (lower level)
        print("Test 1: Open with h5py...")
        try:
            import h5py
            with h5py.File(test_file, 'r') as f:
                print("  ✓ SUCCESS with h5py")
                print(f"  Keys: {list(f.keys())[:5]}")
        except Exception as e:
            print(f"  ✗ FAILED with h5py: {type(e).__name__}: {e}")
        print()
        
        # Try netCDF4 (higher level)
        print("Test 2: Open with netCDF4...")
        try:
            import netCDF4
            with netCDF4.Dataset(test_file, 'r') as ds:
                print("  ✓ SUCCESS with netCDF4")
                print(f"  Dimensions: {list(ds.dimensions.keys())[:5]}")
                print(f"  Variables: {list(ds.variables.keys())[:5]}")
        except Exception as e:
            print(f"  ✗ FAILED with netCDF4: {type(e).__name__}: {e}")
        print()
    else:
        print(f"File not found: {test_file}")
        print()

print("=" * 60)
print("Recommendations:")
print("=" * 60)
print()
print("If h5py works but netCDF4 fails:")
print("  - NetCDF4 library incompatibility")
print("  - Try: conda install -c conda-forge netcdf4")
print()
print("If both h5py and netCDF4 fail:")
print("  - HDF5 library incompatibility")
print("  - Try: conda install -c conda-forge hdf5 h5py netcdf4")
print()
print("If file opens on Mac but not HPC:")
print("  - Different HDF5 versions between systems")
print("  - Re-download file directly on HPC system")
print("  - Or rebuild conda environment with matching library versions")
