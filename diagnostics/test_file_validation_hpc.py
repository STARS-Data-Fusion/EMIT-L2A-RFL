#!/usr/bin/env python
"""
Test file validation with HPC-specific workarounds.

This version handles HPC-specific issues, particularly HDF5 file locking
on network filesystems (NFS, Lustre, GPFS).

IMPORTANT: For best results on HPC systems, set the environment variable
BEFORE running the script:

    export HDF5_USE_FILE_LOCKING=FALSE
    python test_file_validation_hpc.py

Or in a single command:

    HDF5_USE_FILE_LOCKING=FALSE python test_file_validation_hpc.py

If you must set it within Python, this script attempts to do so, but
it may not work if netCDF4 has already been imported elsewhere.
"""

import os
import sys
from pathlib import Path

# Check if environment variable is already set
env_var_set = os.environ.get('HDF5_USE_FILE_LOCKING', '').upper() == 'FALSE'

if not env_var_set:
    print("⚠ WARNING: HDF5_USE_FILE_LOCKING is not set to FALSE")
    print("   This may cause issues on network filesystems (NFS/Lustre)")
    print("   Attempting to set it now, but it's better to set it before running:")
    print("   $ export HDF5_USE_FILE_LOCKING=FALSE")
    print("   $ python test_file_validation_hpc.py")
    print()
    
    # Try to set it (may be too late if modules already loaded)
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    print("Setting HDF5_USE_FILE_LOCKING=FALSE...")
else:
    print("✓ HDF5_USE_FILE_LOCKING is already set to FALSE")

# Check if netCDF4 is already loaded
if 'netCDF4' in sys.modules:
    print("⚠ WARNING: netCDF4 was already imported before this script ran")
    print("   The HDF5_USE_FILE_LOCKING setting may not take effect")
    print()

# Now import EMITL2ARFL
from EMITL2ARFL import validate_NetCDF_file

filename = "~/data/EMIT_L2A_RFL_001_20230129T004447_2302816_003.nc"

# WORKAROUND 2: Ensure path is fully resolved
filepath = Path(os.path.expanduser(filename)).resolve()
print(f"\nValidating file: {filepath}")
print(f"File exists: {filepath.exists()}")
if filepath.exists():
    print(f"File size: {filepath.stat().st_size:,} bytes")
else:
    print(f"ERROR: File does not exist at {filepath}")
    print(f"Please update the 'filename' variable in this script to point to a valid NetCDF file")
    sys.exit(1)

print("\nAttempting validation...")
try:
    validate_NetCDF_file(str(filepath))
    print(f"✓ SUCCESS: Valid NetCDF file")
except Exception as e:
    print(f"✗ FAILED: Validation error")
    print(f"  Error type: {type(e).__name__}")
    print(f"  Error message: {e}")
    print()
    
    # Additional debugging - try direct access
    print("Attempting direct netCDF4 access for debugging...")
    try:
        import netCDF4
        print(f"  netCDF4 version: {netCDF4.__version__}")
        print(f"  Attempting to open file...")
        with netCDF4.Dataset(str(filepath), 'r') as ds:
            print(f"  ✓ Direct netCDF4 access WORKS!")
            print(f"    Dimensions: {len(ds.dimensions)}, Variables: {len(ds.variables)}")
            print()
            print("  This suggests the issue is in the validate_NetCDF_file function,")
            print("  not with HDF5 file locking.")
    except Exception as e2:
        print(f"  ✗ Direct netCDF4 access also FAILS: {e2}")
        print()
        print("  RECOMMENDATION:")
        print("  1. Make sure HDF5_USE_FILE_LOCKING=FALSE is set BEFORE running Python:")
        print("     $ export HDF5_USE_FILE_LOCKING=FALSE")
        print("     $ python test_file_validation_hpc.py")
        print()
        print("  2. Or run with the environment variable in the command:")
        print("     $ HDF5_USE_FILE_LOCKING=FALSE python test_file_validation_hpc.py")
        print()
        print("  3. Add to your ~/.bashrc for permanent effect:")
        print("     export HDF5_USE_FILE_LOCKING=FALSE")
    
    sys.exit(1)
