#!/usr/bin/env python
"""
Test that the package automatically handles HDF5 file locking on HPC systems.

This script verifies that:
1. The package sets HDF5_USE_FILE_LOCKING=FALSE automatically
2. NetCDF files can be opened without manual environment configuration
3. The fix happens transparently without user intervention
"""

import sys
import os

def test_automatic_fix():
    """Test that the package sets HDF5_USE_FILE_LOCKING automatically."""
    
    print("="*60)
    print("Testing EMITL2ARFL Automatic HPC Fix")
    print("="*60)
    print()
    
    # Check initial state
    initial_value = os.environ.get('HDF5_USE_FILE_LOCKING')
    print(f"1. Initial HDF5_USE_FILE_LOCKING: {initial_value or '<not set>'}")
    
    # Verify netCDF4 is not loaded yet
    netcdf4_preloaded = 'netCDF4' in sys.modules
    print(f"2. netCDF4 pre-loaded: {netcdf4_preloaded}")
    
    if netcdf4_preloaded:
        print("   ⚠ WARNING: netCDF4 was already imported!")
        print("   The automatic fix may not work correctly.")
        print()
    
    # Import EMITL2ARFL
    print("3. Importing EMITL2ARFL...")
    from EMITL2ARFL import validate_NetCDF_file
    print("   ✓ Import successful")
    
    # Check if environment variable was set
    final_value = os.environ.get('HDF5_USE_FILE_LOCKING')
    print(f"4. Final HDF5_USE_FILE_LOCKING: {final_value}")
    
    # Verify netCDF4 is now loaded
    netcdf4_loaded = 'netCDF4' in sys.modules
    print(f"5. netCDF4 loaded: {netcdf4_loaded}")
    print()
    
    # Test with a file if provided
    if len(sys.argv) > 1:
        test_file = sys.argv[1]
        print(f"6. Testing file validation: {test_file}")
        try:
            validate_NetCDF_file(test_file)
            print("   ✓ File validation successful!")
        except Exception as e:
            print(f"   ✗ File validation failed: {e}")
            return False
        print()
    
    # Summary
    print("="*60)
    print("Summary:")
    print("="*60)
    
    if initial_value and initial_value.upper() != 'FALSE':
        print(f"✓ User had set HDF5_USE_FILE_LOCKING={initial_value}")
        print("  Package respected user's setting")
        success = True
    elif initial_value and initial_value.upper() == 'FALSE':
        print("✓ User had already set HDF5_USE_FILE_LOCKING=FALSE")
        print("  Package maintained the setting")
        success = True
    elif final_value and final_value.upper() == 'FALSE':
        print("✓ Package automatically set HDF5_USE_FILE_LOCKING=FALSE")
        print("  No user configuration was required!")
        success = True
    else:
        print("✗ HDF5_USE_FILE_LOCKING was not properly set")
        print("  Expected: FALSE")
        print(f"  Got: {final_value}")
        success = False
    
    print()
    if success:
        print("🎉 SUCCESS: Automatic HPC fix is working!")
        print()
        print("This package will work on HPC systems without requiring users")
        print("to manually set environment variables. Just import and use!")
        return True
    else:
        print("✗ FAILURE: Automatic fix did not work as expected")
        return False

if __name__ == "__main__":
    success = test_automatic_fix()
    sys.exit(0 if success else 1)
