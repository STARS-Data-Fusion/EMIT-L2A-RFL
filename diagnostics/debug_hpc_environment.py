#!/usr/bin/env python
"""
Diagnostic script to debug HPC environment issues with NetCDF file access.

This script downloads a sample EMIT L2A RFL file and runs comprehensive diagnostics
to identify any issues with the HPC environment, file access, or library configuration.

Usage:
    # Download a sample file and run diagnostics (recommended for first-time use):
    python debug_hpc_environment.py
    
    # Run diagnostics on an existing file:
    python debug_hpc_environment.py <path_to_netcdf_file>
    
Examples:
    # Download and test (will authenticate with NASA Earthdata):
    python debug_hpc_environment.py
    
    # Test an existing file:
    python debug_hpc_environment.py ~/data/EMIT_L2A_RFL_001_20230129T004447_2302816_003.nc
    
The script will:
    1. Download a sample EMIT L2A RFL granule to ~/data (if no file provided)
    2. Check your environment configuration (library versions, environment variables)
    3. Test file access with different methods (direct netCDF4, with workarounds, via EMITL2ARFL)
    4. Provide specific recommendations for any issues found
"""

import os
import sys
from pathlib import Path
from datetime import datetime

def print_section(title):
    """Print a formatted section header."""
    print(f"\n{'='*60}")
    print(f" {title}")
    print('='*60)

def check_file_details(filename):
    """Check basic file system details."""
    print_section("FILE SYSTEM DETAILS")
    
    filepath = Path(os.path.expanduser(filename))
    print(f"Original path: {filename}")
    print(f"Expanded path: {filepath}")
    print(f"Absolute path: {filepath.resolve()}")
    print(f"File exists: {filepath.exists()}")
    
    if filepath.exists():
        stat = filepath.stat()
        print(f"File size: {stat.st_size:,} bytes")
        print(f"File mode: {oct(stat.st_mode)}")
        print(f"Is readable: {os.access(filepath, os.R_OK)}")
        print(f"Is writable: {os.access(filepath, os.W_OK)}")
        
        # Check if it's on a network filesystem
        try:
            import subprocess
            result = subprocess.run(['df', '-T', str(filepath)], 
                                  capture_output=True, text=True, timeout=5)
            if result.returncode == 0:
                print(f"\nFilesystem info:")
                print(result.stdout)
        except Exception as e:
            print(f"Could not check filesystem type: {e}")

def check_library_versions():
    """Check versions of critical libraries."""
    print_section("LIBRARY VERSIONS")
    
    try:
        import netCDF4
        print(f"netCDF4 version: {netCDF4.__version__}")
        print(f"netCDF4 file: {netCDF4.__file__}")
    except ImportError as e:
        print(f"ERROR: Cannot import netCDF4: {e}")
        return False
    
    try:
        import h5py
        print(f"h5py version: {h5py.__version__}")
        print(f"h5py file: {h5py.__file__}")
        print(f"HDF5 version: {h5py.version.hdf5_version}")
    except ImportError:
        print("h5py not available (optional)")
    
    try:
        import numpy
        print(f"numpy version: {numpy.__version__}")
    except ImportError:
        print("numpy not available")
    
    return True

def check_environment_variables():
    """Check relevant environment variables."""
    print_section("ENVIRONMENT VARIABLES")
    
    env_vars = [
        'HDF5_USE_FILE_LOCKING',
        'HDF5_PLUGIN_PATH',
        'NETCDF_PLUGIN_DIR',
        'LD_LIBRARY_PATH',
        'DYLD_LIBRARY_PATH',
        'PYTHONPATH',
    ]
    
    for var in env_vars:
        value = os.environ.get(var, '<not set>')
        print(f"{var}: {value}")

def test_direct_netcdf4_access(filename):
    """Test direct netCDF4.Dataset access."""
    print_section("DIRECT NETCDF4 ACCESS TEST")
    
    import netCDF4
    filepath = Path(os.path.expanduser(filename)).resolve()
    
    print(f"Attempting to open: {filepath}")
    
    try:
        # Try opening with explicit path conversion
        print("\nTest 1: Opening with resolved absolute path...")
        with netCDF4.Dataset(str(filepath), 'r') as ds:
            print(f"✓ SUCCESS: File opened")
            print(f"  Dimensions: {len(ds.dimensions)}")
            print(f"  Variables: {len(ds.variables)}")
            if ds.dimensions:
                print(f"  Dimension names: {list(ds.dimensions.keys())[:5]}")
            if ds.variables:
                print(f"  Variable names: {list(ds.variables.keys())[:5]}")
        return True
        
    except Exception as e:
        print(f"✗ FAILED: {type(e).__name__}: {e}")
        print(f"\nError details:")
        print(f"  Error type: {type(e)}")
        print(f"  Error args: {e.args}")
        
        # Check for specific error patterns
        error_str = str(e)
        if "errno -101" in error_str.lower():
            print(f"\n⚠ HDF Error -101 detected!")
            print(f"  This typically indicates:")
            print(f"  - File corruption")
            print(f"  - Network filesystem caching issues")
            print(f"  - HDF5 library incompatibility")
        
        return False

def test_with_file_locking_disabled(filename):
    """Test with HDF5 file locking disabled."""
    print_section("TEST WITH FILE LOCKING DISABLED")
    
    # Save original value
    original_value = os.environ.get('HDF5_USE_FILE_LOCKING')
    
    try:
        # Disable file locking
        os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
        print("Set HDF5_USE_FILE_LOCKING=FALSE")
        
        # Force reload of netCDF4 to pick up new environment
        import netCDF4
        if 'netCDF4' in sys.modules:
            print("Reloading netCDF4 module...")
            import importlib
            importlib.reload(netCDF4)
        
        filepath = Path(os.path.expanduser(filename)).resolve()
        print(f"Attempting to open: {filepath}")
        
        with netCDF4.Dataset(str(filepath), 'r') as ds:
            print(f"✓ SUCCESS with HDF5_USE_FILE_LOCKING=FALSE")
            return True
            
    except Exception as e:
        print(f"✗ FAILED even with file locking disabled: {e}")
        return False
        
    finally:
        # Restore original value
        if original_value is None:
            os.environ.pop('HDF5_USE_FILE_LOCKING', None)
        else:
            os.environ['HDF5_USE_FILE_LOCKING'] = original_value

def test_validate_function(filename):
    """Test the actual validate_NetCDF_file function."""
    print_section("VALIDATE_NETCDF_FILE FUNCTION TEST")
    
    try:
        from EMITL2ARFL import validate_NetCDF_file
        print(f"Attempting validation with validate_NetCDF_file()...")
        
        validate_NetCDF_file(filename)
        print(f"✓ SUCCESS: File passed validation")
        return True
        
    except Exception as e:
        print(f"✗ FAILED: {type(e).__name__}: {e}")
        
        # Show the full traceback
        import traceback
        print(f"\nFull traceback:")
        traceback.print_exc()
        return False

def download_sample_granule():
    """Download a sample EMIT L2A RFL granule for testing."""
    print_section("DOWNLOADING SAMPLE GRANULE")
    
    download_dir = Path.home() / "data"
    download_dir.mkdir(exist_ok=True)
    
    print(f"Download directory: {download_dir}")
    print()
    
    try:
        print("Authenticating with NASA Earthdata...")
        import earthaccess
        earthaccess.login()
        print("✓ Authentication successful")
        print()
        
        # Search for a specific known granule (small file from early 2023)
        print("Searching for sample granule...")
        print("  Looking for: EMIT_L2A_RFL_001_20230129T004447_2302816_003")
        
        from EMITL2ARFL import search_EMIT_L2A_RFL_granules
        
        # First try to search for the specific granule by orbit and scene
        # This granule is: orbit=2302816, scene=003
        granules = search_EMIT_L2A_RFL_granules(
            start_UTC=datetime(2023, 1, 29),
            end_UTC=datetime(2023, 1, 30),
            orbit=2302816,
            scene=3
        )
        
        if not granules:
            print("✗ Specific granule not found")
            print("  Trying to find any granule from January 2023...")
            granules = search_EMIT_L2A_RFL_granules(
                start_UTC=datetime(2023, 1, 20),
                end_UTC=datetime(2023, 1, 31)
            )
            if granules:
                print(f"  Found {len(granules)} granules, using first one")
            else:
                print("✗ No granules found")
                return None
        else:
            print(f"✓ Found granule")
        
        print()
        print("Downloading granule (this may take a few minutes)...")
        print(f"  Granule ID: {granules[0]['umm']['GranuleUR']}")
        print()
        print("Note: File validation is disabled during download to prevent")
        print("      automatic deletion. The diagnostic tests will validate")
        print("      the file after download completes.")
        print()
        
        from EMITL2ARFL import retrieve_EMIT_L2A_RFL_granule
        
        granule = retrieve_EMIT_L2A_RFL_granule(
            remote_granule=granules[0],
            download_directory=str(download_dir),
            max_retries=1,  # Only try once since validation is disabled
            skip_validation=True  # Don't validate during download - we'll test it after
        )
        
        print("✓ Download successful")
        print(f"  Reflectance file: {granule.reflectance_filename}")
        print(f"  File will remain in: {download_dir}")
        print()
        
        return granule.reflectance_filename
        
    except ImportError as e:
        print(f"✗ Cannot import required packages: {e}")
        print("  Make sure earthaccess and EMITL2ARFL are installed")
        return None
    except Exception as e:
        print(f"✗ Download failed: {e}")
        import traceback
        traceback.print_exc()
        return None

def main():
    """Main diagnostic routine."""
    print("="*60)
    print(" HPC NetCDF Environment Diagnostic Tool")
    print("="*60)
    print()
    print("This tool will:")
    print("  1. Download a sample EMIT L2A RFL file (if not provided)")
    print("  2. Check your environment configuration")
    print("  3. Test file access with different methods")
    print("  4. Provide specific recommendations")
    print()
    
    # Determine which file to test
    if len(sys.argv) >= 2:
        filename = sys.argv[1]
        print(f"Using provided file: {filename}")
        print()
    else:
        print("No file provided, will download a sample granule...")
        print()
        filename = download_sample_granule()
        
        if filename is None:
            print()
            print("="*60)
            print("Could not download sample file.")
            print("Please provide a file path manually:")
            print("  python debug_hpc_environment.py <path_to_netcdf_file>")
            print("="*60)
            sys.exit(1)
    
    # Run all diagnostic checks
    check_file_details(filename)
    
    if not check_library_versions():
        print("\nERROR: Required libraries not available. Exiting.")
        sys.exit(1)
    
    check_environment_variables()
    
    # Try different access methods
    test1 = test_direct_netcdf4_access(filename)
    test2 = test_with_file_locking_disabled(filename)
    test3 = test_validate_function(filename)
    
    # Summary
    print_section("SUMMARY")
    print(f"Direct netCDF4 access: {'✓ PASS' if test1 else '✗ FAIL'}")
    print(f"With file locking disabled: {'✓ PASS' if test2 else '✗ FAIL'}")
    print(f"validate_NetCDF_file function: {'✓ PASS' if test3 else '✗ FAIL'}")
    
    print()
    print("Note: As of the latest version, EMITL2ARFL automatically sets")
    print("HDF5_USE_FILE_LOCKING=FALSE during import. The package should")
    print("work without any manual environment configuration.")
    
    if test1 and test2 and test3:
        print()
        print("="*60)
        print("🎉 ALL TESTS PASSED!")
        print("="*60)
        print()
        print("Your environment is properly configured and the file is valid.")
        print("You should be able to use EMITL2ARFL without any issues.")
        return
    
    if not test1 and test2:
        print()
        print("\n⚠ RECOMMENDATION:")
        print("  The file opens successfully with HDF5_USE_FILE_LOCKING=FALSE")
        print("  This is a common issue on NFS/Lustre filesystems.")
        print()
        print("  The EMITL2ARFL package should handle this automatically.")
        print("  If you're still experiencing issues:")
        print()
        print("  1. Make sure you have the latest version of EMITL2ARFL")
        print("  2. Import EMITL2ARFL before any other packages that use netCDF4")
        print("  3. Or set the environment variable manually:")
        print("     export HDF5_USE_FILE_LOCKING=FALSE")
        print()
        print("  For more help, see: HPC_QUICKSTART.md")
    
    if not test1 and not test2:
        print("\n⚠ RECOMMENDATION:")
        print("  The file cannot be opened even with file locking disabled.")
        print("  This suggests:")
        print("    1. File may be corrupted")
        print("    2. HDF5/NetCDF library incompatibility")
        print("    3. Filesystem-level issues")
        print("  Try:")
        print("    - Re-downloading the file")
        print("    - Checking md5/checksum integrity")
        print("    - Verifying HDF5 library versions match")

if __name__ == "__main__":
    main()
