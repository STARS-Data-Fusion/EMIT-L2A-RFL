#!/usr/bin/env python3
"""
Diagnostic script to download a single EMIT file and check all dependencies.
"""

import os
import sys
import subprocess
import hashlib
from pathlib import Path

def print_section(title):
    """Print a section header."""
    print(f"\n{'='*80}")
    print(f"  {title}")
    print(f"{'='*80}\n")

def check_python_environment():
    """Check Python environment details."""
    print_section("Python Environment")
    print(f"Python version: {sys.version}")
    print(f"Python executable: {sys.executable}")
    print(f"Platform: {sys.platform}")
    
    # Check conda environment
    conda_env = os.environ.get('CONDA_DEFAULT_ENV', 'Not in conda')
    print(f"Conda environment: {conda_env}")
    print(f"Virtual environment: {sys.prefix}")

def check_dependencies():
    """Check all required dependencies and their versions."""
    print_section("Dependency Versions")
    
    dependencies = [
        'netCDF4',
        'h5py',
        'earthaccess',
        'numpy',
        'requests'
    ]
    
    for dep in dependencies:
        try:
            module = __import__(dep)
            version = getattr(module, '__version__', 'unknown')
            print(f"✓ {dep}: {version}")
            
            # For netCDF4, check the underlying HDF5 library
            if dep == 'netCDF4':
                import netCDF4
                print(f"  - NetCDF library version: {netCDF4.__netcdf4libversion__}")
                print(f"  - HDF5 library version: {netCDF4.__hdf5libversion__}")
        except ImportError as e:
            print(f"✗ {dep}: NOT INSTALLED - {e}")
        except Exception as e:
            print(f"✗ {dep}: ERROR - {e}")

def check_filesystem(path):
    """Check filesystem properties."""
    print_section(f"Filesystem Check: {path}")
    
    path_obj = Path(path)
    
    # Check if path exists and is writable
    if not path_obj.exists():
        print(f"Creating directory: {path}")
        path_obj.mkdir(parents=True, exist_ok=True)
    
    print(f"Path exists: {path_obj.exists()}")
    print(f"Is directory: {path_obj.is_dir()}")
    print(f"Is writable: {os.access(path, os.W_OK)}")
    
    # Get filesystem type (Unix-like systems)
    try:
        if sys.platform != 'win32':
            result = subprocess.run(['df', '-T', str(path)], 
                                  capture_output=True, text=True, timeout=5)
            print(f"\nFilesystem info:\n{result.stdout}")
    except Exception as e:
        print(f"Could not determine filesystem type: {e}")
    
    # Check disk space
    try:
        stat = os.statvfs(path)
        free_space_gb = (stat.f_bavail * stat.f_frsize) / (1024**3)
        print(f"Free space: {free_space_gb:.2f} GB")
    except Exception as e:
        print(f"Could not determine free space: {e}")

def compute_checksum(filepath):
    """Compute MD5 checksum of a file."""
    md5 = hashlib.md5()
    try:
        with open(filepath, 'rb') as f:
            for chunk in iter(lambda: f.read(8192), b''):
                md5.update(chunk)
        return md5.hexdigest()
    except Exception as e:
        return f"ERROR: {e}"

def test_netcdf_operations(filepath):
    """Test various NetCDF operations on a file."""
    print_section(f"NetCDF Operations Test: {Path(filepath).name}")
    
    import netCDF4
    
    # Test 1: Check if file exists and is readable
    print(f"1. File exists: {os.path.exists(filepath)}")
    print(f"   File size: {os.path.getsize(filepath) / (1024**2):.2f} MB")
    print(f"   File checksum: {compute_checksum(filepath)}")
    
    # Test 2: Try to open with netCDF4
    print("\n2. Opening with netCDF4...")
    try:
        with netCDF4.Dataset(filepath, 'r') as ds:
            print("   ✓ Successfully opened")
            print(f"   - Format: {ds.data_model}")
            print(f"   - Dimensions: {list(ds.dimensions.keys())}")
            print(f"   - Variables: {list(ds.variables.keys())[:5]}...")  # First 5
    except Exception as e:
        print(f"   ✗ Failed to open: {type(e).__name__}: {e}")
        return False
    
    # Test 3: Try with h5py
    print("\n3. Opening with h5py...")
    try:
        import h5py
        with h5py.File(filepath, 'r') as f:
            print("   ✓ Successfully opened")
            print(f"   - Keys: {list(f.keys())[:5]}...")  # First 5
    except Exception as e:
        print(f"   ✗ Failed to open: {type(e).__name__}: {e}")
    
    # Test 4: Check file integrity with ncdump (if available)
    print("\n4. Checking with ncdump...")
    try:
        result = subprocess.run(['ncdump', '-h', filepath], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            print("   ✓ ncdump succeeded")
            # Print first few lines
            lines = result.stdout.split('\n')[:10]
            for line in lines:
                print(f"   {line}")
        else:
            print(f"   ✗ ncdump failed: {result.stderr}")
    except FileNotFoundError:
        print("   - ncdump not available (not installed)")
    except Exception as e:
        print(f"   - ncdump error: {e}")
    
    return True

def download_single_file(download_dir):
    """Download a single EMIT granule for testing."""
    print_section("Downloading Single EMIT File")
    
    import earthaccess
    
    # Authenticate
    print("Authenticating with NASA Earthdata...")
    try:
        auth = earthaccess.login()
        print("✓ Authentication successful")
    except Exception as e:
        print(f"✗ Authentication failed: {e}")
        return None
    
    # Search for a single granule
    print("\nSearching for a recent granule...")
    try:
        results = earthaccess.search_data(
            short_name='EMITL2ARFL',
            temporal=('2023-01-29', '2023-01-30'),
            count=1
        )
        
        if not results:
            print("✗ No granules found")
            return None
        
        granule = results[0]
        print(f"✓ Found granule: {granule['umm']['GranuleUR']}")
        
    except Exception as e:
        print(f"✗ Search failed: {e}")
        return None
    
    # Download with single thread
    print(f"\nDownloading to: {download_dir}")
    print("Using single-threaded download (threads=1)...")
    
    try:
        downloaded = earthaccess.download(
            granules=[granule],
            local_path=download_dir,
            threads=1  # Single-threaded
        )
        
        if downloaded:
            filepath = downloaded[0]
            print(f"✓ Download complete: {filepath}")
            
            # Sync filesystem
            print("\nSyncing filesystem...")
            os.sync()
            
            # Wait a moment
            import time
            time.sleep(1)
            
            return filepath
        else:
            print("✗ Download failed - no files returned")
            return None
            
    except Exception as e:
        print(f"✗ Download failed: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return None

def main():
    """Run all diagnostic checks."""
    print("\n" + "="*80)
    print("  EMIT NetCDF Download Diagnostic")
    print("="*80)
    
    # Get download directory from command line or use default
    if len(sys.argv) > 1:
        download_dir = sys.argv[1]
    else:
        download_dir = '/tmp/emit_diagnostic'
    
    print(f"\nDownload directory: {download_dir}")
    
    # Run checks
    check_python_environment()
    check_dependencies()
    check_filesystem(download_dir)
    
    # Download a file
    filepath = download_single_file(download_dir)
    
    if filepath:
        # Test NetCDF operations
        test_netcdf_operations(filepath)
        
        print_section("Summary")
        print(f"✓ File downloaded successfully: {filepath}")
        print(f"✓ File size: {os.path.getsize(filepath) / (1024**2):.2f} MB")
        
        # Try to open with EMITL2ARFL
        print("\nTesting with EMITL2ARFL package...")
        try:
            from EMITL2ARFL.validate_NetCDF_file import validate_NetCDF_file
            if validate_NetCDF_file(filepath):
                print("✓ EMITL2ARFL validation PASSED")
            else:
                print("✗ EMITL2ARFL validation FAILED")
        except Exception as e:
            print(f"✗ EMITL2ARFL validation error: {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()
    else:
        print_section("Summary")
        print("✗ Download failed")
    
    print("\n" + "="*80)
    print("  Diagnostic Complete")
    print("="*80 + "\n")

if __name__ == '__main__':
    main()
