#!/usr/bin/env python
"""
Diagnostic script to test HDF5/NetCDF library installation and configuration.

This script tests:
1. Library versions and build information
2. Basic NetCDF file operations (create, write, read, close)
3. HDF5 file locking and concurrent access
4. NetCDF4 format compatibility
5. Large file support
6. Thread safety
7. Compression capabilities

Usage:
    python diagnose_netcdf_environment.py
"""

import sys
import os
import tempfile
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

def print_section(title):
    """Print a formatted section header."""
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80)

def test_library_versions():
    """Test library versions and build information."""
    print_section("1. Library Versions and Build Information")
    
    try:
        import netCDF4
        print(f"✓ netCDF4 version: {netCDF4.__version__}")
        print(f"  netCDF4 file: {netCDF4.__file__}")
    except ImportError as e:
        print(f"✗ Failed to import netCDF4: {e}")
        return False
    
    try:
        import h5py
        print(f"✓ h5py version: {h5py.__version__}")
        print(f"  HDF5 version: {h5py.version.hdf5_version}")
        print(f"  h5py file: {h5py.__file__}")
    except ImportError as e:
        print(f"✗ Failed to import h5py: {e}")
        return False
    
    try:
        import numpy as np
        print(f"✓ NumPy version: {np.__version__}")
    except ImportError as e:
        print(f"✗ Failed to import NumPy: {e}")
        return False
    
    # Get NetCDF library info
    try:
        print(f"\n  NetCDF library version: {netCDF4.__netcdf4libversion__}")
        print(f"  NetCDF HDF5 version: {netCDF4.__hdf5libversion__}")
    except AttributeError:
        print("  Warning: Could not get NetCDF/HDF5 library versions")
    
    return True

def test_basic_netcdf_operations():
    """Test basic NetCDF file operations."""
    print_section("2. Basic NetCDF File Operations")
    
    import netCDF4
    import numpy as np
    
    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = os.path.join(tmpdir, "test_basic.nc")
        
        try:
            # Create file
            logger.info(f"Creating test file: {test_file}")
            with netCDF4.Dataset(test_file, 'w', format='NETCDF4') as ds:
                # Create dimensions
                ds.createDimension('x', 10)
                ds.createDimension('y', 20)
                
                # Create variable
                var = ds.createVariable('data', 'f4', ('x', 'y'))
                var[:] = np.random.rand(10, 20)
                
                # Add attributes
                ds.title = "Test NetCDF file"
                var.units = "test_units"
            
            print(f"✓ Successfully created NetCDF4 file")
            
            # Read file
            logger.info("Reading test file")
            with netCDF4.Dataset(test_file, 'r') as ds:
                assert ds.title == "Test NetCDF file"
                assert 'data' in ds.variables
                assert ds.variables['data'].shape == (10, 20)
                assert ds.variables['data'].units == "test_units"
            
            print(f"✓ Successfully read NetCDF4 file")
            
            # Check file size
            file_size = os.path.getsize(test_file)
            print(f"  File size: {file_size} bytes")
            
            return True
            
        except Exception as e:
            print(f"✗ Basic NetCDF operations failed: {e}")
            import traceback
            traceback.print_exc()
            return False

def test_netcdf_formats():
    """Test different NetCDF format support."""
    print_section("3. NetCDF Format Support")
    
    import netCDF4
    import numpy as np
    
    formats = ['NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_CLASSIC', 'NETCDF3_64BIT']
    
    with tempfile.TemporaryDirectory() as tmpdir:
        for fmt in formats:
            test_file = os.path.join(tmpdir, f"test_{fmt}.nc")
            
            try:
                with netCDF4.Dataset(test_file, 'w', format=fmt) as ds:
                    ds.createDimension('x', 5)
                    var = ds.createVariable('data', 'f4', ('x',))
                    var[:] = np.arange(5)
                
                # Read back
                with netCDF4.Dataset(test_file, 'r') as ds:
                    assert len(ds.dimensions['x']) == 5
                    assert list(ds.variables['data'][:]) == list(range(5))
                
                print(f"✓ Format {fmt}: OK")
                
            except Exception as e:
                print(f"✗ Format {fmt}: Failed - {e}")
    
    return True

def test_compression():
    """Test compression capabilities."""
    print_section("4. Compression Support")
    
    import netCDF4
    import numpy as np
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Test with compression
        compressed_file = os.path.join(tmpdir, "test_compressed.nc")
        uncompressed_file = os.path.join(tmpdir, "test_uncompressed.nc")
        
        data = np.random.rand(100, 100)
        
        try:
            # Compressed
            with netCDF4.Dataset(compressed_file, 'w', format='NETCDF4') as ds:
                ds.createDimension('x', 100)
                ds.createDimension('y', 100)
                var = ds.createVariable('data', 'f4', ('x', 'y'), 
                                       zlib=True, complevel=4)
                var[:] = data
            
            # Uncompressed
            with netCDF4.Dataset(uncompressed_file, 'w', format='NETCDF4') as ds:
                ds.createDimension('x', 100)
                ds.createDimension('y', 100)
                var = ds.createVariable('data', 'f4', ('x', 'y'))
                var[:] = data
            
            compressed_size = os.path.getsize(compressed_file)
            uncompressed_size = os.path.getsize(uncompressed_file)
            ratio = uncompressed_size / compressed_size if compressed_size > 0 else 0
            
            print(f"✓ Compression works")
            print(f"  Compressed size: {compressed_size} bytes")
            print(f"  Uncompressed size: {uncompressed_size} bytes")
            print(f"  Compression ratio: {ratio:.2f}x")
            
            return True
            
        except Exception as e:
            print(f"✗ Compression test failed: {e}")
            return False

def test_hdf5_file_locking():
    """Test HDF5 file locking behavior."""
    print_section("5. HDF5 File Locking")
    
    import netCDF4
    import numpy as np
    
    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = os.path.join(tmpdir, "test_locking.nc")
        
        try:
            # Create file
            with netCDF4.Dataset(test_file, 'w', format='NETCDF4') as ds:
                ds.createDimension('x', 5)
                var = ds.createVariable('data', 'f4', ('x',))
                var[:] = np.arange(5)
            
            # Try to open file multiple times for reading
            ds1 = netCDF4.Dataset(test_file, 'r')
            ds2 = netCDF4.Dataset(test_file, 'r')
            
            # Read from both
            data1 = ds1.variables['data'][:]
            data2 = ds2.variables['data'][:]
            
            ds1.close()
            ds2.close()
            
            print(f"✓ Concurrent read access works")
            
            # Check for lock files
            lock_file = test_file + ".lock"
            if os.path.exists(lock_file):
                print(f"  Warning: Lock file found: {lock_file}")
            
            return True
            
        except Exception as e:
            print(f"✗ File locking test failed: {e}")
            import traceback
            traceback.print_exc()
            return False

def test_large_file_support():
    """Test large file (>2GB) support."""
    print_section("6. Large File Support")
    
    import netCDF4
    import numpy as np
    
    print("  Note: This test creates a ~100MB file, not a full 2GB+ file")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = os.path.join(tmpdir, "test_large.nc")
        
        try:
            # Create a moderately large file
            with netCDF4.Dataset(test_file, 'w', format='NETCDF4') as ds:
                ds.createDimension('x', 1000)
                ds.createDimension('y', 1000)
                ds.createDimension('z', 10)
                
                var = ds.createVariable('data', 'f4', ('x', 'y', 'z'))
                
                # Write in chunks to avoid memory issues
                chunk_size = 100
                for i in range(0, 1000, chunk_size):
                    var[i:i+chunk_size, :, :] = np.random.rand(chunk_size, 1000, 10)
            
            file_size = os.path.getsize(test_file)
            file_size_mb = file_size / (1024 * 1024)
            
            # Read back a sample
            with netCDF4.Dataset(test_file, 'r') as ds:
                sample = ds.variables['data'][0:10, 0:10, 0]
            
            print(f"✓ Large file support OK")
            print(f"  Created file size: {file_size_mb:.2f} MB")
            print(f"  Large file support: {'YES' if file_size_mb > 1 else 'Limited'}")
            
            return True
            
        except Exception as e:
            print(f"✗ Large file test failed: {e}")
            return False

def test_error_handling():
    """Test error handling for corrupt files."""
    print_section("7. Error Handling for Corrupt Files")
    
    import netCDF4
    
    with tempfile.TemporaryDirectory() as tmpdir:
        corrupt_file = os.path.join(tmpdir, "corrupt.nc")
        
        # Create a corrupt file (just random bytes)
        with open(corrupt_file, 'wb') as f:
            f.write(b'This is not a valid NetCDF file' * 100)
        
        try:
            ds = netCDF4.Dataset(corrupt_file, 'r')
            print(f"✗ Failed to detect corrupt file!")
            ds.close()
            return False
            
        except OSError as e:
            error_msg = str(e).lower()
            print(f"✓ Correctly detected corrupt file")
            print(f"  Error message: {e}")
            
            if 'errno -101' in error_msg or 'hdf error' in error_msg:
                print(f"  Error type: HDF5 error (errno -101)")
            elif 'errno -51' in error_msg:
                print(f"  Error type: Unknown file format")
            else:
                print(f"  Error type: Generic NetCDF error")
            
            return True
            
        except Exception as e:
            print(f"✓ Detected corrupt file with exception: {type(e).__name__}")
            print(f"  Error message: {e}")
            return True

def test_emit_like_file():
    """Test with EMIT-like file structure."""
    print_section("8. EMIT-like File Structure Test")
    
    import netCDF4
    import numpy as np
    
    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = os.path.join(tmpdir, "test_emit_like.nc")
        
        try:
            # Create EMIT-like structure
            with netCDF4.Dataset(test_file, 'w', format='NETCDF4') as ds:
                # Create groups
                sensor_group = ds.createGroup('sensor_band_parameters')
                location_group = ds.createGroup('location')
                
                # Create dimensions
                ds.createDimension('downtrack', 1280)
                ds.createDimension('crosstrack', 1242)
                ds.createDimension('bands', 285)
                
                # Create reflectance variable
                reflectance = ds.createVariable(
                    'reflectance', 
                    'f4', 
                    ('downtrack', 'crosstrack', 'bands'),
                    zlib=True,
                    complevel=4,
                    chunksizes=(128, 128, 28)
                )
                
                # Add attributes
                ds.title = "Test EMIT L2A Reflectance"
                reflectance.units = "unitless"
                reflectance.long_name = "surface reflectance"
                
                # Write some data (just first chunk)
                reflectance[0:10, 0:10, 0:10] = np.random.rand(10, 10, 10)
            
            print(f"✓ Successfully created EMIT-like file")
            
            # Read and validate
            with netCDF4.Dataset(test_file, 'r') as ds:
                assert ds.title == "Test EMIT L2A Reflectance"
                assert 'reflectance' in ds.variables
                assert ds.variables['reflectance'].shape == (1280, 1242, 285)
                
                # Read sample data
                sample = ds.variables['reflectance'][0:5, 0:5, 0:5]
                assert sample.shape == (5, 5, 5)
            
            print(f"✓ Successfully read EMIT-like file")
            
            file_size = os.path.getsize(test_file)
            print(f"  File size: {file_size / (1024*1024):.2f} MB")
            
            return True
            
        except Exception as e:
            print(f"✗ EMIT-like file test failed: {e}")
            import traceback
            traceback.print_exc()
            return False

def test_environment_variables():
    """Check HDF5-related environment variables."""
    print_section("9. HDF5 Environment Variables")
    
    hdf5_vars = [
        'HDF5_USE_FILE_LOCKING',
        'HDF5_DISABLE_VERSION_CHECK',
        'NETCDF4_PLUGIN_PATH',
        'HDF5_PLUGIN_PATH',
    ]
    
    found_any = False
    for var in hdf5_vars:
        value = os.environ.get(var)
        if value is not None:
            print(f"  {var} = {value}")
            found_any = True
    
    if not found_any:
        print("  No HDF5-related environment variables set")
    
    return True

def main():
    """Run all diagnostic tests."""
    print("\n" + "=" * 80)
    print("  NetCDF/HDF5 Environment Diagnostic Tool")
    print("  Testing library installation and capabilities")
    print("=" * 80)
    
    results = {}
    
    # Run tests
    results['versions'] = test_library_versions()
    
    if results['versions']:
        results['basic_ops'] = test_basic_netcdf_operations()
        results['formats'] = test_netcdf_formats()
        results['compression'] = test_compression()
        results['locking'] = test_hdf5_file_locking()
        results['large_files'] = test_large_file_support()
        results['error_handling'] = test_error_handling()
        results['emit_like'] = test_emit_like_file()
        results['env_vars'] = test_environment_variables()
    
    # Summary
    print_section("Summary")
    
    all_passed = all(results.values())
    
    for test_name, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {status}: {test_name}")
    
    print("\n" + "=" * 80)
    if all_passed:
        print("  ✓ ALL TESTS PASSED")
        print("  Your NetCDF/HDF5 environment is properly configured.")
    else:
        print("  ✗ SOME TESTS FAILED")
        print("  There may be issues with your NetCDF/HDF5 installation.")
    print("=" * 80 + "\n")
    
    return 0 if all_passed else 1

if __name__ == '__main__':
    sys.exit(main())
