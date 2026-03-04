import sys
from os.path import expanduser

from EMITL2ARFL import validate_NetCDF_file

if len(sys.argv) < 2:
    print("Usage: python test_file_validation.py <netcdf_file>")
    sys.exit(1)

filename = sys.argv[1]

try:
    validate_NetCDF_file(filename)
    print(f"✓ Valid NetCDF file: {filename}")
except Exception as e:
    print(f"✗ Invalid NetCDF file: {filename} - {type(e).__name__}: {e}")

