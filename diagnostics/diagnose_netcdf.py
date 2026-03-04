#!/usr/bin/env python3
"""
Standalone diagnostic tool for NetCDF file issues on HPC systems.

Usage:
    python diagnose_netcdf.py <file_or_directory>
    python diagnose_netcdf.py ~/data/EMIT_L2A_RFL/EMIT_L2A_RFL_001_20220813T232430_2222515_007/
"""

from EMITL2ARFL.diagnose_netcdf_issues import diagnose_netcdf_file, diagnose_directory
import sys
from pathlib import Path


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ['-h', '--help']:
        print("Usage: python diagnose_netcdf.py <file_or_directory>")
        print("\nExamples:")
        print("  python diagnose_netcdf.py data.nc")
        print("  python diagnose_netcdf.py ~/data/EMIT_download/")
        sys.exit(0 if len(sys.argv) > 1 else 1)
    
    path = Path(sys.argv[1]).expanduser()
    
    if not path.exists():
        print(f"Error: Path does not exist: {path}")
        sys.exit(1)
    
    if path.is_file():
        result = diagnose_netcdf_file(path, verbose=True)
        sys.exit(0 if result['valid_netcdf'] else 1)
    elif path.is_dir():
        results = diagnose_directory(path, pattern="*.nc", verbose=True)
        invalid_count = sum(1 for r in results if not r['valid_netcdf'])
        sys.exit(0 if invalid_count == 0 else 1)


if __name__ == "__main__":
    main()
