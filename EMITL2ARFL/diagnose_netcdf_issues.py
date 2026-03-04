"""
Diagnostic tool for troubleshooting NetCDF file issues on HPC systems.

This script helps identify and diagnose NetCDF file corruption issues,
particularly those related to HDF errors on HPC/network filesystems.
"""

import logging
import sys
from pathlib import Path
from typing import List, Union
import netCDF4

from .validate_NetCDF_file import validate_NetCDF_file
from .file_utils import compute_file_checksum, verify_file_readable
from .exceptions import NetCDFValidationError

logger = logging.getLogger(__name__)


def diagnose_netcdf_file(filepath: Union[str, Path], verbose: bool = True) -> dict:
    """
    Perform comprehensive diagnostics on a NetCDF file.
    
    Parameters
    ----------
    filepath : str or Path
        Path to the NetCDF file to diagnose
    verbose : bool, optional
        If True, print detailed diagnostic information. Defaults to True.
    
    Returns
    -------
    dict
        Dictionary containing diagnostic results with keys:
        - 'exists': bool - File exists
        - 'size': int - File size in bytes
        - 'readable': bool - File can be opened for reading
        - 'valid_netcdf': bool - File passes NetCDF validation
        - 'dimensions': dict - NetCDF dimensions (if readable)
        - 'variables': list - NetCDF variable names (if readable)
        - 'error': str - Error message if validation failed
        - 'recommendation': str - Recommended action
    """
    filepath = Path(filepath)
    result = {
        'filepath': str(filepath),
        'exists': False,
        'size': 0,
        'readable': False,
        'valid_netcdf': False,
        'dimensions': {},
        'variables': [],
        'error': None,
        'recommendation': None
    }
    
    # Check existence
    if not filepath.exists():
        result['error'] = "File does not exist"
        result['recommendation'] = "Check file path or download the file"
        if verbose:
            print(f"❌ File does not exist: {filepath}")
        return result
    
    result['exists'] = True
    
    # Check size
    try:
        result['size'] = filepath.stat().st_size
        if verbose:
            size_mb = result['size'] / (1024 * 1024)
            print(f"📁 File size: {size_mb:.2f} MB ({result['size']} bytes)")
    except Exception as e:
        result['error'] = f"Cannot access file stats: {e}"
        if verbose:
            print(f"❌ {result['error']}")
        return result
    
    # Check if empty
    if result['size'] == 0:
        result['error'] = "File is empty (0 bytes)"
        result['recommendation'] = "Delete and re-download the file"
        if verbose:
            print(f"❌ File is empty")
        return result
    
    # Check basic readability
    result['readable'] = verify_file_readable(filepath)
    if not result['readable']:
        result['error'] = "File cannot be opened for reading"
        result['recommendation'] = "Check file permissions or filesystem issues"
        if verbose:
            print(f"❌ File cannot be read (permission or I/O error)")
        return result
    
    if verbose:
        print(f"✅ File is readable")
    
    # Try NetCDF validation
    try:
        validate_NetCDF_file(filepath, file_type="NetCDF")
        result['valid_netcdf'] = True
        if verbose:
            print(f"✅ File passes NetCDF validation")
    except NetCDFValidationError as e:
        result['error'] = str(e)
        
        # Check for specific error patterns
        error_str = str(e).lower()
        if "errno -101" in error_str or "hdf error" in error_str:
            result['recommendation'] = (
                "HDF error detected - file is likely corrupted during download. "
                "Delete the file and re-download. "
                "On HPC systems, consider increasing retry delays."
            )
            if verbose:
                print(f"❌ HDF/NetCDF format error (likely corruption)")
                print(f"   Error: {e}")
        else:
            result['recommendation'] = "File validation failed - see error message"
            if verbose:
                print(f"❌ NetCDF validation failed")
                print(f"   Error: {e}")
        
        return result
    
    # Try to read NetCDF structure
    try:
        with netCDF4.Dataset(filepath, 'r') as ds:
            result['dimensions'] = {name: dim.size for name, dim in ds.dimensions.items()}
            result['variables'] = list(ds.variables.keys())
            
            if verbose:
                print(f"✅ NetCDF structure:")
                print(f"   Dimensions: {list(result['dimensions'].keys())}")
                print(f"   Variables: {len(result['variables'])} found")
                
                # Show a few sample variables
                sample_vars = result['variables'][:5]
                for var in sample_vars:
                    var_obj = ds.variables[var]
                    print(f"      - {var}: shape={var_obj.shape}, dtype={var_obj.dtype}")
                
                if len(result['variables']) > 5:
                    print(f"      ... and {len(result['variables']) - 5} more variables")
    
    except Exception as e:
        result['error'] = f"Cannot read NetCDF structure: {e}"
        result['recommendation'] = "File may be corrupted - delete and re-download"
        if verbose:
            print(f"❌ Cannot read NetCDF structure")
            print(f"   Error: {e}")
        return result
    
    if verbose:
        print(f"\n✅ File appears to be valid and complete")
    
    return result


def diagnose_directory(directory: Union[str, Path], pattern: str = "*.nc", verbose: bool = True) -> List[dict]:
    """
    Diagnose all NetCDF files in a directory.
    
    Parameters
    ----------
    directory : str or Path
        Path to directory containing NetCDF files
    pattern : str, optional
        File pattern to match. Defaults to "*.nc".
    verbose : bool, optional
        If True, print detailed diagnostic information. Defaults to True.
    
    Returns
    -------
    List[dict]
        List of diagnostic results for each file
    """
    directory = Path(directory)
    
    if not directory.exists():
        logger.error(f"Directory does not exist: {directory}")
        return []
    
    files = list(directory.glob(pattern))
    
    if not files:
        logger.warning(f"No files matching '{pattern}' found in {directory}")
        return []
    
    if verbose:
        print(f"\n{'='*70}")
        print(f"Diagnosing {len(files)} NetCDF files in: {directory}")
        print(f"{'='*70}\n")
    
    results = []
    for i, filepath in enumerate(files, 1):
        if verbose:
            print(f"\n[{i}/{len(files)}] Diagnosing: {filepath.name}")
            print("-" * 70)
        
        result = diagnose_netcdf_file(filepath, verbose=verbose)
        results.append(result)
    
    # Summary
    if verbose:
        print(f"\n{'='*70}")
        print("SUMMARY")
        print(f"{'='*70}")
        
        valid_count = sum(1 for r in results if r['valid_netcdf'])
        invalid_count = len(results) - valid_count
        
        print(f"Total files: {len(results)}")
        print(f"Valid: {valid_count}")
        print(f"Invalid: {invalid_count}")
        
        if invalid_count > 0:
            print(f"\nInvalid files:")
            for r in results:
                if not r['valid_netcdf']:
                    print(f"  - {Path(r['filepath']).name}")
                    if r['error']:
                        print(f"    Error: {r['error'][:100]}")
                    if r['recommendation']:
                        print(f"    Action: {r['recommendation']}")
    
    return results


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Diagnose NetCDF file issues on HPC systems"
    )
    parser.add_argument(
        "path",
        help="Path to NetCDF file or directory containing NetCDF files"
    )
    parser.add_argument(
        "--pattern",
        default="*.nc",
        help="File pattern to match (default: *.nc)"
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress detailed output"
    )
    
    args = parser.parse_args()
    
    path = Path(args.path)
    
    if path.is_file():
        result = diagnose_netcdf_file(path, verbose=not args.quiet)
        sys.exit(0 if result['valid_netcdf'] else 1)
    elif path.is_dir():
        results = diagnose_directory(path, pattern=args.pattern, verbose=not args.quiet)
        invalid_count = sum(1 for r in results if not r['valid_netcdf'])
        sys.exit(0 if invalid_count == 0 else 1)
    else:
        print(f"Error: Path does not exist: {path}", file=sys.stderr)
        sys.exit(1)
