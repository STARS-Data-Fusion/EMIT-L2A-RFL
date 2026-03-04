# HPC NetCDF Error Fixes - Summary

## Problem
Users were experiencing HDF error -101 when opening NetCDF files on HPC machines, indicating file corruption during download due to:
- Network interruptions
- Filesystem buffering issues (NFS, Lustre, etc.)
- Incomplete file writes
- Storage system congestion

## Solution Overview
Added robust retry logic, filesystem synchronization, file stability checks, and diagnostic tools specifically designed for HPC/network filesystem environments.

## Changes Made

### 1. Enhanced `retrieve_EMIT_L2A_RFL_granule.py`
- **Added `retry_delay` parameter** (default: 2.0 seconds)
  - Allows configurable wait time between retry attempts
  - Uses exponential backoff (2x increase each retry)
- **Added filesystem synchronization**
  - Calls `sync()` before and after file operations
  - Ensures writes complete on network filesystems
- **Added file stability checks**
  - Waits for files to stop changing after download
  - Prevents premature validation of incomplete files
- **Improved file removal**
  - Uses `safe_file_remove()` with retry logic
  - Handles network filesystem delays

### 2. Enhanced `validate_NetCDF_file.py`
- **Improved HDF error detection**
  - Specifically identifies errno -101 HDF errors
  - Provides actionable error messages
- **Added optional integrity checking**
  - Can verify file data is accessible
  - Helps identify subtle corruption issues
- **Better error messages**
  - Includes file sizes and specific error codes
  - Provides recommendations for fixing issues

### 3. New `file_utils.py` Module
Utilities for robust file operations on HPC systems:
- `compute_file_checksum()` - Calculate file checksums
- `wait_for_file_stability()` - Wait for file writes to complete
- `safe_file_remove()` - Safely remove files with retries
- `verify_file_readable()` - Check basic file readability

### 4. New `diagnose_netcdf_issues.py` Module
Comprehensive diagnostic tools:
- `diagnose_netcdf_file()` - Detailed diagnostics for single files
- `diagnose_directory()` - Batch diagnostics for directories
- Provides actionable recommendations for fixing issues
- Can be used as standalone diagnostic tool

### 5. New `diagnose_netcdf.py` Script
Standalone command-line tool for diagnosing NetCDF issues:
```bash
python diagnose_netcdf.py <file_or_directory>
```

### 6. New `HPC_TROUBLESHOOTING.md` Documentation
Comprehensive guide including:
- Problem explanation
- Configuration examples
- Best practices for HPC environments
- Troubleshooting steps
- Technical details

## Usage Examples

### Basic Usage (backward compatible)
```python
from EMITL2ARFL import retrieve_EMIT_L2A_RFL_granule

# Works exactly as before with improved reliability
granule = retrieve_EMIT_L2A_RFL_granule(
    remote_granule=granule,
    download_directory="~/data/EMIT_download"
)
```

### HPC-Optimized Usage
```python
# Increase retries and delays for unreliable networks
granule = retrieve_EMIT_L2A_RFL_granule(
    remote_granule=granule,
    download_directory="~/data/EMIT_download",
    max_retries=5,        # Try up to 5 times (default: 3)
    retry_delay=5.0       # Wait 5 seconds, then 10, 20, 40... (default: 2.0)
)
```

### Diagnose Existing Files
```python
from EMITL2ARFL import diagnose_netcdf_file, diagnose_directory

# Check a single file
result = diagnose_netcdf_file("path/to/file.nc", verbose=True)
print(f"Valid: {result['valid_netcdf']}")

# Check all files in directory
results = diagnose_directory("path/to/directory", verbose=True)
```

### Command-Line Diagnostics
```bash
# Diagnose a single file
python diagnose_netcdf.py ~/data/EMIT_download/EMIT_L2A_RFL_001_20220813T232430_2222515_007.nc

# Diagnose all files in directory
python diagnose_netcdf.py ~/data/EMIT_download/
```

## Benefits

1. **Automatic Recovery**: Files corrupted during download are automatically detected and re-downloaded
2. **Exponential Backoff**: Progressive delays give network/filesystem time to stabilize
3. **Filesystem Sync**: Ensures writes complete on network filesystems
4. **Stability Checks**: Prevents premature validation of incomplete files
5. **Backward Compatible**: Existing code works without changes, with improved reliability
6. **Diagnostic Tools**: Easy identification of corrupted files
7. **Clear Documentation**: HPC-specific guidance and troubleshooting

## Testing

All existing tests pass:
```bash
make test  # 10/10 tests pass
```

New functionality verified:
- ✅ Module imports work correctly
- ✅ New utilities are accessible
- ✅ Diagnostic script runs
- ✅ Backward compatibility maintained

## Files Modified

1. `EMITL2ARFL/retrieve_EMIT_L2A_RFL_granule.py` - Enhanced retry logic
2. `EMITL2ARFL/validate_NetCDF_file.py` - Improved error detection
3. `EMITL2ARFL/EMITL2ARFL.py` - Added new module exports

## Files Created

1. `EMITL2ARFL/file_utils.py` - HPC file operation utilities
2. `EMITL2ARFL/diagnose_netcdf_issues.py` - Diagnostic tools
3. `diagnose_netcdf.py` - Standalone diagnostic script
4. `HPC_TROUBLESHOOTING.md` - Comprehensive documentation
5. `CHANGES_SUMMARY.md` - This file

## Next Steps for Users

1. **Re-install the package** (if using HPC):
   ```bash
   conda activate EMITL2ARFL
   pip install -e .
   ```

2. **Diagnose existing corrupted files**:
   ```bash
   python diagnose_netcdf.py ~/data/EMIT_download/
   ```

3. **Delete corrupted files** (based on diagnostic output)

4. **Re-run downloads** with enhanced retry parameters:
   ```python
   generate_EMIT_L2A_RFL_timeseries(
       ...,
       max_retries=5,
       retry_delay=5.0
   )
   ```

5. **Read HPC_TROUBLESHOOTING.md** for detailed guidance

## Technical Notes

- Uses Python's `os.sync()` to force filesystem writes
- Implements exponential backoff: `delay * (2 ** (retry - 1))`
- File stability checked by monitoring size changes over time
- Compatible with NFS, Lustre, and other network filesystems
- All changes are backward compatible (new parameters are optional)
