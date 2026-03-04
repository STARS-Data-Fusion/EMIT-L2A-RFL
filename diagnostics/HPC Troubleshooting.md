# Troubleshooting NetCDF/HDF Errors on HPC Systems

## Common Issue: HDF Error -101

If you're seeing errors like this on HPC systems:

```
[Errno -101] NetCDF: HDF error: '/path/to/file.nc'
ValueError: Reflectance file is not a valid NetCDF file
```

This indicates file corruption, typically caused by:
- Network interruptions during download
- Filesystem buffering issues (NFS, Lustre, etc.)
- Incomplete file writes on network filesystems
- Storage system congestion

## Solutions

### 1. Increase Retry Parameters

The `retrieve_EMIT_L2A_RFL_granule()` function now supports configurable retry behavior:

```python
from EMITL2ARFL import retrieve_EMIT_L2A_RFL_granule

# Increase retries and delays for HPC environments
granule = retrieve_EMIT_L2A_RFL_granule(
    remote_granule=granule,
    download_directory="~/data/EMIT_download",
    max_retries=5,        # Try up to 5 times (default: 3)
    retry_delay=5.0       # Wait 5 seconds between attempts (default: 2.0)
)
```

The function uses exponential backoff, so with `retry_delay=5.0`:
- 1st retry: wait 5 seconds
- 2nd retry: wait 10 seconds
- 3rd retry: wait 20 seconds
- etc.

### 2. Diagnose Existing Files

Use the diagnostic tool to check which files are corrupted:

```bash
# Check a single file
python diagnose_netcdf.py ~/data/EMIT_download/EMIT_L2A_RFL_001_20220813T232430_2222515_007.nc

# Check all files in a directory
python diagnose_netcdf.py ~/data/EMIT_download/
```

Or from Python:

```python
from EMITL2ARFL import diagnose_netcdf_file, diagnose_directory

# Diagnose a single file
result = diagnose_netcdf_file("path/to/file.nc", verbose=True)

# Diagnose all files in a directory
results = diagnose_directory("path/to/directory", verbose=True)
```

### 3. Clean Up Corrupted Files

Delete corrupted files before re-downloading:

```bash
# Find and remove corrupted files
find ~/data/EMIT_download -name "*.nc" -type f -exec python -c "
from EMITL2ARFL import validate_NetCDF_file
from EMITL2ARFL.exceptions import NetCDFValidationError
import sys
try:
    validate_NetCDF_file(sys.argv[1])
except NetCDFValidationError:
    print('Corrupted:', sys.argv[1])
" {} \;
```

Or use the diagnostic script's output to identify files to delete.

### 4. HPC-Specific Best Practices

**For Slurm/PBS Job Scripts:**

```bash
#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --mem=32GB

# Load your environment
module load python/3.10
conda activate EMITL2ARFL

# Set longer retry parameters for HPC
export EMIT_MAX_RETRIES=5
export EMIT_RETRY_DELAY=10.0

# Run your processing script
python generate_kings_canyon_timeseries.py
```

**In Your Python Script:**

```python
import os

# Read environment variables if set
max_retries = int(os.environ.get('EMIT_MAX_RETRIES', 3))
retry_delay = float(os.environ.get('EMIT_RETRY_DELAY', 2.0))

# Use when retrieving granules
filenames = generate_EMIT_L2A_RFL_timeseries(
    start_date_UTC=start_date_UTC,
    end_date_UTC=end_date_UTC,
    geometry=grid,
    output_directory=output_directory,
    max_retries=max_retries,
    retry_delay=retry_delay
)
```

### 5. Filesystem-Specific Issues

**NFS/Lustre Filesystems:**
- Files may not immediately appear after download
- The new code includes `sync()` calls and stability checks
- Consider using local scratch space for downloads, then copying to network storage

**Example:**

```python
import shutil
from pathlib import Path

# Download to local scratch
local_scratch = Path("/tmp/emit_download")
local_scratch.mkdir(exist_ok=True)

granule = retrieve_EMIT_L2A_RFL_granule(
    remote_granule=granule,
    download_directory=str(local_scratch),
    max_retries=5
)

# Copy to permanent storage after successful download
permanent_storage = Path("~/data/EMIT_download").expanduser()
permanent_storage.mkdir(exist_ok=True)
shutil.copytree(local_scratch, permanent_storage, dirs_exist_ok=True)
```

## New Features for HPC Reliability

The following improvements have been added specifically for HPC environments:

1. **Exponential Backoff**: Automatically increases wait time between retries
2. **Filesystem Sync**: Forces filesystem synchronization after file operations
3. **File Stability Checks**: Waits for files to stop changing before validation
4. **Safe File Removal**: Robust file deletion with retries
5. **Detailed HDF Error Detection**: Specifically identifies HDF/NetCDF format errors
6. **Diagnostic Tools**: Built-in tools to identify and troubleshoot corrupted files

## Still Having Issues?

If you continue to experience problems:

1. **Check disk space**: Ensure sufficient space for large NetCDF files (1-2 GB each)
2. **Check filesystem health**: Run filesystem checks (`df -h`, `lfs df` for Lustre)
3. **Try different storage**: Some HPC filesystems are more reliable than others
4. **Download during off-peak hours**: Reduce network congestion
5. **Contact HPC support**: May indicate systemic storage issues

## Technical Details

### What Changed

**In `retrieve_EMIT_L2A_RFL_granule.py`:**
- Added `retry_delay` parameter with exponential backoff
- Added `sync()` calls to force filesystem writes
- Added file stability checks using `wait_for_file_stability()`
- Replaced basic file removal with `safe_file_remove()`

**In `validate_NetCDF_file.py`:**
- Enhanced error detection for HDF errors (errno -101)
- Added optional integrity checking
- Improved error messages with specific recommendations

**New modules:**
- `file_utils.py`: Utilities for robust file operations on network filesystems
- `diagnose_netcdf_issues.py`: Comprehensive diagnostic tools
- `diagnose_netcdf.py`: Standalone diagnostic script

### How It Works

1. **Initial validation**: Checks if existing files are valid
2. **Retry loop**: If validation fails:
   - Waits with exponential backoff
   - Safely removes corrupted files
   - Forces filesystem sync
   - Re-downloads files
   - Waits for files to stabilize
   - Re-validates
3. **Failure handling**: After max retries, reports which files remain corrupted

The improved retry logic significantly reduces corruption issues on HPC systems while maintaining backward compatibility.
