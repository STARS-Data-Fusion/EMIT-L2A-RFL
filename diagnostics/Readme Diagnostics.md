# Diagnostic Tools for HDF/NetCDF Issues

This directory contains several diagnostic scripts to help troubleshoot HDF5 error -101 and other NetCDF-related issues on HPC systems.

## 🆕 Recommended: All-in-One Diagnostic Tool

### `debug_hpc_environment.py` (New!)
**Purpose:** Comprehensive diagnostic with automatic file download - best for first-time setup and troubleshooting

**What it does:**
- **Automatically downloads** a sample EMIT L2A RFL file if you don't have one
- Tests environment configuration (libraries, variables, filesystem)
- Tests file access with multiple methods
- Identifies HPC-specific issues (file locking, network filesystems)
- Provides specific, actionable recommendations
- Works on both local machines and HPC systems

**Usage:**
```bash
# Download sample file and run full diagnostics (recommended)
conda activate EMITL2ARFL
python debug_hpc_environment.py

# Or test your own file
python debug_hpc_environment.py ~/data/your_file.nc
```

**When to use:**
- First-time setup on HPC system
- After installing/updating EMITL2ARFL
- When experiencing "HDF Error -101"
- Before running large processing jobs
- Comparing local vs HPC environments

---

## Other Available Diagnostic Tools

### 1. `diagnose_netcdf_environment.py`
**Purpose:** Test HDF5/NetCDF library installation and configuration

**What it tests:**
- Library versions (netCDF4, HDF5, h5py, NumPy)
- Basic NetCDF file operations (create, read, write)
- Different NetCDF format support (NETCDF4, NETCDF3, etc.)
- Compression capabilities
- HDF5 file locking behavior
- Large file support (>2GB)
- Error handling for corrupt files
- EMIT-like file structure creation/reading
- HDF5 environment variables

**Usage:**
```bash
# On Mac/local
conda activate EMITL2ARFL
python diagnose_netcdf_environment.py

# On HPC
python diagnose_netcdf_environment.py
```

**Expected output:** All tests should pass if libraries are correctly installed.

### 2. `diagnose_single_download.py`
**Purpose:** Comprehensive diagnostic for downloading and validating a single EMIT granule

**What it tests:**
- NASA Earthdata authentication
- Granule search functionality
- File download to specified directory
- NetCDF file validation
- Environment information (Python, packages, system)
- File integrity checks

**Usage:**
```bash
python diagnose_single_download.py
```

**Options:**
- Prompts for download directory (default: `/tmp/emit_test`)
- Tests a single granule from 2024-08-15
- Provides detailed logging of each step

### 3. `test_single_granule.py`
**Purpose:** Simple test script for downloading and validating a single granule

**What it tests:**
- Basic download functionality
- File validation
- Cache clearing on corruption

**Usage:**
```bash
python test_single_granule.py
```

### 4. `clean_cache.py`
**Purpose:** Manual cache cleanup utility

**What it does:**
- Removes all cached EMIT granule files
- Validates each file before removal (optional)
- Reports corrupted files found

**Usage:**
```bash
python clean_cache.py
```

## Running Diagnostics on HPC

### Quick Test Sequence

1. **Test HDF5/NetCDF libraries:**
   ```bash
   python diagnose_netcdf_environment.py > hpc_netcdf_test.log 2>&1
   ```
   
   Look for:
   - ✓ ALL TESTS PASSED
   - Any failures in error handling or file operations

2. **Test single granule download:**
   ```bash
   python diagnose_single_download.py
   ```
   
   When prompted, use `/tmp` for the download directory (or another local scratch space).
   
   Look for:
   - Authentication success
   - Download completion
   - Validation passing
   - No HDF error -101

3. **Compare results:**
   - If Mac tests pass but HPC tests fail → HPC filesystem or threading issue
   - If both fail → library or authentication issue
   - If validation fails on HPC but not Mac → likely multi-threaded download issue

## Common Issues and Solutions

### Issue: HDF error -101 on HPC

**Symptoms:**
- Files download but fail validation with "NetCDF: HDF error"
- Error occurs consistently on first download attempt
- File sizes are correct but HDF format is corrupt

**Cause:**
- Multi-threaded downloads (8+ threads) on network filesystems
- Network filesystem buffering can't handle concurrent writes

**Solution:**
- Default to single-threaded downloads: `threads=1` (already set as default)
- Download to local scratch space (not network filesystem): `/tmp` or `/scratch`

**Test:**
```python
from EMITL2ARFL import retrieve_EMIT_L2A_RFL_granule
import earthaccess

earthaccess.login()
results = earthaccess.search_data(
    short_name='EMITL2ARFL',
    temporal=('2024-08-15', '2024-08-15')
)

# Single-threaded (default, safe for HPC)
granule = retrieve_EMIT_L2A_RFL_granule(
    remote_granule=results[0],
    download_directory='/tmp/emit_test'
)

# Or explicitly set threads=1
granule = retrieve_EMIT_L2A_RFL_granule(
    remote_granule=results[0],
    download_directory='/tmp/emit_test',
    threads=1  # Single-threaded for HPC safety
)
```

### Issue: Persistent corruption even after retries

**Symptoms:**
- All retry attempts fail with same error
- Validation catches corruption immediately
- Cache clearing doesn't help

**Cause:**
- Using multi-threaded downloads on network filesystem
- Corrupted files being reused from cache

**Solution:**
1. Ensure `threads=1` (now default)
2. Clear cache manually if needed:
   ```bash
   python clean_cache.py
   ```
3. Download to local storage (not network filesystem)

### Issue: Authentication failures

**Symptoms:**
- Cannot login to NASA Earthdata
- No search results returned

**Solution:**
```bash
# Create/update .netrc file
echo "machine urs.earthdata.nasa.gov login YOUR_USERNAME password YOUR_PASSWORD" > ~/.netrc
chmod 600 ~/.netrc
```

## Thread Configuration

The package now defaults to **single-threaded downloads** (`threads=1`) for maximum reliability on HPC/network filesystems.

**When to use single-threaded (threads=1):**
- ✓ HPC systems with network filesystems (NFS, Lustre, GPFS)
- ✓ When experiencing HDF error -101
- ✓ When reliability is more important than speed
- ✓ Default behavior (safest option)

**When to use multi-threaded (threads=8):**
- ✓ Local systems with local storage (SSD, HDD)
- ✓ When download speed is critical
- ✓ After confirming no corruption issues

**Example:**
```python
# HPC (default) - single-threaded for reliability
granule = retrieve_EMIT_L2A_RFL_granule(remote_granule)

# Local system - multi-threaded for speed
granule = retrieve_EMIT_L2A_RFL_granule(remote_granule, threads=8)
```

## Expected Test Results

### Mac/Local System (Working)
```
✓ ALL TESTS PASSED
Your NetCDF/HDF5 environment is properly configured.
```

### HPC System (Should Pass After Fix)
```
✓ ALL TESTS PASSED
Your NetCDF/HDF5 environment is properly configured.

# With single-threaded downloads:
[INFO] Downloading granule files (attempt 1/3, threads=1)...
[INFO] Downloaded file validated successfully
```

## Getting Help

If diagnostics reveal persistent issues:

1. **Share diagnostic output:**
   - Output from `diagnose_netcdf_environment.py`
   - Logs from failed downloads
   - HPC system details (filesystem type, scratch space location)

2. **Check HPC documentation:**
   - Recommended scratch space locations
   - Known filesystem limitations
   - Python/conda module recommendations

3. **Contact HPC support:**
   - Share the diagnostic results
   - Ask about HDF5 file locking settings
   - Request guidance on optimal download locations
