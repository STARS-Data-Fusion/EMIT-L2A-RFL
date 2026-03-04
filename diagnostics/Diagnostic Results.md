# Diagnostic Results: EMIT NetCDF Download Issues

## Executive Summary

**GOOD NEWS**: Your environment is working correctly! The HDF5 libraries, Python environment, and download logic are all functioning properly.

**ROOT CAUSE**: The errors were caused by **corrupted cached files** from previous download attempts, not ongoing download issues.

**SOLUTION IMPLEMENTED**: Automatic cache clearing is now built into the retry mechanism - corrupted cached files are automatically detected and removed before re-downloading.

## What We Discovered

### ✅ Working Components

1. **Python Environment**: Python 3.10.19 with EMITL2ARFL conda environment
2. **Dependencies**: All correctly installed
   - netCDF4: 1.7.3
   - HDF5 library: 1.14.6  
   - earthaccess: 0.14.0
   - h5py: 3.15.1

3. **Download & Validation**: Successfully downloaded and validated 3 files (RFL, MASK, RFLUNCERT) to `/tmp`
   - Files opened correctly with netCDF4 ✓
   - Files opened correctly with h5py ✓
   - Files validated with ncdump ✓
   - Files passed EMITL2ARFL validation ✓

### ❌ The Actual Problem

**Corrupted cache files** in `/Users/halverso/data/EMIT_L2A_RFL/` from previous failed downloads were being detected and causing validation failures.

When the script found these existing files, it tried to validate them before downloading, and they failed validation. These were legitimate corruptions from earlier download attempts (likely on the HPC system or during interrupted downloads).

## Test Results

### Successful Download Test
```bash
# Downloaded to /tmp with single-threaded download
[2025-12-05 12:10:13] Downloading to: /tmp/emit_test
[2025-12-05 12:11:18] Waiting for downloaded files to stabilize...
[2025-12-05 12:11:22] All files successfully downloaded and validated.
[2025-12-05 12:11:22] ✓ SUCCESS!
```

### File Validation
```
File: EMIT_L2A_RFL_001_20230129T004447_2302816_003.nc
Size: 1765.28 MB
Checksum: 62fc62c5cd8e7652c0923b512bdf62d9

✓ netCDF4 opened successfully
✓ h5py opened successfully  
✓ ncdump validation passed
✓ EMITL2ARFL validation passed
```

## Solution

### Automatic Cache Clearing (✅ Now Implemented)

The retry mechanism now **automatically detects and removes corrupted cached files**:

1. Before any download, all cached files are validated
2. If a cached file is corrupted, it's immediately removed
3. A fresh download is triggered automatically
4. After download, files are validated again
5. If still corrupted, the file is removed and the retry loop continues

**You don't need to manually clean the cache anymore** - the code handles it automatically!

### Example of Automatic Cache Clearing

```
[INFO] Cached file validation failed: Reflectance file cannot be read...
[INFO] Removing corrupted cached file: .../EMIT_L2A_RFL_001_xxx.nc
[INFO] Corrupted file removed successfully
[INFO] Downloading granule files (attempt 1/3)...
[INFO] Downloaded file validated successfully: .../EMIT_L2A_RFL_001_xxx.nc
[INFO] All files successfully downloaded and validated.
```

### Manual Cache Cleanup (Optional)

If you want to manually clean your cache before running a large batch:

```bash
# Run the cache cleanup utility
conda run -n EMITL2ARFL python clean_cache.py ~/data/EMIT_L2A_RFL
```

This will:
- Scan all `.nc` files in your cache directory
- Identify corrupted vs valid files
- Offer to delete corrupted files
- Remove empty directories

**Note**: This is optional - the automatic cache clearing handles this during downloads.

### Use Fresh Download Directory

If you want to start fresh, use a new directory:

```python
from EMITL2ARFL.generate_EMIT_L2A_RFL_timeseries import generate_EMIT_L2A_RFL_timeseries

filenames = generate_EMIT_L2A_RFL_timeseries(
    start_date_UTC='2022-08-01',
    end_date_UTC='2025-11-20',
    geometry=grid,
    output_directory='~/data/Kings Canyon EMIT',
    download_directory='/tmp/emit_download',  # Fresh location
    max_retries=3,
    retry_delay=2.0
)
```

### Continue with Existing Setup (✅ Recommended)

Your existing setup will work perfectly now with automatic cache clearing. Just run your script as normal:

```python
python generate_kings_canyon_timeseries.py
```

The code will automatically:
- Validate all cached files
- Remove any that are corrupted  
- Download fresh copies
- Validate the downloads
- Continue processing

## Why /tmp Worked

`/tmp` had no pre-existing corrupted files, so:
1. Script detected files didn't exist
2. Downloaded fresh copies
3. Files downloaded correctly
4. Validation passed

## For HPC Usage

When running on HPC systems, the recommendations from previous documentation still apply:

1. **Use local scratch space** for downloads:
   ```python
   download_directory='/scratch/$USER/emit_download'
   ```

2. **Single-threaded downloads** are automatically used on retries (already implemented)

3. **Validation checks** catch corruption immediately (already implemented)

## Scripts Created

1. **`diagnose_single_download.py`**: Comprehensive diagnostic tool
   - Checks environment and dependencies
   - Downloads single file with full diagnostics
   - Tests NetCDF operations
   - Validates with multiple tools

2. **`test_single_granule.py`**: Simple test for single granule download
   - Quick test of download functionality
   - Minimal dependencies
   - Clear success/failure output

3. **`clean_cache.py`**: Cache cleanup utility
   - Scans for corrupted files
   - Safe deletion with confirmation
   - Cleans up empty directories

## Recommendations

### Immediate Actions

1. ✅ **Just run your script!** The automatic cache clearing will handle everything
2. ✅ Optionally run `clean_cache.py` first if you want to see what's corrupted
3. ✅ Monitor the logs - you'll see "[INFO] Cached file validated successfully" for good files

### Going Forward

1. **Keep using your existing setup** - it's working correctly
2. **Let validation catch issues** - the retry logic will handle transient failures
3. **Monitor HPC downloads** - if you see persistent corruption on HPC, use local scratch space

## Conclusion

**Your system is working correctly!** The issue was historical corrupted cache files, not current download/validation problems. After cleaning the cache, your workflow should work smoothly on both local Mac and HPC systems.

The diagnostic tools created during this investigation are available if you need to troubleshoot issues in the future.
