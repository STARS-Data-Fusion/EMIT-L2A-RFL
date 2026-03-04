# Package Update: Automatic HDF5 File Locking Fix

## What Changed

The EMITL2ARFL package now **automatically handles the HDF5 file locking issue** on HPC systems. No user configuration is required!

### Before This Update

Users had to manually set the environment variable:
```bash
export HDF5_USE_FILE_LOCKING=FALSE
python your_script.py
```

### After This Update

Just import and use the package:
```python
from EMITL2ARFL import validate_NetCDF_file

# Works automatically on HPC systems!
validate_NetCDF_file("~/data/file.nc")
```

## How It Works

When you import EMITL2ARFL, the package automatically sets `HDF5_USE_FILE_LOCKING=FALSE` **before** importing netCDF4 or any other HDF5-dependent libraries. This happens at package initialization time, ensuring the setting takes effect.

## Backwards Compatibility

- ✅ If you've already set `HDF5_USE_FILE_LOCKING` in your environment, your setting is respected
- ✅ All existing scripts continue to work without modification
- ✅ Works on both local machines and HPC systems
- ✅ No breaking changes

## For HPC Users

You can now simplify your workflow:

### Old Way (still works)
```bash
export HDF5_USE_FILE_LOCKING=FALSE
conda activate EMITL2ARFL
python generate_kings_canyon_timeseries.py
```

### New Way (simpler)
```bash
conda activate EMITL2ARFL
python generate_kings_canyon_timeseries.py  # Just works!
```

### For SLURM Jobs

You can remove the environment variable from job scripts if you want:

```bash
#!/bin/bash
#SBATCH --job-name=emit_processing

# No longer strictly necessary (but doesn't hurt to keep it)
# export HDF5_USE_FILE_LOCKING=FALSE

conda activate EMITL2ARFL
python generate_kings_canyon_timeseries.py
```

## Technical Details

The fix is implemented in `/EMITL2ARFL/__init__.py`:

```python
import os
import sys

# Set HDF5_USE_FILE_LOCKING before any netCDF4/HDF5 imports
if 'HDF5_USE_FILE_LOCKING' not in os.environ:
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

# Check if netCDF4 was already imported
if 'netCDF4' in sys.modules:
    import warnings
    warnings.warn(
        "netCDF4 was already imported before EMITL2ARFL initialization...",
        RuntimeWarning
    )
```

## Edge Cases

### If you import netCDF4 before EMITL2ARFL

```python
import netCDF4  # Imported first
from EMITL2ARFL import validate_NetCDF_file  # Shows warning
```

You'll see a warning that the fix may not work. Solution: Import EMITL2ARFL first.

### If you want to force file locking ON

Set it explicitly before importing:
```python
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'TRUE'
from EMITL2ARFL import validate_NetCDF_file  # Will respect your setting
```

## Verification

Test that it's working:

```python
import os
print("Before:", os.environ.get('HDF5_USE_FILE_LOCKING'))
# Output: Before: None

from EMITL2ARFL import validate_NetCDF_file
print("After:", os.environ.get('HDF5_USE_FILE_LOCKING'))
# Output: After: FALSE
```

## Recommendation

While the package now handles this automatically, it's still good practice to set the environment variable in your HPC environment (`~/.bashrc`) as a defense-in-depth measure, especially if you use other packages that access HDF5/NetCDF files outside of EMITL2ARFL.

## Migration Guide

### If you were using manual workarounds:

1. ✅ Your existing scripts continue to work
2. ✅ You can remove manual environment variable setting if you want
3. ✅ Job scripts can be simplified
4. ✅ No code changes required

### If you were using test_hpc.sh or other wrappers:

They still work, but are no longer strictly necessary for scripts that only use EMITL2ARFL.

## Summary

🎉 **The HPC file locking issue is now fixed automatically!**

- No configuration needed
- Works transparently
- Backwards compatible
- One less thing to worry about when deploying to HPC systems

Just `pip install` or `conda install` the updated package and everything works.
