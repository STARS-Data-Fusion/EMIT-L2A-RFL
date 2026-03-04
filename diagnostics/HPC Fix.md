# HPC Testing Instructions

## What You Told Me

- **Local machine:** test_file_validation.py works ✓
- **HPC machine:** 
  - test_file_validation.py fails with HDF error ✗
  - debug_hpc_environment.py works ✓
  - iPython manual testing works ✓

## Root Cause

The issue is that `HDF5_USE_FILE_LOCKING=FALSE` must be set **at the shell level BEFORE Python starts**. 

- `debug_hpc_environment.py` works because it sets the variable, imports netCDF4 directly, then imports EMITL2ARFL
- `test_file_validation.py` and `test_file_validation_hpc.py` fail because EMITL2ARFL imports netCDF4 during package initialization, which happens before the script can set the environment variable
- iPython works because your shell configuration already has the variable set

## The Fix for HPC

On your HPC machine, run:

```bash
# 1. Set the environment variable at shell level
export HDF5_USE_FILE_LOCKING=FALSE

# 2. Activate environment
conda activate EMITL2ARFL

# 3. Run your script
python test_file_validation.py
```

Or in one command:
```bash
HDF5_USE_FILE_LOCKING=FALSE python test_file_validation.py
```

## Make It Permanent

Add this to your `~/.bashrc` on the HPC system:

```bash
echo 'export HDF5_USE_FILE_LOCKING=FALSE' >> ~/.bashrc
source ~/.bashrc
```

## Testing on HPC

### Option 1: Use the wrapper script (easiest)
```bash
conda activate EMITL2ARFL
./test_hpc.sh
```

### Option 2: Set manually
```bash
export HDF5_USE_FILE_LOCKING=FALSE
conda activate EMITL2ARFL
python test_file_validation.py
```

### Option 3: One-liner
```bash
HDF5_USE_FILE_LOCKING=FALSE python test_file_validation.py
```

## Verification

Check if the variable is set:
```bash
echo $HDF5_USE_FILE_LOCKING  # Should print: FALSE
```

Check in Python:
```bash
python -c "import os; print(os.environ.get('HDF5_USE_FILE_LOCKING'))"
```

## For Your Processing Scripts

Update your SLURM/PBS job scripts to include:

```bash
#!/bin/bash
#SBATCH --job-name=emit_processing
#SBATCH --time=04:00:00
#SBATCH --mem=32GB

# CRITICAL: Set this before running Python
export HDF5_USE_FILE_LOCKING=FALSE

# Activate environment
conda activate EMITL2ARFL

# Run your processing
python generate_kings_canyon_timeseries.py
```

## Why Setting It in Python Doesn't Work

This approach fails on HPC:

```python
# test_file_validation_hpc.py
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'  # TOO LATE!
from EMITL2ARFL import validate_NetCDF_file      # netCDF4 already initialized
```

By the time you set the variable in Python, the EMITL2ARFL package has already imported netCDF4, and the HDF5 library has already initialized with file locking enabled.

**The environment variable must be set at the shell level before Python starts.**

## Summary for HPC Use

1. ✅ Add `export HDF5_USE_FILE_LOCKING=FALSE` to `~/.bashrc`
2. ✅ Source your bashrc: `source ~/.bashrc`
3. ✅ Verify: `echo $HDF5_USE_FILE_LOCKING`
4. ✅ Run scripts normally: `python test_file_validation.py`

That's it! All your scripts will now work on the HPC system.
