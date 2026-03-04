# Quick Start Guide for HPC Users

If you're experiencing HDF/NetCDF errors on an HPC system (but the same file works in iPython), follow these steps:

## The Problem

**Symptom:** Your script fails with an HDF error, but the file opens fine when you manually test it in iPython:

```python
# This works in iPython:
import netCDF4
ds = netCDF4.Dataset("~/data/file.nc", 'r')  # ✓ Works!

# But this fails in a script:
python test_file_validation.py  # ✗ HDF Error -101
```

**Root Cause:** HDF5 file locking doesn't work on network filesystems (NFS/Lustre/GPFS).

## The Solution

Set this environment variable **BEFORE** running Python:

```bash
export HDF5_USE_FILE_LOCKING=FALSE
```

### Step-by-Step Fix

1. **Add to your shell configuration** (`~/.bashrc` or `~/.bash_profile`):
   ```bash
   echo 'export HDF5_USE_FILE_LOCKING=FALSE' >> ~/.bashrc
   source ~/.bashrc
   ```

2. **Activate your environment:**
   ```bash
   conda activate EMITL2ARFL
   ```

3. **Run your scripts:**
   ```bash
   python test_file_validation.py  # Should now work!
   ```

### Quick Test

Use the provided wrapper script to test:

```bash
conda activate EMITL2ARFL
./test_hpc.sh
```

If this succeeds, your environment is properly configured.

## Diagnostic Tool

If you're still having issues, run the diagnostic:

```bash
python debug_hpc_environment.py ~/data/EMIT_L2A_RFL_*.nc
```

This will tell you exactly what's wrong and provide specific recommendations.

## For SLURM Jobs

Add to the top of your job script:

```bash
#!/bin/bash
#SBATCH --job-name=emit_process
#SBATCH --time=04:00:00

# Critical for network filesystems
export HDF5_USE_FILE_LOCKING=FALSE

conda activate EMITL2ARFL
python your_script.py
```

## Why This Happens

- **Network filesystems** (NFS, Lustre, GPFS) don't support proper file locking
- **HDF5 library** tries to use file locking by default for data integrity
- **iPython** might have the environment variable already set from your shell config
- **Scripts** start fresh without the environment variable

## Permanent Fix

Make the environment variable automatic for your conda environment:

```bash
conda activate EMITL2ARFL

# Create activation script
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
echo 'export HDF5_USE_FILE_LOCKING=FALSE' > $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh

# Optional: Create deactivation script
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
echo 'unset HDF5_USE_FILE_LOCKING' > $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh
```

Now the environment variable is automatically set whenever you activate the environment.

## Verify It's Working

```bash
# Check the environment variable is set
echo $HDF5_USE_FILE_LOCKING  # Should print: FALSE

# Check in Python
python -c "import os; print(os.environ.get('HDF5_USE_FILE_LOCKING'))"  # Should print: FALSE
```

## Common Mistakes

❌ **Setting it in Python code:**
```python
# This is often TOO LATE!
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
from EMITL2ARFL import validate_NetCDF_file
```

✅ **Setting it at shell level:**
```bash
# This is the RIGHT way!
export HDF5_USE_FILE_LOCKING=FALSE
python my_script.py
```

## Still Not Working?

1. **Run the diagnostic:** `python debug_hpc_environment.py <file>`
2. **Check the file exists:** `ls -lh ~/data/*.nc`
3. **Verify environment:** `conda env list`
4. **Check disk space:** `df -h ~/data`
5. **See full guide:** [HPC_SETUP.md](HPC_SETUP.md)

## Summary

```bash
# Quick fix (one-time):
export HDF5_USE_FILE_LOCKING=FALSE
conda activate EMITL2ARFL
python test_file_validation.py

# Permanent fix:
echo 'export HDF5_USE_FILE_LOCKING=FALSE' >> ~/.bashrc
source ~/.bashrc
```
