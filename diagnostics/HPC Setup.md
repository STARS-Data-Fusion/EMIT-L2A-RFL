# HPC Environment Setup Guide

## Common Issue: Script Fails but iPython Works

If you can open NetCDF files in iPython but your scripts fail with HDF errors, this is almost always caused by **HDF5 file locking on network filesystems**.

### Why This Happens

- **Network filesystems** (NFS, Lustre, GPFS) often don't support proper file locking
- **HDF5 library** tries to use file locking by default for data safety
- **iPython** may have different environment variables or initialization order
- **Scripts** import libraries in a specific order that may trigger locking issues

### Solution 1: Disable File Locking (Recommended)

**CRITICAL:** The environment variable MUST be set BEFORE Python starts, not within Python code. Setting it in Python code is often too late because the HDF5 library initializes when netCDF4 is first imported.

#### Option A: Set Environment Variable at Shell Level (Most Reliable)

Add this to your `~/.bashrc` or `~/.bash_profile` (or `~/.zshrc` for zsh):

```bash
export HDF5_USE_FILE_LOCKING=FALSE
```

Then reload your shell:
```bash
source ~/.bashrc
```

#### Option B: Set Before Running Script (For Testing)

Set it in the same command line BEFORE running Python:
```bash
HDF5_USE_FILE_LOCKING=FALSE python test_file_validation.py
```

**Important:** This will NOT work:
```bash
# WRONG - too late!
python -c "import os; os.environ['HDF5_USE_FILE_LOCKING']='FALSE'; exec(open('test_file_validation.py').read())"
```

#### Option C: Use the HPC Test Wrapper (Easiest)

```bash
conda activate EMITL2ARFL
./test_hpc.sh
```

This wrapper script ensures the environment variable is set at the shell level before Python starts.

#### Option D: Add to Job Scripts
```bash
#!/bin/bash
#SBATCH --job-name=emit_process
#SBATCH --time=04:00:00
#SBATCH --mem=32GB

# Set environment variable
export HDF5_USE_FILE_LOCKING=FALSE

# Activate environment
conda activate EMITL2ARFL

# Run script
python generate_kings_canyon_timeseries.py
```

For PBS/Torque jobs:
```bash
#!/bin/bash
#PBS -N emit_process
#PBS -l walltime=04:00:00
#PBS -l mem=32GB

export HDF5_USE_FILE_LOCKING=FALSE
conda activate EMITL2ARFL
cd $PBS_O_WORKDIR
python generate_kings_canyon_timeseries.py
```

### Solution 2: Use the HPC Test Wrapper

The [test_hpc.sh](test_hpc.sh) wrapper ensures proper environment setup:

```bash
conda activate EMITL2ARFL
./test_hpc.sh
```

This script:
- Sets `HDF5_USE_FILE_LOCKING=FALSE` at shell level (before Python starts)
- Verifies the conda environment is activated
- Runs the test and provides clear diagnostics

### Solution 3: Diagnose Your Specific Issue

Run the diagnostic script to identify the exact problem:

```bash
python debug_hpc_environment.py ~/data/EMIT_L2A_RFL_001_20230129T004447_2302816_003.nc
```

This will:
- Check file permissions and filesystem type
- Verify library versions
- Test different access methods
- Identify if file locking is the issue
- Provide specific recommendations

### Why iPython Works But Scripts Don't

When you start iPython, you might:
1. Have `HDF5_USE_FILE_LOCKING=FALSE` already set in your `.bashrc`
2. Import packages in a different order
3. Have Setting in Python Code Often Doesn't Work

You might see code that tries to set the environment variable in Python:

```python
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
```

**This often fails on HPC systems because:**

1. The HDF5 library initializes when netCDF4 is first imported
2. If any module has already imported netCDF4 (directly or indirectly), it's too late
3. The EMITL2ARFL package imports netCDF4 during its initialization
4. By the time your script sets the variable, netCDF4 is already loaded

**The solution:** Set the environment variable at the **shell level** before Python starts, not within Python code.

### Why the environment variable set from a previous session
4. Be using a different conda environment

To check your iPython environment:

```python
import os
print(os.environ.get('HDF5_USE_FILE_LOCKING', 'not set'))
```

### Verification Steps

After applying the fix (make sure conda environment is activated first!):

```bash
conda activate EMITL2ARFL
```

1. **Test the original script:**
   ```bash
   python test_file_validation.py
   ```

2. **Test with the HPC-aware version:**
   ```bash
   python test_file_validation_hpc.py
   ```

3. **Verify in iPython:**
   ```bash
   ipython
   ```
   ```python
   import os
   os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
   from EMITL2ARFL import validate_NetCDF_file
   validate_NetCDF_file("~/data/EMIT_L2A_RFL_001_20230129T004447_2302816_003.nc")
   ```

### Other Potential Issues

If file locking isn't the problem, check:

1. **File corruption:**
   ```bash
   md5sum ~/data/EMIT_L2A_RFL_001_20230129T004447_2302816_003.nc
   ```
   Compare with the checksum from NASA (if available).

2. **Library version mismatch:**
   ```bash
   python -c "import netCDF4; print(netCDF4.__version__)"
   python -c "import h5py; print(h5py.version.hdf5_version)"
   ```

3. **Filesystem issues:**
   ```bash
   df -h ~/data  # Check disk space
   df -T ~/data  # Check filesystem type
   ```

4. **File permissions:**
   ```bash
   ls -lh ~/data/EMIT_L2A_RFL_001_20230129T004447_2302816_003.nc
   ```

### System-Specific Notes

#### For Lustre Filesystems
Lustre often has issues with HDF5 file locking. Always use:
```bash
export HDF5_USE_FILE_LOCKING=FALSE
```

#### For NFS Filesystems
Similar to Lustre. You might also need:
```bash
export HDF5_USE_FILE_LOCKING=FALSE
export NETCDF4_USE_ZLIB=FALSE  # If you see compression errors
```

#### For GPFS Filesystems
Usually works better, but still recommended:
```bash
export HDF5_USE_FILE_LOCKING=FALSE
```

### Adding to Your Conda Environment

To make this permanent for your EMITL2ARFL environment:

```bash
conda activate EMITL2ARFL

# Create activation script
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
echo 'export HDF5_USE_FILE_LOCKING=FALSE' > $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh

# Create deactivation script (optional, to clean up)
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
echo 'unset HDF5_USE_FILE_LOCKING' > $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh
```

Now the environment variable will be set automatically whenever you activate the environment.

### Quick Reference

| Problem | Solution |
|---------|----------|
| Script fails, iPython works | Set `HDF5_USE_FILE_LOCKING=FALSE` |
**First, activate the conda environment:**
```bash
conda activate EMITL2ARFL
```

**Then run all three test approaches:**

```bash
# 1. Original test (should now work)
python test_file_validation.py

# 2. HPC-aware test
python test_file_validation_hpc.py

# 3. Full diagnostic
python debug_hpc_environment.py ~/data/EMIT_L2A_RFL_001_20230129T004447_2302816_003.nc
```

All three should succeed after applying the fix.

### Common Mistake: Wrong Python Environment

If you see `ModuleNotFoundError: No module named 'netCDF4'`, you're using the wrong Python environment. Make sure to:

1. **Activate the conda environment first:**
   ```bash
   conda activate EMITL2ARFL
   ```

2. **Verify you're using the right Python:**
   ```bash
   which python  # Should show path to EMITL2ARFL environment
   ```

3. **Or use the full path:**
   ```bash
   # Find your conda environment path
   conda env list
   
   # Use full path (example)
   ~/miniconda3/envs/EMITL2ARFL/bin/python test_file_validation_hpc.py
   ```
python test_file_validation_hpc.py

# 3. Full diagnostic
python debug_hpc_environment.py ~/data/EMIT_L2A_RFL_001_20230129T004447_2302816_003.nc
```

All three should succeed after applying the fix.
