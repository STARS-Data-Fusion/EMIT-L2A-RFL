# Production Workflow Guide for HPC

Now that you've confirmed the fix works on your HPC machine, here's how to apply it to your actual EMIT data processing workflows.

## Confirmed Working Solution

✅ `./test_hpc.sh` passed on HPC machine

This confirms that setting `HDF5_USE_FILE_LOCKING=FALSE` at the shell level before Python starts resolves the issue.

## Apply to Your Processing Workflows

### For Interactive Work

Add to your `~/.bashrc` on the HPC machine (if not already done):

```bash
export HDF5_USE_FILE_LOCKING=FALSE
```

Then reload:
```bash
source ~/.bashrc
```

Now all your Python scripts will work normally:
```bash
conda activate EMITL2ARFL
python generate_kings_canyon_timeseries.py
```

### For Batch Jobs (SLURM)

Update your job scripts to include the environment variable:

```bash
#!/bin/bash
#SBATCH --job-name=emit_timeseries
#SBATCH --time=08:00:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=4

# CRITICAL: Set this for network filesystems
export HDF5_USE_FILE_LOCKING=FALSE

# Load environment
conda activate EMITL2ARFL

# Run processing
python generate_kings_canyon_timeseries.py
```

### For PBS/Torque Jobs

```bash
#!/bin/bash
#PBS -N emit_timeseries
#PBS -l walltime=08:00:00
#PBS -l mem=64GB
#PBS -l ncpus=4

# CRITICAL: Set this for network filesystems
export HDF5_USE_FILE_LOCKING=FALSE

# Load environment
conda activate EMITL2ARFL

# Change to working directory
cd $PBS_O_WORKDIR

# Run processing
python generate_kings_canyon_timeseries.py
```

## Verify Your Setup

Before running large processing jobs, verify the environment:

```bash
# Check environment variable
echo $HDF5_USE_FILE_LOCKING  # Should print: FALSE

# Test with a single file
python test_file_validation.py

# Or use the wrapper
./test_hpc.sh
```

## Example: Processing Upper Kings Canyon Data

```bash
#!/bin/bash
#SBATCH --job-name=upper_kings
#SBATCH --time=12:00:00
#SBATCH --mem=64GB
#SBATCH --output=upper_kings_%j.out
#SBATCH --error=upper_kings_%j.err

# Set environment variable
export HDF5_USE_FILE_LOCKING=FALSE

# Activate conda environment
conda activate EMITL2ARFL

# Run the processing
python generate_kings_canyon_timeseries.py

echo "Processing complete!"
```

Submit with:
```bash
sbatch process_upper_kings.sh
```

## For Jupyter Notebooks on HPC

If using Jupyter on the HPC system, set the environment variable before starting Jupyter:

```bash
export HDF5_USE_FILE_LOCKING=FALSE
conda activate EMITL2ARFL
jupyter lab --no-browser --port=8888
```

Or add it to your Jupyter kernel's environment. Edit:
```
~/.local/share/jupyter/kernels/emitl2arfl/kernel.json
```

Add to the env section:
```json
{
  "argv": [...],
  "display_name": "EMITL2ARFL",
  "language": "python",
  "env": {
    "HDF5_USE_FILE_LOCKING": "FALSE"
  }
}
```

## What You Can Remove/Ignore

Now that you have the proper fix, you can:

1. **Ignore `test_file_validation_hpc.py`** - It was an attempt to set the variable in Python, which doesn't work reliably
2. **Use standard scripts** - All your existing scripts will work once the environment variable is set at the shell level
3. **Focus on `test_hpc.sh`** - Use this as a quick test whenever you're unsure if the environment is properly configured

## Quick Reference Commands

```bash
# One-time setup (if not in ~/.bashrc)
export HDF5_USE_FILE_LOCKING=FALSE

# Verify setup
echo $HDF5_USE_FILE_LOCKING  # Should print: FALSE
./test_hpc.sh                # Should pass all tests

# Run your scripts normally
conda activate EMITL2ARFL
python generate_kings_canyon_timeseries.py
python generate_EMIT_L2A_RFL_timeseries.py
# etc.
```

## Files You Need on HPC

Essential files to upload to your HPC system:

1. ✅ **test_hpc.sh** - For quick environment verification
2. ✅ **debug_hpc_environment.py** - For troubleshooting if issues arise
3. ✅ **HPC_QUICKSTART.md** or **HPC_FIX.md** - For reference

Optional (for deeper troubleshooting):
- HPC_SETUP.md
- test_file_validation_hpc.py

## Success Indicators

You'll know everything is working when:

- ✅ `./test_hpc.sh` passes
- ✅ `python test_file_validation.py` succeeds
- ✅ Your timeseries generation scripts complete without HDF errors
- ✅ No more `errno -101` or "NetCDF: HDF error" messages

## Next Steps

1. **Add to .bashrc** (if not already done):
   ```bash
   echo 'export HDF5_USE_FILE_LOCKING=FALSE' >> ~/.bashrc
   source ~/.bashrc
   ```

2. **Update your job scripts** to include `export HDF5_USE_FILE_LOCKING=FALSE`

3. **Test your actual workflows:**
   ```bash
   python generate_kings_canyon_timeseries.py
   ```

4. **Monitor your jobs** to ensure they complete successfully

You should now be able to run all your EMIT data processing workflows on the HPC system without any HDF/NetCDF errors! 🎉
