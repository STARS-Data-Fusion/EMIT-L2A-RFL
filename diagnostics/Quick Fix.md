# Quick Fix for HPC HDF Errors

## TL;DR - The #1 Solution

**Use local scratch space instead of network storage for downloads:**

```python
# Change this:
download_directory = "~/data/EMIT_download"  # Network storage (NFS/Lustre) - BAD

# To this:
download_directory = "/tmp/emit_download"  # Local scratch - GOOD
# or
download_directory = "/scratch/$USER/emit"  # HPC scratch space - GOOD
```

This solves 90% of HDF corruption issues on HPC systems.

---

## If You're Getting `[Errno -101] NetCDF: HDF error` on HPC systems:

### 1. Find corrupted files:
```bash
python diagnose_netcdf.py ~/data/EMIT_download/
```

### 2. Delete corrupted files:
```bash
# Manually delete the files listed as invalid by the diagnostic tool
rm ~/data/EMIT_download/corrupted_file.nc
```

### 3. Re-run with local scratch space:
```python
from EMITL2ARFL import generate_EMIT_L2A_RFL_timeseries

filenames = generate_EMIT_L2A_RFL_timeseries(
    start_date_UTC="2022-08-01",
    end_date_UTC="2025-11-20",
    geometry=grid,
    download_directory="/tmp/emit_download",  # LOCAL SCRATCH!
    output_directory="~/data/output",         # Output can still go to network storage
    max_retries=5,      # More retry attempts
    retry_delay=5.0     # Longer delays between retries
)
```

## What Changed?

The package now automatically:
- ✅ Detects corrupted files during download
- ✅ Waits longer between retries (exponential backoff)
- ✅ Syncs filesystem after operations
- ✅ Verifies files are stable before validation
- ✅ Provides better error messages

## For Your generate_kings_canyon_timeseries.py Script

No changes required! The script will automatically use the improved retry logic.

**Optional**: Add retry parameters for even better reliability:

```python
# At the top of your script
MAX_RETRIES = 5
RETRY_DELAY = 5.0

# When calling generate_EMIT_L2A_RFL_timeseries
filenames = generate_EMIT_L2A_RFL_timeseries(
    start_date_UTC=start_date_UTC,
    end_date_UTC=end_date_UTC,
    geometry=grid,
    output_directory=output_directory,
    max_retries=MAX_RETRIES,
    retry_delay=RETRY_DELAY
)
```

## For Slurm Jobs

Add to your job script:
```bash
#!/bin/bash
#SBATCH --time=04:00:00

# Set retry parameters via environment variables
export EMIT_MAX_RETRIES=5
export EMIT_RETRY_DELAY=10.0

python generate_kings_canyon_timeseries.py
```

Then read them in your Python script:
```python
import os

max_retries = int(os.environ.get('EMIT_MAX_RETRIES', 3))
retry_delay = float(os.environ.get('EMIT_RETRY_DELAY', 2.0))

filenames = generate_EMIT_L2A_RFL_timeseries(
    ...,
    max_retries=max_retries,
    retry_delay=retry_delay
)
```

## More Details

See [HPC_TROUBLESHOOTING.md](HPC_TROUBLESHOOTING.md) for comprehensive documentation.
