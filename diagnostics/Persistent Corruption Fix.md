# Fixing Persistent NetCDF Corruption on HPC Systems

## Problem: Files Keep Getting Corrupted

If you're seeing errors like this even after multiple retries:

```
Failed to download valid files after 3 attempts. Missing or corrupted files: [...]
[Errno -101] NetCDF: HDF error
```

This indicates a deeper issue with your HPC environment's storage or network.

## Root Causes

1. **Network Filesystem Issues**: NFS, Lustre, GPFS can buffer writes incorrectly
2. **Multi-threaded Download Conflicts**: Parallel downloads can corrupt files on some filesystems
3. **Storage System Problems**: Disk errors, RAID issues, or storage congestion
4. **Network Instability**: Packet loss during large file transfers

## Solutions (Try in Order)

### 1. Use Local Scratch Space (RECOMMENDED)

Network filesystems are often the problem. Download to local disk first:

```python
import shutil
from pathlib import Path

# Download to local scratch (fast, reliable)
local_scratch = "/tmp/emit_download"  # or "/scratch/$USER/emit"

filenames = generate_EMIT_L2A_RFL_timeseries(
    start_date_UTC=start_date_UTC,
    end_date_UTC=end_date_UTC,
    geometry=grid,
    download_directory=local_scratch,  # Local disk!
    output_directory=output_directory,  # Can still be network storage
    max_retries=5,
    retry_delay=5.0
)

# Files are processed directly from local scratch
# Output files go to network storage which is fine for writing
```

### 2. Force Single-Threaded Downloads

The new version automatically uses single-threaded downloads on retries, but you can force it from the start by modifying your local installation temporarily.

Add this to your script before downloading:

```python
import earthaccess

# Monkey-patch to force single-threaded downloads
_original_download = earthaccess.download

def _safe_download(granules, local_path=None, provider=None, threads=1, **kwargs):
    """Force single-threaded downloads"""
    return _original_download(granules, local_path, provider, threads=1, **kwargs)

earthaccess.download = _safe_download
```

### 3. Extreme Retry Settings

If you must use network storage:

```python
filenames = generate_EMIT_L2A_RFL_timeseries(
    start_date_UTC=start_date_UTC,
    end_date_UTC=end_date_UTC,
    geometry=grid,
    download_directory=download_directory,
    output_directory=output_directory,
    max_retries=10,      # Many more attempts
    retry_delay=15.0     # Much longer delays (15s, 30s, 60s, 120s...)
)
```

### 4. Download Outside Peak Hours

HPC storage systems are often congested during business hours:

```bash
# Schedule your job for off-peak hours
#SBATCH --begin=20:00  # Start at 8 PM
#SBATCH --time=12:00:00
```

### 5. Check Your HPC Storage

Run diagnostics:

```bash
# Check filesystem health
df -h /path/to/your/download/directory

# For Lustre
lfs df -h

# Check for errors in system logs
dmesg | grep -i error

# Test write speed
dd if=/dev/zero of=/path/to/test bs=1M count=1000
rm /path/to/test
```

### 6. Use Alternative Download Methods

Download files manually using wget or curl first, then process them:

```bash
# Get file URLs from earthaccess
python -c "
import earthaccess
earthaccess.login()
results = earthaccess.search_data(
    short_name='EMITL2ARFL',
    temporal=('2022-08-13', '2022-08-13')
)
for r in results:
    for url in r.data_links():
        print(url)
" > urls.txt

# Download with wget (more reliable on some systems)
wget --no-check-certificate -i urls.txt -P /path/to/download/
```

Then process the pre-downloaded files.

### 7. Skip Validation (LAST RESORT)

**WARNING**: Only use this if you're desperate. Files may still be corrupted and cause crashes later.

```python
from EMITL2ARFL import retrieve_EMIT_L2A_RFL_granule

# This will attempt to use files even if they're corrupted
granule = retrieve_EMIT_L2A_RFL_granule(
    remote_granule=remote_granule,
    download_directory=download_directory,
    skip_validation=True  # DANGEROUS!
)
```

## HPC-Specific Configurations

### For Slurm on Lustre

```bash
#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=4

# Use local scratch
export TMPDIR=/local/scratch/$SLURM_JOB_ID
mkdir -p $TMPDIR

# Set download directory
export EMIT_DOWNLOAD_DIR=$TMPDIR/emit

# Run script
conda activate EMITL2ARFL
python generate_kings_canyon_timeseries.py

# Copy results to permanent storage
cp -r $TMPDIR/output/* /gpfs/home/$USER/results/
```

### For PBS on NFS

```bash
#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l mem=64gb

# Use node-local scratch
cd $PBS_O_WORKDIR
export EMIT_DOWNLOAD_DIR=/scratch/$USER/emit_$$
mkdir -p $EMIT_DOWNLOAD_DIR

python generate_kings_canyon_timeseries.py

# Clean up
rm -rf $EMIT_DOWNLOAD_DIR
```

## Checking if Local Scratch Helped

Compare download reliability:

```python
import time

# Test network storage
start = time.time()
try:
    granule = retrieve_EMIT_L2A_RFL_granule(
        orbit=22251, scene=7,
        download_directory="/gpfs/home/user/data"  # Network
    )
    print(f"Network: SUCCESS in {time.time()-start:.1f}s")
except Exception as e:
    print(f"Network: FAILED - {e}")

# Test local scratch
start = time.time()
try:
    granule = retrieve_EMIT_L2A_RFL_granule(
        orbit=22251, scene=7,
        download_directory="/tmp/emit"  # Local
    )
    print(f"Local: SUCCESS in {time.time()-start:.1f}s")
except Exception as e:
    print(f"Local: FAILED - {e}")
```

## When to Contact Support

Contact your HPC support team if:

1. Files corrupt even on local scratch storage
2. `dmesg` shows disk errors
3. Downloads work from login nodes but not compute nodes
4. Problem started recently (may indicate hardware failure)
5. Other users report similar issues

Include this information in your ticket:

```bash
# Gather diagnostic info
uname -a
df -h
lfs df -h  # If using Lustre
dmesg | tail -100
cat /proc/meminfo | head -20

# Include error messages from your Python script
python generate_kings_canyon_timeseries.py 2>&1 | tee error.log
```

## Summary

**Quick Fix**: Use local scratch space for downloads:
```python
download_directory="/tmp/emit_download"  # or /scratch/$USER/emit
```

This solves the problem in 90% of cases.
