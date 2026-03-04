# Installing EMITL2ARFL on HPC Systems

This guide ensures proper HDF5/NetCDF library configuration on HPC systems.

## Problem

HPC systems often have system-wide HDF5/NetCDF libraries that conflict with Python packages. This causes `HDF Error -101` and other read failures.

## Solution: Use Conda (Recommended)

Conda-forge builds include compatible, bundled HDF5 libraries that don't conflict with system libraries.

### 1. Create Clean Conda Environment

```bash
# Create new environment
conda create -n EMITL2ARFL python=3.10

# Activate it
conda activate EMITL2ARFL

# Install HDF5/NetCDF from conda-forge FIRST
conda install -c conda-forge hdf5 h5py netcdf4

# Now install EMITL2ARFL (will use conda's HDF5)
pip install git+https://github.com/STARS-Data-Fusion/EMIT-L2A-RFL.git

# Or for local development:
cd /path/to/EMIT-L2A-RFL
pip install -e .
```

### 2. Set Environment Variable

Add to your `~/.bashrc` or `~/.config/fish/config.fish`:

**For bash:**
```bash
export HDF5_USE_FILE_LOCKING=FALSE
```

**For fish:**
```fish
set -Ux HDF5_USE_FILE_LOCKING FALSE
```

Then log out and back in.

## Alternative: Pip-Only Installation (Not Recommended for HPC)

If you must use pip-only:

```bash
# Create virtual environment
python3 -m venv ~/venvs/EMITL2ARFL
source ~/venvs/EMITL2ARFL/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install with explicit HDF5 from PyPI
pip install h5py netCDF4
pip install git+https://github.com/STARS-Data-Fusion/EMIT-L2A-RFL.git
```

**Important:** Pip-installed netCDF4 may still link against system HDF5, causing conflicts.

## Verify Installation

```bash
python check_hdf5_version.py
```

Should show:
- `h5py` and `netCDF4` from your conda/venv (not system paths)
- `HDF5_USE_FILE_LOCKING: FALSE`

## Test with Sample File

```bash
# Download test file
python -c "
import earthaccess
from EMITL2ARFL import search_EMIT_L2A_RFL_granules, retrieve_EMIT_L2A_RFL_granule

earthaccess.login()
granules = search_EMIT_L2A_RFL_granules(
    start_UTC='2023-01-29',
    end_UTC='2023-01-30',
    orbit=2302816,
    scene=3
)

granule = retrieve_EMIT_L2A_RFL_granule(
    remote_granule=granules[0],
    download_directory='/tmp/emit_test'
)
print(f'Success: {granule.reflectance_filename}')
"
```

## Troubleshooting

### Still Getting HDF Error -101?

1. Check library sources:
   ```bash
   python -c "import netCDF4; print(netCDF4.__file__)"
   python -c "import h5py; print(h5py.__file__)"
   ```
   
   Should show paths in your conda env or venv, NOT `/usr/lib` or `/lib64`.

2. Check HDF5 version compatibility:
   ```bash
   python check_hdf5_version.py /path/to/file.nc
   ```

3. If h5py works but netCDF4 fails:
   ```bash
   conda remove netcdf4
   conda install -c conda-forge netcdf4
   ```

4. Nuclear option - rebuild environment:
   ```bash
   conda deactivate
   conda env remove -n EMITL2ARFL
   # Start over from step 1
   ```

### Files Downloaded on Mac Don't Work on HPC?

This is normal - HDF5 files can be version-specific. **Solution:**

- Download files directly on the HPC system where they'll be used
- Don't copy HDF5/NetCDF files between systems with different library versions
- If you must share files, use systems with matching HDF5 versions

## For System Administrators

If managing shared HPC environments:

```bash
# Load compatible modules before installing
module purge
module load python/3.10
module load hdf5/1.14
module load netcdf/4.9

# Or provide conda/mamba for users
module load conda
```
