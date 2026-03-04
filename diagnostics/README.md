# EMIT L2A RFL Diagnostics

This folder contains diagnostic tools, test scripts, and troubleshooting documentation for the EMITL2ARFL package, particularly focused on HPC system compatibility and HDF5/NetCDF issues.

## Quick Start

### Installation Verification

```bash
# Check if your installation is properly configured
python verify_installation.py

# Check HDF5/NetCDF versions and test a file
python check_hdf5_version.py /path/to/file.nc
```

### Environment Setup

```bash
# Run with proper environment settings (bash)
bash test_hpc.sh

# Verify HPC setup
bash verify_hpc_setup.sh
```

## Documentation

### Installation Guides

- **[Install HPC.md](Install%20HPC.md)** - Complete guide for installing on HPC systems
- **[HPC Setup.md](HPC%20Setup.md)** - Detailed setup instructions with examples
- **[HPC Quickstart.md](HPC%20Quickstart.md)** - Quick reference for HPC users

### Troubleshooting

- **[HPC Troubleshooting.md](HPC%20Troubleshooting.md)** - Common issues and solutions
- **[HPC Fix.md](HPC%20Fix.md)** - Explanation of file locking issue and fixes
- **[Quick Fix.md](Quick%20Fix.md)** - Fast solutions for common problems

### Technical Details

- **[Package Update.md](Package%20Update.md)** - Automatic fix implementation details
- **[Cache Clearing Implementation.md](Cache%20Clearing%20Implementation.md)** - Cache management strategy
- **[Persistent Corruption Fix.md](Persistent%20Corruption%20Fix.md)** - Solutions for corruption issues
- **[Changes Summary.md](Changes%20Summary.md)** - Summary of diagnostic changes
- **[Diagnostic Results.md](Diagnostic%20Results.md)** - Test results documentation
- **[Readme Diagnostics.md](Readme%20Diagnostics.md)** - Overview of diagnostic tools

## Diagnostic Scripts

### Environment Verification

- **verify_installation.py** - Verify HDF5/NetCDF library configuration
- **check_hdf5_version.py** - Display library versions and test file access
- **debug_hpc_environment.py** - Comprehensive environment diagnostics with automatic file download

### File Testing

- **test_file_validation.py** - Test NetCDF file validation
- **test_single_granule.py** - Test downloading and validating a single granule
- **diagnose_netcdf.py** - Diagnose NetCDF file issues
- **diagnose_netcdf_environment.py** - Check NetCDF environment configuration
- **diagnose_single_download.py** - Test single file download

### HPC-Specific Tools

- **test_hpc.sh** - Test script for HPC systems (sets environment before Python)
- **test_hpc_env.sh** - Alternative HPC test wrapper
- **test_local_copy.sh** - Compare network vs local filesystem performance
- **run_with_env.sh** - Run Python with proper HDF5 environment settings
- **verify_hpc_setup.sh** - Verify HPC setup is correct
- **test_file_validation_hpc.py** - HPC-specific file validation test
- **test_automatic_hpc_fix.py** - Test automatic environment variable setting

### Utilities

- **clean_cache.py** - Clear cached downloaded files
- **download_to_scratch.py** - Download to local scratch before copying to network storage

## Common Issues

### HDF Error -101

**Symptom:** `[Errno -101] NetCDF: HDF error` when opening files on HPC systems

**Cause:** Conflict between system HDF5 libraries and Python packages

**Solution:**
```bash
# Rebuild environment with conda-forge packages
mamba create -n EMITL2ARFL -c conda-forge python=3.10 hdf5 h5py netcdf4
mamba activate EMITL2ARFL
cd /path/to/EMIT-L2A-RFL
pip install -e .

# Set environment variable
set -Ux HDF5_USE_FILE_LOCKING FALSE  # fish
export HDF5_USE_FILE_LOCKING=FALSE   # bash (add to ~/.bashrc)
```

See [HPC Troubleshooting.md](HPC%20Troubleshooting.md) for more solutions.

### File Works Locally But Not on HPC

**Cause:** HDF5 version incompatibility between systems

**Solution:** Download files directly on the HPC system where they'll be used. Don't copy HDF5/NetCDF files between systems with different library versions.

### Files Appear Corrupted After Download

**Cause:** Network filesystem caching/buffering issues

**Solution:** Use `download_to_scratch.py` to download to local scratch first:
```bash
python download_to_scratch.py
```

## Running Tests

```bash
# Verify installation
python verify_installation.py

# Test with a specific file
python check_hdf5_version.py /path/to/file.nc

# Run comprehensive diagnostics
python debug_hpc_environment.py

# Test download and validation
python test_single_granule.py
```

## Support

For additional help:
1. Run `python verify_installation.py` to check your setup
2. Review [HPC Troubleshooting.md](HPC%20Troubleshooting.md)
3. Check [Install HPC.md](Install%20HPC.md) for installation instructions
4. Open an issue at https://github.com/STARS-Data-Fusion/EMIT-L2A-RFL/issues
