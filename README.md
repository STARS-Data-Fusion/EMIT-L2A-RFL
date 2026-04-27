# EMIT L2A Estimated Surface Reflectance and Uncertainty and Masks 60 m Search and Download Utility

This Python package is a tool for searching, downloading, and reading hypspectral surface reflectance from the Earth Surface Mineral Dust Source Investigation (EMIT) L2A RFL data product.

Gregory H. Halverson (they/them)<br>
[gregory.h.halverson@jpl.nasa.gov](mailto:gregory.h.halverson@jpl.nasa.gov)<br>
NASA Jet Propulsion Laboratory 329G<br>

## Installation

### Standard Installation (pip)

```bash
pip install EMITL2ARFL
```

### HPC Systems with HDF5 Issues

If you encounter HDF5 errors (like `HDF Error -101`) on HPC systems, the issue is typically caused by conflicts between system HDF5 libraries and pip-installed packages. The solution is to use conda-forge packages which bundle compatible HDF5 libraries.

**Using the makefile (recommended):**

```bash
git clone https://github.com/STARS-Data-Fusion/EMIT-L2A-RFL.git
cd EMIT-L2A-RFL
make environment  # Creates mamba environment with conda-forge HDF5
mamba activate EMITL2ARFL
make install
```

**Manual setup:**

```bash
# Create environment with HDF5 from conda-forge
mamba create -n EMITL2ARFL -c conda-forge python=3.10 hdf5 h5py netcdf4
mamba activate EMITL2ARFL

# Install package
pip install EMITL2ARFL
# Or for development: pip install -e .
```

**Set environment variable (HPC systems only):**

```bash
# For fish shell
set -Ux HDF5_USE_FILE_LOCKING FALSE

# For bash shell (add to ~/.bashrc)
export HDF5_USE_FILE_LOCKING=FALSE
```

See [diagnostics/Install HPC.md](diagnostics/Install%20HPC.md) for detailed HPC installation instructions and troubleshooting.

## Caveats

### Full Granule Download Requirement

The current version of this package **requires downloading entire EMIT L2A RFL granules** in order to obtain spatial subsets. While the package provides functions to extract data for a specific geometry or region of interest, the spatial subsetting occurs **after** the complete granule files (reflectance, mask, and uncertainty) have been downloaded to the local filesystem.

This means:
- **Full files are downloaded**: Each granule (~300-500 MB) is downloaded in its entirety, regardless of subset size
- **Local processing**: Spatial subsetting is performed locally after download using the geolocation information
- **Bandwidth impact**: For small areas of interest, this results in downloading significantly more data than necessary

**Workaround considerations:**
- Cache full granules locally if processing multiple overlapping geometries, as the same granule can be reused
- Consider the download bandwidth constraints when planning time series processing for large date ranges
- Granules that do not intersect your geometry are skipped automatically through the search and filtering process

Future versions of this package may support server-side subsetting via OGC Web Coverage Service (WCS) endpoints or cloud-native data formats (Zarr, Cloud-Optimized GeoTIFF) if such capabilities become available through NASA LPDAAC.

## Future Work: Server-Side Subsetting

Several promising approaches could enable efficient remote subsetting without downloading full granules:

### 1. HTTP Range Requests with Kerchunk (Most Feasible)
Use the `kerchunk` library to create virtual Zarr indices of remote NetCDF files, allowing byte-range requests for specific data chunks. This approach:
- Requires no infrastructure changes at NASA LPDAAC
- Leverages existing `earthaccess` authentication
- Downloads only the bytes needed for the requested subset
- Works immediately if EMIT data servers support HTTP Range requests (standard for most NASA data systems)

**Implementation**: Add optional `use_remote=True` parameter to retrieval functions that bypasses full download in favor of remote chunked access.

### 2. OGC Web Coverage Service (WCS)
If NASA LPDAAC exposes WCS endpoints for EMIT data, the package could:
- Submit server-side subsetting requests with geometry parameters
- Receive pre-subset data directly from the server
- Reduce bandwidth by orders of magnitude for small AOIs

**Implementation**: Create `retrieve_via_wcs.py` module using `owslib` library to construct and execute WCS GetCoverage requests.

### 3. Cloud-Native Data Formats
If EMIT data becomes available in cloud-optimized formats:
- **Zarr**: Native chunked format with efficient remote access via `xarray` and `fsspec`
- **Cloud-Optimized GeoTIFF (COG)**: GDAL can read remote tiles via HTTP Range requests
- Both formats enable subsetting at the storage level without custom protocols

**Implementation**: Add format detection and route requests to appropriate readers based on available data sources.

### 4. OpenDAP Protocol
If EMIT granules are served via OpenDAP (OPeNDAP) servers:
- Use `xarray.open_dataset()` with constraint expressions
- Request subsets by coordinate ranges directly in the URL
- Standard protocol supported across NASA data centers

**Implementation**: Modify `retrieve_EMIT_L2A_RFL_granule()` to accept OpenDAP URLs and use `xarray` remote access instead of download.

### 5. Enhanced Earthaccess Integration
Monitor and contribute to the `earthaccess` library for:
- Native subsetting support in `earthaccess.download()`
- Automatic detection and utilization of OGC services
- Chunked download capabilities for large files

**Community involvement**: File feature requests and collaborate with earthaccess maintainers on subsetting capabilities.

### Development Priority

The **HTTP Range + Kerchunk approach** (Option 1) is recommended as the first implementation target because:
- It has the fewest external dependencies
- No waiting for NASA infrastructure changes
- Can be implemented and tested immediately
- Provides graceful fallback to full download if Range requests are unavailable

Contributions and pull requests implementing any of these approaches are welcome!

## Diagnostics

Diagnostic tools and troubleshooting guides are available in the [diagnostics/](diagnostics/) folder:

- **[Install HPC.md](diagnostics/Install%20HPC.md)** - Complete HPC installation guide
- **[HPC Troubleshooting.md](diagnostics/HPC%20Troubleshooting.md)** - Common HPC issues and solutions
- **verify_installation.py** - Check if HDF5/NetCDF libraries are properly configured
- **check_hdf5_version.py** - Display library versions and test file opening
- **debug_hpc_environment.py** - Comprehensive environment diagnostic tool

## References

* Green, R. O., et al. (2023). Earth Surface Mineral Dust Source Investigation (EMIT) L2A Estimated Surface Reflectance and Uncertainty and Masks, Version 1. [Data set]. NASA EOSDIS Land Processes DAAC. [doi:10.5067/EMIT/EMITL2ARFL.001](https://doi.org/10.5067/EMIT/EMITL2ARFL.001)

* Green, R. O., et al. (2024). The Earth Surface Mineral Dust Source Investigation (EMIT) on the International Space Station: In-flight instrument performance and first results. *Remote Sensing of Environment*, 282, 113277. [doi:10.1016/j.rse.2023.113277](https://doi.org/10.1016/j.rse.2023.113277)
