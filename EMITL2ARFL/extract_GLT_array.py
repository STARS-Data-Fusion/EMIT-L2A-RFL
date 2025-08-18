import numpy as np
import xarray as xr
from rasterio.windows import Window

from .constants import *
from .read_netcdf_array import read_netcdf_array

def extract_GLT_array(
        filename: str,
        window: Window = None,
        adjust_indices: bool = False,
        GLT_nodata_value: int = GLT_NODATA_VALUE) -> np.ndarray:
    """
    Extracts the EMIT Geometry Lookup Table (GLT) index pairs from an xarray.Dataset or NetCDF file.

    The GLT is a geospatial mapping array that, for each pixel in the output geographic grid, provides
    the corresponding (row, column) indices in the original satellite swath data. This enables orthorectification:
    the process of transforming raw sensor data into a georeferenced product by mapping each output pixel
    to its correct location in the input data.

    This function:
    - Loads the swath dataset (if a filename is provided).
    - Extracts the 'glt_x' and 'glt_y' arrays, which contain the swath column and row indices for each output pixel.
    - Stacks these into a single array of shape (latitude, longitude, 2), where the last dimension holds (row, column) pairs: (row, column) = (glt_y, glt_x).
    - Replaces any missing values (NaN) with the specified GLT_nodata_value, ensuring all output indices are valid integers.

    Parameters
    ----------
    swath_ds : xr.Dataset | str
        EMIT swath xarray Dataset containing 'glt_x' and 'glt_y' arrays, or a filename to load.
    GLT_nodata_value : int, optional
        Value to use for missing GLT indices (default: GLT_NODATA_VALUE).

    Returns
    -------
    np.ndarray
        Array of shape (latitude, longitude, 2) with GLT index pairs (row, column), dtype=int.
        Missing values are set to GLT_nodata_value.
    """
    # # Step 1: If input is a filename, load the xarray.Dataset
    # if isinstance(swath_dataset, str):
    #     # Local import to avoid circular import
    #     from .emit_xarray import emit_xarray
    #     ds: xr.Dataset = emit_xarray(swath_dataset, ortho=False)
    # else:
    #     ds: xr.Dataset = swath_dataset

    # # Step 2: Extract GLT x (column) and y (row) indices from the dataset
    # # These arrays map each output pixel to its location in the original swath
    # GLT_x: np.ndarray = ds["glt_x"].data  # swath column indices
    # GLT_y: np.ndarray = ds["glt_y"].data  # swath row indices
    GLT_x = read_netcdf_array(
        filename=filename,
        variable="glt_x", 
        group="location",
        window=window
    )

    GLT_y = read_netcdf_array(
        filename=filename,
        variable="glt_y", 
        group="location",
        window=window
    )

    # Step 3: Stack row and column indices into a single array of shape (latitude, longitude, 2)
    # The last dimension holds (row, column) pairs for each output pixel: (row, column) = (glt_y, glt_x)
    GLT_array: np.ndarray = np.nan_to_num(
        np.stack([GLT_y, GLT_x], axis=-1),
        nan=GLT_nodata_value
    ).astype(int)

    # Step 4: Optionally adjust indices if requested
    if adjust_indices and window is not None:
        # window.col_off and window.row_off are the offsets for the subset
        # Only adjust non-nodata values
        mask = (GLT_array[..., 0] != GLT_nodata_value) & (GLT_array[..., 1] != GLT_nodata_value)
        GLT_array[..., 0][mask] -= int(window.row_off)
        GLT_array[..., 1][mask] -= int(window.col_off)

    # Step 5: Return the GLT array, ready for use in geospatial orthorectification
    return GLT_array
