
from typing import Optional
import netCDF4
import numpy as np
from rasterio.windows import Window

def read_netcdf_array(
    filename: str,
    variable: str,
    group: Optional[str] = None,
    window: Optional[Window] = None
) -> np.ndarray:
    """
    Read a variable array from a specified group or the root in a NetCDF file.

    Parameters
    ----------
    filename : str
        Path to the NetCDF file to read from.
    group : Optional[str], default None
        Name of the group within the NetCDF file containing the variable. If None, reads from the root of the NetCDF data structure (global scope).
    variable : str
        Name of the variable to read from the group.
    window : Optional[rasterio.windows.Window], default None
        If provided, must be a rasterio Window object specifying the subset (window) to read.
        The window must have attributes row_off, col_off, height, and width.

    Returns
    -------
    np.ndarray
        The requested array, or the specified subset if window is provided.

    Raises
    ------
    KeyError
        If the group or variable does not exist in the NetCDF file.
    AttributeError
        If the window object does not have the required attributes.
    """
    # Open the NetCDF file for reading
    with netCDF4.Dataset(filename, "r") as ds:
        # Access the specified group and variable, or root if group is None
        if group is None:
            var = ds.variables[variable]
        else:
            var = ds.groups[group].variables[variable]
        if window is not None:
            # Ensure the window object has the required attributes
            # (row_off, col_off, height, width)
            row_off = int(window.row_off)
            col_off = int(window.col_off)
            height = int(window.height)
            width = int(window.width)
            # Read the specified window (subset) from the variable
            arr = var[row_off:row_off+height, col_off:col_off+width]
        else:
            # Read the entire variable array
            arr = var[:]
    # If the array has 3 dimensions, transpose from (rows, cols, bands) to (bands, rows, cols)
    if arr.ndim == 3:
        arr = np.transpose(arr, (2, 0, 1))
    return arr
