
from typing import Optional
import numpy as np
from rasterio.windows import Window
from .read_netcdf_array import read_netcdf_array

def read_elevation(filename: str, window: Optional[Window] = None) -> np.ndarray:
    """
    Read the `elev` array from the `location` group in the reflectance NetCDF file.

    Parameters
    ----------
    filename : str
        Path to the NetCDF file.
    window : Optional[Window], default None
        If provided, only the subset defined by the window will be read.

    Returns
    -------
    np.ndarray
        The elevation array or its subset.
    """
from typing import Optional
import numpy as np
from rasterio.windows import Window
from .read_netcdf_array import read_netcdf_array

def read_elevation(filename: str, window: Optional[Window] = None) -> np.ndarray:
    """
    Read the `elev` array from the `location` group in the reflectance NetCDF file.

    Parameters
    ----------
    filename : str
        Path to the NetCDF file.
    window : Optional[Window], default None
        If provided, only the subset defined by the window will be read.

    Returns
    -------
    np.ndarray
        The elevation array or its subset.
    """
    return read_netcdf_array(
        filename=filename,
        variable="elev",
        group="location",
        window=window
    )
