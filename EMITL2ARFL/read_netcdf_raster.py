import numpy as np
from typing import Optional
from rasterio.windows import Window
import rasters as rt
from rasters import Raster
from .read_netcdf_array import read_netcdf_array
from .read_geolocation import read_geolocation


def read_netcdf_raster(
    filename: str,
    group: str,
    variable: str,
    window: Optional[Window] = None
) -> Raster:
    """
    Read a variable array from a NetCDF file and return as a rasters.Raster object with geolocation.

    Parameters
    ----------
    filename : str
        Path to the NetCDF file to read from.
    group : str
        Name of the group within the NetCDF file containing the variable.
    variable : str
        Name of the variable to read from the group.
    window : Optional[rasterio.windows.Window], default None
        If provided, must be a rasterio Window object specifying the subset (window) to read.

    Returns
    -------
    Raster
        The requested data as a Raster object with geolocation.
    """
    array = read_netcdf_array(
        filename=filename,
        group=group,
        variable=variable,
        window=window
    )
    
    geolocation = read_geolocation(filename)
    raster = Raster(array, geometry=geolocation)

    return raster
