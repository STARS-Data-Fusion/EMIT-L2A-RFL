import numpy as np
from typing import Optional
from rasterio.windows import Window
import rasters as rt
from rasters import Raster, RasterGeometry
from .read_netcdf_array import read_netcdf_array
from .read_geolocation import read_geolocation

def read_netcdf_raster(
    filename: str,
    variable: str,
    group: Optional[str] = None,
    geometry: Optional[RasterGeometry] = None,
    window: Optional[Window] = None,
    resampling: str = "nearest"
) -> Raster:
    """
    Read a variable array from a NetCDF file and return as a rasters.Raster object with geolocation, supporting spatial subsetting.

    Parameters
    ----------
    filename : str
        Path to the NetCDF file to read from.
    variable : str
        Name of the variable to read from the group or root.
    group : Optional[str], default None
        Name of the group within the NetCDF file containing the variable. If None, reads from the root of the NetCDF data structure (global scope).
    geometry : Optional[RasterGeometry], default None
        If provided, specifies a spatial geometry to subset the data. Ignored if window is given.
    window : Optional[rasterio.windows.Window], default None
        If provided, must be a rasterio Window object specifying the spatial subset (window) to read. Takes precedence over geometry.
    resampling : str, default "nearest"
        Resampling method to use if geometry is provided and a reprojection or resampling is needed.

    Returns
    -------
    Raster
        The requested data as a Raster object with geolocation, optionally spatially subsetted.
    """
    # If a window is not provided but a geometry is, compute the window from the geometry
    if window is None and geometry is not None:
        # read the scene geolocation
        geolocation = read_geolocation(filename=filename)

        # calculate the indices window that covers the target geometry
        window = geolocation.window(geometry)

    # Read the data array from the NetCDF file, using the window if provided
    array = read_netcdf_array(
        filename=filename,
        variable=variable,
        group=group,
        window=window
    )

    # Read the geolocation, using the same window for spatial alignment
    geolocation = read_geolocation(filename, window=window)

    # Wrap the data array and geolocation in a Raster object
    raster = Raster(array, geometry=geolocation)

    # If a geometry is provided, reproject or resample the raster to match the geometry
    if geometry is not None:
        raster = raster.to_geometry(geometry, resampling=resampling)

    return raster
