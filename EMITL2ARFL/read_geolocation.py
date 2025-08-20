from rasters import RasterGeolocation
from .read_latitude_array import read_latitude_array
from .read_longitude_array import read_longitude_array

from typing import Optional
from rasterio.windows import Window

def read_geolocation(filename: str, window: Optional[Window] = None) -> RasterGeolocation:
    """
    Reads the latitude and longitude arrays from a NetCDF reflectance file and constructs a RasterGeolocation object.

    This function uses helper functions to extract the latitude and longitude arrays from the specified NetCDF file.
    It then creates and returns a RasterGeolocation object, which is used to represent the geospatial coordinates
    associated with the raster data. This is useful for mapping pixel locations to geographic coordinates.

    Args:
        filename (str): Path to the NetCDF reflectance file containing the 'lat' and 'lon' arrays in the 'location' group.
        window (Optional[rasterio.windows.Window], optional):
            If provided, a rasterio Window object specifying the subset (window) of the latitude and longitude arrays to read.
            The window must have attributes row_off, col_off, height, and width. If None, the entire arrays are read.

    Returns:
        RasterGeolocation: An object containing the longitude (x) and latitude (y) arrays for georeferencing.
    """
    # Read the latitude array from the NetCDF file using the helper function
    lat = read_latitude_array(filename, window=window)
    # Read the longitude array from the NetCDF file using the helper function
    lon = read_longitude_array(filename, window=window)
    # Create a RasterGeolocation object using the longitude and latitude arrays
    geolocation = RasterGeolocation(x=lon, y=lat)

    # Return the RasterGeolocation object for use in geospatial operations
    return geolocation
