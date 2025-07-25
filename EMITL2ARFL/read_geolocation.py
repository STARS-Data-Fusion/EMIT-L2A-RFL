from rasters import RasterGeolocation
from .read_latitude import read_latitude
from .read_longitude import read_longitude

def read_geolocation(filename: str) -> RasterGeolocation:
    """
    Reads the latitude and longitude arrays from a NetCDF reflectance file and constructs a RasterGeolocation object.

    This function uses helper functions to extract the latitude and longitude arrays from the specified NetCDF file.
    It then creates and returns a RasterGeolocation object, which is used to represent the geospatial coordinates
    associated with the raster data. This is useful for mapping pixel locations to geographic coordinates.

    Args:
        filename (str): Path to the NetCDF reflectance file containing the 'lat' and 'lon' arrays in the 'location' group.

    Returns:
        RasterGeolocation: An object containing the longitude (x) and latitude (y) arrays for georeferencing.
    """
    # Read the latitude array from the NetCDF file using the helper function
    lat = read_latitude(filename)
    # Read the longitude array from the NetCDF file using the helper function
    lon = read_longitude(filename)
    # Create a RasterGeolocation object using the longitude and latitude arrays
    geolocation = RasterGeolocation(x=lon, y=lat)

    # Return the RasterGeolocation object for use in geospatial operations
    return geolocation
