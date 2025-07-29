import netCDF4
import numpy as np

def read_longitude(reflectance_filename: str) -> np.ndarray:
    """Read the `lon` array from the `location` group in the reflectance NetCDF file."""
    with netCDF4.Dataset(reflectance_filename, "r") as ds:
        lon = ds.groups["location"].variables["lon"][:]
    return lon
