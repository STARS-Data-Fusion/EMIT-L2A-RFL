import netCDF4
import numpy as np

def read_latitude(reflectance_filename: str) -> np.ndarray:
    """Read the `lat` array from the `location` group in the reflectance NetCDF file."""
    with netCDF4.Dataset(reflectance_filename, "r") as ds:
        lat = ds.groups["location"].variables["lat"][:]
    return lat
