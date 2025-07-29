from os.path import abspath, expanduser
import numpy as np
from rasters import Raster, RasterGeolocation
from .read_netcdf_raster import read_netcdf_raster
from .read_netcdf_array import read_netcdf_array
from .read_latitude_array import read_latitude_array
from .read_longitude_array import read_longitude_array
from .read_geolocation import read_geolocation
from .extract_GLT import extract_GLT
from .GLT import GLT
from .show_netcdf_tree import show_netcdf_tree

class EMITNetCDF:
    """
    Provides convenient access to EMIT NetCDF files, including reading arrays, rasters, geolocation, and GLT.
    """
    def __init__(self, filename: str) -> None:
        """
        Initialize an EMITNetCDF object for a given NetCDF file.

        Args:
            filename (str): Path to the NetCDF file.
        """
        self.filename: str = abspath(expanduser(filename))

    def __repr__(self) -> str:
        """
        Return a string representation of the EMITNetCDF object.
        """
        return f"EMITNetCDF(filename=\"{self.filename}\")"

    def show_tree(self, indent: int = 0) -> str:
        """
        Show the NetCDF file structure as a tree.

        Args:
            indent (int, optional): Indentation level for the tree. Defaults to 0.

        Returns:
            str: String representation of the NetCDF file tree.
        """
        return show_netcdf_tree(self.filename, indent)
    
    tree = property(show_tree)

    @property
    def lat(self) -> np.ndarray:
        """
        Read the latitude array from the NetCDF file.

        Returns:
            np.ndarray: Latitude values.
        """
        return read_latitude_array(self.filename)

    @property
    def lon(self) -> np.ndarray:
        """
        Read the longitude array from the NetCDF file.

        Returns:
            np.ndarray: Longitude values.
        """
        return read_longitude_array(self.filename)

    @property
    def geolocation(self) -> RasterGeolocation:
        """
        Read the geolocation information from the NetCDF file.

        Returns:
            RasterGeolocation: Geolocation object containing spatial info.
        """
        return read_geolocation(self.filename)
    
    def extract_GLT(self) -> GLT:
        """
        Extract the Geometry Lookup Table (GLT) from the NetCDF file.

        Returns:
            GLT: Geometry Lookup Table object.
        """
        return extract_GLT(self.filename)

    GLT = property(extract_GLT)

    def read(self, group: str, variable: str) -> Raster:
        """
        Read a variable as a Raster object from a specified group in the NetCDF file.

        Args:
            group (str): Name of the group in the NetCDF file.
            variable (str): Name of the variable to read.

        Returns:
            Raster: Raster object of the requested variable with geolocation.
        """
        return read_netcdf_raster(
            filename=self.filename,
            group=group,
            variable=variable
        )

    def read_elevation(self) -> Raster:
        """
        Read the elevation raster from the NetCDF file.

        Returns:
            Raster: Elevation raster object.
        """
        return read_netcdf_raster(
            filename=self.filename,
            group="location",
            variable="elev"
        )
    
    elevation = property(read_elevation)

    def read_array(self, group: str, variable: str) -> np.ndarray:
        """
        Read a variable array from a specified group in the NetCDF file.

        Args:
            group (str): Name of the group in the NetCDF file.
            variable (str): Name of the variable to read.

        Returns:
            np.ndarray: Array of the requested variable.
        """
        return read_netcdf_array(
            filename=self.filename,
            group=group,
            variable=variable
        )
