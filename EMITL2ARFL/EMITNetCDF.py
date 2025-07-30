from typing import List, Dict, Optional
from rasterio.windows import Window

from os.path import abspath, expanduser
import numpy as np

import netCDF4

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

    def read(
        self,
        variable: str,
        group: str = None,
        geometry: Optional[RasterGeolocation] = None,
        window: Optional[Window] = None,
        resampling: str = "nearest"
    ) -> Raster:
        """
        Read a variable as a Raster object from a specified group in the NetCDF file, supporting spatial subsetting.

        Args:
            variable (str): Name of the variable to read.
            group (str, optional): Name of the group in the NetCDF file. If None (default),
                the variable is read from the root of the NetCDF data structure (global scope).
            geometry (RasterGeolocation, optional): Spatial geometry for subsetting. Ignored if window is provided.
            window (rasterio.windows.Window, optional): Spatial window for subsetting. If provided, only the subset is read. Takes precedence over geometry.
            resampling (str, optional): Resampling method if geometry is provided. Defaults to "nearest".

        Returns:
            Raster: Raster object of the requested variable with geolocation, optionally spatially subsetted.
        """
        return read_netcdf_raster(
            filename=self.filename,
            variable=variable,
            group=group,
            geometry=geometry,
            window=window,
            resampling=resampling
        )

    def read_elevation(
        self,
        geometry: Optional[RasterGeolocation] = None,
        window: Optional[Window] = None,
        resampling: str = "nearest"
    ) -> Raster:
        """
        Read the elevation raster from the NetCDF file, with optional spatial subsetting.

        Args:
            geometry (RasterGeolocation, optional): Spatial geometry for subsetting. Ignored if window is provided.
            window (rasterio.windows.Window, optional): Spatial window for subsetting. If provided, only the subset is read. Takes precedence over geometry.
            resampling (str, optional): Resampling method if geometry is provided. Defaults to "nearest".

        Returns:
            Raster: Elevation raster object, optionally spatially subsetted.
        """
        return read_netcdf_raster(
            filename=self.filename,
            variable="elev",
            group="location",
            geometry=geometry,
            window=window,
            resampling=resampling
        )
    
    elevation = property(read_elevation)

    def read_array(self, variable: str, group: str) -> np.ndarray:
        """
        Read a variable array from a specified group in the NetCDF file.

        Args:
            variable (str): Name of the variable to read.
            group (str): Name of the group in the NetCDF file.

        Returns:
            np.ndarray: Array of the requested variable.
        """
        return read_netcdf_array(
            filename=self.filename,
            variable=variable,
            group=group
        )

    @property
    def metadata(self) -> Dict[str, str]:
        """
        Return a dictionary of global NetCDF file attributes (ncattrs).

        Returns:
            dict: Dictionary of attribute names and their values.
        """
        with netCDF4.Dataset(self.filename, 'r') as ds:
            return {attr: ds.getncattr(attr) for attr in ds.ncattrs()}
    
    @property
    def groups(self) -> List[str]:
        """
        List all groups in the NetCDF file.

        Returns:
            List[str]: List of group names.
        """
        with netCDF4.Dataset(self.filename, 'r') as ds:
            return list(ds.groups.keys())
        
    def variables(self, group: str = None) -> List[str]:
        """
        List all variables in the NetCDF file or a specific group.

        Args:
            group (str, optional): Name of the group to list variables from. Defaults to None (global scope).

        Returns:
            List[str]: List of variable names.
        """
        with netCDF4.Dataset(self.filename, 'r') as ds:
            if group is None:
                return list(ds.variables.keys())
            else:
                return list(ds.groups[group].variables.keys())