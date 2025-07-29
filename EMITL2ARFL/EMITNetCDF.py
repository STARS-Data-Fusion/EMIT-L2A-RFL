from glob import glob
from os.path import join, abspath, dirname, expanduser
from typing import List
import numpy as np
from rasters import Raster, RasterGeometry, RasterGeolocation
from .read_netcdf_raster import read_netcdf_raster
from .read_netcdf_array import read_netcdf_array
from .read_latitude_array import read_latitude_array
from .read_longitude_array import read_longitude_array
from .read_geolocation import read_geolocation
from .constants import *
from .emit_ortho_raster import emit_ortho_raster
from .quality_mask import quality_mask
from .extract_GLT_array import extract_GLT_array
from .extract_GLT import extract_GLT
from .GLT import GLT
from .show_netcdf_tree import show_netcdf_tree

class EMITNetCDF:
    def __init__(self, filename: str) -> None:
        self.filename: str = abspath(expanduser(filename))

    def __repr__(self) -> str:
        return f"EMITNetCDF(filename=\"{self.filename}\")"

    def show_tree(self, indent: int = 0) -> str:
        return show_netcdf_tree(self.filename, indent)
    
    tree = property(show_tree)

    @property
    def lat(self) -> np.ndarray:
        return read_latitude_array(self.filename)

    @property
    def lon(self) -> np.ndarray:
        return read_longitude_array(self.filename)

    @property
    def geolocation(self) -> RasterGeolocation:
        return read_geolocation(self.filename)
    
    def extract_GLT(self) -> GLT:
        return extract_GLT(self.filename)

    GLT = property(extract_GLT)

    def read(self, group: str, variable: str):
        return read_netcdf_raster(
            filename=self.filename,
            group=group,
            variable=variable
        )

    def read_elevation(self) -> Raster:
        return read_netcdf_raster(
            filename=self.filename,
            group="location",
            variable="elev"
        )
    
    elevation = property(read_elevation)

    def read_array(self, group: str, variable: str) -> np.ndarray:
        return read_netcdf_array(self.filename, group, variable)
