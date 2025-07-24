from glob import glob
from os.path import join, abspath, dirname, expanduser
from typing import List

import numpy as np

import netCDF4

from rasters import Raster, RasterGeometry, RasterGeolocation

from .constants import *
from .emit_ortho_raster import emit_ortho_raster
from .quality_mask import quality_mask
from .extract_GLT_array import extract_GLT_array

class EMITL2ARFL:
    def __init__(self, directory: str):
        self.directory = directory
    
    def __repr__(self) -> str:
        return f"EMITL2ARFL(directory=\"{self.directory}\")"

    @property
    def directory_absolute(self) -> str:
        return abspath(expanduser(self.directory))

    @property
    def files(self):
        return glob(join(self.directory_absolute, "*.nc"))
    
    @property
    def reflectance_filename(self) -> str:
        return glob(join(self.directory_absolute, "*_RFL_*.nc"))[0]
    
    @property
    def mask_filename(self) -> str:
        return glob(join(self.directory_absolute, "*_MASK_*.nc"))[0]
    
    @property
    def uncertainty_filename(self) -> str:
        return glob(join(self.directory_absolute, "*_RFLUNCERT_*.nc"))[0]
    
    @property
    def lat(self) -> np.ndarray:
        # read the `lat` array from the `location` group in the reflectance NetCDF file    
        with netCDF4.Dataset(self.reflectance_filename, "r") as ds:
            lat = ds.groups["location"].variables["lat"][:]

        return lat
    
    @property
    def lon(self) -> np.ndarray:
        # read the `lon` array from the `location` group in the reflectance NetCDF file    
        with netCDF4.Dataset(self.reflectance_filename, "r") as ds:
            lon = ds.groups["location"].variables["lon"][:]

        return lon
    
    @property
    def geolocation(self) -> RasterGeolocation:
        return RasterGeolocation(
            x=self.lon,
            y=self.lat
        )

    def GLT(self) -> Raster:
        return Raster(extract_GLT_array(swath_ds=self.reflectance_filename), geometry=self.geolocation)

    def quality_mask(self, quality_bands: List[int] = QUALITY_BANDS) -> np.ndarray:
        qmask = quality_mask(
            filepath=self.mask_filename,
            quality_bands=quality_bands
        )

        return qmask

    def reflectance(self, geometry: RasterGeometry = None) -> Raster:
        qmask = self.quality_mask()

        raster = emit_ortho_raster(
            filepath=self.reflectance_filename,
            layer_name="reflectance",
            qmask=qmask
        )

        if geometry is not None:
            raster = raster.to_geometry(geometry)
        
        return raster
    