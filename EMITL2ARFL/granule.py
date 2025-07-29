from glob import glob
from os.path import join, abspath, dirname, expanduser
from typing import List

import numpy as np




from rasters import Raster, RasterGeometry, RasterGeolocation
from .read_latitude import read_latitude
from .read_longitude import read_longitude
from .read_geolocation import read_geolocation

from .constants import *
from .emit_ortho_raster import emit_ortho_raster
from .quality_mask import quality_mask
from .extract_GLT_array import extract_GLT_array
from .extract_GLT import extract_GLT
from .GLT import GLT

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
        return read_latitude(self.reflectance_filename)
    
    @property
    def lon(self) -> np.ndarray:
        return read_longitude(self.reflectance_filename)
    
    @property
    def geolocation(self) -> RasterGeolocation:
        return read_geolocation(self.reflectance_filename)

    def GLT(self) -> GLT:
        return extract_GLT(self.reflectance_filename)

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
    