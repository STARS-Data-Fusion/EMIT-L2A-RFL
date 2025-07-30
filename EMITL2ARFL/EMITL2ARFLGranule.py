from os.path import abspath, expanduser
from typing import List
import numpy as np
from .EMITL2ARFLNetCDF import EMITL2ARFLNetCDF
from .EMITL2AMASKNetCDF import EMITL2AMASKNetCDF
from .EMITL2ARFLUNCERTNetCDF import EMITL2ARFLUNCERTNetCDF
from .constants import QUALITY_BANDS
from .quality_mask import quality_mask
from .emit_ortho_raster import emit_ortho_raster
from rasters import Raster, RasterGeometry, RasterGeolocation
from .GLT import GLT

class EMITL2ARFLGranule:
    def __init__(self, reflectance_filename: str, mask_filename: str, uncertainty_filename: str) -> None:
        self.reflectance_filename: str = abspath(expanduser(reflectance_filename))
        self.mask_filename: str = abspath(expanduser(mask_filename))
        self.uncertainty_filename: str = abspath(expanduser(uncertainty_filename))

    def __repr__(self) -> str:
        return (f"EMITL2ARFL(reflectance_filename=\"{self.reflectance_filename}\", "
                f"mask_filename=\"{self.mask_filename}\", "
                f"uncertainty_filename=\"{self.uncertainty_filename}\")")

    @property
    def reflectance_netcdf(self) -> EMITL2ARFLNetCDF:
        return EMITL2ARFLNetCDF(self.reflectance_filename)

    @property
    def mask_netcdf(self) -> EMITL2AMASKNetCDF:
        return EMITL2AMASKNetCDF(self.mask_filename)

    @property
    def uncertainty_netcdf(self) -> EMITL2ARFLUNCERTNetCDF:
        return EMITL2ARFLUNCERTNetCDF(self.uncertainty_filename)

    @property
    def lat(self) -> np.ndarray:
        return self.reflectance_netcdf.lat

    @property
    def lon(self) -> np.ndarray:
        return self.reflectance_netcdf.lon

    @property
    def geolocation(self) -> RasterGeolocation:
        return self.reflectance_netcdf.geolocation

    def extract_GLT(self) -> GLT:
        return self.reflectance_netcdf.extract_GLT()

    GLT = property(extract_GLT)

    def quality_mask(self, quality_bands: List[int] = QUALITY_BANDS) -> np.ndarray:
        qmask: np.ndarray = quality_mask(
            filepath=self.mask_filename,
            quality_bands=quality_bands
        )
        
        return qmask

    def reflectance(self, geometry: RasterGeometry = None) -> Raster:
        qmask: np.ndarray = self.quality_mask()

        raster: Raster = emit_ortho_raster(
            filepath=self.reflectance_filename,
            layer_name="reflectance",
            qmask=qmask
        )

        if geometry is not None:
            raster = raster.to_geometry(geometry)

        return raster
