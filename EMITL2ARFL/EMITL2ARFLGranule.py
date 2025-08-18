from os.path import abspath, expanduser
from typing import List
import numpy as np

from rasterio.windows import Window

from .EMITL2ARFLNetCDF import EMITL2ARFLNetCDF
from .EMITL2AMASKNetCDF import EMITL2AMASKNetCDF
from .EMITL2ARFLUNCERTNetCDF import EMITL2ARFLUNCERTNetCDF
from .constants import QUALITY_BANDS
from .read_qmask import read_qmask
from .emit_ortho_raster import emit_ortho_raster
from rasters import Raster, RasterGeometry, RasterGeolocation, RasterGrid
from .GLT import GeometryLookupTable
from .read_netcdf_raster import read_netcdf_raster

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
    
    @property
    def grid(self) -> RasterGrid:
        return self.reflectance_netcdf.grid

    def extract_GLT(
            self, 
            geometry: RasterGeometry = None,
            window: Window = None
            ) -> GeometryLookupTable:
        return self.reflectance_netcdf.extract_GLT(
            geometry=geometry, 
            window=window
        )

    GLT = property(extract_GLT)

    def quality_mask(
            self, 
            window: Window = None,
            quality_bands: List[int] = QUALITY_BANDS) -> np.ndarray:
        qmask: np.ndarray = read_qmask(
            filepath=self.mask_filename,
            window=window,
            quality_bands=quality_bands
        )
        
        return qmask

    def reflectance(
            self, 
            geometry: RasterGeometry = None,
            window: Window = None) -> Raster:
        return read_netcdf_raster(
            filename=self.reflectance_filename,
            variable="reflectance",
            geometry=geometry
        )
    
        # return emit_ortho_raster(
        #     filename=self.reflectance_filename,
        #     layer_name="reflectance",
        #     geometry=geometry,
        #     grid_window=window
        # )

        # processing_subset = geometry is not None or window is not None

        # scene_grid = self.reflectance_netcdf.grid

        # if processing_subset:
        #     if window is None and geometry is not None:
        #         window = scene_grid.window(geometry)
            
        #     if geometry is None:
        #         geometry = scene_grid.subset(window)

        # qmask: np.ndarray = self.quality_mask(window=window)

        # raster: Raster = emit_ortho_raster(
        #     filename=self.reflectance_filename,
        #     layer_name="reflectance",
        #     # qmask=qmask,
        #     geometry=geometry,
        #     window=window
        # )

        # if processing_subset:
        #     raster = raster.to_geometry(geometry)

        # return raster
