from os.path import abspath, expanduser
from typing import List, Union
import numpy as np

from rasterio.windows import Window

import rasters as rt
from rasters import Raster, RasterGeometry, RasterGeolocation, RasterGrid

from .EMITL2ARFLNetCDF import EMITL2ARFLNetCDF
from .EMITL2AMASKNetCDF import EMITL2AMASKNetCDF
from .EMITL2ARFLUNCERTNetCDF import EMITL2ARFLUNCERTNetCDF
from .constants import QUALITY_BANDS
from .read_qmask import read_qmask
from .emit_ortho_raster import emit_ortho_raster

from .GLT import GeometryLookupTable
from .read_netcdf_raster import read_netcdf_raster
from .read_geolocation import read_geolocation

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
            swath_window: Window = None,
            geometry: RasterGeometry = None,
            quality_bands: List[int] = QUALITY_BANDS) -> Union[Raster, np.ndarray]:
        granule_geolocation = read_geolocation(filename=self.reflectance_filename)
        
        if swath_window is None and geometry is not None:
            swath_window = granule_geolocation.window(geometry)
            subset_geolocation = granule_geolocation[swath_window]

        if geometry is None and swath_window is not None:
            subset_geolocation = granule_geolocation[swath_window]

        if geometry is None and swath_window is None:
            subset_geolocation = granule_geolocation

        qmask: np.ndarray = read_qmask(
            filename=self.mask_filename,
            window=swath_window,
            quality_bands=quality_bands
        )
        
        if geometry is not None:
            qmask = Raster(qmask, geometry=subset_geolocation)
            qmask = qmask.to_geometry(geometry)
            
        qmask = rt.where(np.isnan(qmask), 0, qmask)
        qmask = qmask.astype(bool)
        
        return qmask

    def reflectance(
            self, 
            geometry: RasterGeometry = None,
            swath_window: Window = None,
            qmask: np.ndarray = None) -> Raster:
            # If a window is not provided but a geometry is, compute the window from the geometry
        if swath_window is None and geometry is not None:
            # read the scene geolocation
            geolocation = read_geolocation(filename=self.reflectance_filename)

        # calculate the indices window that covers the target geometry
        swath_window = geolocation.window(geometry)

        qmask = self.quality_mask(swath_window=swath_window)

        result = read_netcdf_raster(
            filename=self.reflectance_filename,
            variable="reflectance",
            geometry=geometry,
            swath_window=swath_window,
            qmask=qmask
        )

        return result
    
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
