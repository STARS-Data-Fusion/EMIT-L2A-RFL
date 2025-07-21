from glob import glob
from os.path import join, abspath, dirname, expanduser

from rasters import Raster, RasterGeometry

from .emit_ortho_raster import emit_ortho_raster
from .quality_mask import quality_mask

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
    
    def quality_mask(self, geometry: RasterGeometry) -> Raster:
        raster = quality_mask(
            filepath=self.mask_filename,
            quality_bands=[0, 1, 2, 3, 4]
        )

        if geometry is not None:
            raster = raster.to_geometry(geometry)
        
        return raster

    def reflectance(self, geometry: RasterGeometry) -> Raster:
        raster = emit_ortho_raster(
            filepath=self.reflectance_filename,
            layer_name="reflectance"
        )

        if geometry is not None:
            raster = raster.to_geometry(geometry)
        
        return raster
    