import numpy as np
from affine import Affine
import rasters as rt

from .constants import *
from .emit_xarray import emit_xarray

def emit_ortho_raster(
        filepath: str, 
        layer_name: str,
        qmask: np.ndarray = None, 
        unpacked_bmask: np.ndarray = None, 
        fill_value: int = FILL_VALUE,
        engine: str = ENGINE) -> rt.Raster:
    """
    Load an EMIT NetCDF data layer and orthorectify it as `rasters.Raster` object.

    Parameters:
    filepath: a filepath to an EMIT netCDF file
    layer_name: the name of the data layer to be orthorectified
    qmask: a numpy array output from the quality_mask function used to mask pixels based on quality flags selected in that function. Any non-orthorectified array with the proper crosstrack and downtrack dimensions can also be used.
    unpacked_bmask: a numpy array from  the band_mask function that can be used to mask band-specific pixels that have been interpolated.

    Returns:
    raster.Raster object containing the orthorectified EMIT data layer
    """
    ortho_ds = emit_xarray(
        filepath=filepath,
        ortho=True,
        qmask=qmask,
        unpacked_bmask=unpacked_bmask,
        fill_value=fill_value,
        engine=engine
    )

    latitude_length, longitude_length, bands = ortho_ds.reflectance.shape
    affine = Affine.from_gdal(*ortho_ds.geotransform)
    grid = rt.RasterGrid.from_affine(affine, longitude_length, latitude_length)
    raster = rt.MultiRaster(np.transpose(np.array(ortho_ds[layer_name]), (2, 0, 1)), geometry=grid)

    return raster
