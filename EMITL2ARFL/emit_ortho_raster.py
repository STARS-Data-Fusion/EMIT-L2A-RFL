import numpy as np
from affine import Affine
import rasters as rt
from rasters import RasterGeometry
from rasterio.windows import Window

from .constants import *
from .emit_xarray import emit_xarray
from .extract_grid import extract_grid
from .GLT import GeometryLookupTable
from .extract_GLT import extract_GLT

def emit_ortho_raster(
        filename: str, 
        layer_name: str,
        geometry: RasterGeometry = None,
        window: Window = None,
        GLT: GeometryLookupTable = None,
        qmask: np.ndarray = None, 
        unpacked_bmask: np.ndarray = None, 
        fill_value: int = FILL_VALUE,
        engine: str = ENGINE,
        GLT_nodata_value: int = GLT_NODATA_VALUE) -> rt.Raster:
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
    processing_subset = geometry is not None or window is not None

    scene_grid = extract_grid(filename=filename)

    if processing_subset:
        if processing_subset:
            if window is None and geometry is not None:
                window = scene_grid.window(geometry)
            
            if geometry is None:
                geometry = scene_grid.subset(window)
    else:
        geometry = scene_grid

    # FIXME need to generate swath window from grid window based on GLT
    if GLT is None:
        GLT = extract_GLT(
            filename=filename,
            window=window,
            GLT_nodata_value=GLT_nodata_value
        )

    swath_window = GLT.swath_window

    ortho_ds = emit_xarray(
        filepath=filename,
        ortho=True,
        swath_window=swath_window,
        qmask=qmask,
        unpacked_bmask=unpacked_bmask,
        fill_value=fill_value,
        engine=engine
    )

    # latitude_length, longitude_length, bands = ortho_ds.reflectance.shape
    # affine = Affine.from_gdal(*ortho_ds.geotransform)
    # grid = rt.RasterGrid.from_affine(affine, longitude_length, latitude_length)
    raster = rt.MultiRaster(np.transpose(np.array(ortho_ds[layer_name]), (2, 0, 1)), geometry=geometry)

    return raster
