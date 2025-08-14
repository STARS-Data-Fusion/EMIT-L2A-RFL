import xarray as xr
from rasterio.windows import Window

from .extract_grid import extract_grid

from .extract_GLT_array import extract_GLT_array
from .GLT import GeometryLookupTable
from .emit_xarray import emit_xarray
from .constants import GLT_NODATA_VALUE

def extract_GLT(
        filename: str,
        window: Window = None,
        GLT_nodata_value: int = GLT_NODATA_VALUE):
    """
    Wrapper function to extract a GLT object from a NetCDF file.
    """
    GLT_array = extract_GLT_array(filename, GLT_nodata_value=GLT_nodata_value, window=window)
    grid = extract_grid(filename, window=window)

    # a GLT is a grid of swath indices that needs to be georeferenced as a grid

    return GeometryLookupTable(GLT_array=GLT_array, geometry=grid)
