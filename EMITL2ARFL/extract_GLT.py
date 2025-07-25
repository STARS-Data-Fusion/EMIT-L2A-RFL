import xarray as xr

from .read_geolocation import read_geolocation
from .extract_GLT_array import extract_GLT_array
from .GLT import GLT
from .emit_xarray import emit_xarray
from .constants import GLT_NODATA_VALUE

def extract_GLT(swath_dataset: xr.Dataset, GLT_nodata_value: int = GLT_NODATA_VALUE):
    """
    Wrapper function to extract a GLT object from a dataset or filename.
    """
    # If input is a filename, load the xarray.Dataset
    if isinstance(swath_dataset, str):
        ds: xr.Dataset = emit_xarray(swath_dataset, ortho=False)
    else:
        ds: xr.Dataset = swath_dataset

    GLT_array = extract_GLT_array(ds, GLT_nodata_value)
    geolocation = read_geolocation(swath_dataset)

    return GLT(GLT_array=GLT_array, geolocation=geolocation)
