import numpy as np
import xarray as xr

from .constants import *

def extract_GLT(swath_ds: xr.Dataset, GLT_nodata_value: int = GLT_NODATA_VALUE) -> np.ndarray:
    """
    Extracts EMIT geometry look-up table (GLT) as a numpy array of size (latitude, longitude, 2) from EMIT swath xarray Dataset containing `glt_x` and `glt_y` arrays.
    If a filename is provided instead of an xarray.Dataset, the function will load the dataset using xr.open_dataset.
    """
    # If input is a string (filename), load the dataset
    if isinstance(swath_ds, str):
        from .emit_xarray import emit_xarray
        ds = emit_xarray(swath_ds, ortho=False)
    else:
        ds = swath_ds

    # extract GLT x indices
    GLT_x = ds["glt_x"].data
    # extract GLT y indices
    GLT_y = ds["glt_y"].data
    # stack GLT index pairs
    GLT_array = np.nan_to_num(np.stack([GLT_x, GLT_y], axis=-1), nan=GLT_nodata_value).astype(int)

    return GLT_array
