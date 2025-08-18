from typing import List
import numpy as np
import xarray as xr
from rasterio.windows import Window

from .constants import *

def read_qmask(
        filepath: str, 
        window: Window = None,
        quality_bands: List[str] = QUALITY_BANDS, 
        engine: str = ENGINE) -> np.ndarray:
    """
    This function builds a single layer mask to apply based on the bands selected from an EMIT L2A Mask file.

    Parameters:
    filepath: an EMIT L2A Mask netCDF file.
    quality_bands: a list of bands (quality flags only) from the mask file that should be used in creation of  mask.

    Returns:
    qmask: a numpy array that can be used with the emit_xarray function to apply a quality mask.
    """
    # Open Dataset
    mask_ds = xr.open_dataset(filepath, engine=ENGINE)

    # Open Sensor band Group
    mask_parameters_ds = xr.open_dataset(
        filepath,
        engine=engine,
        group="sensor_band_parameters"
    )

    # Print Flags used
    flags_used = mask_parameters_ds["mask_bands"].data[quality_bands]

    # Check for data bands and build mask
    if any(x in quality_bands for x in [5, 6]):
        err_str = f"Selected flags include a data band (5 or 6) not just flag bands"
        raise AttributeError(err_str)
    else:
        if window is not None:
            # window: Window(col_off, row_off, width, height)
            row_off, col_off = int(window.row_off), int(window.col_off)
            height, width = int(window.height), int(window.width)
            qmask = np.sum(
                mask_ds["mask"][row_off:row_off+height, col_off:col_off+width, quality_bands].values,
                axis=-1
            )
        else:
            qmask = np.sum(mask_ds["mask"][:, :, quality_bands].values, axis=-1)
        qmask[qmask > 1] = 1

    return qmask
