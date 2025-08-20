import numpy as np
import xarray as xr
from .apply_qmask import apply_qmask

def apply_qmask_xr(dataset: xr.Dataset, qmask: np.ndarray, fill_value: int = -9999) -> xr.Dataset:
    """
    Apply a quality mask to all data variables in an xarray.Dataset.

    Parameters:
    dataset: xarray.Dataset to mask
    qmask: numpy array with the same spatial shape as the dataset, where 1 indicates masked pixels
    fill_value: value to assign to masked pixels (default: -9999)

    Returns:
    Masked xarray.Dataset
    """
    masked_ds = dataset.copy(deep=True)
    for var in masked_ds.data_vars:
        if qmask is not None:
            arr = masked_ds[var].data
            masked_arr = apply_qmask(arr, qmask, fill_value)
            masked_ds[var].data = masked_arr
    return masked_ds
