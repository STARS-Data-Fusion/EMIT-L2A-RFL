import numpy as np
import xarray as xr

def quality_mask(filepath, quality_bands):
    """
    This function builds a single layer mask to apply based on the bands selected from an EMIT L2A Mask file.

    Parameters:
    filepath: an EMIT L2A Mask netCDF file.
    quality_bands: a list of bands (quality flags only) from the mask file that should be used in creation of  mask.

    Returns:
    qmask: a numpy array that can be used with the emit_xarray function to apply a quality mask.
    """
    # Open Dataset
    mask_ds = xr.open_dataset(filepath, engine="h5netcdf")
    # Open Sensor band Group
    mask_parameters_ds = xr.open_dataset(
        filepath, engine="h5netcdf", group="sensor_band_parameters"
    )
    # Print Flags used
    flags_used = mask_parameters_ds["mask_bands"].data[quality_bands]
    print(f"Flags used: {flags_used}")
    # Check for data bands and build mask
    if any(x in quality_bands for x in [5, 6]):
        err_str = f"Selected flags include a data band (5 or 6) not just flag bands"
        raise AttributeError(err_str)
    else:
        qmask = np.sum(mask_ds["mask"][:, :, quality_bands].values, axis=-1)
        qmask[qmask > 1] = 1
    return qmask