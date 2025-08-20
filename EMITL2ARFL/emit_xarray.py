import netCDF4 as nc
import os
from spectral.io import envi
import numpy as np
import math
from skimage import io
import pandas as pd
import geopandas as gpd
import xarray as xr
import rasterio as rio
import rioxarray as rxr
import s3fs
from rioxarray.merge import merge_arrays
from fsspec.implementations.http import HTTPFile
from rasterio.windows import Window
from os.path import join, expanduser

from .constants import *
from .ortho_xr import ortho_xr
from .extract_GLT_array import extract_GLT_array

def emit_xarray(
    filename: str, 
    ortho: bool = False, 
    swath_window: Window = None,
    grid_window: Window = None,
    GLT_array: np.ndarray = None,
    adjusted_GLT: np.ndarray = None,
    qmask: np.ndarray = None, 
    unpacked_bmask: np.ndarray = None, 
    fill_value: int = FILL_VALUE,
    engine: str = ENGINE
) -> xr.Dataset:
    """
    Load an EMIT NetCDF dataset as an xarray.Dataset.

    Parameters:
    filepath: a filepath to an EMIT netCDF file
    ortho: hyperspectral cube is left is swath spatial dimensions if `ortho` is False and orthorectified to a grid using the included GLT if True
    qmask: a numpy array output from the quality_mask function used to mask pixels based on quality flags selected in that function. Any non-orthorectified array with the proper crosstrack and downtrack dimensions can also be used.
    unpacked_bmask: a numpy array from  the band_mask function that can be used to mask band-specific pixels that have been interpolated.

    Returns:
    out_xr: an xarray.Dataset constructed based on the parameters provided.

    """
    # Grab granule filename to check product

    if type(filename) == s3fs.core.S3File:
        granule_id = filename.info()["name"].split("/", -1)[-1].split(".", -1)[0]
    elif type(filename) == HTTPFile:
        granule_id = filename.path.split("/", -1)[-1].split(".", -1)[0]
    else:
        granule_id = os.path.splitext(os.path.basename(filename))[0]

    # Read in Data as Xarray Datasets
    wvl_group = None

    # load swath dataset to xarray with lazy loading
    ds = xr.open_dataset(filename, engine=engine, chunks={})
    if swath_window is not None:
        # Slicing assumes dimensions are [downtrack, crosstrack]
        row_start = int(swath_window.row_off)
        row_end = row_start + int(swath_window.height)
        col_start = int(swath_window.col_off)
        col_end = col_start + int(swath_window.width)
        ds = ds.isel(downtrack=slice(row_start, row_end), crosstrack=slice(col_start, col_end)).load()
    # load location dataset to xarray with lazy loading
    loc = xr.open_dataset(filename, engine=engine, group="location", chunks={})
    if swath_window is not None:
        loc = loc.isel(downtrack=slice(row_start, row_end), crosstrack=slice(col_start, col_end)).load()

    # Check if mineral dataset and read in groups (only ds/loc for minunc)

    if "L2B_MIN_" in granule_id:
        wvl_group = "mineral_metadata"
    elif "L2B_MINUNC" not in granule_id:
        wvl_group = "sensor_band_parameters"

    wvl = None

    if wvl_group:
        wvl = xr.open_dataset(filename, engine=engine, group=wvl_group)

    # Building Flat Dataset from Components
    data_vars = {**ds.variables}

    # Format xarray coordinates based upon emit product (no wvl for mineral uncertainty)
    coords = {
        "downtrack": (["downtrack"], ds.downtrack.data),
        "crosstrack": (["crosstrack"], ds.crosstrack.data),
        **loc.variables,
    }

    product_band_map = {
        "L2B_MIN_": "name",
        "L2A_MASK_": "mask_bands",
        "L1B_OBS_": "observation_bands",
        "L2A_RFL_": "wavelengths",
        "L1B_RAD_": "wavelengths",
        "L2A_RFLUNCERT_": "wavelengths",
    }

    # if band := product_band_map.get(next((k for k in product_band_map.keys() if k in granule_id), 'unknown'), None):
    # coords['bands'] = wvl[band].data

    if wvl:
        coords = {**coords, **wvl.variables}

    # create xarray dataset with coordinates
    out_xr = xr.Dataset(data_vars=data_vars, coords=coords, attrs=ds.attrs)
    # set granule ID attribute
    out_xr.attrs["granule_id"] = granule_id

    if band := product_band_map.get(
        next((k for k in product_band_map.keys() if k in granule_id), "unknown"), None
    ):
        if "minerals" in list(out_xr.dims):
            out_xr = out_xr.swap_dims({"minerals": band})
            out_xr = out_xr.rename({band: "mineral_name"})
        else:
            out_xr = out_xr.swap_dims({"bands": band})

    # for each data layer, apply quality masks and set fill values to NaN
    for var in list(ds.data_vars):
        if qmask is not None:
            out_xr[var].data[qmask == 1] = -9999
        if unpacked_bmask is not None:
            out_xr[var].data[unpacked_bmask == 1] = -9999

    # orthorectify the swath cube if the ortho flag is set
    if ortho is True:
        # If GLT_array is not provided, extract and slice it from the NetCDF file
        if GLT_array is None:
            GLT_array, adjusted_GLT, swath_window = extract_GLT_array(
                filename=filename, 
                window=grid_window,
                adjust_indices=True
            )
        # orthorectify the swath cube to a grid cube
        out_xr = ortho_xr(
            out_xr,
            GLT_array=adjusted_GLT,
            fill_value=fill_value
        )
        # set `Orthorectified` attribute to True
        out_xr.attrs["Orthorectified"] = "True"

    return out_xr
