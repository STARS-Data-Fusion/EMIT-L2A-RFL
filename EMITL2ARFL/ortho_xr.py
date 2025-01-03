import numpy as np
import xarray as xr

from .constants import *
from .extract_GLT import extract_GLT
from .apply_geometry_lookup_table import apply_GLT
from .get_pixel_center_coords import get_pixel_center_coords

def ortho_xr(swath_ds: xr.Dataset, GLT_nodata_value: int = GLT_NODATA_VALUE, fill_value: int = FILL_VALUE) -> xr.Dataset:
    """
    This function uses `apply_GLT` to create an orthorectified xarray dataset.

    Parameters:
    swath_ds: an xarray dataset produced by emit_xarray
    GLT_nodata_value: no data value for the GLT tables, 0 by default
    fill_value: the fill value for EMIT datasets, -9999 by default

    Returns:
    ortho_ds: an orthocorrected xarray dataset.
    """
    # extract GLT
    GLT_array = extract_GLT(swath_ds)

    # List Variables
    var_list = list(swath_ds.data_vars)

    # Remove flat field from data vars - the flat field is only useful with additional information before orthorectification
    if "flat_field_update" in var_list:
        var_list.remove("flat_field_update")

    # Create empty dictionary for orthocorrected data vars
    data_vars = {}

    # Extract Rawspace Dataset Variable Values (Typically Reflectance)
    for var in var_list:
        swath_array = swath_ds[var].data
        swath_dimensions = swath_ds[var].dims

        # Apply GLT to dataset
        out_ds = apply_GLT(swath_array, GLT_array, GLT_nodata_value=GLT_nodata_value)

        # Update variables - Only works for 2 or 3 dimensional arays
        if swath_array.ndim == 2:
            out_ds = out_ds.squeeze()
            data_vars[var] = (["latitude", "longitude"], out_ds)
        else:
            data_vars[var] = (["latitude", "longitude", swath_dimensions[-1]], out_ds)

        del swath_array

    # Calculate Lat and Lon Vectors
    lon, lat = get_pixel_center_coords(swath_ds)  # Reorder this function to make sense in case of multiple variables

    # Apply GLT to elevation
    elev_ds = apply_GLT(swath_ds["elev"].data, GLT_array, fill_value=fill_value)

    # Delete glt_ds - no longer needed
    del GLT_array

    # Create Coordinate Dictionary
    coords = {
        "latitude": (["latitude"], lat),
        "longitude": (["longitude"], lon),
        **swath_ds.coords,
    }  # unpack to add appropriate coordinates

    # Remove Unnecessary Coords
    for key in ["downtrack", "crosstrack", "lat", "lon", "glt_x", "glt_y", "elev"]:
        del coords[key]

    # Add Orthocorrected Elevation
    coords["elev"] = (["latitude", "longitude"], np.squeeze(elev_ds))

    # Build Output xarray Dataset and assign data_vars array attributes
    ortho_ds = xr.Dataset(data_vars=data_vars, coords=coords, attrs=swath_ds.attrs)

    del out_ds

    # Assign Attributes from Original Datasets
    for var in var_list:
        ortho_ds[var].attrs = swath_ds[var].attrs

    ortho_ds.coords["latitude"].attrs = swath_ds["lat"].attrs
    ortho_ds.coords["longitude"].attrs = swath_ds["lon"].attrs
    ortho_ds.coords["elev"].attrs = swath_ds["elev"].attrs

    # Add Spatial Reference in recognizable format
    ortho_ds.rio.write_crs(swath_ds.spatial_ref, inplace=True)

    return ortho_ds
