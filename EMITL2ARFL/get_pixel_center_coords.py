import numpy as np
import xarray as xr
from rasterio.windows import Window

# Function to Calculate the center of pixel Lat and Lon Coordinates of the GLT grid
def get_pixel_center_coords(
        ds: xr.Dataset,
        grid_window: Window = None):
    """
    This function calculates the gridded latitude and longitude pixel centers for the dataset using the geotransform and GLT arrays.

    Parameters:
    ds: an emit dataset opened with emit_xarray function

    Returns:
    x_geo, y_geo: longitude and latitude pixel centers of glt (gridded data)

    """
    # Retrieve GLT
    GT = ds.geotransform
    # Get Shape of GLT
    dim_x = ds.glt_x.shape[1]
    dim_y = ds.glt_y.shape[0]
    
    # Build Arrays containing pixel centers
    x_geo = (GT[0] + 0.5 * GT[1]) + np.arange(dim_x) * GT[1]
    y_geo = (GT[3] + 0.5 * GT[5]) + np.arange(dim_y) * GT[5]

    # Subset arrays if grid_window is provided
    if grid_window is not None:
        x_geo = x_geo[grid_window.col_off : grid_window.col_off + grid_window.width]
        y_geo = y_geo[grid_window.row_off : grid_window.row_off + grid_window.height]

    return x_geo, y_geo
