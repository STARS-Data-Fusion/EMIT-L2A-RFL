import xarray as xr
from rasterio.windows import Window

from rasters import RasterGeometry

from .extract_grid import extract_grid

from .extract_GLT_array import extract_GLT_array
from .GLT import GeometryLookupTable
from .emit_xarray import emit_xarray
from .constants import GLT_NODATA_VALUE

def extract_GLT(
        filename: str,
        geometry: RasterGeometry = None,
        window: Window = None,
        GLT_nodata_value: int = GLT_NODATA_VALUE):
    """
    Wrapper function to extract a GLT object from a NetCDF file.
    """
    processing_subset = geometry is not None or window is not None

    scene_grid = extract_grid(filename=filename)

    if processing_subset:
        if geometry is None:
            geometry = extract_grid(
                filename, 
                window=window
            )

        if window is None:
            window = scene_grid.window(geometry)
            subset_grid = scene_grid.subset(window)
    else:
        geometry = scene_grid

    GLT_array = extract_GLT_array(
        filename, 
        window=window,
        GLT_nodata_value=GLT_nodata_value
    )

    # a GLT is a grid of swath indices that needs to be georeferenced as a grid

    GLT = GeometryLookupTable(GLT_array=GLT_array, geometry=geometry)

    if processing_subset:
        adjusted_GLT = GLT.adjust_indices(window)
        swath_window: Window = GLT.swath_window

        # calculate the row/col shape of the subset defined by the swath_window Window object
        swath_shape = (swath_window.height, swath_window.width)

        # check that the maximum row/col indices in the adjusted GLT are within the swath shape
        if adjusted_GLT is not None:
            if adjusted_GLT.max_row >= swath_shape[0] or adjusted_GLT.max_col >= swath_shape[1]:
                raise ValueError(
                    f"Adjusted GLT indices are out of bounds for the swath shape: "
                    f"max_row={adjusted_GLT.max_row}, max_col={adjusted_GLT.max_col}, "
                    f"swath_shape={swath_shape}"
                )

    else:
        adjusted_GLT = None
        swath_window = None

    return GLT, adjusted_GLT, swath_window
