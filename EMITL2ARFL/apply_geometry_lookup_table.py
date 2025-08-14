from typing import Union

import numpy as np

from .constants import *
from .GLT import GeometryLookupTable

def apply_GLT(
            swath_array: np.ndarray, 
            GLT: np.ndarray, 
            fill_value: int = FILL_VALUE, 
            GLT_nodata_value: int = GLT_NODATA_VALUE) -> np.ndarray:
    """
    Orthorectifies satellite swath data using a Geometry Lookup Table (GLT).

    This function remaps raw satellite data (swath_array) onto a georeferenced output grid using a GLT.
    The GLT provides, for each output pixel (in geographic coordinates), the corresponding (row, column)
    indices in the original swath data. This process corrects for geometric distortions due to sensor
    viewing geometry, terrain, and other effects, producing an orthorectified (geospatially accurate) output.

    Geospatial process overview:
    - The GLT is a 3D array of shape (latitude, longitude, 2), where the last dimension holds the (row, column)
        indices in the swath data for each output pixel location in the georeferenced grid.
    - For each output pixel, if the GLT provides valid indices, the corresponding value(s) from the swath data
        are copied to the output. If the GLT entry is invalid (nodata), the output pixel is set to a fill value.
    - The result is an orthorectified array, spatially aligned to the geographic grid defined by the GLT.

    Parameters:
            swath_array (np.ndarray): Raw satellite data to be orthorectified. Shape: (rows, cols) or (rows, cols, bands).
            GLT_array (np.ndarray): Geometry Lookup Table. Shape: (latitude, longitude, 2), with (row, col) indices.
            fill_value (int, optional): Value to use for unmapped output pixels. Defaults to FILL_VALUE.
            GLT_nodata_value (int, optional): Value in GLT indicating invalid mapping. Defaults to GLT_NODATA_VALUE.

    Returns:
            np.ndarray: Orthorectified data array, shape (latitude, longitude, bands).

    Raises:
            ValueError: If input array dimensions are incompatible or GLT last dimension is not size 2.
    """
    # 1. Validate GLT shape: must be (latitude, longitude, 2) for geospatial mapping
    if GLT.ndim not in [2, 3] or (GLT.ndim == 3 and GLT.shape[-1] != 2):
        raise ValueError("GLT_array must be 2D or 3D with the last dimension of size 2.")

    # 2. Ensure swath_array is 3D for consistent band handling
    if swath_array.ndim == 2:
        # Convert single-band data to shape (rows, cols, 1)
        swath_array = swath_array[:, :, np.newaxis]

    # 3. Prepare output array shape: (latitude, longitude, bands)
    latitude_length, longitude_length = GLT.shape[:2]
    band_length = swath_array.shape[-1]
    ortho_array_shape = (latitude_length, longitude_length, band_length)

    # 4. Initialize output with fill_value (for unmapped pixels)
    ortho_array = np.full(ortho_array_shape, fill_value, dtype=np.float32)

    # 5. Identify valid GLT entries (where both row and col indices are not nodata)
    #    valid_GLT is a 2D boolean mask of shape (latitude, longitude)
    valid_GLT = np.all(GLT != GLT_nodata_value, axis=-1)

    # 6. Convert GLT indices from 1-based to 0-based (Python convention)
    zero_based_indices = GLT - 1

    # 7. For each valid output pixel, copy the corresponding swath data using GLT indices
    #    This remaps swath_array values to their georeferenced locations in ortho_array
    #    - zero_based_indices[..., 1] gives swath column indices
    #    - zero_based_indices[..., 0] gives swath row indices
    #    - valid_GLT mask selects only valid mappings
    col_indices = zero_based_indices[valid_GLT, 1]
    # print(f"col indices shape: {col_indices.shape} min: {np.nanmin(col_indices)} max: {np.nanmax(col_indices)}")
    row_indices = zero_based_indices[valid_GLT, 0]
    # print(f"row indices shape: {row_indices.shape} min: {np.nanmin(row_indices)} max: {np.nanmax(row_indices)}")
    # print(f"swath_array shape: {swath_array.shape}")

    ortho_array[valid_GLT, :] = swath_array[row_indices, col_indices, :]

    # 8. Replace any fill value of -9999 with np.nan for easier downstream analysis
    ortho_array = np.where(ortho_array == -9999, np.nan, ortho_array)

    # 9. Return the orthorectified, geospatially aligned output array
    return ortho_array
