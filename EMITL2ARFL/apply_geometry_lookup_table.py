import numpy as np

from .constants import *

def apply_GLT(
        swath_array: np.ndarray, 
        GLT_array: np.ndarray, 
        fill_value: int = FILL_VALUE, 
        GLT_nodata_value: int = GLT_NODATA_VALUE) -> np.ndarray:
    """
    Applies a Geometry Lookup Table (GLT) to a numpy array representing satellite data, 
    to orthorectify it based on the GLT. This function supports input arrays of 2 or 3 dimensions.

    Parameters:
    - swath_array (np.ndarray): The input satellite data array to be orthorectified. 
                                Can be 2D (single band) or 3D (multiple bands).
    - GLT_array (np.ndarray): The Geometry Lookup Table array, which maps the input array's 
                              pixels to geographic locations. It is a 3-dimensional array in the 
                              shape of (latitude, longitude, 2), with the last dimension 
                              representing (row, column) indices.
    - fill_value (int, optional): The value used to fill the output array wherever the GLT 
                                  does not provide a mapping. Defaults to FILL_VALUE from constants.
    - GLT_nodata_value (int, optional): The value in the GLT_array that indicates no data or 
                                        invalid mapping. Pixels with this value in the GLT are 
                                        filled with `fill_value` in the output array. Defaults to 
                                        GLT_NODATA_VALUE from constants.

    Returns:
    - np.ndarray: A numpy array of the same number of dimensions as `swath_array`, containing the 
                  orthorectified data. The shape of the output array is determined by the dimensions 
                  of the GLT_array and the number of bands in `swath_array`.

    Raises:
    - ValueError: If the dimensions of the input arrays are not compatible or if the GLT_array does 
                  not have the last dimension of size 2.
    """

    # Ensure GLT_array has the correct shape
    if GLT_array.ndim not in [2, 3] or (GLT_array.ndim == 3 and GLT_array.shape[-1] != 2):
        raise ValueError("GLT_array must be 2D or 3D with the last dimension of size 2.")

    # Adjust swath_array dimensions if necessary
    if swath_array.ndim == 2:
        swath_array = swath_array[:, :, np.newaxis]

    # Extract dimensions for the output array
    latitude_length, longitude_length = GLT_array.shape[:2]
    band_length = swath_array.shape[-1]
    ortho_array_shape = (latitude_length, longitude_length, band_length)

    # Initialize the output array
    ortho_array = np.full(ortho_array_shape, fill_value, dtype=np.float32)

    # Identify valid GLT entries
    valid_GLT = np.all(GLT_array != GLT_nodata_value, axis=-1)

    # Adjust GLT indices to zero-based
    zero_based_indices = GLT_array - 1

    # Apply GLT to swath_array
    ortho_array[valid_GLT, :] = swath_array[zero_based_indices[valid_GLT, 1], zero_based_indices[valid_GLT, 0], :]
    ortho_array = np.where(ortho_array == -9999, np.nan, ortho_array)

    return ortho_array
