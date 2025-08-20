import numpy as np

def apply_qmask(array: np.ndarray, qmask: np.ndarray, fill_value: int = -9999) -> np.ndarray:
    """
    Apply a quality mask to a numpy array.

    Parameters:
    arr: numpy array to mask
    qmask: numpy array with the same spatial shape as arr, where 1 indicates masked pixels
    fill_value: value to assign to masked pixels (default: -9999)

    Returns:
    Masked numpy array
    """
    return np.where(qmask == 1, fill_value, array)
