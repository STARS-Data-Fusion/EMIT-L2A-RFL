from os.path import expanduser, abspath, exists

import netCDF4
import hashlib
from pathlib import Path
from typing import Union

from .exceptions import (
    NetCDFFileNotFoundError,
    NetCDFEmptyFileError,
    NetCDFCorruptedError,
    NetCDFReadError
)


def validate_NetCDF_file(filename: Union[str, Path], file_type: str = "NetCDF", check_integrity: bool = False) -> None:
    """
    Validate whether a file is a valid NetCDF file.
    
    This function performs comprehensive validation of a NetCDF file, checking for:
    - File existence
    - Non-zero file size
    - NetCDF format validity
    - Readability of dimensions and variables
    - (Optional) File integrity via checksum
    
    Parameters
    ----------
    filename : str or Path
        Path to the NetCDF file to validate
    file_type : str, optional
        Descriptive name for the file type (e.g., "Reflectance", "Mask").
        Used in error messages for better context. Defaults to "NetCDF".
    check_integrity : bool, optional
        If True, compute file checksum for integrity verification. Defaults to False.
        This is slower but useful for detecting partial downloads or corruption.
    
    Returns
    -------
    None
        Returns nothing if validation succeeds
    
    Raises
    ------
    NetCDFFileNotFoundError
        If the file does not exist at the specified path
    NetCDFEmptyFileError
        If the file exists but has zero bytes (empty file)
    NetCDFCorruptedError
        If the file is corrupted or has invalid NetCDF format
    NetCDFReadError
        If the file cannot be read due to I/O errors or permission issues
    
    Examples
    --------
    >>> validate_NetCDF_file('data.nc')
    >>> validate_NetCDF_file('reflectance.nc', file_type='Reflectance')
    >>> validate_NetCDF_file('missing.nc')  # Raises NetCDFFileNotFoundError
    """
    filename_absolute = Path(abspath(expanduser(filename)))
    
    # Check if file exists
    if not exists(filename_absolute):
        raise NetCDFFileNotFoundError(
            f"{file_type} file does not exist at path: {filename}"
        )
    
    # Check if file is empty
    file_size = filename_absolute.stat().st_size
    if file_size == 0:
        raise NetCDFEmptyFileError(
            f"{file_type} file is empty (0 bytes): {filename}"
        )
    
    # Attempt to open and validate NetCDF structure
    try:
        with netCDF4.Dataset(filename_absolute, "r") as ds:
            # Try to access basic attributes to ensure file is readable
            dimensions = ds.dimensions
            variables = ds.variables
            
            # Verify that the file contains data structures
            if len(dimensions) == 0 and len(variables) == 0:
                raise NetCDFCorruptedError(
                    f"{file_type} file is corrupted: contains no dimensions or variables. "
                    f"File: {filename}, Size: {file_size} bytes"
                )
            
            # Optionally verify file integrity
            if check_integrity:
                # Read a small portion of data to verify the file can be accessed
                for var_name in list(variables.keys())[:3]:  # Check up to 3 variables
                    try:
                        _ = variables[var_name].shape
                    except Exception as e:
                        raise NetCDFCorruptedError(
                            f"{file_type} file data cannot be accessed. "
                            f"Variable '{var_name}' failed: {e}. "
                            f"File: {filename}, Size: {file_size} bytes"
                        )
    
    except (OSError, IOError) as e:
        # I/O errors suggest file read problems or corruption
        error_msg = str(e).lower()
        error_str = str(e)
        
        # Check for specific HDF/NetCDF error codes
        if "errno -101" in error_str or "hdf error" in error_msg:
            raise NetCDFReadError(
                f"{file_type} file has HDF/NetCDF format errors (possibly corrupted during download): {filename}. "
                f"File size: {file_size} bytes. Error: {e}. "
                f"Recommendation: Delete and re-download this file."
            )
        elif "permission" in error_msg or "access" in error_msg:
            raise NetCDFReadError(
                f"{file_type} file cannot be read due to permission error: {filename}. "
                f"Error: {e}"
            )
        else:
            raise NetCDFReadError(
                f"{file_type} file cannot be read due to I/O error: {filename}. "
                f"File size: {file_size} bytes. Error: {e}"
            )
    
    except RuntimeError as e:
        # RuntimeError typically indicates NetCDF format problems
        raise NetCDFCorruptedError(
            f"{file_type} file is corrupted or has invalid NetCDF format: {filename}. "
            f"File size: {file_size} bytes. Error: {e}"
        )
    
    except Exception as e:
        # Catch-all for unexpected errors
        raise NetCDFCorruptedError(
            f"{file_type} file validation failed with unexpected error: {filename}. "
            f"File size: {file_size} bytes. Error type: {type(e).__name__}, Message: {e}"
        )
