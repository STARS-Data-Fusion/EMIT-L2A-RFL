class EMITNotAvailable(Exception):
    """Custom exception to indicate that EMIT data is not available for the specified parameters."""
    pass


class NetCDFValidationError(Exception):
    """Base exception for NetCDF file validation errors."""
    pass


class NetCDFFileNotFoundError(NetCDFValidationError, FileNotFoundError):
    """Exception raised when a NetCDF file does not exist at the specified path."""
    pass


class NetCDFEmptyFileError(NetCDFValidationError):
    """Exception raised when a NetCDF file exists but is empty (0 bytes)."""
    pass


class NetCDFCorruptedError(NetCDFValidationError):
    """Exception raised when a NetCDF file is corrupted or malformed."""
    pass


class NetCDFReadError(NetCDFValidationError):
    """Exception raised when a NetCDF file cannot be read due to I/O errors."""
    pass
