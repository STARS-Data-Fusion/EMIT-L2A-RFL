import netCDF4

def read_dimensions(filename: str) -> dict[str, int]:
    """
    Read the dimensions from a NetCDF file and return a dictionary mapping
    dimension names to their sizes.

    Args:
        filename (str): Path to the NetCDF file.

    Returns:
        dict[str, int]: Dictionary where keys are dimension names and values are sizes.
    """
    with netCDF4.Dataset(filename, 'r') as ds:
        return {dim_name: len(dim) for dim_name, dim in ds.dimensions.items()}
