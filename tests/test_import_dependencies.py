import pytest

# List of dependencies
dependencies = [
    "colored_logging",
    "earthaccess",
    "geopandas",
    "netCDF4",
    "dateutil",
    "rasters",
    "rioxarray",
    "spectral",
    "xarray"
]

# Generate individual test functions for each dependency
@pytest.mark.parametrize("dependency", dependencies)
def test_dependency_import(dependency):
    __import__(dependency)
