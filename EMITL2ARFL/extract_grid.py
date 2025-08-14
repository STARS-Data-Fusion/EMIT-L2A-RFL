import netCDF4
from affine import Affine

from rasters import RasterGrid

def extract_grid(filename: str) -> RasterGrid:
    with netCDF4.Dataset(filename, 'r') as ds:
        geotransform = ds.getncattr("geotransform")
    
    affine = Affine.from_gdal(geotransform)
    grid = RasterGrid.from_affine(affine)

    return grid
