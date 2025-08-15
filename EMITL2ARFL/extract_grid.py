import netCDF4
from affine import Affine
from rasterio.windows import Window

from rasters import RasterGrid

from .read_dimensions import read_dimensions

def extract_grid(
        filename: str,
        window: Window = None
        ) -> RasterGrid:
    with netCDF4.Dataset(filename, 'r') as ds:
        geotransform = ds.getncattr("geotransform")
    
    dimensions = read_dimensions(filename=filename)
    rows = int(dimensions["ortho_y"])
    cols = int(dimensions["ortho_x"])
    affine = Affine.from_gdal(*geotransform)
    
    grid = RasterGrid.from_affine(
        affine=affine,
        rows=rows,
        cols=cols
    )

    if window is not None:
        grid = grid.subset(window)

    return grid
