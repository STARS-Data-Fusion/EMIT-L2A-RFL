from typing import Union
from datetime import date, datetime

import rasters as rt
from rasters import MultiRaster, Point, Polygon, RasterGeometry

from .constants import *
from .exceptions import *
from .search_EMIT_L2A_RFL_granules import search_EMIT_L2A_RFL_granules
from .retrieve_EMIT_L2A_RFL_granule import retrieve_EMIT_L2A_RFL_granule

def retrieve_EMIT_L2A_RFL(
        date_UTC: Union[date, datetime, str],
        geometry: Union[Point, Polygon, RasterGeometry],
        download_directory: str = DOWNLOAD_DIRECTORY) -> MultiRaster:
    search_results = search_EMIT_L2A_RFL_granules(
        start_UTC=date_UTC,
        end_UTC=date_UTC,
        geometry=geometry
    )
    
    if len(search_results) == 0:
        raise EMITNotAvailable(f"No EMIT L2A RFL granules found for date {date_UTC} and specified geometry.")
    
    granules = [
        retrieve_EMIT_L2A_RFL_granule(
            remote_granule=search_result,
            download_directory=download_directory
        ) 
        for search_result 
        in search_results
    ]
    
    subset_cubes = [granule.reflectance(geometry=geometry) for granule in granules]
    
    merged_cube = rt.mosaic(subset_cubes, geometry=geometry)

    return merged_cube
