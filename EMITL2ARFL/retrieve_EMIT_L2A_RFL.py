from typing import Union
from datetime import date, datetime
import gc

import rasters as rt
from rasters import MultiRaster, Point, Polygon, RasterGeometry
import logging

from .constants import *
from .exceptions import *
from .search_EMIT_L2A_RFL_granules import search_EMIT_L2A_RFL_granules
from .retrieve_EMIT_L2A_RFL_granule import retrieve_EMIT_L2A_RFL_granule

logger = logging.getLogger(__name__)

def retrieve_EMIT_L2A_RFL(
        date_UTC: Union[date, datetime, str],
        geometry: Union[Point, Polygon, RasterGeometry],
        download_directory: str = DOWNLOAD_DIRECTORY,
        max_retries: int = 3,
        retry_delay: float = 2.0,
        threads: int = 1) -> MultiRaster:
    search_results = search_EMIT_L2A_RFL_granules(
        start_UTC=date_UTC,
        end_UTC=date_UTC,
        geometry=geometry
    )
    
    if len(search_results) == 0:
        raise EMITNotAvailable(f"No EMIT L2A RFL granules found for date {date_UTC} and specified geometry.")
    
    logger.info(f"found {len(search_results)} granules for date {date_UTC}")
    
    # Process granules one at a time to minimize memory usage
    subset_cubes = []
    for search_result in search_results:
        try:
            granule = retrieve_EMIT_L2A_RFL_granule(
                remote_granule=search_result,
                download_directory=download_directory,
                max_retries=max_retries,
                retry_delay=retry_delay,
                threads=threads
            )
            subset_cube = granule.reflectance(geometry=geometry)
            subset_cubes.append(subset_cube)
            # Clean up granule object immediately
            del granule
        except Exception as e:
            logger.warning(f"Failed to process granule: {e}")
            continue
    
    if not subset_cubes:
        raise EMITNotAvailable(f"No valid EMIT L2A RFL data retrieved for date {date_UTC}")
    
    merged_cube = rt.mosaic(subset_cubes, geometry=geometry)
    
    # Clean up subset cubes before returning
    del subset_cubes
    gc.collect()

    return merged_cube
