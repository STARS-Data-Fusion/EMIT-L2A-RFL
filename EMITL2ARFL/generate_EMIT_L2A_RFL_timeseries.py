from typing import Union, List
import os
from datetime import date
import pandas as pd
from rasters import RasterGeometry
import logging
from os.path import exists, abspath, expanduser
import gc

from .constants import *
from .exceptions import *
from .retrieve_EMIT_L2A_RFL import retrieve_EMIT_L2A_RFL

logger = logging.getLogger(__name__)

def generate_EMIT_L2A_RFL_timeseries(
        start_date_UTC: Union[date, str],
        end_date_UTC: Union[date, str],
        geometry: RasterGeometry,
        output_directory: str,
        download_directory: str = DOWNLOAD_DIRECTORY,
        max_retries: int = 3,
        retry_delay: float = 2.0,
        threads: int = 1) -> List[str]:
    logger.info(f"generating EMIT L2A RFL timeseries from {start_date_UTC} to {end_date_UTC}")
    
    filenames = []
    
    # iterate through dates
    for date_UTC in pd.date_range(start=start_date_UTC, end=end_date_UTC):
        logger.info(f"processing date: {date_UTC}")
        
        # generate output filename
        output_filename = os.path.join(
            output_directory,
            f"EMIT_L2A_RFL_{date_UTC.strftime('%Y%m%d')}.tif"
        )
        
        if exists(abspath(expanduser(output_filename))):
            logger.info(f"output file already exists: {output_filename}")
            filenames.append(output_filename)
            continue
        
        try:
            # retrieve data for date
            merged_cube = retrieve_EMIT_L2A_RFL(
                date_UTC=date_UTC,
                geometry=geometry,
                download_directory=download_directory,
                max_retries=max_retries,
                retry_delay=retry_delay,
                threads=threads
            )
            
            logger.info(f"saving merged cube: {output_filename}")
            # save merged cube to file
            filenames.append(output_filename)
            merged_cube.to_geotiff(output_filename)
            
            # Explicitly clean up to free memory
            del merged_cube
            gc.collect()
            
        except EMITNotAvailable as e:
            logger.info(f"no EMIT granules available for date {date_UTC}")
            continue
    
    return filenames
