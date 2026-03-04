"""
Generate EMIT L2A Reflectance time series for Kings Canyon area in chunks.

This script processes the data in yearly chunks to avoid memory issues on systems
with limited RAM (e.g., 1 GB). Run this script multiple times with different year
parameters, or let it process all years sequentially.
"""

import logging
from os.path import join
import sys

import earthaccess
import geopandas as gpd
import matplotlib.pyplot as plt
import rasters as rt

import colored_logging as cl

from EMITL2ARFL import *


# Configuration parameters
YEARS_TO_PROCESS = [2022, 2023, 2024, 2025]  # Add/remove years as needed
download_directory = "/tmp/EMIT_download"
output_directory = "~/data/Kings Canyon EMIT"

# Load Upper Kings area of interest
logger.info("Loading Upper Kings area of interest...")
gdf = gpd.read_file("upper_kings.kml")
logger.info(f"Loaded geometry: {gdf.geometry[0]}")

# Create UTM bounding box and raster grid
logger.info("Creating UTM bounding box and raster grid...")
bbox_UTM = rt.Polygon(gdf.unary_union).UTM.bbox
logger.info(f"UTM bounding box: {bbox_UTM}")

grid = rt.RasterGrid.from_bbox(bbox_UTM, cell_size=60, crs=bbox_UTM.crs)
logger.info(f"Raster grid: {grid}")

# Log into earthaccess using netrc credentials
logger.info("Logging into earthaccess...")
earthaccess.login(strategy="netrc", persist=True)

# Process each year separately
all_filenames = []
for year in YEARS_TO_PROCESS:
    # Determine date range for this year
    if year == 2022:
        start_date = "2022-08-01"  # EMIT started in August 2022
    else:
        start_date = f"{year}-01-01"
    
    if year == 2025:
        end_date = "2025-11-20"  # Current end date
    else:
        end_date = f"{year}-12-31"
    
    logger.info(f"\n{'='*60}")
    logger.info(f"Processing year {year}: {start_date} to {end_date}")
    logger.info(f"{'='*60}\n")
    
    try:
        # Generate EMIT L2A reflectance time series for this year
        filenames = generate_EMIT_L2A_RFL_timeseries(
            start_date_UTC=start_date,
            end_date_UTC=end_date,
            geometry=grid,
            download_directory=download_directory,
            output_directory=output_directory
        )
        
        logger.info(f"Generated {len(filenames)} files for year {year}")
        all_filenames.extend(filenames)
        
    except Exception as e:
        logger.error(f"Error processing year {year}: {e}")
        logger.error("Continuing with next year...")
        continue

# Summary
logger.info(f"\n{'='*60}")
logger.info(f"PROCESSING COMPLETE")
logger.info(f"{'='*60}")
logger.info(f"Total files generated: {len(all_filenames)}")
logger.info(f"Output directory: {output_directory}")
