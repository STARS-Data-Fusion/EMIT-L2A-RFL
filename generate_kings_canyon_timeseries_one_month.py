"""
Generate EMIT L2A Reflectance time series for Kings Canyon area.

This script downloads and processes EMIT L2A reflectance data for the Upper Kings
area of interest over a specified date range and generates a time series of
reflectance data.
"""

import logging
from os.path import join

import earthaccess
import geopandas as gpd
import matplotlib.pyplot as plt
import rasters as rt

import colored_logging as cl

from EMITL2ARFL import *


# Configuration parameters
start_date_UTC = "2022-08-01"
end_date_UTC = "2022-08-30"
download_directory = "/tmp/EMIT_download"
output_directory = "~/data/Kings_Canyon_EMIT"

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

# Generate EMIT L2A reflectance time series
logger.info("Generating EMIT L2A reflectance time series...")
filenames = generate_EMIT_L2A_RFL_timeseries(
    start_date_UTC=start_date_UTC,
    end_date_UTC=end_date_UTC,
    geometry=grid,
    download_directory=download_directory,
    output_directory=output_directory
)

logger.info(f"Generated {len(filenames)} files:")
for filename in filenames:
    logger.info(f"  {filename}")

logger.info("\nTime series generation complete!")

