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

from EMITL2ARFL import *


def main():
    """Main function to generate Kings Canyon EMIT time series."""
    
    # Configure logging to see info messages
    logging.basicConfig(level=logging.INFO, format='%(name)s - %(levelname)s - %(message)s')
    
    # Configuration parameters
    start_date_UTC = "2022-08-01"
    end_date_UTC = "2025-11-20"
    download_directory = "~/data/EMIT_download"
    output_directory = "~/data/Kings Canyon EMIT"
    
    # Load Upper Kings area of interest
    print("Loading Upper Kings area of interest...")
    gdf = gpd.read_file("upper_kings.kml")
    print(f"Loaded geometry: {gdf.geometry[0]}")
    
    # Create UTM bounding box and raster grid
    print("Creating UTM bounding box and raster grid...")
    bbox_UTM = rt.Polygon(gdf.unary_union).UTM.bbox
    print(f"UTM bounding box: {bbox_UTM}")
    
    grid = rt.RasterGrid.from_bbox(bbox_UTM, cell_size=60, crs=bbox_UTM.crs)
    print(f"Raster grid: {grid}")
    
    # Log into earthaccess using netrc credentials
    print("Logging into earthaccess...")
    earthaccess.login(strategy="netrc", persist=True)
    
    # Generate EMIT L2A reflectance time series
    print("Generating EMIT L2A reflectance time series...")
    filenames = generate_EMIT_L2A_RFL_timeseries(
        start_date_UTC=start_date_UTC,
        end_date_UTC=end_date_UTC,
        geometry=grid,
        output_directory=output_directory
    )
    
    print(f"Generated {len(filenames)} files:")
    for filename in filenames:
        print(f"  {filename}")
    
    # Process and display each file
    print("\nProcessing generated files...")
    for filename in filenames:
        print(f"\nProcessing: {filename}")
        try:
            raster = MultiRaster.open(filename)
            print(f"Successfully opened raster: {raster}")
            # Note: display() function removed as it's typically for Jupyter notebooks
            # You can add specific processing or visualization code here if needed
        except Exception as e:
            print(f"Error processing {filename}: {e}")


if __name__ == "__main__":
    main()



