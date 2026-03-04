#!/usr/bin/env python3
"""
Simple test to download a single granule to /tmp and diagnose any issues.
"""

import logging
from datetime import date
import earthaccess
from EMITL2ARFL.search_EMIT_L2A_RFL_granules import search_EMIT_L2A_RFL_granules
from EMITL2ARFL.retrieve_EMIT_L2A_RFL_granule import retrieve_EMIT_L2A_RFL_granule

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

logger = logging.getLogger(__name__)

def main():
    # Use /tmp for download
    download_dir = '/tmp/emit_test'
    
    logger.info("Logging into earthaccess...")
    earthaccess.login()
    
    # Search for a single granule
    logger.info("Searching for granule from 2023-01-29...")
    granules = search_EMIT_L2A_RFL_granules(
        start_UTC=date(2023, 1, 29),
        end_UTC=date(2023, 1, 30)
    )
    
    if not granules:
        logger.error("No granules found")
        return
    
    logger.info(f"Found {len(granules)} granules")
    logger.info(f"Testing with first granule: {granules[0]['umm']['GranuleUR']}")
    
    # Try to retrieve it
    logger.info(f"Downloading to: {download_dir}")
    try:
        granule = retrieve_EMIT_L2A_RFL_granule(
            remote_granule=granules[0],
            download_directory=download_dir,
            max_retries=3,
            retry_delay=2.0
        )
        
        logger.info("✓ SUCCESS!")
        logger.info(f"Reflectance: {granule.reflectance_filename}")
        logger.info(f"Mask: {granule.mask_filename}")
        logger.info(f"Uncertainty: {granule.uncertainty_filename}")
        
    except Exception as e:
        logger.error(f"✗ FAILED: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    main()
