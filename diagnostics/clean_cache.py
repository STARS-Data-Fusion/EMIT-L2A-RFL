#!/usr/bin/env python3
"""
Clean up corrupted EMIT cache files and validate remaining ones.
"""

import os
import logging
from pathlib import Path
from EMITL2ARFL.validate_NetCDF_file import validate_NetCDF_file

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

logger = logging.getLogger(__name__)

def clean_cache_directory(cache_dir):
    """Clean up corrupted NetCDF files in the cache directory."""
    cache_path = Path(cache_dir).expanduser()
    
    if not cache_path.exists():
        logger.info(f"Cache directory does not exist: {cache_dir}")
        return
    
    logger.info(f"Scanning cache directory: {cache_dir}")
    
    corrupted_files = []
    valid_files = []
    total_size_corrupted = 0
    
    # Find all .nc files
    nc_files = list(cache_path.rglob('*.nc'))
    logger.info(f"Found {len(nc_files)} NetCDF files")
    
    for nc_file in nc_files:
        try:
            validate_NetCDF_file(nc_file, file_type="NetCDF")
            valid_files.append(nc_file)
            logger.info(f"✓ Valid: {nc_file.name}")
        except Exception as e:
            corrupted_files.append(nc_file)
            file_size = nc_file.stat().st_size
            total_size_corrupted += file_size
            logger.warning(f"✗ Corrupted: {nc_file.name} ({file_size / (1024**2):.2f} MB) - {type(e).__name__}")
    
    # Summary
    logger.info(f"\n{'='*80}")
    logger.info(f"Summary:")
    logger.info(f"  Valid files: {len(valid_files)}")
    logger.info(f"  Corrupted files: {len(corrupted_files)}")
    logger.info(f"  Space wasted by corrupted files: {total_size_corrupted / (1024**3):.2f} GB")
    logger.info(f"{'='*80}\n")
    
    if corrupted_files:
        response = input("Delete corrupted files? (yes/no): ")
        if response.lower() in ['yes', 'y']:
            for nc_file in corrupted_files:
                try:
                    logger.info(f"Deleting: {nc_file}")
                    nc_file.unlink()
                    
                    # Also delete parent directory if empty
                    parent = nc_file.parent
                    if parent != cache_path and not list(parent.iterdir()):
                        logger.info(f"Deleting empty directory: {parent}")
                        parent.rmdir()
                except Exception as e:
                    logger.error(f"Could not delete {nc_file}: {e}")
            
            logger.info(f"✓ Cleaned up {len(corrupted_files)} corrupted files")
        else:
            logger.info("Skipped deletion")
    else:
        logger.info("✓ No corrupted files found!")

if __name__ == '__main__':
    import sys
    
    if len(sys.argv) > 1:
        cache_dir = sys.argv[1]
    else:
        cache_dir = '~/data/EMIT_L2A_RFL'
    
    clean_cache_directory(cache_dir)
