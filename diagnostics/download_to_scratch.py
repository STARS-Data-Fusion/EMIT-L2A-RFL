#!/usr/bin/env python
"""
Download EMIT granules to local scratch, then copy to network storage.

This avoids HDF5 file locking issues on network filesystems.
"""

import shutil
from pathlib import Path
import earthaccess

from EMITL2ARFL import (
    search_EMIT_L2A_RFL_granules,
    retrieve_EMIT_L2A_RFL_granule
)

# Configuration
LOCAL_SCRATCH = Path("/tmp/emit_download")  # Local filesystem
NETWORK_STORAGE = Path.home() / "data"      # Network filesystem (NFS/Lustre)

# Create directories
LOCAL_SCRATCH.mkdir(exist_ok=True)
NETWORK_STORAGE.mkdir(exist_ok=True)

# Login
earthaccess.login()

# Search for granule
granules = search_EMIT_L2A_RFL_granules(
    start_UTC="2023-01-29",
    end_UTC="2023-01-30",
    orbit=2302816,
    scene=3
)

if granules:
    print(f"Found granule: {granules[0]['umm']['GranuleUR']}")
    
    # Download to LOCAL scratch
    print(f"Downloading to local scratch: {LOCAL_SCRATCH}")
    granule = retrieve_EMIT_L2A_RFL_granule(
        remote_granule=granules[0],
        download_directory=str(LOCAL_SCRATCH)
    )
    
    print(f"✓ Downloaded and validated on local storage")
    print(f"  Files validated successfully in: {Path(granule.reflectance_filename).parent}")
    
    # Copy to network storage
    src_dir = Path(granule.reflectance_filename).parent
    granule_name = src_dir.name
    dest_dir = NETWORK_STORAGE / granule_name
    
    print(f"\nCopying to network storage: {dest_dir}")
    if dest_dir.exists():
        print(f"  Destination already exists, skipping copy")
    else:
        shutil.copytree(src_dir, dest_dir)
        print(f"✓ Copied to network storage")
    
    print(f"\nFiles available at:")
    print(f"  Local: {src_dir}")
    print(f"  Network: {dest_dir}")
else:
    print("No granules found")
