import posixpath
from os.path import join, expanduser, abspath, exists

import earthaccess

from .constants import *
from .EMITL2ARFLGranule import EMITL2ARFLGranule
from .find_EMIT_L2A_RFL_granule import find_EMIT_L2A_RFL_granule

def retrieve_EMIT_L2A_RFL_granule(
        remote_granule: earthaccess.search.DataGranule = None,
        orbit: int = None,
        scene: int = None, 
        download_directory: str = DOWNLOAD_DIRECTORY) -> EMITL2ARFLGranule:
    """
    Retrieve an EMIT L2A Reflectance granule.

    This function retrieves an EMIT L2A Reflectance granule based on the provided granule, orbit, and scene.
    If the granule is not provided, it searches for the granule using the orbit and scene parameters.
    The granule is then downloaded to the specified directory and wrapped in an EMITL2ARFL object.

    Args:
        granule (earthaccess.search.DataGranule, optional): The granule to retrieve. Defaults to None.
        orbit (int, optional): The orbit number to search for the granule. Defaults to None.
        scene (int, optional): The scene number to search for the granule. Defaults to None.
        download_directory (str, optional): The directory to download the granule files to. Defaults to ".".

    Returns:
        EMITL2ARFL: The retrieved EMIT L2A Reflectance granule wrapped in an EMITL2ARFL object.

    Raises:
        ValueError: If no granule is found for the provided orbit and scene, or if the provided granule is not an EMIT L2A Reflectance collection 1 granule.
    """
    if remote_granule is None and orbit is not None and scene is not None:
        remote_granule = find_EMIT_L2A_RFL_granule(granule=remote_granule, orbit=orbit, scene=scene)
    
    if remote_granule is None:
        raise ValueError("either granule or orbit and scene must be provided")

    base_filenames = [posixpath.basename(URL) for URL in remote_granule.data_links()]

    if not all([exists(filename) for filename in base_filenames]):
        # parse granule ID
        granule_ID = posixpath.splitext(posixpath.basename(remote_granule.data_links()[0]))[0]

        # make sure that this is an EMIT L2A Reflectance collection 1 granule
        if not granule_ID.startswith("EMIT_L2A_RFL_001_"):
            raise ValueError("The provided granule is not an EMIT L2A Reflectance collection 1 granule.")

        # create a subdirectory for the granule
        directory = join(download_directory, granule_ID)
        # download the granule files to the directory
        earthaccess.download(remote_granule.data_links(), local_path=abspath(expanduser(directory)))

    # Use base_filenames to determine the local filenames for RFL, MASK, and RFLUNCERT
    if 'directory' not in locals():
        granule_ID = posixpath.splitext(posixpath.basename(remote_granule.data_links()[0]))[0]
        directory = join(download_directory, granule_ID)

    abs_directory = abspath(expanduser(directory))

    # Map the base filenames to their full local paths
    local_files = [join(abs_directory, fname) for fname in base_filenames]
    reflectance_filename = next((f for f in local_files if '_RFL_' in f and not '_RFLUNCERT_' in f), None)
    mask_filename = next((f for f in local_files if '_MASK_' in f), None)
    uncertainty_filename = next((f for f in local_files if '_RFLUNCERT_' in f), None)

    if not (reflectance_filename and mask_filename and uncertainty_filename):
        raise FileNotFoundError('Could not find all required NetCDF files (RFL, MASK, RFLUNCERT) in the granule directory.')

    # Check that each file exists on disk
    missing_files = [f for f in [reflectance_filename, mask_filename, uncertainty_filename] if not exists(f)]

    if missing_files:
        raise FileNotFoundError(f"The following required files do not exist: {missing_files}")

    local_granule = EMITL2ARFLGranule(
        reflectance_filename=reflectance_filename,
        mask_filename=mask_filename,
        uncertainty_filename=uncertainty_filename
    )

    return local_granule