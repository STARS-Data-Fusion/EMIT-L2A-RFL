import posixpath
from os.path import join, expanduser, abspath, exists

import earthaccess

from .constants import *
from .granule import EMITL2ARFL
from .find_EMIT_L2A_RFL_granule import find_EMIT_L2A_RFL_granule

def retrieve_EMIT_L2A_RFL_granule(
    remote_granule: earthaccess.search.DataGranule = None,
    orbit: int = None,
    scene: int = None, 
    download_directory: str = DOWNLOAD_DIRECTORY) -> EMITL2ARFL:
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

    # wrap the directory in an EMITL2ARFL object
    local_granule = EMITL2ARFL(directory)

    return local_granule