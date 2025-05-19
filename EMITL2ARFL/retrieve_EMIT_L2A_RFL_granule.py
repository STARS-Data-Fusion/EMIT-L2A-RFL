import posixpath
from os.path import join, expanduser, abspath

import earthaccess

from .constants import *
from .granule import EMITL2ARFL
from .search_EMIT_L2A_RFL_granules import search_EMIT_L2A_RFL_granules

def retrieve_EMIT_L2A_RFL_granule(
        granule: earthaccess.search.DataGranule = None,
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
    if granule is None and orbit is not None and scene is not None:
        remote_granules = search_EMIT_L2A_RFL_granules(orbit=orbit, scene=scene)
        
        if len(remote_granules) == 0:
            raise ValueError(f"no EMIT L2A RFL granule found for orbit {orbit} and scene {scene}")
        
        granule = remote_granules[0]
    
    if granule is None:
        raise ValueError("either granule or orbit and scene must be provided")  

    # parse granule ID
    granule_ID = posixpath.splitext(posixpath.basename(granule.data_links()[0]))[0]

    # make sure that this is an EMIT L2A Reflectance collection 1 granule
    if not granule_ID.startswith("EMIT_L2A_RFL_001_"):
        raise ValueError("The provided granule is not an EMIT L2A Reflectance collection 1 granule.")

    # create a subdirectory for the granule
    directory = join(download_directory, granule_ID)
    # download the granule files to the directory
    earthaccess.download(granule.data_links(), local_path=abspath(expanduser(directory)))
    # wrap the directory in an EMITL2ARFL object
    granule = EMITL2ARFL(directory)

    return granule