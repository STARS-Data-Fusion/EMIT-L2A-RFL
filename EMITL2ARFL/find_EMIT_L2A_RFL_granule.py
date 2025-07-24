from .search_EMIT_L2A_RFL_granules import search_EMIT_L2A_RFL_granules
import earthaccess

def find_EMIT_L2A_RFL_granule(granule: earthaccess.search.DataGranule = None, orbit: int = None, scene: int = None):
    """
    Find an EMIT L2A Reflectance granule by granule object or by orbit and scene.

    Args:
        granule (earthaccess.search.DataGranule, optional): The granule to retrieve. Defaults to None.
        orbit (int, optional): The orbit number to search for the granule. Defaults to None.
        scene (int, optional): The scene number to search for the granule. Defaults to None.

    Returns:
        earthaccess.search.DataGranule: The found granule.

    Raises:
        ValueError: If no granule is found for the provided orbit and scene, or if neither granule nor orbit/scene are provided.
    """
    if granule is None and orbit is not None and scene is not None:
        remote_granules = search_EMIT_L2A_RFL_granules(orbit=orbit, scene=scene)
    
        if len(remote_granules) == 0:
            raise ValueError(f"no EMIT L2A RFL granule found for orbit {orbit} and scene {scene}")
    
        granule = remote_granules[0]
    
    if granule is None:
        raise ValueError("either granule or orbit and scene must be provided")
    
    return granule
