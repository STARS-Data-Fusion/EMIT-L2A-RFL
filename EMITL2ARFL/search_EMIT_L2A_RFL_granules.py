from typing import Union, List
from datetime import date, datetime

import earthaccess

from rasters import Point, Polygon, RasterGeometry

from .constants import *
from .search_earthaccess_granules import search_earthaccess_granules

def search_EMIT_L2A_RFL_granules(
        start_UTC: Union[date, datetime, str] = None,
        end_UTC: Union[date, datetime, str] = None,
        geometry: Union[Point, Polygon, RasterGeometry] = None,
        orbit: int = None,
        scene: int = None) -> List[earthaccess.search.DataGranule]:
    """
    Search for EMIT L2A Reflectance granules using specified parameters.

    Parameters:
    start_UTC (Union[date, datetime, str], optional): The start date and time for the search.
    end_UTC (Union[date, datetime, str], optional): The end date and time for the search.
    geometry (Union[Point, Polygon, RasterGeometry], optional): The spatial geometry for the search.
    orbit (int, optional): The orbit number for the search.
    scene (int, optional): The scene number for the search.

    Returns:
    List[earthaccess.search.DataGranule]: A list of found data granules.
    """
    readable_granule_name: str = None
        
    if orbit is not None and scene is not None:
        readable_granule_name=f"*{orbit}_{scene:03d}*"

    granules = search_earthaccess_granules(
        concept_ID=EMIT_L2A_REFLECTANCE_CONCEPT_ID,
        start_UTC=start_UTC,
        end_UTC=end_UTC,
        geometry=geometry,
        readable_granule_name=readable_granule_name
    )

    return granules
