from typing import Union, List
from datetime import date, datetime

import earthaccess

from rasters import Point, Polygon, RasterGeometry

from .generate_earthaccess_query import generate_earthaccess_query

def search_earthaccess_granules(
        concept_ID: str,
        start_UTC: Union[date, datetime, str] = None,
        end_UTC: Union[date, datetime, str] = None,
        geometry: Union[Point, Polygon, RasterGeometry] = None,
        readable_granule_name: str = None) -> List[earthaccess.search.DataGranule]:
    """
    Search for Earth data granules using the earthaccess library.

    Parameters:
    concept_ID (str): The concept ID of the dataset to search.
    start_UTC (Union[date, datetime, str], optional): The start date/time for the search. Defaults to None.
    end_UTC (Union[date, datetime, str], optional): The end date/time for the search. Defaults to None.
    geometry (Union[Point, Polygon, RasterGeometry], optional): The spatial geometry for the search. Defaults to None.
    readable_granule_name (str, optional): The search pattern to filter the search. Defaults to None.

    Returns:
    List[earthaccess.search.DataGranule]: A list of data granules matching the search criteria.
    """
    query = generate_earthaccess_query(
        concept_ID=concept_ID,
        start_UTC=start_UTC,
        end_UTC=end_UTC,
        geometry=geometry,
        readable_granule_name=readable_granule_name
    )

    granules = query.get()
    
    return granules