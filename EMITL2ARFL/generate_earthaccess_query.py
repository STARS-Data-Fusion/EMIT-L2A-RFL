from typing import Union
from datetime import date, datetime

import earthaccess

from rasters import Point, Polygon, RasterGeometry

from .temporally_constrain_earthaccess_query import temporally_constrain_earthaccess_query
from .spatially_constrain_earthaccess_query import spatially_constrain_earthaccess_query

def generate_earthaccess_query(
        concept_ID: str,
        start_UTC: Union[date, datetime, str] = None,
        end_UTC: Union[date, datetime, str] = None,
        geometry: Union[Point, Polygon, RasterGeometry] = None,
        readable_granule_name: str = None) -> earthaccess.search.DataGranules:
    """
    Generate an EarthAccess query with optional temporal and spatial constraints.

    Parameters:
    concept_ID (str): The concept ID to include in the query.
    start_UTC (Union[date, datetime, str], optional): The start date/time for the temporal constraint.
    end_UTC (Union[date, datetime, str], optional): The end date/time for the temporal constraint.
    geometry (Union[Point, Polygon, RasterGeometry], optional): The spatial geometry for the spatial constraint.
    readable_granule_name (str, optional): The search pattern to include in the query.

    Returns:
    earthaccess.search.DataGranules: The generated EarthAccess query with the specified constraints.
    """
    # create blank query
    query = earthaccess.granule_query()
    # include concept ID in query
    query = query.concept_id(concept_ID)

    if start_UTC is not None and end_UTC is not None:
        # include time range in query
        query = temporally_constrain_earthaccess_query(query, start_UTC, end_UTC)

    if geometry is not None:
        # spatially constrain query
        query = spatially_constrain_earthaccess_query(query, geometry)

    if readable_granule_name is not None:
        # include a search pattern
        query = query.readable_granule_name(readable_granule_name)

    return query