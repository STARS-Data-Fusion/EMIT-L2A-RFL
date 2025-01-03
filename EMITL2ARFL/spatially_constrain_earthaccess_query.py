from typing import Union

import earthaccess

from rasters import Point, Polygon, RasterGeometry

__author__ = "Gregory H. Halverson, Evan Davis"

def spatially_constrain_earthaccess_query(
        query: earthaccess.search.DataGranules, 
        geometry: Union[Point, Polygon, RasterGeometry]) -> earthaccess.search.DataGranules:
    """
    Add spatial constraints to an earthaccess query based on the provided geometry.

    Parameters:
    query (earthaccess.search.DataGranules): The initial earthaccess query object.
    geometry (Union[Point, Polygon, RasterGeometry]): The spatial geometry to constrain the query. 
        Can be a Point, Polygon, or RasterGeometry.

    Returns:
    earthaccess.search.DataGranules: The constrained earthaccess query object.
    """
    # Add spatial constraints to the query if provided
    if isinstance(geometry, Point):
        # If the target geometry is a Point, add a point constraint to the query
        query = query.point(geometry.x, geometry.y)
    
    if isinstance(geometry, Polygon):
        # If the target geometry is a Polygon, add a polygon constraint to the query
        ring = geometry.exterior
        
        # Ensure the ring is counter-clockwise
        if not ring.is_ccw:
            ring = ring.reverse()
        
        coordinates = ring.coords
        
        # Add the polygon coordinates to the query
        query = query.polygon(coordinates)
    
    if isinstance(geometry, RasterGeometry):
        # If the target geometry is a RasterGeometry, add a polygon constraint to the query
        ring = geometry.corner_polygon_latlon.exterior
        
        # Ensure the ring is counter-clockwise
        if not ring.is_ccw:
            ring = ring.reverse()
        
        coordinates = ring.coords
        
        # Add the polygon coordinates to the query
        query = query.polygon(coordinates)
    
    return query