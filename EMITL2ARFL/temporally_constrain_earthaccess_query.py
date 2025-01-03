from typing import Union
from datetime import datetime, date
from dateutil import parser

import earthaccess

__author__ = "Gregory H. Halverson, Evan Davis"

def start_of_day(d: Union[date, str]) -> datetime:
    """
    Convert a date or date string to the earliest datetime of that date.

    Args:
        d (Union[date, str]): The input date or date string.

    Returns:
        datetime: The earliest datetime of the input date.
    """
    if isinstance(d, str):
        d = parser.parse(d)

    if isinstance(d, datetime):
        d = d.date()

    date_string = d.strftime("%Y-%m-%d")
    dt = parser.parse(f"{date_string}T00:00:00Z")

    return dt

def end_of_day(d: Union[date, str]) -> datetime:
    """
    Convert a date or date string to the latest datetime of that date.

    Args:
        d (Union[date, str]): The input date or date string.

    Returns:
        datetime: The latest datetime of the input date.
    """
    if isinstance(d, str):
        d = parser.parse(d)

    if isinstance(d, datetime):
        d = d.date()

    date_string = d.strftime("%Y-%m-%d")
    dt = parser.parse(f"{date_string}T23:59:59Z")

    return dt

def temporally_constrain_earthaccess_query(
        query: earthaccess.search.DataGranules,
        start_date: Union[date, str],
        end_date: Union[date, str]) -> earthaccess.search.DataGranules:
    """
    Temporally constrain an earthaccess query to the specified date range.

    Args:
        query (earthaccess.search.DataGranules): The earthaccess query object to be constrained.
        start_date (Union[date, str]): The start date of the constraint.
        end_date (Union[date, str]): The end date of the constraint.

    Returns:
        earthaccess.search.DataGranules: The constrained earthaccess query object.
    """
    return query.temporal(start_of_day(start_date), end_of_day(end_date))
