"""
Common unit conversions module
"""

from math import pi
from datetime import datetime
from time import mktime


def ft2m(val):
    """
    This function converts from ft to m s

    :param val: value in ft
    :returns: vaue in m
    """
    return round(float(val) * 0.3048, 10)


def nm2m(val):
    """
    This function converts from nautical miles to m

    :param val: value in nautical miles
    :returns: vaue in m
    """
    return val * 1852.0


def h2s(val):
    """
    This function converts from hours to m seconds

    :param val: value in hours
    :returns: vaue in seconds
    """
    return val * 3600.0


def kt2ms(val):
    """
    This function converts from kt to m s^-1

    :param val: value in kt
    :returns: vaue in m s^-1
    """

    if val is None:
        return None
    else:
        return round(float(val) * 0.514444, 10)


def lb2kg(val):
    """This function converts from lb to kg.

    :param val: value in lb
    :returns: vaue in kg
    """
    return val * 0.453592


def deg2rad(val):
    """This function converts from decimal degrees to radians.

    :param val: value in decimal degrees
    :returns: vaue in radians
    """
    return val * pi / 180.0


def m2ft(val):
    """This function converts from meters to feets.

    :param val: value in meters
    :returns: value in feets
    """

    return round(float(val) / 0.3048, 10)


def m2nm(val):
    """This function converts from meters to nautical miles.

    :param val: value in meters
    :returns: value in nautical miles
    """
    return val / 1852.0


def s2h(val):
    """This function converts from seconds to hours.

    :param val: value in seconds
    :returns: value in hours
    """
    return val / 3600.0


def ms2kt(val):
    """This function converts from m s^-1 to kt.

    :param val: value in m s^-1
    :returns: value in kt
    """

    if val is None:
        return None
    else:
        return round(float(val) / 0.514444, 10)


def kg2lb(val):
    """This function converts from kg to lb.

    :param val: value in kg
    :returns: value in lb
    """
    return val / 0.453592


def rad2deg(val):
    """This function converts from radians to decimal degrees.

    :param val: value in radians
    :returns: value in decimal degrees
    """
    return val / pi * 180.0


def hp2W(val):
    """This function converts from horsepower to watts.

    :param val: value in horsepower
    :returns: value in watts
    """
    return val * 745.699872


def date2posix(val):
    """This function converts a date format to posix.

    :param val: date in %Y-%m-%d %H:%M:%S format
    :returns: posix time referenced to 01-01-1970 [s]
    """
    return mktime(datetime.strptime(val, "%Y-%m-%d %H:%M:%S").timetuple())


def unix2date(val):
    """This function converts posix to date format.

    :param val: time referenced to 01-01-1970 [s]
    :returns: date in %Y-%m-%d %H:%M:%S format
    """
    return datetime.fromtimestamp(int(val)).strftime("%Y-%m-%d %H:%M:%S")


convertFrom = {
    "unix": unix2date,
    "ft": ft2m,
    "nm": nm2m,
    "h": h2s,
    "kt": kt2ms,
    "lb": lb2kg,
    "deg": deg2rad,
    "date": date2posix,
    "rad": rad2deg,
    "ms": ms2kt,
    "m": m2ft,
    "kg": kg2lb,
    "s": s2h,
}
