"""
Common unit conversions module
"""

from datetime import datetime
from math import pi
from time import mktime

import numpy as np
import xarray as xr
import pandas as pd

from pyBADA import utils

def ft2m(val):
    """
    Convert from feet to meters, preserving metadata for xarray and pandas inputs.

    :param val: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Input value(s) in feet.
    :returns: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Converted value(s) in meters, matching the type and metadata of `val`.
    """
    factor = 0.3048

    if isinstance(val, xr.DataArray):
        return (val * factor).round(10)

    if isinstance(val, (pd.Series, pd.DataFrame)):
        return val.mul(factor).round(10)

    arr = utils.to_numpy(val) * factor
    out = np.round(arr, 10)
    return out.item() if (isinstance(out, np.ndarray) and out.ndim == 0) else out


def nm2m(val):
    """
    Convert from nautical miles to meters, preserving metadata for xarray and pandas inputs.

    :param val: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Input value(s) in nautical miles.
    :returns: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Converted value(s) in meters, matching the type and metadata of `val`.
    """
    factor = 1852.0

    if isinstance(val, xr.DataArray):
        return (val * factor).round(10)

    if isinstance(val, (pd.Series, pd.DataFrame)):
        return val.mul(factor).round(10)

    arr = utils.to_numpy(val) * factor
    out = np.round(arr, 10)
    return out.item() if (isinstance(out, np.ndarray) and out.ndim == 0) else out


def h2s(val):
    """
    Convert from hours to seconds, preserving metadata for xarray and pandas inputs.

    :param val: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Input value(s) in hours.
    :returns: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Converted value(s) in seconds, matching the type and metadata of `val`.
    """
    factor = 3600.0

    if isinstance(val, xr.DataArray):
        return val * factor

    if isinstance(val, (pd.Series, pd.DataFrame)):
        return val.mul(factor)

    arr = utils.to_numpy(val) * factor
    out = np.round(arr, 10)
    return out.item() if (isinstance(out, np.ndarray) and out.ndim == 0) else out


def kt2ms(val):
    """
    Convert from knots to meters per second, preserving metadata for xarray and pandas inputs.

    :param val: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Input value(s) in knots.
    :returns: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Converted value(s) in meters per second, matching the type and metadata of `val`.
    """
    factor = 0.514444

    if isinstance(val, xr.DataArray):
        return (val * factor).round(10)

    if isinstance(val, (pd.Series, pd.DataFrame)):
        return val.mul(factor).round(10)

    arr = utils.to_numpy(val) * factor
    out = np.round(arr, 10)
    return out.item() if (isinstance(out, np.ndarray) and out.ndim == 0) else out


def lb2kg(val):
    """
    Convert from pounds to kilograms, preserving metadata for xarray and pandas inputs.

    :param val: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Input value(s) in pounds.
    :returns: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Converted value(s) in kilograms, matching the type and metadata of `val`.
    """
    factor = 0.453592

    if isinstance(val, xr.DataArray):
        return (val * factor).round(10)

    if isinstance(val, (pd.Series, pd.DataFrame)):
        return val.mul(factor).round(10)

    arr = utils.to_numpy(val) * factor
    out = np.round(arr, 10)
    return out.item() if (isinstance(out, np.ndarray) and out.ndim == 0) else out


def deg2rad(val):
    """This function converts from decimal degrees to radians.

    :param val: value in decimal degrees
    :returns: vaue in radians
    """
    return val * pi / 180.0

def deg2rad(val):
    """
    Convert from decimal degrees to radians, preserving metadata for xarray and pandas inputs.

    :param val: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Input value(s) in decimal degrees.
    :returns: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Converted value(s) in radians, matching the type and metadata of `val`.
    """
    factor = pi / 180.0

    if isinstance(val, xr.DataArray):
        return (val * factor)

    if isinstance(val, (pd.Series, pd.DataFrame)):
        return val.mul(factor)

    arr = utils.to_numpy(val) * factor
    out = np.round(arr, 10)
    return out.item() if (isinstance(out, np.ndarray) and out.ndim == 0) else out


def m2ft(val):
    """
    Convert from meters to feet, preserving metadata for xarray and pandas inputs.

    :param val: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Input value(s) in meters.
    :returns: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Converted value(s) in feet, matching the type and metadata of `val`.
    """
    factor = 1 / 0.3048

    if isinstance(val, xr.DataArray):
        return (val * factor).round(10)

    if isinstance(val, (pd.Series, pd.DataFrame)):
        return val.mul(factor).round(10)

    arr = utils.to_numpy(val) * factor
    out = np.round(arr, 10)
    return out.item() if (isinstance(out, np.ndarray) and out.ndim == 0) else out


def m2nm(val):
    """
    Convert from meters to nautical miles, preserving metadata for xarray and pandas inputs.

    :param val: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Input value(s) in meters.
    :returns: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Converted value(s) in nautical miles, matching the type and metadata of `val`.
    """
    factor = 1 / 1852.0

    if isinstance(val, xr.DataArray):
        return (val * factor)

    if isinstance(val, (pd.Series, pd.DataFrame)):
        return val.mul(factor)

    arr = utils.to_numpy(val) * factor
    out = np.round(arr, 10)
    return out.item() if (isinstance(out, np.ndarray) and out.ndim == 0) else out


def s2h(val):
    """
    Convert from seconds to hours, preserving metadata for xarray and pandas inputs.

    :param val: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Input value(s) in seconds.
    :returns: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Converted value(s) in hours, matching the type and metadata of `val`.
    """
    factor = 1 / 3600.0

    if isinstance(val, xr.DataArray):
        return (val * factor)

    if isinstance(val, (pd.Series, pd.DataFrame)):
        return val.mul(factor)

    arr = utils.to_numpy(val) * factor
    out = np.round(arr, 10)
    return out.item() if (isinstance(out, np.ndarray) and out.ndim == 0) else out


def ms2kt(val):
    """
    Convert from meters per second to knots, preserving metadata for xarray and pandas inputs.

    :param val: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Input value(s) in meters per second.
    :returns: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Converted value(s) in knots, matching the type and metadata of `val`.
    """
    if val is None:
        return None

    factor = 1 / 0.514444

    if isinstance(val, xr.DataArray):
        return (val * factor).round(10)

    if isinstance(val, (pd.Series, pd.DataFrame)):
        return val.mul(factor).round(10)

    arr = utils.to_numpy(val) * factor
    out = np.round(arr, 10)
    return out.item() if (isinstance(out, np.ndarray) and out.ndim == 0) else out


def kg2lb(val):
    """
    Convert from kilograms to pounds, preserving metadata for xarray and pandas inputs.

    :param val: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Input value(s) in kilograms.
    :returns: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Converted value(s) in pounds, matching the type and metadata of `val`.
    """
    factor = 1 / 0.453592

    if isinstance(val, xr.DataArray):
        return (val * factor).round(10)

    if isinstance(val, (pd.Series, pd.DataFrame)):
        return val.mul(factor).round(10)

    arr = utils.to_numpy(val) * factor
    out = np.round(arr, 10)
    return out.item() if (isinstance(out, np.ndarray) and out.ndim == 0) else out


def rad2deg(val):
    """
    Convert from radians to decimal degrees, preserving metadata for xarray and pandas inputs.

    :param val: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Input value(s) in radians.
    :returns: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Converted value(s) in decimal degrees, matching the type and metadata of `val`.
    """
    factor = 180.0 / pi

    if isinstance(val, xr.DataArray):
        return (val * factor)

    if isinstance(val, (pd.Series, pd.DataFrame)):
        return val.mul(factor)

    arr = utils.to_numpy(val) * factor
    out = np.round(arr, 10)
    return out.item() if (isinstance(out, np.ndarray) and out.ndim == 0) else out


def hp2W(val):
    """
    Convert from horsepower to watts, preserving metadata for xarray and pandas inputs.

    :param val: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Input value(s) in horsepower.
    :returns: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Converted value(s) in watts, matching the type and metadata of `val`.
    """
    factor = 745.699872

    if val is None:
        return None

    if isinstance(val, xr.DataArray):
        return (val * factor).round(10)

    if isinstance(val, (pd.Series, pd.DataFrame)):
        return val.mul(factor).round(10)

    arr = utils.to_numpy(val) * factor
    out = np.round(arr, 10)
    return out.item() if (isinstance(out, np.ndarray) and out.ndim == 0) else out


def date2posix(val):
    """
    Convert date(s) to POSIX timestamp in seconds since 1970-01-01.

    :param val: str, datetime.datetime, np.ndarray of strings/datetimes,
                pandas.Series/DataFrame of datetime-like, or xr.DataArray of datetime dtype
    :returns: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        POSIX time in seconds, matching the type and metadata of `val`.
    """
    # xarray DataArray of datetime64
    if isinstance(val, _xr.DataArray):
        # convert datetime64 to seconds since epoch
        out = val.astype('datetime64[s]').astype(int)
        out.attrs.update(units='s', long_name='POSIX timestamp')
        return out

    # pandas Series/DataFrame
    if isinstance(val, (_pd.Series, _pd.DataFrame)):
        ts = _pd.to_datetime(val)
        # view nanoseconds and convert to seconds
        arr = ts.values.astype('datetime64[s]').astype(int)
        if isinstance(val, _pd.Series):
            return _pd.Series(arr, index=val.index, name=val.name)
        return _pd.DataFrame(arr, index=val.index, columns=val.columns)

    # single string or datetime
    if isinstance(val, str):
        dt = datetime.strptime(val, "%Y-%m-%d %H:%M:%S")
        return mktime(dt.timetuple())
    try:
        # handle datetime.datetime
        if isinstance(val, datetime):
            return mktime(val.timetuple())
    except NameError:
        pass

    # numpy array or list of strings/datetimes
    arr = _np.asarray(val)
    # attempt to parse or convert
    try:
        dt64 = arr.astype('datetime64[s]')
        out = dt64.astype(int)
        return out
    except Exception:
        # try string parse fallback
        flat = arr.ravel()
        result = []
        for v in flat:
            if isinstance(v, str):
                dt = datetime.strptime(v, "%Y-%m-%d %H:%M:%S")
                result.append(mktime(dt.timetuple()))
            else:
                # assume datetime
                result.append(mktime(v.timetuple()))
        out = _np.array(result)
        return out.reshape(arr.shape)


def unix2date(val):
    """
    Convert POSIX timestamp(s) in seconds to date string(s) in "%Y-%m-%d %H:%M:%S" format, preserving metadata for xarray and pandas inputs.

    :param val: float, np.ndarray, pandas.Series/DataFrame of floats, or xr.DataArray of ints
        Input POSIX timestamp(s).
    :returns: str, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Converted date string(s) matching the format and metadata of `val`.
    """
    # xarray DataArray
    if isinstance(val, _xr.DataArray):
        # convert seconds to datetime64 then to string
        dt64 = val.astype('datetime64[s]')
        out = dt64.dt.strftime("%Y-%m-%d %H:%M:%S")
        return out

    # pandas Series/DataFrame
    if isinstance(val, (_pd.Series, _pd.DataFrame)):
        ts = _pd.to_datetime(val, unit='s')
        out = ts.dt.strftime("%Y-%m-%d %H:%M:%S")
        return out

    # numpy array or scalar
    arr = _np.asarray(val, dtype=float)
    flat = arr.ravel()
    result = []
    for v in flat:
        dt = datetime.fromtimestamp(int(v))
        result.append(dt.strftime("%Y-%m-%d %H:%M:%S"))
    res_arr = _np.array(result)
    if res_arr.ndim == 0:
        return res_arr.item()
    return res_arr.reshape(arr.shape)


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
