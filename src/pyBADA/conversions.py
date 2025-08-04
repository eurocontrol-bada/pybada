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
    Convert from feet to meters, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param val: Input value(s) in feet (ft). Scalar, numpy array,
                pandas Series/DataFrame, or xarray.DataArray.
    :type val: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Converted value(s) in meters (m), matching input type.
    """
    factor = 0.3048
    arr_val = utils._extract(val)
    core = arr_val * factor
    return utils._wrap(core)


def nm2m(val):
    """
    Convert from nautical miles to meters, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param val: Input value(s) in nautical miles (nm). Scalar, numpy array,
                pandas Series/DataFrame, or xarray.DataArray.
    :type val: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Converted value(s) in meters (m), matching input type.
    """
    factor = 1852.0
    arr_val = utils._extract(val)
    core = arr_val * factor
    return utils._wrap(core)


def h2s(val):
    """
    Convert from hours to seconds, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param val: Input value(s) in hours (h). Scalar, numpy array,
                pandas Series/DataFrame, or xarray.DataArray.
    :type val: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Converted value(s) in seconds (s), matching input type.
    """
    factor = 3600.0
    arr_val = utils._extract(val)
    core = arr_val * factor
    return utils._wrap(core)


def kt2ms(val):
    """
    Convert from knots to meters per second, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param val: Input value(s) in knots (kt). Scalar, numpy array,
                pandas Series/DataFrame, or xarray.DataArray.
    :type val: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Converted value(s) in meters per second (m/s), matching input type.
    """
    factor = 0.514444
    arr_val = utils._extract(val)
    core = arr_val * factor
    return utils._wrap(core)


def lb2kg(val):
    """
    Convert from pounds to kilograms, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param val: Input value(s) in pounds (lb). Scalar, numpy array,
                pandas Series/DataFrame, or xarray.DataArray.
    :type val: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Converted value(s) in kilograms (kg), matching input type.
    """
    factor = 0.453592
    arr_val = utils._extract(val)
    core = arr_val * factor
    return utils._wrap(core)


def deg2rad(val):
    """
    Convert from decimal degrees to radians, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param val: Input value(s) in decimal degrees (°). Scalar, numpy array,
                pandas Series/DataFrame, or xarray.DataArray.
    :type val: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Converted value(s) in radians (rad), matching input type.
    """
    factor = pi / 180.0
    arr_val = utils._extract(val)
    core = arr_val * factor
    return utils._wrap(core)

def rad2deg(val):
    """
    Convert from decimal degrees to radians, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param val: Input value(s) in decimal degrees (°). Scalar, numpy array,
                pandas Series/DataFrame, or xarray.DataArray.
    :type val: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Converted value(s) in radians (rad), matching input type.
    """
    factor = 180 / pi
    arr_val = utils._extract(val)
    core = arr_val * factor
    return utils._wrap(core)


def m2ft(val):
    """
    Convert from meters to feet, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param val: Input value(s) in meters (m). Scalar, numpy array,
                pandas Series/DataFrame, or xarray.DataArray.
    :type val: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Converted value(s) in feet (ft), matching input type.
    """
    factor = 1 / 0.3048
    arr_val = utils._extract(val)
    core = arr_val * factor
    return utils._wrap(core)


def m2nm(val):
    """
    Convert from meters to nautical miles, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param val: Input value(s) in meters (m). Scalar, numpy array,
                pandas Series/DataFrame, or xarray.DataArray.
    :type val: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Converted value(s) in nautical miles (nm), matching input type.
    """
    factor = 1 / 1852.0
    arr_val = utils._extract(val)
    core = arr_val * factor
    return utils._wrap(core)


def s2h(val):
    """
    Convert from seconds to hours, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param val: Input value(s) in seconds (s). Scalar, numpy array,
                pandas Series/DataFrame, or xarray.DataArray.
    :type val: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Converted value(s) in hours (h), matching input type.
    """
    factor = 1 / 3600.0
    arr_val = utils._extract(val)
    core = arr_val * factor
    return utils._wrap(core)


def ms2kt(val):
    """
    Convert from meters per second to knots, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param val: Input value(s) in meters per second (m/s). Scalar, numpy array,
                pandas Series/DataFrame, or xarray.DataArray.
    :type val: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Converted value(s) in knots (kt), matching input type.
    """
    factor = 1 / 0.514444
    arr_val = utils._extract(val)
    core = arr_val * factor
    return utils._wrap(core)


def kg2lb(val):
    """
    Convert from kilograms to pounds, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param val: Input value(s) in kilograms (kg). Scalar, numpy array,
                pandas Series/DataFrame, or xarray.DataArray.
    :type val: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Converted value(s) in pounds (lb), matching input type.
    """
    factor = 1 / 0.453592
    arr_val = utils._extract(val)
    core = arr_val * factor
    return utils._wrap(core)


def hp2W(val):
    """
    Convert from horsepower to watts, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param val: Input value(s) in horsepower (hp). Scalar, numpy array,
                pandas Series/DataFrame, or xarray.DataArray.
    :type val: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Converted value(s) in watts (W), matching input type.
    """
    factor = 745.699872
    arr_val = utils._extract(val)
    core = arr_val * factor
    return utils._wrap(core)


def date2posix(val):
    """
    Convert date(s) to POSIX timestamp in seconds since 1970-01-01, vectorized for
    numpy arrays, pandas Series/DataFrame, and xarray.DataArray.

    :param val: Input date(s). Can be string, datetime, numpy array of strings/datetimes,
                pandas Series/DataFrame of datetime-like, or xarray.DataArray of datetime64.
    :type val: str or datetime or array-like or pandas Series/DataFrame or xarray.DataArray
    :returns: POSIX timestamp(s) in seconds, matching input type.
    """
    # xarray DataArray of datetime64
    if isinstance(val, xr.DataArray):
        out = val.astype('datetime64[s]').astype(int)
        out.attrs.update(units='s', long_name='POSIX timestamp')
        return out

    if isinstance(val, (pd.Series, pd.DataFrame)):
        ts = pd.to_datetime(val)
        arr = ts.values.astype('datetime64[s]').astype(int)
        if isinstance(val, pd.Series):
            return pd.Series(arr, index=val.index, name=val.name)
        return pd.DataFrame(arr, index=val.index, columns=val.columns)

    if isinstance(val, str):
        dt = datetime.strptime(val, "%Y-%m-%d %H:%M:%S")
        return mktime(dt.timetuple())

    try:
        if isinstance(val, datetime):
            return mktime(val.timetuple())
    except NameError:
        pass

    arr = np.asarray(val)
    try:
        dt64 = arr.astype('datetime64[s]')
        return dt64.astype(int)
    except Exception:
        flat = arr.ravel()
        result = []
        for v in flat:
            if isinstance(v, str):
                dt = datetime.strptime(v, "%Y-%m-%d %H:%M:%S")
                result.append(mktime(dt.timetuple()))
            else:
                result.append(mktime(v.timetuple()))
        out = np.array(result)
        return out.reshape(arr.shape)


def unix2date(val):
    """
    Convert POSIX timestamp(s) in seconds to date string(s) in "%Y-%m-%d %H:%M:%S" format, vectorized for
    numpy arrays, pandas Series/DataFrame, and xarray.DataArray.

    :param val: Input POSIX timestamp(s). Float, numpy array, pandas Series/DataFrame of floats,
                or xarray.DataArray of ints.
    :type val: float or array-like or pandas Series/DataFrame or xarray.DataArray
    :returns: Date string(s) in "%Y-%m-%d %H:%M:%S", matching input type.
    """
    if isinstance(val, xr.DataArray):
        dt64 = val.astype('datetime64[s]')
        return dt64.dt.strftime("%Y-%m-%d %H:%M:%S")

    if isinstance(val, (pd.Series, pd.DataFrame)):
        ts = pd.to_datetime(val, unit='s')
        return ts.dt.strftime("%Y-%m-%d %H:%M:%S")

    arr = np.asarray(val, dtype=float)
    flat = arr.ravel()
    result = []
    for v in flat:
        dt = datetime.fromtimestamp(int(v))
        result.append(dt.strftime("%Y-%m-%d %H:%M:%S"))
    res_arr = np.array(result)
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
