"""
Common unit conversions module
"""

from datetime import datetime
from math import pi
from time import mktime

import numpy as np
import pandas as pd
import xarray as xr

from pyBADA import utils


def _unit_converter(factor):
    """
    Returns a function that multiplies its input by `factor`,
    vectorized through utils._extract / utils._wrap.
    """

    def converter(val):
        # pull out the raw array
        arr = utils._extract(val)
        # do the multiplication
        core = arr * factor
        # pass both the result and the original to _wrap
        return utils._wrap(core, original=val)

    return converter


# Linear unit conversions (factor-based)
_factor_map = {
    "ft2m": 0.3048,
    "nm2m": 1852.0,
    "h2s": 3600.0,
    "kt2ms": 0.514444,
    "lb2kg": 0.453592,
    "deg2rad": pi / 180.0,
    "rad2deg": 180.0 / pi,
    "m2ft": 1 / 0.3048,
    "m2nm": 1 / 1852.0,
    "s2h": 1 / 3600.0,
    "ms2kt": 1 / 0.514444,
    "kg2lb": 1 / 0.453592,
    "hp2W": 745.699872,
}

# Dynamically create each converter in the module namespace
for _name, _factor in _factor_map.items():
    globals()[_name] = _unit_converter(_factor)


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
        out = val.astype("datetime64[s]").astype(int)
        out.attrs.update(units="s", long_name="POSIX timestamp")
        return out

    if isinstance(val, (pd.Series, pd.DataFrame)):
        ts = pd.to_datetime(val)
        arr = ts.values.astype("datetime64[s]").astype(int)
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
        dt64 = arr.astype("datetime64[s]")
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
        dt64 = val.astype("datetime64[s]")
        return dt64.dt.strftime("%Y-%m-%d %H:%M:%S")

    if isinstance(val, (pd.Series, pd.DataFrame)):
        ts = pd.to_datetime(val, unit="s")
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
