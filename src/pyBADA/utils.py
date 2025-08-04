import numpy as np
import xarray as xr
import pandas as pd
from decimal import Decimal, ROUND_HALF_UP


def _round_scalar(x, dec):
    """Round a single scalar value using half-up rounding."""
    quant = Decimal('1.' + '0'*dec)
    d = Decimal(str(x))
    return float(d.quantize(quant, rounding=ROUND_HALF_UP))


def proper_round(num, dec=0):
    """Forced half-up rounding, vectorized over numpy arrays, pandas, and xarray.

    :param num: Input scalar, array-like, pandas Series/DataFrame, or xarray DataArray
    :param dec: Number of decimal places
    :returns: Rounded values with half-up rule; preserves infinities.
    """
    # Scalar helper handling infinities
    def _f(v):
        # Preserve infinities
        try:
            if np.isinf(v):
                return v
        except Exception:
            pass
        return _round_scalar(v, dec)

    # xarray DataArray
    if isinstance(num, xr.DataArray):
        arr = num.data
        rounded = np.vectorize(_f)(arr)
        return xr.DataArray(rounded, coords=num.coords, dims=num.dims)

    # pandas Series
    if isinstance(num, pd.Series):
        return num.map(_f)

    # pandas DataFrame
    if isinstance(num, pd.DataFrame):
        return num.applymap(_f)

    # numpy array
    if isinstance(num, np.ndarray):
        return np.vectorize(_f)(num)

    # scalar fallback
    return _f(num)


def to_numpy(x):
    """
    Convert xarray.DataArray, pandas.Series/DataFrame, or any array-like
    to a NumPy array of floats.
    """
    if isinstance(x, xr.DataArray):
        return x.values
    if isinstance(x, (pd.Series, pd.DataFrame)):
        return x.to_numpy(copy=False)
    return np.asarray(x, dtype=float)


def _extract(x):
    if isinstance(x, xr.DataArray):
        return x.data
    if isinstance(x, (pd.Series, pd.DataFrame)):
        return x.values
    return np.asarray(x, dtype=float)

def _wrap(x):
    if isinstance(x, xr.DataArray):
        return xr.DataArray(x, coords=v.coords, dims=v.dims)
    if isinstance(x, pd.Series):
        return pd.Series(x, index=v.index, name=v.name)
    if isinstance(x, pd.DataFrame):
        return pd.DataFrame(x, index=v.index, columns=v.columns)
    return x

def _vectorized_wrapper(core_func, *args):
    """Generic vectorized wrapper for functions with N inputs."""
    # Extract raw arrays
    arrs = [_extract(a) for a in args]
    core = core_func(*arrs)

    # Determine output type from first argument
    first = args[0]
    if isinstance(first, xr.DataArray) and all(isinstance(a, xr.DataArray) for a in args):
        da = xr.DataArray(core, coords=first.coords, dims=first.dims)
        return proper_round(da, 10)

    if isinstance(first, pd.Series) and all(isinstance(a, pd.Series) for a in args):
        return proper_round(pd.Series(core, index=first.index, name=first.name), 10)

    if isinstance(first, pd.DataFrame) and all(isinstance(a, pd.DataFrame) for a in args):
        return proper_round(pd.DataFrame(core, index=first.index, columns=first.columns), 10)

    # numpy array or scalar fallback
    return proper_round(core, 10)

def checkArgument(argument, **kwargs):
    if kwargs.get(argument) is not None:
        return kwargs.get(argument)
    else:
        raise TypeError("Missing " + argument + " argument")