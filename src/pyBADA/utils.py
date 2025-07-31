import numpy as np
import xarray as xr
import pandas as pd
from decimal import Decimal, ROUND_HALF_UP


def proper_round(num, dec=0):
    """
    Implementation of forced and precise half_up rounding
    """
    quant = Decimal('1.' + '0'*dec)
    d = Decimal(str(num))
    return float(d.quantize(quant, rounding=ROUND_HALF_UP))

def to_numpy(x):
    """
    Convert xarray.DataArray, pandas.Series/DataFrame, or any array-like
    to a NumPy array of floats.
    """
    if isinstance(x, xr.DataArray):
        return x.values
    if isinstance(x, (pd.Series, pd.DataFrame)):
        return x.to_numpy()
    return np.asarray(x, dtype=float)

def checkArgument(argument, **kwargs):
    if kwargs.get(argument) is not None:
        return kwargs.get(argument)
    else:
        raise TypeError("Missing " + argument + " argument")