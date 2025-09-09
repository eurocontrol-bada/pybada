import numbers
from decimal import ROUND_HALF_UP, Decimal

import numpy as np
import pandas as pd
import xarray as xr


def _round_scalar(x, dec):
    """Round a single scalar value using half-up rounding."""
    quant = Decimal("1." + "0" * dec)
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
    if isinstance(x, numbers.Real) and not isinstance(
        x, (np.generic, np.ndarray)
    ):
        return float(x)
    if isinstance(x, xr.DataArray):
        return x.data
    if isinstance(x, (pd.Series, pd.DataFrame)):
        return x.values
    return np.asarray(x, dtype=float)


def _broadcast(*arrays):
    """
    Broadcast any number of array-like inputs to a common shape.
    - If *all* inputs are Python scalars, just return them as a tuple.
    - If exactly two inputs, and one is 1-D while the other is N-D (N>1),
      then align the 1-D array to whichever axis of the N-D array has the same length.
    - Otherwise prepend leading singleton dims to match trailing-dims broadcasting.
    """
    # 1) Scalar passthrough
    if all(
        isinstance(a, numbers.Real)
        and not isinstance(a, (np.generic, np.ndarray))
        for a in arrays
    ):
        return tuple(arrays)

    # Convert everything up front
    arrs = [np.asarray(a) for a in arrays]

    # Special‐case exactly two inputs, one 1-D and one N-D
    if len(arrs) == 2:
        a0, a1 = arrs
        # identify which is 1-D and which is higher-D
        if a0.ndim == 1 and a1.ndim > 1:
            arrs[0] = _align_1d_to_nd(a0, a1)
            return np.broadcast_arrays(arrs[0], a1)
        if a1.ndim == 1 and a0.ndim > 1:
            arrs[1] = _align_1d_to_nd(a1, a0)
            return np.broadcast_arrays(a0, arrs[1])

    # Fallback: pad all inputs with leading singleton dims
    max_ndim = max(a.ndim for a in arrs)
    padded = [a.reshape((1,) * (max_ndim - a.ndim) + a.shape) for a in arrs]
    try:
        return np.broadcast_arrays(*padded)
    except ValueError:
        shapes = [a.shape for a in arrs]
        raise ValueError(f"Cannot broadcast input shapes {shapes}")


def _align_1d_to_nd(one_d, nd):
    """
    Take one_d of shape (N,) and an array nd of shape (...),
    find axis i where nd.shape[i] == N, and reshape one_d to
    (1,1,...,N,...,1) so it lines up on that axis.
    """
    N = one_d.shape[0]
    for i, dim in enumerate(nd.shape):
        if dim == N:
            # build a shape of length nd.ndim: all 1s except axis i holds N
            new_shape = tuple(1 if j != i else N for j in range(nd.ndim))
            return one_d.reshape(new_shape)
    # no matching dimension found
    raise ValueError(
        f"Cannot align 1D array of length {N} with ND shape {nd.shape}"
    )


def _wrap(core, original):
    # 1) Plain Python floats
    if isinstance(original, numbers.Real) and not isinstance(
        original, (np.generic, np.ndarray)
    ):
        # core might be a 0-d array or numpy scalar
        return float(np.asarray(core).item())

    # xarray
    if isinstance(original, xr.DataArray):
        return xr.DataArray(
            core,
            coords=original.coords,
            dims=original.dims,
            name=original.name,
            attrs=original.attrs,
        )

    # pandas Series
    if isinstance(original, pd.Series):
        return pd.Series(core, index=original.index, name=original.name)

    # pandas DataFrame
    if isinstance(original, pd.DataFrame):
        return pd.DataFrame(
            core, index=original.index, columns=original.columns
        )

    # fallback: NumPy arrays/scalars
    return core


def _vectorized_wrapper(core_func, *args):
    """Generic vectorized wrapper for functions with N inputs."""
    # Extract raw arrays or scalars
    arrs = [_extract(a) for a in args]
    core = core_func(*arrs)

    # If *all* inputs were plain Python real numbers, return a Python float
    first = args[0]
    if isinstance(first, numbers.Real) and not isinstance(
        first, (np.generic, np.ndarray)
    ):
        return float(np.asarray(core).item())

    # xarray: if every arg was a DataArray, wrap back to DataArray
    if isinstance(first, xr.DataArray) and all(
        isinstance(a, xr.DataArray) for a in args
    ):
        return xr.DataArray(core, coords=first.coords, dims=first.dims)

    # pandas Series
    if isinstance(first, pd.Series) and all(
        isinstance(a, pd.Series) for a in args
    ):
        return pd.Series(core, index=first.index, name=first.name)

    # pandas DataFrame
    if isinstance(first, pd.DataFrame) and all(
        isinstance(a, pd.DataFrame) for a in args
    ):
        return pd.DataFrame(core, index=first.index, columns=first.columns)

    # fallback: NumPy array or scalar (leave as-is)
    return core


def checkArgument(argument, **kwargs):
    if kwargs.get(argument) is not None:
        return kwargs.get(argument)
    else:
        raise TypeError("Missing " + argument + " argument")
