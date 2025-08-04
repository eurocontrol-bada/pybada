"""Atmosphere module."""

from math import pow, sqrt

import numpy as np
import xarray as xr
import pandas as pd

from pyBADA import utils
from pyBADA import constants as const
from pyBADA import conversions as conv


def _theta_core(arr_h, arr_dT):
    """
    Core normalized temperature ratio computation
    with constant temperature above the tropopause.
    """

    return np.where(
        arr_h < const.h_11,
        1 - const.temp_h * arr_h / const.temp_0 + arr_dT / const.temp_0,
        (const.temp_11 + arr_dT) / const.temp_0
    )


def theta(h, DeltaTemp):
    """
    Normalized temperature according to the ISA model, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param h: Altitude in meters (m). Can be scalar, numpy array,
              pandas Series/DataFrame, or xarray.DataArray.
    :param DeltaTemp: Deviation from ISA temperature in Kelvin (K). Same type as h.
    :type h: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :type DeltaTemp: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Normalized temperature [-], rounded with half-up.
    """

    arr_h  = utils._extract(h)
    arr_dT = utils._extract(DeltaTemp)

    core = _theta_core(arr_h, arr_dT)

    return utils._wrap(core)


def _delta_core(arr_h, arr_dT):
    """
    Core normalized pressure computation using numpy arrays.
    Applies the ISA model with tropopause correction.
    Requires arr_h and arr_dT; computes arr_theta internally.

    :param arr_h: Altitude array (np.ndarray or scalar).
    :param arr_dT: ISA temperature deviation array (np.ndarray or scalar).
    :returns: Normalized pressure ratio array (np.ndarray or scalar).
    """

    arr_theta = _theta_core(arr_h, arr_dT)

    exponent = const.g / (const.temp_h * const.R)
    base = arr_theta - (arr_dT / const.temp_0)
    p = np.power(base, exponent)

    strato_factor = np.exp(
        - const.g / (const.R * const.temp_11) * (arr_h - const.h_11)
    )

    return np.where(
        arr_h <= const.h_11,
        p,
        p * strato_factor
    )


def delta(h, DeltaTemp):
    """
    Normalized pressure according to the ISA model, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param h: Altitude in meters (m). Can be scalar, numpy array,
              pandas Series/DataFrame, or xarray.DataArray.
    :param DeltaTemp: Deviation from ISA temperature in Kelvin (K). Same type as h.
    :type h: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :type DeltaTemp: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Normalized pressure [-], rounded with half-up.
    """

    arr_h  = utils._extract(h)
    arr_dT = utils._extract(DeltaTemp)

    core = _delta_core(arr_h, arr_dT)

    return utils._wrap(core)


def _sigma_core(arr_theta, arr_delta):
    """
    Core normalized density computation from precomputed arrays
    """

    # Apply ideal gas law to normalized variables
    return (arr_delta * const.p_0) / (arr_theta * const.temp_0 * const.R) / const.rho_0


def sigma(h=None, DeltaTemp=None, theta=None, delta=None):
    """Normalized density according to ISA, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    You can either provide:
      - `h` (pressure altitude) and `DeltaTemp` (temperature deviation) to compute
        normalized temperature and pressure internally;
      - Or precomputed `theta` (normalized temperature) and `delta`
        (normalized pressure) directly to avoid recomputation.

    :param h: pressure altitude AMSL [m], required if `theta` is None.
    :type h: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :param DeltaTemp: Temperature deviation from ISA at sea level [K], required if `theta` is None.
    :type DeltaTemp: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :param theta: Precomputed normalized temperature [-], optional.
    :type theta: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :param delta: Precomputed normalized pressure [-], optional.
    :type delta: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Normalized air density [-]
    """

    if theta is not None and delta is not None:
        arr_t = utils._extract(theta)
        arr_d = utils._extract(delta)
        core = _sigma_core(arr_t, arr_d)

        return utils._wrap(core)

    if h is None or DeltaTemp is None:
        raise ValueError("Either provide both h & DeltaTemp, or theta & delta.")
        arr_h  = utils._extract(h)
        arr_dT = utils._extract(DeltaTemp)
        arr_t  = _theta_core(arr_h, arr_dT)
        arr_d  = _delta_core(arr_h, arr_dT)
        core   = _sigma_core(arr_t, arr_d)

        return utils._wrap(core)


def _aSound_core(arr_theta):
    """
    Core speed-of-sound computation from normalized temperature arrays
    """
    return np.sqrt(const.Agamma * const.R * arr_theta * const.temp_0)


def aSound(theta):
    """Calculates the speed of sound based on the normalized air temperature.

    :param theta: Normalized temperature [-]. Can be scalar, numpy array,
                  pandas Series/DataFrame, or xarray.DataArray.
    :returns: Speed of sound in meters per second (m/s)
    """

    arr = utils._extract(theta)

    core = _aSound_core(arr)

    return utils._wrap(core)


def _mach2Tas_core(arr_mach, arr_theta):
    """
    Core true airspeed computation
    """

    return arr_mach * _aSound_core(arr_theta)


def _tas2Mach_core(arr_v, arr_theta):
    """
    Core Mach number computation from true airspeed and normalized temperature
    """

    return arr_v / _aSound_core(arr_theta)


def mach2Tas(Mach, theta):
    """Converts Mach number to true airspeed (TAS), vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param Mach: Mach number [-]. Can be scalar, numpy array,
                  pandas Series/DataFrame, or xarray.DataArray.
    :param theta: Normalized air temperature [-]. Same type as Mach.
    :type Mach: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :type theta: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: True airspeed in meters per second (m/s).
    """

    return utils._vectorized_wrapper(_mach2Tas_core, Mach, theta)


def tas2Mach(v, theta):
    """Converts true airspeed (TAS) to Mach number, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param v: True airspeed in meters per second (m/s). Can be scalar, numpy array,
              pandas Series/DataFrame, or xarray.DataArray.
    :param theta: Normalized air temperature [-]. Same type as v.
    :type v: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :type theta: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Mach number [-].
    """

    return utils._vectorized_wrapper(_tas2Mach_core, v, theta)


def _tas2Cas_core(arr_tas, arr_delta, arr_sigma):
    """
    Core calibrated airspeed computation:
    Uses compressibility correction formula.
    """
    rho = arr_sigma * const.rho_0
    p = arr_delta * const.p_0
    A = pow(1 + const.Amu * rho * arr_tas**2 / (2 * p), 1 / const.Amu) - 1
    B = pow(1 + arr_delta * A, const.Amu) - 1
    return sqrt(2 * const.p_0 * B / (const.Amu * const.rho_0))


def _cas2Tas_core(arr_cas, arr_delta, arr_sigma):
    """
    Core true airspeed from calibrated airspeed:
    Inverts compressibility corrections.
    """
    rho = arr_sigma * const.rho_0
    p = arr_delta * const.p_0
    A = pow(1 + const.Amu * const.rho_0 * arr_cas**2 / (2 * const.p_0), 1 / const.Amu) - 1
    B = pow(1 + (1 / arr_delta) * A, const.Amu) - 1
    return sqrt(2 * p * B / (const.Amu * rho))


def tas2Cas(tas, delta, sigma):
    """Converts true airspeed (TAS) to calibrated airspeed (CAS), vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param tas: True airspeed in meters per second (m/s). Can be scalar, numpy array,
                pandas Series/DataFrame, or xarray.DataArray.
    :param delta: Normalized air pressure [-]. Same type as tas.
    :param sigma: Normalized air density [-]. Same type as tas.
    :type tas: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :type delta: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :type sigma: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Calibrated airspeed in meters per second (m/s)
    """

    return utils._vectorized_wrapper(_tas2Cas_core, tas, delta, sigma)


def cas2Tas(cas, delta, sigma):
    """Converts calibrated airspeed (CAS) to true airspeed (TAS), vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param cas: Calibrated airspeed in meters per second (m/s). Can be scalar, numpy array,
                pandas Series/DataFrame, or xarray.DataArray.
    :param delta: Normalized air pressure [-]. Same type as cas.
    :param sigma: Normalized air density [-]. Same type as cas.
    :type cas: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :type delta: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :type sigma: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: True airspeed in meters per second (m/s).
    """

    return utils._vectorized_wrapper(_cas2Tas_core, cas, delta, sigma)


def _mach2Cas_core(arr_mach, arr_theta, arr_delta, arr_sigma):
    """
    Core conversion from Mach to calibrated airspeed:
    Compute TAS then CAS core.
    """
    tas = _mach2Tas_core(arr_mach, arr_theta)
    return _tas2Cas_core(tas, arr_delta, arr_sigma)


def _cas2Mach_core(arr_cas, arr_theta, arr_delta, arr_sigma):
    """
    Core conversion from calibrated airspeed to Mach:
    Compute TAS then Mach core.
    """
    tas = _cas2Tas_core(arr_cas, arr_delta, arr_sigma)
    return _tas2Mach_core(tas, arr_theta)


def mach2Cas(Mach, theta, delta, sigma):
    """Converts Mach number to calibrated airspeed (CAS), vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param Mach: Mach number [-]. Can be scalar, numpy array,
                  pandas Series/DataFrame, or xarray.DataArray.
    :param theta: Normalized air temperature [-]. Same type as Mach.
    :param delta: Normalized air pressure [-]. Same type as Mach.
    :param sigma: Normalized air density [-]. Same type as Mach.
    :type Mach: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :type theta: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :type delta: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :type sigma: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Calibrated airspeed in meters per second (m/s).
    """

    return utils._vectorized_wrapper(_mach2Cas_core, Mach, theta, delta, sigma)


def cas2Mach(cas, theta, delta, sigma):
    """Converts calibrated airspeed (CAS) to Mach number, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param cas: Calibrated airspeed in meters per second (m/s). Can be scalar, numpy array,
                pandas Series/DataFrame, or xarray.DataArray.
    :param theta: Normalized air temperature [-]. Same type as cas.
    :param delta: Normalized air pressure [-]. Same type as cas.
    :param sigma: Normalized air density [-]. Same type as cas.
    :type cas: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :type theta: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :type delta: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :type sigma: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Mach number [-].
    """

    return utils._vectorized_wrapper(_cas2Mach_core, cas, theta, delta, sigma)


def _pressureAltitude_core(arr_p, QNH):
    """
    Core pressure altitude computation using ISA:
    """

    part1 = (const.temp_0/const.temp_h) * (1 - np.power(arr_p/QNH, const.R*const.temp_h/const.g))
    part2 = const.h_11 + (const.R*const.temp_11/const.g) * np.log(const.p_11/arr_p)
    return np.where(arr_p > const.p_11, part1, part2)


def pressureAltitude(pressure, QNH=101325.0):
    """Calculates pressure altitude based on air pressure and reference pressure (QNH),
    vectorized for xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param pressure: Air pressure in Pascals (Pa). Can be scalar, numpy array,
                     pandas Series/DataFrame, or xarray.DataArray.
    :param QNH:     Reference sea-level pressure in Pascals (Pa). Default 101325 Pa.
    :type pressure: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :type QNH:      float
    :returns:       Pressure altitude in meters (m), rounded with half-up, preserving type and metadata.
    """

    arr = utils._extract(pressure)

    core = _pressureAltitude_core(arr, QNH)

    return utils._wrap(core)


def _isaTempDev_core(arr_temp, arr_h):
    """
    Core ISA temperature deviation computation
    """

    return arr_temp - _theta_core(arr_h, np.zeros_like(arr_h)) * const.temp_0


def ISATemperatureDeviation(temperature, pressureAltitude):
    """Calculates deviation from ISA temperature at a specific pressure altitude,
    vectorized for xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param temperature: Air temperature in Kelvin. Can be scalar, numpy array,
                        pandas Series/DataFrame, or xarray.DataArray.
    :param pressureAltitude: Pressure altitude in meters. Same type as temperature.
    :type temperature: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :type pressureAltitude: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: ISA temperature deviation in Kelvin, rounded with half-up.
    """

    arr_temp = utils._extract(temperature)
    arr_h = utils._extract(pressureAltitude)

    core = _isaTempDev_core(arr_temp, arr_h)

    return utils._wrap(core)


def _crossOver_core(arr_cas, arr_mach):
    """
    Core cross-over altitude computation:
    Calculates the altitude where CAS and Mach yield the same TAS.
    """
    # Pressure at transition
    p_trans = const.p_0 * (
        (
            pow(1 + ((const.Agamma - 1.0)/2.0) * ((arr_cas/const.a_0)**2), 1/const.Amu) - 1.0
        ) /
        (
            pow(1 + ((const.Agamma - 1.0)/2.0) * (arr_mach**2), 1/const.Amu) - 1.0
        )
    )
    # Temperature ratio at transition
    theta_trans = np.power(p_trans/const.p_0, (const.temp_h*const.R)/const.g)
    # Altitude branch
    h = np.where(
        p_trans < const.p_11,
        const.h_11 - (const.R*const.temp_11/const.g)*np.log(p_trans/const.p_11),
        (const.temp_0/-const.temp_h)*(theta_trans - 1)
    )
    return h


def crossOver(cas, Mach):
    """Calculates the cross-over altitude where calibrated airspeed (CAS) and Mach number intersect, vectorized for
    xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param cas: Calibrated airspeed in meters per second (m/s). Can be scalar, numpy array,
                pandas Series/DataFrame, or xarray.DataArray.
    :param Mach: Mach number [-]. Same type as cas.
    :type cas: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :type Mach: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: Cross-over altitude in meters (m), rounded with half-up.
    """
    return utils._vectorized_wrapper(_crossOver_core, cas, Mach)


def atmosphereProperties(h, DeltaTemp):
    """
    Calculates atmospheric properties: normalized temperature, pressure,
    and density ratios based on altitude and temperature deviation from ISA.
    Vectorized for xarray.DataArray, pandas Series/DataFrame, and numpy arrays/scalars.

    :param h: Altitude in meters (m). Can be scalar, numpy array,
              pandas Series/DataFrame, or xarray.DataArray.
    :param DeltaTemp: Deviation from ISA temperature in Kelvin (K). Same type as h.
    :type h: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :type DeltaTemp: float or array-like or xarray.DataArray or pandas Series/DataFrame
    :returns: List of [theta_norm, delta_norm, sigma_norm], each matching the type of input.
    """

    arr_h  = utils._extract(h)
    arr_dT = utils._extract(DeltaTemp)

    arr_theta = _theta_core(arr_h, arr_dT)
    arr_delta = _delta_core(arr_h, arr_dT)
    arr_sigma = _sigma_core(arr_theta, arr_delta)

    return [
        utils._wrap(arr_theta),
        utils._wrap(arr_delta),
        utils._wrap(arr_sigma),
    ]


def convertSpeed(v, speedType, theta, delta, sigma):
    """
    Calculates Mach, true airspeed (TAS), and calibrated airspeed (CAS)
    based on input speed and its type. Vectorized for xarray.DataArray,
    pandas Series/DataFrame, and numpy arrays/scalars.

    :param v: Airspeed value, depending on the type provided (M, CAS, TAS) [-, kt, kt].
              Can be scalar, numpy array, pandas Series/DataFrame,
              or xarray.DataArray.
    :param speedType: Type of input speed: "M" (Mach), "CAS", or "TAS".
    :param theta: Normalized air temperature [-]. Same type as v.
    :param delta: Normalized air pressure [-]. Same type as v.
    :param sigma: Normalized air density [-]. Same type as v.
    :returns: [Mach number, CAS (m/s), TAS (m/s)] each matching the type of input.
    """
    # Extract raw numpy arrays
    arr_v     = utils._extract(v)
    arr_theta = utils._extract(theta)
    arr_delta = utils._extract(delta)
    arr_sigma = utils._extract(sigma)

    # Compute core values on raw arrays
    if speedType.upper() == "TAS":
        arr_TAS = conv.kt2ms(arr_v)
        arr_CAS = _tas2Cas_core(arr_TAS, arr_delta, arr_sigma)
        arr_M   = _tas2Mach_core(arr_TAS, arr_theta)

    elif speedType.upper() == "CAS":
        arr_CAS = conv.kt2ms(arr_v)
        arr_TAS = _cas2Tas_core(arr_CAS, arr_delta, arr_sigma)
        arr_M   = _tas2Mach_core(arr_TAS, arr_theta)

    elif speedType.upper() == "M":
        arr_M   = arr_v
        arr_CAS = _mach2Cas_core(arr_M, arr_theta, arr_delta, arr_sigma)
        arr_TAS = _cas2Tas_core(arr_CAS, arr_delta, arr_sigma)

    else:
        raise ValueError(f"Expected speedType 'TAS', 'CAS' or 'M', got: {speedType!r}")

    return [
        utils._wrap(arr_M),
        utils._wrap(arr_CAS),
        utils._wrap(arr_TAS),
    ]