"""Atmosphere module."""

from math import pow, sqrt

import numpy as np
import xarray as xr
import pandas as pd

from pyBADA import utils
from pyBADA import constants as const
from pyBADA import conversions as conv


def theta(h, DeltaTemp):
    """Calculates the normalized temperature according to the International
    Standard Atmosphere (ISA) model.

    :param h: Altitude in meters (m).
    :param DeltaTemp: Deviation from ISA temperature in Kelvin (K).
    :type h: float
    :type DeltaTemp: float
    :returns: Normalized temperature [-]. The function accounts for whether
        the altitude is below or above the tropopause (11,000 m). Below the
        tropopause, it applies the temperature lapse rate. Above the
        tropopause, a constant temperature is assumed.
    """

    # xarray DataArray branch: native vectorized operations preserve metadata
    if isinstance(h, xr.DataArray):
        # align DeltaTemp to h, broadcasting if scalar or array
        if isinstance(DeltaTemp, xr.DataArray):
            dt = DeltaTemp
        else:
            # broadcast scalar or numpy array to match h
            dt = xr.full_like(h, utils.to_numpy(DeltaTemp))

        cond = h < const.h_11
        base = 1 - const.temp_h * h / const.temp_0
        tropo = const.temp_11 / const.temp_0
        out = xr.where(cond, base + dt / const.temp_0, tropo + dt / const.temp_0)
        return out.round(10)

    # pandas Series/DataFrame branch: native vectorized operations preserve index/columns
    if isinstance(h, (pd.Series, pd.DataFrame)):
        # align DeltaTemp
        if isinstance(DeltaTemp, type(h)):
            dt = DeltaTemp
        else:
            arr_dt = utils.to_numpy(DeltaTemp)
            dt = pd.Series(arr_dt, index=h.index) if isinstance(h, pd.Series) else pd.DataFrame(arr_dt, index=h.index, columns=h.columns)

        base = 1 - const.temp_h * h / const.temp_0
        tropo = const.temp_11 / const.temp_0
        result = base.where(h < const.h_11, tropo) + dt / const.temp_0
        return result.round(10)

    # NumPy array or scalar fallback: inline piecewise ISA computation
    arr_h = utils.to_numpy(h)
    arr_dT = utils.to_numpy(DeltaTemp)
    out = np.where(
        arr_h < const.h_11,
        1 - const.temp_h * arr_h / const.temp_0 + arr_dT / const.temp_0,
        (const.temp_11 + arr_dT) / const.temp_0
    )
    out = np.round(out, 10)
    return out.item() if (isinstance(out, np.ndarray) and out.ndim == 0) else out


def delta(h, DeltaTemp):
    """Calculates the normalized pressure according to the ISA model.

    :param h: Altitude in meters (m).
    :param DeltaTemp: Deviation from ISA temperature in Kelvin (K).
    :type h: float
    :type DeltaTemp: float
    :returns: Normalized pressure [-]. The function uses the barometric
        equation for pressure changes below and above the tropopause.
    """

    p = pow(
        (theta(h, DeltaTemp) - DeltaTemp / const.temp_0),
        const.g / (const.temp_h * const.R),
    )

    if h <= const.h_11:
        delta = p
    else:
        delta = p * np.exp(
            -const.g / const.R / const.temp_11 * (h - const.h_11)
        )

    return utils.proper_round(delta, 10)


def sigma(theta, delta):
    """Calculates the normalized air density according to the ISA model.

    :param theta: Normalized temperature [-].
    :param delta: Normalized pressure [-].
    :type theta: float
    :type delta: float
    :returns: Normalized air density [-]. The function uses the ideal gas law
        to relate pressure, temperature, and density.
    """

    return utils.proper_round(
        ((delta * const.p_0) / (theta * const.temp_0 * const.R)) / const.rho_0,
        10,
    )


def aSound(theta):
    """Calculates the speed of sound based on the normalized air temperature.

    :param theta: Normalized temperature [-].
    :type theta: float
    :returns: Speed of sound in meters per second (m/s). The speed of sound
        depends on air temperature and is calculated using the specific heat
        ratio and the gas constant.
    """

    a = sqrt(const.Agamma * const.R * theta * const.temp_0)
    return utils.proper_round(a, 10)


def mach2Tas(Mach, theta):
    """Converts Mach number to true airspeed (TAS).

    :param Mach: Mach number [-].
    :param theta: Normalized air temperature [-].
    :type Mach: float
    :type theta: float
    :returns: True airspeed in meters per second (m/s).
    """

    if Mach == float("inf"):
        tas = float("inf")
    elif Mach == float("-inf"):
        tas = float("-inf")
    else:
        tas = Mach * aSound(theta)

    return tas


def tas2Mach(v, theta):
    """Converts true airspeed (TAS) to Mach number.

    :param v: True airspeed in meters per second (m/s).
    :param theta: Normalized air temperature [-].
    :type v: float
    :type theta: float
    :returns: Mach number [-].
    """

    return v / aSound(theta)


def tas2Cas(tas, delta, sigma):
    """Converts true airspeed (TAS) to calibrated airspeed (CAS).

    :param tas: True airspeed in meters per second (m/s).
    :param sigma: Normalized air density [-].
    :param delta: Normalized air pressure [-].
    :type tas: float
    :type sigma: float
    :type delta: float
    :returns: Calibrated airspeed in meters per second (m/s). The function
        uses a complex formula to account for air compressibility effects at
        high speeds.
    """

    if tas == float("inf"):
        cas = float("inf")
    elif tas == float("-inf"):
        cas = float("-inf")
    else:
        rho = sigma * const.rho_0
        p = delta * const.p_0

        A = pow(1 + const.Amu * rho * tas * tas / (2 * p), 1 / const.Amu) - 1
        B = pow(1 + delta * A, const.Amu) - 1
        cas = sqrt(2 * const.p_0 * B / (const.Amu * const.rho_0))

    return cas


def cas2Tas(cas, delta, sigma):
    """Converts calibrated airspeed (CAS) to true airspeed (TAS).

    :param cas: Calibrated airspeed in meters per second (m/s).
    :param sigma: Normalized air density [-].
    :param delta: Normalized air pressure [-].
    :type cas: float
    :type delta: float
    :type sigma: float
    :returns: True airspeed in meters per second (m/s). This function inverts
        the compressibility adjustments to compute TAS from CAS.
    """

    rho = sigma * const.rho_0
    p = delta * const.p_0

    A = (
        pow(
            1 + const.Amu * const.rho_0 * cas * cas / (2 * const.p_0),
            1 / const.Amu,
        )
        - 1
    )
    B = pow(1 + (1 / delta) * A, const.Amu) - 1
    tas = sqrt(2 * p * B / (const.Amu * rho))

    return utils.proper_round(tas, 10)


def mach2Cas(Mach, theta, delta, sigma):
    """Converts Mach number to calibrated airspeed (CAS).

    :param Mach: Mach number [-].
    :param theta: Normalized air temperature [-].
    :param delta: Normalized air pressure [-].
    :param sigma: Normalized air density [-].
    :type Mach: float
    :type theta: float
    :type delta: float
    :type sigma: float
    :returns: Calibrated airspeed in meters per second (m/s).
    """

    if Mach == float("inf"):
        cas = float("inf")
    elif Mach == float("-inf"):
        cas = float("-inf")
    else:
        tas = mach2Tas(Mach=Mach, theta=theta)
        cas = tas2Cas(tas=tas, delta=delta, sigma=sigma)

    return cas


def cas2Mach(cas, theta, delta, sigma):
    """Converts calibrated airspeed (CAS) to Mach number.

    :param cas: Calibrated airspeed in meters per second (m/s).
    :param theta: Normalized air temperature [-].
    :param delta: Normalized air pressure [-].
    :param sigma: Normalized air density [-].
    :type cas: float
    :type theta: float
    :type delta: float
    :type sigma: float
    :returns: Mach number [-].
    """

    tas = cas2Tas(cas, delta, sigma)
    M = tas2Mach(tas, theta)

    return utils.proper_round(M, 10)


def pressureAltitude(pressure, QNH=101325.0):
    """
    Calculates pressure altitude based on air pressure and reference pressure (QNH).

    :param pressure: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Air pressure in Pascals (Pa).
    :param QNH: float
        Reference sea level pressure in Pascals (Pa). Default is 101325 Pa.
    :returns: float, np.ndarray, pandas.Series, pandas.DataFrame, or xr.DataArray
        Pressure altitude in meters (m), matching the type and metadata of `pressure`.
    """
    # xarray DataArray branch: preserves coords and attrs
    if isinstance(pressure, xr.DataArray):
        p = pressure
        branch1 = (const.temp_0/const.temp_h) * (
            1 - (p/QNH)**(const.R*const.temp_h/const.g)
        )
        branch2 = const.h_11 + (const.R*const.temp_11/const.g) * xr.ufuncs.log(
            const.p_11/p
        )
        out = xr.where(p > const.p_11, branch1, branch2)
        out.attrs.update(units="m", long_name="pressure altitude")
        return out.round(2)

    # pandas Series/DataFrame branch: preserves index/columns
    if isinstance(pressure, (pd.Series, pd.DataFrame)):
        p = pressure.astype(float)
        part1 = (const.temp_0/const.temp_h) * (
            1 - (p/QNH)**(const.R*const.temp_h/const.g)
        )
        part2 = const.h_11 + (const.R*const.temp_11/const.g) * np.log(
            const.p_11/p
        )
        arr_out = np.where(p > const.p_11, part1, part2)
        if isinstance(p, pd.Series):
            out = pd.Series(arr_out, index=p.index, name=p.name)
        else:
            out = pd.DataFrame(arr_out, index=p.index, columns=p.columns)
        return out.round(2)

    # NumPy array or scalar fallback: preserves shape and returns scalar if needed
    arr = np.asarray(pressure, dtype=float)
    part1 = (const.temp_0/const.temp_h) * (
        1 - (arr/QNH)**(const.R*const.temp_h/const.g)
    )
    part2 = const.h_11 + (const.R*const.temp_11/const.g) * np.log(
        const.p_11/arr
    )
    out = np.where(arr > const.p_11, part1, part2)
    out = np.round(out, 2)
    return out.item() if (isinstance(out, np.ndarray) and out.ndim == 0) else out


def ISATemperatureDeviation(temperature, pressureAltitude):
    """Calculates deviation from ISA temperature at a specific pressure
    altitude.

    :param temperature: air temperature (Kelvin)
    :param pressureAltitude: pressure altitude (m)
    :type temperature: float
    :type pressureAltitude: float
    :returns: ISA temperature deviation (Kelvin).
    """

    stdTemperature = theta(h=pressureAltitude, DeltaTemp=0) * const.temp_0
    deltaISATemp = temperature - stdTemperature

    return deltaISATemp


def crossOver(cas, Mach):
    """Calculates the cross-over altitude where calibrated airspeed (CAS) and
    Mach number intersect.

    :param cas: Calibrated airspeed in meters per second (m/s).
    :param Mach: Mach number [-].
    :type cas: float
    :type Mach: float
    :returns: Cross-over altitude in meters (m). The cross-over altitude is
        where CAS and Mach produce the same true airspeed. The function
        calculates pressure and temperature at this altitude based on the
        given Mach number and CAS.
    """

    p_trans = const.p_0 * (
        (
            pow(
                1 + ((const.Agamma - 1.0) / 2.0) * ((cas / const.a_0) ** 2),
                pow(const.Amu, -1),
            )
            - 1.0
        )
        / (
            pow(
                1 + ((const.Agamma - 1.0) / 2.0) * (Mach**2),
                pow(const.Amu, -1),
            )
            - 1.0
        )
    )

    theta_trans = pow(p_trans / const.p_0, (const.temp_h * const.R) / const.g)

    if p_trans < const.p_11:
        crossover = const.h_11 - (const.R * const.temp_11 / const.g) * np.log(
            p_trans / const.p_11
        )
    else:
        crossover = (const.temp_0 / -const.temp_h) * (theta_trans - 1)

    return crossover


def atmosphereProperties(h, DeltaTemp):
    """
    Calculates atmospheric properties: normalized temperature, pressure, and density ratios based on altitude and temperature deviation from ISA.

    :param h: Altitude in meters (m).
    :param DeltaTemp: Deviation from ISA temperature in Kelvin (K).
    :type h: float
    :type DeltaTemp: float
    :returns: Normalized temperature, pressure, and density ratios as a list [-].
    """

    theta_norm = theta(h=h, DeltaTemp=DeltaTemp)
    delta_norm = delta(h=h, DeltaTemp=DeltaTemp)
    sigma_norm = sigma(theta=theta_norm, delta=delta_norm)

    return [theta_norm, delta_norm, sigma_norm]


def convertSpeed(v, speedType, theta, delta, sigma):
    """Calculates Mach, true airspeed (TAS), and calibrated airspeed (CAS)
    based on input speed and its type.

    :param v: Airspeed value, depending on the type provided (M, CAS, TAS) [-,
        kt, kt].
    :param speedType: Type of input speed, which can be one of "M" (Mach),
        "CAS", or "TAS".
    :param theta: Normalized air temperature [-].
    :param delta: Normalized air pressure [-].
    :param sigma: Normalized air density [-].
    :type v: float
    :type speedType: string
    :type theta: float
    :type delta: float
    :type sigma: float
    :returns: A list of [Mach number, CAS in m/s, TAS in m/s].
    """

    if speedType == "TAS":
        TAS = conv.kt2ms(v)
        CAS = tas2Cas(tas=TAS, delta=delta, sigma=sigma)
        M = tas2Mach(v=TAS, theta=theta)

    elif speedType == "CAS":
        CAS = conv.kt2ms(v)
        TAS = cas2Tas(cas=CAS, delta=delta, sigma=sigma)
        M = tas2Mach(v=TAS, theta=theta)

    elif speedType == "M":
        M = v
        CAS = mach2Cas(Mach=M, theta=theta, delta=delta, sigma=sigma)
        TAS = cas2Tas(cas=CAS, delta=delta, sigma=sigma)
    else:
        raise Exception("Expected TAS, CAS or M, received: " + speedType)

    return [M, CAS, TAS]
