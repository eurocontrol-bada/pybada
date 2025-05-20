"""Atmosphere module."""

from math import sqrt, pow
import numpy as np

from pyBADA import constants as const
from pyBADA import conversions as conv


def proper_round(num, dec=0):
    # First, round the number to the specified number of decimal places
    rounded_num = round(num, dec)

    # Check if the result is an integer (no decimal part)
    if rounded_num.is_integer():
        return int(rounded_num)  # Return as an integer

    return rounded_num  # Return as a float if there is a decimal part


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

    if h < const.h_11:
        theta = 1 - const.temp_h * h / const.temp_0 + DeltaTemp / const.temp_0

    else:
        theta = (const.temp_11 + DeltaTemp) / const.temp_0

    return proper_round(theta, 10)


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

    return proper_round(delta, 10)


def sigma(theta, delta):
    """Calculates the normalized air density according to the ISA model.

    :param theta: Normalized temperature [-].
    :param delta: Normalized pressure [-].
    :type theta: float
    :type delta: float
    :returns: Normalized air density [-]. The function uses the ideal gas law
        to relate pressure, temperature, and density.
    """

    return proper_round(
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
    return proper_round(a, 10)


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

    return proper_round(tas, 10)


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

    return proper_round(M, 10)


def pressureAltitude(pressure, QNH=101325.0):
    """Calculates pressure altitude based on normalized pressure and reference
    pressure (QNH).

    :param QNH: Reference pressure in Pascals (Pa), default is standard sea
        level pressure (101325 Pa).
    :param pressure: air pressure (Pa)
    :type pressure: float
    :type QNH: float
    :returns: Pressure altitude in meters (m). The pressure altitude is
        calculated by applying the barometric formula. Below the tropopause,
        the altitude is computed using the standard temperature lapse rate.
        Above the tropopause, it applies an exponential relationship for
        altitude based on pressure ratio.
    """

    if pressure > const.p_11:
        hp = (const.temp_0 / const.temp_h) * (
            1 - pow(pressure / QNH, const.R * const.temp_h / const.g)
        )
    else:
        hp = const.h_11 + const.R * const.temp_11 / const.g * np.log(
            const.p_11 / pressure
        )

    return hp


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
