# -*- coding: utf-8 -*-
"""
pyBADA
Atmosphere module
Developped @EUROCONTROL (EIH)
2024
"""

__author__ = "Henrich Glaser-Opitz"
__copyright__ = "Copyright 2024, EUROCONTROL (EIH)"
__license__ = "BADA Eurocontrol"
__version__ = "1.0.0"
__maintainer__ = "Henrich Glaser-Opitz"
__email__ = "henrich.glaser-opitz@eurocontrol.int"
__status__ = "Development"
__docformat__ = "reStructuredText"


from math import sqrt, pow, exp, log

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
    """This function returns normalised temperature according to the ISA model

    :param h: altitude [m]
    :param DeltaTemp: deviation with respect to ISA [K]
    :type h: float
    :type DeltaTemp: float
    :returns: normalised temperature [-]

    """
    if h < const.h_11:
        theta = 1 - const.temp_h * h / const.temp_0 + DeltaTemp / const.temp_0

    else:
        theta = (const.temp_11 + DeltaTemp) / const.temp_0

    return proper_round(theta, 10)


def delta(h, DeltaTemp):
    """This function returns normalised pressure according to the ISA model

    :param h: altitude [m]
    :param DeltaTemp: deviation with respect to ISA [K]
    :type h: float
    :type DeltaTemp: float
    :returns: normalised pressure [-]

    """
    p = pow(
        (theta(h, DeltaTemp) - DeltaTemp / const.temp_0),
        const.g / (const.temp_h * const.R),
    )

    if h <= const.h_11:
        delta = p
    else:
        delta = p * exp(-const.g / const.R / const.temp_11 * (h - const.h_11))

    return proper_round(delta, 10)


def sigma(theta, delta):
    """This function returns normalised air denstity according to the ISA model

    :param theta: normalised temperature according to the ISA model [-]
    :param delta: normalised pressure according to the ISA model [-]
    :type theta: float
    :type delta: float
    :returns: normalised air density [-]

    """

    return proper_round(
        ((delta * const.p_0) / (theta * const.temp_0 * const.R)) / const.rho_0, 10
    )


def aSound(theta):
    """This function calculates the speed of sound

    :param theta: normalised air temperature [-]
    :type theta: float
    :returns: speed of sound [m s^-1]

    """

    a = sqrt(const.Agamma * const.R * theta * const.temp_0)
    return proper_round(a, 10)


def mach2Tas(Mach, theta):
    """This function converts Mach number to true airspeed

    :param Mach: Mach number [-]
    :param theta: normalised air temperature [-]
    :type theta: float
    :type Mach: float
    :returns: true airspeed [m s^-1]

    """
    if Mach == float("inf"):
        tas = float("inf")
    elif Mach == float("-inf"):
        tas = float("-inf")
    else:
        tas = Mach * aSound(theta)

    return tas


def tas2Mach(v, theta):
    """This function converts true airspeed to Mach

    :param v: true airspeed [m s^-1]
    :param theta: normalised air temperature [-]
    :type v: float
    :type theta: float
    :returns: Mach number [-]

    """
    return v / aSound(theta)


def tas2Cas(tas, delta, sigma):
    """This function converts true airspeed to callibrated airspeed

    :param tas: callibrated airspeed [m s^-1]
    :param sigma: normalised air density [-]
    :param delta: normalised air pressure [-]
    :type sigma: float
    :type delta: float
    :type tas: float
    :returns: callibrated airspeed [m s^-1]

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
    """This function converts callibrated airspeed to true airspeed

    :param cas: callibrated airspeed [m s^-1]
    :param sigma: normalised air density [-]
    :param delta: normalised air pressure [-]
    :type delta: float
    :type sigma: float
    :type cas: float
    :returns: true airspeed [m s^-1]

    """
    rho = sigma * const.rho_0
    p = delta * const.p_0

    A = (
        pow(1 + const.Amu * const.rho_0 * cas * cas / (2 * const.p_0), 1 / const.Amu)
        - 1
    )
    B = pow(1 + (1 / delta) * A, const.Amu) - 1
    tas = sqrt(2 * p * B / (const.Amu * rho))

    return proper_round(tas, 10)


def mach2Cas(Mach, theta, delta, sigma):
    """This function converts Mach to callibrated airspeed

    :param Mach: Mach number [-]
    :param theta: normalised air temperature [-]
    :param delta: normalised air pressure [-]
    :param sigma: normalised air density [-]
    :type Mach: float
    :type theta: float
    :type delta: float
    :type sigma: float
    :type h: float
    :type DeltaTemp: float
    :returns: true airspeed [m s^-1]

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
    """This function converts callibrated airspeed to Mach

    :param cas: callibrated airspeed [m s^-1]
    :param theta: normalised air temperature [-]
    :param delta: normalised air pressure [-]
    :param sigma: normalised air density [-]
    :type cas: float
    :type theta: float
    :type delta: float
    :type sigma: float
    :returns: true airspeed [m s^-1]

    """
    tas = cas2Tas(cas, delta, sigma)
    M = tas2Mach(tas, theta)

    return proper_round(M, 10)


def hp(delta, QNH=101325.0):
    """This function calculates pressure altitude

    :param QNH: reference pressure [Pa]
    :param delta: normalised air pressure [-]
    :type delta: float
    :type QNH: float
    :returns: pressure altitude [m]

    """
    if delta * const.p_0 > const.p_11:
        hp = (const.temp_0 / const.temp_h) * (
            1 - pow(delta * const.p_0 / QNH, const.R * const.temp_h / const.g)
        )
    else:
        hp = const.h_11 - const.R * const.temp_11 / const.g * log(
            delta * const.p_0 / const.p_11
        )

    return hp


def crossOver(cas, Mach):
    """This function calculates cross-over altitude

    :param cas: callibrated airspeed [m s^-1]
    :param Mach: Mach number [-]
    :type cas: float
    :type Mach: float
    :returns: cross-over altitude [m]

    """
    p_trans = const.p_0 * (
        (
            pow(
                1 + ((const.Agamma - 1.0) / 2.0) * ((cas / const.a_0) ** 2),
                pow(const.Amu, -1),
            )
            - 1.0
        )
        / (pow(1 + ((const.Agamma - 1.0) / 2.0) * (Mach**2), pow(const.Amu, -1)) - 1.0)
    )

    theta_trans = pow(p_trans / const.p_0, (const.temp_h * const.R) / const.g)

    if p_trans < const.p_11:
        crossover = const.h_11 - (const.R * const.temp_11 / const.g) * log(
            p_trans / const.p_11
        )
    else:
        crossover = (const.temp_0 / -const.temp_h) * (theta_trans - 1)

    return crossover


def atmosphereProperties(h, DeltaTemp):
    """This function calculates the atmosphere properties in form of a density, temperature and pressure ratio
    based on altitude and deviation from ISA temperature

    :param h: altitude [m]
    :param DeltaTemp: deviation with respect to ISA [K]
    :type h: float
    :type DeltaTemp: float
    :returns: normalised temperature, pressure and density [-]
    """

    theta_norm = theta(h=h, DeltaTemp=DeltaTemp)
    delta_norm = delta(h=h, DeltaTemp=DeltaTemp)
    sigma_norm = sigma(theta=theta_norm, delta=delta_norm)

    return [theta_norm, delta_norm, sigma_norm]


def convertSpeed(v, speedType, theta, delta, sigma):
    """This function calculates the M, TAS and CAS speed based on imput speed and its type

    :param v: airspeed {M,CAS,TAS}[-,kt,kt]
    :param speedType: type of speed as input {M,CAS,TAS}
    :param theta: normalised air temperature [-]
    :param delta: normalised air pressure [-]
    :param sigma: normalised air density [-]
    :type v: float
    :type speedType: string
    :type theta: float
    :type delta: float
    :type sigma: float
    :returns: [M, CAS, TAS] [-, m/s, m/s]
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
