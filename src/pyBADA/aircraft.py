# -*- coding: utf-8 -*-
"""
pyBADA
Generic airplane/helicopter performance module
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


import abc
from math import sqrt, pow, cos, radians, atan, tan, degrees

from pyBADA import constants as const
from pyBADA import conversions as conv
from pyBADA import atmosphere as atm


def checkArgument(argument, **kwargs):
    if kwargs.get(argument) is not None:
        return kwargs.get(argument)
    else:
        raise TypeError("Missing " + argument + " argument")


class BadaFamily(object):
    """This class sets the token for the respected BADA Family."""

    def __init__(self, BADA3=False, BADA4=False, BADAH=False, BADAE=False):
        self.BADA3 = BADA3
        self.BADA4 = BADA4
        self.BADAH = BADAH
        self.BADAE = BADAE


class Airplane(object):
    """This is a generic airplane class based on a three-degrees-of-freedom point mass model (where all the forces
    are applied at the center of gravity).

    .. note::this generic class only implements basic aircraft dynamics
            calculations, aircraft performance and optimisation can be obtained
            from its inherited classes

    """

    __metaclass__ = abc.ABCMeta

    def __init__(self):
        pass

    @staticmethod
    def loadFactor(fi):
        """This function computes the load factor from bank angle

        :param fi: bank angle [deg].
        :type fi: float.
        :returns: load factor [-].
        :rtype: float.

        """
        # rounding implemented to try to minimize errors on small decimal places, like 1.9999999999999 instead of 2.0
        return 1 / round(cos(radians(fi)), 10)

    @staticmethod
    def bankAngle(rateOfTurn, v):
        """This function computes bank angle based on TAS and rate of turn

        :param v: true airspeed TAS [m s^-1].
        :param rateOfTurn: rateOfTurn [deg/s].
        :type v: float.
        :type rateOfTurn: float.
        :returns: bank angle [deg].
        :rtype: float

        """

        ROT = conv.deg2rad(rateOfTurn)

        BA = atan((ROT * v) / const.g)
        return conv.rad2deg(BA)

    @staticmethod
    def rateOfTurn_bankAngle(TAS, bankAngle):
        """This function computes the rate of turn

        :param TAS: true airspeed TAS [m s^-1].
        :param bankAngle: bank angle [deg].
        :type TAS: float.
        :type bankAngle: float.
        :returns: rate of turn [deg s^-1].
        :rtype: float.

        """

        ROT = tan(radians(bankAngle)) * const.g / TAS

        return degrees(ROT)

    @staticmethod
    def rateOfTurn(v, nz=1.0):
        """This function computes the rate of turn

        :param v: true airspeed TAS [m s^-1].
        :param nz: load factor [-].
        :type v: float.
        :type nz: float.
        :returns: rate of turn [rad s^-1].
        :rtype: float.

        """
        return degrees((const.g / v) * sqrt(nz * nz - 1))

    @staticmethod
    def turnRadius(v, nz=1.0):
        """This function computes the turn radius from load factor and speed

        :param v: true airspeed TAS [m s^-1].
        :param nz: load factor [-].
        :type v: float.
        :type nz: float.
        :returns: turn radius [m].
        :rtype: float.

        """
        return (v * v / const.g) * (1 / sqrt(nz * nz - 1))

    @staticmethod
    def turnRadius_bankAngle(v, ba):
        """This function computes the turn radius from bank angle and speed

        :param v: true airspeed [m s^-1].
        :param ba: bank angle [deg].
        :type v: float.
        :type ba: float.
        :returns: turn radius [m].
        :rtype: float.

        """

        return (v * v / const.g) * (1 / tan(conv.deg2rad(ba)))

    @staticmethod
    def GS(tas, gamma, Ws):
        """This function computes the ground speed

        :param tas: true airspeed [m s^-1].
        :param gamma: flight path angle [deg].
        :param Ws: longitudinal wind speed [m s^-1].
        :type tas: float.
        :type gamma: float.
        :type Ws: float.
        :returns: ground speed [m s^-1].
        :rtype: float.

        """
        return tas * cos(radians(gamma)) + Ws

    @staticmethod
    def esf(**kwargs):
        """This function computes the energy share factor

        :param h: altitude [m]
        :param DeltaTemp: deviation with respect to ISA [K]
        :param flightEvolution: character of the flight evolution [constM/constCAS/acc/dec][-]
        :param phase: phase of flight [cl/des][-]
        :param v: constant speed [M][-]
        :type h: float.
        :type DeltaTemp: float.
        :type flightEvolution: str.
        :type phase: str.
        :type v: float.
        :returns: energy share factor [-].
        :rtype: float.
        """

        flightEvolution = checkArgument("flightEvolution", **kwargs)

        if flightEvolution == "acc" or flightEvolution == "dec":
            phase = checkArgument("phase", **kwargs)
            # acceleration in climb or deceleration in descent
            if (flightEvolution == "acc" and phase == "cl") or (
                flightEvolution == "dec" and phase == "des"
            ):
                ESF = 0.3
            # deceleration in climb or acceleration in descent
            elif (flightEvolution == "dec" and phase == "cl") or (
                flightEvolution == "acc" and phase == "des"
            ):
                ESF = 1.7
            else:
                ESF = float("Nan")
        else:
            h = checkArgument("h", **kwargs)

            # constant M above tropopause
            if flightEvolution == "constM" and h > const.h_11:
                ESF = 1

            # constant M below or at tropopause
            elif flightEvolution == "constM" and h <= const.h_11:
                M = checkArgument("M", **kwargs)
                DeltaTemp = checkArgument("DeltaTemp", **kwargs)

                temp = atm.theta(h, DeltaTemp) * const.temp_0
                ESF = 1 / (
                    1
                    + (const.Agamma * const.R * (-const.temp_h) * M * M / (2 * const.g))
                    * ((temp - DeltaTemp) / temp)
                )

            # constant CAS below or at tropopause
            elif flightEvolution == "constCAS" and h <= const.h_11:
                M = checkArgument("M", **kwargs)
                DeltaTemp = checkArgument("DeltaTemp", **kwargs)

                temp = atm.theta(h, DeltaTemp) * const.temp_0
                A = (
                    const.Agamma * const.R * (-const.temp_h) * M * M / (2 * const.g)
                ) * ((temp - DeltaTemp) / temp)
                B = pow(1 + (const.Agamma - 1) * M * M / 2, -1 / (const.Agamma - 1))
                C = pow(1 + (const.Agamma - 1) * M * M / 2, 1 / const.Amu) - 1
                ESF = 1 / (1 + A + B * C)

            # constant CAS above tropopause
            elif flightEvolution == "constCAS" and h > const.h_11:
                M = checkArgument("M", **kwargs)

                ESF = 1 / (
                    1
                    + (pow(1 + (const.Agamma - 1) * M * M / 2, -1 / (const.Agamma - 1)))
                    * (pow(1 + (const.Agamma - 1) * M * M / 2, 1 / const.Amu) - 1)
                )

            # contant TAS
            elif flightEvolution == "constTAS":
                ESF = 1

            else:
                ESF = float("Nan")

        return ESF


class Helicopter(object):
    """This is a generic helicopter class based on a Total-Energy Model (TEM)

    .. note::this generic class only implements basic aircraft dynamics
            calculations, aircraft performance and optimisation can be obtained
            from its inherited classes

    """

    __metaclass__ = abc.ABCMeta

    def __init__(self):
        pass

    @staticmethod
    def loadFactor(fi):
        """This function computes the load factor from bank angle

        :param fi: bank angle [deg].
        :type fi: float.
        :returns: load factor [-].
        :rtype: float.

        """
        return 1 / round(cos(radians(fi)), 10)

    @staticmethod
    def rateOfTurn(v, nz=1.0):
        """This function computes the rate of turn

        :param v: true airspeed TAS [m s^-1].
        :param nz: load factor [-].
        :type v: float.
        :type nz: float.
        :returns: rate of turn [rad s^-1].
        :rtype: float.

        """
        return degrees((const.g / v) * sqrt(nz * nz - 1))

    @staticmethod
    def rateOfTurn_bankAngle(TAS, bankAngle):
        """This function computes the rate of turn

        :param TAS: true airspeed TAS [m s^-1].
        :param bankAngle: bank angle [deg].
        :type TAS: float.
        :type bankAngle: float.
        :returns: rate of turn [deg s^-1].
        :rtype: float.

        """

        ROT = tan(radians(bankAngle)) * const.g / TAS

        return degrees(ROT)

    @staticmethod
    def turnRadius(v, nz=1.0):
        """This function computes the turn radius from load factor and speed

        :param v: true airspeed TAS [m s^-1].
        :param nz: load factor [-].
        :type v: float.
        :type nz: float.
        :returns: turn radius [m].
        :rtype: float.

        """
        return (v * v / const.g) * (1 / sqrt(nz * nz - 1))

    @staticmethod
    def turnRadius_bankAngle(v, ba):
        """This function computes the turn radius from bank angle and speed

        :param v: true airspeed [m s^-1].
        :param ba: bank angle [deg].
        :type v: float.
        :type ba: float.
        :returns: turn radius [m].
        :rtype: float.

        """

        return (v * v / const.g) * (1 / tan(conv.deg2rad(ba)))

    @staticmethod
    def esf(**kwargs):
        """This function computes the energy share factor

        :param h: altitude [m]
        :param DeltaTemp: deviation with respect to ISA [K]
        :param flightEvolution: character of the flight evolution {constTAS,constCAS,acc,dec}[-]
        :param phase: phase of flight {Climb,Descent}[-]
        :param M: constant speed M [-]
        :type h: float.
        :type DeltaTemp: float.
        :type flightEvolution: str.
        :type phase: str.
        :type M: float.
        :returns: energy share factor ESF [-].
        :rtype: float.
        """

        flightEvolution = checkArgument("flightEvolution", **kwargs)

        if flightEvolution == "acc" or flightEvolution == "dec":
            phase = checkArgument("phase", **kwargs)
            # acceleration in climb or deceleration in descent
            if (flightEvolution == "acc" and phase == "Climb") or (
                flightEvolution == "dec" and phase == "Descent"
            ):
                ESF = 0.3
            # deceleration in climb or acceleration in descent
            elif (flightEvolution == "dec" and phase == "Climb") or (
                flightEvolution == "acc" and phase == "Descent"
            ):
                ESF = 1.7
            else:
                ESF = float("Nan")
        else:
            # contant CAS
            if flightEvolution == "constCAS":
                h = checkArgument("h", **kwargs)
                M = checkArgument("M", **kwargs)
                DeltaTemp = checkArgument("DeltaTemp", **kwargs)

                theta = atm.theta(h, DeltaTemp)
                temp = theta * const.temp_0

                A = (
                    const.Agamma * const.R * (-const.temp_h) * M * M / (2 * const.g)
                ) * ((temp - DeltaTemp) / temp)
                B = pow(1 + (const.Agamma - 1) * M * M / 2, -1 / (const.Agamma - 1))
                C = pow(1 + (const.Agamma - 1) * M * M / 2, 1 / const.Amu) - 1
                ESF = 1 / (1 + A + B * C)

            # contant TAS
            elif flightEvolution == "constTAS":
                ESF = 1

            else:
                ESF = float("Nan")

        return ESF

    @staticmethod
    def bankAngle(rateOfTurn, v):
        """This function computes bank angle based on TAS and rate of turn

        :param v: true airspeed TAS [m s^-1].
        :param rateOfTurn: rateOfTurn [deg/s].
        :type v: float.
        :type rateOfTurn: float.
        :returns: bank angle [deg].
        :rtype: float

        """

        ROT = conv.deg2rad(rateOfTurn)

        BA = atan((ROT * v) / const.g)
        return conv.rad2deg(BA)
