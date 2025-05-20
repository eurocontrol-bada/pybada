"""Generic airplane/helicopter performance module."""

from math import sqrt, pow, cos, radians, atan, tan, degrees

from pyBADA import constants as const
from pyBADA import conversions as conv
from pyBADA import atmosphere as atm


def checkArgument(argument, **kwargs):
    if kwargs.get(argument) is not None:
        return kwargs.get(argument)
    else:
        raise TypeError("Missing " + argument + " argument")


class Bada:
    """This class implements the mechanisms applicable across all BADA
    families."""

    def __init__(self):
        pass

    @staticmethod
    def getBADAParameters(df, acName, parameters):
        """Retrieves specified parameters for a given aircraft name from a
        DataFrame.

        :param df: DataFrame containing BADA aircraft data.
        :param acName: Name of the aircraft or list of aircraft names to
            search for.
        :param parameters: List of column names (or a single column name) to
            retrieve.
        :type df: pd.DataFrame
        :type acName: list or str
        :type parameters: list or str
        :returns: A DataFrame containing the specified parameters for the
            given aircraft.
        :rtype: pd.DataFrame
        :raises ValueError: If any of the specified columns or aircraft names
            are not found.
        """

        # Ensure parameters is a list
        if isinstance(parameters, str):
            parameters = [parameters]

        # Ensure acName is a list
        if isinstance(acName, str):
            acName = [acName]

        # Ensure all requested parameters exist in the DataFrame
        missing_cols = [col for col in parameters if col not in df.columns]
        if missing_cols:
            raise ValueError(
                f"The following parameters are not in the DataFrame columns: {missing_cols}"
            )

        # Filter rows where 'acName' matches any of the specified aircraft names
        filtered_df = df[df["acName"].isin(acName)]

        # Check if any rows were found
        if filtered_df.empty:
            raise ValueError(f"No entries found for aircraft(s): {acName}.")
        else:
            # Select the required columns
            result_df = filtered_df[["acName"] + parameters].reset_index(
                drop=True
            )
            return result_df

    @staticmethod
    def loadFactor(fi):
        """Computes the load factor from a given bank angle.

        The load factor is calculated based on the cosine of the bank angle,
        which is expressed in degrees. A small rounding operation is applied
        to avoid precision issues with small decimal places.

        :param fi: Bank angle in degrees.
        :type fi: float
        :returns: The load factor (dimensionless).
        :rtype: float
        """

        return 1 / round(cos(radians(fi)), 10)

    @staticmethod
    def bankAngle(rateOfTurn, v):
        """Computes the bank angle based on true airspeed (TAS) and rate of
        turn.

        :param v: True airspeed (TAS) in meters per second (m/s).
        :param rateOfTurn: Rate of turn in degrees per second (deg/s).
        :type v: float
        :type rateOfTurn: float
        :returns: Bank angle in degrees.
        :rtype: float
        """

        ROT = conv.deg2rad(rateOfTurn)

        BA = atan((ROT * v) / const.g)
        return conv.rad2deg(BA)

    @staticmethod
    def rateOfTurn(v, nz=1.0):
        """Computes the rate of turn based on true airspeed (TAS) and load
        factor.

        :param v: True airspeed (TAS) in meters per second (m/s).
        :param nz: Load factor (default is 1.0), dimensionless.
        :type v: float
        :type nz: float
        :returns: Rate of turn in degrees per second (deg/s).
        :rtype: float
        """

        return degrees((const.g / v) * sqrt(nz * nz - 1))

    @staticmethod
    def rateOfTurn_bankAngle(TAS, bankAngle):
        """Computes the rate of turn based on true airspeed (TAS) and bank
        angle.

        :param TAS: True airspeed (TAS) in meters per second (m/s).
        :param bankAngle: Bank angle in degrees.
        :type TAS: float
        :type bankAngle: float
        :returns: Rate of turn in degrees per second (deg/s).
        :rtype: float
        """

        ROT = tan(radians(bankAngle)) * const.g / TAS

        return degrees(ROT)

    @staticmethod
    def turnRadius(v, nz=1.0):
        """Computes the turn radius based on true airspeed (TAS) and load
        factor.

        :param v: True airspeed (TAS) in meters per second (m/s).
        :param nz: Load factor (default is 1.0), dimensionless.
        :type v: float
        :type nz: float
        :returns: Turn radius in meters.
        :rtype: float
        """

        return (v * v / const.g) * (1 / sqrt(nz * nz - 1))

    @staticmethod
    def turnRadius_bankAngle(v, ba):
        """Computes the turn radius based on true airspeed (TAS) and bank
        angle.

        :param v: True airspeed (TAS) in meters per second (m/s).
        :param ba: Bank angle in degrees.
        :type v: float
        :type ba: float
        :returns: Turn radius in meters.
        :rtype: float
        """

        return (v * v / const.g) * (1 / tan(conv.deg2rad(ba)))

    @staticmethod
    def GS(tas, gamma, Ws):
        """Computes the ground speed based on true airspeed (TAS), flight path
        angle, and wind speed.

        :param tas: True airspeed (TAS) in meters per second (m/s).
        :param gamma: Flight path angle in degrees.
        :param Ws: Longitudinal wind speed in meters per second (m/s).
        :type tas: float
        :type gamma: float
        :type Ws: float
        :returns: Ground speed in meters per second (m/s).
        :rtype: float
        """

        return tas * cos(radians(gamma)) + Ws


class BadaFamily:
    """This class sets the token for the respected BADA Family."""

    def __init__(self, BADA3=False, BADA4=False, BADAH=False, BADAE=False):
        self.BADA3 = BADA3
        self.BADA4 = BADA4
        self.BADAH = BADAH
        self.BADAE = BADAE


class Airplane:
    """This is a generic airplane class based on a three-degrees-of-freedom
    point mass model (where all the forces are applied at the center of
    gravity).

    .. note::this generic class only implements basic aircraft dynamics
            calculations, aircraft performance and optimisation can be obtained
            from its inherited classes
    """

    def __init__(self):
        pass

    @staticmethod
    def esf(**kwargs):
        """Computes the energy share factor based on flight conditions.

        :param h: Altitude in meters.
        :param DeltaTemp: Temperature deviation with respect to ISA in Kelvin.
        :param flightEvolution: Type of flight evolution
            [constM/constCAS/acc/dec].
        :param phase: Phase of flight [cl/des].
        :param v: Constant speed (Mach number).
        :type h: float
        :type DeltaTemp: float
        :type flightEvolution: str
        :type phase: str
        :type v: float
        :returns: Energy share factor (dimensionless).
        :rtype: float
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
                    + (
                        const.Agamma
                        * const.R
                        * (-const.temp_h)
                        * M
                        * M
                        / (2 * const.g)
                    )
                    * ((temp - DeltaTemp) / temp)
                )

            # constant CAS below or at tropopause
            elif flightEvolution == "constCAS" and h <= const.h_11:
                M = checkArgument("M", **kwargs)
                DeltaTemp = checkArgument("DeltaTemp", **kwargs)

                temp = atm.theta(h, DeltaTemp) * const.temp_0
                A = (
                    const.Agamma
                    * const.R
                    * (-const.temp_h)
                    * M
                    * M
                    / (2 * const.g)
                ) * ((temp - DeltaTemp) / temp)
                B = pow(
                    1 + (const.Agamma - 1) * M * M / 2, -1 / (const.Agamma - 1)
                )
                C = pow(1 + (const.Agamma - 1) * M * M / 2, 1 / const.Amu) - 1
                ESF = 1 / (1 + A + B * C)

            # constant CAS above tropopause
            elif flightEvolution == "constCAS" and h > const.h_11:
                M = checkArgument("M", **kwargs)

                ESF = 1 / (
                    1
                    + (
                        pow(
                            1 + (const.Agamma - 1) * M * M / 2,
                            -1 / (const.Agamma - 1),
                        )
                    )
                    * (
                        pow(1 + (const.Agamma - 1) * M * M / 2, 1 / const.Amu)
                        - 1
                    )
                )

            # contant TAS
            elif flightEvolution == "constTAS":
                ESF = 1

            else:
                ESF = float("Nan")

        return ESF


class Helicopter:
    """This is a generic helicopter class based on a Total-Energy Model (TEM)

    .. note::this generic class only implements basic aircraft dynamics
            calculations, aircraft performance and optimisation can be obtained
            from its inherited classes
    """

    def __init__(self):
        pass

    @staticmethod
    def esf(**kwargs):
        """Computes the energy share factor based on flight conditions.

        :param h: Altitude in meters.
        :param DeltaTemp: Temperature deviation with respect to ISA in Kelvin.
        :param flightEvolution: Type of flight evolution
            [constTAS/constCAS/acc/dec].
        :param phase: Phase of flight [Climb/Descent].
        :param v: Constant speed (Mach number).
        :type h: float
        :type DeltaTemp: float
        :type flightEvolution: str
        :type phase: str
        :type v: float
        :returns: Energy share factor (dimensionless).
        :rtype: float
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
                    const.Agamma
                    * const.R
                    * (-const.temp_h)
                    * M
                    * M
                    / (2 * const.g)
                ) * ((temp - DeltaTemp) / temp)
                B = pow(
                    1 + (const.Agamma - 1) * M * M / 2, -1 / (const.Agamma - 1)
                )
                C = pow(1 + (const.Agamma - 1) * M * M / 2, 1 / const.Amu) - 1
                ESF = 1 / (1 + A + B * C)

            # contant TAS
            elif flightEvolution == "constTAS":
                ESF = 1

            else:
                ESF = float("Nan")

        return ESF
