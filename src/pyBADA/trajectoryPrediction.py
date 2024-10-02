# -*- coding: utf-8 -*-
"""
pyBADA
Basic calculations for the Trajectory Prediction (TP)
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


from pyBADA import atmosphere as atm
from pyBADA import conversions as conv
from math import exp


class TrajectoryPrediction:
    """This class implements some basic calculations required for the trajectory prediction (TP)."""

    def __init__(self):
        pass

    @staticmethod
    def getInitialMass(
        AC,
        distance,
        altitude,
        M,
        payload=60,
        fuelReserve=3600,
        flightPlanInitialMass=None,
    ):
        """This function returns calculated (estimated) aircrat initial mass
                derived from the breguet leduc formula

        :param AC: aircraft {BADA3/4/H/E}
        :param distance: distance flown in [m]
        :param altitude: cruising altitude in [m]
        :param M: cruising speed in MACH [-]
        :param payload: payload mass in [%] of maximum payload
        :param fuelReserve: fuel reserve in [s]
        :param flightPlanInitialMass: initial mass defined from the flight plan [kg]
        :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}.
        :type distance: float
        :type altitude: float
        :type M: float
        :type payload: float
        :type fuelReserve: float
        :type flightPlanInitialMass: float
        :returns: initial mass [kg]
                :rtype: float
        """

        def initialMassCalculation(
            AC, distance, altitude, M, payload, fuelReserve, flightPlanInitialMass
        ):
            DeltaTemp = 0
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=altitude, DeltaTemp=DeltaTemp
            )
            TAS = atm.mach2Tas(Mach=M, theta=theta)

            config = "CR"
            flightPhase = "Cruise"
            mass = AC.MREF

            if AC.BADAFamily.BADA3:
                # compute lift coefficient
                CL = AC.CL(tas=TAS, sigma=sigma, mass=mass)
                # compute drag coefficient
                CD = AC.CD(CL=CL, config=config)
                # compute drag force
                Drag = AC.D(tas=TAS, sigma=sigma, CD=CD)
                # compute thrust force and fuel flow
                THR = Drag

                fuelFlow = AC.ff(
                    h=altitude, v=TAS, T=THR, config=config, flightPhase=flightPhase
                )

            elif AC.BADAFamily.BADA4:
                # compute lift coefficient
                CL = AC.CL(M=M, delta=delta, mass=mass)
                # compute drag coefficient
                [HLid, LG] = AC.flightEnvelope.getAeroConfig(config=config)
                CD = AC.CD(M=M, CL=CL, HLid=HLid, LG=LG)
                # compute drag force
                Drag = AC.D(M=M, delta=delta, CD=CD)
                # compute thrust force and fuel flow
                THR = Drag
                CT = AC.CT(Thrust=THR, delta=delta)

                fuelFlow = AC.ff(
                    CT=CT, delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
                )  # [kg/s]

            elif AC.BADAFamily.BADAH:
                # compute Power required for level flight
                Preq = AC.Preq(sigma=sigma, tas=TAS, mass=mass, phi=0)
                Peng_i = Preq
                # Pav_i = AC.Pav(rating="MCNT", theta=theta, delta=delta)  # assume MCNT rating as the limit

                # if Pav_i < Preq:
                # warnings.warn("Power Available is lower than Power Required",UserWarning)

                # compute fuel flow for level flight
                CP = AC.CP(Peng=Preq)
                fuelFlow = AC.ff(delta=delta, CP=CP)  # [kg/s]

            fuelReserveMass = fuelReserve * fuelFlow

            if AC.MPL is not None:
                maximumPayload = AC.MPL
            else:
                maximumPayload = AC.MTOW - AC.OEW - AC.MFL
            payloadMass = (payload / 100) * maximumPayload

            minimumLandingMass = AC.OEW + payloadMass + fuelReserveMass

            # in case of no wind, the ground speed is the same as true airspeed
            GS = TAS

            initialMass = minimumLandingMass * exp(
                (fuelFlow * distance) / (AC.MREF * GS)
            )

            # set Initial Mass from FPL check
            if flightPlanInitialMass is not None:
                initialMass = flightPlanInitialMass

            # envelope check
            initialMass = min(max(initialMass, AC.OEW), AC.MTOW)

            return initialMass

        if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
            if (AC.MMO is not None and AC.MMO >= 1.0) or (
                AC.VMO is not None and AC.VMO >= 400
            ):
                # identified as fighter jet
                initialMass = AC.MREF
            else:
                initialMass = initialMassCalculation(
                    AC=AC,
                    distance=distance,
                    altitude=altitude,
                    M=M,
                    payload=payload,
                    fuelReserve=fuelReserve,
                    flightPlanInitialMass=flightPlanInitialMass,
                )

        elif AC.BADAFamily.BADAH:
            if AC.vne is not None and AC.vne >= 400:
                # identified as fighter
                initialMass = AC.MREF
            else:
                initialMass = initialMassCalculation(
                    AC=AC,
                    distance=distance,
                    altitude=altitude,
                    M=M,
                    payload=payload,
                    fuelReserve=fuelReserve,
                    flightPlanInitialMass=flightPlanInitialMass,
                )

        return initialMass
