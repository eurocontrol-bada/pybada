"""
Basic calculations for the Trajectory Prediction (TP) using BADA
"""

from pyBADA import atmosphere as atm
from math import exp


def cruiseFuelConsumption(AC, altitude, M, deltaTemp):
    """
    Calculate the cruise fuel consumption for an aircraft during cruise flight using BADA.

    :param AC: Aircraft object (instance of Bada3Aircraft, Bada4Aircraft, or BadaHAircraft).
    :param altitude: Altitude in meters.
    :param M: Mach number at cruising altitude.
    :param deltaTemp: Temperature deviation from standard atmosphere.
    :type AC: object
    :type altitude: float
    :type M: float
    :type deltaTemp: float
    :return: Fuel flow in kg/s.
    :rtype: float
    """

    [theta, delta, sigma] = atm.atmosphereProperties(
        h=altitude, DeltaTemp=deltaTemp
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
            h=altitude,
            v=TAS,
            T=THR,
            config=config,
            flightPhase=flightPhase,
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
            CT=CT, delta=delta, theta=theta, M=M, deltaTemp=deltaTemp
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

    return fuelFlow


def breguetLeducInitialMass(
    AC, distance, GS, cruiseFuelFlow, payload, fuelReserve
):
    """Calculate the estimated initial mass required for the aircraft using
    the Breguet Leduc formula.

    :param AC: Aircraft object (instance of Bada3Aircraft, Bada4Aircraft, or
        BadaHAircraft).
    :param distance: Flight distance in meters.
    :param GS: Ground speed in m/s (assumed equal to true airspeed under no-
        wind conditions).
    :param cruiseFuelFlow: Fuel flow rate during cruise in kg/s.
    :param payload: Payload percentage (of the maximum payload mass) to be
        used.
    :param fuelReserve: Fuel reserve time in seconds.
    :type AC: object
    :type distance: float
    :type GS: float
    :type cruiseFuelFlow: float
    :type payload: float
    :type fuelReserve: float
    :return: Initial mass in kg.
    :rtype: float
    """

    fuelReserveMass = fuelReserve * cruiseFuelFlow

    if AC.MPL is not None:
        maximumPayload = AC.MPL
    else:
        maximumPayload = AC.MTOW - AC.OEW - AC.MFL
    payloadMass = (payload / 100) * maximumPayload

    minimumLandingMass = AC.OEW + payloadMass + fuelReserveMass

    initialMass = minimumLandingMass * exp(
        (cruiseFuelFlow * distance) / (AC.MREF * GS)
    )

    return initialMass


def getInitialMass(
    AC,
    distance,
    altitude,
    M,
    payload=60,
    fuelReserve=3600,
    flightPlanInitialMass=None,
    deltaTemp=0,
):
    """Calculates the estimated initial aircraft mass assumig cruise phase,
    combining flight plan data, aircraft envelope constraints, and an
    exponential fuel consumption model inspired by the Breguet Leduc formula.

    :param AC: Aircraft object (instance of Bada3Aircraft, Bada4Aircraft,
        BadaEAircraft, or BadaHAircraft).
    :param distance: Distance to be flown in meters.
    :param altitude: Cruising altitude in meters.
    :param M: Mach number at cruising altitude.
    :param payload: Percentage of the maximum payload mass (default is 60%).
    :param fuelReserve: Fuel reserve time in seconds (default is 3600 seconds,
        or 1 hour).
    :param flightPlanInitialMass: Optional initial mass from a flight plan, in
        kg.
    :param deltaTemp: Temperature deviation from standard atmosphere.
    :type AC: object
    :type distance: float
    :type altitude: float
    :type M: float
    :type payload: float, optional
    :type fuelReserve: float, optional
    :type flightPlanInitialMass: float, optional
    :type deltaTemp: float, optional
    :return: Estimated initial aircraft mass in kg.
    :rtype: float
    """

    # set Initial Mass from FPL check
    if flightPlanInitialMass is not None:
        initialMass = flightPlanInitialMass
    else:
        # in case of no wind, the ground speed is the same as true airspeed
        [theta, delta, sigma] = atm.atmosphereProperties(
            h=altitude, DeltaTemp=deltaTemp
        )
        TAS = atm.mach2Tas(Mach=M, theta=theta)
        GS = TAS

        if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
            if (AC.MMO is not None and AC.MMO >= 1.0) or (
                AC.VMO is not None and AC.VMO >= 400
            ):
                # identified as fighter jet
                initialMass = AC.MREF
            else:
                cruiseFuelFlow = cruiseFuelConsumption(
                    AC=AC, altitude=altitude, M=M, deltaTemp=deltaTemp
                )
                initialMass = breguetLeducInitialMass(
                    AC=AC,
                    distance=distance,
                    GS=GS,
                    cruiseFuelFlow=cruiseFuelFlow,
                    payload=payload,
                    fuelReserve=fuelReserve,
                )

        elif AC.BADAFamily.BADAH:
            if AC.vne is not None and AC.vne >= 400:
                # identified as fighter
                initialMass = AC.MREF
            else:
                cruiseFuelFlow = cruiseFuelConsumption(
                    AC=AC, altitude=altitude, M=M, deltaTemp=deltaTemp
                )
                initialMass = breguetLeducInitialMass(
                    AC=AC,
                    distance=distance,
                    GS=GS,
                    cruiseFuelFlow=cruiseFuelFlow,
                    payload=payload,
                    fuelReserve=fuelReserve,
                )

    # envelope check
    initialMass = min(max(initialMass, AC.OEW), AC.MTOW)

    return initialMass
