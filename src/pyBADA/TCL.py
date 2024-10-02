# -*- coding: utf-8 -*-
"""
pyBADA
Trajectory Computation Light (TCL) for BADAH/BADAE/BADA3/BADA4 aircraft performance module
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

import os
import numpy as np
from pyBADA.geodesic import Vincenty as vincenty
from pyBADA.geodesic import RhumbLine as rhumb
from pyBADA.geodesic import Turn as turn

from dataclasses import dataclass
import importlib.util
import itertools
import warnings

from math import atan, asin, sin, tan, cos, radians, degrees

from pyBADA import atmosphere as atm
from pyBADA import conversions as conv
from pyBADA import constants as const
from pyBADA.flightTrajectory import FlightTrajectory as FT


@dataclass
class target:
    ROCDtarget: float = None
    slopetarget: float = None
    acctarget: float = None
    ESFtarget: float = None


def checkArgument(argument, **kwargs):
    if kwargs.get(argument) is not None:
        return kwargs.get(argument)
    else:
        raise TypeError("Missing " + argument + " argument")


def constantSpeedLevel(
    AC,
    lengthType,
    length,
    speedType,
    v,
    Hp_init,
    m_init,
    DeltaTemp,
    maxRFL=float("Inf"),
    wS=0.0,
    turnMetrics={"rateOfTurn": 0.0, "bankAngle": 0.0, "directionOfTurn": None},
    stepClimb=False,
    Lat=None,
    Lon=None,
    initialHeading={"magnetic": None, "true": None, "constantHeading": None},
    flightPhase="Cruise",
    magneticDeclinationGrid=None,
    **kwargs,
):
    """This function computes time and fuel required by an aircraft to fly given distance at constant speed in level flight

    :param AC: aircraft {BADA3/4/H/E}
    :param lengthType: what kind of length applies {distance, time}.
    :param length: length of a segment to fly - [NM] distance to fly or [s] time to fly
    :param step_length: length of a step of a segment - [NM] distance to fly or [s] time to fly
    :param speedType: what kind of speed is followed {M, CAS, TAS}.
    :param v: what speed is followed - [kt] CAS/TAS speed to follow or [-] MACH speed to follow.
    :param Hp_init: initial pressure altitude [ft].
    :param m_init: initial aircraft mass [kg].
    :param stepClimb: kind of vertical evolution during cruise {constAlt=False, stepClimb=True}.
    :param HpStep: altitude step for the stepClimb.
    :param DeltaTemp: deviation with respect to ISA [K].
    :param maxRFL: maximum cruise altitude [ft].
    :param wS: longitudinal wind speed (TAS) [kt].
    :param turnMetrics: Metrics for turn performance {"rateOfTurn":0.0,"bankAngle":0.0,"directionOfTurn":None} {[deg/s],[deg],[LEFT/RIGHT]}
    :param SOC_init: initial state of charge [%].
    :param config: aircraft default aerodynamic configuration {TO,IC,CR,AP,LD}.
    :param speedBrakes: deployed or not speedbrakes including value to be added to the drag coeffcient {deployed:False,value:0.03} {deployed:[True/False],value:[-]}.
    :param ROCD_min: lower ROCD threshold to identify the climbing capabilities (service ceiling) [ft/min].
    :param Lat: Geographical Latitude [deg]
    :param Lon: Geographical Longitude [deg]
    :param initialHeading: aircraft magnetic heading, true heading and definition of constant heading(ORTHODROME=False, LOXODROME=True) {[deg],[deg],-}
    :param magneticDeclinationGrid: geographical grid of a magnetic declination on Earth [deg]
    :param mass_const: kind of mass canculation {mass_integrated=False, mass_constant=True}.
    :param m_iter: number of iterations for integration loop [-]
    :param flightPhase: aircraft phase of flight {Climb,Cruise,Descent}
    :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}.
    :type lengthType: string.
    :type length: float.
    :type step_length: float.
    :type speedType: string.
    :type v: float.
    :type Hp_init: float.
    :type m_init: float.
    :type stepClimb: boolean.
    :type HpStep: float.
    :type DeltaTemp: float.
    :type maxRFL: float.
    :type wS: float.
    :type turnMetrics: {float,float,string}.
    :type SOC_init: float.
    :type config: string.
    :type speedBrakes: dict{boolean,float}.
    :type ROCD_min: float.
    :type Lat: float.
    :type Lon: float.
    :type initialHeading: {float,float,boolean}.
    :type magneticDeclinationGrid: magneticDeclinationGrid.
    :type mass_const: boolean.
    :type m_iter: integer.
    :type flightPhase: string.
    :returns:
            BADA3: [Hp, TAS, CAS, M, GS, acc, ROCD, ESF, FUEL, FUELCONSUMED, THR, time, dist, slope, mass, config, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, N, s, NM, deg, kg, -, -,deg,deg,deg,deg,deg,deg/s]
            BADA4: [Hp, TAS, CAS, M, GS, acc, ROCD, ESF, FUEL, FUELCONSUMED, THR, P[Pmec, Pbat, Pelc Ibat, Vbat, Vgbat, SOCr, SOC], time, dist, slope, mass, config, HLid, LG, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, N, [W,W,W,A,V,V,%/h,%], s, NM, deg, kg, -, -, -, -,deg,deg,deg,deg,deg,deg/s]
            BADAH: [Hp, TAS, CAS, M, GS, ROCD, ESF, FUEL, FUELCONSUMED, Peng, Preq, Pav, time, dist, slope, mass, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, W, W, W, s, NM, deg, kg, -,deg,deg,deg,deg,deg,deg/s]
            BADAE: [time, dist, Hp, TAS, CAS, M, GS, acc, ROCD, ESF, slope, mass, P[Pmec, Pelc, Pbat, SOCr, SOC, Ibat, Vbat, Vgbat] comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [s, NM, ft, kt, kt, -, kt, m/s^2, ft/min, deg, kg, [W,W,W,%/h,%,A,V,V], -,deg,deg,deg,deg,deg,deg/s]
    :rtype: dict[list[float]}.
    """

    rateOfTurn = turnMetrics["rateOfTurn"]
    bankAngle = turnMetrics["bankAngle"]
    directionOfTurn = turnMetrics["directionOfTurn"]

    turnFlight = True
    if turnMetrics["rateOfTurn"] == 0.0 and turnMetrics["bankAngle"] == 0.0:
        turnFlight = False

    # conversion of Magnetic Heading to True Heading
    if magneticDeclinationGrid is not None:
        magneticDeclination = magneticDeclinationGrid.getMagneticDeclination(
            LAT_target=Lat, LON_target=Lon
        )
    else:
        magneticDeclination = 0

    # retrieve magnetic and true heading inputs
    magneticHeading = initialHeading["magnetic"]
    trueHeading = initialHeading["true"]
    constantHeading = initialHeading["constantHeading"]

    if Lat and Lon and (magneticHeading or trueHeading):
        if trueHeading is not None and magneticHeading is None:
            # fly TRUE Heading
            headingToFly = "TRUE"
            magneticHeading = trueHeading - magneticDeclination
        elif magneticHeading is not None and trueHeading is None:
            # fly MAGNETIC Heading
            if constantHeading == True:
                headingToFly = "MAGNETIC"
                trueHeading = magneticHeading + magneticDeclination
            else:
                raise Exception("Cannot fly non-constant magnetic heading")

        else:
            raise Exception("Undefined Heading value combination")

    # calculation with constant mass (True) or integrated (False)
    mass_const = kwargs.get("mass_const", False)

    # optional parameter to define initial Baterry State of Charge (SOC)
    if AC.BADAFamily.BADAE:
        SOC_init = kwargs.get("SOC_init", 100)
    else:
        SOC_init = None

    # speed brakes application
    if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
        speedBrakes = kwargs.get("speedBrakes", {"deployed": False, "value": 0.03})

    # optional parameter - iteration step length based on the type of aircraft
    if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
        step_length = kwargs.get(
            "step_length", 10
        )  # [NM] or [s] based on the 'lengthType'
    else:
        step_length = kwargs.get(
            "step_length", 100
        )  # [NM] or [s] based on the 'lengthType'

    #  weight iteration constant
    if AC.BADAFamily.BADAE:
        m_iter = kwargs.get("m_iter", 1)  # number of iterations for integration loop[-]
    else:
        m_iter = kwargs.get("m_iter", 2)  # number of iterations for integration loop[-]

    # comment line describing type of trajectory calculation
    if flightPhase != "Cruise":
        levelOffComment = "_levelOff"
    else:
        levelOffComment = ""

    if turnFlight:
        turnComment = "_turn"
    else:
        turnComment = ""

    if constantHeading:
        constHeadingStr = "_const_Heading"
    elif constantHeading is False or constantHeading is None:
        constHeadingStr = ""

    comment = (
        flightPhase
        + levelOffComment
        + turnComment
        + "_const_"
        + speedType
        + constHeadingStr
    )

    if Lat and Lon and (magneticHeading or trueHeading):
        comment = comment + "_" + headingToFly + "_Heading"

    # Altitude step for step cruise
    HpStep = kwargs.get("HpStep", 2000)  # [ft]

    # minimum remaining ROCD to determine cruise ceiling
    if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
        ROCD_min = kwargs.get("ROCD_min", 50)  # [ft/min]
    else:
        if AC.engineType == "PISTON" or AC.engineType == "ELECTRIC":
            ROCD_min = kwargs.get("ROCD_min", 100)  # [ft/min]
        else:
            ROCD_min = kwargs.get("ROCD_min", 300)  # [ft/min]

    # get the default aerodynamic configuration if provided to be used for the whole segment
    config_default = kwargs.get("config", None)
    if config_default is not None:
        if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
            if not (
                config_default == "TO"
                or config_default == "IC"
                or config_default == "CR"
                or config_default == "AP"
                or config_default == "LD"
            ):
                print(
                    "WRONG default configuration set. Available values are: TO/IC/CR/AP/LD. Configuration will be calculated automatically"
                )

    # initialize output parameters
    Hp = []
    TAS = []
    CAS = []
    GS = []
    M = []
    ROCD = []
    esf = []
    FUEL = []
    FUELCONSUMED = []
    Comment = []
    Slope = []
    acc = []
    THR = []
    DRAG = []
    config = []
    HLid = []
    LG = []
    time = []
    dist = []
    mass = []
    BankAngle = []
    ROT = []

    # BADAH specific
    Preq = []
    Peng = []
    Pav = []

    # optional GPS coordiantes and HDG definition
    LAT = []
    LON = []
    HDGMagnetic = []
    HDGTrue = []

    # BADAE specific
    Pmec = []
    Pbat = []
    SOCr = []
    SOC = []
    Pelc = []
    Ibat = []
    Vbat = []
    Vgbat = []

    if AC.BADAFamily.BADAE:
        SOC_i = SOC_init

    # init loop parameters
    totalLength = 0
    Hp_i = Hp_init
    mass_i = m_init
    time_i = 0.0
    dist_i = 0.0
    fuelConsumed_i = 0.0

    Lat_i = Lat
    Lon_i = Lon
    HDGMagnetic_i = magneticHeading
    HDGTrue_i = trueHeading

    while True:
        ## PART 1: compute parameters at end of step that are known without uncertainties:
        ##         atmosphere, speeds, elapsed time

        # atmosphere properties
        H_m = conv.ft2m(Hp_i)  # altitude [m]
        [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
        # aircraft speed
        [M_i, CAS_i, TAS_i] = atm.convertSpeed(
            v=v, speedType=speedType, theta=theta, delta=delta, sigma=sigma
        )
        GS_i = conv.ms2kt(TAS_i) - wS  # ground speed [kt]

        if turnFlight:
            if turnMetrics["bankAngle"] != 0.0:
                # bankAngle is defined
                rateOfTurn = AC.rateOfTurn_bankAngle(TAS=TAS_i, bankAngle=bankAngle)
            else:
                # rateOfTurn is defined
                bankAngle = AC.bankAngle(rateOfTurn=rateOfTurn, v=TAS_i)  # [degrees]

        if lengthType == "distance":
            # step time is: distance differantial divided by ground speed
            step_distance = totalLength - dist_i  # [NM]
            step_time = 3600 * step_distance / GS_i  # [s]
        elif lengthType == "time":
            step_time = totalLength - time_i  # [s]

            if turnFlight:
                step_distance = conv.m2nm(
                    turn.distance(
                        rateOfTurn=rateOfTurn, TAS=TAS_i, timeOfTurn=step_time
                    )
                )  # arcLength during the turn [NM]
            else:
                step_distance = GS_i * step_time / 3600  # [NM]

        # Load factor
        nz = 1 / cos(radians(bankAngle))

        ## PART 2: compute parameters at end of step that are known only with uncertainties
        ##           (due to unknown mass at end of step):
        ##         weight, lift, drag , thrust, fuel flow

        for _ in itertools.repeat(None, m_iter):
            # BADAH or BADAE
            if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
                # compute Power required for level flight
                Preq_i = AC.Preq(sigma=sigma, tas=TAS_i, mass=mass_i, phi=bankAngle)
                Peng_i = Preq_i
                if AC.BADAFamily.BADAH:
                    Pav_i = AC.Pav(
                        rating="MCNT", theta=theta, delta=delta
                    )  # assume MCNT rating as the limit
                elif AC.BADAFamily.BADAE:
                    Pav_i = AC.Pav(
                        rating="MCNT", SOC=SOC_i
                    )  # assume MCNT rating as the limit

                if Pav_i < Preq_i:
                    warnings.warn(
                        "Power Available is lower than Power Required", UserWarning
                    )

                # BADAH
                if AC.BADAFamily.BADAH:
                    # compute fuel flow for level flight
                    CP = AC.CP(Peng=Preq_i)
                    FUEL_i = AC.ff(delta=delta, CP=CP)  # [kg/s]

                # BADAE
                elif AC.BADAFamily.BADAE:
                    Pbat_i = AC.Pbat(Preq=Preq_i, SOC=SOC_i)
                    SOCr_i = AC.SOCrate(Preq=Preq_i, SOC=SOC_i)

                    # debug data
                    Pelc_i = Preq_i / AC.eta
                    Ibat_i = AC.Ibat(P=Pelc_i, SOC=SOC_i)
                    Vbat_i = AC.Vbat(I=Ibat_i, SOC=SOC_i)
                    Vgbat_i = AC.Vocbat(SOC=SOC_i) - AC.R0bat(SOC=SOC_i) * Ibat_i

            # BADA4
            elif AC.BADAFamily.BADA4:
                # aircraft configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=flightPhase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=flightPhase,
                        previousConfig=config[-1],
                        currentConfig=config_i,
                    )

                [HLid_i, LG_i] = AC.flightEnvelope.getAeroConfig(config=config_i)

                # compute lift coefficient
                CL = AC.CL(M=M_i, delta=delta, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(M=M_i, CL=CL, HLid=HLid_i, LG=LG_i, speedBrakes=speedBrakes)
                # compute drag force
                Drag = AC.D(M=M_i, delta=delta, CD=CD)
                # compute thrust force and fuel flow
                THR_i = Drag
                CT = AC.CT(Thrust=THR_i, delta=delta)
                FUEL_i = AC.ff(
                    CT=CT, delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )  # [kg/s]

            # BADA3
            elif AC.BADAFamily.BADA3:
                # aircraft configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=flightPhase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=flightPhase,
                        previousConfig=config[-1],
                        currentConfig=config_i,
                    )

                # compute lift coefficient
                CL = AC.CL(tas=TAS_i, sigma=sigma, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(CL=CL, config=config_i, speedBrakes=speedBrakes)
                # compute drag force
                Drag = AC.D(tas=TAS_i, sigma=sigma, CD=CD)
                # compute thrust force and fuel flow
                THR_i = Drag

                if flightPhase == "Descent":
                    FUEL_i = AC.ff(
                        flightPhase=flightPhase,
                        v=TAS_i,
                        h=H_m,
                        T=THR_i,
                        config=config_i,
                        adapted=True,
                    )
                else:
                    FUEL_i = AC.ff(
                        flightPhase=flightPhase,
                        v=TAS_i,
                        h=H_m,
                        T=THR_i,
                        config=config_i,
                        adapted=False,
                    )

            # compute fuel burn over current step
            if not time:
                break
            else:
                time_i = time[-1] + step_time
                dist_i = dist[-1] + step_distance

                if (
                    Lat_i is not None
                    and Lon_i is not None
                    and (HDGMagnetic_i is not None or HDGTrue is not None)
                ):

                    if headingToFly == "TRUE":
                        if not turnFlight:
                            if not constantHeading:
                                # fly ORTHODROME
                                (Lat_i, Lon_i, HDGTrue_i) = (
                                    vincenty.destinationPoint_finalBearing(
                                        LAT_init=LAT[-1],
                                        LON_init=LON[-1],
                                        distance=conv.nm2m(step_distance),
                                        bearing=HDGTrue[-1],
                                    )
                                )

                                if magneticDeclinationGrid is not None:
                                    HDGMagnetic_i = (
                                        HDGTrue_i
                                        - magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGMagnetic_i = HDGTrue_i

                            elif constantHeading:
                                # fly LOXODROME
                                (Lat_i, Lon_i) = rhumb.destinationPoint(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearing=HDGTrue[-1],
                                    distance=conv.nm2m(step_distance),
                                )

                                if magneticDeclinationGrid is not None:
                                    HDGMagnetic_i = (
                                        HDGTrue_i
                                        - magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGMagnetic_i = HDGTrue_i

                        else:
                            # calculate the turn
                            (Lat_i, Lon_i, HDGTrue_i) = (
                                turn.destinationPoint_finalBearing(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearingInit=HDGTrue[-1],
                                    TAS=TAS_i,
                                    rateOfTurn=rateOfTurn,
                                    timeOfTurn=step_time,
                                    directionOfTurn=directionOfTurn,
                                )
                            )

                            if magneticDeclinationGrid is not None:
                                HDGMagnetic_i = (
                                    HDGTrue_i
                                    - magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGMagnetic_i = HDGTrue_i

                    elif headingToFly == "MAGNETIC":
                        if not turnFlight:
                            if constantHeading:
                                (Lat_i, Lon_i) = rhumb.destinationPoint(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearing=HDGTrue[-1],
                                    distance=conv.nm2m(step_distance),
                                )

                                if magneticDeclinationGrid is not None:
                                    HDGTrue_i = (
                                        HDGMagnetic_i
                                        + magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGTrue_i = HDGMagnetic_i

                        else:
                            # calculate the turn
                            (Lat_i, Lon_i, HDGTrue_i) = (
                                turn.destinationPoint_finalBearing(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearingInit=HDGTrue[-1],
                                    TAS=TAS_i,
                                    rateOfTurn=rateOfTurn,
                                    timeOfTurn=step_time,
                                    directionOfTurn=directionOfTurn,
                                )
                            )

                            if magneticDeclinationGrid is not None:
                                HDGMagnetic_i = (
                                    HDGTrue_i
                                    - magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGMagnetic_i = HDGTrue_i

                # BADAE
                if AC.BADAFamily.BADAE:
                    # Average SOC rate over step is the mean of initial and final ones
                    step_SOCr = (SOCr[-1] + SOCr_i) / 2  # [%/h]
                    # SOC change is: average SOC rate multiplied by step time
                    step_SOC = step_SOCr * step_time / 3600  # [%]
                    # Update SOC estimate at end of step
                    SOC_i = SOC[-1] - step_SOC  # [%]
                    # update of aircraft mass estimate at end of step - mass is not changing for ELECTRIC engine (no fuel is consumed)
                    mass_i = mass[-1]  # [kg]
                else:
                    # average fuel flow over step is the mean of initial and final ones
                    step_FUEL = (FUEL[-1] + FUEL_i) / 2  # [kg/s]
                    # fuel burnt is: average fuel flow multiplied by step time
                    step_mass = step_FUEL * step_time  # [kg]
                    # update of aircraft mass estimate at end of step
                    if not mass_const:
                        mass_i = mass[-1] - step_mass  # [kg]
                        fuelConsumed_i = step_FUEL * step_time
                    fuelConsumed_i = FUELCONSUMED[-1] + step_FUEL * step_time

        ## PART 3: store information about end of step point

        # point data
        Hp.append(Hp_i)
        TAS.append(conv.ms2kt(TAS_i))
        CAS.append(conv.ms2kt(CAS_i))
        GS.append(GS_i)
        M.append(M_i)
        ROCD.append(0.0)
        Comment.append(comment)
        mass.append(mass_i)
        esf.append(0.0)

        Slope.append(0.0)  # slope (since ROCD=0, then slope=0)
        acc.append(0.0)

        BankAngle.append(bankAngle)
        ROT.append(rateOfTurn)

        time.append(time_i)
        dist.append(dist_i)

        if (
            Lat_i is not None
            and Lon_i is not None
            and (HDGMagnetic_i is not None or HDGTrue_i is not None)
        ):
            LAT.append(Lat_i)
            LON.append(Lon_i)
            HDGMagnetic.append(HDGMagnetic_i)
            HDGTrue.append(HDGTrue_i)

        # everything except electric BADAE
        if not AC.BADAFamily.BADAE:
            FUEL.append(FUEL_i)
            FUELCONSUMED.append(fuelConsumed_i)

        # BADAH
        if AC.BADAFamily.BADAH:
            Preq.append(Preq_i)
            Peng.append(Peng_i)
            Pav.append(Pav_i)

        # BADAE
        elif AC.BADAFamily.BADAE:
            Pmec.append(Peng_i)
            Pbat.append(Pbat_i)
            SOCr.append(SOCr_i)
            SOC.append(SOC_i)
            Pelc.append(Pelc_i)
            Ibat.append(Ibat_i)
            Vbat.append(Vbat_i)
            Vgbat.append(Vgbat_i)

        # BADA3 & BADA4
        elif AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
            THR.append(THR_i)
            DRAG.append(Drag)
            config.append(config_i)

        # BADA4
        if AC.BADAFamily.BADA4:
            HLid.append(HLid_i)
            LG.append(LG_i)

        ## PART 4: check whether a climb step should be performed (applicable for BADA3 and BADA4)
        if (
            stepClimb
            and totalLength != 0
            and (AC.BADAFamily.BADA4 or AC.BADAFamily.BADA3)
        ):
            # determine atmosphere properties at upper cruise altitude
            nextHp = min(Hp_i + HpStep, maxRFL)
            H_m = conv.ft2m(nextHp)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)

            # aircraft speed at upper cruise altitude
            [M_up, CAS_up, TAS_up] = atm.convertSpeed(
                v=v, speedType=speedType, theta=theta, delta=delta, sigma=sigma
            )

            # Determine fuel flow at upper cruise altitude:
            # BADA3
            if AC.BADAFamily.BADA3:
                # aircraft configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=flightPhase,
                        v=CAS_up,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=flightPhase,
                        previousConfig=config[-1],
                        currentConfig=config_i,
                    )

                # compute lift coefficient
                CL = AC.CL(tas=TAS_up, sigma=sigma, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(CL=CL, config=config_i, speedBrakes=speedBrakes)
                # compute drag force
                Drag = AC.D(tas=TAS_up, sigma=sigma, CD=CD)
                # compute thrust force and fuel flow
                THR_up = Drag
                FUEL_up = AC.ff(
                    flightPhase="Cruise",
                    v=TAS_up,
                    h=H_m,
                    T=THR_up,
                    config=config_i,
                    adapted=False,
                )

            # BADA4
            elif AC.BADAFamily.BADA4:
                # aircraft configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=flightPhase,
                        v=CAS_up,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=flightPhase,
                        previousConfig=config[-1],
                        currentConfig=config_i,
                    )

                [HLid_i, LG_i] = AC.flightEnvelope.getAeroConfig(config=config_i)

                # compute lift coefficient
                CL = AC.CL(M=M_up, delta=delta, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(M=M_up, CL=CL, HLid=HLid_i, LG=LG_i, speedBrakes=speedBrakes)
                # compute drag force
                Drag = AC.D(M=M_up, delta=delta, CD=CD)
                # compute thrust force and fuel flow
                THR_up = Drag
                CT = AC.CT(Thrust=THR_up, delta=delta)
                FUEL_up = AC.ff(
                    CT=CT, delta=delta, theta=theta, M=M_up, DeltaTemp=DeltaTemp
                )  # [kg/s]

            # Compare specific range at current and upper cruise altitudes
            if (TAS_up / FUEL_up) > (TAS_i / FUEL_i):
                # Check that available ROCD at upper cruise altitude allows a climb step
                if AC.BADAFamily.BADA4:
                    THR_CL = AC.Thrust(
                        rating="MCMB",
                        delta=delta,
                        theta=theta,
                        M=M_up,
                        DeltaTemp=DeltaTemp,
                    )  # MCMB Thrust
                elif AC.BADAFamily.BADA3:
                    THR_CL = AC.Thrust(
                        rating="MCMB",
                        v=TAS,
                        h=H_m,
                        DeltaTemp=DeltaTemp,
                        config=config_i,
                    )  # MCMB Thrust

                flightEvolution = "const" + speedType
                ESF_i = AC.esf(
                    h=H_m, flightEvolution=flightEvolution, M=M_up, DeltaTemp=DeltaTemp
                )
                temp_const = (theta * const.temp_0) / (
                    theta * const.temp_0 - DeltaTemp
                )  # T/T-dT
                ROCD_up = (
                    conv.m2ft(
                        (1 / temp_const)
                        * (THR_CL - Drag)
                        * TAS_up
                        * ESF_i
                        / (mass_i * const.g)
                    )
                    * 60
                )  # [ft/min]

                if ROCD_up >= ROCD_min:
                    # Compute climb step
                    if speedType == "M":
                        speed = M_up
                    elif speedType == "TAS":
                        speed = conv.ms2kt(TAS_up)
                    elif speedType == "CAS":
                        speed = conv.ms2kt(CAS_up)

                    if Lat and Lon and magneticHeading:
                        flightTrajectory_CL = constantSpeedRating(
                            AC=AC,
                            speedType=speedType,
                            v=speed,
                            Hp_init=Hp_i,
                            Hp_final=nextHp,
                            Hp_step=HpStep,
                            m_init=mass_i,
                            wS=wS,
                            DeltaTemp=DeltaTemp,
                            Lat=LAT[-1],
                            Lon=LON[-1],
                            initialHeading={
                                "magnetic": HDGMagnetic[-1],
                                "true": HDGTrue[-1],
                                "constantHeading": constantHeading,
                            },
                            turnMetrics=turnMetrics,
                            magneticDeclinationGrid=magneticDeclinationGrid,
                        )
                    else:
                        flightTrajectory_CL = constantSpeedRating(
                            AC=AC,
                            speedType=speedType,
                            v=speed,
                            Hp_init=Hp_i,
                            Hp_final=nextHp,
                            Hp_step=HpStep,
                            m_init=mass_i,
                            wS=wS,
                            DeltaTemp=DeltaTemp,
                            turnMetrics=turnMetrics,
                        )

                    # Avoid a step just before the Top of Descent, which can cause stability issues in the main cruise+descent loop
                    if lengthType == "distance":
                        length_CL = flightTrajectory_CL["dist"].iloc[-1]
                    elif lengthType == "time":
                        length_CL = flightTrajectory_CL["time"].iloc[-1]

                    if (totalLength + length_CL) < length - step_length:
                        # add stepClimb segment data to previously calculated data - combine segments
                        Hp.extend(flightTrajectory_CL["Hp"])
                        TAS.extend(flightTrajectory_CL["TAS"])
                        CAS.extend(flightTrajectory_CL["CAS"])
                        GS.extend(flightTrajectory_CL["GS"])
                        M.extend(flightTrajectory_CL["M"])
                        esf.extend(flightTrajectory_CL["ESF"])
                        ROCD.extend(flightTrajectory_CL["ROCD"])
                        Slope.extend(flightTrajectory_CL["slope"])
                        acc.extend(flightTrajectory_CL["acc"])
                        THR.extend(flightTrajectory_CL["THR"])
                        DRAG.extend(flightTrajectory_CL["DRAG"])
                        config.extend(flightTrajectory_CL["config"])
                        FUEL.extend(flightTrajectory_CL["FUEL"])
                        mass.extend(flightTrajectory_CL["mass"])
                        BankAngle.extend(flightTrajectory_CL["BankAngle"])
                        ROT.extend(flightTrajectory_CL["ROT"])

                        comment_CL = flightTrajectory_CL["comment"]
                        Comment.extend([com + "_stepClimb" for com in comment_CL])

                        # BADA4
                        if AC.BADAFamily.BADA4:
                            HLid.extend(flightTrajectory_CL["HLid"])
                            LG.extend(flightTrajectory_CL["LG"])

                        time_seg1 = time
                        for k in flightTrajectory_CL["time"]:
                            time.append(time_seg1[-1] + k)

                        dist_seg1 = dist
                        for k in flightTrajectory_CL["dist"]:
                            dist.append(dist_seg1[-1] + k)

                        FUELCONSUMED_seg1 = FUELCONSUMED
                        for k in flightTrajectory_CL["FUELCONSUMED"]:
                            FUELCONSUMED.append(FUELCONSUMED_seg1[-1] + k)

                        if Lat and Lon and magneticHeading:
                            LAT.extend(flightTrajectory_CL["LAT"])
                            LON.extend(flightTrajectory_CL["LON"])
                            HDGMagnetic.extend(flightTrajectory_CL["HDGMagnetic"])
                            HDGTrue.extend(flightTrajectory_CL["HDGTrue"])

                        # Compute cruise fuel at upper altitude
                        # BADA3
                        if AC.BADAFamily.BADA3:
                            # aircraft configuration
                            if config_default is None:
                                config_i = AC.flightEnvelope.getConfig(
                                    h=H_m,
                                    phase=flightPhase,
                                    v=CAS_up,
                                    mass=mass[-1],
                                    DeltaTemp=DeltaTemp,
                                )
                            else:
                                config_i = config_default

                            # ensure continuity of configuration change within the segment
                            if config:
                                config_i = (
                                    AC.flightEnvelope.checkConfigurationContinuity(
                                        phase=flightPhase,
                                        previousConfig=config[-1],
                                        currentConfig=config_i,
                                    )
                                )

                            # compute lift coefficient
                            CL = AC.CL(tas=TAS_up, sigma=sigma, mass=mass[-1], nz=nz)
                            # compute drag coefficient
                            CD = AC.CD(CL=CL, config=config_i, speedBrakes=speedBrakes)
                            # compute drag force
                            Drag = AC.D(tas=TAS_up, sigma=sigma, CD=CD)
                            # compute thrust force and fuel flow
                            THR_up = Drag
                            FUEL_up = AC.ff(
                                flightPhase="Cruise",
                                v=TAS_up,
                                h=H_m,
                                T=THR_up,
                                config=config_i,
                                adapted=False,
                            )

                        # BADA4
                        elif AC.BADAFamily.BADA4:
                            # aircraft configuration
                            if config_default is None:
                                config_i = AC.flightEnvelope.getConfig(
                                    h=H_m,
                                    phase=flightPhase,
                                    v=CAS_up,
                                    mass=mass[-1],
                                    DeltaTemp=DeltaTemp,
                                )
                            else:
                                config_i = config_default

                            # ensure continuity of configuration change within the segment
                            if config:
                                config_i = (
                                    AC.flightEnvelope.checkConfigurationContinuity(
                                        phase=flightPhase,
                                        previousConfig=config[-1],
                                        currentConfig=config_i,
                                    )
                                )

                            [HLid_i, LG_i] = AC.flightEnvelope.getAeroConfig(
                                config=config_i
                            )

                            # compute lift coefficient
                            CL = AC.CL(M=M_up, delta=delta, mass=mass[-1], nz=nz)
                            # compute drag coefficient
                            CD = AC.CD(
                                M=M_up,
                                CL=CL,
                                HLid=HLid_i,
                                LG=LG_i,
                                speedBrakes=speedBrakes,
                            )
                            # compute drag force
                            Drag = AC.D(M=M_up, delta=delta, CD=CD)
                            # compute thrust force and fuel flow
                            THR_up = Drag
                            CT = AC.CT(Thrust=THR_up, delta=delta)
                            FUEL_up = AC.ff(
                                CT=CT,
                                delta=delta,
                                theta=theta,
                                M=M_up,
                                DeltaTemp=DeltaTemp,
                            )  # [kg/s]

                        Hp.append(Hp[-1])
                        TAS.append(TAS[-1])
                        CAS.append(CAS[-1])
                        GS.append(GS[-1])
                        esf.append(esf[-1])
                        M.append(M[-1])
                        ROCD.append(0)
                        Comment.append(comment)
                        Slope.append(0)
                        acc.append(0)
                        THR.append(THR_up)
                        DRAG.append(Drag)
                        config.append(config_i)
                        FUEL.append(FUEL_up)
                        FUELCONSUMED.append(FUELCONSUMED[-1])
                        mass.append(mass[-1])
                        BankAngle.append(BankAngle[-1])
                        ROT.append(ROT[-1])
                        time.append(time[-1])
                        dist.append(dist[-1])

                        if Lat and Lon and magneticHeading:
                            LAT.append(LAT[-1])
                            LON.append(LON[-1])
                            HDGMagnetic.append(HDGMagnetic[-1])
                            HDGTrue.append(HDGTrue[-1])

                        # BADA4
                        if AC.BADAFamily.BADA4:
                            HLid.append(HLid_i)
                            LG.append(LG_i)

                        Hp_i = Hp[-1]

                        if lengthType == "distance":
                            totalLength = dist[-1]
                        elif lengthType == "time":
                            totalLength = time[-1]

        if totalLength + step_length < length:
            totalLength += step_length
        elif totalLength < length:
            totalLength = length
        else:
            break

    flightData = {
        "Hp": Hp,
        "TAS": TAS,
        "CAS": CAS,
        "GS": GS,
        "M": M,
        "acc": acc,
        "ROCD": ROCD,
        "ESF": esf,
        "Pmec": Pmec,
        "Pelc": Pelc,
        "Pbat": Pbat,
        "SOCr": SOCr,
        "SOC": SOC,
        "Ibat": Ibat,
        "Vbat": Vbat,
        "Vgbat": Vgbat,
        "FUEL": FUEL,
        "FUELCONSUMED": FUELCONSUMED,
        "Preq": Preq,
        "Peng": Peng,
        "Pav": Pav,
        "THR": THR,
        "DRAG": DRAG,
        "time": time,
        "dist": dist,
        "slope": Slope,
        "mass": mass,
        "config": config,
        "HLid": HLid,
        "LG": LG,
        "LAT": LAT,
        "LON": LON,
        "HDGTrue": HDGTrue,
        "HDGMagnetic": HDGMagnetic,
        "BankAngle": BankAngle,
        "ROT": ROT,
        "comment": Comment,
    }

    flightTrajectory = FT.createFlightTrajectoryDataframe(flightData)
    return flightTrajectory


def constantSpeedROCD(
    AC,
    speedType,
    v,
    Hp_init,
    Hp_final,
    ROCDtarget,
    m_init,
    DeltaTemp,
    wS=0.0,
    turnMetrics={"rateOfTurn": 0.0, "bankAngle": 0.0, "directionOfTurn": None},
    Lat=None,
    Lon=None,
    initialHeading={"magnetic": None, "true": None, "constantHeading": None},
    reducedPower=None,
    directionOfTurn=None,
    magneticDeclinationGrid=None,
    **kwargs,
):
    """This function computes time and fuel required by an aircraft to perform a climb/descent from Hp_init to Hp_final at constant speed and constant rate of climb/descent

    :param AC: aircraft {BADA3/4/H/E}
    :param speedType: what kind of speed is followed {M, CAS, TAS}.
    :param v: what kind of speed is followed - [kt] CAS/TAS speed to follow or [-] MACH speed to follow.
    :param Hp_init: initial pressure altitude [ft].
    :param Hp_final: final pressure altitude [ft].
    :param ROCDtarget: Rate of climb/descent to be followed [ft/min].
    :param m_init: initial aircraft mass [kg].
    :param DeltaTemp: deviation with respect to ISA [K].
    :param wS: longitudinal wind speed (TAS) [kt].
    :param turnMetrics: Metrics for turn performance {"rateOfTurn":0.0,"bankAngle":0.0,"directionOfTurn":None} {[deg/s],[deg],[LEFT/RIGHT]}
    :param Hp_step: length of an altitude step of a segment [ft].
    :param SOC_init: initial state of charge [%].
    :param config: aircraft default aerodynamic configuration {TO,IC,CR,AP,LD}.
    :param speedBrakes: deployed or not speedbrakes including value to be added to the drag coeffcient {deployed:False,value:0.03} {deployed:[True/False],value:[-]}.
    :param ROCD_min: lower ROCD threshold to identify the climbing capabilities (service ceiling) [ft/min].
    :param Lat: Geographical Latitude [deg]
    :param Lon: Geographical Longitude [deg]
    :param initialHeading: aircraft magnetic heading, true heading and definition of constant heading(ORTHODROME=False, LOXODROME=True) {[deg],[deg],-}
    :param magneticDeclinationGrid: geographical grid of a magnetic declination on Earth [deg]
    :param mass_const: kind of mass canculation {mass_integrated=False, mass_constant=True}.
    :param m_iter: number of iterations for integration loop [-]
    :param reducedPower: reduction of Power during the climb {True/False}
    :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}.
    :type speedType: string.
    :type v: float.
    :type Hp_init: float.
    :type Hp_final: float.
    :type ROCDtarget: float.
    :type m_init: float.
    :type DeltaTemp: float.
    :type wS: float.
    :type turnMetrics: {float,float,string}.
    :type Hp_step: float.
    :type SOC_init: float.
    :type config: string.
    :type speedBrakes: dict{boolean,float}.
    :type ROCD_min: float.
    :type Lat: float.
    :type Lon: float.
    :type initialHeading: {float,float,boolean}.
    :type magneticDeclinationGrid: magneticDeclinationGrid.
    :type mass_const: boolean.
    :type m_iter: integer.
    :type reducedPower: boolean.
    :returns:
            BADA3: [Hp, TAS, CAS, M, GS, acc, ROCD, ESF, FUEL, FUELCONSUMED, THR, time, dist, slope, mass, config, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, N, s, NM, deg, kg, -, -,deg,deg,deg,deg,deg,deg/s]
            BADA4: [Hp, TAS, CAS, M, GS, acc, ROCD, ESF, FUEL, FUELCONSUMED, THR, P[Pmec, Pbat, Pelc Ibat, Vbat, Vgbat, SOCr, SOC], time, dist, slope, mass, config, HLid, LG, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, N, [W,W,W,A,V,V,%/h,%], s, NM, deg, kg, -, -, -, -,deg,deg,deg,deg,deg,deg/s]
            BADAH: [Hp, TAS, CAS, M, GS, ROCD, ESF, FUEL, FUELCONSUMED, Peng, Preq, Pav, time, dist, slope, mass, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, W, W, W, s, NM, deg, kg, -,deg,deg,deg,deg,deg,deg/s]
            BADAE: [time, dist, Hp, TAS, CAS, M, GS, acc, ROCD, ESF, slope, mass, P[Pmec, Pelc, Pbat, SOCr, SOC, Ibat, Vbat, Vgbat] comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [s, NM, ft, kt, kt, -, kt, m/s^2, ft/min, deg, kg, [W,W,W,%/h,%,A,V,V], -,deg,deg,deg,deg,deg,deg/s]
    :rtype: dict[list[float]}.
    """

    rateOfTurn = turnMetrics["rateOfTurn"]
    bankAngle = turnMetrics["bankAngle"]
    directionOfTurn = turnMetrics["directionOfTurn"]

    turnFlight = True
    if turnMetrics["rateOfTurn"] == 0.0 and turnMetrics["bankAngle"] == 0.0:
        turnFlight = False

    # conversion of Magnetic Heading to True Heading
    if magneticDeclinationGrid is not None:
        magneticDeclination = magneticDeclinationGrid.getMagneticDeclination(
            LAT_target=Lat, LON_target=Lon
        )
    else:
        magneticDeclination = 0

    # retrieve magnetic and true heading inputs
    magneticHeading = initialHeading["magnetic"]
    trueHeading = initialHeading["true"]
    constantHeading = initialHeading["constantHeading"]

    if Lat and Lon and (magneticHeading or trueHeading):
        if trueHeading is not None and magneticHeading is None:
            # fly TRUE Heading
            headingToFly = "TRUE"
            magneticHeading = trueHeading - magneticDeclination
        elif magneticHeading is not None and trueHeading is None:
            # fly MAGNETIC Heading
            if constantHeading == True:
                headingToFly = "MAGNETIC"
                trueHeading = magneticHeading + magneticDeclination
            else:
                raise Exception("Cannot fly non-constant magnetic heading")

        else:
            raise Exception("Undefined Heading value combination")

    # calculation with constant mass (True) or integrated (False)
    mass_const = kwargs.get("mass_const", False)

    # optional parameter to define initial Baterry State of Charge (SOC)
    if AC.BADAFamily.BADAE:
        SOC_init = kwargs.get("SOC_init", 100)
    else:
        SOC_init = None

    # speed brakes application
    if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
        speedBrakes = kwargs.get("speedBrakes", {"deployed": False, "value": 0.03})

    # optional parameter - iteration step for altitude loop
    if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
        Hp_step = kwargs.get("Hp_step", 500)  # [ft]
    else:
        # NB: it must be a multiple of 1000ft so that interrupted climbs end on a regular cruise altitude.
        Hp_step = kwargs.get("Hp_step", 1000)  # [ft]

    # minimum remaining ROCD to determine cruise ceiling
    if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
        ROCD_min = kwargs.get("ROCD_min", 50)  # [ft/min]
    else:
        if AC.engineType == "PISTON" or AC.engineType == "ELECTRIC":
            ROCD_min = kwargs.get("ROCD_min", 100)  # [ft/min]
        else:
            ROCD_min = kwargs.get("ROCD_min", 300)  # [ft/min]

    # determine if  vertical evolution over the segment is CLIMB or DESCENT
    if Hp_init < Hp_final:
        phase = "Climb"
    else:
        phase = "Descent"
        Hp_step = -Hp_step

    # phase of flight defined in case of no altitude step available, where Hp_init = Hp_final
    # phase = kwargs.get('phase', None)

    # check the consistency of ROCD and climb/descent phase of flight
    # if incosistent, change the sign on ROCD target value
    if phase == "Climb" and ROCDtarget < 0:
        ROCDtarget = abs(ROCDtarget)
        print("ROCDtarget for Climb should be positive")
    elif phase == "Descent" and ROCDtarget > 0:
        ROCDtarget = ROCDtarget * (-1)
        print("ROCDtarget for Descent should be negative")

    if turnFlight:
        turnComment = "_turn"
    else:
        turnComment = ""

    if constantHeading:
        constHeadingStr = "_const_Heading"
    elif constantHeading is False or constantHeading is None:
        constHeadingStr = ""

    # comment line describing type of trajectory calculation
    comment = phase + turnComment + "_const_ROCD_" + speedType + constHeadingStr

    if Lat and Lon and (magneticHeading or trueHeading):
        comment = comment + "_" + headingToFly + "_Heading"

    # get the default aerodynamic configuration if provided to be used for the whole segment
    config_default = kwargs.get("config", None)
    if config_default is not None:
        if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
            if not (
                config_default == "TO"
                or config_default == "IC"
                or config_default == "CR"
                or config_default == "AP"
                or config_default == "LD"
            ):
                print(
                    "WRONG default configuration set. Available values are: TO/IC/CR/AP/LD. Configuration will be calculated automatically"
                )

    #  weight iteration constant
    m_iter = kwargs.get("m_iter", 5)  # number of iterations for integration loop[-]

    # convert ROCD to IS units
    ROCDisu = conv.ft2m(ROCDtarget) / 60

    # initialize output parameters
    Hp = []
    TAS = []
    CAS = []
    GS = []
    M = []
    ROCD = []
    esf = []
    FUEL = []
    FUELCONSUMED = []
    time = [0]
    dist = [0]
    mass = [m_init]
    Comment = []
    Slope = []
    acc = []
    THR = []
    DRAG = []
    config = []
    HLid = []
    LG = []
    BankAngle = []
    ROT = []

    if not AC.BADAFamily.BADAE:
        FUELCONSUMED = [0]

    # BADAH specific
    Peng = []
    Preq = []
    Pav = []

    # optional GPS coordiantes and HDG definition
    if Lat and Lon and (magneticHeading or trueHeading):
        LAT = [Lat]
        LON = [Lon]
        HDGMagnetic = [magneticHeading]
        HDGTrue = [trueHeading]

    else:
        LAT = []
        LON = []
        HDGMagnetic = []
        HDGTrue = []

    # BADAE specific
    Pmec = []
    Pbat = []
    SOCr = []
    SOC = [SOC_init]
    Pelc = []
    Ibat = []
    Vbat = []
    Vgbat = []

    # init loop parameters
    Hp_i = Hp_init

    while True:
        ## PART 1: compute parameters at end of step that are known without uncertainties:
        ##         atmosphere, speeds, ESF

        # atmosphere properties
        H_m = conv.ft2m(Hp_i)  # altitude [m]
        [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
        temp_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)

        # aircraft speed
        [M_i, CAS_i, TAS_i] = atm.convertSpeed(
            v=v, speedType=speedType, theta=theta, delta=delta, sigma=sigma
        )

        # compute Energy Share Factor (ESF)
        ESF_i = AC.esf(
            h=H_m, M=M_i, DeltaTemp=DeltaTemp, flightEvolution=("const" + speedType)
        )

        if turnFlight:
            if turnMetrics["bankAngle"] != 0.0:
                # bankAngle is defined
                rateOfTurn = AC.rateOfTurn_bankAngle(TAS=TAS_i, bankAngle=bankAngle)
            else:
                # rateOfTurn is defined
                bankAngle = AC.bankAngle(rateOfTurn=rateOfTurn, v=TAS_i)  # [degrees]

        # Load factor
        nz = 1 / cos(radians(bankAngle))

        ## PART 2: compute parameters at end of step that are known only with uncertainties
        ##         (due to unknown mass at end of step):
        ##         weight, lift, drag , thrust, fuel flow

        mass_i = mass[-1]
        for _ in itertools.repeat(None, m_iter):
            # BADAH or BADAE
            if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
                # compute Power required for level flight
                Preq_i = AC.Preq(sigma=sigma, tas=TAS_i, mass=mass_i, phi=bankAngle)
                # Compute power required for target ROCD
                Preq_target_i = AC.Peng_target(
                    temp=theta * const.temp_0,
                    DeltaTemp=DeltaTemp,
                    ROCD=ROCDisu,
                    mass=mass_i,
                    Preq=Preq_i,
                    ESF=ESF_i,
                )

            # BADA4
            elif AC.BADAFamily.BADA4:
                # aircraft aerodynamic configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=phase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=phase, previousConfig=config[-1], currentConfig=config_i
                    )

                [HLid_i, LG_i] = AC.flightEnvelope.getAeroConfig(config=config_i)

                # compute lift coefficient
                CL = AC.CL(M=M_i, delta=delta, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(M=M_i, CL=CL, HLid=HLid_i, LG=LG_i, speedBrakes=speedBrakes)
                # compute drag force
                Drag = AC.D(M=M_i, delta=delta, CD=CD)
                # compute thrust force
                THR_i = (
                    ROCDisu * mass_i * const.g * temp_const / (TAS_i * ESF_i) + Drag
                )  # [N]

            # BADA3
            elif AC.BADAFamily.BADA3:
                # aircraft aerodynamic configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=phase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=phase,
                        previousConfig=config[-1],
                        currentConfig=config_i,
                    )

                # compute lift coefficient
                CL = AC.CL(tas=TAS_i, sigma=sigma, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(CL=CL, config=config_i, speedBrakes=speedBrakes)
                # compute drag force
                Drag = AC.D(tas=TAS_i, sigma=sigma, CD=CD)
                # compute thrust force
                THR_i = (
                    ROCDisu * mass_i * const.g * temp_const / (TAS_i * ESF_i) + Drag
                )  # [N]

            # check that required thrust/power fits in the avialable thrust/power envelope,
            # recompute ROCD if necessary and compute fuel flow

            # BADAH
            if AC.BADAFamily.BADAH:
                Pmin = 0.1 * AC.P0  # No minimum power model: assume 10% torque
                Pav_i = AC.Pav(
                    rating="MTKF", theta=theta, delta=delta
                )  # assume MTKF rating as the limit
                Pmax = Pav_i

                if Preq_target_i < Pmin:
                    Preq_target_i = Pmin
                    ROCD_i = (
                        conv.m2ft(
                            AC.ROCD(
                                Peng=Preq_target_i,
                                Preq=Preq_i,
                                mass=mass_i,
                                ESF=ESF_i,
                                theta=theta,
                                DeltaTemp=DeltaTemp,
                            )
                        )
                        * 60
                    )
                elif Preq_target_i > Pmax:
                    Preq_target_i = Pmax
                    ROCD_i = (
                        conv.m2ft(
                            AC.ROCD(
                                Peng=Preq_target_i,
                                Preq=Preq_i,
                                mass=mass_i,
                                ESF=ESF_i,
                                theta=theta,
                                DeltaTemp=DeltaTemp,
                            )
                        )
                        * 60
                    )
                else:
                    ROCD_i = ROCDtarget

                # compute fuel flow for level flight
                CP = AC.CP(Peng=Preq_target_i)
                FUEL_i = AC.ff(delta=delta, CP=CP)  # [kg/s]

            # BADAE
            elif AC.BADAFamily.BADAE:
                Pmin = 0.1 * AC.P0  # No minimum power model: assume 10% torque
                Pav_i = AC.Pav(
                    rating="MTKF", SOC=SOC[-1]
                )  # assume MTKF rating as the limit
                Pmax = Pav_i

                if Preq_target_i < Pmin:
                    Preq_target_i = Pmin
                    ROCD_i = (
                        conv.m2ft(
                            AC.ROCD(
                                Peng=Preq_target_i,
                                Preq=Preq_i,
                                mass=mass_i,
                                ESF=ESF_i,
                                theta=theta,
                                DeltaTemp=DeltaTemp,
                            )
                        )
                        * 60
                    )
                elif Preq_target_i > Pmax:
                    Preq_target_i = Pmax
                    ROCD_i = (
                        conv.m2ft(
                            AC.ROCD(
                                Peng=Preq_target_i,
                                Preq=Preq_i,
                                mass=mass_i,
                                ESF=ESF_i,
                                theta=theta,
                                DeltaTemp=DeltaTemp,
                            )
                        )
                        * 60
                    )
                else:
                    ROCD_i = ROCDtarget

                Pbat_i = AC.Pbat(Preq=Preq_target_i, SOC=SOC[-1])
                SOCr_i = AC.SOCrate(Preq=Preq_target_i, SOC=SOC[-1])

                # debug data
                Pelc_i = Preq_target_i / AC.eta
                Ibat_i = AC.Ibat(P=Pelc_i, SOC=SOC[-1])
                Vbat_i = AC.Vbat(I=Ibat_i, SOC=SOC[-1])
                Vgbat_i = AC.Vocbat(SOC=SOC[-1]) - AC.R0bat(SOC=SOC[-1]) * Ibat_i

            # BADA4
            elif AC.BADAFamily.BADA4:
                THR_min = AC.Thrust(
                    rating="LIDL", delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )  # IDLE Thrust
                FUEL_min = AC.ff(
                    rating="LIDL", delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )  # IDLE Fuel Flow
                THR_max = AC.Thrust(
                    rating="MCMB", delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )  # MCMB Thrust
                FUEL_max = AC.ff(
                    rating="MCMB", delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )  # MCMB Fuel Flow
                if THR_i < THR_min:
                    THR_i = THR_min
                    FUEL_i = FUEL_min
                    ROCD_i = (
                        conv.m2ft(
                            (1 / temp_const)
                            * (THR_i - Drag)
                            * TAS_i
                            * ESF_i
                            / (mass_i * const.g)
                        )
                        * 60
                    )
                elif THR_i > THR_max:
                    THR_i = THR_max
                    FUEL_i = FUEL_max
                    ROCD_i = (
                        conv.m2ft(
                            (1 / temp_const)
                            * (THR_i - Drag)
                            * TAS_i
                            * ESF_i
                            / (mass_i * const.g)
                        )
                        * 60
                    )
                else:
                    CT = AC.CT(Thrust=THR_i, delta=delta)
                    FUEL_i = AC.ff(
                        CT=CT, delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                    )  # [kg/s]
                    ROCD_i = ROCDtarget

            # BADA3
            elif AC.BADAFamily.BADA3:
                THR_min = AC.Thrust(
                    rating="LIDL", v=TAS_i, h=H_m, config="CR", DeltaTemp=DeltaTemp
                )  # IDLE Thrust
                FUEL_min = AC.ff(
                    flightPhase="Descent",
                    v=TAS_i,
                    h=H_m,
                    T=THR_min,
                    config="CR",
                    adapted=False,
                )  # IDLE Fuel Flow
                THR_max = AC.Thrust(
                    rating="MCMB", v=TAS_i, h=H_m, DeltaTemp=DeltaTemp, config="CR"
                )  # MCMB Thrust
                FUEL_max = AC.ff(
                    flightPhase="Climb",
                    v=TAS_i,
                    h=H_m,
                    T=THR_max,
                    config="CR",
                    adapted=False,
                )  # MCMB Fuel Flow

                if THR_i < THR_min:
                    print("below minimum")
                    THR_i = THR_min
                    FUEL_i = FUEL_min
                    ROCD_i = (
                        conv.m2ft(
                            (1 / temp_const)
                            * (THR_i - Drag)
                            * TAS_i
                            * ESF_i
                            / (mass_i * const.g)
                        )
                        * 60
                    )
                elif THR_i > THR_max:
                    THR_i = THR_max
                    FUEL_i = FUEL_max
                    ROCD_i = (
                        conv.m2ft(
                            (1 / temp_const)
                            * (THR_i - Drag)
                            * TAS_i
                            * ESF_i
                            / (mass_i * const.g)
                        )
                        * 60
                    )
                else:
                    FUEL_i = AC.ff(
                        v=TAS_i, h=H_m, T=THR_i, config=config_i, adapted=True
                    )
                    ROCD_i = ROCDtarget

            # Compute elapsed time and fuel burn over current step
            if Hp_i == Hp_init:
                break
            else:
                # Average ROCD over step is the mean of initial and final ones
                step_ROCD = (ROCD[-1] + ROCD_i) / 2  # [ft/min]
                # Step time is: altitude differential divided by average ROCD
                step_time = 60 * (Hp_i - Hp[-1]) / step_ROCD  # [s]

                # BADAE
                if AC.BADAFamily.BADAE:
                    # Average SOC rate over step is the mean of initial and final ones
                    step_SOCr = (SOCr[-1] + SOCr_i) / 2  # [%/h]
                    # SOC change is: average SOC rate multiplied by step time
                    step_SOC = step_SOCr * step_time / 3600  # [%]
                    # Update SOC estimate at end of step
                    SOC_i = SOC[-1] - step_SOC  # [%]
                    # update of aircraft mass estimate at end of step - mass is not changing for ELECTRIC engine (no fuel is consumed)
                    mass_i = mass[-1]  # [kg]

                else:
                    # average fuel flow over step is the mean of initial and final ones
                    step_FUEL = (FUEL[-1] + FUEL_i) / 2  # [kg/s]
                    # fuel burnt is: average fuel flow multiplied by step time
                    step_mass = step_FUEL * step_time  # [kg]
                    # update of aircraft mass estimate at end of step
                    if not mass_const:
                        mass_i = mass[-1] - step_mass  # [kg]
                        fuelConsumed_i = step_FUEL * step_time
                    fuelConsumed_i = FUELCONSUMED[-1] + step_FUEL * step_time

        ## PART 3: store information about end of step point
        writeOutputData = True
        if phase == "Climb" and ROCD_i < 0:
            warnings.warn(
                "Value ROCD = "
                + str(ROCD_i)
                + " [ft/min] is negative at the altitude "
                + str(Hp_i)
                + " [ft].",
                UserWarning,
            )
            go_on = False
            writeOutputData = False

        elif phase == "Climb" and ROCD_i < ROCD_min:
            warnings.warn(
                "Value ROCD = "
                + str(ROCD_i)
                + " [ft/min] exceeds the service ceiling limit defined by minimum ROCD = "
                + str(ROCD_min)
                + " [ft/min] at the altitude "
                + str(Hp_i)
                + " [ft].",
                UserWarning,
            )
            go_on = False
            writeOutputData = True

        if writeOutputData:
            # point data
            Hp.append(Hp_i)
            TAS.append(conv.ms2kt(TAS_i))
            CAS.append(conv.ms2kt(CAS_i))
            M.append(M_i)
            ROCD.append(ROCD_i)
            esf.append(ESF_i)
            Comment.append(comment)

            # everything except electric BADAE
            if not AC.BADAFamily.BADAE:
                FUEL.append(FUEL_i)

            # BADAH
            if AC.BADAFamily.BADAH:
                Peng.append(Preq_target_i)
                Preq.append(Preq_i)
                Pav.append(Pav_i)

            # BADAE
            elif AC.BADAFamily.BADAE:
                Pmec.append(Preq_target_i)
                Pbat.append(Pbat_i)
                SOCr.append(SOCr_i)
                Pelc.append(Pelc_i)
                Ibat.append(Ibat_i)
                Vbat.append(Vbat_i)
                Vgbat.append(Vgbat_i)

            # BADA3 & BADA4
            elif AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
                THR.append(THR_i)
                DRAG.append(Drag)
                config.append(config_i)

            # BADA4
            if AC.BADAFamily.BADA4:
                HLid.append(HLid_i)
                LG.append(LG_i)

            # calculation of the slope
            if TAS_i == 0:
                gamma_i = 90 * np.sign(ROCD_i)
            else:
                if AC.BADAFamily.BADAE:
                    gamma_i = degrees(atan(conv.ft2m(ROCD_i) * temp_const / 60 / TAS_i))
                else:
                    # using SIN assumes the TAS to be in the direction of the aircraft axis, not ground plane. Which means, this should be mathematically the correct equation for all the aircraft
                    gamma_i = degrees(asin(conv.ft2m(ROCD_i) * temp_const / 60 / TAS_i))

            Slope.append(gamma_i)

            # ground speed can be calcualted as TAS projected on the ground minus wind
            GS_i = cos(radians(gamma_i)) * TAS_i - wS
            GS.append(conv.ms2kt(GS_i))

            # integrated data
            if Hp_i != Hp_init:
                if AC.BADAFamily.BADAE:
                    SOC.append(SOC_i)

                mass.append(mass_i)
                time.append(time[-1] + step_time)

                # everything except electric BADAE
                if not AC.BADAFamily.BADAE:
                    FUELCONSUMED.append(fuelConsumed_i)

                # Average TAS over step is the mean of initial and final ones
                step_TAS = (TAS[-2] + TAS[-1]) / 2  # [kt]
                # Average slope over the step
                step_gamma = radians((Slope[-2] + Slope[-1]) / 2)  # radians
                # Average ground speed over step
                # since this is not level flight, TAS speed should be projected on the ground, then GS can be calculated applying the wind speed
                step_TAS_projected = cos(step_gamma) * step_TAS
                step_GS = step_TAS_projected - wS  # [kt]
                # Step distance is: average GS multiplied by step time

                if turnFlight:
                    step_distance = conv.m2nm(
                        turn.distance(
                            rateOfTurn=rateOfTurn, TAS=TAS_i, timeOfTurn=step_time
                        )
                    )  # arcLength during the turn [NM]
                else:
                    step_distance = step_GS * step_time / 3600  # [NM]

                # Distance at end of step is distance at start of step plus step distance
                dist.append(dist[-1] + step_distance)

                # add GPS calculation
                if Lat and Lon and (magneticHeading or trueHeading):
                    if headingToFly == "TRUE":
                        if not turnFlight:
                            if not constantHeading:
                                # fly ORTHODROME
                                (Lat_i, Lon_i, HDGTrue_i) = (
                                    vincenty.destinationPoint_finalBearing(
                                        LAT_init=LAT[-1],
                                        LON_init=LON[-1],
                                        distance=conv.nm2m(step_distance),
                                        bearing=HDGTrue[-1],
                                    )
                                )

                                if magneticDeclinationGrid is not None:
                                    HDGMagnetic_i = (
                                        HDGTrue_i
                                        - magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGMagnetic_i = HDGTrue_i

                            elif constantHeading:
                                # fly LOXODROME
                                (Lat_i, Lon_i) = rhumb.destinationPoint(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearing=HDGTrue[-1],
                                    distance=conv.nm2m(step_distance),
                                )
                                HDGTrue_i = HDGTrue[-1]

                                if magneticDeclinationGrid is not None:
                                    HDGMagnetic_i = (
                                        HDGTrue_i
                                        - magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGMagnetic_i = HDGTrue_i

                        else:
                            # calculate the turn
                            (Lat_i, Lon_i, HDGTrue_i) = (
                                turn.destinationPoint_finalBearing(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearingInit=HDGTrue[-1],
                                    TAS=TAS_i,
                                    rateOfTurn=rateOfTurn,
                                    timeOfTurn=step_time,
                                    directionOfTurn=directionOfTurn,
                                )
                            )

                            if magneticDeclinationGrid is not None:
                                HDGMagnetic_i = (
                                    HDGTrue_i
                                    - magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGMagnetic_i = HDGTrue_i

                    elif headingToFly == "MAGNETIC":
                        if not turnFlight:
                            if constantHeading:
                                (Lat_i, Lon_i) = rhumb.destinationPoint(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearing=HDGTrue[-1],
                                    distance=conv.nm2m(step_distance),
                                )
                                HDGMagnetic_i = HDGMagnetic[-1]

                                if magneticDeclinationGrid is not None:
                                    HDGTrue_i = (
                                        HDGMagnetic_i
                                        + magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGTrue_i = HDGMagnetic_i

                        else:
                            # calculate the turn
                            (Lat_i, Lon_i, HDGTrue_i) = (
                                turn.destinationPoint_finalBearing(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearingInit=HDGTrue[-1],
                                    TAS=TAS_i,
                                    rateOfTurn=rateOfTurn,
                                    timeOfTurn=step_time,
                                    directionOfTurn=directionOfTurn,
                                )
                            )

                            if magneticDeclinationGrid is not None:
                                HDGMagnetic_i = (
                                    HDGTrue_i
                                    - magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGMagnetic_i = HDGTrue_i

                    LAT.append(Lat_i)
                    LON.append(Lon_i)
                    HDGMagnetic.append(HDGMagnetic_i)
                    HDGTrue.append(HDGTrue_i)

            acc.append(0.0)
            BankAngle.append(bankAngle)
            ROT.append(rateOfTurn)

            # Determine end altitude of next step
            Hp_next = Hp_i + Hp_step

            if phase == "Climb":
                if Hp_next < Hp_final:
                    Hp_i = Hp_next - (Hp_i % Hp_step)
                # remaining altitude step would cross over the final altitude
                elif Hp_i < Hp_final:
                    Hp_i = Hp_final
                else:
                    break
            else:
                if Hp_next > Hp_final:
                    Hp_i = Hp_next - (Hp_i % Hp_step)
                # remaining altitude step would cross over the final altitude
                elif Hp_i > Hp_final:
                    Hp_i = Hp_final
                else:
                    break

    flightData = {
        "Hp": Hp,
        "TAS": TAS,
        "CAS": CAS,
        "GS": GS,
        "M": M,
        "acc": acc,
        "ROCD": ROCD,
        "ESF": esf,
        "FUEL": FUEL,
        "Pmec": Pmec,
        "Pelc": Pelc,
        "Pbat": Pbat,
        "SOCr": SOCr,
        "SOC": SOC,
        "Ibat": Ibat,
        "Vbat": Vbat,
        "Vgbat": Vgbat,
        "FUELCONSUMED": FUELCONSUMED,
        "Preq": Preq,
        "Peng": Peng,
        "Pav": Pav,
        "THR": THR,
        "DRAG": DRAG,
        "time": time,
        "dist": dist,
        "slope": Slope,
        "mass": mass,
        "config": config,
        "HLid": HLid,
        "LG": LG,
        "LAT": LAT,
        "LON": LON,
        "HDGTrue": HDGTrue,
        "HDGMagnetic": HDGMagnetic,
        "BankAngle": BankAngle,
        "ROT": ROT,
        "comment": Comment,
    }

    flightTrajectory = FT.createFlightTrajectoryDataframe(flightData)
    return flightTrajectory


def constantSpeedROCD_time(
    AC,
    length,
    speedType,
    v,
    Hp_init,
    ROCDtarget,
    m_init,
    DeltaTemp,
    wS=0.0,
    turnMetrics={"rateOfTurn": 0.0, "bankAngle": 0.0, "directionOfTurn": None},
    Lat=None,
    Lon=None,
    initialHeading={"magnetic": None, "true": None, "constantHeading": None},
    reducedPower=None,
    directionOfTurn=None,
    magneticDeclinationGrid=None,
    **kwargs,
):
    """This function computes time and fuel required by an aircraft to perform a climb/descent from Hp_init for set amount of time at constant speed and constant rate of climb/descent

    :param AC: aircraft {BADA3/4/H/E}
    :param speedType: what kind of speed is followed {M, CAS, TAS}.
    :param length: length of a segment to fly [s]
    :param step_length: length of a step of a segment - [s]
    :param v: what kind of speed is followed - [kt] CAS/TAS speed to follow or [-] MACH speed to follow.
    :param Hp_init: initial pressure altitude [ft].
    :param ROCDtarget: Rate of climb/descent to be followed [ft/min].
    :param m_init: initial aircraft mass [kg].
    :param DeltaTemp: deviation with respect to ISA [K].
    :param wS: longitudinal wind speed (TAS) [kt].
    :param turnMetrics: Metrics for turn performance {"rateOfTurn":0.0,"bankAngle":0.0,"directionOfTurn":None} {[deg/s],[deg],[LEFT/RIGHT]}
    :param SOC_init: initial state of charge [%].
    :param config: aircraft default aerodynamic configuration {TO,IC,CR,AP,LD}.
    :param speedBrakes: deployed or not speedbrakes including value to be added to the drag coeffcient {deployed:False,value:0.03} {deployed:[True/False],value:[-]}.
    :param ROCD_min: lower ROCD threshold to identify the climbing capabilities (service ceiling) [ft/min].
    :param Lat: Geographical Latitude [deg]
    :param Lon: Geographical Longitude [deg]
    :param initialHeading: aircraft magnetic heading, true heading and definition of constant heading(ORTHODROME=False, LOXODROME=True) {[deg],[deg],-}
    :param magneticDeclinationGrid: geographical grid of a magnetic declination on Earth [deg]
    :param mass_const: kind of mass canculation {mass_integrated=False, mass_constant=True}.
    :param m_iter: number of iterations for integration loop [-]
    :param reducedPower: reduction of Power during the climb {True/False}
    :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}.
    :type speedType: string.
    :type length: float.
    :type step_length: float.
    :type v: float.
    :type Hp_init: float.
    :type ROCDtarget: float.
    :type m_init: float.
    :type DeltaTemp: float.
    :type wS: float.
    :type turnMetrics: {float,float,string}.
    :type SOC_init: float.
    :type config: string.
    :type speedBrakes: dict{boolean,float}.
    :type ROCD_min: float.
    :type Lat: float.
    :type Lon: float.
    :type initialHeading: {float,float,boolean}.
    :type magneticDeclinationGrid: magneticDeclinationGrid.
    :type mass_const: boolean.
    :type m_iter: integer.
    :type reducedPower: boolean.
    :returns:
            BADA3: [Hp, TAS, CAS, M, GS, acc, ROCD, ESF, FUEL, FUELCONSUMED, THR, time, dist, slope, mass, config, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, N, s, NM, deg, kg, -, -,deg,deg,deg,deg,deg,deg/s]
            BADA4: [Hp, TAS, CAS, M, GS, acc, ROCD, ESF, FUEL, FUELCONSUMED, THR, P[Pmec, Pbat, Pelc Ibat, Vbat, Vgbat, SOCr, SOC], time, dist, slope, mass, config, HLid, LG, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, N, [W,W,W,A,V,V,%/h,%], s, NM, deg, kg, -, -, -, -,deg,deg,deg,deg,deg,deg/s]
            BADAH: [Hp, TAS, CAS, M, GS, ROCD, ESF, FUEL, FUELCONSUMED, Peng, Preq, Pav, time, dist, slope, mass, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, W, W, W, s, NM, deg, kg, -,deg,deg,deg,deg,deg,deg/s]
            BADAE: [time, dist, Hp, TAS, CAS, M, GS, acc, ROCD, ESF, slope, mass, P[Pmec, Pelc, Pbat, SOCr, SOC, Ibat, Vbat, Vgbat] comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [s, NM, ft, kt, kt, -, kt, m/s^2, ft/min, deg, kg, [W,W,W,%/h,%,A,V,V], -,deg,deg,deg,deg,deg,deg/s]
    :rtype: dict[list[float]}.
    """

    rateOfTurn = turnMetrics["rateOfTurn"]
    bankAngle = turnMetrics["bankAngle"]
    directionOfTurn = turnMetrics["directionOfTurn"]

    turnFlight = True
    if turnMetrics["rateOfTurn"] == 0.0 and turnMetrics["bankAngle"] == 0.0:
        turnFlight = False

    # conversion of Magnetic Heading to True Heading
    if magneticDeclinationGrid is not None:
        magneticDeclination = magneticDeclinationGrid.getMagneticDeclination(
            LAT_target=Lat, LON_target=Lon
        )
    else:
        magneticDeclination = 0

    # retrieve magnetic and true heading inputs
    magneticHeading = initialHeading["magnetic"]
    trueHeading = initialHeading["true"]
    constantHeading = initialHeading["constantHeading"]

    if Lat and Lon and (magneticHeading or trueHeading):
        if trueHeading is not None and magneticHeading is None:
            # fly TRUE Heading
            headingToFly = "TRUE"
            magneticHeading = trueHeading - magneticDeclination
        elif magneticHeading is not None and trueHeading is None:
            # fly MAGNETIC Heading
            if constantHeading == True:
                headingToFly = "MAGNETIC"
                trueHeading = magneticHeading + magneticDeclination
            else:
                raise Exception("Cannot fly non-constant magnetic heading")

        else:
            raise Exception("Undefined Heading value combination")

    # calculation with constant mass (True) or integrated (False)
    mass_const = kwargs.get("mass_const", False)

    # optional parameter to define initial Baterry State of Charge (SOC)
    if AC.BADAFamily.BADAE:
        SOC_init = kwargs.get("SOC_init", 100)
    else:
        SOC_init = None

    # speed brakes application
    if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
        speedBrakes = kwargs.get("speedBrakes", {"deployed": False, "value": 0.03})

    # step size in [s]
    step_length = kwargs.get("step_length", 1)

    # minimum remaining ROCD to determine cruise ceiling
    if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
        ROCD_min = kwargs.get("ROCD_min", 50)  # [ft/min]
    else:
        if AC.engineType == "PISTON" or AC.engineType == "ELECTRIC":
            ROCD_min = kwargs.get("ROCD_min", 100)  # [ft/min]
        else:
            ROCD_min = kwargs.get("ROCD_min", 300)  # [ft/min]

    # check the consistency of ROCD and climb/descent phase of flight
    if ROCDtarget < 0:
        phase = "Descent"
    elif ROCDtarget > 0:
        phase = "Climb"
    else:
        print("ROCDtarget should be different from 0")

    if turnFlight:
        turnComment = "_turn"
    else:
        turnComment = ""

    if constantHeading:
        constHeadingStr = "_const_Heading"
    elif constantHeading is False or constantHeading is None:
        constHeadingStr = ""

    # comment line describing type of trajectory calculation
    comment = phase + turnComment + "_const_ROCD_" + speedType + constHeadingStr

    if Lat and Lon and (magneticHeading or trueHeading):
        comment = comment + "_" + headingToFly + "_Heading"

    # get the default aerodynamic configuration if provided to be used for the whole segment
    config_default = kwargs.get("config", None)
    if config_default is not None:
        if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
            if not (
                config_default == "TO"
                or config_default == "IC"
                or config_default == "CR"
                or config_default == "AP"
                or config_default == "LD"
            ):
                print(
                    "WRONG default configuration set. Available values are: TO/IC/CR/AP/LD. Configuration will be calculated automatically"
                )

    #  weight iteration constant
    m_iter = kwargs.get("m_iter", 5)  # number of iterations for integration loop[-]

    # convert ROCD to IS units
    ROCDisu = conv.ft2m(ROCDtarget) / 60

    # initialize output parameters
    Hp = [Hp_init]
    TAS = []
    CAS = []
    GS = []
    M = []
    ROCD = []
    esf = []
    FUEL = []
    FUELCONSUMED = []
    time = [0]
    dist = [0]
    mass = [m_init]
    Comment = []
    Slope = []
    acc = []
    THR = []
    DRAG = []
    config = []
    HLid = []
    LG = []
    BankAngle = []
    ROT = []

    if not AC.BADAFamily.BADAE:
        FUELCONSUMED = [0]

    # BADAH specific
    Peng = []
    Preq = []
    Pav = []

    # optional GPS coordiantes and HDG definition
    if Lat and Lon and (magneticHeading or trueHeading):
        LAT = [Lat]
        LON = [Lon]
        HDGMagnetic = [magneticHeading]
        HDGTrue = [trueHeading]
    else:
        LAT = []
        LON = []
        HDGMagnetic = []
        HDGTrue = []

    # BADAE specific
    Pmec = []
    Pbat = []
    SOCr = []
    SOC = [SOC_init]
    Pelc = []
    Ibat = []
    Vbat = []
    Vgbat = []

    # init loop parameters
    length_loop = 0
    time_i = time[-1]
    go_on = True

    while go_on:
        ## PART 1: compute parameters at end of step that are known without uncertainties:
        ##         atmosphere, speeds, ESF

        ## PART 2: compute parameters at end of step that are known only with uncertainties
        ##           (due to unknown mass at end of step):
        ##         weight, lift, drag , thrust, fuel flow

        mass_i = mass[-1]
        Hp_i = Hp[-1]

        for _ in itertools.repeat(None, m_iter):
            # atmosphere properties
            H_m = conv.ft2m(Hp_i)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
            temp_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)

            # aircraft speed
            [M_i, CAS_i, TAS_i] = atm.convertSpeed(
                v=v, speedType=speedType, theta=theta, delta=delta, sigma=sigma
            )

            # compute Energy Share Factor (ESF)
            ESF_i = AC.esf(
                h=H_m, M=M_i, DeltaTemp=DeltaTemp, flightEvolution=("const" + speedType)
            )

            if turnFlight:
                if turnMetrics["bankAngle"] != 0.0:
                    # bankAngle is defined
                    rateOfTurn = AC.rateOfTurn_bankAngle(TAS=TAS_i, bankAngle=bankAngle)
                else:
                    # rateOfTurn is defined
                    bankAngle = AC.bankAngle(
                        rateOfTurn=rateOfTurn, v=TAS_i
                    )  # [degrees]

            # Load factor
            nz = 1 / cos(radians(bankAngle))

            step_time = length_loop - time[-1]

            # BADAH or BADAE
            if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
                # compute Power required for level flight
                Preq_i = AC.Preq(sigma=sigma, tas=TAS_i, mass=mass_i, phi=bankAngle)
                # Compute power required for target ROCD
                Preq_target_i = AC.Peng_target(
                    temp=theta * const.temp_0,
                    DeltaTemp=DeltaTemp,
                    ROCD=ROCDisu,
                    mass=mass_i,
                    Preq=Preq_i,
                    ESF=ESF_i,
                )

            # BADA4
            elif AC.BADAFamily.BADA4:
                # aircraft aerodynamic configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=phase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=phase, previousConfig=config[-1], currentConfig=config_i
                    )

                [HLid_i, LG_i] = AC.flightEnvelope.getAeroConfig(config=config_i)

                # compute lift coefficient
                CL = AC.CL(M=M_i, delta=delta, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(M=M_i, CL=CL, HLid=HLid_i, LG=LG_i, speedBrakes=speedBrakes)
                # compute drag force
                Drag = AC.D(M=M_i, delta=delta, CD=CD)
                # compute thrust force
                THR_i = (
                    ROCDisu * mass_i * const.g * temp_const / (TAS_i * ESF_i) + Drag
                )  # [N]

            # BADA3
            elif AC.BADAFamily.BADA3:
                # aircraft aerodynamic configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=phase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=phase,
                        previousConfig=config[-1],
                        currentConfig=config_i,
                    )

                # compute lift coefficient
                CL = AC.CL(tas=TAS_i, sigma=sigma, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(CL=CL, config=config_i, speedBrakes=speedBrakes)
                # compute drag force
                Drag = AC.D(tas=TAS_i, sigma=sigma, CD=CD)
                # compute thrust force
                THR_i = (
                    ROCDisu * mass_i * const.g * temp_const / (TAS_i * ESF_i) + Drag
                )  # [N]

            # check that required thrust/power fits in the avialable thrust/power envelope,
            # recompute ROCD if necessary and compute fuel flow

            # BADAH
            if AC.BADAFamily.BADAH:
                Pmin = 0.1 * AC.P0  # No minimum power model: assume 10% torque
                Pav_i = AC.Pav(
                    rating="MTKF", theta=theta, delta=delta
                )  # assume MTKF rating as the limit
                Pmax = Pav_i

                if Preq_target_i < Pmin:
                    Preq_target_i = Pmin
                    ROCD_i = (
                        conv.m2ft(
                            AC.ROCD(
                                Peng=Preq_target_i,
                                Preq=Preq_i,
                                mass=mass_i,
                                ESF=ESF_i,
                                theta=theta,
                                DeltaTemp=DeltaTemp,
                            )
                        )
                        * 60
                    )
                elif Preq_target_i > Pmax:
                    Preq_target_i = Pmax
                    ROCD_i = (
                        conv.m2ft(
                            AC.ROCD(
                                Peng=Preq_target_i,
                                Preq=Preq_i,
                                mass=mass_i,
                                ESF=ESF_i,
                                theta=theta,
                                DeltaTemp=DeltaTemp,
                            )
                        )
                        * 60
                    )
                else:
                    ROCD_i = ROCDtarget

                # compute fuel flow for level flight
                CP = AC.CP(Peng=Preq_target_i)
                FUEL_i = AC.ff(delta=delta, CP=CP)  # [kg/s]

            # BADAE
            elif AC.BADAFamily.BADAE:
                Pmin = 0.1 * AC.P0  # No minimum power model: assume 10% torque
                Pav_i = AC.Pav(
                    rating="MTKF", theta=theta, delta=delta
                )  # assume MTKF rating as the limit
                Pmax = Pav_i

                if Preq_target_i < Pmin:
                    Preq_target_i = Pmin
                    ROCD_i = (
                        conv.m2ft(
                            AC.ROCD(
                                Peng=Preq_target_i,
                                Preq=Preq_i,
                                mass=mass_i,
                                ESF=ESF_i,
                                theta=theta,
                                DeltaTemp=DeltaTemp,
                            )
                        )
                        * 60
                    )
                elif Preq_target_i > Pmax:
                    Preq_target_i = Pmax
                    ROCD_i = (
                        conv.m2ft(
                            AC.ROCD(
                                Peng=Preq_target_i,
                                Preq=Preq_i,
                                mass=mass_i,
                                ESF=ESF_i,
                                theta=theta,
                                DeltaTemp=DeltaTemp,
                            )
                        )
                        * 60
                    )
                else:
                    ROCD_i = ROCDtarget

                Pbat_i = AC.Pbat(Preq=Preq_target_i, SOC=SOC[-1])
                SOCr_i = AC.SOCrate(Preq=Preq_target_i, SOC=SOC[-1])

                # debug data
                Pelc_i = Preq_target_i / AC.eta
                Ibat_i = AC.Ibat(P=Pelc_i, SOC=SOC[-1])
                Vbat_i = AC.Vbat(I=Ibat_i, SOC=SOC[-1])
                Vgbat_i = AC.Vocbat(SOC=SOC[-1]) - AC.R0bat(SOC=SOC[-1]) * Ibat_i

            # BADA4
            elif AC.BADAFamily.BADA4:
                THR_min = AC.Thrust(
                    rating="LIDL", delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )  # IDLE Thrust
                FUEL_min = AC.ff(
                    rating="LIDL", delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )  # IDLE Fuel Flow
                THR_max = AC.Thrust(
                    rating="MCMB", delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )  # MCMB Thrust
                FUEL_max = AC.ff(
                    rating="MCMB", delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )  # MCMB Fuel Flow
                if THR_i < THR_min:
                    THR_i = THR_min
                    FUEL_i = FUEL_min
                    ROCD_i = (
                        conv.m2ft(
                            (1 / temp_const)
                            * (THR_i - Drag)
                            * TAS_i
                            * ESF_i
                            / (mass_i * const.g)
                        )
                        * 60
                    )
                elif THR_i > THR_max:
                    THR_i = THR_max
                    FUEL_i = FUEL_max
                    ROCD_i = (
                        conv.m2ft(
                            (1 / temp_const)
                            * (THR_i - Drag)
                            * TAS_i
                            * ESF_i
                            / (mass_i * const.g)
                        )
                        * 60
                    )
                else:
                    CT = AC.CT(Thrust=THR_i, delta=delta)
                    FUEL_i = AC.ff(
                        CT=CT, delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                    )  # [kg/s]
                    ROCD_i = ROCDtarget

            # BADA3
            elif AC.BADAFamily.BADA3:
                THR_min = AC.Thrust(
                    rating="LIDL", v=TAS_i, h=H_m, config="CR", DeltaTemp=DeltaTemp
                )  # IDLE Thrust
                FUEL_min = AC.ff(
                    flightPhase="Descent",
                    v=TAS_i,
                    h=H_m,
                    T=THR_min,
                    config="CR",
                    adapted=False,
                )  # IDLE Fuel Flow
                THR_max = AC.Thrust(
                    rating="MCMB", v=TAS_i, h=H_m, DeltaTemp=DeltaTemp, config="CR"
                )  # MCMB Thrust
                FUEL_max = AC.ff(
                    flightPhase="Climb",
                    v=TAS_i,
                    h=H_m,
                    T=THR_max,
                    config="CR",
                    adapted=False,
                )  # MCMB Fuel Flow

                if THR_i < THR_min:
                    THR_i = THR_min
                    FUEL_i = FUEL_min
                    ROCD_i = (
                        conv.m2ft(
                            (1 / temp_const)
                            * (THR_i - Drag)
                            * TAS_i
                            * ESF_i
                            / (mass_i * const.g)
                        )
                        * 60
                    )
                elif THR_i > THR_max:
                    THR_i = THR_max
                    FUEL_i = FUEL_max
                    ROCD_i = (
                        conv.m2ft(
                            (1 / temp_const)
                            * (THR_i - Drag)
                            * TAS_i
                            * ESF_i
                            / (mass_i * const.g)
                        )
                        * 60
                    )
                else:
                    FUEL_i = AC.ff(
                        v=TAS_i, h=H_m, T=THR_i, config=config_i, adapted=True
                    )
                    ROCD_i = ROCDtarget

            # Compute elapsed time, altitude and fuel burn over current step
            if length_loop == 0:
                # no need to loop for first point: initial m/Hp already known
                break
            else:
                # Average ROCD over step is the mean of initial and final ones
                step_ROCD = (ROCD[-1] + ROCD_i) / 2  # [ft/min]
                # Altitude differential is: average ROCD multiplied by step time
                step_Hp = step_ROCD * step_time / 60  # [ft]
                # Update altitude estimate at end of step
                Hp_i = Hp[-1] + step_Hp  # [ft]

                # BADAE
                if AC.BADAFamily.BADAE:
                    # Average SOC rate over step is the mean of initial and final ones
                    step_SOCr = (SOCr[-1] + SOCr_i) / 2  # [%/h]
                    # SOC change is: average SOC rate multiplied by step time
                    step_SOC = step_SOCr * step_time / 3600  # [%]
                    # Update SOC estimate at end of step
                    SOC_i = SOC[-1] - step_SOC  # [%]
                    # update of aircraft mass estimate at end of step - mass is not changing for ELECTRIC engine (no fuel is consumed)
                    mass_i = mass[-1]  # [kg]

                else:
                    # average fuel flow over step is the mean of initial and final ones
                    step_FUEL = (FUEL[-1] + FUEL_i) / 2  # [kg/s]
                    # fuel burnt is: average fuel flow multiplied by step time
                    step_mass = step_FUEL * step_time  # [kg]
                    # update of aircraft mass estimate at end of step
                    if not mass_const:
                        mass_i = mass[-1] - step_mass  # [kg]
                        fuelConsumed_i = step_FUEL * step_time
                    fuelConsumed_i = FUELCONSUMED[-1] + step_FUEL * step_time

        ## PART 3: store information about end of step point
        writeOutputData = True
        if phase == "Climb" and ROCD_i < 0:
            warnings.warn(
                "Value ROCD = "
                + str(ROCD_i)
                + " [ft/min] is negative at the altitude "
                + str(Hp_i)
                + " [ft].",
                UserWarning,
            )
            go_on = False
            writeOutputData = False

        elif phase == "Climb" and ROCD_i < ROCD_min:
            warnings.warn(
                "Value ROCD = "
                + str(ROCD_i)
                + " [ft/min] exceeds the service ceiling limit defined by minimum ROCD = "
                + str(ROCD_min)
                + " [ft/min] at the altitude "
                + str(Hp_i)
                + " [ft].",
                UserWarning,
            )
            go_on = False
            writeOutputData = True

        if writeOutputData:
            # point data
            TAS.append(conv.ms2kt(TAS_i))
            CAS.append(conv.ms2kt(CAS_i))
            M.append(M_i)
            ROCD.append(ROCD_i)
            esf.append(ESF_i)
            Comment.append(comment)

            # everything except electric BADAE
            if not AC.BADAFamily.BADAE:
                FUEL.append(FUEL_i)

            # BADAH
            if AC.BADAFamily.BADAH:
                Peng.append(Preq_target_i)
                Preq.append(Preq_i)
                Pav.append(Pav_i)

            # BADAE
            elif AC.BADAFamily.BADAE:
                Pmec.append(Preq_target_i)
                Pbat.append(Pbat_i)
                SOCr.append(SOCr_i)
                Pelc.append(Pelc_i)
                Ibat.append(Ibat_i)
                Vbat.append(Vbat_i)
                Vgbat.append(Vgbat_i)

            # BADA3 & BADA4
            elif AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
                THR.append(THR_i)
                DRAG.append(Drag)
                config.append(config_i)

            # BADA4
            if AC.BADAFamily.BADA4:
                HLid.append(HLid_i)
                LG.append(LG_i)

            # calculation of the slope
            if TAS_i == 0:
                gamma_i = 90 * np.sign(ROCD_i)
            else:
                [theta, delta, sigma] = atm.atmosphereProperties(
                    h=conv.ft2m(Hp_i), DeltaTemp=DeltaTemp
                )
                temp_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)
                if AC.BADAFamily.BADAE:
                    gamma_i = degrees(atan(conv.ft2m(ROCD_i) * temp_const / 60 / TAS_i))
                else:
                    # using SIN assumes the TAS to be in the direction of the aircraft axis, not ground plane. Which means, this should be mathematically the correct equation for all the aircraft
                    gamma_i = degrees(asin(conv.ft2m(ROCD_i) * temp_const / 60 / TAS_i))

            # ground speed can be calcualted as TAS projected on the ground minus wind
            GS_i = cos(radians(gamma_i)) * TAS_i - wS
            GS.append(conv.ms2kt(GS_i))

            Slope.append(gamma_i)
            acc.append(0.0)

            BankAngle.append(bankAngle)
            ROT.append(rateOfTurn)

            # integrated data
            if length_loop != 0:
                if AC.BADAFamily.BADAE:
                    SOC.append(SOC_i)

                Hp.append(Hp_i)
                mass.append(mass_i)
                time.append(time[-1] + step_time)

                # everything except electric BADAE
                if not AC.BADAFamily.BADAE:
                    FUELCONSUMED.append(fuelConsumed_i)

                # Average TAS over step is the mean of initial and final ones
                step_TAS = (TAS[-2] + TAS[-1]) / 2  # [kt]
                # Average slope over the step
                step_gamma = radians((Slope[-2] + Slope[-1]) / 2)  # radians
                # Average ground speed over step
                # since this is not level flight, TAS speed should be projected on the ground, then GS can be calculated applying the wind speed
                step_TAS_projected = cos(step_gamma) * step_TAS
                step_GS = step_TAS_projected - wS  # [kt]
                # Step distance is: average GS multiplied by step time
                if turnFlight:
                    step_distance = conv.m2nm(
                        turn.distance(
                            rateOfTurn=rateOfTurn, TAS=TAS_i, timeOfTurn=step_time
                        )
                    )  # arcLength during the turn [NM]
                else:
                    step_distance = step_GS * step_time / 3600  # [NM]
                # Distance at end of step is distance at start of step plus step distance
                dist.append(dist[-1] + step_distance)

                # add GPS calculation
                if Lat and Lon and (magneticHeading or trueHeading):
                    if headingToFly == "TRUE":
                        if not turnFlight:
                            if not constantHeading:
                                # fly ORTHODROME
                                (Lat_i, Lon_i, HDGTrue_i) = (
                                    vincenty.destinationPoint_finalBearing(
                                        LAT_init=LAT[-1],
                                        LON_init=LON[-1],
                                        distance=conv.nm2m(step_distance),
                                        bearing=HDGTrue[-1],
                                    )
                                )

                                if magneticDeclinationGrid is not None:
                                    HDGMagnetic_i = (
                                        HDGTrue_i
                                        - magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGMagnetic_i = HDGTrue_i

                            elif constantHeading:
                                # fly LOXODROME
                                (Lat_i, Lon_i) = rhumb.destinationPoint(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearing=HDGTrue[-1],
                                    distance=conv.nm2m(step_distance),
                                )
                                HDGTrue_i = HDGTrue[-1]

                                if magneticDeclinationGrid is not None:
                                    HDGMagnetic_i = (
                                        HDGTrue_i
                                        - magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGMagnetic_i = HDGTrue_i

                        else:
                            # calculate the turn
                            (Lat_i, Lon_i, HDGTrue_i) = (
                                turn.destinationPoint_finalBearing(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearingInit=HDGTrue[-1],
                                    TAS=TAS_i,
                                    rateOfTurn=rateOfTurn,
                                    timeOfTurn=step_time,
                                    directionOfTurn=directionOfTurn,
                                )
                            )

                            if magneticDeclinationGrid is not None:
                                HDGMagnetic_i = (
                                    HDGTrue_i
                                    - magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGMagnetic_i = HDGTrue_i

                    elif headingToFly == "MAGNETIC":
                        if not turnFlight:
                            if constantHeading:
                                (Lat_i, Lon_i) = rhumb.destinationPoint(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearing=HDGTrue[-1],
                                    distance=conv.nm2m(step_distance),
                                )
                                HDGMagnetic_i = HDGMagnetic[-1]

                                if magneticDeclinationGrid is not None:
                                    HDGTrue_i = (
                                        HDGMagnetic_i
                                        + magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGTrue_i = HDGMagnetic_i

                        else:
                            # calculate the turn
                            (Lat_i, Lon_i, HDGTrue_i) = (
                                turn.destinationPoint_finalBearing(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearingInit=HDGTrue[-1],
                                    TAS=TAS_i,
                                    rateOfTurn=rateOfTurn,
                                    timeOfTurn=step_time,
                                    directionOfTurn=directionOfTurn,
                                )
                            )

                            if magneticDeclinationGrid is not None:
                                HDGMagnetic_i = (
                                    HDGTrue_i
                                    - magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGMagnetic_i = HDGTrue_i

                    LAT.append(Lat_i)
                    LON.append(Lon_i)
                    HDGMagnetic.append(HDGMagnetic_i)
                    HDGTrue.append(HDGTrue_i)

            if length_loop + step_length < length:
                length_loop += step_length
            elif length_loop < length:
                length_loop = length
            else:
                go_on = False

    flightData = {
        "Hp": Hp,
        "TAS": TAS,
        "CAS": CAS,
        "GS": GS,
        "M": M,
        "acc": acc,
        "ROCD": ROCD,
        "ESF": esf,
        "FUEL": FUEL,
        "Pmec": Pmec,
        "Pelc": Pelc,
        "Pbat": Pbat,
        "SOCr": SOCr,
        "SOC": SOC,
        "Ibat": Ibat,
        "Vbat": Vbat,
        "Vgbat": Vgbat,
        "FUELCONSUMED": FUELCONSUMED,
        "Preq": Preq,
        "Peng": Peng,
        "Pav": Pav,
        "THR": THR,
        "DRAG": DRAG,
        "time": time,
        "dist": dist,
        "slope": Slope,
        "mass": mass,
        "config": config,
        "HLid": HLid,
        "LG": LG,
        "LAT": LAT,
        "LON": LON,
        "HDGTrue": HDGTrue,
        "HDGMagnetic": HDGMagnetic,
        "BankAngle": BankAngle,
        "ROT": ROT,
        "comment": Comment,
    }

    flightTrajectory = FT.createFlightTrajectoryDataframe(flightData)
    return flightTrajectory


def constantSpeedSlope(
    AC,
    speedType,
    v,
    Hp_init,
    Hp_final,
    slopetarget,
    m_init,
    DeltaTemp,
    wS=0.0,
    turnMetrics={"rateOfTurn": 0.0, "bankAngle": 0.0, "directionOfTurn": None},
    Lat=None,
    Lon=None,
    initialHeading={"magnetic": None, "true": None, "constantHeading": None},
    reducedPower=None,
    directionOfTurn=None,
    magneticDeclinationGrid=None,
    **kwargs,
):
    """This function computes time and fuel required by an aircraft to perform a climb/descent from Hp_init to Hp_final at constant speed and constant slope

    :param AC: aircraft {BADA3/4/H/E}
    :param speedType: what kind of speed is followed {M, CAS, TAS}.
    :param v: what kind of speed is followed - [kt] CAS/TAS speed to follow or [-] MACH speed to follow.
    :param Hp_init: initial pressure altitude [ft].
    :param Hp_final: final pressure altitude [ft].
    :param slopetarget: slope to be followed [deg].
    :param m_init: initial aircraft mass [kg].
    :param DeltaTemp: deviation with respect to ISA [K].
    :param wS: longitudinal wind speed (TAS) [kt].
    :param turnMetrics: Metrics for turn performance {"rateOfTurn":0.0,"bankAngle":0.0,"directionOfTurn":None} {[deg/s],[deg],[LEFT/RIGHT]}
    :param Hp_step: length of an altitude step of a segment [ft].
    :param SOC_init: initial state of charge [%].
    :param config: aircraft default aerodynamic configuration {TO,IC,CR,AP,LD}.
    :param speedBrakes: deployed or not speedbrakes including value to be added to the drag coeffcient {deployed:False,value:0.03} {deployed:[True/False],value:[-]}.
    :param ROCD_min: lower ROCD threshold to identify the climbing capabilities (service ceiling) [ft/min].
    :param Lat: Geographical Latitude [deg]
    :param Lon: Geographical Longitude [deg]
    :param initialHeading: aircraft magnetic heading, true heading and definition of constant heading(ORTHODROME=False, LOXODROME=True) {[deg],[deg],-}
    :param magneticDeclinationGrid: geographical grid of a magnetic declination on Earth [deg]
    :param mass_const: kind of mass canculation {mass_integrated=False, mass_constant=True}.
    :param m_iter: number of iterations for integration loop [-]
    :param reducedPower: reduction of Power during the climb {True/False}
    :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}.
    :type speedType: string.
    :type v: float.
    :type Hp_init: float.
    :type Hp_final: float.
    :type slopetarget: float.
    :type m_init: float.
    :type DeltaTemp: float.
    :type wS: float.
    :type turnMetrics: {float,float,string}.
    :type Hp_step: float.
    :type SOC_init: float.
    :type config: string.
    :type speedBrakes: dict{boolean,float}.
    :type ROCD_min: float.
    :type Lat: float.
    :type Lon: float.
    :type initialHeading: {float,float,boolean}.
    :type magneticDeclinationGrid: magneticDeclinationGrid.
    :type mass_const: boolean.
    :type m_iter: integer.
    :type reducedPower: boolean.
    :returns:
            BADA3: [Hp, TAS, CAS, M, GS, acc, ROCD, ESF, FUEL, FUELCONSUMED, THR, time, dist, slope, mass, config, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, N, s, NM, deg, kg, -, -,deg,deg,deg,deg,deg,deg/s]
            BADA4: [Hp, TAS, CAS, M, GS, acc, ROCD, ESF, FUEL, FUELCONSUMED, THR, P[Pmec, Pbat, Pelc Ibat, Vbat, Vgbat, SOCr, SOC], time, dist, slope, mass, config, HLid, LG, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, N, [W,W,W,A,V,V,%/h,%], s, NM, deg, kg, -, -, -, -,deg,deg,deg,deg,deg,deg/s]
            BADAH: [Hp, TAS, CAS, M, GS, ROCD, ESF, FUEL, FUELCONSUMED, Peng, Preq, Pav, time, dist, slope, mass, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, W, W, W, s, NM, deg, kg, -,deg,deg,deg,deg,deg,deg/s]
            BADAE: [time, dist, Hp, TAS, CAS, M, GS, acc, ROCD, ESF, slope, mass, P[Pmec, Pelc, Pbat, SOCr, SOC, Ibat, Vbat, Vgbat] comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [s, NM, ft, kt, kt, -, kt, m/s^2, ft/min, deg, kg, [W,W,W,%/h,%,A,V,V], -,deg,deg,deg,deg,deg,deg/s]
    :rtype: dict[list[float]}.
    """

    rateOfTurn = turnMetrics["rateOfTurn"]
    bankAngle = turnMetrics["bankAngle"]
    directionOfTurn = turnMetrics["directionOfTurn"]

    turnFlight = True
    if turnMetrics["rateOfTurn"] == 0.0 and turnMetrics["bankAngle"] == 0.0:
        turnFlight = False

    # conversion of Magnetic Heading to True Heading
    if magneticDeclinationGrid is not None:
        magneticDeclination = magneticDeclinationGrid.getMagneticDeclination(
            LAT_target=Lat, LON_target=Lon
        )
    else:
        magneticDeclination = 0

    # retrieve magnetic and true heading inputs
    magneticHeading = initialHeading["magnetic"]
    trueHeading = initialHeading["true"]
    constantHeading = initialHeading["constantHeading"]

    if Lat and Lon and (magneticHeading or trueHeading):
        if trueHeading is not None and magneticHeading is None:
            # fly TRUE Heading
            headingToFly = "TRUE"
            magneticHeading = trueHeading - magneticDeclination
        elif magneticHeading is not None and trueHeading is None:
            # fly MAGNETIC Heading
            if constantHeading == True:
                headingToFly = "MAGNETIC"
                trueHeading = magneticHeading + magneticDeclination
            else:
                raise Exception("Cannot fly non-constant magnetic heading")

        else:
            raise Exception("Undefined Heading value combination")

    # calculation with constant mass (True) or integrated (False)
    mass_const = kwargs.get("mass_const", False)

    # optional parameter to define initial Baterry State of Charge (SOC)
    if AC.BADAFamily.BADAE:
        SOC_init = kwargs.get("SOC_init", 100)
    else:
        SOC_init = None

    # speed brakes application
    if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
        speedBrakes = kwargs.get("speedBrakes", {"deployed": False, "value": 0.03})

    # optional parameter - iteration step for altitude loop
    if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
        Hp_step = kwargs.get("Hp_step", 500)  # [ft]
    else:
        # NB: it must be a multiple of 1000ft so that interrupted climbs end on a regular cruise altitude.
        Hp_step = kwargs.get("Hp_step", 1000)  # [ft]

    # minimum remaining ROCD to determine cruise ceiling
    if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
        ROCD_min = kwargs.get("ROCD_min", 50)  # [ft/min]
    else:
        if AC.engineType == "PISTON" or AC.engineType == "ELECTRIC":
            ROCD_min = kwargs.get("ROCD_min", 100)  # [ft/min]
        else:
            ROCD_min = kwargs.get("ROCD_min", 300)  # [ft/min]

    # determine if  vertical evolution over the segment is CLIMB or DESCENT
    if Hp_init < Hp_final:
        phase = "Climb"
    else:
        phase = "Descent"
        Hp_step = -Hp_step

    # check the consistency of SLOPE and climb/descent phase of flight
    # if incosistent, change the sign on slope target value
    if phase == "Climb" and slopetarget < 0:
        slopetarget = abs(slopetarget)
        print("Slopetarget for Climb should be positive")
    elif phase == "Descent" and slopetarget > 0:
        slopetarget = slopetarget * (-1)
        print("Slopetarget for Descent should be negative")

    if turnFlight:
        turnComment = "_turn"
    else:
        turnComment = ""

    if constantHeading:
        constHeadingStr = "_const_Heading"
    elif constantHeading is False or constantHeading is None:
        constHeadingStr = ""

    # comment line describing type of trajectory calculation
    comment = phase + turnComment + "_const_Slope_" + speedType + constHeadingStr

    if Lat and Lon and (magneticHeading or trueHeading):
        comment = comment + "_" + headingToFly + "_Heading"

    # get the default aerodynamic configuration if provided to be used for the whole segment
    config_default = kwargs.get("config", None)
    if config_default is not None:
        if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
            if not (
                config_default == "TO"
                or config_default == "IC"
                or config_default == "CR"
                or config_default == "AP"
                or config_default == "LD"
            ):
                print(
                    "WRONG default configuration set. Available values are: TO/IC/CR/AP/LD. Configuration will be calculated automatically"
                )

    #  weight iteration constant
    m_iter = kwargs.get("m_iter", 5)  # number of iterations for integration loop[-]

    # initialize output parameters
    Hp = []
    TAS = []
    CAS = []
    GS = []
    M = []
    ROCD = []
    esf = []
    FUEL = []
    FUELCONSUMED = []
    time = [0]
    dist = [0]
    mass = [m_init]
    Comment = []
    Slope = []
    acc = []
    THR = []
    DRAG = []
    config = []
    HLid = []
    LG = []
    BankAngle = []
    ROT = []

    if not AC.BADAFamily.BADAE:
        FUELCONSUMED = [0]

    # BADAH specific
    Peng = []
    Preq = []
    Pav = []

    # optional GPS coordiantes and HDG definition
    if Lat and Lon and (magneticHeading or trueHeading):
        LAT = [Lat]
        LON = [Lon]
        HDGMagnetic = [magneticHeading]
        HDGTrue = [trueHeading]
    else:
        LAT = []
        LON = []
        HDGMagnetic = []
        HDGTrue = []

    # BADAE specific
    Pmec = []
    Pbat = []
    SOCr = []
    SOC = [SOC_init]
    Pelc = []
    Ibat = []
    Vbat = []
    Vgbat = []

    # init loop parameters
    Hp_i = Hp_init
    go_on = True

    while go_on:
        ## PART 1: compute parameters at end of step that are known without uncertainties:
        ##         atmosphere, speeds, ESF

        # atmosphere properties
        H_m = conv.ft2m(Hp_i)  # altitude [m]
        [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
        temp_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)

        # aircraft speed
        [M_i, CAS_i, TAS_i] = atm.convertSpeed(
            v=v, speedType=speedType, theta=theta, delta=delta, sigma=sigma
        )

        if turnFlight:
            if turnMetrics["bankAngle"] != 0.0:
                # bankAngle is defined
                rateOfTurn = AC.rateOfTurn_bankAngle(TAS=TAS_i, bankAngle=bankAngle)
            else:
                # rateOfTurn is defined
                bankAngle = AC.bankAngle(rateOfTurn=rateOfTurn, v=TAS_i)  # [degrees]

        # Load factor
        nz = 1 / cos(radians(bankAngle))

        # compute Energy Share Factor (ESF)
        ESF_i = AC.esf(
            h=H_m, M=M_i, DeltaTemp=DeltaTemp, flightEvolution=("const" + speedType)
        )

        if AC.BADAFamily.BADAE:
            # special case for BADAE, in future it may apply also for BADAH
            ROCDisu = tan(conv.deg2rad(slopetarget)) * TAS_i * (1 / temp_const)
        else:
            ROCDisu = sin(conv.deg2rad(slopetarget)) * TAS_i * (1 / temp_const)

        ## PART 2: compute parameters at end of step that are known only with uncertainties
        ##           (due to unknown mass at end of step):
        ##         weight, lift, drag , thrust, fuel flow

        mass_i = mass[-1]
        for _ in itertools.repeat(None, m_iter):
            # BADAH or BADAE
            if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
                # compute Power required for level flight
                Preq_i = AC.Preq(sigma=sigma, tas=TAS_i, mass=mass_i, phi=bankAngle)
                # Compute power required for target ROCD
                Preq_target_i = AC.Peng_target(
                    temp=theta * const.temp_0,
                    DeltaTemp=DeltaTemp,
                    ROCD=ROCDisu,
                    mass=mass_i,
                    Preq=Preq_i,
                    ESF=ESF_i,
                )

            # BADA4
            elif AC.BADAFamily.BADA4:
                # aircraft aerodynamic configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=phase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=phase, previousConfig=config[-1], currentConfig=config_i
                    )

                [HLid_i, LG_i] = AC.flightEnvelope.getAeroConfig(config=config_i)

                # compute lift coefficient
                CL = AC.CL(M=M_i, delta=delta, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(M=M_i, CL=CL, HLid=HLid_i, LG=LG_i, speedBrakes=speedBrakes)
                # compute drag force
                Drag = AC.D(M=M_i, delta=delta, CD=CD)
                # compute thrust force
                THR_i = (
                    ROCDisu * mass_i * const.g * temp_const / (TAS_i * ESF_i) + Drag
                )  # [N]

            # BADA3
            elif AC.BADAFamily.BADA3:
                # aircraft aerodynamic configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=phase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=phase,
                        previousConfig=config[-1],
                        currentConfig=config_i,
                    )

                # compute lift coefficient
                CL = AC.CL(tas=TAS_i, sigma=sigma, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(CL=CL, config=config_i, speedBrakes=speedBrakes)
                # compute drag force
                Drag = AC.D(tas=TAS_i, sigma=sigma, CD=CD)
                # compute thrust force
                THR_i = (
                    ROCDisu * mass_i * const.g * temp_const / (TAS_i * ESF_i) + Drag
                )  # [N]

            # check that required thrust/power fits in the avialable thrust/power envelope,
            # recompute ROCD if necessary and compute fuel flow

            # BADAH
            if AC.BADAFamily.BADAH:
                Pmin = 0.1 * AC.P0  # No minimum power model: assume 10% torque
                Pav_i = AC.Pav(
                    rating="MTKF", theta=theta, delta=delta
                )  # assume MTKF rating as the limit
                Pmax = Pav_i

                if Preq_target_i < Pmin:
                    Preq_target_i = Pmin
                    ROCD_i = (
                        conv.m2ft(
                            AC.ROCD(
                                Peng=Preq_target_i,
                                Preq=Preq_i,
                                mass=mass_i,
                                ESF=ESF_i,
                                theta=theta,
                                DeltaTemp=DeltaTemp,
                            )
                        )
                        * 60
                    )
                elif Preq_target_i > Pmax:
                    Preq_target_i = Pmax
                    ROCD_i = (
                        conv.m2ft(
                            AC.ROCD(
                                Peng=Preq_target_i,
                                Preq=Preq_i,
                                mass=mass_i,
                                ESF=ESF_i,
                                theta=theta,
                                DeltaTemp=DeltaTemp,
                            )
                        )
                        * 60
                    )
                else:
                    ROCD_i = conv.m2ft(ROCDisu) * 60

                # compute fuel flow for level flight
                CP = AC.CP(Peng=Preq_target_i)
                FUEL_i = AC.ff(delta=delta, CP=CP)  # [kg/s]

            # BADAE
            elif AC.BADAFamily.BADAE:
                Pmin = 0.1 * AC.P0  # No minimum power model: assume 10% torque
                Pav_i = AC.Pav(
                    rating="MTKF", SOC=SOC[-1]
                )  # assume MTKF rating as the limit
                Pmax = Pav_i

                if Preq_target_i < Pmin:
                    Preq_target_i = Pmin
                    ROCD_i = (
                        conv.m2ft(
                            AC.ROCD(
                                Peng=Preq_target_i,
                                Preq=Preq_i,
                                mass=mass_i,
                                ESF=ESF_i,
                                theta=theta,
                                DeltaTemp=DeltaTemp,
                            )
                        )
                        * 60
                    )
                elif Preq_target_i > Pmax:
                    Preq_target_i = Pmax
                    ROCD_i = (
                        conv.m2ft(
                            AC.ROCD(
                                Peng=Preq_target_i,
                                Preq=Preq_i,
                                mass=mass_i,
                                ESF=ESF_i,
                                theta=theta,
                                DeltaTemp=DeltaTemp,
                            )
                        )
                        * 60
                    )
                else:
                    ROCD_i = conv.m2ft(ROCDisu) * 60

                Pbat_i = AC.Pbat(Preq=Preq_target_i, SOC=SOC[-1])
                SOCr_i = AC.SOCrate(Preq=Preq_target_i, SOC=SOC[-1])

                # debug data
                Pelc_i = Preq_target_i / AC.eta
                Ibat_i = AC.Ibat(P=Pelc_i, SOC=SOC[-1])
                Vbat_i = AC.Vbat(I=Ibat_i, SOC=SOC[-1])
                Vgbat_i = AC.Vocbat(SOC=SOC[-1]) - AC.R0bat(SOC=SOC[-1]) * Ibat_i

            # BADA4
            elif AC.BADAFamily.BADA4:
                THR_min = AC.Thrust(
                    rating="LIDL", delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )  # IDLE Thrust
                FUEL_min = AC.ff(
                    rating="LIDL", delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )  # IDLE Fuel Flow
                THR_max = AC.Thrust(
                    rating="MCMB", delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )  # MCMB Thrust
                FUEL_max = AC.ff(
                    rating="MCMB", delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )  # MCMB Fuel Flow

                if THR_i < THR_min:
                    THR_i = THR_min
                    FUEL_i = FUEL_min
                    ROCD_i = (
                        conv.m2ft(
                            (1 / temp_const)
                            * (THR_i - Drag)
                            * TAS_i
                            * ESF_i
                            / (mass_i * const.g)
                        )
                        * 60
                    )
                elif THR_i > THR_max:
                    THR_i = THR_max
                    FUEL_i = FUEL_max
                    ROCD_i = (
                        conv.m2ft(
                            (1 / temp_const)
                            * (THR_i - Drag)
                            * TAS_i
                            * ESF_i
                            / (mass_i * const.g)
                        )
                        * 60
                    )
                else:
                    CT = AC.CT(Thrust=THR_i, delta=delta)
                    FUEL_i = AC.ff(
                        CT=CT, delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                    )  # [kg/s]
                    ROCD_i = conv.m2ft(ROCDisu) * 60

            # BADA3
            elif AC.BADAFamily.BADA3:
                THR_min = AC.Thrust(
                    rating="LIDL", v=TAS_i, h=H_m, config="CR", DeltaTemp=DeltaTemp
                )  # IDLE Thrust
                FUEL_min = AC.ff(
                    flightPhase="Descent",
                    v=TAS_i,
                    h=H_m,
                    T=THR_min,
                    config="CR",
                    adapted=False,
                )  # IDLE Fuel Flow
                THR_max = AC.Thrust(
                    rating="MCMB", v=TAS_i, h=H_m, DeltaTemp=DeltaTemp, config="CR"
                )  # MCMB Thrust
                FUEL_max = AC.ff(
                    flightPhase="Climb",
                    v=TAS_i,
                    h=H_m,
                    T=THR_max,
                    config="CR",
                    adapted=False,
                )  # MCMB Fuel Flow

                if THR_i < THR_min:
                    THR_i = THR_min
                    FUEL_i = FUEL_min
                    ROCD_i = (
                        conv.m2ft(
                            (1 / temp_const)
                            * (THR_i - Drag)
                            * TAS_i
                            * ESF_i
                            / (mass_i * const.g)
                        )
                        * 60
                    )
                elif THR_i > THR_max:
                    THR_i = THR_max
                    FUEL_i = FUEL_max
                    ROCD_i = (
                        conv.m2ft(
                            (1 / temp_const)
                            * (THR_i - Drag)
                            * TAS_i
                            * ESF_i
                            / (mass_i * const.g)
                        )
                        * 60
                    )
                else:
                    FUEL_i = AC.ff(
                        v=TAS_i, h=H_m, T=THR_i, config=config_i, adapted=True
                    )
                    ROCD_i = conv.m2ft(ROCDisu) * 60

            # Compute elapsed time and fuel burn over current step
            if Hp_i == Hp_init:
                break
            else:
                # Average ROCD over step is the mean of initial and final ones
                step_ROCD = (ROCD[-1] + ROCD_i) / 2  # [ft/min]
                # Step time is: altitude differential divided by average ROCD
                step_time = 60 * (Hp_i - Hp[-1]) / ROCD_i  # [s]

                # BADAE
                if AC.BADAFamily.BADAE:
                    # Average SOC rate over step is the mean of initial and final ones
                    step_SOCr = (SOCr[-1] + SOCr_i) / 2  # [%/h]
                    # SOC change is: average SOC rate multiplied by step time
                    step_SOC = step_SOCr * step_time / 3600  # [%]
                    # Update SOC estimate at end of step
                    SOC_i = SOC[-1] - step_SOC  # [%]
                    # update of aircraft mass estimate at end of step - mass is not changing for ELECTRIC engine (no fuel is consumed)
                    mass_i = mass[-1]  # [kg]

                else:
                    # average fuel flow over step is the mean of initial and final ones
                    step_FUEL = (FUEL[-1] + FUEL_i) / 2  # [kg/s]
                    # fuel burnt is: average fuel flow multiplied by step time
                    step_mass = step_FUEL * step_time  # [kg]
                    # update of aircraft mass estimate at end of step
                    if not mass_const:
                        mass_i = mass[-1] - step_mass  # [kg]
                        fuelConsumed_i = step_FUEL * step_time
                    fuelConsumed_i = FUELCONSUMED[-1] + step_FUEL * step_time

        ## PART 3: store information about end of step point
        writeOutputData = True
        if phase == "Climb" and ROCD_i < 0:
            warnings.warn(
                "Value ROCD = "
                + str(ROCD_i)
                + " [ft/min] is negative at the altitude "
                + str(Hp_i)
                + " [ft].",
                UserWarning,
            )
            go_on = False
            writeOutputData = False

        elif phase == "Climb" and ROCD_i < ROCD_min:
            warnings.warn(
                "Value ROCD = "
                + str(ROCD_i)
                + " [ft/min] exceeds the service ceiling limit defined by minimum ROCD = "
                + str(ROCD_min)
                + " [ft/min] at the altitude "
                + str(Hp_i)
                + " [ft].",
                UserWarning,
            )
            go_on = False
            writeOutputData = True

        if writeOutputData:
            # point data
            Hp.append(Hp_i)
            TAS.append(conv.ms2kt(TAS_i))
            CAS.append(conv.ms2kt(CAS_i))
            M.append(M_i)
            ROCD.append(ROCD_i)
            esf.append(ESF_i)
            Comment.append(comment)

            # everything except electric BADAE
            if not AC.BADAFamily.BADAE:
                FUEL.append(FUEL_i)

            # BADAH
            if AC.BADAFamily.BADAH:
                Peng.append(Preq_target_i)
                Preq.append(Preq_i)
                Pav.append(Pav_i)

            # BADAE
            elif AC.BADAFamily.BADAE:
                Pmec.append(Preq_target_i)
                Pbat.append(Pbat_i)
                SOCr.append(SOCr_i)
                Pelc.append(Pelc_i)
                Ibat.append(Ibat_i)
                Vbat.append(Vbat_i)
                Vgbat.append(Vgbat_i)

            # BADA3 & BADA4
            elif AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
                THR.append(THR_i)
                DRAG.append(Drag)
                config.append(config_i)

            # BADA4
            if AC.BADAFamily.BADA4:
                HLid.append(HLid_i)
                LG.append(LG_i)

            # calculation of the slope
            if TAS_i == 0:
                gamma_i = 90 * np.sign(ROCD_i)
            else:
                if AC.BADAFamily.BADAE:
                    gamma_i = degrees(atan(conv.ft2m(ROCD_i) * temp_const / 60 / TAS_i))
                else:
                    # using SIN assumes the TAS to be in the direction of the aircraft axis, not ground plane. Which means, this should be mathematically the correct equation for all the aircraft
                    gamma_i = degrees(asin(conv.ft2m(ROCD_i) * temp_const / 60 / TAS_i))

            # ground speed can be calcualted as TAS projected on the ground minus wind
            GS_i = cos(radians(gamma_i)) * TAS_i - wS
            GS.append(conv.ms2kt(GS_i))

            Slope.append(gamma_i)
            acc.append(0.0)

            BankAngle.append(bankAngle)
            ROT.append(rateOfTurn)

            # integrated data
            if Hp_i != Hp_init:
                if AC.BADAFamily.BADAE:
                    SOC.append(SOC_i)

                mass.append(mass_i)
                time.append(time[-1] + step_time)

                # everything except electric BADAE
                if not AC.BADAFamily.BADAE:
                    FUELCONSUMED.append(fuelConsumed_i)

                # Average TAS over step is the mean of initial and final ones
                step_TAS = (TAS[-2] + TAS[-1]) / 2  # [kt]
                # Average slope over the step
                step_gamma = radians((Slope[-2] + Slope[-1]) / 2)  # radians
                # Average ground speed over step
                # since this is not level flight, TAS speed should be projected on the ground, then GS can be calculated applying the wind speed
                step_TAS_projected = cos(step_gamma) * step_TAS
                step_GS = step_TAS_projected - wS  # [kt]
                # Step distance is: average GS multiplied by step time
                if turnFlight:
                    step_distance = conv.m2nm(
                        turn.distance(
                            rateOfTurn=rateOfTurn, TAS=TAS_i, timeOfTurn=step_time
                        )
                    )  # arcLength during the turn [NM]
                else:
                    step_distance = step_GS * step_time / 3600  # [NM]
                # Distance at end of step is distance at start of step plus step distance
                dist.append(dist[-1] + step_distance)

                # add GPS calculation
                if Lat and Lon and (magneticHeading or trueHeading):
                    if headingToFly == "TRUE":
                        if not turnFlight:
                            if not constantHeading:
                                # fly ORTHODROME
                                (Lat_i, Lon_i, HDGTrue_i) = (
                                    vincenty.destinationPoint_finalBearing(
                                        LAT_init=LAT[-1],
                                        LON_init=LON[-1],
                                        distance=conv.nm2m(step_distance),
                                        bearing=HDGTrue[-1],
                                    )
                                )

                                if magneticDeclinationGrid is not None:
                                    HDGMagnetic_i = (
                                        HDGTrue_i
                                        - magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGMagnetic_i = HDGTrue_i

                            elif constantHeading:
                                # fly LOXODROME
                                (Lat_i, Lon_i) = rhumb.destinationPoint(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearing=HDGTrue[-1],
                                    distance=conv.nm2m(step_distance),
                                )
                                HDGTrue_i = HDGTrue[-1]

                                if magneticDeclinationGrid is not None:
                                    HDGMagnetic_i = (
                                        HDGTrue_i
                                        - magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGMagnetic_i = HDGTrue_i

                        else:
                            # calculate the turn
                            (Lat_i, Lon_i, HDGTrue_i) = (
                                turn.destinationPoint_finalBearing(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearingInit=HDGTrue[-1],
                                    TAS=TAS_i,
                                    rateOfTurn=rateOfTurn,
                                    timeOfTurn=step_time,
                                    directionOfTurn=directionOfTurn,
                                )
                            )

                            if magneticDeclinationGrid is not None:
                                HDGMagnetic_i = (
                                    HDGTrue_i
                                    - magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGMagnetic_i = HDGTrue_i

                    elif headingToFly == "MAGNETIC":
                        if not turnFlight:
                            if constantHeading:
                                (Lat_i, Lon_i) = rhumb.destinationPoint(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearing=HDGTrue[-1],
                                    distance=conv.nm2m(step_distance),
                                )
                                HDGMagnetic_i = HDGMagnetic[-1]

                                if magneticDeclinationGrid is not None:
                                    HDGTrue_i = (
                                        HDGMagnetic_i
                                        + magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGTrue_i = HDGMagnetic_i

                        else:
                            # calculate the turn
                            (Lat_i, Lon_i, HDGTrue_i) = (
                                turn.destinationPoint_finalBearing(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearingInit=HDGTrue[-1],
                                    TAS=TAS_i,
                                    rateOfTurn=rateOfTurn,
                                    timeOfTurn=step_time,
                                    directionOfTurn=directionOfTurn,
                                )
                            )

                            if magneticDeclinationGrid is not None:
                                HDGMagnetic_i = (
                                    HDGTrue_i
                                    - magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGMagnetic_i = HDGTrue_i

                    LAT.append(Lat_i)
                    LON.append(Lon_i)
                    HDGMagnetic.append(HDGMagnetic_i)
                    HDGTrue.append(HDGTrue_i)

            # Determine end altitude of next step
            Hp_next = Hp_i + Hp_step

            if phase == "Climb":
                if Hp_next < Hp_final:
                    Hp_i = Hp_next - (Hp_i % Hp_step)
                # remaining altitude step would cross over the final altitude
                elif Hp_i < Hp_final:
                    Hp_i = Hp_final
                else:
                    go_on = False
            else:
                if Hp_next > Hp_final:
                    Hp_i = Hp_next - (Hp_i % Hp_step)
                # remaining altitude step would cross over the final altitude
                elif Hp_i > Hp_final:
                    Hp_i = Hp_final
                else:
                    go_on = False

    flightData = {
        "Hp": Hp,
        "TAS": TAS,
        "CAS": CAS,
        "GS": GS,
        "M": M,
        "acc": acc,
        "ROCD": ROCD,
        "ESF": esf,
        "FUEL": FUEL,
        "Pmec": Pmec,
        "Pelc": Pelc,
        "Pbat": Pbat,
        "SOCr": SOCr,
        "SOC": SOC,
        "Ibat": Ibat,
        "Vbat": Vbat,
        "Vgbat": Vgbat,
        "FUELCONSUMED": FUELCONSUMED,
        "Preq": Preq,
        "Peng": Peng,
        "Pav": Pav,
        "THR": THR,
        "DRAG": DRAG,
        "time": time,
        "dist": dist,
        "slope": Slope,
        "mass": mass,
        "config": config,
        "HLid": HLid,
        "LG": LG,
        "LAT": LAT,
        "LON": LON,
        "HDGTrue": HDGTrue,
        "HDGMagnetic": HDGMagnetic,
        "BankAngle": BankAngle,
        "ROT": ROT,
        "comment": Comment,
    }

    flightTrajectory = FT.createFlightTrajectoryDataframe(flightData)
    return flightTrajectory


def constantSpeedSlope_time(
    AC,
    length,
    speedType,
    v,
    Hp_init,
    slopetarget,
    m_init,
    DeltaTemp,
    wS=0.0,
    turnMetrics={"rateOfTurn": 0.0, "bankAngle": 0.0, "directionOfTurn": None},
    Lat=None,
    Lon=None,
    initialHeading={"magnetic": None, "true": None, "constantHeading": None},
    reducedPower=None,
    directionOfTurn=None,
    magneticDeclinationGrid=None,
    **kwargs,
):
    """This function computes time and fuel required by an aircraft to perform a climb/descent from Hp_init for set amount of time at constant speed and constant slope

    :param AC: aircraft {BADA3/4/H/E}
    :param speedType: what kind of speed is followed {M, CAS, TAS}.
    :param length: length of a segment to fly [s]
    :param step_length: length of a step of a segment - [s]
    :param v: what kind of speed is followed - [kt] CAS/TAS speed to follow or [-] MACH speed to follow.
    :param Hp_init: initial pressure altitude [ft].
    :param slopetarget: slope to be followed [deg].
    :param m_init: initial aircraft mass [kg].
    :param DeltaTemp: deviation with respect to ISA [K].
    :param wS: longitudinal wind speed (TAS) [kt].
    :param turnMetrics: Metrics for turn performance {"rateOfTurn":0.0,"bankAngle":0.0,"directionOfTurn":None} {[deg/s],[deg],[LEFT/RIGHT]}
    :param Hp_step: length of an altitude step of a segment [ft].
    :param SOC_init: initial state of charge [%].
    :param config: aircraft default aerodynamic configuration {TO,IC,CR,AP,LD}.
    :param speedBrakes: deployed or not speedbrakes including value to be added to the drag coeffcient {deployed:False,value:0.03} {deployed:[True/False],value:[-]}.
    :param ROCD_min: lower ROCD threshold to identify the climbing capabilities (service ceiling) [ft/min].
    :param Lat: Geographical Latitude [deg]
    :param Lon: Geographical Longitude [deg]
    :param initialHeading: aircraft magnetic heading, true heading and definition of constant heading(ORTHODROME=False, LOXODROME=True) {[deg],[deg],-}
    :param magneticDeclinationGrid: geographical grid of a magnetic declination on Earth [deg]
    :param mass_const: kind of mass canculation {mass_integrated=False, mass_constant=True}.
    :param m_iter: number of iterations for integration loop [-]
    :param reducedPower: reduction of Power during the climb {True/False}
    :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}.
    :type speedType: string.
    :type length: float.
    :type step_length: float.
    :type v: float.
    :type Hp_init: float.
    :type slopetarget: float.
    :type m_init: float.
    :type DeltaTemp: float.
    :type wS: float.
    :type turnMetrics: {float,float,string}.
    :type Hp_step: float.
    :type SOC_init: float.
    :type config: string.
    :type speedBrakes: dict{boolean,float}.
    :type ROCD_min: float.
    :type Lat: float.
    :type Lon: float.
    :type initialHeading: {float,float,boolean}.
    :type magneticDeclinationGrid: magneticDeclinationGrid.
    :type mass_const: boolean.
    :type m_iter: integer.
    :returns:
            BADA3: [Hp, TAS, CAS, M, GS, acc, ROCD, ESF, FUEL, FUELCONSUMED, THR, time, dist, slope, mass, config, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, N, s, NM, deg, kg, -, -,deg,deg,deg,deg,deg,deg/s]
            BADA4: [Hp, TAS, CAS, M, GS, acc, ROCD, ESF, FUEL, FUELCONSUMED, THR, P[Pmec, Pbat, Pelc Ibat, Vbat, Vgbat, SOCr, SOC], time, dist, slope, mass, config, HLid, LG, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, N, [W,W,W,A,V,V,%/h,%], s, NM, deg, kg, -, -, -, -,deg,deg,deg,deg,deg,deg/s]
            BADAH: [Hp, TAS, CAS, M, GS, ROCD, ESF, FUEL, FUELCONSUMED, Peng, Preq, Pav, time, dist, slope, mass, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, W, W, W, s, NM, deg, kg, -,deg,deg,deg,deg,deg,deg/s]
            BADAE: [time, dist, Hp, TAS, CAS, M, GS, acc, ROCD, ESF, slope, mass, P[Pmec, Pelc, Pbat, SOCr, SOC, Ibat, Vbat, Vgbat] comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [s, NM, ft, kt, kt, -, kt, m/s^2, ft/min, deg, kg, [W,W,W,%/h,%,A,V,V], -,deg,deg,deg,deg,deg,deg/s]
    :rtype: dict[list[float]}.
    """

    rateOfTurn = turnMetrics["rateOfTurn"]
    bankAngle = turnMetrics["bankAngle"]
    directionOfTurn = turnMetrics["directionOfTurn"]

    turnFlight = True
    if turnMetrics["rateOfTurn"] == 0.0 and turnMetrics["bankAngle"] == 0.0:
        turnFlight = False

    # conversion of Magnetic Heading to True Heading
    if magneticDeclinationGrid is not None:
        magneticDeclination = magneticDeclinationGrid.getMagneticDeclination(
            LAT_target=Lat, LON_target=Lon
        )
    else:
        magneticDeclination = 0

    # retrieve magnetic and true heading inputs
    magneticHeading = initialHeading["magnetic"]
    trueHeading = initialHeading["true"]
    constantHeading = initialHeading["constantHeading"]

    if Lat and Lon and (magneticHeading or trueHeading):
        if trueHeading is not None and magneticHeading is None:
            # fly TRUE Heading
            headingToFly = "TRUE"
            magneticHeading = trueHeading - magneticDeclination
        elif magneticHeading is not None and trueHeading is None:
            # fly MAGNETIC Heading
            if constantHeading == True:
                headingToFly = "MAGNETIC"
                trueHeading = magneticHeading + magneticDeclination
            else:
                raise Exception("Cannot fly non-constant magnetic heading")

        else:
            raise Exception("Undefined Heading value combination")

    # calculation with constant mass (True) or integrated (False)
    mass_const = kwargs.get("mass_const", False)

    # optional parameter to define initial Baterry State of Charge (SOC)
    if AC.BADAFamily.BADAE:
        SOC_init = kwargs.get("SOC_init", 100)
    else:
        SOC_init = None

    # speed brakes application
    if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
        speedBrakes = kwargs.get("speedBrakes", {"deployed": False, "value": 0.03})

    # step size in [s]
    step_length = kwargs.get("step_length", 1)

    # minimum remaining ROCD to determine cruise ceiling
    if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
        ROCD_min = kwargs.get("ROCD_min", 50)  # [ft/min]
    else:
        if AC.engineType == "PISTON" or AC.engineType == "ELECTRIC":
            ROCD_min = kwargs.get("ROCD_min", 100)  # [ft/min]
        else:
            ROCD_min = kwargs.get("ROCD_min", 300)  # [ft/min]

    # check the consistency of ROCD and climb/descent phase of flight
    if slopetarget < 0:
        phase = "Descent"
    elif slopetarget > 0:
        phase = "Climb"
    else:
        print("Slopetarget should be different from 0")

    if turnFlight:
        turnComment = "_turn"
    else:
        turnComment = ""

    if constantHeading:
        constHeadingStr = "_const_Heading"
    elif constantHeading is False or constantHeading is None:
        constHeadingStr = ""

    # comment line describing type of trajectory calculation
    comment = phase + turnComment + "_const_Slope_" + speedType + constHeadingStr

    if Lat and Lon and (magneticHeading or trueHeading):
        comment = comment + "_" + headingToFly + "_Heading"

    # get the default aerodynamic configuration if provided to be used for the whole segment
    config_default = kwargs.get("config", None)
    if config_default is not None:
        if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
            if not (
                config_default == "TO"
                or config_default == "IC"
                or config_default == "CR"
                or config_default == "AP"
                or config_default == "LD"
            ):
                print(
                    "WRONG default configuration set. Available values are: TO/IC/CR/AP/LD. Configuration will be calculated automatically"
                )

    #  weight iteration constant
    m_iter = kwargs.get("m_iter", 5)  # number of iterations for integration loop[-]

    # initialize output parameters
    Hp = [Hp_init]
    TAS = []
    CAS = []
    GS = []
    M = []
    ROCD = []
    esf = []
    FUEL = []
    FUELCONSUMED = []
    time = [0]
    dist = [0]
    mass = [m_init]
    Comment = []
    Slope = []
    acc = []
    THR = []
    DRAG = []
    config = []
    HLid = []
    LG = []
    BankAngle = []
    ROT = []

    if not AC.BADAFamily.BADAE:
        FUELCONSUMED = [0]

    # BADAH specific
    Peng = []
    Preq = []
    Pav = []

    # optional GPS coordiantes and HDG definition
    if Lat and Lon and (magneticHeading or trueHeading):
        LAT = [Lat]
        LON = [Lon]
        HDGMagnetic = [magneticHeading]
        HDGTrue = [trueHeading]
    else:
        LAT = []
        LON = []
        HDGMagnetic = []
        HDGTrue = []

    # BADAE specific
    Pmec = []
    Pbat = []
    SOCr = []
    SOC = [SOC_init]
    Pelc = []
    Ibat = []
    Vbat = []
    Vgbat = []

    # init loop parameters
    length_loop = 0
    time_i = time[-1]
    go_on = True

    while go_on:
        ## PART 1: compute parameters at end of step that are known without uncertainties:
        ##         atmosphere, speeds, ESF

        ## PART 2: compute parameters at end of step that are known only with uncertainties
        ##           (due to unknown mass at end of step):
        ##         weight, lift, drag , thrust, fuel flow

        mass_i = mass[-1]
        Hp_i = Hp[-1]
        for _ in itertools.repeat(None, m_iter):
            # atmosphere properties
            H_m = conv.ft2m(Hp_i)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
            temp_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)

            # aircraft speed
            [M_i, CAS_i, TAS_i] = atm.convertSpeed(
                v=v, speedType=speedType, theta=theta, delta=delta, sigma=sigma
            )

            if turnFlight:
                if turnMetrics["bankAngle"] != 0.0:
                    # bankAngle is defined
                    rateOfTurn = AC.rateOfTurn_bankAngle(TAS=TAS_i, bankAngle=bankAngle)
                else:
                    # rateOfTurn is defined
                    bankAngle = AC.bankAngle(
                        rateOfTurn=rateOfTurn, v=TAS_i
                    )  # [degrees]

            # Load factor
            nz = 1 / cos(radians(bankAngle))

            # compute Energy Share Factor (ESF)
            ESF_i = AC.esf(
                h=H_m, M=M_i, DeltaTemp=DeltaTemp, flightEvolution=("const" + speedType)
            )

            step_time = length_loop - time[-1]

            # Compute required ROCD
            if AC.BADAFamily.BADAE:
                # special case for BADAE, in future it may apply also for BADAH
                ROCDisu = tan(conv.deg2rad(slopetarget)) * TAS_i * (1 / temp_const)
            else:
                ROCDisu = sin(conv.deg2rad(slopetarget)) * TAS_i * (1 / temp_const)

            # BADAH or BADAE
            if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
                # compute Power required for level flight
                Preq_i = AC.Preq(sigma=sigma, tas=TAS_i, mass=mass_i, phi=bankAngle)
                # Compute power required for target ROCD
                Preq_target_i = AC.Peng_target(
                    temp=theta * const.temp_0,
                    DeltaTemp=DeltaTemp,
                    ROCD=ROCDisu,
                    mass=mass_i,
                    Preq=Preq_i,
                    ESF=ESF_i,
                )

            # BADA4
            elif AC.BADAFamily.BADA4:
                # aircraft aerodynamic configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=phase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=phase, previousConfig=config[-1], currentConfig=config_i
                    )

                [HLid_i, LG_i] = AC.flightEnvelope.getAeroConfig(config=config_i)

                # compute lift coefficient
                CL = AC.CL(M=M_i, delta=delta, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(M=M_i, CL=CL, HLid=HLid_i, LG=LG_i, speedBrakes=speedBrakes)
                # compute drag force
                Drag = AC.D(M=M_i, delta=delta, CD=CD)
                # compute thrust force
                THR_i = (
                    ROCDisu * mass_i * const.g * temp_const / (TAS_i * ESF_i) + Drag
                )  # [N]

            # BADA3
            elif AC.BADAFamily.BADA3:
                # aircraft aerodynamic configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=phase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=phase,
                        previousConfig=config[-1],
                        currentConfig=config_i,
                    )

                # compute lift coefficient
                CL = AC.CL(tas=TAS_i, sigma=sigma, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(CL=CL, config=config_i, speedBrakes=speedBrakes)
                # compute drag force
                Drag = AC.D(tas=TAS_i, sigma=sigma, CD=CD)
                # compute thrust force
                THR_i = (
                    ROCDisu * mass_i * const.g * temp_const / (TAS_i * ESF_i) + Drag
                )  # [N]

            # check that required thrust/power fits in the avialable thrust/power envelope,
            # recompute ROCD if necessary and compute fuel flow

            # BADAH
            if AC.BADAFamily.BADAH:
                Pmin = 0.1 * AC.P0  # No minimum power model: assume 10% torque
                Pav_i = AC.Pav(
                    rating="MTKF", theta=theta, delta=delta
                )  # assume MTKF rating as the limit
                Pmax = Pav_i

                if Preq_target_i < Pmin:
                    Preq_target_i = Pmin
                    ROCD_i = (
                        conv.m2ft(
                            AC.ROCD(
                                Peng=Preq_target_i,
                                Preq=Preq_i,
                                mass=mass_i,
                                ESF=ESF_i,
                                theta=theta,
                                DeltaTemp=DeltaTemp,
                            )
                        )
                        * 60
                    )
                elif Preq_target_i > Pmax:
                    Preq_target_i = Pmax
                    ROCD_i = (
                        conv.m2ft(
                            AC.ROCD(
                                Peng=Preq_target_i,
                                Preq=Preq_i,
                                mass=mass_i,
                                ESF=ESF_i,
                                theta=theta,
                                DeltaTemp=DeltaTemp,
                            )
                        )
                        * 60
                    )
                else:
                    ROCD_i = conv.m2ft(ROCDisu) * 60

                # compute fuel flow for level flight
                CP = AC.CP(Peng=Preq_target_i)
                FUEL_i = AC.ff(delta=delta, CP=CP)  # [kg/s]

            # BADAE
            elif AC.BADAFamily.BADAE:
                Pmin = 0.1 * AC.P0  # No minimum power model: assume 10% torque
                Pav_i = AC.Pav(
                    rating="MTKF", SOC=SOC[-1]
                )  # assume MTKF rating as the limit
                Pmax = Pav_i

                if Preq_target_i < Pmin:
                    Preq_target_i = Pmin
                    ROCD_i = (
                        conv.m2ft(
                            AC.ROCD(
                                Peng=Preq_target_i,
                                Preq=Preq_i,
                                mass=mass_i,
                                ESF=ESF_i,
                                theta=theta,
                                DeltaTemp=DeltaTemp,
                            )
                        )
                        * 60
                    )
                elif Preq_target_i > Pmax:
                    Preq_target_i = Pmax
                    ROCD_i = (
                        conv.m2ft(
                            AC.ROCD(
                                Peng=Preq_target_i,
                                Preq=Preq_i,
                                mass=mass_i,
                                ESF=ESF_i,
                                theta=theta,
                                DeltaTemp=DeltaTemp,
                            )
                        )
                        * 60
                    )
                else:
                    ROCD_i = conv.m2ft(ROCDisu) * 60

                Pbat_i = AC.Pbat(Preq=Preq_target_i, SOC=SOC[-1])
                SOCr_i = AC.SOCrate(Preq=Preq_target_i, SOC=SOC[-1])

                # debug data
                Pelc_i = Preq_target_i / AC.eta
                Ibat_i = AC.Ibat(P=Pelc_i, SOC=SOC[-1])
                Vbat_i = AC.Vbat(I=Ibat_i, SOC=SOC[-1])
                Vgbat_i = AC.Vocbat(SOC=SOC[-1]) - AC.R0bat(SOC=SOC[-1]) * Ibat_i

            # BADA4
            elif AC.BADAFamily.BADA4:
                THR_min = AC.Thrust(
                    rating="LIDL", delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )  # IDLE Thrust
                FUEL_min = AC.ff(
                    rating="LIDL", delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )  # IDLE Fuel Flow
                THR_max = AC.Thrust(
                    rating="MCMB", delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )  # MCMB Thrust
                FUEL_max = AC.ff(
                    rating="MCMB", delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )  # MCMB Fuel Flow

                if THR_i < THR_min:
                    THR_i = THR_min
                    FUEL_i = FUEL_min
                    ROCD_i = (
                        conv.m2ft(
                            (1 / temp_const)
                            * (THR_i - Drag)
                            * TAS_i
                            * ESF_i
                            / (mass_i * const.g)
                        )
                        * 60
                    )
                elif THR_i > THR_max:
                    THR_i = THR_max
                    FUEL_i = FUEL_max
                    ROCD_i = (
                        conv.m2ft(
                            (1 / temp_const)
                            * (THR_i - Drag)
                            * TAS_i
                            * ESF_i
                            / (mass_i * const.g)
                        )
                        * 60
                    )
                else:
                    CT = AC.CT(Thrust=THR_i, delta=delta)
                    FUEL_i = AC.ff(
                        CT=CT, delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                    )  # [kg/s]
                    ROCD_i = conv.m2ft(ROCDisu) * 60

            # BADA3
            elif AC.BADAFamily.BADA3:
                THR_min = AC.Thrust(
                    rating="LIDL", v=TAS_i, h=H_m, config="CR", DeltaTemp=DeltaTemp
                )  # IDLE Thrust
                FUEL_min = AC.ff(
                    flightPhase="Descent",
                    v=TAS_i,
                    h=H_m,
                    T=THR_min,
                    config="CR",
                    adapted=False,
                )  # IDLE Fuel Flow
                THR_max = AC.Thrust(
                    rating="MCMB", v=TAS_i, h=H_m, DeltaTemp=DeltaTemp, config="CR"
                )  # MCMB Thrust
                FUEL_max = AC.ff(
                    flightPhase="Climb",
                    v=TAS_i,
                    h=H_m,
                    T=THR_max,
                    config="CR",
                    adapted=False,
                )  # MCMB Fuel Flow

                if THR_i < THR_min:
                    THR_i = THR_min
                    FUEL_i = FUEL_min
                    ROCD_i = (
                        conv.m2ft(
                            (1 / temp_const)
                            * (THR_i - Drag)
                            * TAS_i
                            * ESF_i
                            / (mass_i * const.g)
                        )
                        * 60
                    )
                elif THR_i > THR_max:
                    THR_i = THR_max
                    FUEL_i = FUEL_max
                    ROCD_i = (
                        conv.m2ft(
                            (1 / temp_const)
                            * (THR_i - Drag)
                            * TAS_i
                            * ESF_i
                            / (mass_i * const.g)
                        )
                        * 60
                    )
                else:
                    FUEL_i = AC.ff(
                        v=TAS_i, h=H_m, T=THR_i, config=config_i, adapted=True
                    )
                    ROCD_i = conv.m2ft(ROCDisu) * 60

            # Compute elapsed time and fuel burn over current step
            if length_loop == 0:
                # no need to loop for first point: initial m/Hp already known
                break
            else:
                # Average ROCD over step is the mean of initial and final ones
                step_ROCD = (ROCD[-1] + ROCD_i) / 2  # [ft/min]
                # Altitude differential is: average ROCD multiplied by step time
                step_Hp = step_ROCD * step_time / 60  # [ft]
                # Update altitude estimate at end of step
                Hp_i = Hp[-1] + step_Hp  # [ft]

                # BADAE
                if AC.BADAFamily.BADAE:
                    # Average SOC rate over step is the mean of initial and final ones
                    step_SOCr = (SOCr[-1] + SOCr_i) / 2  # [%/h]
                    # SOC change is: average SOC rate multiplied by step time
                    step_SOC = step_SOCr * step_time / 3600  # [%]
                    # Update SOC estimate at end of step
                    SOC_i = SOC[-1] - step_SOC  # [%]
                    # update of aircraft mass estimate at end of step - mass is not changing for ELECTRIC engine (no fuel is consumed)
                    mass_i = mass[-1]  # [kg]

                else:
                    # average fuel flow over step is the mean of initial and final ones
                    step_FUEL = (FUEL[-1] + FUEL_i) / 2  # [kg/s]
                    # fuel burnt is: average fuel flow multiplied by step time
                    step_mass = step_FUEL * step_time  # [kg]
                    # update of aircraft mass estimate at end of step
                    if not mass_const:
                        mass_i = mass[-1] - step_mass  # [kg]
                        fuelConsumed_i = step_FUEL * step_time
                    fuelConsumed_i = FUELCONSUMED[-1] + step_FUEL * step_time

        ## PART 3: store information about end of step point
        writeOutputData = True
        if phase == "Climb" and ROCD_i < 0:
            warnings.warn(
                "Value ROCD = "
                + str(ROCD_i)
                + " [ft/min] is negative at the altitude "
                + str(Hp_i)
                + " [ft].",
                UserWarning,
            )
            go_on = False
            writeOutputData = False

        elif phase == "Climb" and ROCD_i < ROCD_min:
            warnings.warn(
                "Value ROCD = "
                + str(ROCD_i)
                + " [ft/min] exceeds the service ceiling limit defined by minimum ROCD = "
                + str(ROCD_min)
                + " [ft/min] at the altitude "
                + str(Hp_i)
                + " [ft].",
                UserWarning,
            )
            go_on = False
            writeOutputData = True

        if writeOutputData:
            # point data
            TAS.append(conv.ms2kt(TAS_i))
            CAS.append(conv.ms2kt(CAS_i))
            M.append(M_i)
            ROCD.append(ROCD_i)
            esf.append(ESF_i)
            Comment.append(comment)

            # everything except electric BADAE
            if not AC.BADAFamily.BADAE:
                FUEL.append(FUEL_i)

            # BADAH
            if AC.BADAFamily.BADAH:
                Peng.append(Preq_target_i)
                Preq.append(Preq_i)
                Pav.append(Pav_i)

            # BADAE
            elif AC.BADAFamily.BADAE:
                Pmec.append(Preq_target_i)
                Pbat.append(Pbat_i)
                SOCr.append(SOCr_i)
                Pelc.append(Pelc_i)
                Ibat.append(Ibat_i)
                Vbat.append(Vbat_i)
                Vgbat.append(Vgbat_i)

            # BADA3 & BADA4
            elif AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
                THR.append(THR_i)
                DRAG.append(Drag)
                config.append(config_i)

            # BADA4
            if AC.BADAFamily.BADA4:
                HLid.append(HLid_i)
                LG.append(LG_i)

            # calculation of the slope
            if TAS_i == 0:
                gamma_i = 90 * np.sign(ROCD_i)
            else:
                [theta, delta, sigma] = atm.atmosphereProperties(
                    h=conv.ft2m(Hp_i), DeltaTemp=DeltaTemp
                )
                temp_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)
                if AC.BADAFamily.BADAE:
                    gamma_i = degrees(atan(conv.ft2m(ROCD_i) * temp_const / 60 / TAS_i))
                else:
                    # using SIN assumes the TAS to be in the direction of the aircraft axis, not ground plane. Which means, this should be mathematically the correct equation for all the aircraft
                    gamma_i = degrees(asin(conv.ft2m(ROCD_i) * temp_const / 60 / TAS_i))

            # ground speed can be calcualted as TAS projected on the ground minus wind
            GS_i = cos(radians(gamma_i)) * TAS_i - wS
            GS.append(conv.ms2kt(GS_i))

            Slope.append(gamma_i)
            acc.append(0.0)

            BankAngle.append(bankAngle)
            ROT.append(rateOfTurn)

            # integrated data
            if length_loop != 0:
                if AC.BADAFamily.BADAE:
                    SOC.append(SOC_i)

                Hp.append(Hp_i)
                mass.append(mass_i)
                time.append(time[-1] + step_time)

                # everything except electric BADAE
                if not AC.BADAFamily.BADAE:
                    FUELCONSUMED.append(fuelConsumed_i)

                # Average TAS over step is the mean of initial and final ones
                step_TAS = (TAS[-2] + TAS[-1]) / 2  # [kt]
                # Average slope over the step
                step_gamma = radians((Slope[-2] + Slope[-1]) / 2)  # radians
                # Average ground speed over step
                # since this is not level flight, TAS speed should be projected on the ground, then GS can be calculated applying the wind speed
                step_TAS_projected = cos(step_gamma) * step_TAS
                step_GS = step_TAS_projected - wS  # [kt]
                # Step distance is: average GS multiplied by step time
                if turnFlight:
                    step_distance = conv.m2nm(
                        turn.distance(
                            rateOfTurn=rateOfTurn, TAS=TAS_i, timeOfTurn=step_time
                        )
                    )  # arcLength during the turn [NM]
                else:
                    step_distance = step_GS * step_time / 3600  # [NM]
                # Distance at end of step is distance at start of step plus step distance
                dist.append(dist[-1] + step_distance)

                # add GPS calculation
                if Lat and Lon and (magneticHeading or trueHeading):
                    if headingToFly == "TRUE":
                        if not turnFlight:
                            if not constantHeading:
                                # fly ORTHODROME
                                (Lat_i, Lon_i, HDGTrue_i) = (
                                    vincenty.destinationPoint_finalBearing(
                                        LAT_init=LAT[-1],
                                        LON_init=LON[-1],
                                        distance=conv.nm2m(step_distance),
                                        bearing=HDGTrue[-1],
                                    )
                                )

                                if magneticDeclinationGrid is not None:
                                    HDGMagnetic_i = (
                                        HDGTrue_i
                                        - magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGMagnetic_i = HDGTrue_i

                            elif constantHeading:
                                # fly LOXODROME
                                (Lat_i, Lon_i) = rhumb.destinationPoint(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearing=HDGTrue[-1],
                                    distance=conv.nm2m(step_distance),
                                )
                                HDGTrue_i = HDGTrue[-1]

                                if magneticDeclinationGrid is not None:
                                    HDGMagnetic_i = (
                                        HDGTrue_i
                                        - magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGMagnetic_i = HDGTrue_i

                        else:
                            # calculate the turn
                            (Lat_i, Lon_i, HDGTrue_i) = (
                                turn.destinationPoint_finalBearing(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearingInit=HDGTrue[-1],
                                    TAS=TAS_i,
                                    rateOfTurn=rateOfTurn,
                                    timeOfTurn=step_time,
                                    directionOfTurn=directionOfTurn,
                                )
                            )

                            if magneticDeclinationGrid is not None:
                                HDGMagnetic_i = (
                                    HDGTrue_i
                                    - magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGMagnetic_i = HDGTrue_i

                    elif headingToFly == "MAGNETIC":
                        if not turnFlight:
                            if constantHeading:
                                (Lat_i, Lon_i) = rhumb.destinationPoint(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearing=HDGTrue[-1],
                                    distance=conv.nm2m(step_distance),
                                )
                                HDGMagnetic_i = HDGMagnetic[-1]

                                if magneticDeclinationGrid is not None:
                                    HDGTrue_i = (
                                        HDGMagnetic_i
                                        + magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGTrue_i = HDGMagnetic_i

                        else:
                            # calculate the turn
                            (Lat_i, Lon_i, HDGTrue_i) = (
                                turn.destinationPoint_finalBearing(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearingInit=HDGTrue[-1],
                                    TAS=TAS_i,
                                    rateOfTurn=rateOfTurn,
                                    timeOfTurn=step_time,
                                    directionOfTurn=directionOfTurn,
                                )
                            )

                            if magneticDeclinationGrid is not None:
                                HDGMagnetic_i = (
                                    HDGTrue_i
                                    - magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGMagnetic_i = HDGTrue_i

                    LAT.append(Lat_i)
                    LON.append(Lon_i)
                    HDGMagnetic.append(HDGMagnetic_i)
                    HDGTrue.append(HDGTrue_i)

            if length_loop + step_length < length:
                length_loop += step_length
            elif length_loop < length:
                length_loop = length
            else:
                go_on = False

    flightData = {
        "Hp": Hp,
        "TAS": TAS,
        "CAS": CAS,
        "GS": GS,
        "M": M,
        "acc": acc,
        "ROCD": ROCD,
        "ESF": esf,
        "FUEL": FUEL,
        "Pmec": Pmec,
        "Pelc": Pelc,
        "Pbat": Pbat,
        "SOCr": SOCr,
        "SOC": SOC,
        "Ibat": Ibat,
        "Vbat": Vbat,
        "Vgbat": Vgbat,
        "FUELCONSUMED": FUELCONSUMED,
        "Preq": Preq,
        "Peng": Peng,
        "Pav": Pav,
        "THR": THR,
        "DRAG": DRAG,
        "time": time,
        "dist": dist,
        "slope": Slope,
        "mass": mass,
        "config": config,
        "HLid": HLid,
        "LG": LG,
        "LAT": LAT,
        "LON": LON,
        "HDGTrue": HDGTrue,
        "HDGMagnetic": HDGMagnetic,
        "BankAngle": BankAngle,
        "ROT": ROT,
        "comment": Comment,
    }

    flightTrajectory = FT.createFlightTrajectoryDataframe(flightData)
    return flightTrajectory


def constantSpeedRating(
    AC,
    speedType,
    v,
    Hp_init,
    Hp_final,
    m_init,
    DeltaTemp,
    wS=0.0,
    turnMetrics={"rateOfTurn": 0.0, "bankAngle": 0.0, "directionOfTurn": None},
    Lat=None,
    Lon=None,
    initialHeading={"magnetic": None, "true": None, "constantHeading": None},
    reducedPower=None,
    directionOfTurn=None,
    expedite=False,
    magneticDeclinationGrid=None,
    initRating=None,
    **kwargs,
):
    """This function computes time and fuel required by an aircraft to perform a climb/descent from Hp_init to Hp_final at constant speed and constant engine rating

    :param AC: aircraft {BADA3/4/H/E}
    :param speedType: what kind of speed is followed {M, CAS, TAS}.
    :param v: what kind of speed is followed - [kt] CAS/TAS speed to follow or [-] MACH speed to follow.
    :param Hp_init: initial pressure altitude [ft].
    :param Hp_final: final pressure altitude [ft].
    :param m_init: initial aircraft mass [kg].
    :param DeltaTemp: deviation with respect to ISA [K].
    :param wS: longitudinal wind speed (TAS) [kt].
    :param turnMetrics: Metrics for turn performance {"rateOfTurn":0.0,"bankAngle":0.0,"directionOfTurn":None} {[deg/s],[deg],[LEFT/RIGHT]}
    :param Hp_step: length of an altitude step of a segment [ft].
    :param SOC_init: initial state of charge [%].
    :param config: aircraft default aerodynamic configuration {TO,IC,CR,AP,LD}.
    :param speedBrakes: deployed or not speedbrakes including value to be added to the drag coeffcient {deployed:False,value:0.03} {deployed:[True/False],value:[-]}.
    :param ROCD_min: lower ROCD threshold to identify the climbing capabilities (service ceiling) [ft/min].
    :param Lat: Geographical Latitude [deg]
    :param Lon: Geographical Longitude [deg]
    :param initialHeading: aircraft magnetic heading, true heading and definition of constant heading(ORTHODROME=False, LOXODROME=True) {[deg],[deg],-}
    :param magneticDeclinationGrid: geographical grid of a magnetic declination on Earth [deg]
    :param mass_const: kind of mass canculation {mass_integrated=False, mass_constant=True}.
    :param m_iter: number of iterations for integration loop [-]
    :param reducedPower: reduction of Power during the climb {True/False}
    :param initRating: default rating settings
    :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}.
    :type speedType: string.
    :type v: float.
    :type Hp_init: float.
    :type Hp_final: float.
    :type m_init: float.
    :type DeltaTemp: float.
    :type wS: float.
    :type turnMetrics: {float,float,string}.
    :type Hp_step: float.
    :type SOC_init: float.
    :type config: string.
    :type speedBrakes: dict{boolean,float}.
    :type ROCD_min: float.
    :type Lat: float.
    :type Lon: float.
    :type initialHeading: {float,float,boolean}.
    :type magneticDeclinationGrid: magneticDeclinationGrid.
    :type mass_const: boolean.
    :type m_iter: integer.
    :type reducedPower: boolean.
    :type initRating: string.
    :returns:
            BADA3: [Hp, TAS, CAS, M, GS, acc, ROCD, ESF, FUEL, FUELCONSUMED, THR, time, dist, slope, mass, config, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, N, s, NM, deg, kg, -, -,deg,deg,deg,deg,deg,deg/s]
            BADA4: [Hp, TAS, CAS, M, GS, acc, ROCD, ESF, FUEL, FUELCONSUMED, THR, P[Pmec, Pbat, Pelc Ibat, Vbat, Vgbat, SOCr, SOC], time, dist, slope, mass, config, HLid, LG, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, N, [W,W,W,A,V,V,%/h,%], s, NM, deg, kg, -, -, -, -,deg,deg,deg,deg,deg,deg/s]
            BADAH: [Hp, TAS, CAS, M, GS, ROCD, ESF, FUEL, FUELCONSUMED, Peng, Preq, Pav, time, dist, slope, mass, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, W, W, W, s, NM, deg, kg, -,deg,deg,deg,deg,deg,deg/s]
            BADAE: [time, dist, Hp, TAS, CAS, M, GS, acc, ROCD, ESF, slope, mass, P[Pmec, Pelc, Pbat, SOCr, SOC, Ibat, Vbat, Vgbat] comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [s, NM, ft, kt, kt, -, kt, m/s^2, ft/min, deg, kg, [W,W,W,%/h,%,A,V,V], -,deg,deg,deg,deg,deg,deg/s]
    :rtype: dict[list[float]}.
    """

    rateOfTurn = turnMetrics["rateOfTurn"]
    bankAngle = turnMetrics["bankAngle"]
    directionOfTurn = turnMetrics["directionOfTurn"]

    turnFlight = True
    if turnMetrics["rateOfTurn"] == 0.0 and turnMetrics["bankAngle"] == 0.0:
        turnFlight = False

    # conversion of Magnetic Heading to True Heading
    if magneticDeclinationGrid is not None:
        magneticDeclination = magneticDeclinationGrid.getMagneticDeclination(
            LAT_target=Lat, LON_target=Lon
        )
    else:
        magneticDeclination = 0

    # retrieve magnetic and true heading inputs
    magneticHeading = initialHeading["magnetic"]
    trueHeading = initialHeading["true"]
    constantHeading = initialHeading["constantHeading"]

    if Lat and Lon and (magneticHeading or trueHeading):
        if trueHeading is not None and magneticHeading is None:
            # fly TRUE Heading
            headingToFly = "TRUE"
            magneticHeading = trueHeading - magneticDeclination
        elif magneticHeading is not None and trueHeading is None:
            # fly MAGNETIC Heading
            if constantHeading == True:
                headingToFly = "MAGNETIC"
                trueHeading = magneticHeading + magneticDeclination
            else:
                raise Exception("Cannot fly non-constant magnetic heading")

        else:
            raise Exception("Undefined Heading value combination")

    # calculation with constant mass (True) or integrated (False)
    mass_const = kwargs.get("mass_const", False)

    # optional parameter to define initial Baterry State of Charge (SOC)
    if AC.BADAFamily.BADAE:
        SOC_init = kwargs.get("SOC_init", 100)
    else:
        SOC_init = None

    # speed brakes application
    if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
        speedBrakes = kwargs.get("speedBrakes", {"deployed": False, "value": 0.03})

    # optional parameter - iteration step for altitude loop
    if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
        Hp_step = kwargs.get("Hp_step", 500)  # [ft]
    else:
        # NB: it must be a multiple of 1000ft so that interrupted climbs end on a regular cruise altitude.
        Hp_step = kwargs.get("Hp_step", 1000)  # [ft]

    # minimum remaining ROCD to determine cruise ceiling
    if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
        ROCD_min = kwargs.get("ROCD_min", 50)  # [ft/min]
    else:
        if AC.engineType == "PISTON" or AC.engineType == "ELECTRIC":
            ROCD_min = kwargs.get("ROCD_min", 100)  # [ft/min]
        else:
            ROCD_min = kwargs.get("ROCD_min", 300)  # [ft/min]

    # determine if  vertical evolution over the segment is CLIMB or DESCENT
    # and associate engine rating  and altitude iteration direction

    if Hp_init < Hp_final:
        phase = "Climb"
    else:
        phase = "Descent"
        Hp_step = -Hp_step

    if initRating is None:
        if phase == "Climb":
            if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
                if v == 0:
                    rating = "MTKF"
                else:
                    rating = "MCNT"
            else:
                rating = "MCMB"
        elif phase == "Descent":
            if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
                if v == 0:
                    rating = "UNKNOWN"
                else:
                    rating = "UNKNOWN"
            else:
                rating = "LIDL"
    else:
        rating = initRating

    if turnFlight:
        turnComment = "_turn"
    else:
        turnComment = ""

    if constantHeading:
        constHeadingStr = "_const_Heading"
    elif constantHeading is False or constantHeading is None:
        constHeadingStr = ""

    # comment line describing type of trajectory calculation
    comment = (
        phase + turnComment + "_const_" + speedType + "_" + rating + constHeadingStr
    )

    if Lat and Lon and (magneticHeading or trueHeading):
        comment = comment + "_" + headingToFly + "_Heading"

    if expedite:
        comment = comment + "_expedite"

    # get the default aerodynamic configuration if provided to be used for the whole segment
    config_default = kwargs.get("config", None)
    if config_default is not None:
        if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
            if not (
                config_default == "TO"
                or config_default == "IC"
                or config_default == "CR"
                or config_default == "AP"
                or config_default == "LD"
            ):
                print(
                    "WRONG default configuration set. Available values are: TO/IC/CR/AP/LD. Configuration will be calculated automatically"
                )

    #  weight iteration constant
    m_iter = kwargs.get("m_iter", 5)  # number of iterations for integration loop[-]

    # The thrust_fuel method for BADA 3 models applies the cruise fuel correction
    # whenever the thrust is adapted, instead of only in cruise: this correction
    # needs to be reverted when thrust is adapted for constant ROC/slope.

    # cruise_correction = 1/f(5)

    # initialize output parameters
    Hp = []
    TAS = []
    CAS = []
    GS = []
    M = []
    ROCD = []
    esf = []
    FUEL = []
    FUELCONSUMED = []
    time = [0]
    dist = [0]
    mass = [m_init]
    Comment = []
    Slope = []
    acc = []
    THR = []
    DRAG = []
    config = []
    HLid = []
    LG = []
    BankAngle = []
    ROT = []

    if not AC.BADAFamily.BADAE:
        FUELCONSUMED = [0]

    # BADAH specific
    Peng = []
    Preq = []
    Pav = []

    # optional GPS coordiantes and HDG definition
    if Lat and Lon and (magneticHeading or trueHeading):
        LAT = [Lat]
        LON = [Lon]
        HDGMagnetic = [magneticHeading]
        HDGTrue = [trueHeading]
    else:
        LAT = []
        LON = []
        HDGMagnetic = []
        HDGTrue = []

    # BADAE specific
    Pmec = []
    Pbat = []
    SOCr = []
    SOC = [SOC_init]
    Pelc = []
    Ibat = []
    Vbat = []
    Vgbat = []

    # init loop parameters
    Hp_i = Hp_init
    go_on = True

    while go_on:
        ## PART 1: compute parameters at end of step that are known without uncertainties:
        ##         atmosphere, speeds, thrust. fuel flow, ESF
        # atmosphere properties
        H_m = conv.ft2m(Hp_i)  # altitude [m]
        [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
        temp_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)

        # aircraft speed
        [M_i, CAS_i, TAS_i] = atm.convertSpeed(
            v=v, speedType=speedType, theta=theta, delta=delta, sigma=sigma
        )

        if turnFlight:
            if turnMetrics["bankAngle"] != 0.0:
                # bankAngle is defined
                rateOfTurn = AC.rateOfTurn_bankAngle(TAS=TAS_i, bankAngle=bankAngle)
            else:
                # rateOfTurn is defined
                bankAngle = AC.bankAngle(rateOfTurn=rateOfTurn, v=TAS_i)  # [degrees]

        # Load factor
        nz = 1 / cos(radians(bankAngle))

        # compute Energy Share Factor (ESF)
        ESF_i = AC.esf(
            h=H_m, M=M_i, DeltaTemp=DeltaTemp, flightEvolution=("const" + speedType)
        )

        mass_i = mass[-1]

        # BADAH
        if AC.BADAFamily.BADAH:
            # compute available power
            if rating == "UNKNOWN":
                Preq_target_i = 0.1 * AC.P0  # No minimum power model: assume 10% torque
            else:
                Preq_target_i = AC.Pav(rating=rating, theta=theta, delta=delta)

            Pav_i = AC.Pav(rating="MTKF", theta=theta, delta=delta)

            # compute fuel flow for level flight
            CP = AC.CP(Peng=Preq_target_i)
            FUEL_i = AC.ff(delta=delta, CP=CP)  # [kg/s]

        # BADAE
        elif AC.BADAFamily.BADAE:
            # compute available power
            if rating == "UNKNOWN":
                Preq_target_i = 0.1 * AC.P0  # No minimum power model: assume 10% torque
            else:
                Preq_target_i = AC.Pav(rating=rating, SOC=SOC[-1])

                Pav_i = AC.Pav(rating=rating, SOC=SOC[-1])

            Pbat_i = AC.Pbat(Preq=Preq_target_i, SOC=SOC[-1])
            SOCr_i = AC.SOCrate(Preq=Preq_target_i, SOC=SOC[-1])

            # debug data
            Pelc_i = Preq_target_i / AC.eta
            Ibat_i = AC.Ibat(P=Pelc_i, SOC=SOC[-1])
            Vbat_i = AC.Vbat(I=Ibat_i, SOC=SOC[-1])
            Vgbat_i = AC.Vocbat(SOC=SOC[-1]) - AC.R0bat(SOC=SOC[-1]) * Ibat_i

        # BADA4
        elif AC.BADAFamily.BADA4:
            # compute thrust force and fuel flow
            THR_i = AC.Thrust(
                rating=rating, delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
            )  # [N]
            CT = AC.CT(Thrust=THR_i, delta=delta)
            FUEL_i = AC.ff(CT=CT, delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp)

        # BADA3
        elif AC.BADAFamily.BADA3:
            # aircraft aerodynamic configuration
            if config_default is None:
                config_i = AC.flightEnvelope.getConfig(
                    h=H_m,
                    phase=phase,
                    v=CAS_i,
                    mass=mass_i,
                    DeltaTemp=DeltaTemp,
                )
            else:
                config_i = config_default

            # ensure continuity of configuration change within the segment
            if config:
                config_i = AC.flightEnvelope.checkConfigurationContinuity(
                    phase=phase, previousConfig=config[-1], currentConfig=config_i
                )

            # compute thrust force and fuel flow
            THR_i = AC.Thrust(
                rating=rating, v=TAS_i, h=H_m, config=config_i, DeltaTemp=DeltaTemp
            )
            FUEL_i = AC.ff(
                flightPhase=phase,
                v=TAS_i,
                h=H_m,
                T=THR_i,
                config=config_i,
                adapted=True,
            )

        if Hp_i != Hp_init:  # exclude first point: initial m already known
            # BADAE
            if AC.BADAFamily.BADAE:
                # Average SOC rate over step is the mean of initial and final ones
                step_SOCr = (SOCr[-1] + SOCr_i) / 2  # [%/h]
            else:
                # average fuel flow over step is the mean of initial and final ones
                step_FUEL = (FUEL[-1] + FUEL_i) / 2  # [kg/s]

        ## PART 2: compute parameters at end of step that are known only with uncertainties
        ##           (due to unknown mass at end of step):
        ##         weight, lift, drag , ROCD

        for _ in itertools.repeat(None, m_iter):
            # BADAH or BADAE
            if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
                # compute Power required
                Preq_i = AC.Preq(sigma=sigma, tas=TAS_i, mass=mass_i, phi=bankAngle)
                # compute ROCD
                ROCD_i = (
                    conv.m2ft(
                        AC.ROCD(
                            Peng=Preq_target_i,
                            Preq=Preq_i,
                            mass=mass_i,
                            ESF=ESF_i,
                            theta=theta,
                            DeltaTemp=DeltaTemp,
                        )
                    )
                    * 60
                )

            # BADA4
            elif AC.BADAFamily.BADA4:
                # aircraft aerodynamic configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=phase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=phase, previousConfig=config[-1], currentConfig=config_i
                    )

                [HLid_i, LG_i] = AC.flightEnvelope.getAeroConfig(config=config_i)

                # compute lift coefficient
                CL = AC.CL(M=M_i, delta=delta, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(M=M_i, CL=CL, HLid=HLid_i, LG=LG_i, speedBrakes=speedBrakes)
                # compute drag force
                Drag = AC.D(M=M_i, delta=delta, CD=CD)
                # compute ROCD
                ROCD_i = (
                    conv.m2ft(
                        (1 / temp_const)
                        * (THR_i - Drag)
                        * TAS_i
                        * ESF_i
                        / (mass_i * const.g)
                    )
                    * 60
                )

            # BADA3
            elif AC.BADAFamily.BADA3:
                # aircraft aerodynamic configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=phase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=phase,
                        previousConfig=config[-1],
                        currentConfig=config_i,
                    )

                # compute lift coefficient
                CL = AC.CL(tas=TAS_i, sigma=sigma, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(
                    CL=CL, config=config_i, expedite=expedite, speedBrakes=speedBrakes
                )
                # compute drag force
                Drag = AC.D(tas=TAS_i, sigma=sigma, CD=CD)
                # compute ROCD
                ROCD_i = (
                    conv.m2ft(
                        AC.ROCD(
                            T=THR_i,
                            D=Drag,
                            v=TAS_i,
                            mass=mass_i,
                            ESF=ESF_i,
                            h=H_m,
                            DeltaTemp=DeltaTemp,
                            reducedPower=reducedPower,
                        )
                    )
                    * 60
                )

            # Compute elapsed time and fuel burn over current step
            if Hp_i == Hp_init:
                break
            else:
                # Average ROCD over step is the mean of initial and final ones
                step_ROCD = (ROCD[-1] + ROCD_i) / 2  # [ft/min]
                # Step time is: altitude differential divided by average ROCD
                step_time = 60 * (Hp_i - Hp[-1]) / step_ROCD  # [s]

                # BADAE
                if AC.BADAFamily.BADAE:
                    # SOC change is: average SOC rate multiplied by step time
                    step_SOC = step_SOCr * step_time / 3600  # [%]
                    # Update SOC estimate at end of step
                    SOC_i = SOC[-1] - step_SOC  # [%]
                    # update of aircraft mass estimate at end of step - mass is not changing for ELECTRIC engine (no fuel is consumed)
                    mass_i = mass[-1]  # [kg]

                else:
                    # fuel burnt is: average fuel flow multiplied by step time
                    step_mass = step_FUEL * step_time  # [kg]
                    # update of aircraft mass estimate at end of step
                    if not mass_const:
                        mass_i = mass[-1] - step_mass  # [kg]
                        fuelConsumed_i = step_FUEL * step_time
                    fuelConsumed_i = FUELCONSUMED[-1] + step_FUEL * step_time

        ## PART 3: store information about end of step point

        writeOutputData = True
        if phase == "Climb" and ROCD_i < 0:
            warnings.warn(
                "Value ROCD = "
                + str(ROCD_i)
                + " [ft/min] is negative at the altitude "
                + str(Hp_i)
                + " [ft].",
                UserWarning,
            )
            go_on = False
            writeOutputData = False

        elif phase == "Climb" and ROCD_i < ROCD_min:
            warnings.warn(
                "Value ROCD = "
                + str(ROCD_i)
                + " [ft/min] exceeds the service ceiling limit defined by minimum ROCD = "
                + str(ROCD_min)
                + " [ft/min] at the altitude "
                + str(Hp_i)
                + " [ft].",
                UserWarning,
            )
            go_on = False
            writeOutputData = True

        if writeOutputData:
            # point data
            Hp.append(Hp_i)
            TAS.append(conv.ms2kt(TAS_i))
            CAS.append(conv.ms2kt(CAS_i))
            M.append(M_i)
            ROCD.append(ROCD_i)
            esf.append(ESF_i)
            Comment.append(comment)

            # everything except electric BADAE
            if not AC.BADAFamily.BADAE:
                FUEL.append(FUEL_i)

            # BADAH
            if AC.BADAFamily.BADAH:
                Peng.append(Preq_target_i)
                Preq.append(Preq_i)
                Pav.append(Pav_i)

            # BADAE
            elif AC.BADAFamily.BADAE:
                Pmec.append(Preq_target_i)
                Pbat.append(Pbat_i)
                SOCr.append(SOCr_i)
                Pelc.append(Pelc_i)
                Ibat.append(Ibat_i)
                Vbat.append(Vbat_i)
                Vgbat.append(Vgbat_i)

            # BADA3 & BADA4
            elif AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
                THR.append(THR_i)
                DRAG.append(Drag)
                config.append(config_i)

            # BADA4
            if AC.BADAFamily.BADA4:
                HLid.append(HLid_i)
                LG.append(LG_i)

            # calculation of the slope
            if TAS_i == 0:
                gamma_i = 90 * np.sign(ROCD_i)
            else:
                if AC.BADAFamily.BADAE:
                    gamma_i = degrees(atan(conv.ft2m(ROCD_i) * temp_const / 60 / TAS_i))
                else:
                    # using SIN assumes the TAS to be in the direction of the aircraft axis, not ground plane. Which means, this should be mathematically the correct equation for all the aircraft
                    gamma_i = degrees(asin(conv.ft2m(ROCD_i) * temp_const / 60 / TAS_i))

            # ground speed can be calcualted as TAS projected on the ground minus wind
            GS_i = cos(radians(gamma_i)) * TAS_i - wS
            GS.append(conv.ms2kt(GS_i))

            Slope.append(gamma_i)
            acc.append(0.0)

            BankAngle.append(bankAngle)
            ROT.append(rateOfTurn)

            # integrated data
            if Hp_i != Hp_init:
                if AC.BADAFamily.BADAE:
                    SOC.append(SOC_i)

                mass.append(mass_i)
                time.append(time[-1] + step_time)

                # everything except electric BADAE
                if not AC.BADAFamily.BADAE:
                    FUELCONSUMED.append(fuelConsumed_i)

                # Average TAS over step is the mean of initial and final ones
                step_TAS = (TAS[-2] + TAS[-1]) / 2  # [kt]
                # Average slope over the step
                step_gamma = radians((Slope[-2] + Slope[-1]) / 2)  # radians
                # Average ground speed over step
                # since this is not level flight, TAS speed should be projected on the ground, then GS can be calculated applying the wind speed
                step_TAS_projected = cos(step_gamma) * step_TAS
                step_GS = step_TAS_projected - wS  # [kt]
                # Step distance is: average GS multiplied by step time
                if turnFlight:
                    step_distance = conv.m2nm(
                        turn.distance(
                            rateOfTurn=rateOfTurn, TAS=TAS_i, timeOfTurn=step_time
                        )
                    )  # arcLength during the turn [NM]
                else:
                    step_distance = step_GS * step_time / 3600  # [NM]
                # Distance at end of step is distance at start of step plus step distance
                dist.append(dist[-1] + step_distance)

                # add GPS calculation
                if Lat and Lon and (magneticHeading or trueHeading):
                    if headingToFly == "TRUE":
                        if not turnFlight:
                            if not constantHeading:
                                # fly ORTHODROME
                                (Lat_i, Lon_i, HDGTrue_i) = (
                                    vincenty.destinationPoint_finalBearing(
                                        LAT_init=LAT[-1],
                                        LON_init=LON[-1],
                                        distance=conv.nm2m(step_distance),
                                        bearing=HDGTrue[-1],
                                    )
                                )

                                if magneticDeclinationGrid is not None:
                                    HDGMagnetic_i = (
                                        HDGTrue_i
                                        - magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGMagnetic_i = HDGTrue_i

                            elif constantHeading:
                                # fly LOXODROME
                                (Lat_i, Lon_i) = rhumb.destinationPoint(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearing=HDGTrue[-1],
                                    distance=conv.nm2m(step_distance),
                                )
                                HDGTrue_i = HDGTrue[-1]

                                if magneticDeclinationGrid is not None:
                                    HDGMagnetic_i = (
                                        HDGTrue_i
                                        - magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGMagnetic_i = HDGTrue_i

                        else:
                            # calculate the turn
                            (Lat_i, Lon_i, HDGTrue_i) = (
                                turn.destinationPoint_finalBearing(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearingInit=HDGTrue[-1],
                                    TAS=TAS_i,
                                    rateOfTurn=rateOfTurn,
                                    timeOfTurn=step_time,
                                    directionOfTurn=directionOfTurn,
                                )
                            )

                            if magneticDeclinationGrid is not None:
                                HDGMagnetic_i = (
                                    HDGTrue_i
                                    - magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGMagnetic_i = HDGTrue_i

                    elif headingToFly == "MAGNETIC":
                        if not turnFlight:
                            if constantHeading:
                                (Lat_i, Lon_i) = rhumb.destinationPoint(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearing=HDGTrue[-1],
                                    distance=conv.nm2m(step_distance),
                                )
                                HDGMagnetic_i = HDGMagnetic[-1]

                                if magneticDeclinationGrid is not None:
                                    HDGTrue_i = (
                                        HDGMagnetic_i
                                        + magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGTrue_i = HDGMagnetic_i

                        else:
                            # calculate the turn
                            (Lat_i, Lon_i, HDGTrue_i) = (
                                turn.destinationPoint_finalBearing(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearingInit=HDGTrue[-1],
                                    TAS=TAS_i,
                                    rateOfTurn=rateOfTurn,
                                    timeOfTurn=step_time,
                                    directionOfTurn=directionOfTurn,
                                )
                            )

                            if magneticDeclinationGrid is not None:
                                HDGMagnetic_i = (
                                    HDGTrue_i
                                    - magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGMagnetic_i = HDGTrue_i

                    LAT.append(Lat_i)
                    LON.append(Lon_i)
                    HDGMagnetic.append(HDGMagnetic_i)
                    HDGTrue.append(HDGTrue_i)

            # Determine end altitude of next step
            Hp_next = Hp_i + Hp_step

            if phase == "Climb":
                if Hp_next < Hp_final:
                    Hp_i = Hp_next - (Hp_i % Hp_step)
                # remaining altitude step would cross over the final altitude
                elif Hp_i < Hp_final:
                    Hp_i = Hp_final
                else:
                    go_on = False
            else:
                if Hp_next > Hp_final:
                    Hp_i = Hp_next - (Hp_i % Hp_step)
                # remaining altitude step would cross over the final altitude
                elif Hp_i > Hp_final:
                    Hp_i = Hp_final
                else:
                    go_on = False

    flightData = {
        "Hp": Hp,
        "TAS": TAS,
        "CAS": CAS,
        "GS": GS,
        "M": M,
        "acc": acc,
        "ROCD": ROCD,
        "ESF": esf,
        "FUEL": FUEL,
        "Pmec": Pmec,
        "Pelc": Pelc,
        "Pbat": Pbat,
        "SOCr": SOCr,
        "SOC": SOC,
        "Ibat": Ibat,
        "Vbat": Vbat,
        "Vgbat": Vgbat,
        "FUELCONSUMED": FUELCONSUMED,
        "Preq": Preq,
        "Peng": Peng,
        "Pav": Pav,
        "THR": THR,
        "DRAG": DRAG,
        "time": time,
        "dist": dist,
        "slope": Slope,
        "mass": mass,
        "config": config,
        "HLid": HLid,
        "LG": LG,
        "LAT": LAT,
        "LON": LON,
        "HDGTrue": HDGTrue,
        "HDGMagnetic": HDGMagnetic,
        "BankAngle": BankAngle,
        "ROT": ROT,
        "comment": Comment,
    }

    flightTrajectory = FT.createFlightTrajectoryDataframe(flightData)
    return flightTrajectory


def constantSpeedRating_time(
    AC,
    length,
    speedType,
    v,
    Hp_init,
    phase,
    m_init,
    DeltaTemp,
    wS=0.0,
    turnMetrics={"rateOfTurn": 0.0, "bankAngle": 0.0, "directionOfTurn": None},
    Lat=None,
    Lon=None,
    initialHeading={"magnetic": None, "true": None, "constantHeading": None},
    reducedPower=None,
    directionOfTurn=None,
    expedite=False,
    magneticDeclinationGrid=None,
    initRating=None,
    **kwargs,
):
    """This function computes time and fuel required by an aircraft to perform a climb/descent from Hp_init for set amount of time at constant speed and constant engine rating

    :param AC: aircraft {BADA3/4/H/E}
    :param speedType: what kind of speed is followed {M, CAS, TAS}.
    :param length: length of a segment to fly [s]
    :param step_length: length of a step of a segment - [s]
    :param v: what kind of speed is followed - [kt] CAS/TAS speed to follow or [-] MACH speed to follow.
    :param Hp_init: initial pressure altitude [ft].
    :param phase: phase of flight {Climb, Descent}
    :param m_init: initial aircraft mass [kg].
    :param DeltaTemp: deviation with respect to ISA [K].
    :param wS: longitudinal wind speed (TAS) [kt].
    :param turnMetrics: Metrics for turn performance {"rateOfTurn":0.0,"bankAngle":0.0,"directionOfTurn":None} {[deg/s],[deg],[LEFT/RIGHT]}
    :param SOC_init: initial state of charge [%].
    :param config: aircraft default aerodynamic configuration {TO,IC,CR,AP,LD}.
    :param speedBrakes: deployed or not speedbrakes including value to be added to the drag coeffcient {deployed:False,value:0.03} {deployed:[True/False],value:[-]}.
    :param ROCD_min: lower ROCD threshold to identify the climbing capabilities (service ceiling) [ft/min].
    :param Lat: Geographical Latitude [deg]
    :param Lon: Geographical Longitude [deg]
    :param initialHeading: aircraft magnetic heading, true heading and definition of constant heading(ORTHODROME=False, LOXODROME=True) {[deg],[deg],-}
    :param magneticDeclinationGrid: geographical grid of a magnetic declination on Earth [deg]
    :param mass_const: kind of mass canculation {mass_integrated=False, mass_constant=True}.
    :param m_iter: number of iterations for integration loop [-]
    :param reducedPower: reduction of Power during the climb {True/False}
    :param initRating: default rating settings
    :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}.
    :type speedType: string.
    :type length: float.
    :type step_length: float.
    :type v: float.
    :type Hp_init: float.
    :type phase: string.
    :type m_init: float.
    :type DeltaTemp: float.
    :type wS: float.
    :type turnMetrics: {float,float,string}.
    :type SOC_init: float.
    :type config: string.
    :type speedBrakes: dict{boolean,float}.
    :type ROCD_min: float.
    :type Lat: float.
    :type Lon: float.
    :type initialHeading: {float,float,boolean}.
    :type magneticDeclinationGrid: magneticDeclinationGrid.
    :type mass_const: boolean.
    :type m_iter: integer.
    :type reducedPower: boolean.
    :type initRating: string
    :returns:
            BADA3: [Hp, TAS, CAS, M, GS, acc, ROCD, ESF, FUEL, FUELCONSUMED, THR, time, dist, slope, mass, config, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, N, s, NM, deg, kg, -, -,deg,deg,deg,deg,deg,deg/s]
            BADA4: [Hp, TAS, CAS, M, GS, acc, ROCD, ESF, FUEL, FUELCONSUMED, THR, P[Pmec, Pbat, Pelc Ibat, Vbat, Vgbat, SOCr, SOC], time, dist, slope, mass, config, HLid, LG, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, N, [W,W,W,A,V,V,%/h,%], s, NM, deg, kg, -, -, -, -,deg,deg,deg,deg,deg,deg/s]
            BADAH: [Hp, TAS, CAS, M, GS, ROCD, ESF, FUEL, FUELCONSUMED, Peng, Preq, Pav, time, dist, slope, mass, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, W, W, W, s, NM, deg, kg, -,deg,deg,deg,deg,deg,deg/s]
            BADAE: [time, dist, Hp, TAS, CAS, M, GS, acc, ROCD, ESF, slope, mass, P[Pmec, Pelc, Pbat, SOCr, SOC, Ibat, Vbat, Vgbat] comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [s, NM, ft, kt, kt, -, kt, m/s^2, ft/min, deg, kg, [W,W,W,%/h,%,A,V,V], -,deg,deg,deg,deg,deg,deg/s]
    :rtype: dict[list[float]}.
    """

    rateOfTurn = turnMetrics["rateOfTurn"]
    bankAngle = turnMetrics["bankAngle"]
    directionOfTurn = turnMetrics["directionOfTurn"]

    turnFlight = True
    if turnMetrics["rateOfTurn"] == 0.0 and turnMetrics["bankAngle"] == 0.0:
        turnFlight = False

    # conversion of Magnetic Heading to True Heading
    if magneticDeclinationGrid is not None:
        magneticDeclination = magneticDeclinationGrid.getMagneticDeclination(
            LAT_target=Lat, LON_target=Lon
        )
    else:
        magneticDeclination = 0

    # retrieve magnetic and true heading inputs
    magneticHeading = initialHeading["magnetic"]
    trueHeading = initialHeading["true"]
    constantHeading = initialHeading["constantHeading"]

    if Lat and Lon and (magneticHeading or trueHeading):
        if trueHeading is not None and magneticHeading is None:
            # fly TRUE Heading
            headingToFly = "TRUE"
            magneticHeading = trueHeading - magneticDeclination
        elif magneticHeading is not None and trueHeading is None:
            # fly MAGNETIC Heading
            if constantHeading == True:
                headingToFly = "MAGNETIC"
                trueHeading = magneticHeading + magneticDeclination
            else:
                raise Exception("Cannot fly non-constant magnetic heading")

        else:
            raise Exception("Undefined Heading value combination")

    # calculation with constant mass (True) or integrated (False)
    mass_const = kwargs.get("mass_const", False)

    # optional parameter to define initial Baterry State of Charge (SOC)
    if AC.BADAFamily.BADAE:
        SOC_init = kwargs.get("SOC_init", 100)
    else:
        SOC_init = None

    # speed brakes application
    if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
        speedBrakes = kwargs.get("speedBrakes", {"deployed": False, "value": 0.03})

    # step size in [s]
    step_length = kwargs.get("step_length", 1)

    # minimum remaining ROCD to determine cruise ceiling
    if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
        ROCD_min = kwargs.get("ROCD_min", 50)  # [ft/min]
    else:
        if AC.engineType == "PISTON" or AC.engineType == "ELECTRIC":
            ROCD_min = kwargs.get("ROCD_min", 100)  # [ft/min]
        else:
            ROCD_min = kwargs.get("ROCD_min", 300)  # [ft/min]

    # determine if  vertical evolution over the segment is CLIMB or DESCENT
    # and associate engine rating and altitude iteration direction

    if initRating is None:
        if phase == "Climb":
            if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
                if v == 0:
                    rating = "MTKF"
                else:
                    rating = "MCNT"
            else:
                rating = "MCMB"
        elif phase == "Descent":
            if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
                if v == 0:
                    rating = "UNKNOWN"
                else:
                    rating = "UNKNOWN"
            else:
                rating = "LIDL"
        else:
            raise Exception("Phase definition is wrong! It should be Climb or Descent")
    else:
        rating = initRating

    if turnFlight:
        turnComment = "_turn"
    else:
        turnComment = ""

    if constantHeading:
        constHeadingStr = "_const_Heading"
    elif constantHeading is False or constantHeading is None:
        constHeadingStr = ""

    # comment line describing type of trajectory calculation
    comment = (
        phase + turnComment + "_const_" + speedType + "_" + rating + constHeadingStr
    )

    if Lat and Lon and (magneticHeading or trueHeading):
        comment = comment + "_" + headingToFly + "_Heading"

    if expedite:
        comment = comment + "_expedite"

    # get the default aerodynamic configuration if provided to be used for the whole segment
    config_default = kwargs.get("config", None)
    if config_default is not None:
        if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
            if not (
                config_default == "TO"
                or config_default == "IC"
                or config_default == "CR"
                or config_default == "AP"
                or config_default == "LD"
            ):
                print(
                    "WRONG default configuration set. Available values are: TO/IC/CR/AP/LD. Configuration will be calculated automatically"
                )

    #  weight iteration constant
    m_iter = kwargs.get("m_iter", 5)  # number of iterations for integration loop[-]

    # The thrust_fuel method for BADA 3 models applies the cruise fuel correction
    # whenever the thrust is adapted, instead of only in cruise: this correction
    # needs to be reverted when thrust is adapted for constant ROC/slope.

    # cruise_correction = 1/f(5)

    # initialize output parameters
    Hp = [Hp_init]
    TAS = []
    CAS = []
    GS = []
    M = []
    ROCD = []
    esf = []
    FUEL = []
    FUELCONSUMED = []
    time = [0]
    dist = [0]
    mass = [m_init]
    Comment = []
    Slope = []
    acc = []
    THR = []
    DRAG = []
    config = []
    HLid = []
    LG = []
    BankAngle = []
    ROT = []

    if not AC.BADAFamily.BADAE:
        FUELCONSUMED = [0]

    # BADAH specific
    Peng = []
    Preq = []
    Pav = []

    # optional GPS coordiantes and HDG definition
    if Lat and Lon and (magneticHeading or trueHeading):
        LAT = [Lat]
        LON = [Lon]
        HDGMagnetic = [magneticHeading]
        HDGTrue = [trueHeading]
    else:
        LAT = []
        LON = []
        HDGMagnetic = []
        HDGTrue = []

    # BADAE specific
    Pmec = []
    Pbat = []
    SOCr = []
    SOC = [SOC_init]
    Pelc = []
    Ibat = []
    Vbat = []
    Vgbat = []

    # init loop parameters
    length_loop = 0
    time_i = time[-1]
    go_on = True

    while go_on:
        ## PART 1: compute parameters at end of step that are known without uncertainties:
        ##         atmosphere, speeds, thrust. fuel flow, ESF

        ## PART 2: compute parameters at end of step that are known only with uncertainties
        ##           (due to unknown mass at end of step):
        ##         weight, lift, drag , ROCD

        mass_i = mass[-1]
        Hp_i = Hp[-1]
        for _ in itertools.repeat(None, m_iter):
            # atmosphere properties
            H_m = conv.ft2m(Hp_i)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
            temp_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)

            # aircraft speed
            [M_i, CAS_i, TAS_i] = atm.convertSpeed(
                v=v, speedType=speedType, theta=theta, delta=delta, sigma=sigma
            )

            if turnFlight:
                if turnMetrics["bankAngle"] != 0.0:
                    # bankAngle is defined
                    rateOfTurn = AC.rateOfTurn_bankAngle(TAS=TAS_i, bankAngle=bankAngle)
                else:
                    # rateOfTurn is defined
                    bankAngle = AC.bankAngle(
                        rateOfTurn=rateOfTurn, v=TAS_i
                    )  # [degrees]

            # Load factor
            nz = 1 / cos(radians(bankAngle))

            # compute Energy Share Factor (ESF)
            ESF_i = AC.esf(
                h=H_m, M=M_i, DeltaTemp=DeltaTemp, flightEvolution=("const" + speedType)
            )

            step_time = length_loop - time[-1]

            # BADAH
            if AC.BADAFamily.BADAH:
                # compute available power
                if rating == "UNKNOWN":
                    Preq_target_i = (
                        0.1 * AC.P0
                    )  # No minimum power model: assume 10% torque
                else:
                    Preq_target_i = AC.Pav(rating=rating, theta=theta, delta=delta)

                Pav_i = AC.Pav(rating="MTKF", theta=theta, delta=delta)

                # compute fuel flow for level flight
                CP = AC.CP(Peng=Preq_target_i)
                FUEL_i = AC.ff(delta=delta, CP=CP)  # [kg/s]

            # BADAE
            elif AC.BADAFamily.BADAE:
                # compute available power
                if rating == "UNKNOWN":
                    Preq_target_i = (
                        0.1 * AC.P0
                    )  # No minimum power model: assume 10% torque
                else:
                    Preq_target_i = AC.Pav(rating=rating, SOC=SOC[-1])

                Pbat_i = AC.Pbat(Preq=Preq_target_i, SOC=SOC[-1])
                SOCr_i = AC.SOCrate(Preq=Preq_target_i, SOC=SOC[-1])

                # debug data
                Pelc_i = Preq_target_i / AC.eta
                Ibat_i = AC.Ibat(P=Pelc_i, SOC=SOC[-1])
                Vbat_i = AC.Vbat(I=Ibat_i, SOC=SOC[-1])
                Vgbat_i = AC.Vocbat(SOC=SOC[-1]) - AC.R0bat(SOC=SOC[-1]) * Ibat_i

            # BADA4
            elif AC.BADAFamily.BADA4:
                # compute thrust force and fuel flow
                THR_i = AC.Thrust(
                    rating=rating, delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )  # [N]
                CT = AC.CT(Thrust=THR_i, delta=delta)
                FUEL_i = AC.ff(
                    CT=CT, delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                )

            # BADA3
            elif AC.BADAFamily.BADA3:
                # aircraft aerodynamic configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=phase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=phase,
                        previousConfig=config[-1],
                        currentConfig=config_i,
                    )

                # compute thrust force and fuel flow
                THR_i = AC.Thrust(
                    rating=rating, v=TAS_i, h=H_m, config=config_i, DeltaTemp=DeltaTemp
                )
                FUEL_i = AC.ff(
                    flightPhase=phase,
                    v=TAS_i,
                    h=H_m,
                    T=THR_i,
                    config=config_i,
                    adapted=False,
                )  # MCMB(climb) or IDLE(descent)

            # BADAH or BADAE
            if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
                # compute Power required
                Preq_i = AC.Preq(sigma=sigma, tas=TAS_i, mass=mass_i, phi=bankAngle)
                # compute ROCD
                ROCD_i = (
                    conv.m2ft(
                        AC.ROCD(
                            Peng=Preq_target_i,
                            Preq=Preq_i,
                            mass=mass_i,
                            ESF=ESF_i,
                            theta=theta,
                            DeltaTemp=DeltaTemp,
                        )
                    )
                    * 60
                )

            # BADA4
            elif AC.BADAFamily.BADA4:
                # aircraft aerodynamic configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=phase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=phase, previousConfig=config[-1], currentConfig=config_i
                    )

                [HLid_i, LG_i] = AC.flightEnvelope.getAeroConfig(config=config_i)

                # compute lift coefficient
                CL = AC.CL(M=M_i, delta=delta, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(M=M_i, CL=CL, HLid=HLid_i, LG=LG_i, speedBrakes=speedBrakes)
                # compute drag force
                Drag = AC.D(M=M_i, delta=delta, CD=CD)

                ROCD_i = (
                    conv.m2ft(
                        (1 / temp_const)
                        * (THR_i - Drag)
                        * TAS_i
                        * ESF_i
                        / (mass_i * const.g)
                    )
                    * 60
                )

            # BADA3
            elif AC.BADAFamily.BADA3:
                # aircraft aerodynamic configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=phase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=phase,
                        previousConfig=config[-1],
                        currentConfig=config_i,
                    )

                # compute lift coefficient
                CL = AC.CL(tas=TAS_i, sigma=sigma, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(
                    CL=CL, config=config_i, expedite=expedite, speedBrakes=speedBrakes
                )
                # compute drag force
                Drag = AC.D(tas=TAS_i, sigma=sigma, CD=CD)
                # compute ROCD
                ROCD_i = (
                    conv.m2ft(
                        AC.ROCD(
                            T=THR_i,
                            D=Drag,
                            v=TAS_i,
                            mass=mass_i,
                            ESF=ESF_i,
                            h=H_m,
                            DeltaTemp=DeltaTemp,
                            reducedPower=reducedPower,
                        )
                    )
                    * 60
                )

            # Compute elapsed time and fuel burn over current step
            if length_loop == 0:
                # no need to loop for first point: initial m/Hp already known
                break
            else:
                # Average ROCD over step is the mean of initial and final ones
                step_ROCD = (ROCD[-1] + ROCD_i) / 2  # [ft/min]
                # Altitude differential is: average ROCD multiplied by step time
                step_Hp = step_ROCD * step_time / 60  # [ft]
                # Update altitude estimate at end of step
                Hp_i = Hp[-1] + step_Hp  # [ft]

                # BADAE
                if AC.BADAFamily.BADAE:
                    # Average SOC rate over step is the mean of initial and final ones
                    step_SOCr = (SOCr[-1] + SOCr_i) / 2  # [%/h]
                    # SOC change is: average SOC rate multiplied by step time
                    step_SOC = step_SOCr * step_time / 3600  # [%]
                    # Update SOC estimate at end of step
                    SOC_i = SOC[-1] - step_SOC  # [%]
                    # update of aircraft mass estimate at end of step - mass is not changing for ELECTRIC engine (no fuel is consumed)
                    mass_i = mass[-1]  # [kg]
                else:
                    # average fuel flow over step is the mean of initial and final ones
                    step_FUEL = (FUEL[-1] + FUEL_i) / 2  # [kg/s]
                    # fuel burnt is: average fuel flow multiplied by step time
                    step_mass = step_FUEL * step_time  # [kg]
                    # update of aircraft mass estimate at end of step
                    if not mass_const:
                        mass_i = mass[-1] - step_mass  # [kg]
                        fuelConsumed_i = step_FUEL * step_time
                    fuelConsumed_i = FUELCONSUMED[-1] + step_FUEL * step_time

        ## PART 3: store information about end of step point

        writeOutputData = True
        if phase == "Climb" and ROCD_i < 0:
            warnings.warn(
                "Value ROCD = "
                + str(ROCD_i)
                + " [ft/min] is negative at the altitude "
                + str(Hp_i)
                + " [ft].",
                UserWarning,
            )
            go_on = False
            writeOutputData = False

        elif phase == "Climb" and ROCD_i < ROCD_min:
            warnings.warn(
                "Value ROCD = "
                + str(ROCD_i)
                + " [ft/min] exceeds the service ceiling limit defined by minimum ROCD = "
                + str(ROCD_min)
                + " [ft/min] at the altitude "
                + str(Hp_i)
                + " [ft].",
                UserWarning,
            )
            go_on = False
            writeOutputData = True

        if writeOutputData:
            # point data
            TAS.append(conv.ms2kt(TAS_i))
            CAS.append(conv.ms2kt(CAS_i))
            M.append(M_i)
            ROCD.append(ROCD_i)
            esf.append(ESF_i)
            Comment.append(comment)

            # everything except electric BADAE
            if not AC.BADAFamily.BADAE:
                FUEL.append(FUEL_i)

            # BADAH
            if AC.BADAFamily.BADAH:
                Peng.append(Preq_target_i)
                Preq.append(Preq_i)
                Pav.append(Pav_i)

            # BADAE
            elif AC.BADAFamily.BADAE:
                Pmec.append(Preq_target_i)
                Pbat.append(Pbat_i)
                SOCr.append(SOCr_i)
                Pelc.append(Pelc_i)
                Ibat.append(Ibat_i)
                Vbat.append(Vbat_i)
                Vgbat.append(Vgbat_i)

            # BADA3 & BADA4
            elif AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
                THR.append(THR_i)
                DRAG.append(Drag)
                config.append(config_i)

            # BADA4
            if AC.BADAFamily.BADA4:
                HLid.append(HLid_i)
                LG.append(LG_i)

            # calculation of the slope
            if TAS_i == 0:
                gamma_i = 90 * np.sign(ROCD_i)
            else:
                [theta, delta, sigma] = atm.atmosphereProperties(
                    h=conv.ft2m(Hp_i), DeltaTemp=DeltaTemp
                )
                temp_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)
                if AC.BADAFamily.BADAE:
                    gamma_i = degrees(atan(conv.ft2m(ROCD_i) * temp_const / 60 / TAS_i))
                else:
                    # using SIN assumes the TAS to be in the direction of the aircraft axis, not ground plane. Which means, this should be mathematically the correct equation for all the aircraft
                    gamma_i = degrees(asin(conv.ft2m(ROCD_i) * temp_const / 60 / TAS_i))

            # ground speed can be calcualted as TAS projected on the ground minus wind
            GS_i = cos(radians(gamma_i)) * TAS_i - wS
            GS.append(conv.ms2kt(GS_i))

            Slope.append(gamma_i)
            acc.append(0.0)

            BankAngle.append(bankAngle)
            ROT.append(rateOfTurn)

            # integrated data
            if length_loop != 0:
                if AC.BADAFamily.BADAE:
                    SOC.append(SOC_i)

                Hp.append(Hp_i)
                mass.append(mass_i)
                time.append(time[-1] + step_time)

                # everything except electric BADAE
                if not AC.BADAFamily.BADAE:
                    FUELCONSUMED.append(fuelConsumed_i)

                # Average TAS over step is the mean of initial and final ones
                step_TAS = (TAS[-2] + TAS[-1]) / 2  # [kt]
                # Average slope over the step
                step_gamma = radians((Slope[-2] + Slope[-1]) / 2)  # radians
                # Average ground speed over step
                # since this is not level flight, TAS speed should be projected on the ground, then GS can be calculated applying the wind speed
                step_TAS_projected = cos(step_gamma) * step_TAS
                step_GS = step_TAS_projected - wS  # [kt]
                # Step distance is: average GS multiplied by step time
                if turnFlight:
                    step_distance = conv.m2nm(
                        turn.distance(
                            rateOfTurn=rateOfTurn, TAS=TAS_i, timeOfTurn=step_time
                        )
                    )  # arcLength during the turn [NM]
                else:
                    step_distance = step_GS * step_time / 3600  # [NM]
                # Distance at end of step is distance at start of step plus step distance
                dist.append(dist[-1] + step_distance)

                # add GPS calculation
                if Lat and Lon and (magneticHeading or trueHeading):
                    if headingToFly == "TRUE":
                        if not turnFlight:
                            if not constantHeading:
                                # fly ORTHODROME
                                (Lat_i, Lon_i, HDGTrue_i) = (
                                    vincenty.destinationPoint_finalBearing(
                                        LAT_init=LAT[-1],
                                        LON_init=LON[-1],
                                        distance=conv.nm2m(step_distance),
                                        bearing=HDGTrue[-1],
                                    )
                                )

                                if magneticDeclinationGrid is not None:
                                    HDGMagnetic_i = (
                                        HDGTrue_i
                                        - magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGMagnetic_i = HDGTrue_i

                            elif constantHeading:
                                # fly LOXODROME
                                (Lat_i, Lon_i) = rhumb.destinationPoint(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearing=HDGTrue[-1],
                                    distance=conv.nm2m(step_distance),
                                )
                                HDGTrue_i = HDGTrue[-1]

                                if magneticDeclinationGrid is not None:
                                    HDGMagnetic_i = (
                                        HDGTrue_i
                                        - magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGMagnetic_i = HDGTrue_i

                        else:
                            # calculate the turn
                            (Lat_i, Lon_i, HDGTrue_i) = (
                                turn.destinationPoint_finalBearing(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearingInit=HDGTrue[-1],
                                    TAS=TAS_i,
                                    rateOfTurn=rateOfTurn,
                                    timeOfTurn=step_time,
                                    directionOfTurn=directionOfTurn,
                                )
                            )

                            if magneticDeclinationGrid is not None:
                                HDGMagnetic_i = (
                                    HDGTrue_i
                                    - magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGMagnetic_i = HDGTrue_i

                    elif headingToFly == "MAGNETIC":
                        if not turnFlight:
                            if constantHeading:
                                (Lat_i, Lon_i) = rhumb.destinationPoint(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearing=HDGTrue[-1],
                                    distance=conv.nm2m(step_distance),
                                )
                                HDGMagnetic_i = HDGMagnetic[-1]

                                if magneticDeclinationGrid is not None:
                                    HDGTrue_i = (
                                        HDGMagnetic_i
                                        + magneticDeclinationGrid.getMagneticDeclination(
                                            LAT_target=Lat_i, LON_target=Lon_i
                                        )
                                    )
                                else:
                                    magneticDeclination = 0
                                    HDGTrue_i = HDGMagnetic_i

                        else:
                            # calculate the turn
                            (Lat_i, Lon_i, HDGTrue_i) = (
                                turn.destinationPoint_finalBearing(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    bearingInit=HDGTrue[-1],
                                    TAS=TAS_i,
                                    rateOfTurn=rateOfTurn,
                                    timeOfTurn=step_time,
                                    directionOfTurn=directionOfTurn,
                                )
                            )

                            if magneticDeclinationGrid is not None:
                                HDGMagnetic_i = (
                                    HDGTrue_i
                                    - magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGMagnetic_i = HDGTrue_i

                    LAT.append(Lat_i)
                    LON.append(Lon_i)
                    HDGMagnetic.append(HDGMagnetic_i)
                    HDGTrue.append(HDGTrue_i)

            if length_loop + step_length < length:
                length_loop += step_length
            elif length_loop < length:
                length_loop = length
            else:
                go_on = False

    flightData = {
        "Hp": Hp,
        "TAS": TAS,
        "CAS": CAS,
        "GS": GS,
        "M": M,
        "acc": acc,
        "ROCD": ROCD,
        "ESF": esf,
        "FUEL": FUEL,
        "Pmec": Pmec,
        "Pelc": Pelc,
        "Pbat": Pbat,
        "SOCr": SOCr,
        "SOC": SOC,
        "Ibat": Ibat,
        "Vbat": Vbat,
        "Vgbat": Vgbat,
        "FUELCONSUMED": FUELCONSUMED,
        "Preq": Preq,
        "Peng": Peng,
        "Pav": Pav,
        "THR": THR,
        "DRAG": DRAG,
        "time": time,
        "dist": dist,
        "slope": Slope,
        "mass": mass,
        "config": config,
        "HLid": HLid,
        "LG": LG,
        "LAT": LAT,
        "LON": LON,
        "HDGTrue": HDGTrue,
        "HDGMagnetic": HDGMagnetic,
        "BankAngle": BankAngle,
        "ROT": ROT,
        "comment": Comment,
    }

    flightTrajectory = FT.createFlightTrajectoryDataframe(flightData)
    return flightTrajectory


def accDec(
    AC,
    speedType,
    v_init,
    v_final,
    phase,
    Hp_init,
    m_init,
    DeltaTemp,
    wS=0.0,
    turnMetrics={"rateOfTurn": 0.0, "bankAngle": 0.0, "directionOfTurn": None},
    control=None,
    Lat=None,
    Lon=None,
    initialHeading={"magnetic": None, "true": None, "constantHeading": None},
    reducedPower=None,
    magneticDeclinationGrid=None,
    **kwargs,
):
    """This function computes time and fuel required by an aircraft to perform an acceleration/deceleration from v_init to v_final in climb cruise or descent

    .. note::
            The control law used during the segment depends on the targets provided in input parameter 'control':
            - ROCD/slope+ESF:  law is ROCD/slope+ESF
            - ROCD/slope+acc:  law is ROCD/slope+acc
            - ROCD/slope only: law is rating+ROCD/slope
            - ESF only:        law is rating+ESF
            - acc only:        law is rating+acc
            - Neither:         law is rating+default ESF

    :param AC: aircraft {BADA3/4/H/E}
    :param speedType: what kind of speed is followed {M, CAS, TAS}.
    :param v_init: initial speed to follow - [kt] CAS/TAS speed to follow or [-] MACH speed to follow.
    :param v_final: final speed to follow - [kt] CAS/TAS speed to follow or [-] MACH speed to follow.
    :param phase: vertical evolution {Climb, Descent, Cruise}
    :param control: structure containing a combination of the following targets:
            :param ROCDtarget: Rate of climb/descent to be followed [ft/min].
            :param slopetarget: slope (flight path angle) to be followed [deg].
            :param acctarget: acceleration to be followed [m/s^2].
            :param ESFtarget: Energy Share Factor to be followed [-].
    :param maxRating: rating to be used as a limit on the maximum thrust/power [-].
    :param Hp_init: initial pressure altitude [ft].
    :param m_init: initial aircraft mass [kg].
    :param DeltaTemp: deviation with respect to ISA [K].
    :param wS: longitudinal wind speed (TAS) [kt].
    :param turnMetrics: Metrics for turn performance {"rateOfTurn":0.0,"bankAngle":0.0,"directionOfTurn":None} {[deg/s],[deg],[LEFT/RIGHT]}
    :param SOC_init: initial state of charge [%].
    :param config: aircraft default aerodynamic configuration {TO,IC,CR,AP,LD}.
    :param speedBrakes: deployed or not speedbrakes including value to be added to the drag coeffcient {deployed:False,value:0.03} {deployed:[True/False],value:[-]}.
    :param speed_step: step of the speed for the speed iteration loop [-] for M, [kt] for TAS or CAS
    :param Lat: Geographical Latitude [deg]
    :param Lon: Geographical Longitude [deg]
    :param initialHeading: aircraft magnetic heading, true heading and definition of constant heading(ORTHODROME=False, LOXODROME=True) {[deg],[deg],-}
    :param magneticDeclinationGrid: geographical grid of a magnetic declination on Earth [deg]
    :param mass_const: kind of mass canculation {mass_integrated=False, mass_constant=True}.
    :param m_iter: number of iterations for integration loop [-]
    :param reducedPower: reduction of Power during the climb {True/False}
    :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}.
    :type speedType: string.
    :type v_init: float.
    :type v_final: float.
    :type phase: string.
    :type control: structure.
            :type ROCDtarget: float.
            :type slopetarget: float.
            :type acctarget: float.
            :type ESFtarget: float.
    :type maxRating: float.
    :type Hp_init: float.
    :type m_init: float.
    :type DeltaTemp: float.
    :type wS: float.
    :type turnMetrics: {float,float,string}.
    :type SOC_init: float.
    :type config: string.
    :type speedBrakes: dict{boolean,float}.
    :type speed_step: float.
    :type Lat: float.
    :type Lon: float.
    :type initialHeading: {float,float,boolean}.
    :type magneticDeclinationGrid: magneticDeclinationGrid.
    :type mass_const: boolean.
    :type m_iter: integer.
    :type reducedPower: boolean.
    :returns:
            BADA3: [Hp, TAS, CAS, M, GS, acc, ROCD, ESF, FUEL, FUELCONSUMED, THR, time, dist, slope, mass, config, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, N, s, NM, deg, kg, -, -,deg,deg,deg,deg,deg,deg/s]
            BADA4: [Hp, TAS, CAS, M, GS, acc, ROCD, ESF, FUEL, FUELCONSUMED, THR, P[Pmec, Pbat, Pelc Ibat, Vbat, Vgbat, SOCr, SOC], time, dist, slope, mass, config, HLid, LG, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, N, [W,W,W,A,V,V,%/h,%], s, NM, deg, kg, -, -, -, -,deg,deg,deg,deg,deg,deg/s]
            BADAH: [Hp, TAS, CAS, M, GS, ROCD, ESF, FUEL, FUELCONSUMED, Peng, Preq, Pav, time, dist, slope, mass, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, W, W, W, s, NM, deg, kg, -,deg,deg,deg,deg,deg,deg/s]
            BADAE: [time, dist, Hp, TAS, CAS, M, GS, acc, ROCD, ESF, slope, mass, P[Pmec, Pelc, Pbat, SOCr, SOC, Ibat, Vbat, Vgbat] comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [s, NM, ft, kt, kt, -, kt, m/s^2, ft/min, deg, kg, [W,W,W,%/h,%,A,V,V], -,deg,deg,deg,deg,deg,deg/s]
    :rtype: dict[list[float]}.
    """

    rateOfTurn = turnMetrics["rateOfTurn"]
    bankAngle = turnMetrics["bankAngle"]
    directionOfTurn = turnMetrics["directionOfTurn"]

    turnFlight = True
    if turnMetrics["rateOfTurn"] == 0.0 and turnMetrics["bankAngle"] == 0.0:
        turnFlight = False

    # conversion of Magnetic Heading to True Heading
    if magneticDeclinationGrid is not None:
        magneticDeclination = magneticDeclinationGrid.getMagneticDeclination(
            LAT_target=Lat, LON_target=Lon
        )
    else:
        magneticDeclination = 0

    # retrieve magnetic and true heading inputs
    magneticHeading = initialHeading["magnetic"]
    trueHeading = initialHeading["true"]
    constantHeading = initialHeading["constantHeading"]

    if Lat and Lon and (magneticHeading or trueHeading):
        if trueHeading is not None and magneticHeading is None:
            # fly TRUE Heading
            headingToFly = "TRUE"
            magneticHeading = trueHeading - magneticDeclination
        elif magneticHeading is not None and trueHeading is None:
            # fly MAGNETIC Heading
            if constantHeading == True:
                headingToFly = "MAGNETIC"
                trueHeading = magneticHeading + magneticDeclination
            else:
                raise Exception("Cannot fly non-constant magnetic heading")

        else:
            raise Exception("Undefined Heading value combination")

    # calculation with constant mass (True) or integrated (False)
    mass_const = kwargs.get("mass_const", False)

    # optional parameter to define initial Baterry State of Charge (SOC)
    if AC.BADAFamily.BADAE:
        SOC_init = kwargs.get("SOC_init", 100)
    else:
        SOC_init = None

    # speed brakes application
    if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
        speedBrakes = kwargs.get("speedBrakes", {"deployed": False, "value": 0.03})

    # iteratin step of speed loop
    if speedType == "M":
        speed_step = kwargs.get("speed_step", 0.01)  # [-] Mach increment
    elif speedType == "CAS" or speedType == "TAS":
        speed_step = kwargs.get("speed_step", 5.0)  # [kt] CAS/TAS increment

    # number of iteration of mass/altitude loop
    # BADAE
    if AC.BADAFamily.BADAE:
        m_iter = kwargs.get("m_iter", 5)  # number of iterations for integration loop[-]
    # BADA3 or BADA4 or BADAH
    else:
        m_iter = kwargs.get(
            "m_iter", 10
        )  # number of iterations for integration loop[-]

    # Determine if speed evolution over the segment is acceleration or deceleration
    # and associated speed iteration direction
    if v_init < v_final:
        speedEvol = "acc"
    else:
        speedEvol = "dec"
        speed_step = -speed_step

    if control is None:
        # create empty control target
        control = target()

    # check the consistency of SLOPE/ROCD and climb/descent phase of flight
    # if incosistent, change the sign on slope/ROCD target value
    if phase == "Climb":
        if control.slopetarget is not None and control.slopetarget < 0:
            control.slopetarget = abs(control.slopetarget)
            print("Slopetarget for Climb should be positive")
        if control.ROCDtarget is not None and control.ROCDtarget < 0:
            control.ROCDtarget = abs(control.ROCDtarget)
            print("ROCDtarget for Climb should be positive")
    elif phase == "Descent":
        if control.slopetarget is not None and control.slopetarget > 0:
            control.slopetarget = control.slopetarget * (-1)
            print("Slopetarget for Descent should be negative")
        if control.ROCDtarget is not None and control.ROCDtarget > 0:
            control.ROCDtarget = control.ROCDtarget * (-1)
            print("ROCDtarget for Descent should be negative")

    # check the consistency of acc/dec and ESF
    if phase == "Cruise":
        if control.ESFtarget is not None and control.ESFtarget != 0:
            control.ESFtarget = 0
    elif phase == "Climb":
        if (
            control.ESFtarget is not None
            and speedEvol == "acc"
            and control.ESFtarget > 1
        ):
            print("ESFtarget for acceleration in Climb should be < 1")
        if (
            control.ESFtarget is not None
            and speedEvol == "dec"
            and control.ESFtarget < 1
        ):
            print("ESFtarget for deceleration in Climb should be > 1")
    elif phase == "Descent":
        if (
            control.ESFtarget is not None
            and speedEvol == "acc"
            and control.ESFtarget < 1
        ):
            print("ESFtarget for acceleration in Descent should be > 1")
        if (
            control.ESFtarget is not None
            and speedEvol == "dec"
            and control.ESFtarget > 1
        ):
            print("ESFtarget for deceleration in Descent should be < 1")

    # check the consistency of acctarget and acc/dec
    if speedEvol == "acc":
        if control.acctarget is not None and control.acctarget < 0:
            control.acctarget = abs(control.acctarget)
            print("Acctarget in acceleration should be > 1")
    elif speedEvol == "dec":
        if control.acctarget is not None and control.acctarget > 0:
            control.acctarget = control.acctarget * (-1)
            print("Acctarget in deceleration should be < 1")

    if control.ROCDtarget is not None and control.slopetarget is not None:
        print("Both ROCD and SLOPE target provided, priority given to SLOPE")

    # comment line describing type of trajectory calculation
    controlComment = ""
    if control is not None:
        if control.ROCDtarget is not None:
            controlComment += "_" + "ROCDtarget"
        if control.slopetarget is not None:
            controlComment += "_" + "slopetarget"
        if control.acctarget is not None:
            controlComment += "_" + "acctarget"
        if control.ESFtarget is not None:
            controlComment += "_" + "ESFtarget"

    if turnFlight:
        turnComment = "_turn"
    else:
        turnComment = ""

    if constantHeading:
        constHeadingStr = "_const_Heading"
    elif constantHeading is False or constantHeading is None:
        constHeadingStr = ""

    # comment line describing type of trajectory calculation
    comment = (
        phase
        + turnComment
        + "_"
        + speedEvol
        + "_"
        + speedType
        + controlComment
        + constHeadingStr
    )

    if Lat and Lon and (magneticHeading or trueHeading):
        comment = comment + "_" + headingToFly + "_Heading"

    # compute Energy Share Factor
    if control.ESFtarget is not None:
        ESFc = control.ESFtarget
    elif control.ROCDtarget is not None or control.slopetarget is not None:
        ESFc = None
    elif control.acctarget is not None:
        ESFc = None

        # update ROCD target if phase is set to "Cruise"
        if phase == "Cruise":
            control.ROCDtarget = 0
    else:
        # Neither ROCD/slope nor ESF/acc provided means control is rating+default ESF
        if (phase == "Climb" and speedEvol == "acc") or (
            phase == "Descent" and speedEvol == "dec"
        ):
            ESFc = 0.3
        elif (phase == "Climb" and speedEvol == "dec") or (
            phase == "Descent" and speedEvol == "acc"
        ):
            ESFc = 1.7
        elif phase == "Cruise":
            ESFc = 0

    # convert target values to ISU units
    if control.ROCDtarget is not None:
        ROCDtargetisu = conv.ft2m(control.ROCDtarget) / 60  # [m/s]
    if control.slopetarget is not None:
        slopetargetisu = conv.deg2rad(control.slopetarget)

    # Determine max engine rating if not provided
    if "maxRating" not in kwargs:
        if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
            maxRating = "MTKF"
        else:
            maxRating = "MCMB"
    else:
        maxRating = checkArgument("maxRating", **kwargs)

    # Determine engine rating
    if (control.ROCDtarget is not None or control.slopetarget is not None) and (
        control.ESFtarget is not None or control.acctarget is not None
    ):
        rating = None
    else:
        if phase == "Climb" or (phase == "Cruise" and speedEvol == "acc"):
            rating = maxRating
        elif phase == "Descent" or (phase == "Cruise" and speedEvol == "dec"):
            if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
                rating = "UNKNOWN"  # TBD: No minimum power model
            else:
                rating = "LIDL"

    # get the default aerodynamic configuration if provided to be used for the whole segment
    config_default = kwargs.get("config", None)
    if config_default is not None:
        if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
            if not (
                config_default == "TO"
                or config_default == "IC"
                or config_default == "CR"
                or config_default == "AP"
                or config_default == "LD"
            ):
                print(
                    "WRONG default configuration set. Available values are: TO/IC/CR/AP/LD. Configuration will be calculated automatically"
                )

    # initialize output parameters
    Hp = [Hp_init]
    TAS = []
    CAS = []
    GS = []
    M = []
    ROCD = []
    esf = []
    FUEL = []
    FUELCONSUMED = []
    time = [0]
    dist = [0]
    mass = [m_init]
    ESF = []
    Comment = []
    check = []  # TEM consistency check result
    Slope = []
    acc = []
    THR = []
    DRAG = []
    config = []
    HLid = []
    LG = []
    BankAngle = []
    ROT = []

    if not AC.BADAFamily.BADAE:
        FUELCONSUMED = [0]

    # BADAH specific
    Peng = []
    Preq = []
    Pav = []

    # optional GPS coordiantes and HDG definition
    if Lat and Lon and (magneticHeading or trueHeading):
        LAT = [Lat]
        LON = [Lon]
        HDGMagnetic = [magneticHeading]
        HDGTrue = [trueHeading]
    else:
        LAT = []
        LON = []
        HDGMagnetic = []
        HDGTrue = []

    # BADAE specific
    Pmec = []
    Pbat = []
    SOCr = []
    SOC = [SOC_init]
    Pelc = []
    Ibat = []
    Vbat = []
    Vgbat = []

    # init loop parameters
    dVdtisu = []

    # initialize loop parameters: speed at end of step and loop termination
    v_i = v_init
    go_on = True

    while go_on:
        ## PART 1: compute parameters at end of step that are known without uncertainties:
        ##         *none*

        ## PART 2: compute parameters at end of step that are known only with uncertainties
        ##           (due to unknown mass and altitude at end of step):

        # Initialize loop parameters: aircraft mass (resp. altitude) at end of step is
        # first estimated as equal to aircraft mass (resp. altitude) at start of step

        mass_i = mass[-1]
        Hp_i = Hp[-1]
        for _ in itertools.repeat(None, m_iter):
            # atmosphere properties
            H_m = conv.ft2m(Hp_i)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
            temp_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)

            # aircraft speed
            [M_i, CAS_i, TAS_i] = atm.convertSpeed(
                v=v_i, speedType=speedType, theta=theta, delta=delta, sigma=sigma
            )

            if turnFlight:
                if turnMetrics["bankAngle"] != 0.0:
                    # bankAngle is defined
                    rateOfTurn = AC.rateOfTurn_bankAngle(TAS=TAS_i, bankAngle=bankAngle)
                else:
                    # rateOfTurn is defined
                    bankAngle = AC.bankAngle(
                        rateOfTurn=rateOfTurn, v=TAS_i
                    )  # [degrees]

            # Load factor
            nz = 1 / cos(radians(bankAngle))

            # compute ROCD target (if any) on this step
            if control.slopetarget is not None:
                # compute target ROCD corresponding to target slope

                if AC.BADAFamily.BADAE:
                    # special case for BADAE, in future it may apply also for BADAH
                    dh_dt_i = TAS_i * tan(slopetargetisu)
                else:
                    dh_dt_i = TAS_i * sin(slopetargetisu)

                ROCDtargetisu_i = dh_dt_i * (1 / temp_const)
            elif control.ROCDtarget is not None:
                ROCDtargetisu_i = ROCDtargetisu
                dh_dt_i = ROCDtargetisu_i * temp_const

            # BADAH or BADAE
            if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
                # compute Power required
                Preq_i = AC.Preq(sigma=sigma, tas=TAS_i, mass=mass_i, phi=bankAngle)

                # compute engine power
                if rating is None:
                    # compute power required for the manoeuver
                    if ESFc is not None:
                        P_i = dh_dt_i * mass_i * const.g / ESFc + Preq_i  # [W]
                    elif control.acctarget is not None:
                        P_i = (
                            dh_dt_i * mass_i * const.g
                            + mass_i * TAS_i * control.acctarget
                            + Preq_i
                        )  # [W]
                    else:
                        print("Error: neither ESF nor acc target provided")

                    # Check that required thrust/power fits in the available thrust/power envelope,
                    # recompute ROCD if necessary and compute fuel coefficient accordingly
                    Pmin = 0.1 * AC.P0  # No minimum power model: assume 10% torque

                    if AC.BADAFamily.BADAH:
                        Pmax = AC.Pav(rating=maxRating, theta=theta, delta=delta)
                        Pav_i = AC.Pav(rating=maxRating, theta=theta, delta=delta)
                    elif AC.BADAFamily.BADAE:
                        Pmax = AC.Pav(rating=maxRating, SOC=SOC[-1])
                        Pav_i = AC.Pav(rating=maxRating, SOC=SOC[-1])

                    if P_i < Pmin:
                        P_i = Pmin

                        if ESFc is not None:
                            ROCD_i = (
                                conv.m2ft(
                                    AC.ROCD(
                                        Peng=P_i,
                                        Preq=Preq_i,
                                        mass=mass_i,
                                        ESF=ESFc,
                                        theta=theta,
                                        DeltaTemp=DeltaTemp,
                                    )
                                )
                                * 60
                            )
                        elif control.acctarget is not None:
                            ROCD_i = (
                                conv.m2ft(
                                    (P_i - mass_i * TAS_i * control.acctarget - Preq_i)
                                    / (mass_i * const.g * temp_const)
                                )
                                * 60
                            )

                    elif P_i > Pmax:
                        P_i = Pmax

                        if ESFc is not None:
                            ROCD_i = (
                                conv.m2ft(
                                    AC.ROCD(
                                        Peng=P_i,
                                        Preq=Preq_i,
                                        mass=mass_i,
                                        ESF=ESFc,
                                        theta=theta,
                                        DeltaTemp=DeltaTemp,
                                    )
                                )
                                * 60
                            )
                        elif control.acctarget is not None:
                            ROCD_i = (
                                conv.m2ft(
                                    (P_i - mass_i * TAS_i * control.acctarget - Preq_i)
                                    / (mass_i * const.g * temp_const)
                                )
                                * 60
                            )
                    else:
                        ROCD_i = control.ROCDtarget

                else:
                    # Compute available power
                    if rating == "UNKNOWN":
                        P_i = 0.1 * AC.P0  # No minimum power model: assume 10% torque
                        Pav_i = AC.Pav(rating=maxRating, theta=theta, delta=delta)
                    else:
                        if AC.BADAFamily.BADAH:
                            P_i = AC.Pav(rating=rating, theta=theta, delta=delta)
                            Pav_i = AC.Pav(rating=rating, theta=theta, delta=delta)
                        elif AC.BADAFamily.BADAE:
                            P_i = AC.Pav(rating=rating, SOC=SOC[-1])
                            Pav_i = AC.Pav(rating=rating, SOC=SOC[-1])

                # Compute excess power
                Pe_i = P_i - Preq_i  # [W]

                # BADAH
                if AC.BADAFamily.BADAH:
                    # compute fuel flow for level flight
                    CP = AC.CP(Peng=P_i)
                    FUEL_i = AC.ff(delta=delta, CP=CP)  # [kg/s]

                # BADAE
                elif AC.BADAFamily.BADAE:
                    Pbat_i = AC.Pbat(Preq=P_i, SOC=SOC[-1])
                    SOCr_i = AC.SOCrate(Preq=P_i, SOC=SOC[-1])

                    # debug data
                    Pelc_i = P_i / AC.eta
                    Ibat_i = AC.Ibat(P=Pelc_i, SOC=SOC[-1])
                    Vbat_i = AC.Vbat(I=Ibat_i, SOC=SOC[-1])
                    Vgbat_i = AC.Vocbat(SOC=SOC[-1]) - AC.R0bat(SOC=SOC[-1]) * Ibat_i

            # BADA4
            elif AC.BADAFamily.BADA4:
                # aircraft aerodynamic configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=phase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=phase, previousConfig=config[-1], currentConfig=config_i
                    )

                [HLid_i, LG_i] = AC.flightEnvelope.getAeroConfig(config=config_i)

                # compute lift coefficient
                CL = AC.CL(M=M_i, delta=delta, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(M=M_i, CL=CL, HLid=HLid_i, LG=LG_i, speedBrakes=speedBrakes)
                # compute drag force
                Drag = AC.D(M=M_i, delta=delta, CD=CD)

                # compute thrust and fuel flow
                if rating is None:
                    # compute thrust force required for the manoeuver
                    if ESFc is not None:
                        THR_i = (
                            dh_dt_i * mass_i * const.g / (TAS_i * ESFc) + Drag
                        )  # [N]
                    elif control.acctarget is not None:
                        THR_i = (
                            dh_dt_i * mass_i * const.g / TAS_i
                            + mass_i * control.acctarget
                            + Drag
                        )  # [N]
                    else:
                        print("Error: neither ESF nor acc target provided")

                    # Check that required thrust fits in the available thrust envelope,
                    # recompute ROCD if necessary and compute fuel flow accordingly
                    THR_min = AC.Thrust(
                        rating="LIDL",
                        delta=delta,
                        theta=theta,
                        M=M_i,
                        DeltaTemp=DeltaTemp,
                    )  # IDLE Thrust
                    FUEL_min = AC.ff(
                        rating="LIDL",
                        delta=delta,
                        theta=theta,
                        M=M_i,
                        DeltaTemp=DeltaTemp,
                    )  # IDLE Fuel Flow
                    THR_max = AC.Thrust(
                        rating="MCMB",
                        delta=delta,
                        theta=theta,
                        M=M_i,
                        DeltaTemp=DeltaTemp,
                    )  # MCMB Thrust
                    FUEL_max = AC.ff(
                        rating="MCMB",
                        delta=delta,
                        theta=theta,
                        M=M_i,
                        DeltaTemp=DeltaTemp,
                    )  # MCMB Fuel Flow

                    if THR_i < THR_min:
                        THR_i = THR_min
                        FUEL_i = FUEL_min
                    elif THR_i > THR_max:
                        THR_i = THR_max
                        FUEL_i = FUEL_max
                    else:
                        CT = AC.CT(Thrust=THR_i, delta=delta)
                        FUEL_i = AC.ff(
                            CT=CT, delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                        )
                else:
                    THR_i = AC.Thrust(
                        rating=rating,
                        delta=delta,
                        theta=theta,
                        M=M_i,
                        DeltaTemp=DeltaTemp,
                    )  # [N]
                    CT = AC.CT(Thrust=THR_i, delta=delta)
                    FUEL_i = AC.ff(
                        CT=CT, delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                    )

                # compute excess power
                Pe_i = (THR_i - Drag) * TAS_i  # [kg*m^2/s^3]

            # BADA3
            elif AC.BADAFamily.BADA3:
                # aircraft aerodynamic configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=phase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=phase,
                        previousConfig=config[-1],
                        currentConfig=config_i,
                    )

                # compute lift coefficient
                CL = AC.CL(tas=TAS_i, sigma=sigma, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(CL=CL, config=config_i, speedBrakes=speedBrakes)
                # compute drag force
                Drag = AC.D(tas=TAS_i, sigma=sigma, CD=CD)

                # compute thrust and fuel flow
                if rating is None:
                    # compute thrust force required for the manoeuver
                    if ESFc is not None:
                        THR_i = (
                            dh_dt_i * mass_i * const.g / (TAS_i * ESFc) + Drag
                        )  # [N]
                    elif control.acctarget is not None:
                        THR_i = (
                            dh_dt_i * mass_i * const.g / TAS_i
                            + mass_i * control.acctarget
                            + Drag
                        )  # [N]
                    else:
                        print("Error: neither ESF nor acc target provided")

                    # Check that required thrust fits in the available thrust envelope,
                    # recompute ROCD if necessary and compute fuel flow accordingly

                    THR_min = AC.Thrust(
                        rating="LIDL", v=TAS_i, h=H_m, config="CR", DeltaTemp=DeltaTemp
                    )  # IDLE Thrust
                    FUEL_min = AC.ff(
                        flightPhase="Descent",
                        v=TAS_i,
                        h=H_m,
                        T=THR_min,
                        config="CR",
                        adapted=False,
                    )  # IDLE Fuel Flow
                    THR_max = AC.Thrust(
                        rating="MCMB", v=TAS_i, h=H_m, config="CR", DeltaTemp=DeltaTemp
                    )  # MCMB Thrust
                    FUEL_max = AC.ff(
                        flightPhase="Climb",
                        v=TAS_i,
                        h=H_m,
                        T=THR_max,
                        config="CR",
                        adapted=False,
                    )  # MCMB Fuel Flow

                    if THR_i < THR_min:
                        THR_i = THR_min
                        FUEL_i = FUEL_min
                    elif THR_i > THR_max:
                        THR_i = THR_max
                        FUEL_i = FUEL_max
                    else:
                        FUEL_i = AC.ff(
                            v=TAS_i, h=H_m, T=THR_i, config=config_i, adapted=True
                        )
                else:
                    THR_i = AC.Thrust(
                        rating=rating,
                        v=TAS_i,
                        h=H_m,
                        config=config_i,
                        DeltaTemp=DeltaTemp,
                    )
                    if rating == "MCMB" or rating == "MTKF":
                        FUEL_i = AC.ff(
                            flightPhase="Climb",
                            v=TAS_i,
                            h=H_m,
                            T=THR_i,
                            config=config_i,
                        )
                    elif rating == "MCRZ":
                        FUEL_i = AC.ff(
                            flightPhase="Cruise",
                            v=TAS_i,
                            h=H_m,
                            T=THR_i,
                            config=config_i,
                        )
                    elif rating == "LIDL":
                        FUEL_i = AC.ff(
                            flightPhase="Descent",
                            v=TAS_i,
                            h=H_m,
                            T=THR_i,
                            config=config_i,
                        )

                # compute excess power
                Pe_i = (THR_i - Drag) * TAS_i  # [kg*m^2/s^3]

            if ESFc is not None:
                ESF_i = ESFc
                # compute power dedicated to climb
                PC_i = Pe_i * ESF_i  # [kg*m^2/s^3]
                # compute ROCD
                dhdtisu = PC_i / (mass_i * const.g)  # [m/s]
                ROCDisu = dhdtisu * 1 / temp_const  # [m/s]
                ROCD_i = conv.m2ft(ROCDisu) * 60  # [ft/min]
            elif control.acctarget is not None:
                # compute power required for acc/dec rate
                Pa_i = mass_i * TAS_i * control.acctarget
                # check that required power fits in the available power envelope
                if abs(Pa_i) > abs(Pe_i):
                    Pa_i = Pe_i
                # compute power dedicated to climb
                PC_i = Pe_i - Pa_i  # [kg*m^2/s^3]

                if Pe_i != 0:
                    ESF_i = PC_i / Pe_i
                else:
                    # ESF_i = float("Inf")
                    ESF_i = float(0)

                # compute ROCD
                dhdtisu = PC_i / (mass_i * const.g)  # [m/s]
                ROCDisu = dhdtisu * 1 / temp_const  # [m/s]
                ROCD_i = conv.m2ft(ROCDisu) * 60  # [ft/min]
            elif control.slopetarget is not None or control.ROCDtarget is not None:
                dhdtisu = dh_dt_i  # [m/s]
                ROCDisu = dh_dt_i * 1 / temp_const  # [m/s]
                ROCD_i = conv.m2ft(ROCDisu) * 60  # [ft/min]
                PC_i = dh_dt_i * (mass_i * const.g)  # [kg*m^2/s^3]

                if Pe_i != 0:
                    ESF_i = PC_i / Pe_i
                else:
                    ESF_i = float("Inf")
            else:
                print("Error: unexpected combination of targets")

            # compute acceleration
            if TAS_i == 0:
                dVdtisu_i = (Pe_i - PC_i) / (mass_i * (TAS_i + 0.5))  # [m/s^2]
            else:
                dVdtisu_i = (Pe_i - PC_i) / (mass_i * TAS_i)  # [m/s^2]

            # Compute elapsed time, altitude and fuel burn over current step
            if v_i == v_init:
                # no need to loop for first point: initial m/Hp already known
                break
            else:
                # Average acceleration over step is the mean of initial and final ones
                step_dVdtisu = (dVdtisu[-1] + dVdtisu_i) / 2  # [m/s^2]
                # Step time is: TAS differential divided by average acceleration
                step_time = (TAS_i - conv.kt2ms(TAS[-1])) / step_dVdtisu
                # Average ROCD over step is the mean of initial and final ones
                step_ROCD = (ROCD[-1] + ROCD_i) / 2  # [ft/min]
                # Altitude differential is: average ROCD multiplied by step time
                step_Hp = step_ROCD * step_time / 60  # [ft]
                # Update altitude estimate at end of step
                Hp_i = Hp[-1] + step_Hp  # [ft]

                # BADAE
                if AC.BADAFamily.BADAE:
                    # Average SOC rate over step is the mean of initial and final ones
                    step_SOCr = (SOCr[-1] + SOCr_i) / 2  # [%/h]
                    # SOC change is: average SOC rate multiplied by step time
                    step_SOC = step_SOCr * step_time / 3600  # [%]
                    # Update SOC estimate at end of step
                    SOC_i = SOC[-1] - step_SOC  # [%]
                    # update of aircraft mass estimate at end of step - mass is not changing for ELECTRIC engine (no fuel is consumed)
                    mass_i = mass[-1]  # [kg]

                else:
                    # Average fuel flow over step is the mean of initial and final ones
                    step_FUEL = (FUEL[-1] + FUEL_i) / 2  # [kg/s]
                    # Fuel burnt is: average fuel flow multiplied by step time
                    step_mass = step_FUEL * step_time  # [kg]
                    # Update aircraft mass estimate at end of step
                    if not mass_const:
                        mass_i = mass[-1] - step_mass  # [kg]
                        fuelConsumed_i = step_FUEL * step_time
                    fuelConsumed_i = FUELCONSUMED[-1] + step_FUEL * step_time

        ## PART 3: store information about end of step point

        # point data
        TAS.append(conv.ms2kt(TAS_i))
        CAS.append(conv.ms2kt(CAS_i))
        M.append(M_i)
        dVdtisu.append(dVdtisu_i)
        ROCD.append(ROCD_i)
        esf.append(ESF_i)
        Comment.append(comment)

        # everything except electric BADAE
        if not AC.BADAFamily.BADAE:
            FUEL.append(FUEL_i)

        # BADAH
        if AC.BADAFamily.BADAH:
            Peng.append(P_i)
            Preq.append(Preq_i)
            Pav.append(Pav_i)

        # BADAE
        elif AC.BADAFamily.BADAE:
            Pmec.append(P_i)
            Pbat.append(Pbat_i)
            SOCr.append(SOCr_i)
            Pelc.append(Pelc_i)
            Ibat.append(Ibat_i)
            Vbat.append(Vbat_i)
            Vgbat.append(Vgbat_i)

        # BADA3 & BADA4
        elif AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
            THR.append(THR_i)
            DRAG.append(Drag)
            config.append(config_i)

        # BADA4
        if AC.BADAFamily.BADA4:
            HLid.append(HLid_i)
            LG.append(LG_i)

        # TEM consistency check
        # BADAH or BADAE
        if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
            check.append(
                P_i - Preq_i - mass_i * const.g * dhdtisu - mass_i * TAS_i * dVdtisu_i
            )

        # BADA3 or BADA4
        elif AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
            check.append(
                (THR_i - Drag) * TAS_i
                - mass_i * const.g * dhdtisu
                - mass_i * TAS_i * dVdtisu_i
            )

        # calculation of the slope
        if TAS_i == 0:
            gamma_i = 90 * np.sign(ROCD_i)
        else:
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=conv.ft2m(Hp_i), DeltaTemp=DeltaTemp
            )
            temp_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)
            if AC.BADAFamily.BADAE:
                gamma_i = degrees(atan(conv.ft2m(ROCD_i) * temp_const / 60 / TAS_i))
            else:
                # using SIN assumes the TAS to be in the direction of the aircraft axis, not ground plane. Which means, this should be mathematically the correct equation for all the aircraft
                gamma_i = degrees(asin(conv.ft2m(ROCD_i) * temp_const / 60 / TAS_i))

        # ground speed can be calcualted as TAS projected on the ground minus wind
        GS_i = cos(radians(gamma_i)) * TAS_i - wS
        GS.append(conv.ms2kt(GS_i))

        Slope.append(gamma_i)
        BankAngle.append(bankAngle)
        ROT.append(rateOfTurn)

        # integrated data
        if v_i != v_init:  # exclude first point: initial t/d/m already known
            if AC.BADAFamily.BADAE:
                SOC.append(SOC_i)

            # everything except electric BADAE
            if not AC.BADAFamily.BADAE:
                FUELCONSUMED.append(fuelConsumed_i)

            # Altitude at end of step has been termined in PART 2
            Hp.append(Hp_i)
            # Aircraft mass at end of step has been termined in PART 2
            mass.append(mass_i)
            # Time at end of step is time at start of step plus step time
            time.append(time[-1] + step_time)

            # Average TAS over step is the mean of initial and final ones
            step_TAS = (TAS[-2] + TAS[-1]) / 2  # [kt]
            # Average slope over the step
            step_gamma = radians((Slope[-2] + Slope[-1]) / 2)  # radians
            # Average ground speed over step
            # since this is not level flight, TAS speed should be projected on the ground, then GS can be calculated applying the wind speed
            step_TAS_projected = cos(step_gamma) * step_TAS
            step_GS = step_TAS_projected - wS  # [kt]
            # Step distance is: average GS multiplied by step time
            if turnFlight:
                step_distance = conv.m2nm(
                    turn.distance(
                        rateOfTurn=rateOfTurn, TAS=TAS_i, timeOfTurn=step_time
                    )
                )  # arcLength during the turn [NM]
            else:
                step_distance = step_GS * step_time / 3600  # [NM]
            # Distance at end of step is distance at start of step plus step distance
            dist.append(dist[-1] + step_distance)

            # add GPS calculation
            if Lat and Lon and (magneticHeading or trueHeading):
                if headingToFly == "TRUE":
                    if not turnFlight:
                        if not constantHeading:
                            # fly ORTHODROME
                            (Lat_i, Lon_i, HDGTrue_i) = (
                                vincenty.destinationPoint_finalBearing(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    distance=conv.nm2m(step_distance),
                                    bearing=HDGTrue[-1],
                                )
                            )

                            if magneticDeclinationGrid is not None:
                                HDGMagnetic_i = (
                                    HDGTrue_i
                                    - magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGMagnetic_i = HDGTrue_i

                        elif constantHeading:
                            # fly LOXODROME
                            (Lat_i, Lon_i) = rhumb.destinationPoint(
                                LAT_init=LAT[-1],
                                LON_init=LON[-1],
                                bearing=HDGTrue[-1],
                                distance=conv.nm2m(step_distance),
                            )
                            HDGTrue_i = HDGTrue[-1]

                            if magneticDeclinationGrid is not None:
                                HDGMagnetic_i = (
                                    HDGTrue_i
                                    - magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGMagnetic_i = HDGTrue_i

                    else:
                        # calculate the turn
                        (Lat_i, Lon_i, HDGTrue_i) = turn.destinationPoint_finalBearing(
                            LAT_init=LAT[-1],
                            LON_init=LON[-1],
                            bearingInit=HDGTrue[-1],
                            TAS=TAS_i,
                            rateOfTurn=rateOfTurn,
                            timeOfTurn=step_time,
                            directionOfTurn=directionOfTurn,
                        )

                        if magneticDeclinationGrid is not None:
                            HDGMagnetic_i = (
                                HDGTrue_i
                                - magneticDeclinationGrid.getMagneticDeclination(
                                    LAT_target=Lat_i, LON_target=Lon_i
                                )
                            )
                        else:
                            magneticDeclination = 0
                            HDGMagnetic_i = HDGTrue_i

                elif headingToFly == "MAGNETIC":
                    if not turnFlight:
                        if constantHeading:
                            (Lat_i, Lon_i) = rhumb.destinationPoint(
                                LAT_init=LAT[-1],
                                LON_init=LON[-1],
                                bearing=HDGTrue[-1],
                                distance=conv.nm2m(step_distance),
                            )
                            HDGMagnetic_i = HDGMagnetic[-1]

                            if magneticDeclinationGrid is not None:
                                HDGTrue_i = (
                                    HDGMagnetic_i
                                    + magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGTrue_i = HDGMagnetic_i

                    else:
                        # calculate the turn
                        (Lat_i, Lon_i, HDGTrue_i) = turn.destinationPoint_finalBearing(
                            LAT_init=LAT[-1],
                            LON_init=LON[-1],
                            bearingInit=HDGTrue[-1],
                            TAS=TAS_i,
                            rateOfTurn=rateOfTurn,
                            timeOfTurn=step_time,
                            directionOfTurn=directionOfTurn,
                        )

                        if magneticDeclinationGrid is not None:
                            HDGMagnetic_i = (
                                HDGTrue_i
                                - magneticDeclinationGrid.getMagneticDeclination(
                                    LAT_target=Lat_i, LON_target=Lon_i
                                )
                            )
                        else:
                            magneticDeclination = 0
                            HDGMagnetic_i = HDGTrue_i

                LAT.append(Lat_i)
                LON.append(Lon_i)
                HDGMagnetic.append(HDGMagnetic_i)
                HDGTrue.append(HDGTrue_i)

        # Determine end speed of next step
        v_next = v_i + speed_step

        if speedEvol == "acc":
            if v_next < v_final:
                v_i = v_next
            elif v_i < v_final:
                v_i = v_final
            else:
                go_on = False
        else:
            if v_next > v_final:
                v_i = v_next
            elif v_i > v_final:
                v_i = v_final
            else:
                go_on = False

    flightData = {
        "Hp": Hp,
        "TAS": TAS,
        "CAS": CAS,
        "GS": GS,
        "M": M,
        "acc": dVdtisu,
        "ROCD": ROCD,
        "ESF": esf,
        "FUEL": FUEL,
        "Pmec": Pmec,
        "Pelc": Pelc,
        "Pbat": Pbat,
        "SOCr": SOCr,
        "SOC": SOC,
        "Ibat": Ibat,
        "Vbat": Vbat,
        "Vgbat": Vgbat,
        "FUELCONSUMED": FUELCONSUMED,
        "Preq": Preq,
        "Peng": Peng,
        "Pav": Pav,
        "THR": THR,
        "DRAG": DRAG,
        "time": time,
        "dist": dist,
        "slope": Slope,
        "mass": mass,
        "config": config,
        "HLid": HLid,
        "LG": LG,
        "LAT": LAT,
        "LON": LON,
        "HDGTrue": HDGTrue,
        "HDGMagnetic": HDGMagnetic,
        "BankAngle": BankAngle,
        "ROT": ROT,
        "comment": Comment,
    }

    flightTrajectory = FT.createFlightTrajectoryDataframe(flightData)
    return flightTrajectory


def accDec_time(
    AC,
    length,
    speedType,
    v_init,
    speedEvol,
    phase,
    Hp_init,
    m_init,
    DeltaTemp,
    wS=0.0,
    turnMetrics={"rateOfTurn": 0.0, "bankAngle": 0.0, "directionOfTurn": None},
    control=None,
    Lat=None,
    Lon=None,
    initialHeading={"magnetic": None, "true": None, "constantHeading": None},
    reducedPower=None,
    magneticDeclinationGrid=None,
    **kwargs,
):
    """This function computes time and fuel required by an aircraft to perform an acceleration/deceleration from v_init for set amount of time in climb cruise or descent

    .. note::
            The control law used during the segment depends on the targets provided in input parameter 'control':
            - ROCD/slope+ESF:  law is ROCD/slope+ESF
            - ROCD/slope+acc:  law is ROCD/slope+acc
            - ROCD/slope only: law is rating+ROCD/slope
            - ESF only:        law is rating+ESF
            - acc only:        law is rating+acc
            - Neither:         law is rating+default ESF

    :param AC: aircraft {BADA3/4/H/E}
    :param speedType: what kind of speed is followed {M, CAS, TAS}.
    :param length: length of a segment to fly [s]
    :param step_length: length of a step of a segment - [s]
    :param v_init: initial speed to follow - [kt] CAS/TAS speed to follow or [-] MACH speed to follow.
    :param speedEvol: speed evolution {acc, dec}
    :param phase: vertical evolution {Climb, Descent, Cruise}
    :param control: structure containing a combination of the following targets:
            :param ROCDtarget: Rate of climb/descent to be followed [ft/min].
            :param slopetarget: slope (flight path angle) to be followed [deg].
            :param acctarget: acceleration to be followed [m/s^2].
            :param ESFtarget: Energy Share Factor to be followed [-].
    :param maxRating: rating to be used as a limit on the maximum thrust/power [-].
    :param Hp_init: initial pressure altitude [ft].
    :param m_init: initial aircraft mass [kg].
    :param DeltaTemp: deviation with respect to ISA [K].
    :param wS: longitudinal wind speed (TAS) [kt].
    :param turnMetrics: Metrics for turn performance {"rateOfTurn":0.0,"bankAngle":0.0,"directionOfTurn":None} {[deg/s],[deg],[LEFT/RIGHT]}
    :param SOC_init: initial state of charge [%].
    :param config: aircraft default aerodynamic configuration {TO,IC,CR,AP,LD}.
    :param speedBrakes: deployed or not speedbrakes including value to be added to the drag coeffcient {deployed:False,value:0.03} {deployed:[True/False],value:[-]}.
    :param Lat: Geographical Latitude [deg]
    :param Lon: Geographical Longitude [deg]
    :param initialHeading: aircraft magnetic heading, true heading and definition of constant heading(ORTHODROME=False, LOXODROME=True) {[deg],[deg],-}
    :param magneticDeclinationGrid: geographical grid of a magnetic declination on Earth [deg]
    :param mass_const: kind of mass canculation {mass_integrated=False, mass_constant=True}.
    :param m_iter: number of iterations for integration loop [-]
    :param reducedPower: reduction of Power during the climb {True/False}
    :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}.
    :type speedType: string.
    :type length: float.
    :type step_length: float.
    :type v_init: float.
    :type speedEvol: string.
    :type phase: string.
    :type control: structure.
            :type ROCDtarget: float.
            :type slopetarget: float.
            :type acctarget: float.
            :type ESFtarget: float.
    :type maxRating: float.
    :type Hp_init: float.
    :type m_init: float.
    :type DeltaTemp: float.
    :type wS: float.
    :type turnMetrics: {float,float,string}.
    :type SOC_init: float.
    :type config: string.
    :type speedBrakes: dict{boolean,float}.
    :type Lat: float.
    :type Lon: float.
    :type initialHeading: {float,float,boolean}.
    :type magneticDeclinationGrid: magneticDeclinationGrid.
    :type mass_const: boolean.
    :type m_iter: integer.
    :type reducedPower: boolean.
    :returns:
            BADA3: [Hp, TAS, CAS, M, GS, acc, ROCD, ESF, FUEL, FUELCONSUMED, THR, time, dist, slope, mass, config, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, N, s, NM, deg, kg, -, -,deg,deg,deg,deg,deg,deg/s]
            BADA4: [Hp, TAS, CAS, M, GS, acc, ROCD, ESF, FUEL, FUELCONSUMED, THR, P[Pmec, Pbat, Pelc Ibat, Vbat, Vgbat, SOCr, SOC], time, dist, slope, mass, config, HLid, LG, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, N, [W,W,W,A,V,V,%/h,%], s, NM, deg, kg, -, -, -, -,deg,deg,deg,deg,deg,deg/s]
            BADAH: [Hp, TAS, CAS, M, GS, ROCD, ESF, FUEL, FUELCONSUMED, Peng, Preq, Pav, time, dist, slope, mass, comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [ft, kt, kt, -, kt, m/s^2, ft/min, kg/s, kg, W, W, W, s, NM, deg, kg, -,deg,deg,deg,deg,deg,deg/s]
            BADAE: [time, dist, Hp, TAS, CAS, M, GS, acc, ROCD, ESF, slope, mass, P[Pmec, Pelc, Pbat, SOCr, SOC, Ibat, Vbat, Vgbat] comment, LAT, LON, HDGMagnetic, HDGTrue, bank angle, ROT] [s, NM, ft, kt, kt, -, kt, m/s^2, ft/min, deg, kg, [W,W,W,%/h,%,A,V,V], -,deg,deg,deg,deg,deg,deg/s]
    :rtype: dict[list[float]}.
    """

    rateOfTurn = turnMetrics["rateOfTurn"]
    bankAngle = turnMetrics["bankAngle"]
    directionOfTurn = turnMetrics["directionOfTurn"]

    turnFlight = True
    if turnMetrics["rateOfTurn"] == 0.0 and turnMetrics["bankAngle"] == 0.0:
        turnFlight = False

    # conversion of Magnetic Heading to True Heading
    if magneticDeclinationGrid is not None:
        magneticDeclination = magneticDeclinationGrid.getMagneticDeclination(
            LAT_target=Lat, LON_target=Lon
        )
    else:
        magneticDeclination = 0

    # retrieve magnetic and true heading inputs
    magneticHeading = initialHeading["magnetic"]
    trueHeading = initialHeading["true"]
    constantHeading = initialHeading["constantHeading"]

    if Lat and Lon and (magneticHeading or trueHeading):
        if trueHeading is not None and magneticHeading is None:
            # fly TRUE Heading
            headingToFly = "TRUE"
            magneticHeading = trueHeading - magneticDeclination
        elif magneticHeading is not None and trueHeading is None:
            # fly MAGNETIC Heading
            if constantHeading == True:
                headingToFly = "MAGNETIC"
                trueHeading = magneticHeading + magneticDeclination
            else:
                raise Exception("Cannot fly non-constant magnetic heading")

        else:
            raise Exception("Undefined Heading value combination")

    # calculation with constant mass (True) or integrated (False)
    mass_const = kwargs.get("mass_const", False)

    # optional parameter to define initial Baterry State of Charge (SOC)
    if AC.BADAFamily.BADAE:
        SOC_init = kwargs.get("SOC_init", 100)
    else:
        SOC_init = None

    # speed brakes application
    if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
        speedBrakes = kwargs.get("speedBrakes", {"deployed": False, "value": 0.03})

    # step size in [s]
    step_length = kwargs.get("step_length", 1)

    # number of iteration of mass/altitude loop
    # BADAE
    if AC.BADAFamily.BADAE:
        m_iter = kwargs.get("m_iter", 5)  # number of iterations for integration loop[-]
    # BADA3 or BADA4 or BADAH
    else:
        m_iter = kwargs.get(
            "m_iter", 10
        )  # number of iterations for integration loop[-]

    if control is None:
        # create empty control target
        control = target()

    # check the consistency of SLOPE/ROCD and climb/descent phase of flight
    # if incosistent, change the sign on slope/ROCD target value
    if phase == "Climb":
        if control.slopetarget is not None and control.slopetarget < 0:
            control.slopetarget = abs(control.slopetarget)
            print("Slopetarget for Climb should be positive")
        if control.ROCDtarget is not None and control.ROCDtarget < 0:
            control.ROCDtarget = abs(control.ROCDtarget)
            print("ROCDtarget for Climb should be positive")
    elif phase == "Descent":
        if control.slopetarget is not None and control.slopetarget > 0:
            control.slopetarget = control.slopetarget * (-1)
            print("Slopetarget for Descent should be negative")
        if control.ROCDtarget is not None and control.ROCDtarget < 0:
            control.ROCDtarget = control.ROCDtarget * (-1)
            print("ROCDtarget for Descent should be negative")

    # check the consistency of acc/dec and ESF
    if phase == "Cruise":
        if control.ESFtarget is not None and control.ESFtarget != 0:
            control.ESFtarget = 0
    elif phase == "Climb":
        if (
            control.ESFtarget is not None
            and speedEvol == "acc"
            and control.ESFtarget > 1
        ):
            print("ESFtarget for acceleration in Climb should be < 1")
        if (
            control.ESFtarget is not None
            and speedEvol == "dec"
            and control.ESFtarget < 1
        ):
            print("ESFtarget for deceleration in Climb should be > 1")
    elif phase == "Descent":
        if (
            control.ESFtarget is not None
            and speedEvol == "acc"
            and control.ESFtarget < 1
        ):
            print("ESFtarget for acceleration in Descent should be > 1")
        if (
            control.ESFtarget is not None
            and speedEvol == "dec"
            and control.ESFtarget > 1
        ):
            print("ESFtarget for deceleration in Descent should be < 1")

    # check the consistency of acctarget and acc/dec
    if speedEvol == "acc":
        if control.acctarget is not None and control.acctarget < 0:
            control.acctarget = abs(control.acctarget)
            print("Acctarget in acceleration should be > 1")
    elif speedEvol == "dec":
        if control.acctarget is not None and control.acctarget > 0:
            control.acctarget = control.acctarget * (-1)
            print("Acctarget in deceleration should be < 1")

    if control.ROCDtarget is not None and control.slopetarget is not None:
        print("Both ROCD and SLOPE target provided, priority given to SLOPE")

    # comment line describing type of trajectory calculation
    controlComment = ""
    if control is not None:
        if control.ROCDtarget is not None:
            controlComment += "_" + "ROCDtarget"
        if control.slopetarget is not None:
            controlComment += "_" + "slopetarget"
        if control.acctarget is not None:
            controlComment += "_" + "acctarget"
        if control.ESFtarget is not None:
            controlComment += "_" + "ESFtarget"

    if turnFlight:
        turnComment = "_turn"
    else:
        turnComment = ""

    if constantHeading:
        constHeadingStr = "_const_Heading"
    elif constantHeading is False or constantHeading is None:
        constHeadingStr = ""

    # comment line describing type of trajectory calculation
    comment = (
        phase
        + turnComment
        + "_"
        + speedEvol
        + "_"
        + speedType
        + controlComment
        + "_"
        + constHeadingStr
    )

    if Lat and Lon and (magneticHeading or trueHeading):
        comment = comment + "_" + headingToFly + "_Heading"

    # compute Energy Share Factor
    if control.ESFtarget is not None:
        ESFc = control.ESFtarget
    elif control.ROCDtarget is not None or control.slopetarget is not None:
        ESFc = None
    elif control.acctarget is not None:
        ESFc = None

        # update ROCD target if phase is set to "Cruise"
        if phase == "Cruise":
            control.ROCDtarget = 0
    else:
        # Neither ROCD/slope nor ESF/acc provided means control is rating+default ESF
        if (phase == "Climb" and speedEvol == "acc") or (
            phase == "Descent" and speedEvol == "dec"
        ):
            ESFc = 0.3
        elif (phase == "Climb" and speedEvol == "dec") or (
            phase == "Descent" and speedEvol == "acc"
        ):
            ESFc = 1.7
        elif phase == "Cruise":
            ESFc = 0

    # convert target values to ISU units
    if control.ROCDtarget is not None:
        ROCDtargetisu = conv.ft2m(control.ROCDtarget) / 60  # [m/s]
    if control.slopetarget is not None:
        slopetargetisu = conv.deg2rad(control.slopetarget)

    # Determine max engine rating if not provided
    if "maxRating" not in kwargs:
        if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
            maxRating = "MTKF"
        else:
            maxRating = "MCMB"
    else:
        maxRating = checkArgument("maxRating", **kwargs)

    # Determine engine rating
    if (control.ROCDtarget is not None or control.slopetarget is not None) and (
        control.ESFtarget is not None or control.acctarget is not None
    ):
        rating = None
    else:
        if phase == "Climb" or (phase == "Cruise" and speedEvol == "acc"):
            rating = maxRating
        elif phase == "Descent" or (phase == "Cruise" and speedEvol == "dec"):
            if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
                rating = "UNKNOWN"  # TBD: No minimum power model
            else:
                rating = "LIDL"

    # get the default aerodynamic configuration if provided to be used for the whole segment
    config_default = kwargs.get("config", None)
    if config_default is not None:
        if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
            if not (
                config_default == "TO"
                or config_default == "IC"
                or config_default == "CR"
                or config_default == "AP"
                or config_default == "LD"
            ):
                print(
                    "WRONG default configuration set. Available values are: TO/IC/CR/AP/LD. Configuration will be calculated automatically"
                )

    # initialize output parameters
    [theta_init, delta_init, sigma_init] = atm.atmosphereProperties(
        h=conv.ft2m(Hp_init), DeltaTemp=DeltaTemp
    )
    [M_init, CAS_init, TAS_init] = atm.convertSpeed(
        v=v_init,
        speedType=speedType,
        theta=theta_init,
        delta=delta_init,
        sigma=sigma_init,
    )

    Hp = [Hp_init]
    TAS = [conv.ms2kt(TAS_init)]
    CAS = [conv.ms2kt(CAS_init)]
    GS = []
    M = [M_init]
    ROCD = []
    esf = []
    FUEL = []
    FUELCONSUMED = []
    time = [0]
    dist = [0]
    mass = [m_init]
    ESF = []
    Comment = []
    check = []  # TEM consistency check result
    Slope = []
    acc = []
    THR = []
    DRAG = []
    config = []
    HLid = []
    LG = []
    BankAngle = []
    ROT = []

    if not AC.BADAFamily.BADAE:
        FUELCONSUMED = [0]

    # BADAH specific
    Peng = []
    Preq = []
    Pav = []

    # optional GPS coordiantes and HDG definition
    if Lat and Lon and (magneticHeading or trueHeading):
        LAT = [Lat]
        LON = [Lon]
        HDGMagnetic = [magneticHeading]
        HDGTrue = [trueHeading]
    else:
        LAT = []
        LON = []
        HDGMagnetic = []
        HDGTrue = []

    # BADAE specific
    Pmec = []
    Pbat = []
    SOCr = []
    SOC = [SOC_init]
    Pelc = []
    Ibat = []
    Vbat = []
    Vgbat = []

    # init loop parameters
    dVdtisu = []

    # initialize loop parameters: speed at end of step and loop termination
    length_loop = 0
    time_i = time[-1]
    go_on = True

    while go_on:
        ## PART 1: compute parameters at end of step that are known without uncertainties:
        ##         *none*

        ## PART 2: compute parameters at end of step that are known only with uncertainties
        ##           (due to unknown mass and altitude at end of step):

        # Initialize loop parameters: aircraft mass (resp. altitude) at end of step is
        # first estimated as equal to aircraft mass (resp. altitude) at start of step

        mass_i = mass[-1]
        Hp_i = Hp[-1]
        v_i = TAS[-1]

        for _ in itertools.repeat(None, m_iter):
            # atmosphere properties
            H_m = conv.ft2m(Hp_i)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
            temp_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)

            step_time = length_loop - time[-1]

            # aircraft speed
            [M_i, CAS_i, TAS_i] = atm.convertSpeed(
                v=v_i, speedType="TAS", theta=theta, delta=delta, sigma=sigma
            )

            if turnFlight:
                if turnMetrics["bankAngle"] != 0.0:
                    # bankAngle is defined
                    rateOfTurn = AC.rateOfTurn_bankAngle(TAS=TAS_i, bankAngle=bankAngle)
                else:
                    # rateOfTurn is defined
                    bankAngle = AC.bankAngle(
                        rateOfTurn=rateOfTurn, v=TAS_i
                    )  # [degrees]

            # Load factor
            nz = 1 / cos(radians(bankAngle))

            # compute ROCD target (if any) on this step
            if control.slopetarget is not None:
                # compute target ROCD corresponding to target slope
                if AC.BADAFamily.BADAE:
                    # special case for BADAE, in future it may apply also for BADAH
                    dh_dt_i = TAS_i * tan(slopetargetisu)
                else:
                    dh_dt_i = TAS_i * sin(slopetargetisu)

                ROCDtargetisu_i = dh_dt_i * (1 / temp_const)
            elif control.ROCDtarget is not None:
                ROCDtargetisu_i = ROCDtargetisu
                dh_dt_i = ROCDtargetisu_i * temp_const

            # BADAH or BADAE
            if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
                # compute Power required
                Preq_i = AC.Preq(sigma=sigma, tas=TAS_i, mass=mass_i, phi=bankAngle)

                # compute engine power
                if rating is None:
                    # compiute power required for the manoeuver
                    if ESFc is not None:
                        P_i = dh_dt_i * mass_i * const.g / ESFc + Preq_i  # [W]
                    elif control.acctarget is not None:
                        P_i = (
                            dh_dt_i * mass_i * const.g
                            + mass_i * TAS_i * control.acctarget
                            + Preq_i
                        )  # [W]
                    else:
                        print("Error: neither ESF nor acc target provided")

                    # Check that required thrust/power fits in the available thrust/power envelope,
                    # recompute ROCD if necessary and compute fuel coefficient accordingly
                    Pmin = 0.1 * AC.P0  # No minimum power model: assume 10% torque

                    if AC.BADAFamily.BADAH:
                        Pmax = AC.Pav(rating=maxRating, theta=theta, delta=delta)
                        Pav_i = AC.Pav(rating=maxRating, theta=theta, delta=delta)
                    elif AC.BADAFamily.BADAE:
                        Pmax = AC.Pav(rating=maxRating, SOC=SOC[-1])
                        Pav_i = AC.Pav(rating=maxRating, SOC=SOC[-1])

                    if P_i < Pmin:
                        P_i = Pmin
                        if ESFc is not None:

                            ROCD_i = (
                                conv.m2ft(
                                    AC.ROCD(
                                        Peng=P_i,
                                        Preq=Preq_i,
                                        mass=mass_i,
                                        ESF=ESFc,
                                        theta=theta,
                                        DeltaTemp=DeltaTemp,
                                    )
                                )
                                * 60
                            )
                        elif control.acctarget is not None:
                            ROCD_i = (
                                conv.m2ft(
                                    (P_i - mass_i * TAS_i * control.acctarget - Preq_i)
                                    / (mass_i * const.g * temp_const)
                                )
                                * 60
                            )
                    elif P_i > Pmax:
                        P_i = Pmax

                        if ESFc is not None:
                            ROCD_i = (
                                conv.m2ft(
                                    AC.ROCD(
                                        Peng=P_i,
                                        Preq=Preq_i,
                                        mass=mass_i,
                                        ESF=ESFc,
                                        theta=theta,
                                        DeltaTemp=DeltaTemp,
                                    )
                                )
                                * 60
                            )
                        elif control.acctarget is not None:
                            ROCD_i = (
                                conv.m2ft(
                                    (P_i - mass_i * TAS_i * control.acctarget - Preq_i)
                                    / (mass_i * const.g * temp_const)
                                )
                                * 60
                            )
                    else:
                        ROCD_i = control.ROCDtarget

                else:
                    # Compute available power
                    if rating == "UNKNOWN":
                        P_i = 0.1 * AC.P0  # No minimum power model: assume 10% torque
                        Pav_i = AC.Pav(rating=maxRating, theta=theta, delta=delta)
                    else:
                        if AC.BADAFamily.BADAH:
                            P_i = AC.Pav(rating=rating, theta=theta, delta=delta)
                            Pav_i = AC.Pav(rating=rating, theta=theta, delta=delta)
                        elif AC.BADAFamily.BADAE:
                            P_i = AC.Pav(rating=rating, SOC=SOC[-1])
                            Pav_i = AC.Pav(rating=rating, SOC=SOC[-1])

                # Compute excess power
                Pe_i = P_i - Preq_i  # [W]

                # BADAH
                if AC.BADAFamily.BADAH:
                    # compute fuel flow for level flight
                    CP = AC.CP(Peng=P_i)
                    FUEL_i = AC.ff(delta=delta, CP=CP)  # [kg/s]

                # BADAE
                elif AC.BADAFamily.BADAE:
                    Pbat_i = AC.Pbat(Preq=P_i, SOC=SOC[-1])
                    SOCr_i = AC.SOCrate(Preq=P_i, SOC=SOC[-1])

                    # debug data
                    Pelc_i = P_i / AC.eta
                    Ibat_i = AC.Ibat(P=Pelc_i, SOC=SOC[-1])
                    Vbat_i = AC.Vbat(I=Ibat_i, SOC=SOC[-1])
                    Vgbat_i = AC.Vocbat(SOC=SOC[-1]) - AC.R0bat(SOC=SOC[-1]) * Ibat_i

            # BADA4
            elif AC.BADAFamily.BADA4:
                # aircraft aerodynamic configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=phase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=phase, previousConfig=config[-1], currentConfig=config_i
                    )

                [HLid_i, LG_i] = AC.flightEnvelope.getAeroConfig(config=config_i)

                # compute lift coefficient
                CL = AC.CL(M=M_i, delta=delta, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(M=M_i, CL=CL, HLid=HLid_i, LG=LG_i, speedBrakes=speedBrakes)
                # compute drag force
                Drag = AC.D(M=M_i, delta=delta, CD=CD)

                # compute thrust and fuel flow
                if rating is None:
                    # compute thrust force required for the manoeuver
                    if ESFc is not None:
                        THR_i = (
                            dh_dt_i * mass_i * const.g / (TAS_i * ESFc) + Drag
                        )  # [N]
                    elif control.acctarget is not None:
                        THR_i = (
                            dh_dt_i * mass_i * const.g / TAS_i
                            + mass_i * control.acctarget
                            + Drag
                        )  # [N]
                    else:
                        print("Error: neither ESF nor acc target provided")

                    # Check that required thrust fits in the available thrust envelope,
                    # recompute ROCD if necessary and compute fuel flow accordingly
                    THR_min = AC.Thrust(
                        rating="LIDL",
                        delta=delta,
                        theta=theta,
                        M=M_i,
                        DeltaTemp=DeltaTemp,
                    )  # IDLE Thrust
                    FUEL_min = AC.ff(
                        rating="LIDL",
                        delta=delta,
                        theta=theta,
                        M=M_i,
                        DeltaTemp=DeltaTemp,
                    )  # IDLE Fuel Flow
                    THR_max = AC.Thrust(
                        rating="MCMB",
                        delta=delta,
                        theta=theta,
                        M=M_i,
                        DeltaTemp=DeltaTemp,
                    )  # MCMB Thrust
                    FUEL_max = AC.ff(
                        rating="MCMB",
                        delta=delta,
                        theta=theta,
                        M=M_i,
                        DeltaTemp=DeltaTemp,
                    )  # MCMB Fuel Flow

                    if THR_i < THR_min:
                        THR_i = THR_min
                        FUEL_i = FUEL_min
                    elif THR_i > THR_max:
                        THR_i = THR_max
                        FUEL_i = FUEL_max
                    else:
                        CT = AC.CT(Thrust=THR_i, delta=delta)
                        FUEL_i = AC.ff(
                            CT=CT, delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                        )
                else:
                    THR_i = AC.Thrust(
                        rating=rating,
                        delta=delta,
                        theta=theta,
                        M=M_i,
                        DeltaTemp=DeltaTemp,
                    )  # [N]
                    CT = AC.CT(Thrust=THR_i, delta=delta)
                    FUEL_i = AC.ff(
                        rating="LIDL",
                        delta=delta,
                        theta=theta,
                        M=M_i,
                        DeltaTemp=DeltaTemp,
                    )
                    # FUEL_i = AC.ff(
                    # CT=CT, delta=delta, theta=theta, M=M_i, DeltaTemp=DeltaTemp
                    # )

                # compute excess power
                Pe_i = (THR_i - Drag) * TAS_i  # [kg*m^2/s^3]

            # BADA3
            elif AC.BADAFamily.BADA3:
                # aircraft aerodynamic configuration
                if config_default is None:
                    config_i = AC.flightEnvelope.getConfig(
                        h=H_m,
                        phase=phase,
                        v=CAS_i,
                        mass=mass_i,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    config_i = config_default

                # ensure continuity of configuration change within the segment
                if config:
                    config_i = AC.flightEnvelope.checkConfigurationContinuity(
                        phase=phase,
                        previousConfig=config[-1],
                        currentConfig=config_i,
                    )

                # compute lift coefficient
                CL = AC.CL(tas=TAS_i, sigma=sigma, mass=mass_i, nz=nz)
                # compute drag coefficient
                CD = AC.CD(CL=CL, config=config_i, speedBrakes=speedBrakes)
                # compute drag force
                Drag = AC.D(tas=TAS_i, sigma=sigma, CD=CD)

                # compute thrust and fuel flow
                if rating is None:
                    # compute thrust force required for the manoeuver
                    if (ESFc) is not None:
                        THR_i = (
                            dh_dt_i * mass_i * const.g / (TAS_i * ESFc) + Drag
                        )  # [N]
                    elif control.acctarget is not None:
                        THR_i = (
                            dh_dt_i * mass_i * const.g / TAS_i
                            + mass_i * control.acctarget
                            + Drag
                        )  # [N]
                    else:
                        print("Error: neither ESF nor acc target provided")

                    # Check that required thrust fits in the available thrust envelope,
                    # recompute ROCD if necessary and compute fuel flow accordingly

                    THR_min = AC.Thrust(
                        rating="LIDL", v=TAS_i, h=H_m, config="CR", DeltaTemp=DeltaTemp
                    )  # IDLE Thrust
                    FUEL_min = AC.ff(
                        flightPhase="Descent",
                        v=TAS_i,
                        h=H_m,
                        T=THR_min,
                        config="CR",
                        adapted=False,
                    )  # IDLE Fuel Flow
                    THR_max = AC.Thrust(
                        rating="MCMB", v=TAS_i, h=H_m, config="CR", DeltaTemp=DeltaTemp
                    )  # MCMB Thrust
                    FUEL_max = AC.ff(
                        flightPhase="Climb",
                        v=TAS_i,
                        h=H_m,
                        T=THR_max,
                        config="CR",
                        adapted=False,
                    )  # MCMB Fuel Flow

                    if THR_i < THR_min:
                        THR_i = THR_min
                        FUEL_i = FUEL_min
                    elif THR_i > THR_max:
                        THR_i = THR_max
                        FUEL_i = FUEL_max
                    else:
                        FUEL_i = AC.ff(
                            v=TAS_i, h=H_m, T=THR_i, config=config_i, adapted=True
                        )
                else:
                    THR_i = AC.Thrust(
                        rating=rating,
                        v=TAS_i,
                        h=H_m,
                        config=config_i,
                        DeltaTemp=DeltaTemp,
                    )
                    if rating == "MCMB" or rating == "MTKF":
                        FUEL_i = AC.ff(
                            flightPhase="Climb",
                            v=TAS_i,
                            h=H_m,
                            T=THR_i,
                            config=config_i,
                        )
                    elif rating == "MCRZ":
                        FUEL_i = AC.ff(
                            flightPhase="Cruise",
                            v=TAS_i,
                            h=H_m,
                            T=THR_i,
                            config=config_i,
                        )
                    elif rating == "LIDL":
                        FUEL_i = AC.ff(
                            flightPhase="Descent",
                            v=TAS_i,
                            h=H_m,
                            T=THR_i,
                            config=config_i,
                        )

                # compute excess power
                Pe_i = (THR_i - Drag) * TAS_i  # [kg*m^2/s^3]

            if ESFc is not None:
                ESF_i = ESFc
                # compute power dedicated to climb
                PC_i = Pe_i * ESF_i  # [kg*m^2/s^3]
                # compute ROCD
                dhdtisu = PC_i / (mass_i * const.g)  # [m/s]
                ROCDisu = dhdtisu * 1 / temp_const  # [m/s]
                ROCD_i = conv.m2ft(ROCDisu) * 60  # [ft/min]
            elif control.acctarget is not None:
                # compute power required for acc/dec rate
                Pa_i = mass_i * TAS_i * control.acctarget
                # check that required power fits in the available power envelope
                if abs(Pa_i) > abs(Pe_i):
                    Pa_i = Pe_i
                # compute power dedicated to climb
                PC_i = Pe_i - Pa_i  # [kg*m^2/s^3]

                if Pe_i != 0:
                    ESF_i = PC_i / Pe_i
                else:
                    ESF_i = float("Inf")
                # compute ROCD
                dhdtisu = PC_i / (mass_i * const.g)  # [m/s]
                ROCDisu = dhdtisu * 1 / temp_const  # [m/s]
                ROCD_i = conv.m2ft(ROCDisu) * 60  # [ft/min]
            elif control.slopetarget is not None or control.ROCDtarget is not None:
                dhdtisu = dh_dt_i  # [m/s]
                ROCDisu = dh_dt_i * 1 / temp_const  # [m/s]
                ROCD_i = conv.m2ft(ROCDisu) * 60  # [ft/min]
                PC_i = dh_dt_i * (mass_i * const.g)  # [kg*m^2/s^3]

                if Pe_i != 0:
                    ESF_i = PC_i / Pe_i
                else:
                    ESF_i = float("Inf")
            else:
                print("Error: unexpected combination of targets")

            # compute acceleration
            if TAS_i == 0:
                dVdtisu_i = (Pe_i - PC_i) / (mass_i * (TAS_i + 0.5))  # [m/s^2]
            else:
                dVdtisu_i = (Pe_i - PC_i) / (mass_i * TAS_i)  # [m/s^2]

            if length_loop == 0:
                # no need to loop for first point: initial m/Hp already known
                break
            else:
                # Average acceleration over step is the mean of initial and final ones
                step_dVdtisu = (dVdtisu[-1] + dVdtisu_i) / 2  # [m/s^2]
                # step speed is step acceleration multiplied by time step
                step_speed = conv.ms2kt(step_dVdtisu * step_time)
                # acceleration should be always defined in terms of TAS speed, so the next speed will be always defined as TAS
                v_i = TAS[-1] + step_speed
                # Average ROCD over step is the mean of initial and final ones
                step_ROCD = (ROCD[-1] + ROCD_i) / 2  # [ft/min]
                # Altitude differential is: average ROCD multiplied by step time
                step_Hp = step_ROCD * step_time / 60  # [ft]
                # Update altitude estimate at end of step
                Hp_i = Hp[-1] + step_Hp  # [ft]

                # BADAE
                if AC.BADAFamily.BADAE:
                    # Average SOC rate over step is the mean of initial and final ones
                    step_SOCr = (SOCr[-1] + SOCr_i) / 2  # [%/h]
                    # SOC change is: average SOC rate multiplied by step time
                    step_SOC = step_SOCr * step_time / 3600  # [%]
                    # Update SOC estimate at end of step
                    SOC_i = SOC[-1] - step_SOC  # [%]
                    # update of aircraft mass estimate at end of step - mass is not changing for ELECTRIC engine (no fuel is consumed)
                    mass_i = mass[-1]  # [kg]

                else:
                    # Average fuel flow over step is the mean of initial and final ones
                    step_FUEL = (FUEL[-1] + FUEL_i) / 2  # [kg/s]
                    # Fuel burnt is: average fuel flow multiplied by step time
                    step_mass = step_FUEL * step_time  # [kg]
                    # Update aircraft mass estimate at end of step
                    if not mass_const:
                        mass_i = mass[-1] - step_mass  # [kg]
                        fuelConsumed_i = step_FUEL * step_time
                    fuelConsumed_i = FUELCONSUMED[-1] + step_FUEL * step_time

        ## PART 3: store information about end of step point

        # point data
        dVdtisu.append(dVdtisu_i)
        ROCD.append(ROCD_i)
        esf.append(ESF_i)
        Comment.append(comment)

        # everything except electric BADAE
        if not AC.BADAFamily.BADAE:
            FUEL.append(FUEL_i)

        # BADAH
        if AC.BADAFamily.BADAH:
            Peng.append(P_i)
            Preq.append(Preq_i)
            Pav.append(Pav_i)

        # BADAE
        elif AC.BADAFamily.BADAE:
            Pmec.append(P_i)
            Pbat.append(Pbat_i)
            SOCr.append(SOCr_i)
            Pelc.append(Pelc_i)
            Ibat.append(Ibat_i)
            Vbat.append(Vbat_i)
            Vgbat.append(Vgbat_i)

        # BADA3 & BADA4
        elif AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
            THR.append(THR_i)
            DRAG.append(Drag)
            config.append(config_i)

        # BADA4
        if AC.BADAFamily.BADA4:
            HLid.append(HLid_i)
            LG.append(LG_i)

        # TEM consistency check
        # BADAH or BADAE
        if AC.BADAFamily.BADAH or AC.BADAFamily.BADAE:
            check.append(
                P_i - Preq_i - mass_i * const.g * dhdtisu - mass_i * TAS_i * dVdtisu_i
            )

        # BADA3 or BADA4
        elif AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
            check.append(
                (THR_i - Drag) * TAS_i
                - mass_i * const.g * dhdtisu
                - mass_i * TAS_i * dVdtisu_i
            )

        # calculation of the slope
        if TAS_i == 0:
            gamma_i = 90 * np.sign(ROCD_i)
        else:
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=conv.ft2m(Hp_i), DeltaTemp=DeltaTemp
            )
            temp_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)
            if AC.BADAFamily.BADAE:
                gamma_i = degrees(atan(conv.ft2m(ROCD_i) * temp_const / 60 / TAS_i))
            else:
                # using SIN assumes the TAS to be in the direction of the aircraft axis, not ground plane. Which means, this should be mathematically the correct equation for all the aircraft
                gamma_i = degrees(asin(conv.ft2m(ROCD_i) * temp_const / 60 / TAS_i))

        # ground speed can be calcualted as TAS projected on the ground minus wind
        GS_i = cos(radians(gamma_i)) * TAS_i - wS
        GS.append(conv.ms2kt(GS_i))

        Slope.append(gamma_i)
        BankAngle.append(bankAngle)
        ROT.append(rateOfTurn)

        # integrated data
        if length_loop != 0:  # exclude first point: initial t/d/m already known
            if AC.BADAFamily.BADAE:
                SOC.append(SOC_i)

            # everything except electric BADAE
            if not AC.BADAFamily.BADAE:
                FUELCONSUMED.append(fuelConsumed_i)

            # speed at the end of step
            TAS.append(conv.ms2kt(TAS_i))
            CAS.append(conv.ms2kt(CAS_i))
            M.append(M_i)
            # Altitude at end of step has been termined in PART 2
            Hp.append(Hp_i)
            # Aircraft mass at end of step has been termined in PART 2
            mass.append(mass_i)
            # Time at end of step is time at start of step plus step time
            time.append(time[-1] + step_time)

            # Average TAS over step is the mean of initial and final ones
            step_TAS = (TAS[-2] + TAS[-1]) / 2  # [kt]
            # Average slope over the step
            step_gamma = radians((Slope[-2] + Slope[-1]) / 2)  # radians
            # Average ground speed over step
            # since this is not level flight, TAS speed should be projected on the ground, then GS can be calculated applying the wind speed
            step_TAS_projected = cos(step_gamma) * step_TAS
            step_GS = step_TAS_projected - wS  # [kt]
            # Step distance is: average GS multiplied by step time
            if turnFlight:
                step_distance = conv.m2nm(
                    turn.distance(
                        rateOfTurn=rateOfTurn, TAS=TAS_i, timeOfTurn=step_time
                    )
                )  # arcLength during the turn [NM]
            else:
                step_distance = step_GS * step_time / 3600  # [NM]
            # Distance at end of step is distance at start of step plus step distance
            dist.append(dist[-1] + step_distance)

            # add GPS calculation
            if Lat and Lon and (magneticHeading or trueHeading):
                if headingToFly == "TRUE":
                    if not turnFlight:
                        if not constantHeading:
                            # fly ORTHODROME
                            (Lat_i, Lon_i, HDGTrue_i) = (
                                vincenty.destinationPoint_finalBearing(
                                    LAT_init=LAT[-1],
                                    LON_init=LON[-1],
                                    distance=conv.nm2m(step_distance),
                                    bearing=HDGTrue[-1],
                                )
                            )

                            if magneticDeclinationGrid is not None:
                                HDGMagnetic_i = (
                                    HDGTrue_i
                                    - magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGMagnetic_i = HDGTrue_i

                        elif constantHeading:
                            # fly LOXODROME
                            (Lat_i, Lon_i) = rhumb.destinationPoint(
                                LAT_init=LAT[-1],
                                LON_init=LON[-1],
                                bearing=HDGTrue[-1],
                                distance=conv.nm2m(step_distance),
                            )
                            HDGTrue_i = HDGTrue[-1]

                            if magneticDeclinationGrid is not None:
                                HDGMagnetic_i = (
                                    HDGTrue_i
                                    - magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGMagnetic_i = HDGTrue_i

                    else:
                        # calculate the turn
                        (Lat_i, Lon_i, HDGTrue_i) = turn.destinationPoint_finalBearing(
                            LAT_init=LAT[-1],
                            LON_init=LON[-1],
                            bearingInit=HDGTrue[-1],
                            TAS=TAS_i,
                            rateOfTurn=rateOfTurn,
                            timeOfTurn=step_time,
                            directionOfTurn=directionOfTurn,
                        )

                        if magneticDeclinationGrid is not None:
                            HDGMagnetic_i = (
                                HDGTrue_i
                                - magneticDeclinationGrid.getMagneticDeclination(
                                    LAT_target=Lat_i, LON_target=Lon_i
                                )
                            )
                        else:
                            magneticDeclination = 0
                            HDGMagnetic_i = HDGTrue_i

                elif headingToFly == "MAGNETIC":
                    if not turnFlight:
                        if constantHeading:
                            (Lat_i, Lon_i) = rhumb.destinationPoint(
                                LAT_init=LAT[-1],
                                LON_init=LON[-1],
                                bearing=HDGTrue[-1],
                                distance=conv.nm2m(step_distance),
                            )
                            HDGMagnetic_i = HDGMagnetic[-1]

                            if magneticDeclinationGrid is not None:
                                HDGTrue_i = (
                                    HDGMagnetic_i
                                    + magneticDeclinationGrid.getMagneticDeclination(
                                        LAT_target=Lat_i, LON_target=Lon_i
                                    )
                                )
                            else:
                                magneticDeclination = 0
                                HDGTrue_i = HDGMagnetic_i

                    else:
                        # calculate the turn
                        (Lat_i, Lon_i, HDGTrue_i) = turn.destinationPoint_finalBearing(
                            LAT_init=LAT[-1],
                            LON_init=LON[-1],
                            bearingInit=HDGTrue[-1],
                            TAS=TAS_i,
                            rateOfTurn=rateOfTurn,
                            timeOfTurn=step_time,
                            directionOfTurn=directionOfTurn,
                        )

                        if magneticDeclinationGrid is not None:
                            HDGMagnetic_i = (
                                HDGTrue_i
                                - magneticDeclinationGrid.getMagneticDeclination(
                                    LAT_target=Lat_i, LON_target=Lon_i
                                )
                            )
                        else:
                            magneticDeclination = 0
                            HDGMagnetic_i = HDGTrue_i

                LAT.append(Lat_i)
                LON.append(Lon_i)
                HDGMagnetic.append(HDGMagnetic_i)
                HDGTrue.append(HDGTrue_i)

        if length_loop + step_length < length:
            length_loop += step_length
        elif length_loop < length:
            length_loop = length
        else:
            go_on = False

    flightData = {
        "Hp": Hp,
        "TAS": TAS,
        "CAS": CAS,
        "GS": GS,
        "M": M,
        "acc": dVdtisu,
        "ROCD": ROCD,
        "ESF": esf,
        "FUEL": FUEL,
        "Pmec": Pmec,
        "Pelc": Pelc,
        "Pbat": Pbat,
        "SOCr": SOCr,
        "SOC": SOC,
        "Ibat": Ibat,
        "Vbat": Vbat,
        "Vgbat": Vgbat,
        "FUELCONSUMED": FUELCONSUMED,
        "Preq": Preq,
        "Peng": Peng,
        "Pav": Pav,
        "THR": THR,
        "DRAG": DRAG,
        "time": time,
        "dist": dist,
        "slope": Slope,
        "mass": mass,
        "config": config,
        "HLid": HLid,
        "LG": LG,
        "LAT": LAT,
        "LON": LON,
        "HDGTrue": HDGTrue,
        "HDGMagnetic": HDGMagnetic,
        "BankAngle": BankAngle,
        "ROT": ROT,
        "comment": Comment,
    }

    flightTrajectory = FT.createFlightTrajectoryDataframe(flightData)
    return flightTrajectory
