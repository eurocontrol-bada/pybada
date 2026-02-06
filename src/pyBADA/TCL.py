"""Trajectory Computation Light (TCL) for BADAH/BADA3/BADA4 aircraft
performance module."""

import numpy as np

from pyBADA import atmosphere as atm
from pyBADA import constants as const
from pyBADA import conversions as conv
from pyBADA import myTypes, trajectorySegments
from pyBADA.flightTrajectory import FlightTrajectory as FT
from pyBADA.myTypes import (
    AccelerationLevelKind,
    ArrivalProfileType,
    CalculationType,
    ClimbType,
    CruiseSpeedType,
    CruiseType,
    DepartureProfileType,
    DescentType,
    FlightPhase,
    HClimbType,
    HDescentType,
    IntegrationType,
    SpeedType,
    TakeOffProcedureBADA,
    TakeOffProcedureNADP1,
    TakeOffProcedureNADP2,
)


def climbDescentRate(
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo,
    controlTarget: myTypes.ControlTarget,
    speedBrakes: myTypes.SpeedBrakes | None = None,
) -> FT:
    """Calculates the flight trajectory in climb or descent at a specific ROCD.

    :param AC: Aircraft model instance containing performance and aerodynamic data.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing initial, final, and step values.
    :param speed: Speed configuration including the speed type and initial value.
    :param mass: Current aircraft mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param controlTarget: Target objective for the segment, specifically target ROCD.
    :param speedBrakes: Configuration for speed brake state and deployment value.

    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo
    :type controlTarget: myTypes.ControlTarget
    :type speedBrakes: myTypes.SpeedBrakes | None

    :returns: A Flight Trajectory object containing the computed results.
    :rtype: FT
    """

    trajectory = FT()

    if speedBrakes is None:
        speedBrakes = myTypes.SpeedBrakes(deployed=0, value=0)
    else:
        if speedBrakes.value is None:
            speedBrakesValue = 0.0
        else:
            speedBrakesValue = abs(min(speedBrakes.value, 100) * (0.03 / 100))

        speedBrakes = {
            "deployed": speedBrakes.deployed,
            "value": speedBrakesValue,
        }

    flightTrajectory = trajectorySegments.constantSpeedROCD(
        AC=AC,
        speedType=speed.speedType,
        v=speed.initSpeed,
        ROCDtarget=controlTarget.ROCDtarget,
        Hp_init=pressureAltitude.initPressureAltitude,
        Hp_final=pressureAltitude.finalPressureAltitude,
        Hp_step=pressureAltitude.stepPressureAltitude,
        m_init=mass,
        wS=meteo.wS,
        deltaTemp=meteo.deltaTemp,
        speedBrakes=speedBrakes,
        calculationType=calculationType,
    )
    trajectory.append(AC, flightTrajectory)

    return trajectory


def climbDescentSlope(
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo,
    controlTarget: myTypes.ControlTarget,
    speedBrakes: myTypes.SpeedBrakes | None = None,
) -> FT:
    """Calculates the flight trajectory in climb or descent at a specific slope.

    :param AC: Aircraft model instance containing performance and aerodynamic data.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing initial, final, and step values.
    :param speed: Speed configuration including the speed type and initial value.
    :param mass: Current aircraft mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param controlTarget: Target objective for the segment, specifically the target slope.
    :param speedBrakes: Configuration for speed brake state and deployment value.

    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo
    :type controlTarget: myTypes.ControlTarget
    :type speedBrakes: myTypes.SpeedBrakes | None

    :returns: A Flight Trajectory object containing the computed results.
    :rtype: FT
    """

    trajectory = FT()

    if speedBrakes is None:
        speedBrakes = myTypes.SpeedBrakes(deployed=0, value=0)
    else:
        if speedBrakes.value is None:
            speedBrakesValue = 0.0
        else:
            speedBrakesValue = abs(min(speedBrakes.value, 100) * (0.03 / 100))

        speedBrakes = {
            "deployed": speedBrakes.deployed,
            "value": speedBrakesValue,
        }

    flightTrajectory = trajectorySegments.constantSpeedSlope(
        AC=AC,
        speedType=speed.speedType,
        v=speed.initSpeed,
        slopetarget=controlTarget.slopetarget,
        Hp_init=pressureAltitude.initPressureAltitude,
        Hp_final=pressureAltitude.finalPressureAltitude,
        Hp_step=pressureAltitude.stepPressureAltitude,
        m_init=mass,
        wS=meteo.wS,
        deltaTemp=meteo.deltaTemp,
        speedBrakes=speedBrakes,
        calculationType=calculationType,
    )
    trajectory.append(AC, flightTrajectory)

    return trajectory


def climbDescentAccelerationDeceleration(
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    flightPhase: myTypes.FlightPhase,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo,
    controlTarget: myTypes.ControlTarget,
    speedBrakes: myTypes.SpeedBrakes | None = None,
) -> FT:
    """Calculates the flight trajectory for acceleration or deceleration during climb or descent.

    :param AC: Aircraft model instance containing performance and aerodynamic data.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing the initial pressure altitude.
    :param flightPhase: The current phase of flight (CLIMB or DESCENT).
    :param speed: Speed configuration including initial, final, and step speed values.
    :param mass: Current aircraft mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param controlTarget: Target parameters governing the acceleration/deceleration.
    :param speedBrakes: Configuration for speed brake state and deployment value.

    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type flightPhase: myTypes.FlightPhase
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo
    :type controlTarget: myTypes.ControlTarget
    :type speedBrakes: myTypes.SpeedBrakes | None

    :returns: A Flight Trajectory object containing the computed results.
    :rtype: FT
    """

    trajectory = FT()

    if speedBrakes is None:
        speedBrakes = myTypes.SpeedBrakes(deployed=0, value=0)
    else:
        if speedBrakes.value is None:
            speedBrakesValue = 0.0
        else:
            speedBrakesValue = abs(min(speedBrakes.value, 100) * (0.03 / 100))

        speedBrakes = {
            "deployed": speedBrakes.deployed,
            "value": speedBrakesValue,
        }

    if calculationType == calculationType.POINT:
        stepSpeed = abs(speed.finalSpeed - speed.initSpeed)
    else:
        stepSpeed = speed.stepSpeed

    if flightPhase == FlightPhase.CLIMB:
        phase = "Climb"
    elif flightPhase == FlightPhase.DESCENT:
        phase = "Descent"

    flightTrajectory = trajectorySegments.accDec(
        AC=AC,
        speedType=speed.speedType,
        v_init=speed.initSpeed,
        v_final=speed.finalSpeed,
        speed_step=stepSpeed,
        phase=phase,
        Hp_init=pressureAltitude.initPressureAltitude,
        m_init=mass,
        deltaTemp=meteo.deltaTemp,
        wS=meteo.wS,
        control=controlTarget,
        speedBrakes=speedBrakes,
        calculationType=calculationType,
    )
    trajectory.append(AC, flightTrajectory)

    return trajectory


def apcDescentEmergency(
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo,
    controlTarget: myTypes.ControlTarget,
    speedBrakes: myTypes.SpeedBrakes,
) -> FT:
    """Calculates an emergency descent trajectory at maximum operating speeds (VMO/MMO).

    :param AC: Aircraft model instance containing VMO, MMO, and performance limits.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing initial, final, and step values.
    :param speed: Speed configuration (Note: Function defaults to AC.VMO/AC.MMO).
    :param mass: Current aircraft mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param controlTarget: Target parameters for the descent segment.
    :param speedBrakes: Configuration for speed brake state and deployment value.

    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo
    :type controlTarget: myTypes.ControlTarget
    :type speedBrakes: myTypes.SpeedBrakes

    :returns: A Flight Trajectory object representing the high-speed emergency descent.
    :rtype: FT
    """

    trajectory = FT()

    initPressureAltitude = pressureAltitude.initPressureAltitude
    finalPressureAltitude = pressureAltitude.finalPressureAltitude
    stepPressureAltitude = pressureAltitude.stepPressureAltitude

    # descent at max speed
    # calculate crossover altitude [ft]
    crossoverAltitude = conv.m2ft(atm.crossOver(conv.kt2ms(AC.VMO), AC.MMO))

    # trajectory above crossover altitude
    if (
        initPressureAltitude > crossoverAltitude
        and finalPressureAltitude > crossoverAltitude
    ):
        flightTrajectory = trajectorySegments.constantSpeedRating(
            AC=AC,
            speedType="M",
            v=AC.MMO,
            Hp_init=pressureAltitude.initPressureAltitude,
            Hp_final=pressureAltitude.finalPressureAltitude,
            Hp_step=pressureAltitude.stepPressureAltitude,
            m_init=mass,
            deltaTemp=meteo.deltaTemp,
            wS=meteo.wS,
            calculationType=calculationType,
            expedite=True,
        )
        trajectory.append(AC, flightTrajectory)

    # trajectory above and below crossover altitude
    elif (
        initPressureAltitude > crossoverAltitude
        and finalPressureAltitude < crossoverAltitude
    ):
        flightTrajectory = trajectorySegments.constantSpeedRating(
            AC=AC,
            speedType="M",
            v=AC.MMO,
            Hp_init=pressureAltitude.initPressureAltitude,
            Hp_final=pressureAltitude.finalPressureAltitude,
            Hp_step=pressureAltitude.stepPressureAltitude,
            m_init=mass,
            deltaTemp=meteo.deltaTemp,
            wS=meteo.wS,
            calculationType=calculationType,
            expedite=True,
        )
        trajectory.append(AC, flightTrajectory)
        (CAS_current, Hp_current, mass_current, config_current) = (
            trajectory.getFinalValue(AC, ["CAS", "Hp", "mass", "config"])
        )
        if calculationType == "POINT" and Hp_current != finalPressureAltitude:
            trajectory.removeLines(AC, numberOfLines=1)

        flightTrajectory = trajectorySegments.constantSpeedRating(
            AC=AC,
            speedType="CAS",
            v=AC.VMO,
            Hp_init=pressureAltitude.initPressureAltitude,
            Hp_final=pressureAltitude.finalPressureAltitude,
            Hp_step=pressureAltitude.stepPressureAltitude,
            m_init=mass_current,
            deltaTemp=meteo.deltaTemp,
            wS=meteo.wS,
            calculationType=calculationType,
            expedite=True,
        )
        trajectory.append(AC, flightTrajectory)

    # trajectory below crossover altitude
    else:
        flightTrajectory = trajectorySegments.constantSpeedRating(
            AC=AC,
            speedType="CAS",
            v=AC.VMO,
            Hp_init=pressureAltitude.initPressureAltitude,
            Hp_final=pressureAltitude.finalPressureAltitude,
            Hp_step=pressureAltitude.stepPressureAltitude,
            m_init=mass,
            deltaTemp=meteo.deltaTemp,
            wS=meteo.wS,
            calculationType=calculationType,
            expedite=True,
        )
        trajectory.append(AC, flightTrajectory)

    return trajectory


def hpcDescentEmergency(
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo,
    controlTarget: myTypes.ControlTarget,
) -> FT:
    """Calculates an emergency descent trajectory for helicopters at VNE.

    :param AC: Aircraft model instance containing VNE and helicopter performance data.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing initial, final, and step values.
    :param speed: Speed configuration (Note: Function defaults to AC.vne).
    :param mass: Current aircraft mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param controlTarget: Target parameters for the descent segment.

    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo
    :type controlTarget: myTypes.ControlTarget

    :returns: A Flight Trajectory object representing the high-speed helicopter descent.
    :rtype: FT
    """

    trajectory = FT()

    initPressureAltitude = pressureAltitude.initPressureAltitude
    finalPressureAltitude = pressureAltitude.finalPressureAltitude
    stepPressureAltitude = pressureAltitude.stepPressureAltitude

    # descent at max speed
    flightTrajectory = trajectorySegments.constantSpeedRating(
        AC=AC,
        speedType="CAS",
        v=AC.vne,
        Hp_init=pressureAltitude.initPressureAltitude,
        Hp_final=pressureAltitude.finalPressureAltitude,
        Hp_step=pressureAltitude.stepPressureAltitude,
        m_init=mass,
        deltaTemp=meteo.deltaTemp,
        wS=meteo.wS,
        calculationType=calculationType,
    )
    trajectory.append(AC, flightTrajectory)

    return trajectory


def apcClimbCasMach(
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo,
    controlTarget: myTypes.ControlTarget,
    speedBrakes: myTypes.SpeedBrakes,
    casMachSpeedSchedule: myTypes.CASMACHSpeedSchedule,
    takeOffProcedure: myTypes.TakeOffProcedure,
    departureProfile: myTypes.DepartureProfile,
    reducedPower: bool,
) -> FT:
    """Calculates a complete climb trajectory using a CAS/Mach speed schedule.

    This function simulates a multi-segment climb profile by determining
    threshold altitudes for configuration changes and speed transitions. It
    incorporates BADA speed schedules or custom user inputs, manages the
    transition at the crossover altitude, and accounts for specific takeoff
    procedures (e.g., NADP1, NADP2) and engine-specific climb limits.

    :param AC: Aircraft model instance containing performance limits and engine data.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing initial, final, and step values.
    :param speed: Speed configuration including acceleration level types and steps.
    :param mass: Current aircraft mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param controlTarget: Target parameters for the climb segments.
    :param speedBrakes: Configuration for speed brake state and deployment value.
    :param casMachSpeedSchedule: Object defining CAS below/above FL100 and target Mach.
    :param takeOffProcedure: The specific takeoff/climb procedure (BADA, NADP1, or NADP2).
    :param departureProfile: Profile type definition (e.g., constant CAS below 10,000ft).
    :param reducedPower: Boolean flag to indicate if reduced thrust/power should be used.

    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo
    :type controlTarget: myTypes.ControlTarget
    :type speedBrakes: myTypes.SpeedBrakes
    :type casMachSpeedSchedule: myTypes.CASMACHSpeedSchedule
    :type takeOffProcedure: myTypes.TakeOffProcedure
    :type departureProfile: myTypes.DepartureProfile
    :type reducedPower: bool

    :returns: A Flight Trajectory object containing the concatenated climb segments.
    :rtype: FT
    """

    trajectory = FT()

    massCurrent = mass

    # setting the current aero configuration
    config_current = None

    stepPressureAltitude = pressureAltitude.stepPressureAltitude

    # determine the BADA speed schedule [m/s and Mach]
    CASbelowFL100 = casMachSpeedSchedule.CASbelowFL100
    CASaboveFL100 = casMachSpeedSchedule.CASaboveFL100
    Mach = casMachSpeedSchedule.Mach

    # revert to default BADA speed schedule if not defined
    if CASbelowFL100 is None:
        [Vcl1, Vcl2_temp, Mcl_temp] = AC.flightEnvelope.getSpeedSchedule(
            phase="Climb"
        )
    else:
        Vcl1 = conv.kt2ms(CASbelowFL100)

    if CASaboveFL100 is None:
        [Vcl1_temp, Vcl2, Mcl_temp] = AC.flightEnvelope.getSpeedSchedule(
            phase="Climb"
        )
    else:
        Vcl2 = conv.kt2ms(CASaboveFL100)

    if Mach is None:
        [Vcl1_temp, Vcl2_temp, Mcl] = AC.flightEnvelope.getSpeedSchedule(
            phase="Climb"
        )
    else:
        Mcl = Mach

    speedSchedule = [Vcl1, Vcl2, Mcl]

    # calculate crossover altitude [ft]
    crossoverAltitude = conv.m2ft(atm.crossOver(Vcl2, Mcl))

    # calculate tropopause altitude [ft]
    tropopauseAltitude = conv.m2ft(const.h_11)

    # Collect special altitudes that need to be included.
    specialAltitudes = [crossoverAltitude, tropopauseAltitude]
    configChangesGPFAltitudes = [400, 2000]

    if (
        departureProfile.departureProfileType
        == DepartureProfileType.constCASbelow100
    ):
        procedureAltitudeList = [10000]
    else:
        match takeOffProcedure:
            case TakeOffProcedureBADA():
                if AC.engineType == "JET":
                    procedureAltitudeList = [
                        1500,
                        3000,
                        4000,
                        5000,
                        6000,
                        10000,
                    ]

                elif AC.engineType in ("TURBOPROP", "PISTON", "ELECTRIC"):
                    procedureAltitudeList = [500, 1000, 1500, 10000]

                else:
                    procedureAltitudeList = []  # Handle unknown engines safely
            case TakeOffProcedureNADP1(NADP1Threshold=threshold):
                procedureAltitudeList = [threshold, 10000]
            case TakeOffProcedureNADP2(
                NADP2Threshold1=threshold1, NADP2Threshold2=threshold2
            ):
                procedureAltitudeList = [threshold1, threshold2, 10000]

    # Merge the arrays, remove duplicates, and sort.
    altitudeClimbThresholdList = np.sort(
        np.unique(
            np.concatenate(
                (
                    procedureAltitudeList,
                    configChangesGPFAltitudes,
                    specialAltitudes,
                    [pressureAltitude.finalPressureAltitude],
                )
            )
        )
    )

    # Now filter the list: remove altitudes below initPressureAltitude
    # and above finalPressureAltitude.
    altitudeClimbThresholdList = altitudeClimbThresholdList[
        (altitudeClimbThresholdList >= pressureAltitude.initPressureAltitude)
        & (
            altitudeClimbThresholdList
            <= pressureAltitude.finalPressureAltitude
        )
    ]

    # set current altitude
    Hp_current = pressureAltitude.initPressureAltitude

    if Hp_current < crossoverAltitude:
        if (
            departureProfile.departureProfileType
            == DepartureProfileType.constCASbelow100
        ):
            CAS_current = conv.ms2kt(Vcl1)

        else:
            # set current speed
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=conv.ft2m(Hp_current), deltaTemp=meteo.deltaTemp
            )

            match takeOffProcedure:
                case TakeOffProcedureBADA():
                    CAS_current = conv.ms2kt(
                        AC.ARPM.climbSpeed(
                            h=conv.ft2m(Hp_current),
                            mass=massCurrent,
                            theta=theta,
                            delta=delta,
                            deltaTemp=meteo.deltaTemp,
                            speedSchedule_default=speedSchedule,
                        )[0]
                    )
                case TakeOffProcedureNADP1(NADP1Threshold=threshold):
                    if threshold == 0.0:
                        CAS_current = conv.ms2kt(
                            AC.ARPM.climbSpeed(
                                h=conv.ft2m(Hp_current),
                                mass=massCurrent,
                                theta=theta,
                                delta=delta,
                                deltaTemp=meteo.deltaTemp,
                                speedSchedule_default=speedSchedule,
                                procedure=takeOffProcedure.type,
                            )[0]
                        )
                    else:
                        CAS_current = conv.ms2kt(
                            AC.ARPM.climbSpeed(
                                h=conv.ft2m(Hp_current),
                                mass=massCurrent,
                                theta=theta,
                                delta=delta,
                                deltaTemp=meteo.deltaTemp,
                                speedSchedule_default=speedSchedule,
                                procedure=takeOffProcedure.type,
                                NADP1_ALT=threshold,
                            )[0]
                        )
                case TakeOffProcedureNADP2(
                    NADP2Threshold1=threshold1, NADP2Threshold2=threshold2
                ):
                    if threshold1 == 0.0 or threshold2 == 0.0:
                        CAS_current = conv.ms2kt(
                            AC.ARPM.climbSpeed(
                                h=conv.ft2m(Hp_current),
                                mass=massCurrent,
                                theta=theta,
                                delta=delta,
                                deltaTemp=meteo.deltaTemp,
                                speedSchedule_default=speedSchedule,
                                procedure=takeOffProcedure.type,
                            )[0]
                        )
                    else:
                        CAS_current = conv.ms2kt(
                            AC.ARPM.climbSpeed(
                                h=conv.ft2m(Hp_current),
                                mass=massCurrent,
                                theta=theta,
                                delta=delta,
                                deltaTemp=meteo.deltaTemp,
                                speedSchedule_default=speedSchedule,
                                procedure=takeOffProcedure.type,
                                NADP2_ALT=[threshold1, threshold2],
                            )[0]
                        )
    else:
        M_current = Mcl

    for Hp_next in altitudeClimbThresholdList:
        # set current CAS speed [kt]
        [theta, delta, sigma] = atm.atmosphereProperties(
            h=conv.ft2m(Hp_next), deltaTemp=meteo.deltaTemp
        )

        if (
            departureProfile.departureProfileType
            == DepartureProfileType.constCASbelow100
            and Hp_next <= 10000
        ):
            CAS_final = CAS_current
        else:
            match takeOffProcedure:
                case TakeOffProcedureBADA():
                    CAS_final = conv.ms2kt(
                        AC.ARPM.climbSpeed(
                            h=conv.ft2m(Hp_next),
                            mass=massCurrent,
                            theta=theta,
                            delta=delta,
                            deltaTemp=meteo.deltaTemp,
                            speedSchedule_default=speedSchedule,
                        )[0]
                    )
                case TakeOffProcedureNADP1(NADP1Threshold=threshold):
                    if threshold == 0.0:
                        CAS_final = conv.ms2kt(
                            AC.ARPM.climbSpeed(
                                h=conv.ft2m(Hp_next),
                                mass=massCurrent,
                                theta=theta,
                                delta=delta,
                                deltaTemp=meteo.deltaTemp,
                                speedSchedule_default=speedSchedule,
                                procedure=takeOffProcedure.type,
                            )[0]
                        )
                    else:
                        CAS_final = conv.ms2kt(
                            AC.ARPM.climbSpeed(
                                h=conv.ft2m(Hp_next),
                                mass=massCurrent,
                                theta=theta,
                                delta=delta,
                                deltaTemp=meteo.deltaTemp,
                                speedSchedule_default=speedSchedule,
                                procedure=takeOffProcedure.type,
                                NADP1_ALT=threshold,
                            )[0]
                        )
                case TakeOffProcedureNADP2(
                    NADP2Threshold1=threshold1, NADP2Threshold2=threshold2
                ):
                    if threshold1 == 0.0 or threshold2 == 0.0:
                        CAS_final = conv.ms2kt(
                            AC.ARPM.climbSpeed(
                                h=conv.ft2m(Hp_next),
                                mass=massCurrent,
                                theta=theta,
                                delta=delta,
                                deltaTemp=meteo.deltaTemp,
                                speedSchedule_default=speedSchedule,
                                procedure=takeOffProcedure.type,
                            )[0]
                        )
                    else:
                        CAS_final = conv.ms2kt(
                            AC.ARPM.climbSpeed(
                                h=conv.ft2m(Hp_next),
                                mass=massCurrent,
                                theta=theta,
                                delta=delta,
                                deltaTemp=meteo.deltaTemp,
                                speedSchedule_default=speedSchedule,
                                procedure=takeOffProcedure.type,
                                NADP2_ALT=[threshold1, threshold2],
                            )[0]
                        )

        if Hp_next < crossoverAltitude:
            Hp_final = Hp_next

            # check if the next calculated speed is smaller then my current speed, to avoid deceleration during climb
            # if CAS_final < CAS_current:
            # CAS_final = CAS_current

            # how much altitude do I need to Accelerate to next threshold altitude speed?

            if speed.accelerationLevelKind == AccelerationLevelKind.AT:
                # climb to set altitude and accelerate then to reach the next threshold
                flightTrajectory = trajectorySegments.constantSpeedRating(
                    AC=AC,
                    speedType=SpeedType.CAS,
                    v=CAS_current,
                    Hp_init=Hp_current,
                    Hp_final=Hp_final,
                    Hp_step=stepPressureAltitude,
                    m_init=massCurrent,
                    wS=meteo.wS,
                    deltaTemp=meteo.deltaTemp,
                    reducedPower=reducedPower,
                    calculationType=calculationType,
                )
                trajectory.append(AC, flightTrajectory)
                CAS_current, Hp_current, massCurrent, config_current = (
                    trajectory.getFinalValue(
                        AC, ["CAS", "Hp", "mass", "config"]
                    )
                )

                if (
                    calculationType == CalculationType.POINT
                    and Hp_current != pressureAltitude.finalPressureAltitude
                ):
                    trajectory.removeLines(AC, numberOfLines=1)

                if (
                    abs(CAS_current - CAS_final) > 0.3
                ):  # preventing speed jumps due to small accuracy issues
                    flightTrajectory = trajectorySegments.accDec(
                        AC=AC,
                        speedType=SpeedType.CAS,
                        v_init=CAS_current,
                        v_final=CAS_final,
                        speed_step=speed.stepSpeed,
                        Hp_init=Hp_current,
                        control=controlTarget,
                        phase="Cruise",
                        config=config_current,
                        m_init=massCurrent,
                        wS=meteo.wS,
                        deltaTemp=meteo.deltaTemp,
                        calculationType=calculationType,
                    )
                    trajectory.append(AC, flightTrajectory)
                    CAS_current, Hp_current, massCurrent, config_current = (
                        trajectory.getFinalValue(
                            AC, ["CAS", "Hp", "mass", "config"]
                        )
                    )
                    if (
                        calculationType == CalculationType.POINT
                        and Hp_current
                        != pressureAltitude.finalPressureAltitude
                    ):
                        trajectory.removeLines(AC, numberOfLines=1)

            elif speed.accelerationLevelKind == AccelerationLevelKind.BEFORE:
                Hp_final = Hp_current
                diff = Hp_next - Hp_current
                while abs(diff) > 1:
                    traj = FT()
                    flightTrajectory = trajectorySegments.accDec(
                        AC=AC,
                        speedType=SpeedType.CAS,
                        v_init=CAS_current,
                        v_final=CAS_final,
                        speed_step=speed.stepSpeed,
                        Hp_init=Hp_final,
                        control=controlTarget,
                        phase="Climb",
                        config=config_current,
                        m_init=massCurrent,
                        wS=meteo.wS,
                        deltaTemp=meteo.deltaTemp,
                        calculationType=calculationType,
                    )
                    traj.append(AC, flightTrajectory)
                    CAS, Hp_end = traj.getFinalValue(AC, ["CAS", "Hp"])
                    diff = Hp_next - Hp_end
                    Hp_final = Hp_final + diff

                # climb to set altitude and accelerate then to reach the next threshold
                flightTrajectory = trajectorySegments.constantSpeedRating(
                    AC=AC,
                    speedType=SpeedType.CAS,
                    v=CAS_current,
                    Hp_init=Hp_current,
                    Hp_final=Hp_final,
                    Hp_step=stepPressureAltitude,
                    m_init=massCurrent,
                    wS=meteo.wS,
                    deltaTemp=meteo.deltaTemp,
                    reducedPower=reducedPower,
                    calculationType=calculationType,
                )
                trajectory.append(AC, flightTrajectory)
                CAS_current, Hp_current, massCurrent, config_current = (
                    trajectory.getFinalValue(
                        AC, ["CAS", "Hp", "mass", "config"]
                    )
                )
                if (
                    calculationType == CalculationType.POINT
                    and Hp_current != pressureAltitude.finalPressureAltitude
                ):
                    trajectory.removeLines(AC, numberOfLines=1)

                if (
                    abs(CAS_current - CAS_final) > 0.3
                ):  # preventing speed jumps due to small accuracy issues
                    flightTrajectory = trajectorySegments.accDec(
                        AC=AC,
                        speedType=SpeedType.CAS,
                        v_init=CAS_current,
                        v_final=CAS_final,
                        speed_step=speed.stepSpeed,
                        Hp_init=Hp_current,
                        control=controlTarget,
                        phase="Climb",
                        config=config_current,
                        m_init=massCurrent,
                        wS=meteo.wS,
                        deltaTemp=meteo.deltaTemp,
                        calculationType=calculationType,
                    )
                    trajectory.append(AC, flightTrajectory)
                    CAS_current, Hp_current, massCurrent, config_current = (
                        trajectory.getFinalValue(
                            AC, ["CAS", "Hp", "mass", "config"]
                        )
                    )
                    if abs(Hp_current - Hp_next) < 1:
                        trajectory.overwriteLastValue(AC, "Hp", Hp_next)
                        Hp_current = Hp_next
                    if (
                        calculationType == CalculationType.POINT
                        and Hp_current
                        != pressureAltitude.finalPressureAltitude
                    ):
                        trajectory.removeLines(AC, numberOfLines=1)

            elif speed.accelerationLevelKind == AccelerationLevelKind.AFTER:
                # climb to set altitude and Accelerate then to reach the next threshold
                flightTrajectory = trajectorySegments.constantSpeedRating(
                    AC=AC,
                    speedType=SpeedType.CAS,
                    v=CAS_current,
                    Hp_init=Hp_current,
                    Hp_final=Hp_next,
                    Hp_step=stepPressureAltitude,
                    m_init=massCurrent,
                    wS=meteo.wS,
                    deltaTemp=meteo.deltaTemp,
                    reducedPower=reducedPower,
                    calculationType=calculationType,
                )
                trajectory.append(AC, flightTrajectory)
                CAS_current, Hp_current, massCurrent, config_current = (
                    trajectory.getFinalValue(
                        AC, ["CAS", "Hp", "mass", "config"]
                    )
                )
                if (
                    calculationType == CalculationType.POINT
                    and Hp_current != pressureAltitude.finalPressureAltitude
                ):
                    trajectory.removeLines(AC, numberOfLines=1)

                if (
                    abs(CAS_current - CAS_final) > 0.3
                ):  # preventing speed jumps due to small accuracy issues
                    flightTrajectory = trajectorySegments.accDec(
                        AC=AC,
                        speedType=SpeedType.CAS,
                        v_init=CAS_current,
                        v_final=CAS_final,
                        speed_step=speed.stepSpeed,
                        Hp_init=Hp_current,
                        control=controlTarget,
                        phase="Climb",
                        config=config_current,
                        m_init=massCurrent,
                        wS=meteo.wS,
                        deltaTemp=meteo.deltaTemp,
                        calculationType=calculationType,
                    )
                    trajectory.append(AC, flightTrajectory)
                    CAS_current, Hp_current, massCurrent, config_current = (
                        trajectory.getFinalValue(
                            AC, ["CAS", "Hp", "mass", "config"]
                        )
                    )
                    if (
                        calculationType == CalculationType.POINT
                        and Hp_current
                        != pressureAltitude.finalPressureAltitude
                    ):
                        trajectory.removeLines(AC, numberOfLines=1)

        elif Hp_next == crossoverAltitude:
            # check if the next calculated speed is smaller then my current speed, to avoid deceleration during climb
            # if CAS_final < CAS_current:
            # CAS_final = CAS_current

            flightTrajectory = trajectorySegments.constantSpeedRating(
                AC=AC,
                speedType=SpeedType.CAS,
                v=CAS_final,
                Hp_init=Hp_current,
                Hp_final=Hp_next,
                Hp_step=stepPressureAltitude,
                m_init=massCurrent,
                wS=meteo.wS,
                deltaTemp=meteo.deltaTemp,
                reducedPower=reducedPower,
                calculationType=calculationType,
            )

            trajectory.append(AC, flightTrajectory)
            CAS_current, Hp_current, massCurrent, M_current, config_current = (
                trajectory.getFinalValue(
                    AC, ["CAS", "Hp", "mass", "M", "config"]
                )
            )
            if (
                calculationType == CalculationType.POINT
                and Hp_current != pressureAltitude.finalPressureAltitude
            ):
                trajectory.removeLines(AC, numberOfLines=1)

            # check if the calculation stopped before the final altitude goal
            if Hp_current < Hp_next:
                break
        else:
            flightTrajectory = trajectorySegments.constantSpeedRating(
                AC=AC,
                speedType=SpeedType.M,
                v=M_current,
                Hp_init=Hp_current,
                Hp_final=Hp_next,
                Hp_step=stepPressureAltitude,
                m_init=massCurrent,
                wS=meteo.wS,
                deltaTemp=meteo.deltaTemp,
                reducedPower=reducedPower,
                calculationType=calculationType,
            )

            trajectory.append(AC, flightTrajectory)
            CAS_current, Hp_current, massCurrent, config_current = (
                trajectory.getFinalValue(AC, ["CAS", "Hp", "mass", "config"])
            )
            if (
                calculationType == CalculationType.POINT
                and Hp_current != pressureAltitude.finalPressureAltitude
            ):
                trajectory.removeLines(AC, numberOfLines=1)

            # check if the calculation stopped before the final altitude goal
            if Hp_current < Hp_next:
                break

    return trajectory


def apcDescentCasMach(
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo,
    controlTarget: myTypes.ControlTarget,
    speedBrakes: myTypes.SpeedBrakes,
    casMachSpeedSchedule: myTypes.CASMACHSpeedSchedule,
    arrivalProfile: myTypes.ArrivalProfile,
) -> FT:
    """Calculates a complete descent trajectory using a CAS/Mach speed schedule.

    This function simulates a multi-segment descent profile. It manages the
    transition from constant Mach at high altitudes to constant CAS below the
    crossover altitude. It accounts for speed restrictions (e.g., FL100),
    engine-specific threshold altitudes, and various acceleration/deceleration
    logics (AT, BEFORE, AFTER) as defined in the speed configuration.

    :param AC: Aircraft model instance containing performance limits and engine data.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing initial, final, and step values.
    :param speed: Speed configuration including deceleration level types and steps.
    :param mass: Current aircraft mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param controlTarget: Target parameters governing the descent segments.
    :param speedBrakes: Configuration for speed brake state and deployment value.
    :param casMachSpeedSchedule: Object defining target Mach and CAS values for the descent.
    :param arrivalProfile: Profile type definition (e.g., expedite or constant CAS).

    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo
    :type controlTarget: myTypes.ControlTarget
    :type speedBrakes: myTypes.SpeedBrakes
    :type casMachSpeedSchedule: myTypes.CASMACHSpeedSchedule
    :type arrivalProfile: myTypes.ArrivalProfile

    :returns: A Flight Trajectory object containing the concatenated descent segments.
    :rtype: FT
    """

    trajectory = FT()

    mass_current = mass

    # setting the current aero configuration
    config_current = None

    # determine the BADA speed schedule [m/s and Mach]
    CASbelowFL100 = casMachSpeedSchedule.CASbelowFL100
    CASaboveFL100 = casMachSpeedSchedule.CASaboveFL100
    Mach = casMachSpeedSchedule.Mach

    if CASbelowFL100 is None and CASaboveFL100 is None and Mach is None:
        # if there is no definition of speed schedule, take the default one from BADA
        [Vdes1, Vdes2, Mdes] = AC.flightEnvelope.getSpeedSchedule(
            phase="Descent"
        )  # BADA Descent speed schedule
        speedSchedule = AC.flightEnvelope.getSpeedSchedule(phase="Descent")
    else:
        Vdes1 = conv.kt2ms(CASbelowFL100)
        Vdes2 = conv.kt2ms(CASaboveFL100)
        Mdes = Mach
        speedSchedule = [Vdes1, Vdes2, Mdes]

    # arrival profile
    if arrivalProfile.arrivalProfileType == ArrivalProfileType.expedite:
        expedite = True
    else:
        expedite = False

    # calculate crossover altitude [ft]
    crossoverAltitude = conv.m2ft(atm.crossOver(Vdes2, Mdes))
    # calculate tropopause altitude [ft]
    tropopauseAltitude = conv.m2ft(const.h_11)

    # Collect special altitudes that need to be included.
    specialAltitudes = [crossoverAltitude, tropopauseAltitude]
    configChangesGPFAltitudes = [3000, 8000]

    if (
        arrivalProfile.arrivalProfileType
        == ArrivalProfileType.constCASbelow100
    ):
        procedeureAltitudeList = [10000]
    else:
        if AC.engineType == "JET" or AC.engineType == "TURBOPROP":
            procedeureAltitudeList = [1000, 1500, 2000, 3000, 6000, 10000]
        elif AC.engineType == "PISTON" or AC.engineType == "ELECTRIC":
            procedeureAltitudeList = [500, 1000, 1500, 10000]

    # Merge the arrays, remove duplicates, and sort in descending order.
    altitudeDescentThresholdList = np.sort(
        np.unique(
            np.concatenate(
                (
                    procedeureAltitudeList,
                    configChangesGPFAltitudes,
                    specialAltitudes,
                    [pressureAltitude.finalPressureAltitude],
                )
            )
        )
    )[::-1]

    # Now filter the list: remove altitudes below initPressureAltitude
    # and above finalPressureAltitude.
    altitudeDescentThresholdList = altitudeDescentThresholdList[
        (altitudeDescentThresholdList <= pressureAltitude.initPressureAltitude)
        & (
            altitudeDescentThresholdList
            >= pressureAltitude.finalPressureAltitude
        )
    ]

    # set current altitude
    Hp_current = pressureAltitude.initPressureAltitude

    # set current speed
    [theta, delta, sigma] = atm.atmosphereProperties(
        h=conv.ft2m(Hp_current), deltaTemp=meteo.deltaTemp
    )

    CAS_current = conv.ms2kt(
        AC.ARPM.descentSpeed(
            h=conv.ft2m(Hp_current),
            mass=mass_current,
            theta=theta,
            delta=delta,
            deltaTemp=meteo.deltaTemp,
            speedSchedule_default=speedSchedule,
        )[0]
    )
    M_current = atm.cas2Mach(
        cas=conv.kt2ms(CAS_current), theta=theta, delta=delta, sigma=sigma
    )

    if (
        arrivalProfile.arrivalProfileType
        == ArrivalProfileType.constCASbelow100
    ):
        CAS_current = conv.ms2kt(Vdes1)

    for Hp_next in altitudeDescentThresholdList:
        [theta, delta, sigma] = atm.atmosphereProperties(
            h=conv.ft2m(Hp_next), deltaTemp=meteo.deltaTemp
        )
        CAS_final = conv.ms2kt(
            AC.ARPM.descentSpeed(
                h=conv.ft2m(Hp_next - 1),
                mass=mass_current,
                theta=theta,
                delta=delta,
                deltaTemp=meteo.deltaTemp,
                speedSchedule_default=speedSchedule,
            )[0]
        )

        if Hp_current > crossoverAltitude:
            # descent to final altitude - no deceleration, because flying constant M
            flightTrajectory = trajectorySegments.constantSpeedRating(
                AC=AC,
                speedType="M",
                v=M_current,
                Hp_init=Hp_current,
                Hp_final=Hp_next,
                Hp_step=pressureAltitude.stepPressureAltitude,
                m_init=mass_current,
                wS=meteo.wS,
                deltaTemp=meteo.deltaTemp,
                expedite=expedite,
                calculationType=calculationType,
            )
            trajectory.append(AC, flightTrajectory)
            (CAS_current, Hp_current, mass_current, config_current) = (
                trajectory.getFinalValue(AC, ["CAS", "Hp", "mass", "config"])
            )
            if (
                calculationType == "POINT"
                and Hp_current != pressureAltitude.finalPressureAltitude
            ):
                trajectory.removeLines(AC, numberOfLines=1)
        else:
            if (
                arrivalProfile == "constCASbelow100"
                or arrivalProfile == "expedite"
            ):
                CAS_final = conv.ms2kt(Vdes1)

            if speed.accelerationLevelKind == AccelerationLevelKind.AT:
                flightTrajectory = trajectorySegments.constantSpeedRating(
                    AC=AC,
                    speedType="CAS",
                    v=CAS_current,
                    Hp_init=Hp_current,
                    Hp_final=Hp_next,
                    Hp_step=pressureAltitude.stepPressureAltitude,
                    m_init=mass_current,
                    wS=meteo.wS,
                    deltaTemp=meteo.deltaTemp,
                    expedite=expedite,
                    calculationType=calculationType,
                )
                trajectory.append(AC, flightTrajectory)
                (
                    CAS_current,
                    Hp_current,
                    mass_current,
                    config_current,
                ) = trajectory.getFinalValue(
                    AC, ["CAS", "Hp", "mass", "config"]
                )
                if (
                    calculationType == "POINT"
                    and Hp_current != pressureAltitude.finalPressureAltitude
                ):
                    trajectory.removeLines(AC, numberOfLines=1)

                if (
                    abs(CAS_current - CAS_final) > 0.3
                    and Hp_next != pressureAltitude.finalPressureAltitude
                ):
                    flightTrajectory = trajectorySegments.accDec(
                        AC=AC,
                        speedType="CAS",
                        v_init=CAS_current,
                        v_final=CAS_final,
                        speed_step=speed.stepSpeed,
                        Hp_init=Hp_current,
                        control=controlTarget,
                        phase="Cruise",
                        m_init=mass_current,
                        wS=meteo.wS,
                        deltaTemp=meteo.deltaTemp,
                        config=config_current,
                        calculationType=calculationType,
                    )
                    trajectory.append(AC, flightTrajectory)
                    (
                        CAS_current,
                        Hp_current,
                        mass_current,
                        config_current,
                    ) = trajectory.getFinalValue(
                        AC, ["CAS", "Hp", "mass", "config"]
                    )
                    if (
                        calculationType == "POINT"
                        and Hp_current
                        != pressureAltitude.finalPressureAltitude
                    ):
                        trajectory.removeLines(AC, numberOfLines=1)

            elif speed.accelerationLevelKind == AccelerationLevelKind.AFTER:
                flightTrajectory = trajectorySegments.constantSpeedRating(
                    AC=AC,
                    speedType="CAS",
                    v=CAS_current,
                    Hp_init=Hp_current,
                    Hp_final=Hp_next,
                    Hp_step=pressureAltitude.stepPressureAltitude,
                    m_init=mass_current,
                    wS=meteo.wS,
                    deltaTemp=meteo.deltaTemp,
                    expedite=expedite,
                    calculationType=calculationType,
                )
                trajectory.append(AC, flightTrajectory)
                (CAS_current, Hp_current, mass_current, config_current) = (
                    trajectory.getFinalValue(
                        AC, ["CAS", "Hp", "mass", "config"]
                    )
                )
                if (
                    calculationType == "POINT"
                    and Hp_current != pressureAltitude.finalPressureAltitude
                ):
                    trajectory.removeLines(AC, numberOfLines=1)

                if (
                    abs(CAS_current - CAS_final) > 0.3
                    and Hp_next != pressureAltitude.finalPressureAltitude
                ):
                    flightTrajectory = trajectorySegments.accDec(
                        AC=AC,
                        speedType="CAS",
                        v_init=CAS_current,
                        v_final=CAS_final,
                        speed_step=speed.stepSpeed,
                        Hp_init=Hp_current,
                        control=controlTarget,
                        phase="Descent",
                        m_init=mass_current,
                        wS=meteo.wS,
                        deltaTemp=meteo.deltaTemp,
                        config=config_current,
                        calculationType=calculationType,
                    )
                    trajectory.append(AC, flightTrajectory)
                    (
                        CAS_current,
                        Hp_current,
                        mass_current,
                        config_current,
                    ) = trajectory.getFinalValue(
                        AC, ["CAS", "Hp", "mass", "config"]
                    )
                    if (
                        calculationType == "POINT"
                        and Hp_current
                        != pressureAltitude.finalPressureAltitude
                    ):
                        trajectory.removeLines(AC, numberOfLines=1)

            elif speed.accelerationLevelKind == AccelerationLevelKind.BEFORE:
                Hp_final = Hp_next
                # determine how much altitude is needed to decelerate, than descent to that altitude and decelerate.
                if (
                    abs(CAS_current - CAS_final) > 0.3
                    and Hp_next != pressureAltitude.finalPressureAltitude
                ):
                    Hp_final = Hp_current
                    diff = Hp_next - Hp_current
                    while abs(diff) > 1:
                        tempTrajectory = FT()
                        flightTrajectory = trajectorySegments.accDec(
                            AC=AC,
                            speedType="CAS",
                            v_init=CAS_current,
                            v_final=CAS_final,
                            speed_step=speed.stepSpeed,
                            Hp_init=Hp_final,
                            control=controlTarget,
                            phase="Descent",
                            config=config_current,
                            m_init=mass_current,
                            wS=meteo.wS,
                            deltaTemp=meteo.deltaTemp,
                            calculationType=calculationType,
                        )
                        tempTrajectory.append(AC, flightTrajectory)
                        CAS, Hp_end = tempTrajectory.getFinalValue(
                            AC, ["CAS", "Hp"]
                        )
                        diff = Hp_next - Hp_end
                        Hp_final = Hp_final + diff

                # descent to final altitude before deceleration
                flightTrajectory = trajectorySegments.constantSpeedRating(
                    AC=AC,
                    speedType="CAS",
                    v=CAS_current,
                    Hp_init=Hp_current,
                    Hp_final=Hp_final,
                    Hp_step=pressureAltitude.stepPressureAltitude,
                    m_init=mass_current,
                    wS=meteo.wS,
                    deltaTemp=meteo.deltaTemp,
                    expedite=expedite,
                    calculationType=calculationType,
                )
                trajectory.append(AC, flightTrajectory)
                (CAS_current, Hp_current, mass_current, config_current) = (
                    trajectory.getFinalValue(
                        AC, ["CAS", "Hp", "mass", "config"]
                    )
                )
                if (
                    calculationType == "POINT"
                    and Hp_current != pressureAltitude.finalPressureAltitude
                ):
                    trajectory.removeLines(AC, numberOfLines=1)

                # decelerate
                if (
                    abs(CAS_current - CAS_final) > 0.3
                    and Hp_next != pressureAltitude.finalPressureAltitude
                ):
                    flightTrajectory = trajectorySegments.accDec(
                        AC=AC,
                        speedType="CAS",
                        v_init=CAS_current,
                        v_final=CAS_final,
                        speed_step=speed.stepSpeed,
                        Hp_init=Hp_current,
                        control=controlTarget,
                        phase="Descent",
                        m_init=mass_current,
                        wS=meteo.wS,
                        deltaTemp=meteo.deltaTemp,
                        config=config_current,
                        calculationType=calculationType,
                    )
                    trajectory.append(AC, flightTrajectory)
                    (CAS_current, Hp_current, mass_current, config_current) = (
                        trajectory.getFinalValue(
                            AC, ["CAS", "Hp", "mass", "config"]
                        )
                    )
                    if abs(Hp_current - Hp_next) < 1:
                        trajectory.overwriteLastValue(AC, "Hp", Hp_next)
                        Hp_current = Hp_next
                    if (
                        calculationType == "POINT"
                        and Hp_current
                        != pressureAltitude.finalPressureAltitude
                    ):
                        trajectory.removeLines(AC, numberOfLines=1)

    return trajectory


def hpcClimbARPM(
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo,
    controlTarget: myTypes.ControlTarget,
    rating: myTypes.HRating,
) -> FT:
    """Calculates a helicopter climb trajectory using the ARPM  model.

    :param AC: Aircraft (Helicopter) model instance with ARPM procedure methods.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing initial, final, and step values.
    :param speed: Speed configuration (Note: Initial TAS is derived from ARPM).
    :param mass: Current helicopter mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param controlTarget: Target parameters for the climb segment.
    :param rating: The helicopter power rating (e.g., 'ARPM', 'TOP', 'MCP').

    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo
    :type controlTarget: myTypes.ControlTarget
    :type rating: myTypes.HRating

    :returns: A Flight Trajectory object representing the helicopter climb profile.
    :rtype: FT
    """

    trajectory = FT()

    mass_current = mass
    Hp_current = pressureAltitude.initPressureAltitude

    #  create altitude climb list
    altitudeClimbThresholdList = [pressureAltitude.finalPressureAltitude]

    if (
        5 > pressureAltitude.initPressureAltitude
        and 5 < pressureAltitude.finalPressureAltitude
    ):
        altitudeClimbThresholdList.append(5)
        altitudeClimbThresholdList.append(
            pressureAltitude.initPressureAltitude
            + pressureAltitude.stepPressureAltitude
        )
        altitudeClimbThresholdList.sort()

    [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = AC.ARPM.ARPMProcedure(
        phase="Climb",
        h=conv.ft2m(Hp_current),
        deltaTemp=meteo.deltaTemp,
        mass=mass_current,
        rating=rating,
    )
    TAS_current = conv.ms2kt(tas)
    ROCD_current = conv.m2ft(ROCD) * 60

    Hp_speed = Hp_current  # Altitude to calculate the speed for the whole climb segment, which will not change at this point

    for Hp_next in altitudeClimbThresholdList:
        if Hp_next == 5:
            flightTrajectory = trajectorySegments.constantSpeedROCD(
                AC=AC,
                speedType="TAS",
                v=TAS_current,
                ROCDtarget=ROCD_current,
                Hp_init=Hp_current,
                Hp_final=5,
                Hp_step=pressureAltitude.stepPressureAltitude,
                m_init=mass_current,
                wS=meteo.wS,
                deltaTemp=meteo.deltaTemp,
                calculationType=calculationType,
            )
            trajectory.append(AC, flightTrajectory)
        else:
            if Hp_current == 5:
                Hp_speed = Hp_current + 0.1

                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = (
                    AC.ARPM.ARPMProcedure(
                        phase="Climb",
                        h=conv.ft2m(Hp_speed),
                        deltaTemp=meteo.deltaTemp,
                        mass=mass_current,
                        rating=rating,
                    )
                )
                TAS_current = conv.ms2kt(tas)
                ROCD_current = conv.m2ft(ROCD) * 60

            if rating == "ARPM":
                flightTrajectory = trajectorySegments.constantSpeedROCD(
                    AC=AC,
                    speedType="TAS",
                    v=TAS_current,
                    ROCDtarget=ROCD_current,
                    Hp_init=Hp_current,
                    Hp_final=Hp_next,
                    Hp_step=pressureAltitude.stepPressureAltitude,
                    m_init=mass_current,
                    wS=meteo.wS,
                    deltaTemp=meteo.deltaTemp,
                    calculationType=calculationType,
                )
                trajectory.append(AC, flightTrajectory)
            else:
                flightTrajectory = trajectorySegments.constantSpeedRating(
                    AC=AC,
                    speedType="TAS",
                    v=TAS_current,
                    Hp_init=Hp_current,
                    Hp_final=Hp_next,
                    Hp_step=pressureAltitude.stepPressureAltitude,
                    m_init=mass_current,
                    wS=meteo.wS,
                    deltaTemp=meteo.deltaTemp,
                    calculationType=calculationType,
                    initRating=rating,
                )
                trajectory.append(AC, flightTrajectory)

        CAS_current, TAS_current, Hp_current, mass_current = (
            trajectory.getFinalValue(AC, ["CAS", "TAS", "Hp", "mass"])
        )

        if (
            calculationType == CalculationType.POINT
            and Hp_current != pressureAltitude.finalPressureAltitude
        ):
            trajectory.removeLines(AC, numberOfLines=1)

    return trajectory


def hpcDescentARPM(
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo,
    controlTarget: myTypes.ControlTarget,
) -> FT:
    """Calculates a helicopter descent and landing trajectory using the ARPM model.

    This function simulates a rotorcraft descent profile by identifying key
    altitude thresholds (500ft, 150ft, and 5ft) that trigger updates to the
    ARPM performance parameters. It handles the transition from forward
    descent to the final landing phase, automatically calculating the required
    TAS and ROCD based on the ARPM 'Descent' phase logic.

    :param AC: Aircraft (Helicopter) model instance with ARPM procedure methods.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing initial, final, and step values.
    :param speed: Speed configuration (Note: TAS is derived from ARPM).
    :param mass: Current helicopter mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param controlTarget: Target parameters for the descent segment.

    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo
    :type controlTarget: myTypes.ControlTarget

    :returns: A Flight Trajectory object representing the helicopter descent profile.
    :rtype: FT
    """

    trajectory = FT()

    mass_current = mass
    Hp_current = pressureAltitude.initPressureAltitude

    altitudeDescentThresholdList = [pressureAltitude.finalPressureAltitude]

    if (
        500 > pressureAltitude.finalPressureAltitude
        and 500 < pressureAltitude.initPressureAltitude
    ):
        altitudeDescentThresholdList.append(500)
    if (
        150 > pressureAltitude.finalPressureAltitude
        and 150 < pressureAltitude.initPressureAltitude
    ):
        altitudeDescentThresholdList.append(150)
    if (
        5 > pressureAltitude.finalPressureAltitude
        and 5 < pressureAltitude.initPressureAltitude
    ):
        altitudeDescentThresholdList.append(5)

    altitudeDescentThresholdList.sort(reverse=True)

    [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = AC.ARPM.ARPMProcedure(
        phase="Descent",
        h=conv.ft2m(Hp_current),
        deltaTemp=meteo.deltaTemp,
        mass=mass_current,
        rating="ARPM",
    )
    TAS_current = conv.ms2kt(tas)
    ROCD_current = conv.m2ft(ROCD) * 60

    Hp_speed = Hp_current  # Altitude to calculate the speed for the whole climb segment, which will not change at this point
    for Hp_next in altitudeDescentThresholdList:
        if Hp_next == 500:
            flightTrajectory = trajectorySegments.constantSpeedROCD(
                AC=AC,
                speedType=SpeedType.TAS,
                v=TAS_current,
                ROCDtarget=ROCD_current,
                Hp_init=Hp_current,
                Hp_final=Hp_next,
                Hp_step=pressureAltitude.stepPressureAltitude,
                m_init=mass_current,
                wS=meteo.wS,
                deltaTemp=meteo.deltaTemp,
                calculationType=calculationType,
            )
            trajectory.append(AC, flightTrajectory)

        elif Hp_next == 150 or Hp_next == 5:
            Hp_speed = Hp_current + 0.1
            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = (
                AC.ARPM.ARPMProcedure(
                    phase="Descent",
                    h=conv.ft2m(Hp_next),
                    deltaTemp=meteo.deltaTemp,
                    mass=mass_current,
                    rating="ARPM",
                )
            )
            ROCD_current = conv.m2ft(ROCD) * 60

            flightTrajectory = trajectorySegments.constantSpeedROCD(
                AC=AC,
                speedType=SpeedType.TAS,
                v=TAS_current,
                ROCDtarget=ROCD_current,
                Hp_init=Hp_current,
                Hp_final=Hp_next,
                Hp_step=pressureAltitude.stepPressureAltitude,
                m_init=mass_current,
                wS=meteo.wS,
                deltaTemp=meteo.deltaTemp,
                calculationType=calculationType,
            )
            trajectory.append(AC, flightTrajectory)

        else:  # vertical landing
            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = (
                AC.ARPM.ARPMProcedure(
                    phase="Descent",
                    h=conv.ft2m(0),
                    deltaTemp=meteo.deltaTemp,
                    mass=mass_current,
                    rating="ARPM",
                )
            )
            TAS_current = conv.ms2kt(tas)
            ROCD_current = conv.m2ft(ROCD) * 60

            flightTrajectory = trajectorySegments.constantSpeedROCD(
                AC=AC,
                speedType=SpeedType.TAS,
                v=TAS_current,
                ROCDtarget=ROCD_current,
                Hp_init=Hp_current,
                Hp_final=Hp_next,
                Hp_step=pressureAltitude.stepPressureAltitude,
                m_init=mass_current,
                wS=meteo.wS,
                deltaTemp=meteo.deltaTemp,
                calculationType=calculationType,
            )
            trajectory.append(AC, flightTrajectory)

        CAS_current, TAS_current, Hp_current, mass_current = (
            trajectory.getFinalValue(AC, ["CAS", "TAS", "Hp", "mass"])
        )
        if (
            calculationType == CalculationType.POINT
            and Hp_current != pressureAltitude.finalPressureAltitude
        ):
            trajectory.removeLines(AC, numberOfLines=1)

    return trajectory


def apcClimbCalculation(
    climbType: myTypes.ClimbType,
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo | None = None,
    controlTarget: myTypes.ControlTarget | None = None,
    speedBrakes: myTypes.SpeedBrakes | None = None,
    climbCASMACHProfileConfiguration: myTypes.ClimbCASMACHProfileConfiguration
    | None = None,
) -> FT:
    """Dispatches and calculates the aircraft climb trajectory based on the specified climb type.

    This function serves as a high-level interface that selects the appropriate
    climb logic—be it a constant rate of climb, a speed transition (acceleration),
    or a standard CAS/Mach speed schedule. It ensures the correct parameters
    and configurations are passed down to the specialized performance sub-functions.

    :param climbType: The strategy for the climb (RATE, ACCELERATION, or CASMACH).
    :param AC: Aircraft model instance containing performance and engine data.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing initial, final, and step values.
    :param speed: Speed configuration including initial, final, and step speed values.
    :param mass: Current aircraft mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param controlTarget: Target parameters governing the climb (e.g., target ROCD).
    :param speedBrakes: Configuration for speed brake state and deployment value.
    :param climbCASMACHProfileConfiguration: Detailed settings for CAS/Mach climbs, including takeoff procedures.

    :type climbType: myTypes.ClimbType
    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo | None
    :type controlTarget: myTypes.ControlTarget | None
    :type speedBrakes: myTypes.SpeedBrakes | None
    :type climbCASMACHProfileConfiguration: myTypes.ClimbCASMACHProfileConfiguration | None

    :returns: A Flight Trajectory object containing the results of the specific climb calculation.
    :rtype: FT
    """

    trajectory = FT()

    match climbType:
        case ClimbType.RATE:
            trajectory = climbDescentRate(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                speed=speed,
                mass=mass,
                meteo=meteo,
                controlTarget=controlTarget,
                speedBrakes=speedBrakes,
            )
        case ClimbType.ACCELERATION:
            trajectory = climbDescentAccelerationDeceleration(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                flightPhase=FlightPhase.CLIMB,
                speed=speed,
                mass=mass,
                meteo=meteo,
                controlTarget=controlTarget,
                speedBrakes=speedBrakes,
            )
        case ClimbType.CASMACH:
            trajectory = apcClimbCasMach(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                speed=speed,
                mass=mass,
                meteo=meteo,
                controlTarget=controlTarget,
                speedBrakes=speedBrakes,
                casMachSpeedSchedule=climbCASMACHProfileConfiguration.casMachSpeedSchedule,
                takeOffProcedure=climbCASMACHProfileConfiguration.takeOffProcedure,
                departureProfile=climbCASMACHProfileConfiguration.departureProfile,
                reducedPower=climbCASMACHProfileConfiguration.reducedPower,
            )
    return trajectory


def hpcClimbCalculation(
    climbType: myTypes.ClimbType,
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo | None = None,
    controlTarget: myTypes.ControlTarget | None = None,
    hClimbRatingConfiguration: myTypes.HClimbRatingConfiguration | None = None,
) -> FT:
    """Dispatches and calculates the helicopter climb trajectory based on the specified climb type.

    :param climbType: The strategy for the climb (RATE, ACCELERATION, or ARPM).
    :param AC: Aircraft (Helicopter) model instance with ARPM performance data.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing initial, final, and step values.
    :param speed: Speed configuration including initial, final, and step speed values.
    :param mass: Current helicopter mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param controlTarget: Target parameters governing the climb (e.g., target ROCD).
    :param hClimbRatingConfiguration: Configuration for helicopter-specific power ratings.

    :type climbType: myTypes.ClimbType
    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo | None
    :type controlTarget: myTypes.ControlTarget | None
    :type hClimbRatingConfiguration: myTypes.HClimbRatingConfiguration | None

    :returns: A Flight Trajectory object containing the results of the helicopter climb.
    :rtype: FT
    """

    trajectory = FT()

    match climbType:
        case HClimbType.RATE:
            trajectory = climbDescentRate(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                speed=speed,
                mass=mass,
                meteo=meteo,
                controlTarget=controlTarget,
            )
        case HClimbType.ACCELERATION:
            trajectory = climbDescentAccelerationDeceleration(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                flightPhase=FlightPhase.CLIMB,
                speed=speed,
                mass=mass,
                meteo=meteo,
                controlTarget=controlTarget,
            )
        case HClimbType.ARPM:
            trajectory = hpcClimbARPM(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                speed=speed,
                mass=mass,
                meteo=meteo,
                controlTarget=controlTarget,
                rating=hClimbRatingConfiguration.rating,
            )
    return trajectory


def hpcDescentCalculation(
    descentType: myTypes.DescentType,
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo | None = None,
    controlTarget: myTypes.ControlTarget | None = None,
    speedBrakes: myTypes.SpeedBrakes | None = None,
    CASMACHProfileConfiguration: (
        myTypes.DescentCASMACHProfileConfiguration | None
    ) = None,
) -> FT:
    """Dispatches and calculates the helicopter descent trajectory based on the specified type.

    :param descentType: The strategy for the descent (RATE, SLOPE, EMERGENCY, ACCELERATION, or ARPM).
    :param AC: Aircraft (Helicopter) model instance containing VNE and ARPM data.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing initial, final, and step values.
    :param speed: Speed configuration including initial, final, and step speed values.
    :param mass: Current helicopter mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param controlTarget: Target parameters governing the descent (e.g., target ROCD or slope).
    :param speedBrakes: Configuration for speed brake state and deployment value.
    :param CASMACHProfileConfiguration: Detailed settings for fixed-wing style CAS/Mach descent schedules.

    :type descentType: myTypes.DescentType
    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo | None
    :type controlTarget: myTypes.ControlTarget | None
    :type speedBrakes: myTypes.SpeedBrakes | None
    :type CASMACHProfileConfiguration: myTypes.DescentCASMACHProfileConfiguration | None

    :returns: A Flight Trajectory object containing the computed helicopter descent results.
    :rtype: FT
    """

    trajectory = FT()

    match descentType:
        case HDescentType.RATE:
            trajectory = climbDescentRate(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                speed=speed,
                mass=mass,
                meteo=meteo,
                controlTarget=controlTarget,
            )
        case HDescentType.SLOPE:
            trajectory = climbDescentSlope(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                speed=speed,
                mass=mass,
                meteo=meteo,
                controlTarget=controlTarget,
            )
        case HDescentType.EMERGENCY:
            trajectory = hpcDescentEmergency(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                speed=speed,
                mass=mass,
                meteo=meteo,
                controlTarget=controlTarget,
            )
        case HDescentType.ACCELERATION:
            trajectory = climbDescentAccelerationDeceleration(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                flightPhase=FlightPhase.DESCENT,
                speed=speed,
                mass=mass,
                meteo=meteo,
                controlTarget=controlTarget,
            )
        case HDescentType.ARPM:
            trajectory = hpcDescentARPM(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                speed=speed,
                mass=mass,
                meteo=meteo,
                controlTarget=controlTarget,
            )
    return trajectory


def apcDescentCalculation(
    descentType: myTypes.DescentType,
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo | None = None,
    controlTarget: myTypes.ControlTarget | None = None,
    speedBrakes: myTypes.SpeedBrakes | None = None,
    CASMACHProfileConfiguration: (
        myTypes.DescentCASMACHProfileConfiguration | None
    ) = None,
) -> FT:
    """Dispatches and calculates the aircraft descent trajectory based on the specified type.

    This function serves as a high-level interface that selects the appropriate
    descent logic—supporting constant rate (ROCD), fixed gradient (slope),
    high-speed emergency procedures, speed transitions (deceleration), or
    managed CAS/Mach speed schedules. It coordinates the inputs and routes
    them to the specialized performance sub-functions.

    :param descentType: The strategy for the descent (RATE, SLOPE, EMERGENCY, ACCELERATION, or CASMACH).
    :param AC: Aircraft model instance containing performance and engine data.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing initial, final, and step values.
    :param speed: Speed configuration including initial, final, and step speed values.
    :param mass: Current aircraft mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param controlTarget: Target parameters governing the descent (e.g., target ROCD or slope).
    :param speedBrakes: Configuration for speed brake state and deployment value.
    :param CASMACHProfileConfiguration: Detailed settings for managed CAS/Mach descent profiles.

    :type descentType: myTypes.DescentType
    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo | None
    :type controlTarget: myTypes.ControlTarget | None
    :type speedBrakes: myTypes.SpeedBrakes | None
    :type CASMACHProfileConfiguration: myTypes.DescentCASMACHProfileConfiguration | None

    :returns: A Flight Trajectory object containing the results of the specific descent calculation.
    :rtype: FT
    """

    trajectory = FT()

    match descentType:
        case DescentType.RATE:
            trajectory = climbDescentRate(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                speed=speed,
                mass=mass,
                meteo=meteo,
                controlTarget=controlTarget,
                speedBrakes=speedBrakes,
            )
        case DescentType.SLOPE:
            trajectory = climbDescentSlope(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                speed=speed,
                mass=mass,
                meteo=meteo,
                controlTarget=controlTarget,
                speedBrakes=speedBrakes,
            )
        case DescentType.EMERGENCY:
            trajectory = apcDescentEmergency(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                speed=speed,
                mass=mass,
                meteo=meteo,
                controlTarget=controlTarget,
                speedBrakes=speedBrakes,
            )
        case DescentType.ACCELERATION:
            trajectory = climbDescentAccelerationDeceleration(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                flightPhase=FlightPhase.DESCENT,
                speed=speed,
                mass=mass,
                meteo=meteo,
                controlTarget=controlTarget,
                speedBrakes=speedBrakes,
            )
        case DescentType.CASMACH:
            trajectory = apcDescentCasMach(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                speed=speed,
                mass=mass,
                meteo=meteo,
                controlTarget=controlTarget,
                speedBrakes=speedBrakes,
                casMachSpeedSchedule=CASMACHProfileConfiguration.casMachSpeedSchedule,
                arrivalProfile=CASMACHProfileConfiguration.arrivalProfile,
            )
    return trajectory


def hpcLevelConstantSpeed(
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo,
    integrationType: myTypes.IntegrationType,
    cruiseLength: float,
    stepSize: float,
) -> FT:
    """Calculates a level flight trajectory at a constant speed for helicopters.

    This function simulates cruise performance at a fixed altitude. It can
    calculate a single state (POINT) across a range of altitudes or integrate
    fuel burn and time over a specific distance or duration (INTEGRATED). It
    supports automated speed selection for MRC (Maximum Range Cruise),
    LRC (Long Range Cruise), and MEC (Maximum Endurance Cruise).

    :param AC: Aircraft (Helicopter) model instance including optimization modules.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing initial, final, and step values.
    :param speed: Speed configuration; can be a fixed value or an optimization type.
    :param mass: Current helicopter mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param integrationType: Defines if cruiseLength is DISTANCE or TIME.
    :param cruiseLength: The total length of the cruise segment [NM or s].
    :param stepSize: The discretization step for the integrated calculation.

    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo
    :type integrationType: myTypes.IntegrationType
    :type cruiseLength: float
    :type stepSize: float

    :returns: A Flight Trajectory object containing the level flight results.
    :rtype: FT
    """

    trajectory = FT()

    if (
        speed.speedType == CruiseSpeedType.MRC
        or speed.speedType == CruiseSpeedType.LRC
        or speed.speedType == CruiseSpeedType.MEC
    ):
        speedType = SpeedType.TAS
        cruiseSpeed = AC.OPT.getOPTParam(
            speed.speedType,
            pressureAltitude.initPressureAltitude,
            mass,
            meteo.deltaTemp,
        )
    else:
        speedType = speed.speedType
        cruiseSpeed = speed.initSpeed

    if calculationType == CalculationType.POINT:
        Hp_range = np.arange(
            pressureAltitude.initPressureAltitude,
            pressureAltitude.finalPressureAltitude,
            pressureAltitude.stepPressureAltitude,
        )
        Hp_range = np.append(Hp_range, pressureAltitude.finalPressureAltitude)

        for Hp in Hp_range:
            flightTrajectory = trajectorySegments.constantSpeedLevel(
                AC=AC,
                lengthType=IntegrationType.DISTANCE,
                length=1,
                speedType=speedType,
                v=cruiseSpeed,
                Hp_init=Hp,
                m_init=mass,
                wS=meteo.wS,
                deltaTemp=meteo.deltaTemp,
                calculationType=calculationType,
            )
            trajectory.append(AC, flightTrajectory)
            trajectory.removeLines(AC=AC, numberOfLines=1)

    if calculationType == CalculationType.INTEGRATED:
        flightTrajectory = trajectorySegments.constantSpeedLevel(
            AC=AC,
            lengthType=integrationType,
            length=cruiseLength,
            speedType=speedType,
            v=cruiseSpeed,
            Hp_init=pressureAltitude.initPressureAltitude,
            maxRFL=pressureAltitude.finalPressureAltitude,
            HpStep=pressureAltitude.stepPressureAltitude,
            m_init=mass,
            wS=meteo.wS,
            deltaTemp=meteo.deltaTemp,
            calculationType=calculationType,
            step_length=stepSize,
        )
        trajectory.append(AC, flightTrajectory)

    return trajectory


def apcLevelConstantSpeed(
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo,
    integrationType: myTypes.IntegrationType,
    cruiseLength: float,
    stepSize: float,
    costIndex: float,
    stepClimb: bool,
) -> FT:
    """Calculates a level flight trajectory at a constant speed for fixed-wing aircraft.

    This function simulates cruise performance, supporting both point-based
    performance snapshots across an altitude range and integrated flight segments.
    It features advanced speed selection logic, including Maximum Range Cruise (MRC),
    Long Range Cruise (LRC), and Economic (ECON) speeds based on Cost Index.
    Integrated calculations can optionally include step-climb logic to optimize
    flight levels as mass decreases.

    :param AC: Aircraft model instance including performance optimization modules.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing initial, final, and step values.
    :param speed: Speed configuration; supports fixed values or optimized types (MRC, LRC, ECON).
    :param mass: Current aircraft mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param integrationType: Defines if cruiseLength is DISTANCE or TIME.
    :param cruiseLength: The total length of the cruise segment [NM or s].
    :param stepSize: The discretization step for the integrated calculation.
    :param costIndex: The ratio of time-related costs to fuel costs used for ECON speed.
    :param stepClimb: Boolean flag to enable or disable automatic step-climb logic.

    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo
    :type integrationType: myTypes.IntegrationType
    :type cruiseLength: float
    :type stepSize: float
    :type costIndex: float
    :type stepClimb: bool

    :returns: A Flight Trajectory object containing the cruise results.
    :rtype: FT
    """

    trajectory = FT()

    [theta, delta, sigma] = atm.atmosphereProperties(
        h=conv.ft2m(pressureAltitude.initPressureAltitude),
        deltaTemp=meteo.deltaTemp,
    )

    if (
        speed.speedType == CruiseSpeedType.MRC
        or speed.speedType == CruiseSpeedType.LRC
        or speed.speedType == CruiseSpeedType.MEC
    ):
        CW = AC.OPT.CW(mass=mass, delta=delta)
        cruiseSpeed = AC.OPT.getOPTParam(speed.speedType, CW)
        speedType = SpeedType.M
    elif speed.speedType == CruiseSpeedType.ECON:
        CW = AC.OPT.CW(mass=mass, delta=delta)
        CCI = AC.OPT.CCI(theta=theta, delta=delta, cI=costIndex)
        cruiseSpeed = AC.OPT.getOPTParam("ECON", CW, CCI)
        speedType = "M"
    else:
        cruiseSpeed = speed.initSpeed
        speedType = speed.speedType

    if calculationType == CalculationType.POINT:
        Hp_range = np.arange(
            pressureAltitude.initPressureAltitude,
            pressureAltitude.finalPressureAltitude,
            pressureAltitude.stepPressureAltitude,
        )
        Hp_range = np.append(Hp_range, pressureAltitude.finalPressureAltitude)

        for Hp in Hp_range:
            flightTrajectory = trajectorySegments.constantSpeedLevel(
                AC=AC,
                lengthType=IntegrationType.DISTANCE,
                length=1,
                speedType=speedType,
                v=cruiseSpeed,
                Hp_init=Hp,
                m_init=mass,
                wS=meteo.wS,
                deltaTemp=meteo.deltaTemp,
                calculationType=calculationType,
            )
            trajectory.append(AC, flightTrajectory)
            trajectory.removeLines(AC=AC, numberOfLines=1)

    if calculationType == CalculationType.INTEGRATED:
        flightTrajectory = trajectorySegments.constantSpeedLevel(
            AC=AC,
            lengthType=integrationType,
            length=cruiseLength,
            speedType=speedType,
            v=cruiseSpeed,
            Hp_init=pressureAltitude.initPressureAltitude,
            maxRFL=pressureAltitude.finalPressureAltitude,
            HpStep=pressureAltitude.stepPressureAltitude,
            m_init=mass,
            wS=meteo.wS,
            deltaTemp=meteo.deltaTemp,
            calculationType=calculationType,
            step_length=stepSize,
            stepClimb=stepClimb,
        )
        trajectory.append(AC, flightTrajectory)

    return trajectory


def apcLevelAccelerationDeceleration(
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo,
    controlTarget: myTypes.ControlTarget | None = None,
    speedBrakes: myTypes.SpeedBrakes | None = None,
) -> FT:
    """Calculates a level flight trajectory for speed transitions (acceleration or deceleration).

    :param AC: Aircraft model instance containing performance and aerodynamic data.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing the initial pressure altitude.
    :param speed: Speed configuration including initial, final, and step speed values.
    :param mass: Current aircraft mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param controlTarget: Target parameters governing the speed change (e.g., longitudinal acceleration).
    :param speedBrakes: Configuration for speed brake state and deployment value.

    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo
    :type controlTarget: myTypes.ControlTarget | None
    :type speedBrakes: myTypes.SpeedBrakes | None

    :returns: A Flight Trajectory object containing the acceleration/deceleration segment.
    :rtype: FT
    """

    trajectory = FT()

    if speedBrakes is None:
        speedBrakes = myTypes.SpeedBrakes(deployed=0, value=0)
    else:
        if speedBrakes.value is None:
            speedBrakesValue = 0.0
        else:
            speedBrakesValue = abs(min(speedBrakes.value, 100) * (0.03 / 100))

        speedBrakes = {
            "deployed": speedBrakes.deployed,
            "value": speedBrakesValue,
        }

    flightTrajectory = trajectorySegments.accDec(
        AC=AC,
        speedType=speed.speedType,
        v_init=speed.initSpeed,
        v_final=speed.finalSpeed,
        speed_step=speed.stepSpeed,
        phase="Cruise",
        Hp_init=pressureAltitude.initPressureAltitude,
        m_init=mass,
        deltaTemp=meteo.deltaTemp,
        wS=meteo.wS,
        controlTarget=controlTarget,
        calculationType=calculationType,
        speedBrakes=speedBrakes,
    )
    trajectory.append(AC, flightTrajectory)

    return trajectory


def hpcLevelAccelerationDeceleration(
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo,
    controlTarget: myTypes.ControlTarget | None = None,
) -> FT:
    """Calculates a level helicopter trajectory for speed transitions at a constant altitude.

    :param AC: Aircraft (Helicopter) model instance containing performance data.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing the initial pressure altitude.
    :param speed: Speed configuration including initial, final, and step speed values.
    :param mass: Current helicopter mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param controlTarget: Target parameters governing the transition (e.g., specific acceleration).

    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo
    :type controlTarget: myTypes.ControlTarget | None

    :returns: A Flight Trajectory object containing the computed transition segment.
    :rtype: FT
    """

    trajectory = FT()

    flightTrajectory = trajectorySegments.accDec(
        AC=AC,
        speedType=speed.speedType,
        v_init=speed.initSpeed,
        v_final=speed.finalSpeed,
        speed_step=speed.stepSpeed,
        phase="Cruise",
        Hp_init=pressureAltitude.initPressureAltitude,
        m_init=mass,
        deltaTemp=meteo.deltaTemp,
        wS=meteo.wS,
        controlTarget=controlTarget,
        calculationType=calculationType,
    )
    trajectory.append(AC, flightTrajectory)

    return trajectory


def hpcLevelCalculation(
    cruiseType: myTypes.CruiseType,
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo | None = None,
    controlTarget: myTypes.ControlTarget | None = None,
    integrationType: myTypes.IntegrationType | None = None,
    cruiseLength: float | None = None,
    stepSize: float | None = None,
) -> FT:
    """Dispatches and calculates the helicopter level flight trajectory based on cruise type.

    :param cruiseType: The strategy for level flight (CONSTANTSPEED or ACCELERATION).
    :param AC: Aircraft (Helicopter) model instance containing ARPM and OPT data.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing initial, final, and step values.
    :param speed: Speed configuration; fixed values or optimization types.
    :param mass: Current helicopter mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param controlTarget: Target parameters governing speed changes (for acceleration).
    :param integrationType: Defines if cruiseLength is DISTANCE or TIME.
    :param cruiseLength: The total length of the cruise segment [NM or s].
    :param stepSize: The discretization step for integrated calculations.

    :type cruiseType: myTypes.CruiseType
    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo | None
    :type controlTarget: myTypes.ControlTarget | None
    :type integrationType: myTypes.IntegrationType | None
    :type cruiseLength: float | None
    :type stepSize: float | None

    :returns: A Flight Trajectory object containing the results of the level flight calculation.
    :rtype: FT
    """

    trajectory = FT()

    match cruiseType:
        case CruiseType.CONSTANTSPEED:
            trajectory = hpcLevelConstantSpeed(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                speed=speed,
                mass=mass,
                meteo=meteo,
                integrationType=integrationType,
                cruiseLength=cruiseLength,
                stepSize=stepSize,
            )
        case CruiseType.ACCELERATION:
            trajectory = hpcLevelAccelerationDeceleration(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                speed=speed,
                mass=mass,
                meteo=meteo,
                controlTarget=controlTarget,
            )
    return trajectory


def apcLevelCalculation(
    cruiseType: myTypes.CruiseType,
    AC,
    calculationType: myTypes.CalculationType,
    pressureAltitude: myTypes.PressureAltitude,
    speed: myTypes.Speed,
    mass: float,
    meteo: myTypes.Meteo | None = None,
    controlTarget: myTypes.ControlTarget | None = None,
    integrationType: myTypes.IntegrationType | None = None,
    cruiseLength: float | None = None,
    stepSize: float | None = None,
    stepClimb: bool | None = None,
    speedBrakes: myTypes.SpeedBrakes | None = None,
    costIndex: float | None = None,
) -> FT:
    """Dispatches and calculates the fixed-wing aircraft level flight trajectory.

    :param cruiseType: The strategy for level flight (CONSTANTSPEED or ACCELERATION).
    :param AC: Aircraft model instance containing performance and optimization data.
    :param calculationType: The calculation method, either POINT or INTEGRATED.
    :param pressureAltitude: Altitude data containing initial, final, and step values.
    :param speed: Speed configuration; supports fixed values or optimized types (MRC, LRC, ECON).
    :param mass: Current aircraft mass in kilograms [kg].
    :param meteo: Meteorological data (wind speed and temperature deviation).
    :param controlTarget: Target parameters governing speed changes (for acceleration).
    :param integrationType: Defines if cruiseLength is DISTANCE or TIME.
    :param cruiseLength: The total length of the cruise segment [NM or s].
    :param stepSize: The discretization step for integrated calculations.
    :param stepClimb: Boolean flag to enable or disable automatic step-climb logic.
    :param speedBrakes: Configuration for speed brake state and deployment value.
    :param costIndex: The ratio of time-related costs to fuel costs for ECON speed.

    :type cruiseType: myTypes.CruiseType
    :type AC: Aircraft
    :type calculationType: myTypes.CalculationType
    :type pressureAltitude: myTypes.PressureAltitude
    :type speed: myTypes.Speed
    :type mass: float
    :type meteo: myTypes.Meteo | None
    :type controlTarget: myTypes.ControlTarget | None
    :type integrationType: myTypes.IntegrationType | None
    :type cruiseLength: float | None
    :type stepSize: float | None
    :type stepClimb: bool | None
    :type speedBrakes: myTypes.SpeedBrakes | None
    :type costIndex: float | None

    :returns: A Flight Trajectory object containing the results of the level flight calculation.
    :rtype: FT
    """

    trajectory = FT()

    match cruiseType:
        case CruiseType.CONSTANTSPEED:
            trajectory = apcLevelConstantSpeed(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                speed=speed,
                mass=mass,
                meteo=meteo,
                integrationType=integrationType,
                cruiseLength=cruiseLength,
                stepSize=stepSize,
                costIndex=costIndex,
                stepClimb=stepClimb,
            )
        case CruiseType.ACCELERATION:
            trajectory = apcLevelAccelerationDeceleration(
                AC=AC,
                calculationType=calculationType,
                pressureAltitude=pressureAltitude,
                speed=speed,
                mass=mass,
                meteo=meteo,
                controlTarget=controlTarget,
                speedBrakes=speedBrakes,
            )
    return trajectory


def apcFlightEnvelope(
    AC,
    mass: float,
    altitudeStep: float,
    meteo: myTypes.Meteo | None = None,
) -> dict[str, list]:
    """Generates the operational and certified flight envelopes for a fixed-wing aircraft.

    This function calculates the minimum and maximum speed limits across the
    aircraft's entire altitude range (up to HMO). It accounts for stall speeds,
    buffet limits, thrust-limited maximum speeds, and structural limits (VMO/MMO),
    while handling the transition at the crossover altitude.

    :param AC: Aircraft model instance containing flight envelope and limit data.
    :param mass: Current aircraft mass for weight-dependent speed calculations [kg].
    :param altitudeStep: The vertical increment for envelope sampling [ft].
    :param meteo: Meteorological data, specifically temperature deviation (deltaTemp).

    :type AC: Aircraft
    :type mass: float
    :type altitudeStep: float
    :type meteo: myTypes.Meteo | None

    :returns: A dictionary containing lists of speed limits for both envelopes.
    :rtype: dict
    """

    [theta, delta, sigma] = atm.atmosphereProperties(
        h=conv.ft2m(AC.hmo), deltaTemp=0
    )

    VMO = AC.VMO
    if AC.MMO is not None:
        MMO = AC.MMO
    else:
        MMO = atm.cas2Mach(cas=VMO, theta=theta, delta=delta, sigma=sigma)

    operationalFlightEnvelope = []
    certifiedFlightEnvelope = []

    crossoverAltitude = conv.m2ft(
        atm.crossOver(cas=conv.kt2ms(VMO), Mach=MMO)
    )  # [ft]
    tropopauseAltitude = conv.m2ft(const.h_11)  # [ft]

    # Create the basic altitude range including the top altitude.
    altitudeRange = np.arange(0, AC.hmo, altitudeStep)
    altitudeRange = np.append(altitudeRange, AC.hmo)

    if AC.BADAFamilyName == "BADA3":
        maxAltitude = conv.m2ft(
            AC.flightEnvelope.maxAltitude(mass=mass, deltaTemp=meteo.deltaTemp)
        )  # [ft]

        # Collect special altitudes that need to be included if they are less than AC.hmo.
        special_altitudes = [
            alt
            for alt in [crossoverAltitude, tropopauseAltitude, maxAltitude]
            if alt < AC.hmo
        ]

        # Merge the arrays, remove duplicates, and sort.
        altitudeRange = np.sort(
            np.unique(np.concatenate((altitudeRange, special_altitudes)))
        )

        for Hp in altitudeRange:
            alt_m = conv.ft2m(Hp)
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=alt_m, deltaTemp=meteo.deltaTemp
            )

            Vmin_operational = conv.ms2kt(
                AC.flightEnvelope.VMin(
                    h=alt_m, mass=mass, config="CR", deltaTemp=meteo.deltaTemp
                )
            )

            Vmax = AC.flightEnvelope.VMax(h=alt_m, deltaTemp=meteo.deltaTemp)
            buffetLimit = AC.flightEnvelope.lowSpeedBuffetLimit(
                h=alt_m, mass=mass, deltaTemp=meteo.deltaTemp, nz=1.0
            )

            VminCertified = conv.ms2kt(
                AC.flightEnvelope.VMin(
                    h=alt_m,
                    mass=mass,
                    config="CR",
                    deltaTemp=meteo.deltaTemp,
                    nz=1.0,
                    envelopeType="CERTIFIED",
                )
            )

            Vmax_thrustLimited = conv.ms2kt(
                AC.flightEnvelope.Vmax_thrustLimited(
                    h=alt_m,
                    mass=mass,
                    deltaTemp=meteo.deltaTemp,
                    rating="MCRZ",
                    config="CR",
                )
            )

            if Hp < crossoverAltitude:
                VmaxCertified = AC.VMO
                speedType = "CAS"
            else:
                VmaxCertified = AC.MMO
                speedType = "M"

            # aircraft speed
            if VminCertified is None:
                [VminCertified_M, VminCertified_CAS, VminCertified_TAS] = [
                    None,
                    None,
                    None,
                ]
            else:
                [VminCertified_M, VminCertified_CAS, VminCertified_TAS] = (
                    atm.convertSpeed(
                        v=VminCertified,
                        speedType="CAS",
                        theta=theta,
                        delta=delta,
                        sigma=sigma,
                    )
                )
                VminCertified_CAS = conv.ms2kt(VminCertified_CAS)
                VminCertified_TAS = conv.ms2kt(VminCertified_TAS)

            [VmaxCertified_M, VmaxCertified_CAS, VmaxCertified_TAS] = (
                atm.convertSpeed(
                    v=VmaxCertified,
                    speedType=speedType,
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                )
            )
            VmaxCertified_CAS = conv.ms2kt(VmaxCertified_CAS)
            VmaxCertified_TAS = conv.ms2kt(VmaxCertified_TAS)

            # limit the calcuation to where the max certified speed is lower than min certified speed
            if (
                VminCertified_CAS is not None
                and VmaxCertified_CAS is not None
                and VminCertified_CAS > VmaxCertified_CAS
            ):
                break

            if (
                Vmax_thrustLimited is None
                or Hp > maxAltitude
                or (Vmin_operational > Vmax_thrustLimited)
            ):
                [Vmax_M, Vmax_CAS, Vmax_TAS] = [None, None, None]
                [Vmin_M, Vmin_CAS, Vmin_TAS] = [None, None, None]

                operationalFlightEnvelope.append(
                    {
                        "altitude": Hp,
                        "Vmin_CAS": Vmin_CAS,
                        "Vmin_TAS": Vmin_TAS,
                        "Vmin_M": Vmin_M,
                        "Vmax_CAS": Vmax_CAS,
                        "Vmax_TAS": Vmax_TAS,
                        "Vmax_M": Vmax_M,
                    }
                )
                certifiedFlightEnvelope.append(
                    {
                        "altitude": Hp,
                        "Vmin_CAS": VminCertified_CAS,
                        "Vmin_TAS": VminCertified_TAS,
                        "Vmin_M": VminCertified_M,
                        "Vmax_CAS": VmaxCertified_CAS,
                        "Vmax_TAS": VmaxCertified_TAS,
                        "Vmax_M": VmaxCertified_M,
                    }
                )

            else:
                [Vmin_M, Vmin_CAS, Vmin_TAS] = atm.convertSpeed(
                    v=Vmin_operational,
                    speedType="CAS",
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                )
                [Vmax_M, Vmax_CAS, Vmax_TAS] = atm.convertSpeed(
                    v=Vmax_thrustLimited,
                    speedType="CAS",
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                )

                operationalFlightEnvelope.append(
                    {
                        "altitude": Hp,
                        "Vmin_CAS": conv.ms2kt(Vmin_CAS),
                        "Vmin_TAS": conv.ms2kt(Vmin_TAS),
                        "Vmin_M": Vmin_M,
                        "Vmax_CAS": conv.ms2kt(Vmax_CAS),
                        "Vmax_TAS": conv.ms2kt(Vmax_TAS),
                        "Vmax_M": Vmax_M,
                    }
                )
                certifiedFlightEnvelope.append(
                    {
                        "altitude": Hp,
                        "Vmin_CAS": VminCertified_CAS,
                        "Vmin_TAS": VminCertified_TAS,
                        "Vmin_M": VminCertified_M,
                        "Vmax_CAS": VmaxCertified_CAS,
                        "Vmax_TAS": VmaxCertified_TAS,
                        "Vmax_M": VmaxCertified_M,
                    }
                )

    elif AC.BADAFamilyName == "BADA4":
        # Collect special altitudes that need to be included if they are less than AC.hmo
        special_altitudes = [
            alt
            for alt in [crossoverAltitude, tropopauseAltitude]
            if alt < AC.hmo
        ]

        # Merge the arrays, remove duplicates, and sort.
        altitudeRange = np.sort(
            np.unique(np.concatenate((altitudeRange, special_altitudes)))
        )

        for Hp in altitudeRange:
            alt_m = conv.ft2m(Hp)
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=alt_m, deltaTemp=meteo.deltaTemp
            )

            Vmin_operational = conv.ms2kt(
                AC.flightEnvelope.VMin(
                    config="CR", theta=theta, delta=delta, mass=mass
                )
            )
            Vmax = AC.flightEnvelope.VMax(
                h=alt_m,
                HLid=0,
                LG="LGUP",
                delta=delta,
                theta=theta,
                mass=mass,
                nz=1.0,
            )
            VminCertified = conv.ms2kt(
                AC.flightEnvelope.VStall(
                    theta=theta,
                    delta=delta,
                    mass=mass,
                    HLid=0,
                    LG="LGUP",
                    nz=1.0,
                )
            )
            Vmax_thrustLimited = conv.ms2kt(
                AC.flightEnvelope.Vmax_thrustLimited(
                    h=alt_m,
                    mass=mass,
                    deltaTemp=meteo.deltaTemp,
                    rating="MCRZ",
                    config="CR",
                )
            )

            if Hp < crossoverAltitude:
                VmaxCertified = AC.VMO
                speedType = "CAS"
            else:
                VmaxCertified = AC.MMO
                speedType = "M"

            if VminCertified is None or VmaxCertified is None:
                break

            # aircraft speed
            [VminCertified_M, VminCertified_CAS, VminCertified_TAS] = (
                atm.convertSpeed(
                    v=VminCertified,
                    speedType="CAS",
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                )
            )
            [VmaxCertified_M, VmaxCertified_CAS, VmaxCertified_TAS] = (
                atm.convertSpeed(
                    v=VmaxCertified,
                    speedType=speedType,
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                )
            )

            # limit the calcuation to where the max certified speed is lower than min certified speed
            if VminCertified_CAS > VmaxCertified_CAS:
                break

            if (
                Vmax_thrustLimited is None
                or Vmin_operational is None
                or (Vmin_operational > Vmax_thrustLimited)
            ):
                [Vmax_M, Vmax_CAS, Vmax_TAS] = [None, None, None]
                [Vmin_M, Vmin_CAS, Vmin_TAS] = [None, None, None]

                operationalFlightEnvelope.append(
                    {
                        "altitude": Hp,
                        "Vmin_CAS": Vmin_CAS,
                        "Vmin_TAS": Vmin_TAS,
                        "Vmin_M": Vmin_M,
                        "Vmax_CAS": Vmax_CAS,
                        "Vmax_TAS": Vmax_TAS,
                        "Vmax_M": Vmax_M,
                    }
                )
                certifiedFlightEnvelope.append(
                    {
                        "altitude": Hp,
                        "Vmin_CAS": conv.ms2kt(VminCertified_CAS),
                        "Vmin_TAS": conv.ms2kt(VminCertified_TAS),
                        "Vmin_M": VminCertified_M,
                        "Vmax_CAS": conv.ms2kt(VmaxCertified_CAS),
                        "Vmax_TAS": conv.ms2kt(VmaxCertified_TAS),
                        "Vmax_M": VmaxCertified_M,
                    }
                )

            else:
                [Vmin_M, Vmin_CAS, Vmin_TAS] = atm.convertSpeed(
                    v=Vmin_operational,
                    speedType="CAS",
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                )
                [Vmax_M, Vmax_CAS, Vmax_TAS] = atm.convertSpeed(
                    v=Vmax_thrustLimited,
                    speedType="CAS",
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                )

                operationalFlightEnvelope.append(
                    {
                        "altitude": Hp,
                        "Vmin_CAS": conv.ms2kt(Vmin_CAS),
                        "Vmin_TAS": conv.ms2kt(Vmin_TAS),
                        "Vmin_M": Vmin_M,
                        "Vmax_CAS": conv.ms2kt(Vmax_CAS),
                        "Vmax_TAS": conv.ms2kt(Vmax_TAS),
                        "Vmax_M": Vmax_M,
                    }
                )
                certifiedFlightEnvelope.append(
                    {
                        "altitude": Hp,
                        "Vmin_CAS": conv.ms2kt(VminCertified_CAS),
                        "Vmin_TAS": conv.ms2kt(VminCertified_TAS),
                        "Vmin_M": VminCertified_M,
                        "Vmax_CAS": conv.ms2kt(VmaxCertified_CAS),
                        "Vmax_TAS": conv.ms2kt(VmaxCertified_TAS),
                        "Vmax_M": VmaxCertified_M,
                    }
                )

    return {
        "operationalFlightEnvelope": operationalFlightEnvelope,
        "certifiedFlightEnvelope": certifiedFlightEnvelope,
    }


def hpcFlightEnvelope(
    AC,
    mass: float,
    altitudeStep: float,
    meteo: myTypes.Meteo | None = None,
) -> dict[str, list]:
    """Generates the operational and certified flight envelopes for a helicopter.

    This function calculates the minimum and maximum speed boundaries for rotorcraft
    across a range of altitudes. It accounts for certified structural limits (VNE),
    atmospheric crossover altitudes, and power-limited speed envelopes where the
    maximum achievable speed is constrained by engine performance (MCNT rating).

    :param AC: Aircraft (Helicopter) model instance containing envelope data.
    :param mass: Current helicopter mass for performance-limited calculations [kg].
    :param altitudeStep: The vertical increment for envelope sampling [ft].
    :param meteo: Meteorological data, specifically temperature deviation (deltaTemp).

    :type AC: Aircraft
    :type mass: float
    :type altitudeStep: float
    :type meteo: myTypes.Meteo | None

    :returns: A dictionary containing lists of speed limits for both envelopes.
    :rtype: dict
    """

    [theta, delta, sigma] = atm.atmosphereProperties(
        h=conv.ft2m(AC.hmo), deltaTemp=0
    )
    VMO = AC.vne
    MMO = atm.cas2Mach(cas=VMO, theta=theta, delta=delta, sigma=sigma)

    operationalFlightEnvelope = []
    certifiedFlightEnvelope = []

    crossoverAltitude = conv.m2ft(
        atm.crossOver(cas=conv.kt2ms(VMO), Mach=MMO)
    )  # [ft]
    tropopauseAltitude = conv.m2ft(const.h_11)  # [ft]

    # Create the basic altitude range including the top altitude.
    altitudeRange = np.arange(0, AC.hmo, altitudeStep)
    altitudeRange = np.append(altitudeRange, AC.hmo)

    if AC.BADAFamilyName == "BADAH":
        maxAltitude = conv.m2ft(AC.flightEnvelope.maxAltitude())  # [ft]

        # Collect special altitudes that need to be included if they are less than AC.hmo.
        special_altitudes = [
            alt
            for alt in [crossoverAltitude, tropopauseAltitude, maxAltitude]
            if alt < AC.hmo
        ]

        # Merge the arrays, remove duplicates, and sort.
        altitudeRange = np.sort(
            np.unique(np.concatenate((altitudeRange, special_altitudes)))
        )

        for Hp in altitudeRange:
            alt_m = conv.ft2m(Hp)
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=alt_m, deltaTemp=meteo.deltaTemp
            )

            Vmax = AC.flightEnvelope.VMax()
            VminCertified = 0
            [Vmin_operational, Vmax_powerLimited] = (
                AC.flightEnvelope.speedEnvelope_powerLimited(
                    h=alt_m,
                    mass=mass,
                    deltaTemp=meteo.deltaTemp,
                    rating="MCNT",
                    rateOfTurn=0,
                )
            )

            if Hp < crossoverAltitude:
                VmaxCertified = VMO
                speedType = "CAS"
            else:
                VmaxCertified = MMO
                speedType = "M"

            # aircraft speed
            [VminCertified_M, VminCertified_CAS, VminCertified_TAS] = (
                atm.convertSpeed(
                    v=VminCertified,
                    speedType="CAS",
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                )
            )
            [VmaxCertified_M, VmaxCertified_CAS, VmaxCertified_TAS] = (
                atm.convertSpeed(
                    v=VmaxCertified,
                    speedType=speedType,
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                )
            )

            # limit the calcuation to where the max certified speed is lower than min certified speed
            if VminCertified_CAS > VmaxCertified_CAS:
                break

            if (
                Vmax_powerLimited is None
                or Hp > maxAltitude
                or Vmin_operational is None
                or (Vmin_operational > Vmax_powerLimited)
            ):
                [Vmax_M, Vmax_CAS, Vmax_TAS] = [None, None, None]
                [Vmin_M, Vmin_CAS, Vmin_TAS] = [None, None, None]

                operationalFlightEnvelope.append(
                    {
                        "altitude": Hp,
                        "Vmin_CAS": Vmin_CAS,
                        "Vmin_TAS": Vmin_TAS,
                        "Vmin_M": Vmin_M,
                        "Vmax_CAS": Vmax_CAS,
                        "Vmax_TAS": Vmax_TAS,
                        "Vmax_M": Vmax_M,
                    }
                )
                certifiedFlightEnvelope.append(
                    {
                        "altitude": Hp,
                        "Vmin_CAS": conv.ms2kt(VminCertified_CAS),
                        "Vmin_TAS": conv.ms2kt(VminCertified_TAS),
                        "Vmin_M": VminCertified_M,
                        "Vmax_CAS": conv.ms2kt(VmaxCertified_CAS),
                        "Vmax_TAS": conv.ms2kt(VmaxCertified_TAS),
                        "Vmax_M": VmaxCertified_M,
                    }
                )

            else:
                [Vmin_M, Vmin_CAS, Vmin_TAS] = atm.convertSpeed(
                    v=conv.ms2kt(Vmin_operational),
                    speedType="CAS",
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                )
                [Vmax_M, Vmax_CAS, Vmax_TAS] = atm.convertSpeed(
                    v=conv.ms2kt(Vmax_powerLimited),
                    speedType="CAS",
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                )

                operationalFlightEnvelope.append(
                    {
                        "altitude": Hp,
                        "Vmin_CAS": conv.ms2kt(Vmin_CAS),
                        "Vmin_TAS": conv.ms2kt(Vmin_TAS),
                        "Vmin_M": Vmin_M,
                        "Vmax_CAS": conv.ms2kt(Vmax_CAS),
                        "Vmax_TAS": conv.ms2kt(Vmax_TAS),
                        "Vmax_M": Vmax_M,
                    }
                )
                certifiedFlightEnvelope.append(
                    {
                        "altitude": Hp,
                        "Vmin_CAS": conv.ms2kt(VminCertified_CAS),
                        "Vmin_TAS": conv.ms2kt(VminCertified_TAS),
                        "Vmin_M": VminCertified_M,
                        "Vmax_CAS": conv.ms2kt(VmaxCertified_CAS),
                        "Vmax_TAS": conv.ms2kt(VmaxCertified_TAS),
                        "Vmax_M": VmaxCertified_M,
                    }
                )

    elif AC.BADAFamilyName == "BADAE":
        maxAltitude = conv.m2ft(AC.flightEnvelope.maxAltitude())  # [ft]

        # Collect special altitudes that need to be included if they are less than AC.hmo.
        special_altitudes = [
            alt
            for alt in [crossoverAltitude, tropopauseAltitude, maxAltitude]
            if alt < AC.hmo
        ]

        # Merge the arrays, remove duplicates, and sort.
        altitudeRange = np.sort(
            np.unique(np.concatenate((altitudeRange, special_altitudes)))
        )

        for Hp in altitudeRange:
            alt_m = conv.ft2m(Hp)
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=alt_m, deltaTemp=meteo.deltaTemp
            )

            Vmax = AC.flightEnvelope.VMax()
            VminCertified = 0

            # power limited min/max speed is not yet defined for the BADAE, so for now a simplified calculation will be provided
            # [Vmin, Vmax_powerLimited] = AC.flightEnvelope.speedEnvelope_powerLimited(h=alt_m, mass=acMass, deltaTemp=deltaTemp, rating='MCNT', rateOfTurn=0)
            Vmin = 0
            Vmax_powerLimited = Vmax

            if Hp < crossoverAltitude:
                VmaxCertified = VMO
                speedType = "CAS"
            else:
                VmaxCertified = MMO
                speedType = "M"

            # aircraft speed
            [VminCertified_M, VminCertified_CAS, VminCertified_TAS] = (
                atm.convertSpeed(
                    v=VminCertified,
                    speedType="CAS",
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                )
            )
            [VmaxCertified_M, VmaxCertified_CAS, VmaxCertified_TAS] = (
                atm.convertSpeed(
                    v=VmaxCertified,
                    speedType=speedType,
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                )
            )

            if Vmax_powerLimited is None or Hp > maxAltitude:
                [Vmax_M, Vmax_CAS, Vmax_TAS] = [None, None, None]
                [Vmin_M, Vmin_CAS, Vmin_TAS] = [None, None, None]

                operationalFlightEnvelope.append(
                    {
                        "altitude": Hp,
                        "Vmin_CAS": Vmin_CAS,
                        "Vmin_TAS": Vmin_TAS,
                        "Vmin_M": Vmin_M,
                        "Vmax_CAS": Vmax_CAS,
                        "Vmax_TAS": Vmax_TAS,
                        "Vmax_M": Vmax_M,
                    }
                )
                certifiedFlightEnvelope.append(
                    {
                        "altitude": Hp,
                        "Vmin_CAS": conv.ms2kt(VminCertified_CAS),
                        "Vmin_TAS": conv.ms2kt(VminCertified_TAS),
                        "Vmin_M": VminCertified_M,
                        "Vmax_CAS": conv.ms2kt(VmaxCertified_CAS),
                        "Vmax_TAS": conv.ms2kt(VmaxCertified_TAS),
                        "Vmax_M": VmaxCertified_M,
                    }
                )

            else:
                [Vmin_M, Vmin_CAS, Vmin_TAS] = atm.convertSpeed(
                    v=conv.ms2kt(Vmin),
                    speedType="CAS",
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                )
                [Vmax_M, Vmax_CAS, Vmax_TAS] = atm.convertSpeed(
                    v=conv.ms2kt(Vmax_powerLimited),
                    speedType="CAS",
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                )

                operationalFlightEnvelope.append(
                    {
                        "altitude": Hp,
                        "Vmin_CAS": conv.ms2kt(Vmin_CAS),
                        "Vmin_TAS": conv.ms2kt(Vmin_TAS),
                        "Vmin_M": Vmin_M,
                        "Vmax_CAS": conv.ms2kt(Vmax_CAS),
                        "Vmax_TAS": conv.ms2kt(Vmax_TAS),
                        "Vmax_M": Vmax_M,
                    }
                )
                certifiedFlightEnvelope.append(
                    {
                        "altitude": Hp,
                        "Vmin_CAS": conv.ms2kt(VminCertified_CAS),
                        "Vmin_TAS": conv.ms2kt(VminCertified_TAS),
                        "Vmin_M": VminCertified_M,
                        "Vmax_CAS": conv.ms2kt(VmaxCertified_CAS),
                        "Vmax_TAS": conv.ms2kt(VmaxCertified_TAS),
                        "Vmax_M": VmaxCertified_M,
                    }
                )

    return {
        "operationalFlightEnvelope": operationalFlightEnvelope,
        "certifiedFlightEnvelope": certifiedFlightEnvelope,
    }


def flightEnvelope(
    AC,
    mass: float,
    altitudeStep: float,
    meteo: myTypes.Meteo | None = None,
) -> dict[str, list]:
    """Dispatches the flight envelope calculation to the appropriate aircraft or helicopter model.

    :param AC: Aircraft model instance containing the BADA family name and performance data.
    :param mass: Current aircraft mass for weight-dependent speed limits [kg].
    :param altitudeStep: The vertical increment for sampling the envelope [ft].
    :param meteo: Meteorological data, specifically temperature deviation (deltaTemp).

    :type AC: Aircraft
    :type mass: float
    :type altitudeStep: float
    :type meteo: myTypes.Meteo | None

    :returns: A dictionary containing 'operationalFlightEnvelope' and 'certifiedFlightEnvelope' data.
    :rtype: dict
    """

    match AC.BADAFamilyName:
        case "BADA3" | "BADA4":
            flightEnvelopeData = apcFlightEnvelope(
                AC=AC, mass=mass, meteo=meteo, altitudeStep=altitudeStep
            )
        case "BADAH" | "BADAE":
            flightEnvelopeData = hpcFlightEnvelope(
                AC=AC, mass=mass, meteo=meteo, altitudeStep=altitudeStep
            )

    return flightEnvelopeData
