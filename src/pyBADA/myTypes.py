from dataclasses import dataclass, field
from enum import StrEnum
from typing import Literal


class CalculationType(StrEnum):
    POINT = "POINT"
    INTEGRATED = "INTEGRATED"


class FlightPhase(StrEnum):
    CLIMB = "CLIMB"
    DESCENT = "DESCENT"
    CRUISE = "CRUISE"


class AccelerationOption(StrEnum):
    ROCD_ESF = "ROCD_ESF"
    ROCD_ACC = "ROCD_ACC"
    SLOPE_ACC = "SLOPE_ACC"
    SLOPE_ESF = "SLOPE_ESF"
    ROCD_RATING = "ROCD_RATING"
    SLOPE_RATING = "SLOPE_RATING"
    ESF_RATING = "ESF_RATING"
    ACC_RATING = "ACC_RATING"


class SpeedType(StrEnum):
    TAS = "TAS"
    CAS = "CAS"
    M = "M"


class AccelerationLevelKind(StrEnum):
    AT = "AT"
    BEFORE = "BEFORE"
    AFTER = "AFTER"


@dataclass
class Speed:
    speedType: SpeedType | None = None
    accelerationLevelKind: AccelerationLevelKind | None = None
    initSpeed: float | None = None
    finalSpeed: float | None = None
    stepSpeed: float | None = None


@dataclass
class SpeedBrakes:
    deployed: float | None = None
    value: float | None = None


@dataclass
class CASMACHSpeedSchedule:
    CASbelowFL100: float | None = None
    CASaboveFL100: float | None = None
    Mach: float | None = None


@dataclass
class Meteo:
    wS: float = 0.0  # Wind speed [kt], default to 0
    deltaTemp: float = 0.0  # Temperature deviation [K], default to 0


@dataclass
class PressureAltitude:
    initPressureAltitude: float | None = None
    finalPressureAltitude: float | None = None
    stepPressureAltitude: float | None = 1000.0


@dataclass
class ControlTarget:
    ROCDtarget: float | None = None
    slopetarget: float | None = None
    acctarget: float | None = None
    ESFtarget: float | None = None


class TakeOffProcedureType(StrEnum):
    BADA = "BADA"
    NADP1 = "NADP1"
    NADP2 = "NADP2"


@dataclass
class TakeOffProcedureBADA:
    type: Literal[TakeOffProcedureType.BADA] = TakeOffProcedureType.BADA


@dataclass
class TakeOffProcedureNADP1:
    NADP1Threshold: float | None = None
    type: Literal[TakeOffProcedureType.NADP1] = TakeOffProcedureType.NADP1


@dataclass
class TakeOffProcedureNADP2:
    NADP2Threshold1: float | None = None
    NADP2Threshold2: float | None = None
    type: Literal[TakeOffProcedureType.NADP2] = TakeOffProcedureType.NADP2


TakeOffProcedure = (
    TakeOffProcedureBADA | TakeOffProcedureNADP1 | TakeOffProcedureNADP2
)


class DepartureProfileType(StrEnum):
    constCASbelow100 = "constCASbelow100"
    calculatedCAS = "calculatedCAS"


@dataclass
class DepartureProfile:
    departureProfileType: DepartureProfileType | None = (
        DepartureProfileType.calculatedCAS
    )


class ArrivalProfileType(StrEnum):
    constCASbelow100 = "constCASbelow100"
    calculatedCAS = "calculatedCAS"
    expedite = "expedite"


@dataclass
class ArrivalProfile:
    arrivalProfileType: ArrivalProfileType | None = (
        ArrivalProfileType.calculatedCAS
    )


class ClimbType(StrEnum):
    CASMACH = "CASMACH"
    RATE = "RATE"
    ACCELERATION = "ACCELERATION"


class DescentType(StrEnum):
    CASMACH = "CASMACH"
    RATE = "RATE"
    SLOPE = "SLOPE"
    ACCELERATION = "ACCELERATION"
    EMERGENCY = "EMERGENCY"


class CruiseType(StrEnum):
    CONSTANTSPEED = "CONSTANTSPEED"
    ACCELERATION = "ACCELERATION"


class CruiseSpeedType(StrEnum):
    TAS = "TAS"
    CAS = "CAS"
    M = "M"
    MEC = "MEC"
    LRC = "LRC"
    MRC = "MRC"
    ECON = "ECON"


class HDescentType(StrEnum):
    ARPM = "ARPM"
    RATE = "RATE"
    SLOPE = "SLOPE"
    ACCELERATION = "ACCELERATION"
    EMERGENCY = "EMERGENCY"


class HClimbType(StrEnum):
    RATE = "RATE"
    ARPM = "ARPM"
    ACCELERATION = "ACCELERATION"


class HRating(StrEnum):
    MTKF = "MTKF"
    MCNT = "MCNT"
    ARPM = "ARPM"


class IntegrationType(StrEnum):
    ALTITUDE = "ALTITUDE"
    DISTANCE = "DISTANCE"
    TIME = "TIME"


@dataclass
class HClimbRatingConfiguration:
    rating: HRating | None = None


@dataclass
class ClimbCASMACHProfileConfiguration:
    casMachSpeedSchedule: CASMACHSpeedSchedule | None = None
    takeOffProcedure: TakeOffProcedure | None = None
    departureProfile: DepartureProfile | None = None
    reducedPower: bool = False


@dataclass
class DescentCASMACHProfileConfiguration:
    casMachSpeedSchedule: CASMACHSpeedSchedule = field(
        default_factory=CASMACHSpeedSchedule
    )
    arrivalProfile: ArrivalProfile = field(default_factory=ArrivalProfile)
