"""
Aircraft Trajectory Calculation on Earth
========================================

Example of BADA3 and BADA4 trajectory including geodesic calculations using TCL
"""

from dataclasses import dataclass
import matplotlib.pyplot as plt

from pyBADA import atmosphere as atm
from pyBADA import conversions as conv
from pyBADA import TCL as TCL
from pyBADA.flightTrajectory import FlightTrajectory as FT
from pyBADA.magnetic import Grid

from pyBADA.bada4 import Bada4Aircraft
from pyBADA.bada4 import Parser as Bada4Parser


@dataclass
class target:
    ROCDtarget: float = None
    slopetarget: float = None
    acctarget: float = None
    ESFtarget: float = None


badaVersion = "DUMMY"

# allData = Bada3Parser.parseAll(badaVersion=badaVersion)
allData = Bada4Parser.parseAll(badaVersion=badaVersion)
print(allData)

# AC = Bada3Aircraft(badaVersion=badaVersion, acName='J2H', allData=allData)
AC = Bada4Aircraft(
    badaVersion=badaVersion, acName="Dummy-TWIN", allData=allData
)

# get magnetic declination data
magneticDeclinationGrid = Grid()

# create a Flight Trajectory object to store the output from TCL segment calculations
ft = FT()

# default parameters
speedType = "CAS"  # {M, CAS, TAS}
wS = 0  # [kt] wind speed
ba = 0  # [deg] bank angle
DeltaTemp = 0  # [K] delta temperature from ISA

# Initial conditions
m_init = AC.OEW + 0.7 * (AC.MTOW - AC.OEW)  # [kg] initial mass
CAS_init = 170  # [kt] Initial CAS
Hp_RWY = 318.0  # [ft] CDG RWY26R elevation
Lat_init = 48.9982052030771  # CDG RWY26R coordinates
Lon_init = 2.5995367285775965  # CDG RWY26R coordinates

# take-off conditions
[theta, delta, sigma] = atm.atmosphereProperties(
    h=conv.ft2m(Hp_RWY), DeltaTemp=DeltaTemp
)  # atmosphere properties at RWY altitude
[cas_cl1, speedUpdated] = AC.ARPM.climbSpeed(
    h=conv.ft2m(Hp_RWY),
    mass=m_init,
    theta=theta,
    delta=delta,
    DeltaTemp=DeltaTemp,
)  # [m/s] take-off CAS


Hp_CR = 33000  # [ft] CRUISing level

# BADA speed schedule
[Vcl1, Vcl2, Mcl] = AC.flightEnvelope.getSpeedSchedule(
    phase="Climb"
)  # BADA Climb speed schedule
[Vcr1, Vcr2, Mcr] = AC.flightEnvelope.getSpeedSchedule(
    phase="Cruise"
)  # BADA Cruise speed schedule
[Vdes1, Vdes2, Mdes] = AC.flightEnvelope.getSpeedSchedule(
    phase="Descent"
)  # BADA Descent speed schedule

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# CLIMB to threshold altitude 1500ft at take-off speed
# ------------------------------------------------
flightTrajectory = TCL.constantSpeedRating(
    AC=AC,
    speedType="CAS",
    v=conv.ms2kt(cas_cl1),
    Hp_init=Hp_RWY,
    Hp_final=1499,
    m_init=m_init,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=Lat_init,
    Lon=Lon_init,
    initialHeading={"magnetic": None, "true": 50, "constantHeading": False},
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)

# accelerate according to BADA ARPM for below 3000ft
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

[theta, delta, sigma] = atm.atmosphereProperties(
    h=conv.ft2m(2999), DeltaTemp=DeltaTemp
)
[cas_cl2, speedUpdated] = AC.ARPM.climbSpeed(
    h=conv.ft2m(2999),
    mass=m_final,
    theta=theta,
    delta=delta,
    DeltaTemp=DeltaTemp,
)

flightTrajectory = TCL.accDec(
    AC=AC,
    speedType="CAS",
    v_init=CAS_final,
    v_final=conv.ms2kt(cas_cl2),
    Hp_init=Hp,
    control=None,
    phase="Climb",
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)

# CLIMB to threshold altitude 3000ft
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

flightTrajectory = TCL.constantSpeedRating(
    AC=AC,
    speedType="CAS",
    v=CAS_final,
    Hp_init=Hp,
    Hp_final=2999,
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)


# accelerate according to BADA ARPM for below 4000ft
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

[theta, delta, sigma] = atm.atmosphereProperties(
    h=conv.ft2m(3999), DeltaTemp=DeltaTemp
)
[cas_cl3, speedUpdated] = AC.ARPM.climbSpeed(
    h=conv.ft2m(3999),
    mass=m_final,
    theta=theta,
    delta=delta,
    DeltaTemp=DeltaTemp,
)

flightTrajectory = TCL.accDec(
    AC=AC,
    speedType="CAS",
    v_init=CAS_final,
    v_final=conv.ms2kt(cas_cl3),
    Hp_init=Hp,
    control=None,
    phase="Climb",
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)

# CLIMB to threshold altitude 4000ft
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

flightTrajectory = TCL.constantSpeedRating(
    AC=AC,
    speedType="CAS",
    v=CAS_final,
    Hp_init=Hp,
    Hp_final=3999,
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)


# accelerate according to BADA ARPM for below 5000ft
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

[theta, delta, sigma] = atm.atmosphereProperties(
    h=conv.ft2m(4999), DeltaTemp=DeltaTemp
)
[cas_cl4, speedUpdated] = AC.ARPM.climbSpeed(
    h=conv.ft2m(4999),
    mass=m_final,
    theta=theta,
    delta=delta,
    DeltaTemp=DeltaTemp,
)

flightTrajectory = TCL.accDec(
    AC=AC,
    speedType="CAS",
    v_init=CAS_final,
    v_final=conv.ms2kt(cas_cl4),
    Hp_init=Hp,
    control=None,
    phase="Climb",
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)

# CLIMB to threshold altitude 5000ft
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

flightTrajectory = TCL.constantSpeedRating(
    AC=AC,
    speedType="CAS",
    v=CAS_final,
    Hp_init=Hp,
    Hp_final=4999,
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)


# accelerate according to BADA ARPM for below 6000ft
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

[theta, delta, sigma] = atm.atmosphereProperties(
    h=conv.ft2m(5999), DeltaTemp=DeltaTemp
)
[cas_cl5, speedUpdated] = AC.ARPM.climbSpeed(
    h=conv.ft2m(5999),
    mass=m_final,
    theta=theta,
    delta=delta,
    DeltaTemp=DeltaTemp,
)

flightTrajectory = TCL.accDec(
    AC=AC,
    speedType="CAS",
    v_init=CAS_final,
    v_final=conv.ms2kt(cas_cl5),
    Hp_init=Hp,
    control=None,
    phase="Climb",
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)

# CLIMB to threshold altitude 6000ft
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

flightTrajectory = TCL.constantSpeedRating(
    AC=AC,
    speedType="CAS",
    v=CAS_final,
    Hp_init=Hp,
    Hp_final=5999,
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)


# accelerate according to BADA ARPM for below 10000ft
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

[theta, delta, sigma] = atm.atmosphereProperties(
    h=conv.ft2m(9999), DeltaTemp=DeltaTemp
)
[cas_cl6, speedUpdated] = AC.ARPM.climbSpeed(
    h=conv.ft2m(9999),
    mass=m_final,
    theta=theta,
    delta=delta,
    DeltaTemp=DeltaTemp,
)

flightTrajectory = TCL.accDec(
    AC=AC,
    speedType="CAS",
    v_init=CAS_final,
    v_final=conv.ms2kt(cas_cl6),
    Hp_init=Hp,
    control=None,
    phase="Climb",
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)

# CLIMB to threshold altitude 10000ft
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

flightTrajectory = TCL.constantSpeedRating(
    AC=AC,
    speedType="CAS",
    v=CAS_final,
    Hp_init=Hp,
    Hp_final=9999,
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)


# accelerate according to BADA ARPM for above 10000ft and below crossover altitude
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

flightTrajectory = TCL.accDec(
    AC=AC,
    speedType="CAS",
    v_init=CAS_final,
    v_final=conv.ms2kt(Vcl2),
    Hp_init=Hp,
    control=None,
    phase="Climb",
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)


# CLIMB to crossover altitude
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

# calculate the crosover altitude for climb phase
crossoverAltitude = conv.m2ft(atm.crossOver(Vcl2, Mcl))

flightTrajectory = TCL.constantSpeedRating(
    AC=AC,
    speedType="CAS",
    v=CAS_final,
    Hp_init=Hp,
    Hp_final=crossoverAltitude,
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)

# climb at M from crossover altitude
# ------------------------------------------------
# current values
Hp, m_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "LAT", "LON", "HDGTrue"]
)

flightTrajectory = TCL.constantSpeedRating(
    AC=AC,
    speedType="M",
    v=Mcl,
    Hp_init=Hp,
    Hp_final=Hp_CR,
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)

# if not at CR speed -> adapt the speed first (acc/dec)
# ------------------------------------------------
# current values
Hp, m_final, M_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "M", "LAT", "LON", "HDGTrue"]
)

if M_final < Mcr:
    control = target(acctarget=0.5)
    flightTrajectory = TCL.accDec(
        AC=AC,
        speedType="M",
        v_init=M_final,
        v_final=Mcr,
        Hp_init=Hp,
        control=control,
        phase="Cruise",
        m_init=m_final,
        wS=wS,
        bankAngle=ba,
        DeltaTemp=DeltaTemp,
        Lat=LAT_final,
        Lon=LON_final,
        initialHeading={
            "magnetic": None,
            "true": HDGTrue,
            "constantHeading": False,
        },
        magneticDeclinationGrid=magneticDeclinationGrid,
    )
    ft.append(AC, flightTrajectory)

# CRUISE for 200 NM
# ------------------------------------------------
# current values
Hp, m_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "LAT", "LON", "HDGTrue"]
)

flightTrajectory = TCL.constantSpeedLevel(
    AC=AC,
    lengthType="distance",
    length=200,
    speedType="M",
    v=Mcr,
    Hp_init=Hp,
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)

# CRUISE Step for 300 NM
# ------------------------------------------------
# current values
# Hp, m_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
# AC, ["Hp", "mass", "LAT", "LON", "HDGTrue"]
# )

# flightTrajectory = TCL.constantSpeedLevel(
# AC=AC,
# lengthType="distance",
# length=200,
# step_length=50,
# maxRFL=36000,
# speedType="M",
# v=Mcr,
# Hp_init=Hp,
# m_init=m_final,
# stepClimb=True,
# wS=wS,
# bankAngle=ba,
# DeltaTemp=DeltaTemp,
# Lat=LAT_final,
# Lon=LON_final,
# initialHeading={"magnetic":None, "true":HDGTrue, "constantHeading":False},
# magneticDeclinationGrid=magneticDeclinationGrid
# )
# ft.append(AC, flightTrajectory)

# acc/dec to DESCENT speed during the descend
# ------------------------------------------------
# current values
Hp, m_final, M_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "M", "LAT", "LON", "HDGTrue"]
)

flightTrajectory = TCL.accDec(
    AC=AC,
    speedType="M",
    v_init=M_final,
    v_final=Mdes,
    Hp_init=Hp,
    phase="Descent",
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)

# descend to crossover altitude
# ------------------------------------------------
# current values
Hp, m_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "LAT", "LON", "HDGTrue"]
)

# calculate the crosover altitude for descend phase
crossoverAltitude = conv.m2ft(atm.crossOver(Vdes2, Mdes))

flightTrajectory = TCL.constantSpeedRating(
    AC=AC,
    speedType="M",
    v=Mdes,
    Hp_init=Hp,
    Hp_final=crossoverAltitude,
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)

# descend to FL100
# ------------------------------------------------
# current values
Hp, m_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "LAT", "LON", "HDGTrue"]
)

flightTrajectory = TCL.constantSpeedRating(
    AC=AC,
    speedType="CAS",
    v=conv.ms2kt(Vdes2),
    Hp_init=Hp,
    Hp_final=10000,
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)

# decelerate according to BADA ARPM for below FL100
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

# get BADA target speed from BADA ARPM procedure for the altitude bracket below
[theta, delta, sigma] = atm.atmosphereProperties(
    h=conv.ft2m(9999), DeltaTemp=DeltaTemp
)
[cas, speedUpdated] = AC.ARPM.descentSpeed(
    h=conv.ft2m(9999),
    mass=m_final,
    theta=theta,
    delta=delta,
    DeltaTemp=DeltaTemp,
)

flightTrajectory = TCL.accDec(
    AC=AC,
    speedType="CAS",
    v_init=CAS_final,
    v_final=conv.ms2kt(cas),
    Hp_init=Hp,
    phase="Descent",
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)

# descend to 6000ft
# ------------------------------------------------
# current values
Hp, m_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "LAT", "LON", "HDGTrue"]
)

flightTrajectory = TCL.constantSpeedRating(
    AC=AC,
    speedType="CAS",
    v=conv.ms2kt(cas),
    Hp_init=Hp,
    Hp_final=6000,
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)

# decelerate according to BADA ARPM for below 6000
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

# get BADA target speed from BADA ARPM procedure for the altitude bracket below
[theta, delta, sigma] = atm.atmosphereProperties(
    h=conv.ft2m(5999), DeltaTemp=DeltaTemp
)
[cas, speedUpdated] = AC.ARPM.descentSpeed(
    h=conv.ft2m(5999),
    mass=m_final,
    theta=theta,
    delta=delta,
    DeltaTemp=DeltaTemp,
)

flightTrajectory = TCL.accDec(
    AC=AC,
    speedType="CAS",
    v_init=CAS_final,
    v_final=conv.ms2kt(cas),
    Hp_init=Hp,
    phase="Descent",
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)

# descend to 5000ft
# ------------------------------------------------
# current values
Hp, m_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "LAT", "LON", "HDGTrue"]
)

flightTrajectory = TCL.constantSpeedRating(
    AC=AC,
    speedType="CAS",
    v=conv.ms2kt(cas),
    Hp_init=Hp,
    Hp_final=5000,
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)


# descend on ILS with 3deg glideslope to next altitude threshold
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

if AC.BADAFamily.BADA3:
    flightTrajectory = TCL.constantSpeedSlope(
        AC=AC,
        speedType="CAS",
        v=CAS_final,
        Hp_init=Hp,
        Hp_final=3700,
        slopetarget=-3.0,
        config="AP",
        m_init=m_final,
        wS=wS,
        bankAngle=ba,
        DeltaTemp=DeltaTemp,
        Lat=LAT_final,
        Lon=LON_final,
        initialHeading={
            "magnetic": None,
            "true": HDGTrue,
            "constantHeading": False,
        },
        magneticDeclinationGrid=magneticDeclinationGrid,
    )
elif AC.BADAFamily.BADA4:
    flightTrajectory = TCL.constantSpeedSlope(
        AC=AC,
        speedType="CAS",
        v=CAS_final,
        Hp_init=Hp,
        Hp_final=3000,
        slopetarget=-3.0,
        config=None,
        m_init=m_final,
        wS=wS,
        bankAngle=ba,
        DeltaTemp=DeltaTemp,
        Lat=LAT_final,
        Lon=LON_final,
        initialHeading={
            "magnetic": None,
            "true": HDGTrue,
            "constantHeading": False,
        },
        magneticDeclinationGrid=magneticDeclinationGrid,
    )

ft.append(AC, flightTrajectory)


# descend on ILS with 3deg glideslope while decelerating
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

# get BADA target speed from BADA ARPM procedure for the altitude bracket below
[theta, delta, sigma] = atm.atmosphereProperties(
    h=conv.ft2m(2999), DeltaTemp=DeltaTemp
)
[cas, speedUpdated] = AC.ARPM.descentSpeed(
    h=conv.ft2m(2999),
    mass=m_final,
    theta=theta,
    delta=delta,
    DeltaTemp=DeltaTemp,
)

control = target(slopetarget=-3.0)
flightTrajectory = TCL.accDec(
    AC=AC,
    speedType="CAS",
    v_init=CAS_final,
    v_final=conv.ms2kt(cas),
    Hp_init=Hp,
    control=control,
    phase="Descent",
    config="AP",
    speedBrakes={"deployed": True, "value": 0.03},
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)


# descend on ILS with 3deg glideslope to next altitude threshold
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

if Hp > 2000:
    flightTrajectory = TCL.constantSpeedSlope(
        AC=AC,
        speedType="CAS",
        v=CAS_final,
        Hp_init=Hp,
        Hp_final=2000,
        slopetarget=-3.0,
        config=None,
        m_init=m_final,
        wS=wS,
        bankAngle=ba,
        DeltaTemp=DeltaTemp,
        Lat=LAT_final,
        Lon=LON_final,
        initialHeading={
            "magnetic": None,
            "true": HDGTrue,
            "constantHeading": False,
        },
        magneticDeclinationGrid=magneticDeclinationGrid,
    )
    ft.append(AC, flightTrajectory)


# descend on ILS with 3deg glideslope while decelerating
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

# get BADA target speed from BADA ARPM procedure for the altitude bracket below
[theta, delta, sigma] = atm.atmosphereProperties(
    h=conv.ft2m(1999), DeltaTemp=DeltaTemp
)
[cas, speedUpdated] = AC.ARPM.descentSpeed(
    h=conv.ft2m(1999),
    mass=m_final,
    theta=theta,
    delta=delta,
    DeltaTemp=DeltaTemp,
)

control = target(slopetarget=-3.0)
flightTrajectory = TCL.accDec(
    AC=AC,
    speedType="CAS",
    v_init=CAS_final,
    v_final=conv.ms2kt(cas),
    Hp_init=Hp,
    control=control,
    phase="Descent",
    config="LD",
    speedBrakes={"deployed": True, "value": 0.03},
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)


# descend on ILS with 3deg glideslope to next altitude threshold
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

if Hp > 1500:
    flightTrajectory = TCL.constantSpeedSlope(
        AC=AC,
        speedType="CAS",
        v=CAS_final,
        Hp_init=Hp,
        Hp_final=1500,
        slopetarget=-3.0,
        config="LD",
        m_init=m_final,
        wS=wS,
        bankAngle=ba,
        DeltaTemp=DeltaTemp,
        Lat=LAT_final,
        Lon=LON_final,
        initialHeading={
            "magnetic": None,
            "true": HDGTrue,
            "constantHeading": False,
        },
        magneticDeclinationGrid=magneticDeclinationGrid,
    )
    ft.append(AC, flightTrajectory)


# descend on ILS with 3deg glideslope while decelerating
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

# get BADA target speed from BADA ARPM procedure for the altitude bracket below
[theta, delta, sigma] = atm.atmosphereProperties(
    h=conv.ft2m(1499), DeltaTemp=DeltaTemp
)
[cas, speedUpdated] = AC.ARPM.descentSpeed(
    h=conv.ft2m(1499),
    mass=m_final,
    theta=theta,
    delta=delta,
    DeltaTemp=DeltaTemp,
)

control = target(slopetarget=-3.0)
if AC.BADAFamily.BADA3:
    flightTrajectory = TCL.accDec(
        AC=AC,
        speedType="CAS",
        v_init=CAS_final,
        v_final=conv.ms2kt(cas),
        Hp_init=Hp,
        control=control,
        phase="Descent",
        config="LD",
        speedBrakes={"deployed": True, "value": 0.03},
        m_init=m_final,
        wS=wS,
        bankAngle=ba,
        DeltaTemp=DeltaTemp,
        Lat=LAT_final,
        Lon=LON_final,
        initialHeading={
            "magnetic": None,
            "true": HDGTrue,
            "constantHeading": False,
        },
        magneticDeclinationGrid=magneticDeclinationGrid,
    )
elif AC.BADAFamily.BADA4:
    flightTrajectory = TCL.accDec(
        AC=AC,
        speedType="CAS",
        v_init=CAS_final,
        v_final=conv.ms2kt(cas),
        Hp_init=Hp,
        control=control,
        phase="Descent",
        config="LD",
        m_init=m_final,
        wS=wS,
        bankAngle=ba,
        DeltaTemp=DeltaTemp,
        Lat=LAT_final,
        Lon=LON_final,
        initialHeading={
            "magnetic": None,
            "true": HDGTrue,
            "constantHeading": False,
        },
        magneticDeclinationGrid=magneticDeclinationGrid,
    )
ft.append(AC, flightTrajectory)


# descend on ILS with 3deg glideslope to next altitude threshold
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

if Hp > 1000:
    flightTrajectory = TCL.constantSpeedSlope(
        AC=AC,
        speedType="CAS",
        v=CAS_final,
        Hp_init=Hp,
        Hp_final=1000,
        slopetarget=-3.0,
        config=None,
        m_init=m_final,
        wS=wS,
        bankAngle=ba,
        DeltaTemp=DeltaTemp,
        Lat=LAT_final,
        Lon=LON_final,
        initialHeading={
            "magnetic": None,
            "true": HDGTrue,
            "constantHeading": False,
        },
        magneticDeclinationGrid=magneticDeclinationGrid,
    )
    ft.append(AC, flightTrajectory)

# descend on ILS with 3deg glideslope while decelerating
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

# get BADA target speed from BADA ARPM procedure for the altitude bracket below
[theta, delta, sigma] = atm.atmosphereProperties(
    h=conv.ft2m(999), DeltaTemp=DeltaTemp
)
[cas, speedUpdated] = AC.ARPM.descentSpeed(
    h=conv.ft2m(999),
    mass=m_final,
    theta=theta,
    delta=delta,
    DeltaTemp=DeltaTemp,
)

control = target(slopetarget=-3.0)
if AC.BADAFamily.BADA3:
    flightTrajectory = TCL.accDec(
        AC=AC,
        speedType="CAS",
        v_init=CAS_final,
        v_final=conv.ms2kt(cas),
        Hp_init=Hp,
        control=control,
        phase="Descent",
        config=None,
        speedBrakes={"deployed": True, "value": 0.03},
        m_init=m_final,
        wS=wS,
        bankAngle=ba,
        DeltaTemp=DeltaTemp,
        Lat=LAT_final,
        Lon=LON_final,
        initialHeading={
            "magnetic": None,
            "true": HDGTrue,
            "constantHeading": False,
        },
        magneticDeclinationGrid=magneticDeclinationGrid,
    )
elif AC.BADAFamily.BADA4:
    flightTrajectory = TCL.accDec(
        AC=AC,
        speedType="CAS",
        v_init=CAS_final,
        v_final=conv.ms2kt(cas),
        Hp_init=Hp,
        control=control,
        phase="Descent",
        config=None,
        m_init=m_final,
        wS=wS,
        bankAngle=ba,
        DeltaTemp=DeltaTemp,
        Lat=LAT_final,
        Lon=LON_final,
        initialHeading={
            "magnetic": None,
            "true": HDGTrue,
            "constantHeading": False,
        },
        magneticDeclinationGrid=magneticDeclinationGrid,
    )
ft.append(AC, flightTrajectory)


# descend on ILS with 3deg glideslope to next altitude threshold
# ------------------------------------------------
# current values
Hp, m_final, CAS_final, LAT_final, LON_final, HDGTrue = ft.getFinalValue(
    AC, ["Hp", "mass", "CAS", "LAT", "LON", "HDGTrue"]
)

flightTrajectory = TCL.constantSpeedSlope(
    AC=AC,
    speedType="CAS",
    v=CAS_final,
    Hp_init=Hp,
    Hp_final=Hp_RWY,
    slopetarget=-3.0,
    config=None,
    m_init=m_final,
    wS=wS,
    bankAngle=ba,
    DeltaTemp=DeltaTemp,
    Lat=LAT_final,
    Lon=LON_final,
    initialHeading={
        "magnetic": None,
        "true": HDGTrue,
        "constantHeading": False,
    },
    magneticDeclinationGrid=magneticDeclinationGrid,
)
ft.append(AC, flightTrajectory)

# print and plot final trajectory
df = ft.getFT(AC=AC)
print(df)

# Plotting the graph Hp=f(dist)
plt.figure(1, figsize=(8, 6))
plt.plot(df["dist"], df["Hp"], linestyle="-", color="b")
plt.grid(True)
plt.xlabel("Distance [NM]")
plt.ylabel("Pressure Altitude [ft]")
plt.title("Pressure Altitude as a Function of Distance")

# Plot for Calibrated Airspeed (CAS)
plt.figure(2, figsize=(8, 6))
plt.plot(df["dist"], df["CAS"], linestyle="-", color="r")
plt.grid(True)
plt.xlabel("Distance [NM]")
plt.ylabel("CAS [kt]")
plt.title("Calibrated Airspeed (CAS) as a Function of Distance")

# Display the plot
plt.show()


# save the output to a CSV/XLSX file
# ------------------------------------------------
# ft.save2csv(os.path.join(grandParentDir,"flightTrajectory_export"), separator=',')
# ft.save2xlsx(os.path.join(grandParentDir,"flightTrajectory_export"))
# ft.save2kml(os.path.join(grandParentDir,"flightTrajectory_export"))
