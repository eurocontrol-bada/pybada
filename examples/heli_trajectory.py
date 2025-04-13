"""
Helicopter Trajectory Calculation
=================================

Example of BADAH trajectory using TCL
"""

from math import tan, pi
from dataclasses import dataclass
import matplotlib.pyplot as plt

from pyBADA import constants as const
from pyBADA import TCL as TCL
from pyBADA.flightTrajectory import FlightTrajectory as FT

from pyBADA.badaH import Parser as BadaHParser
from pyBADA.badaH import BadaHAircraft


@dataclass
class target:
    ROCDtarget: float = None
    slopetarget: float = None
    acctarget: float = None
    ESFtarget: float = None


# initialization of BADAH
badaVersion = "DUMMY"

allData = BadaHParser.parseAll(badaVersion=badaVersion)
print(allData)

AC = BadaHAircraft(badaVersion=badaVersion, acName="DUMH")

# create a Flight Trajectory object to store the output from TCL segment calculations
ft = FT()

# take-off with const ROCD
speedType = "TAS"  # {M, CAS, TAS}
v = 0  # [kt] CAS/TAS speed to follow or [-] MACH speed to follow
Hp_init = 0  # [ft]
Hp_VTOL = 5  # [ft] upper altitude for vertical take-off and landing
m_init = AC.OEW + 0.7 * (AC.MTOW - AC.OEW)  # [kg] initial mass
wS = 0  # [kt] wind speed
bankAngle = 0  # [deg] bank angle
DeltaTemp = 0  # [K] delta temperature from ISA
Hp_step = 500  # [ft] altitude step
step_length = 10  # iteration step length for cruise
RFL = 3000
maxRFL = 3000

AVT = 0.3 * const.g  # [g] to [m/s^2]
DVT = -0.3 * const.g  # [g] to [m/s^2]

MEC = AC.OPT.getOPTParam("MEC", var_1=0, var_2=m_init, DeltaTemp=DeltaTemp)
LRC = AC.OPT.getOPTParam("LRC", RFL, m_init, DeltaTemp)

flightTrajectory = TCL.constantSpeedROCD(
    AC=AC,
    speedType=speedType,
    v=v,
    Hp_init=Hp_init,
    Hp_final=Hp_VTOL,
    ROCDtarget=100,
    m_init=m_init,
    wS=wS,
    bankAngle=bankAngle,
    DeltaTemp=DeltaTemp,
    Hp_step=Hp_step,
)
ft.append(AC, flightTrajectory)

Hp, m_final = ft.getFinalValue(AC, ["Hp", "mass"])
control = target(ROCDtarget=500, acctarget=AVT)

# acc in climb
flightTrajectory = TCL.accDec(
    AC=AC,
    speedType=speedType,
    v_init=v,
    v_final=MEC,
    Hp_init=Hp,
    phase="Climb",
    control=control,
    maxRating="MTKF",
    m_init=m_final,
    wS=wS,
    bankAngle=bankAngle,
    DeltaTemp=DeltaTemp,
)
ft.append(AC, flightTrajectory)

Hp, m_final, v = ft.getFinalValue(AC, ["Hp", "mass", "TAS"])

# climb const ROCD
flightTrajectory = TCL.constantSpeedROCD(
    AC=AC,
    speedType=speedType,
    v=v,
    Hp_init=Hp,
    Hp_final=3000,
    ROCDtarget=1000,
    m_init=m_final,
    wS=wS,
    bankAngle=bankAngle,
    DeltaTemp=DeltaTemp,
    Hp_step=Hp_step,
)
ft.append(AC, flightTrajectory)

Hp, m_final, v = ft.getFinalValue(AC, ["Hp", "mass", "TAS"])

# acc in cruise
flightTrajectory = TCL.accDec(
    AC=AC,
    speedType=speedType,
    v_init=v,
    v_final=LRC,
    Hp_init=Hp,
    phase="Cruise",
    maxRating="MCNT",
    m_init=m_final,
    wS=wS,
    bankAngle=bankAngle,
    DeltaTemp=DeltaTemp,
)
ft.append(AC, flightTrajectory)

Hp, m_final, v = ft.getFinalValue(AC, ["Hp", "mass", "TAS"])

DEdist = RFL / tan(3 * pi / 180) * 0.3048 / 1852  # [NM]
length = 14.02  # 30 - 3.57 - DEdist

# cruise const TAS
flightTrajectory = TCL.constantSpeedLevel(
    AC=AC,
    lengthType="distance",
    length=length,
    speedType=speedType,
    v=v,
    Hp_init=Hp,
    m_init=m_final,
    maxRFL=maxRFL,
    wS=wS,
    bankAngle=bankAngle,
    DeltaTemp=DeltaTemp,
    step_length=step_length,
)
ft.append(AC, flightTrajectory)

Hp, m_final, v = ft.getFinalValue(AC, ["Hp", "mass", "TAS"])

# descent const ROCD
flightTrajectory = TCL.constantSpeedROCD(
    AC=AC,
    speedType=speedType,
    v=v,
    Hp_init=Hp,
    Hp_final=500,
    ROCDtarget=-500,
    m_init=m_final,
    wS=wS,
    bankAngle=bankAngle,
    DeltaTemp=DeltaTemp,
    Hp_step=Hp_step,
)
ft.append(AC, flightTrajectory)

Hp, m_final, v = ft.getFinalValue(AC, ["Hp", "mass", "TAS"])
control = target(ROCDtarget=-300, ESFtarget=0.3)

# dec in descent const ROCD
flightTrajectory = TCL.accDec(
    AC=AC,
    speedType=speedType,
    v_init=v,
    v_final=MEC,
    Hp_init=Hp,
    phase="Descent",
    control=control,
    m_init=m_final,
    wS=wS,
    bankAngle=bankAngle,
    DeltaTemp=DeltaTemp,
)
ft.append(AC, flightTrajectory)

Hp, m_final, v = ft.getFinalValue(AC, ["Hp", "mass", "TAS"])

# descent const ROCD
flightTrajectory = TCL.constantSpeedROCD(
    AC=AC,
    speedType=speedType,
    v=v,
    Hp_init=Hp,
    Hp_final=150,
    ROCDtarget=-300,
    m_init=m_final,
    wS=wS,
    bankAngle=bankAngle,
    DeltaTemp=DeltaTemp,
    Hp_step=Hp_step,
)
ft.append(AC, flightTrajectory)

Hp, m_final, v = ft.getFinalValue(AC, ["Hp", "mass", "TAS"])
control = target(ROCDtarget=-200, ESFtarget=0.3)

# dec in descent const ROCD
flightTrajectory = TCL.accDec(
    AC=AC,
    speedType=speedType,
    v_init=v,
    v_final=30,
    Hp_init=Hp,
    phase="Descent",
    control=control,
    m_init=m_final,
    wS=wS,
    bankAngle=bankAngle,
    DeltaTemp=DeltaTemp,
)
ft.append(AC, flightTrajectory)

Hp, m_final, v = ft.getFinalValue(AC, ["Hp", "mass", "TAS"])

# descent const ROCD
flightTrajectory = TCL.constantSpeedROCD(
    AC=AC,
    speedType=speedType,
    v=v,
    Hp_init=Hp,
    Hp_final=Hp_VTOL,
    ROCDtarget=-200,
    m_init=m_final,
    wS=wS,
    bankAngle=bankAngle,
    DeltaTemp=DeltaTemp,
    Hp_step=Hp_step,
)
ft.append(AC, flightTrajectory)

Hp, m_final, v = ft.getFinalValue(AC, ["Hp", "mass", "TAS"])
control = target(acctarget=DVT)

# dec in descent const ROCD
flightTrajectory = TCL.accDec(
    AC=AC,
    speedType=speedType,
    v_init=v,
    v_final=0,
    Hp_init=Hp,
    phase="Cruise",
    control=control,
    m_init=m_final,
    wS=wS,
    bankAngle=bankAngle,
    DeltaTemp=DeltaTemp,
)
ft.append(AC, flightTrajectory)

Hp, m_final, v = ft.getFinalValue(AC, ["Hp", "mass", "TAS"])

# descent const ROCD
flightTrajectory = TCL.constantSpeedROCD(
    AC=AC,
    speedType=speedType,
    v=v,
    Hp_init=Hp,
    Hp_final=0,
    ROCDtarget=-100,
    m_init=m_final,
    wS=wS,
    bankAngle=bankAngle,
    DeltaTemp=DeltaTemp,
    Hp_step=Hp_step,
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

# Plot for True Airspeed (TAS)
plt.figure(2, figsize=(8, 6))
plt.plot(df["dist"], df["TAS"], linestyle="-", color="r")
plt.grid(True)
plt.xlabel("Distance [NM]")
plt.ylabel("TAS [kt]")
plt.title("True Airspeed (TAS) as a Function of Distance")

# Display the plot
plt.show()

# save the output to a CSV/XLSX file
# ------------------------------------------------
# ft.save2csv(os.path.join(grandParentDir,"flightTrajectory_export"), separator=',')
# ft.save2xlsx(os.path.join(grandParentDir,"flightTrajectory_export"))
