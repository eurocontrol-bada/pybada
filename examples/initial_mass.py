"""
Initial Mass Calculation
========================

Example calculation of aircraft initial mass

"""

from pyBADA.bada3 import Bada3Aircraft
from pyBADA import trajectoryPrediction as TP
from pyBADA import atmosphere as atm
from pyBADA import conversions as conv


# calculate estimations for the fuel flow, and aircraft initial mass
AC = Bada3Aircraft(badaVersion="DUMMY", acName="J2M")

deltaTemp = 0  # deviation from ISA temperature [K or deg C]
M = 0.7  # Mach number [-]
altitude = conv.ft2m(30000)  # cruise altitude[m]
distance = conv.nm2m(100)  # flown distance [NM]
payload = 80  # payload mass [% of the maximum payload mass]
fuelReserve = 3600  # reserve of the fuel [s]
flightPlanInitialMass = None  # planned initial/takeoff mass [kg]

# fuel flow in cruise
cruiseFuelFlow = TP.cruiseFuelConsumption(
    AC=AC, altitude=altitude, M=M, deltaTemp=deltaTemp
)

# in case of no wind, the ground speed is the same as true airspeed
[theta, delta, sigma] = atm.atmosphereProperties(
    h=altitude, DeltaTemp=deltaTemp
)
TAS = atm.mach2Tas(Mach=M, theta=theta)
GS = TAS

# distance based initial mass
breguetLeducInitialMass = TP.breguetLeducInitialMass(
    AC=AC,
    distance=distance,
    GS=GS,
    cruiseFuelFlow=cruiseFuelFlow,
    payload=payload,
    fuelReserve=fuelReserve,
)

# calculation of initial mass taking into account flight plan data and aircraft flight envelope
initMass = TP.getInitialMass(
    AC=AC,
    distance=distance,
    altitude=altitude,
    M=M,
    payload=payload,
    fuelReserve=fuelReserve,
    flightPlanInitialMass=flightPlanInitialMass,
    deltaTemp=deltaTemp,
)

print(f"cruiseFuelFlow: {cruiseFuelFlow} [kg/s]")
print(f"breguetLeducInitialMass: {breguetLeducInitialMass} [kg]")
print(f"initMass: {initMass} [kg]")
