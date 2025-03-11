"""
Initial Mass Calculation
========================

Example of calculation of aircraft initial mass

"""

from pyBADA.bada3 import Bada3Aircraft
from pyBADA import trajectoryPrediction as TP
from pyBADA import atmosphere as atm
from pyBADA import conversions as conv


# calculate estimations for the fuel flow, and aircraft initial mass
AC = Bada3Aircraft(badaVersion="3.16", acName="A320")

deltaTemp = 0
M = 0.7
altitude = conv.ft2m(30000)
distance = conv.nm2m(100)
payload = 80
fuelReserve = 3600
flightPlanInitialMass = None

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

# distance based initiall mass
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

print("cruiseFuelFlow:", cruiseFuelFlow, "[kg/s]")
print("breguetLeducInitialMass:", breguetLeducInitialMass, "[kg]")
print("initMass:", initMass, "[kg]")
