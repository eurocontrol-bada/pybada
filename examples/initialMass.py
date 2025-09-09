"""
Initial Mass Calculation
========================

Example calculation of aircraft initial mass

"""

from pyBADA import atmosphere as atm
from pyBADA import conversions as conv
from pyBADA import trajectoryPrediction as TP
from pyBADA.bada3 import Bada3Aircraft
from pyBADA.bada4 import Bada4Aircraft
from pyBADA.badaH import BadaHAircraft

# calculate estimations for the fuel flow, and aircraft initial mass
AC = Bada3Aircraft(badaVersion="DUMMY", acName="J2M")

# Common inputs
deltaTemp = 0  # deviation from ISA temperature [K]
M = 0.7  # Mach number [-]
altitude = conv.ft2m(30000)  # cruise altitude [m]
distance = conv.nm2m(100)  # flown distance [m]
payload = 80  # payload mass [% of max payload]
fuelReserve = 3600  # fuel reserve [s]
flightPlanInitMass = None  # planned takeoff mass [kg]

# Precompute atmosphere & speed inputs
theta, delta, sigma = atm.atmosphereProperties(h=altitude, deltaTemp=deltaTemp)
TAS = atm.mach2Tas(Mach=M, theta=theta)
GS = TAS  # assume no wind

# Specify acName per model
for model_name, ACClass, acName in [
    ("BADA3", Bada3Aircraft, "J2M"),
    ("BADA4", Bada4Aircraft, "Dummy-TWIN-plus"),
    ("BADAH", BadaHAircraft, "DUMH"),
]:
    print(f"\n=== {model_name} ===")
    AC = ACClass(badaVersion="DUMMY", acName=acName)

    # 1) Cruise fuel flow [kg/s]
    cruiseFuelFlow = TP.cruiseFuelConsumption(
        AC=AC, altitude=altitude, M=M, deltaTemp=deltaTemp
    )
    print(f"cruiseFuelFlow:          {cruiseFuelFlow:.6f} kg/s")

    # 2) Distance based initial mass via Breguet–Leduc
    breguetInitialMass = TP.breguetLeducInitialMass(
        AC=AC,
        distance=distance,
        GS=GS,
        cruiseFuelFlow=cruiseFuelFlow,
        payload=payload,
        fuelReserve=fuelReserve,
    )
    print(f"breguetLeducInitialMass: {breguetInitialMass:.2f} kg")

    # 3) Initial mass using limited by flight envelope (OEW, MTOW)
    initMass = TP.getInitialMass(
        AC=AC,
        distance=distance,
        altitude=altitude,
        M=M,
        payload=payload,
        fuelReserve=fuelReserve,
        flightPlanInitialMass=flightPlanInitMass,
        deltaTemp=deltaTemp,
    )
    print(f"initMass(limited):       {initMass:.2f} kg")
