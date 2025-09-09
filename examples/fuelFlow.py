"""
Fuel Flow Calculation
========================

Example calculation of aircraft fuel flow in descent

"""

from math import sin

import matplotlib.pyplot as plt
import numpy as np

from pyBADA import atmosphere as atm
from pyBADA import constants as const
from pyBADA import conversions as conv
from pyBADA import utils
from pyBADA.bada3 import Bada3Aircraft
from pyBADA.bada4 import Bada4Aircraft

# create an aircraft
ACList = [
    Bada4Aircraft(badaVersion="DUMMY", acName="Dummy-TWIN"),
    Bada3Aircraft(badaVersion="DUMMY", acName="J2M"),
]

# deviation from ISA temperature
deltaTemp = 0

# definition of altitude range
fl_array = np.arange(0, 401, 10)
altitude_array = conv.ft2m(fl_array * 100)

# --- collect per-aircraft series here ---
series = []  # list of dicts: {"label": str, "alt_ft": list[float], "ff": list[float]}

for AC in ACList:
    # define aircraft mass - here as reference mass
    mass = AC.MREF

    # get the original speed schedule for descent for this aircraft
    [Vdes1, Vdes2, Mdes] = AC.flightEnvelope.getSpeedSchedule(phase="Descent")

    # crossover altitude
    crossAlt = atm.crossOver(cas=Vdes2, Mach=Mdes)

    label = f"({AC.BADAFamilyName}) {AC.acName.strip('_')} (BADA {getattr(AC, 'BADAVersion')}) {mass} kg"

    alt_ft_vals = []
    ff_vals = []

    print(
        f"\n=== {AC.__class__.__name__}  |  BADA {AC.BADAVersion}  |  Mass: {mass:.0f} kg ==="
    )
    print(
        f"{'FL':>4}  {'Alt(ft)':>8}  {'Cfg':>3}  {'M':>6}  {'FF (kg/s)':>10}"
    )
    print("-" * 39)

    for alt in altitude_array:
        # atmosphere properties
        theta, delta, sigma = atm.atmosphereProperties(
            h=alt, deltaTemp=deltaTemp
        )

        # determine the speed acording to BADA ARPM
        [cas, speedUpdated] = AC.ARPM.descentSpeed(
            h=alt, mass=mass, theta=theta, delta=delta, deltaTemp=deltaTemp
        )
        # general speed conversion
        [M, CAS, TAS] = atm.convertSpeed(
            v=conv.ms2kt(cas),
            speedType="CAS",
            theta=theta,
            delta=delta,
            sigma=sigma,
        )

        # determine the aerodynamic configuration if necesary
        config = AC.flightEnvelope.getConfig(
            h=alt, phase="Descent", v=CAS, mass=mass, deltaTemp=deltaTemp
        )

        # calculate Energy Share Factor depending if aircraft is flying constant M or CAS (based on crossover altitude)
        if alt < crossAlt:
            ESF = AC.esf(
                h=alt, flightEvolution="constCAS", M=M, deltaTemp=deltaTemp
            )
        else:
            ESF = AC.esf(
                h=alt, flightEvolution="constM", M=M, deltaTemp=deltaTemp
            )

        # =====
        # BADA4
        # =====
        if AC.BADAFamily.BADA4:
            # =================================================================================
            # for altitudes where aircraft descends on 3degree slope in AP and LD configuration
            # =================================================================================
            if config == "AP" or config == "LD":
                gamma = -3.0
                temp_const = (theta * const.temp_0) / (
                    theta * const.temp_0 - deltaTemp
                )
                ROCD_gamma = sin(conv.deg2rad(gamma)) * TAS * (1 / temp_const)

                n = 1.0
                [HLid, LG] = AC.flightEnvelope.getAeroConfig(config=config)
                CL = AC.CL(M=M, delta=delta, mass=mass, nz=n)
                CD = AC.CD(M=M, CL=CL, HLid=HLid, LG=LG)
                Drag = AC.D(M=M, delta=delta, CD=CD)
                Thrust = (ROCD_gamma * mass * const.g) * temp_const / (
                    ESF * TAS
                ) + Drag
                CT = AC.CT(Thrust=Thrust, delta=delta)
                ff = AC.ff(
                    CT=CT, delta=delta, theta=theta, M=M, deltaTemp=deltaTemp
                )

            # =============================================================
            # for altitudes where aircraft descends in IDLE engine settings
            # =============================================================
            else:
                ff = AC.ff(
                    rating="LIDL",
                    delta=delta,
                    theta=theta,
                    M=M,
                    deltaTemp=deltaTemp,
                )  # [kg/s]

        # =====
        # BADA3
        # =====
        elif AC.BADAFamily.BADA3:
            adaptedThrust = False
            if AC.engineType in ("PISTON", "ELECTRIC"):
                # PISTON  and ELECTRIC uses LIDL throughout the whole descent phase
                config = "CR"
                adaptedThrust = True

            # =================================================================================
            # for altitudes where aircraft descends on 3degree slope in AP and LD configuration
            # =================================================================================
            if config in ("AP", "LD"):
                gamma = -3.0
                temp_const = (theta * const.temp_0) / (
                    theta * const.temp_0 - deltaTemp
                )
                ROCD_gamma = sin(conv.deg2rad(gamma)) * TAS * (1 / temp_const)

                n = 1.0
                CL = AC.CL(sigma=sigma, mass=mass, tas=TAS, nz=n)
                CD = AC.CD(CL=CL, config=config)
                Drag = AC.D(sigma=sigma, tas=TAS, CD=CD)
                Thrust = AC.Thrust(
                    rating="ADAPTED",
                    v=TAS,
                    config=config,
                    h=alt,
                    ROCD=ROCD_gamma,
                    mass=mass,
                    acc=0,
                    deltaTemp=deltaTemp,
                    Drag=Drag,
                )
                ff = AC.ff(
                    flightPhase="Descent",
                    v=TAS,
                    h=alt,
                    T=Thrust,
                    config=config,
                    adapted=adaptedThrust,
                )

            # =============================================================
            # for altitudes where aircraft descends in IDLE engine settings
            # =============================================================
            else:
                ff = AC.ff(
                    v=TAS,
                    h=alt,
                    T=Thrust,
                    flightPhase="Descent",
                    config=config,
                )

        fl = int(utils.proper_round(conv.m2ft(alt) / 100))
        alt_ft = conv.m2ft(alt)
        print(f"{fl:>4d}  {alt_ft:>8.0f}  {config:>3}  {M:>6.3f}  {ff:>10.6f}")

        alt_ft_vals.append(alt_ft)
        ff_vals.append(float(ff))

    series.append(
        {
            "label": label,
            "alt_ft": np.array(alt_ft_vals),
            "ff": np.array(ff_vals),
        }
    )

# --- PLOT: Fuel flow vs Altitude for all aircraft ---
plt.figure(figsize=(8, 5))
for s in series:
    plt.plot(s["alt_ft"], s["ff"], label=s["label"])

plt.xlabel("Altitude (ft)")
plt.ylabel("Fuel Flow (kg/s)")
plt.title("Descent Fuel Flow vs Altitude")
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.show()
