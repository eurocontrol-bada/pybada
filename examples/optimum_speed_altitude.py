"""
Optimisation Calculation
========================

Example of BADA4 and BADAH optimum speed and altitude calculation
"""

from pyBADA import atmosphere as atm
from pyBADA import conversions as conv

from pyBADA.bada4 import Bada4Aircraft
from pyBADA.badaH import BadaHAircraft

badaVersion = "DUMMY"

# initialization of BADA4
AC = Bada4Aircraft(badaVersion=badaVersion, acName="Dummy-TWIN")
print("BADA4 Optimum Speed and Altitude:")

# BADA4
if AC.BADAFamily.BADA4:
    mass = AC.MTOW  # [kg] AC weight
    h = conv.ft2m(33000)  # [m] AC flight altitdue
    DeltaTemp = 0  # [K] temperature deviation from ISA
    cI = 50  # [kg min^-1] cost index
    wS = 0  # [m s^-1] longitudinal wind speed

    [theta, delta, sigma] = atm.atmosphereProperties(
        h=h, DeltaTemp=DeltaTemp
    )  # atmosphere properties

    # Economic Mach Cruise Speed
    econMach = AC.OPT.econMach(
        theta=theta, delta=delta, mass=mass, DeltaTemp=DeltaTemp, cI=cI, wS=wS
    )
    print("EconMach = ", econMach)

    CCI = AC.OPT.CCI(theta=theta, delta=delta, cI=cI)
    CW = AC.OPT.CW(mass=mass, delta=delta)
    econMach_precomputed = AC.OPT.getOPTParam("ECON", CW, CCI)
    print("EconMach_precomputed = ", econMach_precomputed)

    # Maximum Range Cruise (MRC) Mach speed
    MRC = AC.OPT.MRC(
        theta=theta, delta=delta, mass=mass, DeltaTemp=DeltaTemp, wS=wS
    )
    print("MRC = ", MRC)

    MRC_precomputed = AC.OPT.getOPTParam("MRC", CW)
    print("MRC_precomputed = ", MRC_precomputed)

    # Long Range Cruise (LRC) Mach speed
    LRC = AC.OPT.LRC(
        theta=theta, delta=delta, mass=mass, DeltaTemp=DeltaTemp, wS=wS
    )
    print("LRC = ", LRC)

    LRC_precomputed = AC.OPT.getOPTParam("LRC", CW)
    print("LRC_precomputed = ", LRC_precomputed)

    # Maximum Endurance Cruise (MEC) Mach speed
    MEC = AC.OPT.MEC(
        theta=theta, delta=delta, mass=mass, DeltaTemp=DeltaTemp, wS=wS
    )
    print("MEC = ", MEC)

    MEC_precomputed = AC.OPT.getOPTParam("MEC", CW)
    print("MEC_precomputed = ", MEC_precomputed)

    # optimum flight altitude at given M speed
    M = MRC
    optAlt = AC.OPT.optAltitude(M=M, mass=mass, DeltaTemp=DeltaTemp)
    print("optAlt =", optAlt)

    optAlt_precomputed = AC.OPT.getOPTParam("OPTALT", M, mass)
    print("optAlt_precomputed = ", optAlt_precomputed)

# initialization of BADAH
AC = BadaHAircraft(badaVersion=badaVersion, acName="DUMH")
print("\nBADAH Optimum Speed and Altitude:")

# BADAH
if AC.BADAFamily.BADAH:
    mass = 1600  # [kg] AC weight
    Hp = 14000  # [ft] AC flight altitude
    h = conv.ft2m(Hp)  # [m] AC flight altitude
    DeltaTemp = 20  # [K] temperature deviation from ISA
    wS = 0  # [m s^-1] longitudinal wind speed

    [theta, delta, sigma] = atm.atmosphereProperties(
        h=h, DeltaTemp=DeltaTemp
    )  # atmoshpere properties

    # Maximum Range Cruise (MRC) Mach speed
    MRC = AC.OPT.MRC(h=h, mass=mass, DeltaTemp=DeltaTemp, wS=wS)
    # print("MRC = ", conv.ms2kt(MRC))
    print("MRC = ", (MRC))

    MRC_precomputed = AC.OPT.getOPTParam("MRC", Hp, mass, DeltaTemp)
    print("MRC_precomputed = ", conv.kt2ms(MRC_precomputed))

    # Long Range Cruise (LRC) Mach speed
    LRC = AC.OPT.LRC(h=h, mass=mass, DeltaTemp=DeltaTemp, wS=wS)
    print("LRC = ", (LRC))

    LRC_precomputed = AC.OPT.getOPTParam("LRC", Hp, mass, DeltaTemp)
    print("LRC_precomputed = ", conv.kt2ms(LRC_precomputed))

    # Maximum Endurance Cruise (MEC) Mach speed
    MEC = AC.OPT.MEC(h=h, mass=mass, DeltaTemp=DeltaTemp, wS=wS)
    print("MEC = ", (MEC))

    MEC_precomputed = AC.OPT.getOPTParam("MEC", Hp, mass, DeltaTemp)
    print("MEC_precomputed = ", conv.kt2ms(MEC_precomputed))
