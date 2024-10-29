"""
BADA File Parser
================

Example of BADA3, BADA4 and BADAH file parser
"""

from pyBADA.bada3 import Bada3Aircraft
from pyBADA.bada4 import Bada4Aircraft
from pyBADA.badaH import BadaHAircraft

# initialization of BADAH
AC = BadaHAircraft(badaVersion="DUMMY", acName="DUMH")

# BADAH
if AC.BADAFamily.BADAH:
    ICAO = AC.ICAO
    WTC = AC.WTC
    MTOW = AC.MTOW
    print(
        "BADA Family:",
        AC.BADAFamilyName,
        "| BADA Version:",
        AC.BADAVersion,
        "| ICAO:",
        ICAO,
        "| WTC:",
        WTC,
        "| MTOW =",
        MTOW,
    )

# initialization of BADA4
AC = Bada4Aircraft(badaVersion="DUMMY", acName="Dummy-TBP")

# BADA3 or BADA4
if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
    ICAO = AC.ICAO
    WTC = AC.WTC
    VMO = AC.VMO
    MMO = AC.MMO
    MTOW = AC.MTOW
    print(
        "BADA Family:",
        AC.BADAFamilyName,
        "| BADA Version:",
        AC.BADAVersion,
        "| ICAO:",
        ICAO,
        "| WTC:",
        WTC,
        "| VMO =",
        VMO,
        "| MMO =",
        MMO,
        "| MTOW =",
        MTOW,
    )

# initialization of BADA3
AC = Bada3Aircraft(badaVersion="DUMMY", acName="BZJT")

# BADA3 or BADA4
if AC.BADAFamily.BADA3 or AC.BADAFamily.BADA4:
    ICAO = AC.ICAO
    WTC = AC.WTC
    VMO = AC.VMO
    MMO = AC.MMO
    MTOW = AC.MTOW
    print(
        "BADA Family:",
        AC.BADAFamilyName,
        "| BADA Version:",
        AC.BADAVersion,
        "| ICAO:",
        ICAO,
        "| WTC:",
        WTC,
        "| VMO =",
        VMO,
        "| MMO =",
        MMO,
        "| MTOW =",
        MTOW,
    )
