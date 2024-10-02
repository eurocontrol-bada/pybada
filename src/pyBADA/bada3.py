# -*- coding: utf-8 -*-
"""
pyBADA
Generic BADA3 aircraft performance module
Developped @EUROCONTROL (EIH)
2024
"""

__author__ = "Henrich Glaser-Opitz"
__copyright__ = "Copyright 2024, EUROCONTROL (EIH)"
__license__ = "BADA Eurocontrol"
__version__ = "1.0.0"
__maintainer__ = "Henrich Glaser-Opitz"
__email__ = "henrich.glaser-opitz@eurocontrol.int"
__status__ = "Development"
__docformat__ = "reStructuredText"

from math import sqrt, isnan, asin, atan
import numpy as np

import os
from datetime import date
import xml.etree.ElementTree as ET
import pandas as pd

from pyBADA import constants as const
from pyBADA import conversions as conv
from pyBADA import atmosphere as atm
from pyBADA import configuration as configuration
from pyBADA.aircraft import Airplane, BadaFamily


def proper_round(num, dec=0):
    num = str(num)[: str(num).index(".") + dec + 2]
    if num[-1] >= "5":
        return float(num[: -2 - (not dec)] + str(int(num[-2 - (not dec)]) + 1))
    return float(num[:-1])


def checkArgument(argument, **kwargs):
    if kwargs.get(argument) is not None:
        return kwargs.get(argument)
    else:
        raise TypeError("Missing " + argument + " argument")


class Parser(object):
    """This class implements the BADA3 parsing mechanism to parse APF, OPF and GPF BADA3 files."""

    def __init__(self):
        pass

    @staticmethod
    def list_subfolders(folderPath):
        # List all entries in the directory
        entries = os.listdir(folderPath)

        # Filter out entries that are directories
        subfolders = [
            entry for entry in entries if os.path.isdir(os.path.join(folderPath, entry))
        ]

        return subfolders

    @staticmethod
    def parseXML(filePath, badaVersion, acName):
        """This function parses BADA3 xml formatted file

        :param filePath: path to the folder with BADA4 xml formatted file.
        :param acName: name of Aircraft for BADA3 xml formatted file.
        :type filePath: str.
        :type acName: str.
        :raises: IOError
        """

        filename = os.path.join(filePath, "BADA3", badaVersion, acName, acName) + ".xml"

        try:
            tree = ET.parse(filename)
            root = tree.getroot()
        except:
            raise IOError(filename + " not found or in correct format")

        modificationDateOPF = "UNKNOWN"
        modificationDateAPF = "UNKNOWN"

        # Parse general aircraft data
        model = root.find("model").text  # aircraft model
        engineType = root.find("type").text  # engine type
        engines = root.find("engine").text  # engine name

        ICAO = root.find("ICAO").find("designator").text
        WTC = root.find("ICAO").find("WTC").text

        # Parse engine data
        AFCM = root.find("AFCM")  # get AFCM
        PFM = root.find("PFM")  # get PFM
        ALM = root.find("ALM")  # get ALM
        Ground = root.find("Ground")  # get Ground
        ARPM = root.find("ARPM")  # get ARPM

        # AFCM
        S = float(AFCM.find("S").text)
        MREF = float(AFCM.find("mref").text)

        mass = {}
        mass["reference"] = float(AFCM.find("mref").text)

        name = {}
        HLids = []
        d = {}
        CD0 = {}
        CD2 = {}
        Vstall = {}

        for conf in AFCM.findall("Configuration"):
            HLid = int(conf.get("HLid"))
            HLids.append(str(HLid))
            name[HLid] = conf.find("name").text

            d[HLid] = {}
            CD0[HLid] = {}
            CD2[HLid] = {}
            Vstall[HLid] = {}

            LGUP = conf.find("LGUP")
            LGDN = conf.find("LGDN")

            if LGUP is not None:
                DPM = LGUP.find("DPM")

                d[HLid]["LGUP"] = []
                for i in DPM.find("CD").findall("d"):
                    d[HLid]["LGUP"].append(float(i.text))

                CD0[HLid]["LGUP"] = d[HLid]["LGUP"][0]
                CD2[HLid]["LGUP"] = d[HLid]["LGUP"][1]

                BLM = LGUP.find("BLM")

                Vstall[HLid]["LGUP"] = []
                if BLM is not None:  # BLM is not clean
                    Vstall[HLid]["LGUP"] = float(BLM.find("VS").text)

                else:  # BLM is clean
                    BLM = LGUP.find("BLM_clean")

                    Vstall[HLid]["LGUP"] = float(BLM.find("VS").text)

                    CL_clean = BLM.find("CL_clean")

                    Clbo = float(CL_clean.find("Clbo").text)
                    k = float(CL_clean.find("k").text)

            if LGDN is not None:  # Landing gear NOT allowed in clean configuration

                d[HLid]["LGDN"] = []
                for i in LGDN.find("DPM").find("CD").findall("d"):
                    d[HLid]["LGDN"].append(float(i.text))

                CD0[HLid]["LGDN"] = d[HLid]["LGDN"][0]
                CD2[HLid]["LGDN"] = d[HLid]["LGDN"][1]

                if LGDN.find("DPM").find("DeltaCD") is None:
                    DeltaCD = 0.0
                else:
                    DeltaCD = float(LGDN.find("DPM").find("DeltaCD").text)

                Vstall[HLid]["LGDN"] = float(LGDN.find("BLM").find("VS").text)

            elif LGDN is None:
                CD0[HLid]["LGDN"] = 0.0
                CD2[HLid]["LGDN"] = 0.0
                DeltaCD = 0.0

        drone = False
        if Vstall[0]["LGUP"] == 0.0:
            drone = True

        # PFM
        numberOfEngines = float(PFM.find("n_eng").text)

        CT = PFM.find("CT")
        CTc1 = float(CT.find("CTc1").text)
        CTc2 = float(CT.find("CTc2").text)
        CTc3 = float(CT.find("CTc3").text)
        CTc4 = float(CT.find("CTc4").text)
        CTc5 = float(CT.find("CTc5").text)
        Ct = [CTc1, CTc2, CTc3, CTc4, CTc5]

        CTdeslow = float(CT.find("CTdeslow").text)
        CTdeshigh = float(CT.find("CTdeshigh").text)
        CTdesapp = float(CT.find("CTdesapp").text)
        CTdesld = float(CT.find("CTdesld").text)
        HpDes = float(CT.find("Hpdes").text)

        CF = PFM.find("CF")

        Cf1 = float(CF.find("Cf1").text)
        Cf2 = float(CF.find("Cf2").text)
        Cf3 = float(CF.find("Cf3").text)
        Cf4 = float(CF.find("Cf4").text)
        Cfcr = float(CF.find("Cfcr").text)

        CfDes = [Cf3, Cf4]
        CfCrz = float(CF.find("Cfcr").text)
        Cf = [Cf1, Cf2]

        # ALM
        GLM = ALM.find("GLM")
        hmo = float(GLM.find("hmo").text)
        Hmax = float(GLM.find("hmax").text)
        tempGrad = float(GLM.find("temp_grad").text)
        massGrad = float(GLM.find("mass_grad").text)
        mass["mass grad"] = float(GLM.find("mass_grad").text)

        KLM = ALM.find("KLM")
        MMO = float(KLM.find("mmo").text)
        VMO = float(KLM.find("vmo").text)

        DLM = ALM.find("DLM")
        MTOW = float(DLM.find("MTOW").text)
        OEW = float(DLM.find("OEW").text)
        MPL = float(DLM.find("MPL").text)

        mass["minimum"] = float(DLM.find("OEW").text)
        mass["maximum"] = float(DLM.find("MTOW").text)
        mass["max payload"] = float(DLM.find("MPL").text)

        # Ground
        dimensions = Ground.find("Dimensions")
        Runway = Ground.find("Runway")
        TOL = float(Runway.find("TOL").text)
        LDL = float(Runway.find("LDL").text)
        span = float(dimensions.find("span").text)
        length = float(dimensions.find("length").text)

        # ARPM
        aeroConfSchedule = ARPM.find("AeroConfSchedule")

        # all aerodynamic configurations
        aeroConfig = {}
        for conf in aeroConfSchedule.findall("AeroPhase"):
            name = conf.find("name").text
            HLid = int(conf.find("HLid").text)
            LG = "LG" + conf.find("LG").text
            aeroConfig[name] = {"name": name, "HLid": HLid, "LG": LG}

        speedScheduleList = ARPM.find("SpeedScheduleList")
        SpeedSchedule = speedScheduleList.find("SpeedSchedule")

        # all phases of flight
        speedSchedule = {}
        for phaseOfFlight in SpeedSchedule.findall("SpeedPhase"):
            name = phaseOfFlight.find("name").text
            CAS1 = conv.kt2ms(float(phaseOfFlight.find("CAS1").text))
            CAS2 = conv.kt2ms(float(phaseOfFlight.find("CAS2").text))
            M = float(phaseOfFlight.find("M").text)
            speedSchedule[name] = {"CAS1": CAS1, "CAS2": CAS2, "M": M}

        V1 = {}
        V1["cl"] = speedSchedule["Climb"]["CAS1"]
        V1["cr"] = speedSchedule["Cruise"]["CAS1"]
        V1["des"] = speedSchedule["Descent"]["CAS1"]

        V2 = {}
        V2["cl"] = speedSchedule["Climb"]["CAS2"]
        V2["cr"] = speedSchedule["Cruise"]["CAS2"]
        V2["des"] = speedSchedule["Descent"]["CAS2"]

        M = {}
        M["cl"] = speedSchedule["Climb"]["M"]
        M["cr"] = speedSchedule["Cruise"]["M"]
        M["des"] = speedSchedule["Descent"]["M"]

        xmlFiles = True

        # Single row dataframe
        data = {
            "acName": [acName],
            "model": [model],
            "engineType": [engineType],
            "engines": [engines],
            "ICAO": [ICAO],
            "WTC": [WTC],
            "modificationDateOPF": [modificationDateOPF],
            "modificationDateAPF": [modificationDateAPF],
            "S": [S],
            "MREF": [MREF],
            "mass": [mass],
            "name": [name],
            "HLids": [HLids],
            "d": [d],
            "CD0": [CD0],
            "CD2": [CD2],
            "Vstall": [Vstall],
            "Clbo": [Clbo],
            "k": [k],
            "DeltaCD": [DeltaCD],
            "drone": [drone],
            "numberOfEngines": [numberOfEngines],
            "CTc1": [CTc1],
            "CTc2": [CTc2],
            "CTc3": [CTc3],
            "CTc4": [CTc4],
            "CTc5": [CTc5],
            "Ct": [Ct],
            "CTdeslow": [CTdeslow],
            "CTdeshigh": [CTdeshigh],
            "CTdeshigh": [CTdeshigh],
            "CTdesapp": [CTdesapp],
            "CTdesld": [CTdesld],
            "HpDes": [HpDes],
            "Cf1": [Cf1],
            "Cf2": [Cf2],
            "Cf3": [Cf3],
            "Cf4": [Cf4],
            "Cfcr": [Cfcr],
            "CfDes": [CfDes],
            "CfCrz": [CfCrz],
            "Cf": [Cf],
            "hmo": [hmo],
            "Hmax": [Hmax],
            "tempGrad": [tempGrad],
            "massGrad": [massGrad],
            "mass": [mass],
            "MMO": [MMO],
            "VMO": [VMO],
            "MTOW": [MTOW],
            "OEW": [OEW],
            "MPL": [MPL],
            "TOL": [TOL],
            "LDL": [LDL],
            "span": [span],
            "length": [length],
            "span": [span],
            "length": [length],
            "aeroConfig": [aeroConfig],
            "speedSchedule": [speedSchedule],
            "V1": [V1],
            "V2": [V2],
            "M": [M],
            "speedSchedule": [speedSchedule],
            "xmlFiles": [xmlFiles],
        }
        df_single = pd.DataFrame(data)

        return df_single

    @staticmethod
    def findData(f):
        line = f.readline()
        while line is not None and not line.startswith("CD"):
            line = f.readline()

        if line is None:
            return f, None

        line = " ".join(line.split())
        line = line.strip().split(" ")
        return f, line

    @staticmethod
    def parseOPF(filePath, badaVersion, acName):
        """This function parses BADA3 ascii formatted file

        :param filePath: path to the BADA3 ascii formatted file.
        :param acName: ICAO aircraft designation
        :type filePath: str.
        :type acName: str
        :raises: IOError
        """

        filename = (
            os.path.join(
                filePath,
                "BADA3",
                badaVersion,
                acName,
            )
            + ".OPF"
        )

        idx = 0
        with open(filename, "r", encoding="latin-1") as f:
            while True:
                line = f.readline()

                if idx == 13:
                    if "with" in line:
                        engines = line.split("with")[1].split("engines")[0].strip()
                    else:
                        engines = "unknown"
                idx += 1

                if not line:
                    break
                elif "Modification_date" in line:
                    data = line.split(":")[1].strip().split(" ")
                    modificationDateOPF = " ".join([data[0], data[1], data[2]])

                elif "CC====== Actype" in line:
                    f, line = Parser.findData(f=f)
                    if line is None:
                        break
                    ICAO = line[1].replace("_", "")
                    numberOfEngines = int(line[2])
                    engineType = line[4].upper()
                    WTC = line[5]

                elif "CC====== Mass (t)" in line:
                    f, line = Parser.findData(f=f)
                    if line is None:
                        break
                    mass = {}
                    MREF = float(line[1]) * 1000.0
                    mass["reference"] = float(line[1]) * 1000.0
                    mass["minimum"] = float(line[2]) * 1000.0
                    mass["maximum"] = float(line[3]) * 1000.0
                    mass["max payload"] = float(line[4]) * 1000.0
                    mass["mass grad"] = float(line[5])

                    MTOW = mass["maximum"]
                    OEW = mass["minimum"]
                    MPL = mass["max payload"]
                    massGrad = mass["mass grad"]

                elif "CC====== Flight envelope" in line:
                    f, line = Parser.findData(f=f)
                    if line is None:
                        break
                    VMO = float(line[1])
                    MMO = float(line[2])
                    hmo = float(line[3])
                    Hmax = float(line[4])
                    tempGrad = float(line[5])

                elif "CC====== Aerodynamics" in line:
                    f, line = Parser.findData(f=f)
                    if line is None:
                        break
                    ndrst = int(line[1])
                    S = float(line[2])
                    Clbo = float(line[3])
                    k = float(line[4])

                    n = 1
                    Vstall = {}
                    CD0 = {}
                    CD2 = {}
                    HLids = []
                    while n <= ndrst:
                        f, line = Parser.findData(f=f)
                        if line is None:
                            break
                        HLid = line[2]

                        Vstall[HLid] = float(line[-5])
                        CD0[HLid] = float(line[-4])
                        CD2[HLid] = float(line[-3])
                        HLids.append(str(HLid))
                        n += 1

                    drone = False
                    if Vstall["CR"] == 0.0:
                        drone = True

                    iterator = 1
                    while iterator <= 2:
                        f, line = Parser.findData(f=f)
                        if "EXT" in line[2]:
                            CD2["SPOILER_EXT"] = float(line[3])
                        iterator += 1

                    iterator = 1
                    while iterator <= 2:
                        f, line = Parser.findData(f=f)
                        if "DOWN" in line[2]:
                            CD0["GEAR_DOWN"] = float(line[3])
                            CD2["GEAR_DOWN"] = float(line[4])
                        iterator += 1

                    iterator = 1
                    while iterator <= 2:
                        f, line = Parser.findData(f=f)
                        if "ON" in line[2]:
                            CD2["BRAKES_ON"] = float(line[3])
                        iterator += 1

                elif "CC====== Engine Thrust" in line:
                    f, line = Parser.findData(f=f)
                    if line is None:
                        break
                    Ct = [float(i) for i in line[1:-1]]
                    f, line = Parser.findData(f=f)
                    if line is None:
                        break

                    CTdeslow = float(line[1])
                    CTdeshigh = float(line[2])
                    CTdesapp = float(line[4])
                    CTdesld = float(line[5])
                    HpDes = float(line[3])

                    # self.CtDes = {}
                    # self.CtDes["low"] = float(line[1])
                    # self.CtDes["high"] = float(line[2])
                    # self.HpDes = float(line[3])
                    # self.CtDes["app"] = float(line[4])
                    # self.CtDes["lnd"] = float(line[5])
                    f, line = Parser.findData(f=f)
                    if line is None:
                        break

                elif "CC====== Fuel Consumption" in line:
                    f, line = Parser.findData(f=f)
                    if line is None:
                        break
                    Cf = [float(i) for i in line[1:-1]]
                    f, line = Parser.findData(f=f)
                    if line is None:
                        break
                    CfDes = [float(i) for i in line[1:-1]]
                    f, line = Parser.findData(f=f)
                    if line is None:
                        break
                    CfCrz = float(line[1])

                elif "CC====== Ground" in line:
                    f, line = Parser.findData(f=f)
                    if line is None:
                        break
                    TOL = float(line[1])
                    LDL = float(line[2])
                    span = float(line[3])
                    length = float(line[4])

        # Single row dataframe
        data = {
            "acName": [acName],
            "engineType": [engineType],
            "engines": [engines],
            "ICAO": [ICAO],
            "WTC": [WTC],
            "modificationDateOPF": [modificationDateOPF],
            "modificationDateOPF": [modificationDateOPF],
            "S": [S],
            "MREF": [MREF],
            "mass": [mass],
            "HLids": [HLids],
            "CD0": [CD0],
            "CD2": [CD2],
            "Vstall": [Vstall],
            "Clbo": [Clbo],
            "k": [k],
            "drone": [drone],
            "numberOfEngines": [numberOfEngines],
            "Ct": [Ct],
            "CTdeslow": [CTdeslow],
            "CTdeshigh": [CTdeshigh],
            "CTdeshigh": [CTdeshigh],
            "CTdesapp": [CTdesapp],
            "CTdesld": [CTdesld],
            "HpDes": [HpDes],
            "CfDes": [CfDes],
            "CfCrz": [CfCrz],
            "Cf": [Cf],
            "hmo": [hmo],
            "Hmax": [Hmax],
            "tempGrad": [tempGrad],
            "mass": [mass],
            "MMO": [MMO],
            "VMO": [VMO],
            "massGrad": [massGrad],
            "MTOW": [MTOW],
            "OEW": [OEW],
            "MPL": [MPL],
            "TOL": [TOL],
            "LDL": [LDL],
            "span": [span],
            "length": [length],
            "span": [span],
            "length": [length],
        }
        df_single = pd.DataFrame(data)

        return df_single

    @staticmethod
    def parseAPF(filePath, badaVersion, acName):
        """This function parses BADA3 APF ascii formatted file

        :param filePath: path to the BADA3 APF ascii formatted file.
        :param acName: ICAO aircraft designation
        :type filePath: str.
        :type acName: str
        :raises: IOError
        """

        filename = os.path.join(filePath, "BADA3", badaVersion, acName) + ".APF"

        dataLines = list()
        with open(filename, "r", encoding="latin-1") as f:
            while True:
                line = f.readline()

                if line.startswith("CC"):
                    if "Modification_date" in line:
                        data = line.split(":")[1].strip().split(" ")
                        modificationDateAPF = " ".join([data[0], data[1], data[2]])
                if line.startswith("CD"):
                    line = " ".join(line.split())
                    line = line.strip().split(" ")

                    if "LO" in line:
                        line = line[line.index("LO") + 1 :]
                    elif "AV" in line:
                        line = line[line.index("AV") + 1 :]
                    elif "HI" in line:
                        line = line[line.index("HI") + 1 :]

                    dataLines.append(line)
                elif "THE END" in line:
                    break
        dataLines.pop(0)  # remove first line that does not contain usefull data

        # AV - average - line with average data
        AVLine = dataLines[1]
        # reading of V1 parameter from APF file

        V1 = {}
        V1["cl"] = conv.kt2ms(AVLine[0])
        V1["cr"] = conv.kt2ms(AVLine[3])
        V1["des"] = conv.kt2ms(AVLine[8])

        V2 = {}
        V2["cl"] = conv.kt2ms(AVLine[1])
        V2["cr"] = conv.kt2ms(AVLine[4])
        V2["des"] = conv.kt2ms(AVLine[7])

        M = {}
        M["cl"] = float(AVLine[2]) / 100
        M["cr"] = float(AVLine[5]) / 100
        M["des"] = float(AVLine[6]) / 100

        # Single row dataframe
        data = {
            "modificationDateAPF": [modificationDateAPF],
            "V1": [V1],
            "V2": [V2],
            "M": [M],
        }
        df_single = pd.DataFrame(data)

        return df_single

    @staticmethod
    def combineOPF_APF(OPFDataFrame, APFDataFrame):

        # Combine data with GPF data (temporary solution)
        combined_df = pd.concat(
            [OPFDataFrame.reset_index(drop=True), APFDataFrame.reset_index(drop=True)],
            axis=1,
        )

        return combined_df

    @staticmethod
    def readSynonym(filePath, badaVersion):

        filename = os.path.join(filePath, "BADA3", badaVersion, "SYNONYM.NEW")

        # synonym - file name pair dictionary
        synonym_fileName = {}

        if os.path.isfile(filename):
            with open(filename, "r", encoding="latin-1") as f:
                while True:
                    line = f.readline()

                    if not line:
                        break

                    if line.startswith("CD"):
                        line = " ".join(line.split())
                        line = line.strip().split(" ")

                        model = str(line[2])
                        file = str(line[-3])

                        synonym_fileName[model] = file

        return synonym_fileName

    @staticmethod
    def readSynonymXML(filePath, badaVersion):

        filename = os.path.join(filePath, "BADA3", badaVersion, "SYNONYM.xml")

        # synonym - file name pair dictionary
        synonym_fileName = {}

        if os.path.isfile(filename):
            try:
                tree = ET.parse(filename)
                root = tree.getroot()
            except:
                raise IOError(filename + " not found or in correct format")

            for child in root.iter("SYN"):
                code = child.find("code").text
                manufacturer = child.find("manu").text
                file = child.find("file").text
                ICAO = child.find("ICAO").text

                synonym_fileName[code] = file

        return synonym_fileName

    @staticmethod
    def parseSynonym(filePath, badaVersion, acName):

        synonym_fileName = Parser.readSynonym(filePath, badaVersion)

        # if ASCI synonym does not exist, try XML synonym file
        if not synonym_fileName:
            synonym_fileName = Parser.readSynonymXML(filePath, badaVersion)

        if synonym_fileName and acName in synonym_fileName:
            fileName = synonym_fileName[acName]
            return fileName
        else:
            return None

    @staticmethod
    def readGPF(filePath, badaVersion):
        """This function parses BADA3 GPF ascii formatted file

        :param filePath: path to the BADA3 GPF ascii formatted file.
        :type filePath: str.
        :raises: IOError
        """

        filename = os.path.join(filePath, "BADA3", badaVersion, "BADA.GPF")

        GPFparamList = list()

        if os.path.isfile(filename):
            with open(filename, "r", encoding="latin-1") as f:
                while True:
                    line = f.readline()

                    if not line:
                        break

                    if line.startswith("CD"):
                        line = " ".join(line.split())
                        line = line.strip().split(" ")
                        param = {
                            "name": str(line[1]),
                            "value": float(line[5]),
                            "engine": str(line[3]).split(","),
                            "phase": str(line[4]).split(","),
                            "flight": str(line[2]),
                        }

                        GPFparamList.append(param)
        return GPFparamList

    @staticmethod
    def readGPFXML(filePath, badaVersion):
        """This function parses BADA3 GPF xml formatted file

        :param filePath: path to the GPF xml formatted file.
        :type filePath: str.
        :raises: IOError
        """

        filename = os.path.join(filePath, "BADA3", badaVersion, "GPF.xml")

        GPFparamList = list()

        if os.path.isfile(filename):
            try:
                tree = ET.parse(filename)
                root = tree.getroot()
            except:
                raise IOError(filename + " not found or in correct format")

            allEngines = ["JET", "TURBOPROP", "PISTON", "ELECTRIC"]
            allPhases = ["to", "ic", "cl", "cr", "des", "hold", "app", "lnd"]
            allFlights = ["civ", "mil"]

            # Parse general aircraft data
            AccMax = root.find("AccMax")
            GPFparamList.append(
                {
                    "name": "acc_long_max",
                    "value": float(AccMax.find("long").text),
                    "engine": allEngines,
                    "phase": allPhases,
                    "flight": allFlights,
                }
            )
            GPFparamList.append(
                {
                    "name": "acc_norm_max",
                    "value": float(AccMax.find("norm").text),
                    "engine": allEngines,
                    "phase": allPhases,
                    "flight": allFlights,
                }
            )

            AngBank = root.find("AngBank")
            GPFparamList.append(
                {
                    "name": "ang_bank_nom",
                    "value": float(AngBank.find("Nom").find("Civ").find("ToLd").text),
                    "engine": allEngines,
                    "phase": ["to", "ld"],
                    "flight": ["civ"],
                }
            )
            GPFparamList.append(
                {
                    "name": "ang_bank_nom",
                    "value": float(AngBank.find("Nom").find("Civ").find("Others").text),
                    "engine": allEngines,
                    "phase": ["ic", "cl", "cr", "des", "hold", "app"],
                    "flight": ["civ"],
                }
            )
            GPFparamList.append(
                {
                    "name": "ang_bank_nom",
                    "value": float(AngBank.find("Nom").find("Mil").text),
                    "engine": allEngines,
                    "phase": allPhases,
                    "flight": ["mil"],
                }
            )
            GPFparamList.append(
                {
                    "name": "ang_bank_max",
                    "value": float(AngBank.find("Max").find("Civ").find("ToLd").text),
                    "engine": allEngines,
                    "phase": ["to", "ld"],
                    "flight": ["civ"],
                }
            )
            GPFparamList.append(
                {
                    "name": "ang_bank_max",
                    "value": float(AngBank.find("Max").find("Civ").find("Hold").text),
                    "engine": allEngines,
                    "phase": ["hold"],
                    "flight": ["civ"],
                }
            )
            GPFparamList.append(
                {
                    "name": "ang_bank_max",
                    "value": float(AngBank.find("Max").find("Civ").find("Others").text),
                    "engine": allEngines,
                    "phase": ["ic", "cl", "cr", "des", "app"],
                    "flight": ["civ"],
                }
            )
            GPFparamList.append(
                {
                    "name": "ang_bank_max",
                    "value": float(AngBank.find("Max").find("Mil").text),
                    "engine": allEngines,
                    "phase": allPhases,
                    "flight": ["mil"],
                }
            )

            GPFparamList.append(
                {
                    "name": "C_des_exp",
                    "value": float(root.find("CDesExp").text),
                    "engine": allEngines,
                    "phase": allPhases,
                    "flight": allFlights,
                }
            )
            GPFparamList.append(
                {
                    "name": "C_th_to",
                    "value": float(root.find("CThTO").text),
                    "engine": allEngines,
                    "phase": allPhases,
                    "flight": allFlights,
                }
            )
            GPFparamList.append(
                {
                    "name": "C_th_cr",
                    "value": float(root.find("CTcr").text),
                    "engine": allEngines,
                    "phase": allPhases,
                    "flight": allFlights,
                }
            )
            GPFparamList.append(
                {
                    "name": "C_v_min_to",
                    "value": float(root.find("CVminTO").text),
                    "engine": allEngines,
                    "phase": allPhases,
                    "flight": allFlights,
                }
            )
            GPFparamList.append(
                {
                    "name": "C_v_min",
                    "value": float(root.find("CVmin").text),
                    "engine": allEngines,
                    "phase": allPhases,
                    "flight": allFlights,
                }
            )

            HmaxList = {}
            for phase in root.find("HmaxList").findall("HmaxPhase"):

                HmaxList[phase.find("Phase").text] = float(phase.find("Hmax").text)

                if phase.find("Phase").text == "TO":
                    GPFparamList.append(
                        {
                            "name": "H_max_to",
                            "value": float(phase.find("Hmax").text),
                            "engine": allEngines,
                            "phase": ["to"],
                            "flight": allFlights,
                        }
                    )

                elif phase.find("Phase").text == "IC":
                    GPFparamList.append(
                        {
                            "name": "H_max_ic",
                            "value": float(phase.find("Hmax").text),
                            "engine": allEngines,
                            "phase": ["ic"],
                            "flight": allFlights,
                        }
                    )

                elif phase.find("Phase").text == "AP":
                    GPFparamList.append(
                        {
                            "name": "H_max_app",
                            "value": float(phase.find("Hmax").text),
                            "engine": allEngines,
                            "phase": ["app"],
                            "flight": allFlights,
                        }
                    )

                elif phase.find("Phase").text == "LD":
                    GPFparamList.append(
                        {
                            "name": "H_max_ld",
                            "value": float(phase.find("Hmax").text),
                            "engine": allEngines,
                            "phase": ["lnd"],
                            "flight": allFlights,
                        }
                    )

            VdList = {}
            for vdphase in root.find("VdList").findall("VdPhase"):

                Phase = vdphase.find("Phase")
                name = Phase.find("name").text
                index = int(Phase.find("index").text)
                Vd = float(vdphase.find("Vd").text)

                if name not in VdList:
                    VdList[name] = {}

                VdList[name][index] = Vd

                if name == "CL":
                    if index == 1:
                        V_cl_1 = Vd
                        GPFparamList.append(
                            {
                                "name": "V_cl_1",
                                "value": V_cl_1,
                                "engine": allEngines,
                                "phase": ["cl"],
                                "flight": allFlights,
                            }
                        )
                    elif index == 2:
                        V_cl_2 = Vd
                        GPFparamList.append(
                            {
                                "name": "V_cl_2",
                                "value": V_cl_2,
                                "engine": allEngines,
                                "phase": ["cl"],
                                "flight": allFlights,
                            }
                        )
                    elif index == 3:
                        V_cl_3 = Vd
                        GPFparamList.append(
                            {
                                "name": "V_cl_3",
                                "value": V_cl_3,
                                "engine": allEngines,
                                "phase": ["cl"],
                                "flight": allFlights,
                            }
                        )
                    elif index == 4:
                        V_cl_4 = Vd
                        GPFparamList.append(
                            {
                                "name": "V_cl_4",
                                "value": V_cl_4,
                                "engine": allEngines,
                                "phase": ["cl"],
                                "flight": allFlights,
                            }
                        )
                    elif index == 5:
                        V_cl_5 = Vd
                        GPFparamList.append(
                            {
                                "name": "V_cl_5",
                                "value": V_cl_5,
                                "engine": allEngines,
                                "phase": ["cl"],
                                "flight": allFlights,
                            }
                        )
                    elif index == 6:
                        V_cl_6 = Vd
                        GPFparamList.append(
                            {
                                "name": "V_cl_6",
                                "value": V_cl_6,
                                "engine": allEngines,
                                "phase": ["cl"],
                                "flight": allFlights,
                            }
                        )
                    elif index == 7:
                        V_cl_7 = Vd
                        GPFparamList.append(
                            {
                                "name": "V_cl_7",
                                "value": V_cl_7,
                                "engine": allEngines,
                                "phase": ["cl"],
                                "flight": allFlights,
                            }
                        )
                    elif index == 8:
                        V_cl_8 = Vd
                        GPFparamList.append(
                            {
                                "name": "V_cl_8",
                                "value": V_cl_8,
                                "engine": allEngines,
                                "phase": ["cl"],
                                "flight": allFlights,
                            }
                        )

                if name == "DES":
                    if index == 1:
                        V_des_1 = Vd
                        GPFparamList.append(
                            {
                                "name": "V_des_1",
                                "value": V_des_1,
                                "engine": allEngines,
                                "phase": ["des"],
                                "flight": allFlights,
                            }
                        )
                    elif index == 2:
                        V_des_2 = Vd
                        GPFparamList.append(
                            {
                                "name": "V_des_2",
                                "value": V_des_2,
                                "engine": allEngines,
                                "phase": ["des"],
                                "flight": allFlights,
                            }
                        )
                    elif index == 3:
                        V_des_3 = Vd
                        GPFparamList.append(
                            {
                                "name": "V_des_3",
                                "value": V_des_3,
                                "engine": allEngines,
                                "phase": ["des"],
                                "flight": allFlights,
                            }
                        )
                    elif index == 4:
                        V_des_4 = Vd
                        GPFparamList.append(
                            {
                                "name": "V_des_4",
                                "value": V_des_4,
                                "engine": allEngines,
                                "phase": ["des"],
                                "flight": allFlights,
                            }
                        )
                    elif index == 5:
                        V_des_5 = Vd
                        GPFparamList.append(
                            {
                                "name": "V_des_5",
                                "value": V_des_5,
                                "engine": allEngines,
                                "phase": ["des"],
                                "flight": allFlights,
                            }
                        )
                    elif index == 6:
                        V_des_6 = Vd
                        GPFparamList.append(
                            {
                                "name": "V_des_6",
                                "value": V_des_6,
                                "engine": allEngines,
                                "phase": ["des"],
                                "flight": allFlights,
                            }
                        )
                    elif index == 7:
                        V_des_7 = Vd
                        GPFparamList.append(
                            {
                                "name": "V_des_7",
                                "value": V_des_7,
                                "engine": allEngines,
                                "phase": ["des"],
                                "flight": allFlights,
                            }
                        )

            VList = {}
            for vphase in root.find("VList").findall("VPhase"):

                Phase = vphase.find("Phase")
                name = Phase.find("name").text
                index = int(Phase.find("index").text)
                V = float(vphase.find("V").text)

                if name not in VList:
                    VList[name] = {}

                VList[name][index] = V

                if name == "HOLD":
                    if index == 1:
                        V_hold_1 = V
                        GPFparamList.append(
                            {
                                "name": "V_hold_1",
                                "value": V_hold_1,
                                "engine": allEngines,
                                "phase": ["hold"],
                                "flight": allFlights,
                            }
                        )
                    elif index == 2:
                        V_hold_2 = V
                        GPFparamList.append(
                            {
                                "name": "V_hold_2",
                                "value": V_hold_2,
                                "engine": allEngines,
                                "phase": ["hold"],
                                "flight": allFlights,
                            }
                        )
                    elif index == 3:
                        V_hold_3 = V
                        GPFparamList.append(
                            {
                                "name": "V_hold_3",
                                "value": V_hold_3,
                                "engine": allEngines,
                                "phase": ["hold"],
                                "flight": allFlights,
                            }
                        )
                    elif index == 4:
                        V_hold_4 = V
                        GPFparamList.append(
                            {
                                "name": "V_hold_4",
                                "value": V_hold_4,
                                "engine": allEngines,
                                "phase": ["hold"],
                                "flight": allFlights,
                            }
                        )

            V_backtrack = float(root.find("Vground").find("backtrack").text)
            GPFparamList.append(
                {
                    "name": "V_backtrack",
                    "value": V_backtrack,
                    "engine": allEngines,
                    "phase": ["gnd"],
                    "flight": allFlights,
                }
            )
            V_taxi = float(root.find("Vground").find("taxi").text)
            GPFparamList.append(
                {
                    "name": "V_taxi",
                    "value": V_taxi,
                    "engine": allEngines,
                    "phase": ["gnd"],
                    "flight": allFlights,
                }
            )
            V_apron = float(root.find("Vground").find("apron").text)
            GPFparamList.append(
                {
                    "name": "V_apron",
                    "value": V_apron,
                    "engine": allEngines,
                    "phase": ["gnd"],
                    "flight": allFlights,
                }
            )
            V_gate = float(root.find("Vground").find("gate").text)
            GPFparamList.append(
                {
                    "name": "V_gate",
                    "value": V_gate,
                    "engine": allEngines,
                    "phase": ["gnd"],
                    "flight": allFlights,
                }
            )

            CredList = {}
            for CredEng in root.find("CredList").findall("CredEng"):
                EngineType = CredEng.find("EngineType").text
                Cred = float(CredEng.find("Cred").text)

                CredList[EngineType] = Cred

                if EngineType == "JET":
                    GPFparamList.append(
                        {
                            "name": "C_red_jet",
                            "value": Cred,
                            "engine": ["jet"],
                            "phase": ["ic", "cl"],
                            "flight": allFlights,
                        }
                    )
                elif EngineType == "TURBOPROP":
                    GPFparamList.append(
                        {
                            "name": "C_red_turbo",
                            "value": Cred,
                            "engine": ["tbp"],
                            "phase": ["ic", "cl"],
                            "flight": allFlights,
                        }
                    )
                elif EngineType == "PISTON":
                    GPFparamList.append(
                        {
                            "name": "C_red_piston",
                            "value": Cred,
                            "engine": ["pst"],
                            "phase": ["ic", "cl"],
                            "flight": allFlights,
                        }
                    )
                elif EngineType == "ELECTRIC":
                    GPFparamList.append(
                        {
                            "name": "C_red_elec",
                            "value": Cred,
                            "engine": ["elc"],
                            "phase": ["ic", "cl"],
                            "flight": allFlights,
                        }
                    )

        return GPFparamList

    @staticmethod
    def parseGPF(filePath, badaVersion):
        GPFdata = Parser.readGPF(filePath, badaVersion)

        # if ASCI GPF does not exist, try XML GPF file
        if not GPFdata:
            GPFdata = Parser.readGPFXML(filePath, badaVersion)

        # Single row dataframe
        data = {"GPFdata": [GPFdata]}
        df_single = pd.DataFrame(data)

        return df_single

    @staticmethod
    def getGPFValue(GPFdata, name, engine="JET", phase="cr", flight="civ"):
        """This function returns value of the GPF parameter based on the defined features
        like flight, Engine and Phase of flight

        :param Name: name of the GPF parameter.
        :param Engine: type of the engine where this parameter can be applied.
        :param Phase: phase of the flight, where this parameter can be applied.
        :param flight: flight where this parameter can be applied (civ or mil).
        :type Name: str.
        :type Engine: str.
        :type Phase: str.
        :type flight: str.

        """
        # implementation required because 3.16 GPF contains different engine names than 3.15 GPF file
        if engine == "JET":
            engineList = [engine, "jet"]
        if engine == "TURBOPROP":
            engineList = [engine, "turbo", "tbp"]
        if engine == "PISTON":
            engineList = [engine, "piston", "pst"]
        if engine == "ELECTRIC":
            engineList = [engine, "electric", "elc"]

        for param in GPFdata:
            if (
                (param["name"] == name)
                & (any(i in engineList for i in param["engine"]))
                & (phase in param["phase"])
                & (flight in param["flight"])
            ):
                return float(param["value"])
        return None

    @staticmethod
    def combineACDATA_GPF(ACDataFrame, GPFDataframe):
        """This function combines 2 dataframes, the parsed aircraft file
        and parsed GPF file
        """

        # Combine data with GPF data (temporary solution)
        combined_df = pd.concat(
            [ACDataFrame.reset_index(drop=True), GPFDataframe.reset_index(drop=True)],
            axis=1,
        )

        return combined_df

    @staticmethod
    def parseAll(badaVersion, filePath=None):
        """This function parses all BADA3 formatted file and stores
        all data in the final dataframe containing all the BADA data.

        :param filePath: path to the BADA3 formatted file.
        :type filePath: str.
        :raises: IOError
        """

        if filePath == None:
            filePath = configuration.getAircraftPath()
        else:
            filePath = filePath

        # parsing GPF file
        GPFparsedDataframe = Parser.parseGPF(filePath, badaVersion)

        # try to get subfolders, if they exist
        # get names of all the folders in the main BADA model folder to search for XML files
        folderPath = os.path.join(filePath, "BADA3", badaVersion)
        subfolders = Parser.list_subfolders(folderPath)

        if not subfolders:
            # use APF and OPF files
            merged_df = pd.DataFrame()

            # get synonym-filename pairs
            synonym_fileName = Parser.readSynonym(filePath, badaVersion)

            for synonym in synonym_fileName:
                file = synonym_fileName[synonym]

                # parse the original data of a model
                OPFDataFrame = Parser.parseOPF(filePath, badaVersion, file)
                APFDataFrame = Parser.parseAPF(filePath, badaVersion, file)

                df = Parser.combineOPF_APF(OPFDataFrame, APFDataFrame)

                # rename acName in the dateaframe to match the synonym model name
                df.at[0, "acName"] = synonym

                # Combine data with GPF data (temporary solution)
                combined_df = Parser.combineACDATA_GPF(df, GPFparsedDataframe)

                # Merge DataFrames
                merged_df = pd.concat([merged_df, combined_df], ignore_index=True)

            return merged_df

        else:
            # use xml files inside those subfolders
            merged_df = pd.DataFrame()

            # get synonym-filename pairs
            synonym_fileName = Parser.readSynonymXML(filePath, badaVersion)

            for synonym in synonym_fileName:
                file = synonym_fileName[synonym]

                if file in subfolders:
                    # parse the original XML of a model
                    df = Parser.parseXML(filePath, badaVersion, file)

                    # rename acName in the dateaframe to match the synonym model name
                    df.at[0, "acName"] = synonym

                    # Combine data with GPF data (temporary solution)
                    combined_df = Parser.combineACDATA_GPF(df, GPFparsedDataframe)

                    # Merge DataFrames
                    merged_df = pd.concat([merged_df, combined_df], ignore_index=True)

            return merged_df

    @staticmethod
    def getBADAParameters(df, acName, parameters):
        """Retrieves specified parameters for a given aircraft name from the DataFrame.

        :param df: The DataFrame containing aircraft data.
        :param acName: The name of the aircraft to search for
        :param parameters: A list of column names to retrieve or a single column name
        :type df: pandas dataframe.
        :type acName: list[string].
        :type parameters: list[string].
        :return: parameter values: dataframe
        :rtype: dataframe.
        """

        # Ensure parameters is a list
        if isinstance(parameters, str):
            parameters = [parameters]

        # Ensure acName is a list
        if isinstance(acName, str):
            acName = [acName]

        # Ensure all requested parameters exist in the DataFrame
        missing_cols = [col for col in parameters if col not in df.columns]
        if missing_cols:
            raise ValueError(
                f"The following parameters are not in the DataFrame columns: {missing_cols}"
            )

        # Filter rows where 'acName' matches any of the specified aircraft names
        filtered_df = df[df["acName"].isin(acName)]

        # Check if any rows were found
        if filtered_df.empty:
            raise ValueError(f"No entries found for aircraft(s): {acName}.")
        else:
            # Select the required columns
            result_df = filtered_df[["acName"] + parameters].reset_index(drop=True)
            return result_df

    @staticmethod
    def safe_get(df, column_name, default_value=None):
        """Accessing a potentially dropped column from a dataframe"""

        if column_name in df.columns:
            return df[column_name].iloc[0]
        else:
            return default_value


class BADA3(Airplane):
    """This class implements the part of BADA3 performance model that will be used in other classes following the BADA3 manual.

    :param AC: parsed aircraft.
    :type AC: bada3.Parse.
    """

    def __init__(self, AC):
        super().__init__()
        self.AC = AC

    def CL(self, sigma, mass, tas, nz=1.0):
        """This function computes the lift coefficient

        :param tas: true airspeed [m s^-1].
        :param sigma: normalised air density [-].
        :param mass: aircraft mass [kg].
        :param nz: load factor [-].
        :type tas: float.
        :type sigma: float.
        :type mass: float.
        :type nz: float.
        :returns: Lift coefficient [-].
        :rtype: float.
        """

        return 2 * mass * const.g * nz / (sigma * const.rho_0 * tas * tas * self.AC.S)

    def CD(
        self, CL, config, expedite=False, speedBrakes={"deployed": False, "value": 0.03}
    ):
        """This function computes the drag coefficient

        :param CL: Lift coefficient [-].
        :param config: aircraft aerodynamic configuration [CR/IC/TO/AP/LD][-].
        :param expedite: expedite descend factor [-].
        :param speedBrakes: speed brakes used or not [-].
        :type CL: float.
        :type config: str.
        :type expedite: boolean.
        :type speedBrakes: boolean.
        :returns: Drag coefficient [-].
        :rtype: float
        :raises: ValueError
        """

        if self.AC.xmlFiles:
            HLid_CR = self.AC.aeroConfig["CR"]["HLid"]
            LG_CR = self.AC.aeroConfig["CR"]["LG"]
            HLid_AP = self.AC.aeroConfig["AP"]["HLid"]
            LG_AP = self.AC.aeroConfig["AP"]["LG"]
            HLid_LD = self.AC.aeroConfig["LD"]["HLid"]
            LG_LD = self.AC.aeroConfig["LD"]["LG"]

            if (
                self.AC.CD0[HLid_AP][LG_AP] == 0.0
                and self.AC.CD0[HLid_LD][LG_LD] == 0.0
                and self.AC.CD2[HLid_AP][LG_AP] == 0.0
                and self.AC.CD2[HLid_LD][LG_LD] == 0.0
                and self.AC.DeltaCD == 0.0
            ):

                CD = self.AC.CD0[HLid_CR][LG_CR] + self.AC.CD2[HLid_CR][LG_CR] * (
                    CL * CL
                )
            else:
                if config == "CR" or config == "IC" or config == "TO":
                    CD = (
                        self.AC.CD0[HLid_CR][LG_CR]
                        + self.AC.CD2[HLid_CR][LG_CR] * CL**2
                    )
                elif config == "AP":
                    CD = (
                        self.AC.CD0[HLid_AP][LG_AP]
                        + self.AC.CD2[HLid_AP][LG_AP] * CL**2
                    )
                elif config == "LD":
                    CD = (
                        self.AC.CD0[HLid_LD][LG_LD]
                        + self.AC.DeltaCD
                        + self.AC.CD2[HLid_LD][LG_LD] * CL**2
                    )
                else:
                    return float("Nan")
        else:
            if (
                self.AC.CD0["AP"] == 0.0
                and self.AC.CD0["LD"] == 0.0
                and self.AC.CD2["AP"] == 0.0
                and self.AC.CD2["LD"] == 0.0
                and self.AC.CD0["GEAR_DOWN"] == 0.0
            ):
                CD = self.AC.CD0["CR"] + self.AC.CD2["CR"] * (CL * CL)

            else:
                if config == "CR" or config == "IC" or config == "TO":
                    CD = self.AC.CD0["CR"] + self.AC.CD2["CR"] * CL**2
                elif config == "AP":
                    CD = self.AC.CD0[config] + self.AC.CD2[config] * CL**2
                elif config == "LD":
                    CD = (
                        self.AC.CD0[config]
                        + self.AC.CD0["GEAR_DOWN"]
                        + self.AC.CD2[config] * CL**2
                    )
                else:
                    return float("Nan")

        # implementation of a simple speed brakes model
        if speedBrakes["deployed"]:
            if speedBrakes["value"] is not None:
                CD = CD + speedBrakes["value"]
            else:
                CD = CD + 0.03
            return CD

        # expedite descent
        C_des_exp = 1.0
        if expedite:
            C_des_exp = Parser.getGPFValue(self.AC.GPFdata, "C_des_exp", phase="des")
            CD = CD * C_des_exp

        return CD

    def D(self, sigma, tas, CD):
        """This function computes the aerodynamic drag

        :param tas: true airspeed [m s^-1].
        :param sigma: normalised air density [-].
        :param CD: Drag coefficient [-].
        :type tas: float.
        :type sigma: float.
        :type CD: float.
        :returns: Aerodynamic drag [N].
        :rtype: float.
        """

        return 0.5 * sigma * const.rho_0 * tas * tas * self.AC.S * CD

    def L(self, sigma, tas, CL):
        """This function computes the aerodynamic lift

        :param tas: true airspeed [m s^-1].
        :param sigma: normalised air density [-].
        :param CL: Lift coefficient [-].
        :type tas: float.
        :type sigma: float.
        :type CL: float.
        :returns: Aerodynamic lift [N].
        :rtype: float.
        """

        return 0.5 * sigma * const.rho_0 * tas * tas * self.AC.S * CL

    def Thrust(self, h, DeltaTemp, rating, v, config, **kwargs):
        """This function computes the aircraft thrust

        :param rating: engine rating {MCMB,MCRZ,MTKF,LIDL}.
        :param h: altitude [m].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param v: true airspeed (TAS) [m s^-1].
        :param config: aircraft aerodynamic configuration [CR/IC/TO/AP/LD][-].
        :type rating: string.
        :type DeltaTemp: float
        :type h: float.
        :type v: float.
        :type config: string.
        :returns: Thrust [N].
        :rtype: float
        """

        if rating == "MCMB":
            # MCMB
            T = self.TMax(h=h, DeltaTemp=DeltaTemp, rating=rating, v=v)

        elif rating == "MTKF":
            # MTKF
            T = self.TMax(h=h, DeltaTemp=DeltaTemp, rating=rating, v=v)

        elif rating == "MCRZ":
            # MCRZ
            T = self.TMax(h=h, DeltaTemp=DeltaTemp, rating=rating, v=v)

        elif rating == "LIDL":
            # LIDL
            T = self.TDes(h=h, DeltaTemp=DeltaTemp, v=v, config=config)

        elif rating == "ADAPTED":
            # ADAPTED
            ROCD = checkArgument("ROCD", **kwargs)
            mass = checkArgument("mass", **kwargs)
            v = checkArgument("v", **kwargs)
            acc = checkArgument("acc", **kwargs)
            Drag = checkArgument("Drag", **kwargs)
            T = self.TAdapted(
                h=h, DeltaTemp=DeltaTemp, ROCD=ROCD, mass=mass, v=v, acc=acc, Drag=Drag
            )

        else:
            T = float("Nan")

        return T

    def TAdapted(self, h, DeltaTemp, ROCD, mass, v, acc, Drag):
        """This function computes the adapted thrust

        :param h: altitude [m].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param ROCD: rate of climb/descend [m s^-1].
        :param mass: actual aircraft weight [kg]
        :param v: true airspeed (TAS) [m s^-1].
        :param acc: acceleration [m s^-2].
        :param Drag: aerodynamic drag [N].
        :type h: float.
        :type DeltaTemp: float.
        :type ROCD: float.
        :type mass: float.
        :type acc: float.
        :type Drag: float.
        :returns: maximum thrust [N].
        :rtype: float
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        tau_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)
        Tadapted = ROCD * mass * const.g * tau_const / v + mass * acc + Drag

        return Tadapted

    def TMax(self, h, DeltaTemp, rating, v):
        """This function computes the maximum thrust

        :param rating: engine rating {MCMB,MCRZ,MTKF}.
        :param h: altitude [m].
        :param v: true airspeed (TAS) [m s^-1].
        :param DeltaTemp: deviation with respect to ISA [K]
        :type rating: string.
        :type h: float.
        :type v: float.
        :type DeltaTemp: float
        :returns: maximum thrust [N].
        :rtype: float
        """

        acModel = self.AC.engineType

        if acModel == "JET":
            TMaxISA = self.AC.Ct[0] * (
                1
                - (conv.m2ft(h)) / self.AC.Ct[1]
                + self.AC.Ct[2] * (conv.m2ft(h)) * (conv.m2ft(h))
            )

        elif acModel == "TURBOPROP":
            TMaxISA = (self.AC.Ct[0] / conv.ms2kt(v)) * (
                1 - conv.m2ft(h) / self.AC.Ct[1]
            ) + self.AC.Ct[2]

        elif acModel == "PISTON" or acModel == "ELECTRIC":
            TMaxISA = self.AC.Ct[0] * (1 - conv.m2ft(h) / self.AC.Ct[1]) + (
                self.AC.Ct[2] / conv.ms2kt(v)
            )

        else:
            return float("Nan")

        DeltaTempEff = DeltaTemp - self.AC.Ct[3]

        if self.AC.Ct[4] < 0:
            Ctc5 = 0
        else:
            Ctc5 = self.AC.Ct[4]

        DeltaTemp_ = Ctc5 * DeltaTempEff

        if DeltaTemp_ <= 0:
            DeltaTemp_ = 0
        elif DeltaTemp_ > 0.4:
            DeltaTemp_ = 0.4

        TMax = TMaxISA * (1 - DeltaTemp_)

        if rating == "MCMB" or rating == "MTKF":
            return TMax

        elif rating == "MCRZ":
            return TMax * Parser.getGPFValue(self.AC.GPFdata, "C_th_cr", phase="cr")

    def TDes(self, h, DeltaTemp, v, config):
        """This function computes the descent thrust

        :param h: altitude [m].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param config: aircraft aerodynamic configuration [CR/IC/TO/AP/LD][-].
        :param v: true airspeed (TAS) [m s^-1].
        :type h: float.
        :type DeltaTemp: float.
        :type config: string.
        :type v: float.
        :returns: minimum thrust [N].
        :rtype: float
        """

        H_max_app = Parser.getGPFValue(self.AC.GPFdata, "H_max_app", phase="app")

        if (
            self.AC.engineType == "PISTON"
            or self.AC.engineType == "ELECTRIC"
            or self.AC.engineType == "TURBOPROP"
        ):
            TMaxClimb = self.TMax(rating="MCMB", h=h, DeltaTemp=DeltaTemp, v=v)
        elif self.AC.engineType == "JET":
            TMaxClimb = self.TMax(rating="MCMB", h=h, DeltaTemp=DeltaTemp, v=v)

        # non-clean data available -> Hp,des cannot be below Hmax,AP
        HpDes_ = self.AC.HpDes

        if self.AC.xmlFiles:
            [HLid_CR, LG_CR] = self.flightEnvelope.getAeroConfig(config="CR")
            [HLid_AP, LG_AP] = self.flightEnvelope.getAeroConfig(config="AP")
            [HLid_LD, LG_LD] = self.flightEnvelope.getAeroConfig(config="LD")

            if (
                self.AC.CD0[HLid_AP][LG_AP] != 0.0
                and self.AC.CD0[HLid_LD][LG_LD] != 0.0
                and self.AC.CD2[HLid_AP][LG_AP] != 0.0
                and self.AC.CD2[HLid_LD][LG_LD] != 0.0
                and self.AC.DeltaCD != 0.0
            ):

                if HpDes_ < H_max_app:
                    HpDes_ = H_max_app

        else:
            if (
                self.AC.CD0["AP"] != 0.0
                and self.AC.CD0["LD"] != 0.0
                and self.AC.CD2["AP"] != 0.0
                and self.AC.CD2["LD"] != 0.0
                and self.AC.CD0["GEAR_DOWN"] != 0.0
            ):
                if HpDes_ < H_max_app:
                    HpDes_ = H_max_app

        if h > conv.ft2m(HpDes_):
            Tdes = self.AC.CTdeshigh * TMaxClimb
        elif h <= conv.ft2m(HpDes_):
            if config == "CR":
                Tdes = self.AC.CTdeslow * TMaxClimb
            elif config == "AP":
                Tdes = self.AC.CTdesapp * TMaxClimb
            elif config == "LD":
                Tdes = self.AC.CTdesld * TMaxClimb
            else:
                Tdes = float("Nan")

        return Tdes

    def ffnom(self, v, T):
        """This function computes the nominal fuel flow

        :param v: true airspeed (TAS) [m s^-1].
        :param T: Thrust [N].
        :type v: float.
        :type T: float.
        :returns: nominal fuel flow [kg s^-1].
        :rtype: float
        """

        if self.AC.engineType == "JET":
            eta = self.AC.Cf[0] * (1 + conv.ms2kt(v) / self.AC.Cf[1]) / (1000 * 60)
            ffnom = eta * T

        elif self.AC.engineType == "TURBOPROP":
            eta = (
                self.AC.Cf[0]
                * (1 - conv.ms2kt(v) / self.AC.Cf[1])
                * (conv.ms2kt(v) / 1000)
                / (1000 * 60)
            )
            ffnom = eta * T

        elif self.AC.engineType == "PISTON" or self.AC.engineType == "ELECTRIC":
            ffnom = self.AC.Cf[0] / 60

        return ffnom

    def ffMin(self, h):
        """This function computes the minimum fuel flow

        :param h: altitude [m].
        :type h: float.
        :returns: Minimum fuel flow [kg s^-1].
        :rtype: float
        """

        if self.AC.engineType == "JET" or self.AC.engineType == "TURBOPROP":
            ffmin = self.AC.CfDes[0] * (1 - (conv.m2ft(h)) / self.AC.CfDes[1]) / 60
        elif self.AC.engineType == "PISTON" or self.AC.engineType == "ELECTRIC":
            ffmin = self.AC.CfDes[0] / 60  # Cf3 param

        return ffmin

    def ff(self, h, v, T, config=None, flightPhase=None, adapted=False):
        """This function computes the fuel flow based on the flight phase and flight situation

        :param h: altitude [m].
        :param rating: engine rating {MCMB,MCRZ,LIDL}.
        :param v: true airspeed (TAS) [m s^-1].
        :param T: Thrust [N].
        :type h: float.
        :type rating: string.
        :type v: float.
        :type T: float.
        :returns: fuel flow [kg s^-1].
        :rtype: float
        """

        if adapted:
            # adapted thrust
            ffnom = self.ffnom(v=v, T=T)
            ff = max(ffnom, self.ffMin(h=h))
        else:
            if flightPhase == "Climb":
                # climb thrust
                ffnom = self.ffnom(v=v, T=T)
                ff = max(ffnom, self.ffMin(h=h))

            elif flightPhase == "Cruise":
                # cruise thrust
                ffnom = self.ffnom(v=v, T=T)
                ff = ffnom * self.AC.CfCrz

            elif flightPhase == "Descent":
                # descent in IDLE
                if config == "CR":
                    ff = self.ffMin(h=h)
                elif config == "AP" or config == "LD":
                    ffnom = self.ffnom(v=v, T=T)
                    ff = max(ffnom, self.ffMin(h=h))
            else:
                ff = float("Nan")

        return ff

    def reducedPower(self, h, mass, DeltaTemp):
        """This function computes the reduced climb power coefficient

        :param h: altitude [m].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param mass: actual aircraft weight [kg]
        :param hMax: aircraft flight envelope max Altitude [m]
        :type h: float.
        :type DeltaTemp: float.
        :type mass: float.
        :type hMax: float.
        :returns: reduced climb power coefficient [-].
        :rtype: float
        """

        hMax = self.flightEnvelope.maxAltitude(mass=mass, DeltaTemp=DeltaTemp)
        mMax = self.AC.mass["maximum"]
        mMin = self.AC.mass["minimum"]

        ep = 1e-6  # floating point precision
        if (h + ep) < 0.8 * hMax:
            if self.AC.engineType == "JET":
                CRed = Parser.getGPFValue(
                    self.AC.GPFdata, "C_red_jet", engine=self.AC.engineType, phase="cl"
                )
            elif self.AC.engineType == "TURBOPROP":
                CRed = Parser.getGPFValue(
                    self.AC.GPFdata, "C_red_turbo", engine="TURBOPROP", phase="cl"
                )
            elif self.AC.engineType == "PISTON":
                CRed = Parser.getGPFValue(
                    self.AC.GPFdata, "C_red_piston", engine="PISTON", phase="cl"
                )
            elif self.AC.engineType == "ELECTRIC":
                CRed = Parser.getGPFValue(
                    self.AC.GPFdata, "C_red_elec", engine="ELECTRIC", phase="cl"
                )
        else:
            CRed = 0

        CPowRed = 1 - CRed * (mMax - mass) / (mMax - mMin)
        return CPowRed

    def ROCD(self, T, D, v, mass, ESF, h, DeltaTemp, reducedPower=False):
        """This function computes the Rate of Climb or Descent

        :param h: altitude [m].
        :param T: aircraft thrust [N].
        :param D: aircraft drag [N].
        :param v: aircraft true airspeed [TAS] [m s^-1].
        :param mass: actual aircraft mass  [kg].
        :param ESF: energy share factor [-].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param reducedPower: power reduction [-]
        :type h: float.
        :type T: float.
        :type D: float.
        :type v: float.
        :type mass: float.
        :type ESF: float.
        :type DeltaTemp: float.
        :type reducedPower: boolean.
        :returns: rate of climb/descend [m/s].
        :rtype: float
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        temp = theta * const.temp_0

        CPowRed = 1.0
        if reducedPower:
            CPowRed = self.reducedPower(h=h, mass=mass, DeltaTemp=DeltaTemp)

        ROCD = (
            ((temp - DeltaTemp) / temp) * (T - D) * v * ESF * CPowRed / (mass * const.g)
        )

        return ROCD


class FlightEnvelope(BADA3):
    """This class is a BADA3 aircraft subclass and implements the flight envelope caclulations
    following the BADA3 manual.

    :param AC: parsed aircraft.
    :type AC: bada3.Parse.
    """

    def __init__(self, AC):
        super().__init__(AC)

    def maxAltitude(self, mass, DeltaTemp):
        """This function computes the maximum altitude for any given aircraft mass

        :param DeltaTemp: deviation with respect to ISA [K]
        :param mass: actual aircraft weight [kg]
        :type DeltaTemp: float.
        :type mass: float.
        :returns: maximum altitude [m].
        :rtype: float
        """

        Gt = self.AC.tempGrad
        Gw = self.AC.mass["mass grad"]
        Ctc4 = self.AC.Ct[3]
        mMax = self.AC.mass["maximum"]

        if Gw < 0:
            Gw = 0
        if Gt > 0:
            Gt = 0

        var = DeltaTemp - Ctc4
        if var < 0:
            var = 0

        if self.AC.Hmax == 0:
            hMax = self.AC.hmo
        else:
            hMax = min(self.AC.hmo, self.AC.Hmax + Gt * var + Gw * (mMax - mass))

        return conv.ft2m(hMax)

    def VStall(self, mass, config):
        """This function computes the mass correction for stall speed calculation

        :param config: aircraft configuration [CR/IC/TO/AP/LD][-]
        :param mass: aircraft operating mass [kg]
        :type config: string.
        :type mass: string
        :returns: aircraft stall speed [m s^-1]
        :rtype: float
        """

        if self.AC.xmlFiles:
            [HLid, LG] = self.getAeroConfig(config=config)
            vStall = conv.kt2ms(self.AC.Vstall[HLid][LG]) * sqrt(
                mass / self.AC.mass["reference"]
            )
        else:
            vStall = conv.kt2ms(self.AC.Vstall[config]) * sqrt(
                mass / self.AC.mass["reference"]
            )

        return vStall

    def VMin(self, h, mass, config, DeltaTemp, nz=1.2, envelopeType="OPERATIONAL"):
        """This function computes the minimum speed

        :param h: altitude [m].
        :param config: aircraft configuration [CR/IC/TO/AP/LD][-]
        :param mass: aircraft operating mass [kg]
        :param DeltaTemp: deviation with respect to ISA [K]
        :type h: float.
        :type config: string
        :type mass: float
        :type DeltaTemp: float.
        :returns: minimum speed [m s^-1].
        :rtype: float
        """

        if envelopeType == "OPERATIONAL":
            if config == "TO":
                VminStall = Parser.getGPFValue(
                    self.AC.GPFdata, "C_v_min_to", phase="to"
                ) * self.VStall(mass=mass, config=config)
            else:
                VminStall = Parser.getGPFValue(
                    self.AC.GPFdata, "C_v_min"
                ) * self.VStall(mass=mass, config=config)
        elif envelopeType == "CERTIFIED":
            VminStall = self.VStall(mass=mass, config=config)

        if self.AC.Clbo == 0.0 and self.AC.k == 0.0:
            Vmin = VminStall
        else:
            if h < conv.ft2m(15000):
                Vmin = VminStall
            elif h >= conv.ft2m(15000):
                # low speed buffeting limit applies only for JET and TURBOPROP
                if self.AC.engineType == "JET" or self.AC.engineType == "TURBOPROP":
                    [theta, delta, sigma] = atm.atmosphereProperties(
                        h=h, DeltaTemp=DeltaTemp
                    )
                    buffetLimit = self.lowSpeedBuffetLimit(
                        h=h, mass=mass, DeltaTemp=DeltaTemp, nz=nz
                    )
                    if buffetLimit == float("Nan"):
                        Vmin = VminStall
                    else:
                        Vmin = max(
                            VminStall,
                            atm.mach2Cas(
                                buffetLimit, theta=theta, delta=delta, sigma=sigma
                            ),
                        )
                elif self.AC.engineType == "PISTON" or self.AC.engineType == "ELECTRIC":
                    Vmin = VminStall

        return Vmin

    def Vmax_thrustLimited(self, h, mass, DeltaTemp, rating, config):
        """This function computes the maximum CAS speed within the certified flight envelope taking into account the trust limitation.

        :param h: altitude [m].
        :param mass: aircraft operating mass [kg]
        :param DeltaTemp: deviation with respect to ISA [K]
        :param rating: aircraft engine rating [MTKF/MCMB/MCRZ][-]
        :param config: aircraft configuration [TO/IC/CR][-]
        :type h: float.
        :type mass: float
        :type DeltaTemp: float.
        :type config: string
        :type rating: string
        :returns: maximum thrust lmited speed [m s^-1].
        :rtype: float
        """

        [theta, delta, sigma] = atm.atmosphereProperties(h=h, DeltaTemp=DeltaTemp)

        VmaxCertified = self.VMax(h=h, DeltaTemp=DeltaTemp)
        VminCertified = self.VMin(
            h=h,
            mass=mass,
            config=config,
            DeltaTemp=DeltaTemp,
            nz=1.0,
            envelopeType="CERTIFIED",
        )

        maxCASList = []
        DragValue = None
        ThrustValue = None
        CDvalue = None
        CASValue = None
        MValue = None
        for CAS in np.linspace(VminCertified, VmaxCertified, num=200, endpoint=True):
            TAS = atm.cas2Tas(cas=CAS, delta=delta, sigma=sigma)
            M = atm.cas2Mach(cas=CAS, theta=theta, delta=delta, sigma=sigma)
            maxThrust = self.Thrust(
                h=h, DeltaTemp=DeltaTemp, rating=rating, v=TAS, config=config
            )
            CL = self.CL(sigma=sigma, mass=mass, tas=TAS, nz=1.0)
            CD = self.CD(CL=CL, config=config)
            Drag = self.D(sigma=sigma, tas=TAS, CD=CD)

            if maxThrust >= Drag:
                maxCASList.append(CAS)
                DragValue = Drag
                ThrustValue = maxThrust
                CDvalue = CD
                CASValue = conv.ms2kt(CAS)
                MValue = M

        if not maxCASList:
            return None
        else:
            return max(maxCASList)

    def Vx(self, h, mass, DeltaTemp, rating, config):
        """This function computes the best angle of climb CAS speed.

        :param h: altitude [m].
        :param mass: aircraft operating mass [kg]
        :param DeltaTemp: deviation with respect to ISA [K]
        :param rating: aircraft engine rating [MTKF/MCMB/MCRZ][-]
        :param config: aircraft configuration [TO/IC/CR][-]
        :type h: float.
        :type mass: float
        :type DeltaTemp: float.
        :type config: string
        :type rating: string
        :returns: Vx - best angle of climb speed [m s^-1].
        :rtype: float
        """

        [theta, delta, sigma] = atm.atmosphereProperties(h=h, DeltaTemp=DeltaTemp)

        VmaxCertified = self.VMax(h=h, DeltaTemp=DeltaTemp)
        VminCertified = self.VMin(
            h=h,
            mass=mass,
            config=config,
            DeltaTemp=DeltaTemp,
            nz=1.0,
            envelopeType="CERTIFIED",
        )

        excessThrustList = []
        VxList = []

        for CAS in np.linspace(VminCertified, VmaxCertified, num=200, endpoint=True):
            TAS = atm.cas2Tas(cas=CAS, delta=delta, sigma=sigma)
            maxThrust = self.Thrust(
                h=h, DeltaTemp=DeltaTemp, rating=rating, v=TAS, config=config
            )
            CL = self.CL(sigma=sigma, mass=mass, tas=TAS, nz=1.0)
            CD = self.CD(CL=CL, config=config)
            Drag = self.D(sigma=sigma, tas=TAS, CD=CD)

            excessThrustList.append(maxThrust - Drag)
            VxList.append(CAS)

        idx = excessThrustList.index(max(excessThrustList))

        return VxList[idx]

    def VMax(self, h, DeltaTemp):
        """This function computes the maximum speed

        :param h: altitude [m].
        :param DeltaTemp: deviation with respect to ISA [K]
        :type h: float.
        :type DeltaTemp: float.
        :returns: maximum speed [m s^-1].
        :rtype: float
        """

        crossoverAlt = atm.crossOver(cas=conv.kt2ms(self.AC.VMO), Mach=self.AC.MMO)

        if h >= crossoverAlt:
            [theta, delta, sigma] = atm.atmosphereProperties(h=h, DeltaTemp=DeltaTemp)
            VMax = atm.mach2Cas(Mach=self.AC.MMO, theta=theta, delta=delta, sigma=sigma)
        else:
            VMax = conv.kt2ms(self.AC.VMO)

        return VMax

    def lowSpeedBuffetLimit(self, h, mass, DeltaTemp, nz=1.2):
        """This function computes the low speed Buffeting limit using numerical methods by numpy

        :param h: altitude [m].
        :param mass: aircraft mass [kg]
        :param DeltaTemp: deviation with respect to ISA [K]
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :returns: low speed buffet limit as M [-].
        :rtype: float
        """

        p = atm.delta(h, DeltaTemp) * const.p_0

        a1 = self.AC.k
        a2 = -(self.AC.Clbo)
        a3 = (mass * const.g) / (self.AC.S * p * (0.7 / nz))

        coef = [a1, a2, 0, a3]
        roots = np.roots(coef)

        Mb = list()
        for root in roots:
            if root > 0 and not isinstance(root, complex):
                Mb.append(root)
        if not Mb:
            return float("Nan")

        return min(Mb)

    def getConfig(self, phase, h, mass, v, DeltaTemp, hRWY=0.0, nz=1.2):
        """This function returns the aircraft aerodynamic configuration
        based on the aircraft altitude and speed and phase of flight

        :param hRWY: runway elevation AMSL [m].
        :param phase: aircraft phase of flight [cl/cr/des][-].
        :param h: altitude [m].
        :param v: calibrated airspeed (CAS) [m s^-1].
        :param mass: aircraft mass [kg]
        :param DeltaTemp: deviation with respect to ISA [K]
        :type hRWY: float.
        :type phase: string.
        :type h: float.
        :type v: float.
        :type mass: float.
        :type DeltaTemp: float.
        :returns: aircraft aerodynamic configuration [TO/IC/CR/AP/LD][-].
        :rtype: string
        """

        config = None

        # aircraft AGL altitude assuming being close to the RWY [m]
        h_AGL = h - hRWY

        HmaxTO_AGL = (
            conv.ft2m(Parser.getGPFValue(self.AC.GPFdata, "H_max_to", phase="to"))
            - hRWY
        )
        HmaxIC_AGL = (
            conv.ft2m(Parser.getGPFValue(self.AC.GPFdata, "H_max_ic", phase="ic"))
            - hRWY
        )
        HmaxAPP_AGL = (
            conv.ft2m(Parser.getGPFValue(self.AC.GPFdata, "H_max_app", phase="app"))
            - hRWY
        )
        HmaxLD_AGL = (
            conv.ft2m(Parser.getGPFValue(self.AC.GPFdata, "H_max_ld", phase="lnd"))
            - hRWY
        )

        if phase == "Climb" and h_AGL <= HmaxTO_AGL:
            config = "TO"
        elif phase == "Climb" and (h_AGL > HmaxTO_AGL and h_AGL <= HmaxIC_AGL):
            config = "IC"
        else:
            vMinCR = self.VMin(h=h, mass=mass, config="CR", DeltaTemp=DeltaTemp, nz=nz)
            vMinAPP = self.VMin(h=h, mass=mass, config="AP", DeltaTemp=DeltaTemp, nz=nz)
            ep = 1e-6
            if (
                phase == "Descent"
                and (h_AGL + ep) < HmaxLD_AGL
                and (v + ep) < (vMinAPP + conv.kt2ms(10))
            ):
                config = "LD"
            elif (
                phase == "Descent"
                and h_AGL >= HmaxLD_AGL
                and (h_AGL + ep) < HmaxAPP_AGL
                and (v - ep) < (vMinCR + conv.kt2ms(10))
            ) or (
                phase == "Descent"
                and (h_AGL + ep) < HmaxLD_AGL
                and (
                    (v - ep) < (vMinCR + conv.kt2ms(10))
                    and v >= (vMinAPP + conv.kt2ms(10))
                )
            ):
                config = "AP"
            elif (
                (phase == "Climb" and h_AGL > HmaxIC_AGL)
                or phase == "Cruise"
                or (phase == "Descent" and h_AGL >= HmaxAPP_AGL)
                or (
                    phase == "Descent"
                    and (h_AGL + ep) < HmaxAPP_AGL
                    and v >= (vMinCR + conv.kt2ms(10))
                )
            ):
                config = "CR"

        if config is None:
            raise TypeError("Unable to determine aircraft configuration")

        return config

    def getAeroConfig(self, config):
        """This function returns the aircraft aerodynamic configuration
        based on the aerodynamic configuration ID

        :param config: aircraft configuration [CR/IC/TO/AP/LD][-]
        :type config: string
        :returns: aircraft aerodynamic configuration combination of HLID and LG [-].
        :rtype: [float, string]
        """

        HLid = self.AC.aeroConfig[config]["HLid"]
        LG = self.AC.aeroConfig[config]["LG"]

        return [HLid, LG]

    def getHoldSpeed(self, h, theta, delta, sigma, DeltaTemp):
        """This function returns the aircraft holding speed based on the altitude

        :param h: altitude [m].
        :param DeltaTemp: deviation with respect to ISA [K]
        :type h: float.
        :type DeltaTemp: float.
        :returns: aircraft Holding calibrated airspeed (CAS) [m s^-1].
        :rtype: float
        """

        if h <= conv.ft2m(14000):
            vHold = Parser.getGPFValue(self.AC.GPFdata, "V_hold_1", phase="hold")
        elif h > conv.ft2m(14000) and h <= conv.ft2m(20000):
            vHold = Parser.getGPFValue(self.AC.GPFdata, "V_hold_2", phase="hold")
        elif h > conv.ft2m(20000) and h <= conv.ft2m(34000):
            vHold = Parser.getGPFValue(self.AC.GPFdata, "V_hold_3", phase="hold")
        elif h > conv.ft2m(34000):
            MHold = Parser.getGPFValue(self.AC.GPFdata, "V_hold_4", phase="hold")
            vHold = atm.mach2Cas(Mach=M, theta=theta, delta=delta, sigma=sigma)

        return conv.kt2ms(vHold)

    def getGroundMovementSpeed(self, pos):
        """This function returns the aircraft ground movement speed based on postion on the ground

        :param pos: aircraft position on airport ground [backtrack/taxi/apron/gate][-].
        :type pos: string.
        :returns: aircraft ground movement calibrated airspeed (CAS) [m s^-1].
        :rtype: float
        """

        if pos == "backtrack":
            vGround = Parser.getGPFValue(self.AC.GPFdata, "V_backtrack", phase="gnd")
        elif pos == "taxi":
            vGround = Parser.getGPFValue(self.AC.GPFdata, "V_taxi", phase="gnd")
        elif pos == "apron":
            vGround = Parser.getGPFValue(self.AC.GPFdata, "V_apron", phase="gnd")
        elif pos == "gate":
            vGround = Parser.getGPFValue(self.AC.GPFdata, "V_gate", phase="gnd")

        return conv.kt2ms(vGround)

    def getBankAngle(self, phase, flightUnit, value):
        """This function returns the aircraft bank angle based on phase of flight

        :param phase: aircraft phase of flight [to/ic/cl/cr/...][-].
        :param flightUnit: flight unit [civ/mil][-].
        :param value: nominal or maximum value [nom/max][-].
        :type phase: string.
        :type flightUnit: string.
        :type value: string.
        :returns: bank angle [deg]
        :rtype: float
        """

        nomBankAngle = Parser.getGPFValue(
            self.AC.GPFdata, "ang_bank_nom", flightUnit=flightUnit, phase=phase
        )
        maxBankAngle = Parser.getGPFValue(
            self.AC.GPFdata, "ang_bank_max", flightUnit=flightUnit, phase=phase
        )

        if value == "nom":
            return nomBankAngle
        elif value == "max":
            return maxBankAngle

    def isAccOK(self, v1, v2, type="long", flightUnit="civ", deltaTime=1.0):
        """This function checks the limits for longitudinal and normal acceleration

        :param type: logitudinal or normal acceleration [long/norm][-].
        :param flightUnit: flight unit [civ/mil][-].
        :param v1: (long) true airspeed (TAS) at step k-1 [m s^-1].
        :param v1: (norm) vertical airspeed (ROCD) at step k-1 [m s^-1].
        :param v2: (long) true airspeed (TAS) at step k [m s^-1].
        :param v2: (norm) vertical airspeed (ROCD) at step k [m s^-1].
        :param deltaTime: time interval between k and k-1 [s].
        :type type: string.
        :type flightUnit: string.
        :type v1: float.
        :type v2: float.
        :type deltaTime: float.
        :returns: acceleration OK [True] or NOK [False] [-]
        :rtype: boolean
        """

        OK = False

        if flightUnit == "civ":
            if type == "long":
                if (
                    abs(v2 - v1)
                    <= conv.ft2m(Parser.getGPFValue(self.AC.GPFdata, "acc_long_max"))
                    * deltaTime
                ):
                    OK = True

        elif type == "norm":
            if (
                abs(v2 - v1)
                <= conv.ft2m(Parser.getGPFValue(self.AC.GPFdata, "acc_norm_max"))
                * deltaTime
            ):
                OK = True

        # currently undefined for BADA3
        elif flightUnit == "mil":
            OK = True

        return OK

    def getSpeedSchedule(self, phase):
        """This function returns the speed schedule
        based on the phase of flight {Climb, Cruise, Descent}

        :param phase: aircraft phase of flight {Climb, Cruise, Descent}
        :type phase: string
        :returns: speed schedule combination of CAS1, CAS2 and M [m s^-1, m s^-1, -].
        :rtype: [float, float, float]
        """

        if phase == "Climb":
            phase = "cl"
        if phase == "Cruise":
            phase = "cr"
        if phase == "Descent":
            phase = "des"

        CAS1 = self.AC.V1[phase]
        CAS2 = self.AC.V2[phase]
        M = self.AC.M[phase]

        return [CAS1, CAS2, M]

    def checkConfigurationContinuity(self, phase, previousConfig, currentConfig):
        """This function ensures the continuity of the configuration change,
        so, the aerodynamic configuration does not change in the wrong direction based on the phase of the flight

        :param phase: aircraft phase of flight {Climb, Cruise, Descent}
        :param previousConfig: aircraft previous aerodynamic configuration
        :param currentConfig: aircraft current aerodynamic configuration
        :type phase: string
        :type previousConfig: string
        :type currentConfig: string
        :returns: speed new current configuration
        :rtype: string
        """

        newConfig = ""

        # previous configuration is NOT empty/unknown
        if previousConfig is not None:
            if phase == "Descent":
                if currentConfig == "CR" and (
                    previousConfig == "AP" or previousConfig == "LD"
                ):
                    newConfig = previousConfig
                elif currentConfig == "AP" and previousConfig == "LD":
                    newConfig = previousConfig
                else:
                    newConfig = currentConfig

            elif phase == "Climb":
                if currentConfig == "TO" and (
                    previousConfig == "IC" or previousConfig == "CR"
                ):
                    newConfig = previousConfig
                elif currentConfig == "IC" and previousConfig == "CR":
                    newConfig = previousConfig
                else:
                    newConfig = currentConfig

            elif phase == "Cruise":
                newConfig = currentConfig

        # previous configuration is empty/unknown
        else:
            newConfig = currentConfig

        return newConfig


class ARPM:
    """This class is a BADA3 aircraft subclass and implements the Airline Procedure Model (ARPM)
    following the BADA3 user manual.

    :param AC: parsed aircraft.
    :type AC: bada3.Parse.
    """

    def __init__(self, AC):
        self.AC = AC

        self.flightEnvelope = FlightEnvelope(AC)

    def climbSpeed(
        self,
        theta,
        delta,
        mass,
        h,
        DeltaTemp,
        speedSchedule_default=None,
        applyLimits=True,
        config=None,
        procedure="BADA",
        NADP1_ALT=3000,
        NADP2_ALT=[1000, 3000],
    ):
        """This function computes the climb speed schedule CAS speed for any given altitude

        :param h: altitude [m].
        :param mass: aircraft mass [kg].
        :param theta: normalised air temperature [-].
        :param delta: normalised air pressure [-].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param speedSchedule_default: default speed schedule that will overwrite the BADA schedule [Vcl1, Vcl2, Mcl].
        :param applyLimits: apply min/max speed limitation [-].
        :type h: float.
        :type mass: float.
        :type theta: float.
        :type delta: float.
        :type DeltaTemp: float.
        :type speedSchedule_default: [float, float, float].
        :type applyLimits: [boolean].
        :returns: climb calibrated airspeed (CAS) [m s^-1].
        :rtype: float
        """

        phase = "cl"
        acModel = self.AC.engineType
        Cvmin = Parser.getGPFValue(self.AC.GPFdata, "C_v_min", phase=phase)
        CvminTO = Parser.getGPFValue(self.AC.GPFdata, "C_v_min_to", phase="to")
        VstallTO = self.flightEnvelope.VStall(config="TO", mass=mass)
        VstallCR = self.flightEnvelope.VStall(config="CR", mass=mass)

        [Vcl1, Vcl2, Mcl] = self.flightEnvelope.getSpeedSchedule(
            phase=phase
        )  # BADA Climb speed schedule

        if speedSchedule_default is not None:
            Vcl1 = speedSchedule_default[0]
            Vcl2 = speedSchedule_default[1]
            Mcl = speedSchedule_default[2]

        crossOverAlt = atm.crossOver(cas=Vcl2, Mach=Mcl)
        sigma = atm.sigma(theta=theta, delta=delta)

        if procedure == "BADA":
            if acModel == "JET":
                speed = list()
                speed.append(min(Vcl1, conv.kt2ms(250)))
                speed.append(
                    Cvmin * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(self.AC.GPFdata, "V_cl_5", phase=phase)
                    )
                )
                speed.append(
                    Cvmin * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(self.AC.GPFdata, "V_cl_4", phase=phase)
                    )
                )
                speed.append(
                    Cvmin * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(self.AC.GPFdata, "V_cl_3", phase=phase)
                    )
                )
                speed.append(
                    Cvmin * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(self.AC.GPFdata, "V_cl_2", phase=phase)
                    )
                )
                speed.append(
                    Cvmin * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(self.AC.GPFdata, "V_cl_1", phase=phase)
                    )
                )

                n = 1
                while n < len(speed):
                    if speed[n] > speed[n - 1]:
                        speed[n] = speed[n - 1]
                    n = n + 1

                if h < conv.ft2m(1500):
                    cas = speed[5]
                elif h >= conv.ft2m(1500) and h < conv.ft2m(3000):
                    cas = speed[4]
                elif h >= conv.ft2m(3000) and h < conv.ft2m(4000):
                    cas = speed[3]
                elif h >= conv.ft2m(4000) and h < conv.ft2m(5000):
                    cas = speed[2]
                elif h >= conv.ft2m(5000) and h < conv.ft2m(6000):
                    cas = speed[1]
                elif h >= conv.ft2m(6000) and h < conv.ft2m(10000):
                    cas = speed[0]
                elif h >= conv.ft2m(10000) and h < crossOverAlt:
                    cas = Vcl2
                elif h >= crossOverAlt:
                    cas = atm.mach2Cas(Mach=Mcl, theta=theta, delta=delta, sigma=sigma)

            elif acModel == "TURBOPROP" or acModel == "PISTON" or acModel == "ELECTRIC":
                speed = list()
                speed.append(min(Vcl1, conv.kt2ms(250)))
                speed.append(
                    Cvmin * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(
                            self.AC.GPFdata, "V_cl_8", engine="TURBOPROP", phase=phase
                        )
                    )
                )
                speed.append(
                    Cvmin * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(
                            self.AC.GPFdata, "V_cl_7", engine="TURBOPROP", phase=phase
                        )
                    )
                )
                speed.append(
                    Cvmin * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(
                            self.AC.GPFdata, "V_cl_6", engine="TURBOPROP", phase=phase
                        )
                    )
                )

                n = 1
                while n < len(speed):
                    if speed[n] > speed[n - 1]:
                        speed[n] = speed[n - 1]
                    n = n + 1

                if h < conv.ft2m(500):
                    cas = speed[3]
                elif h >= conv.ft2m(500) and h < conv.ft2m(1000):
                    cas = speed[2]
                elif h >= conv.ft2m(1000) and h < conv.ft2m(1500):
                    cas = speed[1]
                elif h >= conv.ft2m(1500) and h < conv.ft2m(10000):
                    cas = speed[0]
                elif h >= conv.ft2m(10000) and h < crossOverAlt:
                    cas = Vcl2
                elif h >= crossOverAlt:
                    cas = atm.mach2Cas(Mach=Mcl, theta=theta, delta=delta, sigma=sigma)

        elif procedure == "NADP1":
            if acModel == "JET":
                speed = list()
                speed.append(min(Vcl1, conv.kt2ms(250)))
                speed.append(
                    CvminTO * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(self.AC.GPFdata, "V_cl_2", phase=phase)
                    )
                )
                n = 1
                while n < len(speed):
                    if speed[n] > speed[n - 1]:
                        speed[n] = speed[n - 1]
                    n = n + 1

                if h < conv.ft2m(NADP1_ALT):
                    cas = speed[1]
                elif h >= conv.ft2m(NADP1_ALT) and h < conv.ft2m(10000):
                    cas = speed[0]
                elif h >= conv.ft2m(10000) and h < crossOverAlt:
                    cas = Vcl2
                elif h >= crossOverAlt:
                    sigma = atm.sigma(theta=theta, delta=delta)
                    cas = atm.mach2Cas(Mach=Mcl, theta=theta, delta=delta, sigma=sigma)

            elif acModel == "TURBOPROP" or acModel == "PISTON":
                speed = list()
                speed.append(min(Vcl1, conv.kt2ms(250)))
                speed.append(
                    CvminTO * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(self.AC.GPFdata, "V_cl_1", phase=phase)
                    )
                )

                n = 1
                while n < len(speed):
                    if speed[n] > speed[n - 1]:
                        speed[n] = speed[n - 1]
                    n = n + 1

                if h < conv.ft2m(NADP1_ALT):
                    cas = speed[1]
                elif h >= conv.ft2m(NADP1_ALT) and h < conv.ft2m(10000):
                    cas = speed[0]
                elif h >= conv.ft2m(10000) and h < crossOverAlt:
                    cas = Vcl2
                elif h >= crossOverAlt:
                    sigma = atm.sigma(theta=theta, delta=delta)
                    cas = atm.mach2Cas(Mach=Mcl, theta=theta, delta=delta, sigma=sigma)

        elif procedure == "NADP2":
            if acModel == "JET":
                speed = list()
                speed.append(min(Vcl1, conv.kt2ms(250)))
                speed.append(
                    Cvmin * VstallCR
                    + conv.kt2ms(
                        Parser.getGPFValue(self.AC.GPFdata, "V_cl_2", phase=phase)
                    )
                )
                speed.append(
                    CvminTO * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(self.AC.GPFdata, "V_cl_2", phase=phase)
                    )
                )

                n = 1
                while n < len(speed):
                    if speed[n] > speed[n - 1]:
                        speed[n] = speed[n - 1]
                    n = n + 1

                if h < conv.ft2m(NADP2_ALT[0]):
                    cas = speed[2]
                elif h >= conv.ft2m(NADP2_ALT[0]) and h < conv.ft2m(NADP2_ALT[1]):
                    cas = speed[1]
                elif h >= conv.ft2m(NADP2_ALT[1]) and h < conv.ft2m(10000):
                    cas = speed[0]
                elif h >= conv.ft2m(10000) and h < crossOverAlt:
                    cas = Vcl2
                elif h >= crossOverAlt:
                    sigma = atm.sigma(theta=theta, delta=delta)
                    cas = atm.mach2Cas(Mach=Mcl, theta=theta, delta=delta, sigma=sigma)

            elif acModel == "TURBOPROP" or acModel == "PISTON":
                speed = list()
                speed.append(min(Vcl1, conv.kt2ms(250)))
                speed.append(
                    Cvmin * VstallCR
                    + conv.kt2ms(
                        Parser.getGPFValue(self.AC.GPFdata, "V_cl_2", phase=phase)
                    )
                )
                speed.append(
                    CvminTO * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(self.AC.GPFdata, "V_cl_1", phase=phase)
                    )
                )

                n = 1
                while n < len(speed):
                    if speed[n] > speed[n - 1]:
                        speed[n] = speed[n - 1]
                    n = n + 1

                if h < conv.ft2m(NADP2_ALT[0]):
                    cas = speed[2]
                elif h >= conv.ft2m(NADP2_ALT[0]) and h < conv.ft2m(NADP2_ALT[1]):
                    cas = speed[1]
                elif h >= conv.ft2m(NADP2_ALT[1]) and h < conv.ft2m(10000):
                    cas = speed[0]
                elif h >= conv.ft2m(10000) and h < crossOverAlt:
                    cas = Vcl2
                elif h >= crossOverAlt:
                    sigma = atm.sigma(theta=theta, delta=delta)
                    cas = atm.mach2Cas(Mach=Mcl, theta=theta, delta=delta, sigma=sigma)

        if applyLimits:
            # check if the speed is within the limits of minimum and maximum speed from the flight envelope, if not, overwrite calculated speed with flight envelope min/max speed
            if config is None:
                config = self.flightEnvelope.getConfig(
                    h=h, phase="Climb", v=cas, mass=mass, DeltaTemp=DeltaTemp
                )
            minSpeed = self.flightEnvelope.VMin(
                h=h, mass=mass, config=config, DeltaTemp=DeltaTemp
            )
            maxSpeed = self.flightEnvelope.VMax(h=h, DeltaTemp=DeltaTemp)

            eps = 1e-6  # float calculation precision
            # empty envelope - keep the original calculated CAS speed
            if maxSpeed < minSpeed:
                if (cas - eps) > maxSpeed and (cas - eps) > minSpeed:
                    return [cas, "V"]
                elif (cas + eps) < minSpeed and (cas + eps) < maxSpeed:
                    return [cas, "v"]
                else:
                    return [cas, "vV"]

            if minSpeed > (cas + eps):
                return [minSpeed, "C"]

            if maxSpeed < (cas - eps):
                return [maxSpeed, "C"]

        return [cas, ""]

    def cruiseSpeed(
        self,
        theta,
        delta,
        mass,
        h,
        DeltaTemp,
        speedSchedule_default=None,
        applyLimits=True,
        config=None,
    ):
        """This function computes the cruise speed schedule CAS speed for any given altitude

        :param h: altitude [m].
        :param mass: aircraft mass [kg].
        :param theta: normalised air temperature [-].
        :param delta: normalised air pressure [-].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param speedSchedule_default: default speed schedule that will overwrite the BADA schedule [Vcr1, Vcr2, Mcr].
        :param applyLimits: apply min/max speed limitation [-].
        :type h: float.
        :type mass: float.
        :type theta: float.
        :type delta: float.
        :type DeltaTemp: float.
        :type speedSchedule_default: [float, float, float].
        :type applyLimits: [boolean].
        :returns: climb calibrated airspeed (CAS) [m s^-1].
        :rtype: float
        """

        phase = "cr"
        acModel = self.AC.engineType

        [Vcr1, Vcr2, Mcr] = self.flightEnvelope.getSpeedSchedule(
            phase=phase
        )  # BADA Cruise speed schedule

        if speedSchedule_default is not None:
            Vcr1 = speedSchedule_default[0]
            Vcr2 = speedSchedule_default[1]
            Mcr = speedSchedule_default[2]

        crossOverAlt = atm.crossOver(cas=Vcr2, Mach=Mcr)
        sigma = atm.sigma(theta=theta, delta=delta)

        if acModel == "JET":
            if h < conv.ft2m(3000):
                cas = min(Vcr1, conv.kt2ms(170))
            elif h >= conv.ft2m(3000) and h < conv.ft2m(6000):
                cas = min(Vcr1, conv.kt2ms(220))
            elif h >= conv.ft2m(6000) and h < conv.ft2m(14000):
                cas = min(Vcr1, conv.kt2ms(250))
            elif h >= conv.ft2m(14000) and h < crossOverAlt:
                cas = Vcr2
            elif h >= crossOverAlt:
                cas = atm.mach2Cas(Mach=Mcr, theta=theta, delta=delta, sigma=sigma)

        elif acModel == "TURBOPROP" or acModel == "PISTON" or acModel == "ELECTRIC":
            if h < conv.ft2m(3000):
                cas = min(Vcr1, conv.kt2ms(150))
            elif h >= conv.ft2m(3000) and h < conv.ft2m(6000):
                cas = min(Vcr1, conv.kt2ms(180))
            elif h >= conv.ft2m(6000) and h < conv.ft2m(10000):
                cas = min(Vcr1, conv.kt2ms(250))
            elif h >= conv.ft2m(10000) and h < crossOverAlt:
                cas = Vcr2
            elif h >= crossOverAlt:
                cas = atm.mach2Cas(Mach=Mcr, theta=theta, delta=delta, sigma=sigma)

        if applyLimits:
            # check if the speed is within the limits of minimum and maximum speed from the flight envelope, if not, overwrite calculated speed with flight envelope min/max speed
            if config is None:
                config = self.flightEnvelope.getConfig(
                    h=h, phase="Cruise", v=cas, mass=mass, DeltaTemp=DeltaTemp
                )

            minSpeed = self.flightEnvelope.VMin(
                h=h, mass=mass, config=config, DeltaTemp=DeltaTemp
            )
            maxSpeed = self.flightEnvelope.VMax(h=h, DeltaTemp=DeltaTemp)

            eps = 1e-6  # float calculation precision
            # empty envelope - keep the original calculated CAS speed
            if maxSpeed < minSpeed:
                if (cas - eps) > maxSpeed and (cas - eps) > minSpeed:
                    return [cas, "V"]
                elif (cas + eps) < minSpeed and (cas + eps) < maxSpeed:
                    return [cas, "v"]
                else:
                    return [cas, "vV"]

            if minSpeed > (cas + eps):
                return [minSpeed, "C"]

            if maxSpeed < (cas - eps):
                return [maxSpeed, "C"]

        return [cas, ""]

    def descentSpeed(
        self,
        theta,
        delta,
        mass,
        h,
        DeltaTemp,
        speedSchedule_default=None,
        applyLimits=True,
        config=None,
    ):
        """This function computes the descent speed schedule CAS speed for any given altitude

        :param h: altitude [m].
        :param mass: aircraft mass [kg].
        :param theta: normalised air temperature [-].
        :param delta: normalised air pressure [-].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param speedSchedule_default: default speed schedule that will overwrite the BADA schedule [Vdes1, Vdes2, Mdes].
        :param applyLimits: apply min/max speed limitation [-].
        :type h: float.
        :type mass: float.
        :type theta: float.
        :type delta: float.
        :type DeltaTemp: float.
        :type speedSchedule_default: [float, float, float].
        :type applyLimits: [boolean].
        :returns: climb calibrated airspeed (CAS) [m s^-1].
        :rtype: float
        """

        phase = "des"
        acModel = self.AC.engineType
        Cvmin = Parser.getGPFValue(self.AC.GPFdata, "C_v_min", phase=phase)
        VstallDES = self.flightEnvelope.VStall(config="LD", mass=mass)

        [Vdes1, Vdes2, Mdes] = self.flightEnvelope.getSpeedSchedule(
            phase=phase
        )  # BADA Descent speed schedule

        if speedSchedule_default is not None:
            Vdes1 = speedSchedule_default[0]
            Vdes2 = speedSchedule_default[1]
            Mdes = speedSchedule_default[2]

        crossOverAlt = atm.crossOver(cas=Vdes2, Mach=Mdes)
        sigma = atm.sigma(theta=theta, delta=delta)

        if acModel == "JET" or acModel == "TURBOPROP":
            speed = list()
            speed.append(min(Vdes1, conv.kt2ms(220)))
            speed.append(
                Cvmin * VstallDES
                + conv.kt2ms(
                    Parser.getGPFValue(self.AC.GPFdata, "V_des_4", phase=phase)
                )
            )
            speed.append(
                Cvmin * VstallDES
                + conv.kt2ms(
                    Parser.getGPFValue(self.AC.GPFdata, "V_des_3", phase=phase)
                )
            )
            speed.append(
                Cvmin * VstallDES
                + conv.kt2ms(
                    Parser.getGPFValue(self.AC.GPFdata, "V_des_2", phase=phase)
                )
            )
            speed.append(
                Cvmin * VstallDES
                + conv.kt2ms(
                    Parser.getGPFValue(self.AC.GPFdata, "V_des_1", phase=phase)
                )
            )

            n = 1
            while n < len(speed):
                if speed[n] > speed[n - 1]:
                    speed[n] = speed[n - 1]
                n = n + 1

            if h < conv.ft2m(1000):
                cas = speed[4]
            elif h >= conv.ft2m(1000) and h < conv.ft2m(1500):
                cas = speed[3]
            elif h >= conv.ft2m(1500) and h < conv.ft2m(2000):
                cas = speed[2]
            elif h >= conv.ft2m(2000) and h < conv.ft2m(3000):
                cas = speed[1]
            elif h >= conv.ft2m(3000) and h < conv.ft2m(6000):
                cas = speed[0]
            elif h >= conv.ft2m(6000) and h < conv.ft2m(10000):
                cas = min(Vdes1, conv.kt2ms(250))
            elif h >= conv.ft2m(10000) and h < crossOverAlt:
                cas = Vdes2
            elif h >= crossOverAlt:
                cas = atm.mach2Cas(Mach=Mdes, theta=theta, delta=delta, sigma=sigma)

        elif acModel == "PISTON" or acModel == "ELECTRIC":
            speed = list()
            speed.append(Vdes1)
            speed.append(
                Cvmin * VstallDES
                + conv.kt2ms(
                    Parser.getGPFValue(
                        self.AC.GPFdata, "V_des_7", engine="PISTON", phase=phase
                    )
                )
            )
            speed.append(
                Cvmin * VstallDES
                + conv.kt2ms(
                    Parser.getGPFValue(
                        self.AC.GPFdata, "V_des_6", engine="PISTON", phase=phase
                    )
                )
            )
            speed.append(
                Cvmin * VstallDES
                + conv.kt2ms(
                    Parser.getGPFValue(
                        self.AC.GPFdata, "V_des_5", engine="PISTON", phase=phase
                    )
                )
            )

            n = 1
            while n < len(speed):
                if speed[n] > speed[n - 1]:
                    speed[n] = speed[n - 1]
                n = n + 1

            if h < conv.ft2m(500):
                cas = speed[3]
            elif h >= conv.ft2m(500) and h < conv.ft2m(1000):
                cas = speed[2]
            elif h >= conv.ft2m(1000) and h < conv.ft2m(1500):
                cas = speed[1]
            elif h >= conv.ft2m(1500) and h < conv.ft2m(10000):
                cas = speed[0]
            elif h >= conv.ft2m(10000) and h < crossOverAlt:
                cas = Vdes2
            elif h >= crossOverAlt:
                cas = atm.mach2Cas(Mach=Mdes, theta=theta, delta=delta, sigma=sigma)

        if applyLimits:
            # check if the speed is within the limits of minimum and maximum speed from the flight envelope, if not, overwrite calculated speed with flight envelope min/max speed
            if config is None:
                config = self.flightEnvelope.getConfig(
                    h=h, phase="Descent", v=cas, mass=mass, DeltaTemp=DeltaTemp
                )
            minSpeed = self.flightEnvelope.VMin(
                h=h, mass=mass, config=config, DeltaTemp=DeltaTemp
            )

            maxSpeed = self.flightEnvelope.VMax(h=h, DeltaTemp=DeltaTemp)

            eps = 1e-6  # float calculation precision
            # empty envelope - keep the original calculated CAS speed
            if maxSpeed < minSpeed:
                if (cas - eps) > maxSpeed and (cas - eps) > minSpeed:
                    return [cas, "V"]
                elif (cas + eps) < minSpeed and (cas + eps) < maxSpeed:
                    return [cas, "v"]
                else:
                    return [cas, "vV"]

            if minSpeed > (cas + eps):
                return [minSpeed, "C"]

            if maxSpeed < (cas - eps):
                return [maxSpeed, "C"]

        return [cas, ""]


class PTD(BADA3):
    """This class implements the PTD file creator for BADA3 aircraft following BADA3 manual.

    :param AC: parsed aircraft.
    :type AC: bada3.Parse.
    """

    def __init__(self, AC):
        super().__init__(AC)

        self.flightEnvelope = FlightEnvelope(AC)
        self.ARPM = ARPM(AC)

    def create(self, DeltaTemp, saveToPath):
        """This function creates the BADA3 PTD file

        :param saveToPath: path to directory where PTF should be stored [-]
        :param DeltaTemp: deviation from ISA temperature [K]
        :type saveToPath: string.
        :type DeltaTemp: float.
        :returns: NONE
        """

        # 3 different mass levels [kg]
        if 1.2 * self.AC.mass["minimum"] > self.AC.mass["reference"]:
            massLow = self.AC.mass["minimum"]
        else:
            massLow = 1.2 * self.AC.mass["minimum"]

        massList = [massLow, self.AC.mass["reference"], self.AC.mass["maximum"]]
        max_alt_ft = self.AC.hmo

        # original PTD altitude list
        altitudeList = list(range(0, 2000, 500))
        altitudeList.extend(range(2000, 4000, 1000))

        if int(max_alt_ft) < 30000:
            altitudeList.extend(range(4000, int(max_alt_ft), 2000))
            altitudeList.append(int(max_alt_ft))
        else:
            altitudeList.extend(range(4000, 30000, 2000))
            altitudeList.extend(range(29000, int(max_alt_ft), 2000))
            altitudeList.append(int(max_alt_ft))

        CLList = []
        for mass in massList:
            CLList.append(
                self.PTD_climb(
                    mass=mass, altitudeList=altitudeList, DeltaTemp=DeltaTemp
                )
            )
        DESList_med = self.PTD_descent(
            mass=self.AC.mass["reference"],
            altitudeList=altitudeList,
            DeltaTemp=DeltaTemp,
        )

        self.save2PTD(
            saveToPath=saveToPath,
            CLList_low=CLList[0],
            CLList_med=CLList[1],
            CLList_high=CLList[2],
            DESList_med=DESList_med,
            DeltaTemp=DeltaTemp,
        )

    def save2PTD(
        self, saveToPath, CLList_low, CLList_med, CLList_high, DESList_med, DeltaTemp
    ):
        """This function saves data to PTD file

        :param saveToPath: path to directory where PTD should be stored [-]
        :param CLList_low: list of PTD data in CLIMB at low aircraft mass [-].
        :param CLList_med: list of PTD data in CLIMB at medium aircraft mass [-].
        :param CLList_high: list of PTD data in CLIMB at high aircraft mass [-].
        :param DESList_med: list of PTD data in DESCENT at medium aircraft mass [-].
        :param DeltaTemp: deviation from ISA temperature [K]
        :type saveToPath: string.
        :type CLList_low: list.
        :type CLList_med: list.
        :type CLList_high: list.
        :type DESList_med: list.
        :type DeltaTemp: float.
        :returns: NONE
        """

        def Nan2Zero(list):
            # replace NAN values by 0 for printing purposes
            for k in range(len(list)):
                for m in range(len(list[k])):
                    if isinstance(list[k][m], float):
                        if isnan(list[k][m]):
                            list[k][m] = 0
            return list

        newpath = saveToPath
        if not os.path.exists(newpath):
            os.makedirs(newpath)

        if DeltaTemp == 0.0:
            ISA = ""
        elif DeltaTemp > 0.0:
            ISA = "+" + str(int(DeltaTemp))
        elif DeltaTemp < 0.0:
            ISA = str(int(DeltaTemp))

        acName = self.AC.acName

        while len(acName) < 6:
            acName = acName + "_"
        filename = saveToPath + acName + "_ISA" + ISA + ".PTD"

        file = open(filename, "w")
        file.write("BADA PERFORMANCE FILE RESULTS\n")
        file = open(filename, "a")
        file.write("=============================\n=============================\n\n")
        file.write("Low mass CLIMBS\n")
        file.write("===============\n\n")
        file.write(
            " FL[-] T[K] p[Pa] rho[kg/m3] a[m/s] TAS[kt] CAS[kt]    M[-] mass[kg] Thrust[N] Drag[N] Fuel[kgm] ESF[-] ROC[fpm] TDC[N]  PWC[-]\n"
        )

        # replace NAN values by 0 for printing purposes
        CLList_low = Nan2Zero(CLList_low)
        CLList_med = Nan2Zero(CLList_med)
        CLList_high = Nan2Zero(CLList_high)
        DESList_med = Nan2Zero(DESList_med)

        for k in range(0, len(CLList_low[0])):
            file.write(
                "%6d %3.0f %6.0f %7.3f %7.0f %8.2f %8.2f %7.2f %6.0f %9.0f %9.0f %7.1f %7.2f %7.0f %8.0f %7.2f \n"
                % (
                    CLList_low[0][k],
                    CLList_low[1][k],
                    CLList_low[2][k],
                    CLList_low[3][k],
                    CLList_low[4][k],
                    CLList_low[5][k],
                    CLList_low[6][k],
                    CLList_low[7][k],
                    CLList_low[8][k],
                    CLList_low[9][k],
                    CLList_low[10][k],
                    CLList_low[11][k],
                    CLList_low[12][k],
                    CLList_low[13][k],
                    CLList_low[14][k],
                    CLList_low[15][k],
                )
            )

        file.write("\n\nMedium mass CLIMBS\n")
        file.write("==================\n\n")
        file.write(
            " FL[-] T[K] p[Pa] rho[kg/m3] a[m/s] TAS[kt] CAS[kt]    M[-] mass[kg] Thrust[N] Drag[N] Fuel[kgm] ESF[-] ROC[fpm] TDC[N]  PWC[-]\n"
        )

        for k in range(0, len(CLList_med[0])):
            file.write(
                "%6d %3.0f %6.0f %7.3f %7.0f %8.2f %8.2f %7.2f %6.0f %9.0f %9.0f %7.1f %7.2f %7.0f %8.0f %7.2f \n"
                % (
                    CLList_med[0][k],
                    CLList_med[1][k],
                    CLList_med[2][k],
                    CLList_med[3][k],
                    CLList_med[4][k],
                    CLList_med[5][k],
                    CLList_med[6][k],
                    CLList_med[7][k],
                    CLList_med[8][k],
                    CLList_med[9][k],
                    CLList_med[10][k],
                    CLList_med[11][k],
                    CLList_med[12][k],
                    CLList_med[13][k],
                    CLList_med[14][k],
                    CLList_med[15][k],
                )
            )

        file.write("\n\nHigh mass CLIMBS\n")
        file.write("================\n\n")
        file.write(
            " FL[-] T[K] p[Pa] rho[kg/m3] a[m/s] TAS[kt] CAS[kt]    M[-] mass[kg] Thrust[N] Drag[N] Fuel[kgm] ESF[-] ROC[fpm] TDC[N]  PWC[-]\n"
        )

        for k in range(0, len(CLList_high[0])):
            file.write(
                "%6d %3.0f %6.0f %7.3f %7.0f %8.2f %8.2f %7.2f %6.0f %9.0f %9.0f %7.1f %7.2f %7.0f %8.0f %7.2f \n"
                % (
                    CLList_high[0][k],
                    CLList_high[1][k],
                    CLList_high[2][k],
                    CLList_high[3][k],
                    CLList_high[4][k],
                    CLList_high[5][k],
                    CLList_high[6][k],
                    CLList_high[7][k],
                    CLList_high[8][k],
                    CLList_high[9][k],
                    CLList_high[10][k],
                    CLList_high[11][k],
                    CLList_high[12][k],
                    CLList_high[13][k],
                    CLList_high[14][k],
                    CLList_high[15][k],
                )
            )

        file.write("\nMedium mass DESCENTS\n")
        file.write("====================\n\n")
        file.write(
            " FL[-] T[K] p[Pa] rho[kg/m3] a[m/s] TAS[kt] CAS[kt]    M[-] mass[kg] Thrust[N] Drag[N] Fuel[kgm] ESF[-] ROD[fpm] TDC[N] gammaTAS[deg]\n"
        )

        for k in range(0, len(DESList_med[0])):
            file.write(
                "%6d %3.0f %6.0f %7.3f %7.0f %8.2f %8.2f %7.2f %6.0f %9.0f %9.0f %7.1f %7.2f %7.0f %8.0f %8.2f \n"
                % (
                    DESList_med[0][k],
                    DESList_med[1][k],
                    DESList_med[2][k],
                    DESList_med[3][k],
                    DESList_med[4][k],
                    DESList_med[5][k],
                    DESList_med[6][k],
                    DESList_med[7][k],
                    DESList_med[8][k],
                    DESList_med[9][k],
                    DESList_med[10][k],
                    DESList_med[11][k],
                    DESList_med[12][k],
                    DESList_med[13][k],
                    DESList_med[14][k],
                    DESList_med[15][k],
                )
            )

        file.write("\nTDC stands for (Thrust - Drag) * Cred\n")

    def PTD_climb(self, mass, altitudeList, DeltaTemp):
        """This function calculates the BADA3 PTD data in CLIMB

        :param mass: aircraft mass [kg]
        :param altitudeList: aircraft altitude list [ft]
        :param DeltaTemp: deviation from ISA temperature [K]
        :type mass: float.
        :type altitudeList: list of int.
        :type DeltaTemp: float.
        :returns: list of PTD CLIMB data [-]
        :rtype: list
        """

        FL_complet = []
        T_complet = []
        p_complet = []
        rho_complet = []
        a_complet = []
        TAS_complet = []
        CAS_complet = []
        M_complet = []
        mass_complet = []
        Thrust_complet = []
        Drag_complet = []
        ff_comlet = []
        ESF_complet = []
        ROCD_complet = []
        TDC_complet = []
        PWC_complet = []

        phase = "cl"

        Vcl1 = self.AC.V1[phase]
        Vcl2 = self.AC.V2[phase]
        Mcl = self.AC.M[phase]

        Vcl1 = min(Vcl1, conv.kt2ms(250))
        crossAlt = atm.crossOver(cas=Vcl2, Mach=Mcl)

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
            [cas, speedUpdated] = self.ARPM.climbSpeed(
                theta=theta,
                delta=delta,
                h=H_m,
                mass=mass,
                DeltaTemp=DeltaTemp,
                speedSchedule_default=[Vcl1, Vcl2, Mcl],
                applyLimits=False,
            )
            tas = atm.cas2Tas(cas=cas, delta=delta, sigma=sigma)
            M = atm.tas2Mach(v=tas, theta=theta)
            a = atm.aSound(theta=theta)
            FL = h / 100

            config = self.flightEnvelope.getConfig(
                h=H_m, phase="Climb", v=cas, mass=mass, DeltaTemp=DeltaTemp
            )

            Thrust = self.Thrust(
                rating="MCMB", v=tas, h=H_m, config=config, DeltaTemp=DeltaTemp
            )
            ff = self.ff(flightPhase="Climb", v=tas, h=H_m, T=Thrust) * 60

            CL = self.CL(tas=tas, sigma=sigma, mass=mass)

            CD = self.CD(CL=CL, config=config)
            Drag = self.D(tas=tas, sigma=sigma, CD=CD)

            CPowRed = self.reducedPower(h=H_m, mass=mass, DeltaTemp=DeltaTemp)
            TDC = (Thrust - Drag) * CPowRed

            if H_m < crossAlt:
                ESF = self.esf(
                    h=H_m, flightEvolution="constCAS", M=M, DeltaTemp=DeltaTemp
                )
            else:
                ESF = self.esf(
                    h=H_m, flightEvolution="constM", M=M, DeltaTemp=DeltaTemp
                )

            ROCD = (
                conv.m2ft(
                    self.ROCD(
                        h=H_m,
                        T=Thrust,
                        D=Drag,
                        v=tas,
                        mass=mass,
                        ESF=ESF,
                        DeltaTemp=DeltaTemp,
                        reducedPower=True,
                    )
                )
                * 60
            )

            FL_complet.append(proper_round(FL))
            T_complet.append(theta * const.temp_0)
            p_complet.append(delta * const.p_0)
            rho_complet.append(sigma * const.rho_0)
            a_complet.append(a)
            TAS_complet.append(conv.ms2kt(tas))
            CAS_complet.append(conv.ms2kt(cas))
            M_complet.append(M)
            mass_complet.append(mass)
            Thrust_complet.append(Thrust)
            Drag_complet.append(Drag)
            ff_comlet.append(ff)
            ESF_complet.append(ESF)
            ROCD_complet.append(ROCD)
            TDC_complet.append(TDC)
            PWC_complet.append(CPowRed)

        CLList = [
            FL_complet,
            T_complet,
            p_complet,
            rho_complet,
            a_complet,
            TAS_complet,
            CAS_complet,
            M_complet,
            mass_complet,
            Thrust_complet,
            Drag_complet,
            ff_comlet,
            ESF_complet,
            ROCD_complet,
            TDC_complet,
            PWC_complet,
        ]

        return CLList

    def PTD_descent(self, mass, altitudeList, DeltaTemp):
        """This function calculates the BADA3 PTD data in DESCENT

        :param mass: aircraft mass [kg]
        :param altitudeList: list of aircraft maximum altitude [ft]
        :param DeltaTemp: deviation from ISA temperature [K]
        :type mass: float.
        :type altitudeList: list of int.
        :type DeltaTemp: float.
        :returns: list of PTD DESCENT data [-]
        :rtype: list
        """

        FL_complet = []
        T_complet = []
        p_complet = []
        rho_complet = []
        a_complet = []
        TAS_complet = []
        CAS_complet = []
        M_complet = []
        mass_complet = []
        Thrust_complet = []
        Drag_complet = []
        ff_comlet = []
        ESF_complet = []
        ROCD_complet = []
        TDC_complet = []
        gamma_complet = []

        phase = "des"

        Vdes1 = self.AC.V1[phase]
        Vdes2 = self.AC.V2[phase]
        Mdes = self.AC.M[phase]

        Vdes1 = min(Vdes1, conv.kt2ms(250))
        crossAlt = atm.crossOver(cas=Vdes2, Mach=Mdes)

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
            [cas, speedUpdated] = self.ARPM.descentSpeed(
                theta=theta,
                delta=delta,
                h=H_m,
                mass=mass,
                DeltaTemp=DeltaTemp,
                speedSchedule_default=[Vdes1, Vdes2, Mdes],
                applyLimits=False,
            )
            tas = atm.cas2Tas(cas=cas, delta=delta, sigma=sigma)
            M = atm.tas2Mach(v=tas, theta=theta)
            a = atm.aSound(theta=theta)
            FL = h / 100

            CL = self.CL(tas=tas, sigma=sigma, mass=mass)
            config = self.flightEnvelope.getConfig(
                h=H_m, phase="Descent", v=cas, mass=mass, DeltaTemp=DeltaTemp
            )
            CD = self.CD(CL=CL, config=config)
            Drag = self.D(tas=tas, sigma=sigma, CD=CD)

            if self.AC.engineType == "PISTON" or self.AC.engineType == "ELECTRIC":
                # PISTON  and ELECTRIC uses LIDL throughout the whole descent phase
                Thrust = self.Thrust(
                    rating="LIDL", v=tas, h=H_m, config="CR", DeltaTemp=DeltaTemp
                )
                ff = (
                    self.ff(
                        flightPhase="Descent",
                        v=tas,
                        h=H_m,
                        T=Thrust,
                        config="CR",
                        adapted=False,
                    )
                    * 60
                )
            else:
                Thrust = self.Thrust(
                    rating="LIDL", v=tas, h=H_m, config=config, DeltaTemp=DeltaTemp
                )
                ff = (
                    self.ff(
                        flightPhase="Descent",
                        v=tas,
                        h=H_m,
                        T=Thrust,
                        config=config,
                        adapted=False,
                    )
                    * 60
                )

            CPowRed = 1.0
            TDC = (Thrust - Drag) * CPowRed

            if H_m < crossAlt:
                ESF = self.esf(
                    h=H_m, flightEvolution="constCAS", M=M, DeltaTemp=DeltaTemp
                )
            else:
                ESF = self.esf(
                    h=H_m, flightEvolution="constM", M=M, DeltaTemp=DeltaTemp
                )

            ROCD = (
                conv.m2ft(
                    self.ROCD(
                        h=H_m,
                        T=Thrust,
                        D=Drag,
                        v=tas,
                        mass=mass,
                        ESF=ESF,
                        DeltaTemp=DeltaTemp,
                    )
                )
                * 60
            )

            tau_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)
            dhdt = (conv.ft2m(ROCD / 60)) * tau_const

            if self.AC.drone:
                gamma = conv.rad2deg(atan(dhdt / tas))
            else:
                gamma = conv.rad2deg(asin(dhdt / tas))

            FL_complet.append(proper_round(FL))
            T_complet.append(theta * const.temp_0)
            p_complet.append(delta * const.p_0)
            rho_complet.append(sigma * const.rho_0)
            a_complet.append(a)
            TAS_complet.append(conv.ms2kt(tas))
            CAS_complet.append(conv.ms2kt(cas))
            M_complet.append(M)
            mass_complet.append(mass)
            Thrust_complet.append(Thrust)
            Drag_complet.append(Drag)
            ff_comlet.append(ff)
            ESF_complet.append(ESF)
            ROCD_complet.append(-1 * ROCD)
            TDC_complet.append(TDC)
            gamma_complet.append(gamma)

        DESList = [
            FL_complet,
            T_complet,
            p_complet,
            rho_complet,
            a_complet,
            TAS_complet,
            CAS_complet,
            M_complet,
            mass_complet,
            Thrust_complet,
            Drag_complet,
            ff_comlet,
            ESF_complet,
            ROCD_complet,
            TDC_complet,
            gamma_complet,
        ]

        return DESList


class PTF(BADA3):
    """This class implements the PTD file creator for BADA3 aircraft following BADA3 manual.

    :param AC: parsed aircraft.
    :type AC: bada3.Parse.
    """

    def __init__(self, AC):
        super().__init__(AC)

        self.flightEnvelope = FlightEnvelope(AC)
        self.ARPM = ARPM(AC)

    def create(self, DeltaTemp, saveToPath):
        """This function creates the BADA3 PTF file

        :param saveToPath: path to directory where PTF should be stored [-]
        :param DeltaTemp: deviation from ISA temperature [K]
        :type saveToPath: string.
        :type DeltaTemp: float.
        :returns: NONE
        """

        # 3 different mass levels [kg]
        if 1.2 * self.AC.mass["minimum"] > self.AC.mass["reference"]:
            massLow = self.AC.mass["minimum"]
        else:
            massLow = 1.2 * self.AC.mass["minimum"]

        massList = [massLow, self.AC.mass["reference"], self.AC.mass["maximum"]]
        max_alt_ft = self.AC.hmo

        # original PTF altitude list
        altitudeList = list(range(0, 2000, 500))
        altitudeList.extend(range(2000, 4000, 1000))

        if int(max_alt_ft) < 30000:
            altitudeList.extend(range(4000, int(max_alt_ft), 2000))
            altitudeList.append(max_alt_ft)
        else:
            altitudeList.extend(range(4000, 30000, 2000))
            altitudeList.extend(range(29000, int(max_alt_ft), 2000))
            altitudeList.append(max_alt_ft)

        CRList = self.PTF_cruise(
            massList=massList, altitudeList=altitudeList, DeltaTemp=DeltaTemp
        )
        CLList = self.PTF_climb(
            massList=massList, altitudeList=altitudeList, DeltaTemp=DeltaTemp
        )
        DESList = self.PTF_descent(
            massList=massList, altitudeList=altitudeList, DeltaTemp=DeltaTemp
        )

        self.save2PTF(
            saveToPath=saveToPath,
            altitudeList=altitudeList,
            massList=massList,
            CRList=CRList,
            CLList=CLList,
            DESList=DESList,
            DeltaTemp=DeltaTemp,
        )

    def save2PTF(
        self, saveToPath, CRList, CLList, DESList, altitudeList, massList, DeltaTemp
    ):
        """This function saves data to PTF file

        :param saveToPath: path to directory where PTF should be stored [-]
        :param altitudeList: aircraft altitude list [ft]
        :param massList: aircraft mass list [kg]
        :param CRList: list of PTF data in CRUISE [-].
        :param CLList: list of PTF data in CLIMB [-].
        :param DESList: list of PTF data in DESCENT [-].
        :param DeltaTemp: deviation from ISA temperature [K]
        :type saveToPath: string.
        :type altitudeList: list of int.
        :type massList: list of int.
        :type CRList: list.
        :type CLList: list.
        :type DESList: list.
        :type DeltaTemp: float.
        :returns: NONE
        """

        def Nan2Zero(list):
            # replace NAN values by 0 for printing purposes
            for k in range(len(list)):
                for m in range(len(list[k])):
                    if isinstance(list[k][m], float):
                        if isnan(list[k][m]):
                            list[k][m] = 0
            return list

        newpath = saveToPath
        if not os.path.exists(newpath):
            os.makedirs(newpath)

        if DeltaTemp == 0.0:
            ISA = ""
        elif DeltaTemp > 0.0:
            ISA = "+" + str(int(DeltaTemp))
        elif DeltaTemp < 0.0:
            ISA = str(int(DeltaTemp))

        acName = self.AC.acName

        while len(acName) < 6:
            acName = acName + "_"
        filename = saveToPath + acName + "_ISA" + ISA + ".PTF"

        V1cl = min(250, conv.ms2kt(self.AC.V1["cl"]))
        V2cl = conv.ms2kt(self.AC.V2["cl"])
        Mcl = self.AC.M["cl"]
        V1des = min(250, conv.ms2kt(self.AC.V1["des"]))
        V2des = conv.ms2kt(self.AC.V2["des"])
        Mdes = self.AC.M["des"]
        V1cr = min(250, conv.ms2kt(self.AC.V1["cr"]))
        V2cr = conv.ms2kt(self.AC.V2["cr"])
        Mcr = self.AC.M["cr"]

        today = date.today()
        d3 = today.strftime("%b %d %Y")
        OPFModDate = self.AC.modificationDateOPF
        APFModDate = self.AC.modificationDateAPF

        file = open(filename, "w")
        file.write(
            "BADA PERFORMANCE FILE                                        %s\n\n" % (d3)
        )
        file = open(filename, "a")
        file.write("AC/Type: %s\n" % (acName))
        file.write(
            "                              Source OPF File:               %s\n"
            % (OPFModDate)
        )
        file.write(
            "                              Source APF file:               %s\n\n"
            % (APFModDate)
        )
        file.write(
            " Speeds:   CAS(LO/HI)  Mach   Mass Levels [kg]         Temperature:  ISA%s\n"
            % (ISA)
        )
        file.write(
            " climb   - %3d/%3d     %4.2f   low     -  %.0f\n"
            % (V1cl, V2cl, Mcl, massList[0])
        )
        file.write(
            " cruise  - %3d/%3d     %4.2f   nominal -  %-6.0f        Max Alt. [ft]:%7d\n"
            % (V1cr, V2cr, Mcr, massList[1], altitudeList[-1])
        )
        file.write(
            " descent - %3d/%3d     %4.2f   high    -  %0.f\n"
            % (V1des, V2des, Mdes, massList[2])
        )
        file.write(
            "==========================================================================================\n"
        )
        file.write(
            " FL |          CRUISE           |               CLIMB               |       DESCENT       \n"
        )
        file.write(
            "    |  TAS          fuel        |  TAS          ROCD         fuel   |  TAS  ROCD    fuel  \n"
        )
        file.write(
            "    | [kts]       [kg/min]      | [kts]        [fpm]       [kg/min] | [kts] [fpm] [kg/min]\n"
        )
        file.write(
            "    |          lo   nom    hi   |         lo    nom    hi    nom    |        nom    nom   \n"
        )
        file.write(
            "==========================================================================================\n"
        )

        # replace NAN values by 0 for printing purposes
        CLList = Nan2Zero(CLList)
        DESList = Nan2Zero(DESList)

        for k in range(0, len(altitudeList)):
            FL = proper_round(altitudeList[k] / 100)
            if FL < 30:
                file.write(
                    "%3.0f |                           |  %3.0f   %5.0f %5.0f %5.0f   %5.1f  |  %3.0f  %5.0f  %5.1f  \n"
                    % (
                        FL,
                        CLList[0][k],
                        CLList[1][k],
                        CLList[2][k],
                        CLList[3][k],
                        CLList[4][k],
                        DESList[0][k],
                        DESList[1][k],
                        DESList[2][k],
                    )
                )
            else:
                file.write(
                    "%3.0f |  %3.0f   %5.1f %5.1f %5.1f  |  %3.0f   %5.0f %5.0f %5.0f   %5.1f  |  %3.0f  %5.0f  %5.1f  \n"
                    % (
                        FL,
                        CRList[0][k],
                        CRList[1][k],
                        CRList[2][k],
                        CRList[3][k],
                        CLList[0][k],
                        CLList[1][k],
                        CLList[2][k],
                        CLList[3][k],
                        CLList[4][k],
                        DESList[0][k],
                        DESList[1][k],
                        DESList[2][k],
                    )
                )
            file.write(
                "    |                           |                                   | \n"
            )

        file.write(
            "==========================================================================================\n"
        )

    def PTF_cruise(self, massList, altitudeList, DeltaTemp):
        """This function calculates the BADA3 PTF data in CRUISE

        :param massList: list of aircraft mass [kg]
        :param altitudeList: aircraft altitude list [ft]
        :param DeltaTemp: deviation from ISA temperature [K]
        :type massList: list.
        :type altitudeList: list of int.
        :type DeltaTemp: float.
        :returns: list of PTF CRUISE data [-]
        :rtype: list
        """

        TAS_CR_complet = []
        FF_CR_LO_complet = []
        FF_CR_NOM_complet = []
        FF_CR_HI_complet = []

        phase = "cr"
        massNominal = massList[1]

        Vcr1 = self.AC.V1[phase]
        Vcr2 = self.AC.V2[phase]
        Mcr = self.AC.M[phase]

        Vcr1 = min(Vcr1, conv.kt2ms(250))

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
            [cas, speedUpdated] = self.ARPM.cruiseSpeed(
                theta=theta,
                delta=delta,
                h=H_m,
                mass=massNominal,
                DeltaTemp=DeltaTemp,
                speedSchedule_default=[Vcr1, Vcr2, Mcr],
                applyLimits=False,
            )
            tas_nominal = atm.cas2Tas(cas=cas, delta=delta, sigma=sigma)
            FL = h / 100
            ff = []
            for mass in massList:
                [cas, speedUpdated] = self.ARPM.cruiseSpeed(
                    theta=theta,
                    delta=delta,
                    h=H_m,
                    mass=mass,
                    DeltaTemp=DeltaTemp,
                    speedSchedule_default=[Vcr1, Vcr2, Mcr],
                    applyLimits=False,
                )
                tas = atm.cas2Tas(cas=cas, delta=delta, sigma=sigma)
                CL = self.CL(tas=tas, sigma=sigma, mass=mass)
                CD = self.CD(CL=CL, config="CR")
                Drag = self.D(tas=tas, sigma=sigma, CD=CD)
                Thrust = Drag
                ff.append(self.ff(flightPhase="Cruise", v=tas, h=H_m, T=Thrust) * 60)

            TAS_CR_complet.append(conv.ms2kt(tas_nominal))
            FF_CR_LO_complet.append(ff[0])
            FF_CR_NOM_complet.append(ff[1])
            FF_CR_HI_complet.append(ff[2])

        CRList = [TAS_CR_complet, FF_CR_LO_complet, FF_CR_NOM_complet, FF_CR_HI_complet]

        return CRList

    def PTF_climb(self, massList, altitudeList, DeltaTemp):
        """This function calculates the BADA3 PTF data in CLIMB

        :param massList: list of aircraft mass [kg]
        :param altitudeList: aircraft altitude list [ft]
        :param DeltaTemp: deviation from ISA temperature [K]
        :type massList: list.
        :type altitudeList: list of int.
        :type DeltaTemp: float.
        :returns: list of PTF CLIMB data [-]
        :rtype: list
        """

        TAS_CL_complet = []
        ROCD_CL_LO_complet = []
        ROCD_CL_NOM_complet = []
        ROCD_CL_HI_complet = []
        FF_CL_NOM_complet = []

        phase = "cl"
        massNominal = massList[1]

        Vcl1 = self.AC.V1[phase]
        Vcl2 = self.AC.V2[phase]
        Mcl = self.AC.M[phase]

        Vcl1 = min(Vcl1, conv.kt2ms(250))
        crossAlt = atm.crossOver(cas=Vcl2, Mach=Mcl)

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
            FL = h / 100

            ROC = []
            tas_list = []
            ff_list = []
            for mass in massList:
                [cas, speedUpdated] = self.ARPM.climbSpeed(
                    theta=theta,
                    delta=delta,
                    h=H_m,
                    mass=mass,
                    DeltaTemp=DeltaTemp,
                    speedSchedule_default=[Vcl1, Vcl2, Mcl],
                    applyLimits=False,
                )
                tas = atm.cas2Tas(cas=cas, delta=delta, sigma=sigma)
                M = atm.tas2Mach(v=tas, theta=theta)
                CL = self.CL(tas=tas, sigma=sigma, mass=mass)
                config = self.flightEnvelope.getConfig(
                    h=H_m, phase="Climb", v=cas, mass=massNominal, DeltaTemp=DeltaTemp
                )
                CD = self.CD(CL=CL, config=config)
                Drag = self.D(tas=tas, sigma=sigma, CD=CD)
                Thrust = self.Thrust(
                    rating="MCMB", v=tas, h=H_m, config=config, DeltaTemp=DeltaTemp
                )
                ff = self.ff(flightPhase="Climb", v=tas, h=H_m, T=Thrust) * 60

                if H_m < crossAlt:
                    ESF = self.esf(
                        h=H_m, flightEvolution="constCAS", M=M, DeltaTemp=DeltaTemp
                    )
                else:
                    ESF = self.esf(
                        h=H_m, flightEvolution="constM", M=M, DeltaTemp=DeltaTemp
                    )

                # I think this should use all config, not just for nominal weight
                ROC_val = (
                    conv.m2ft(
                        self.ROCD(
                            h=H_m,
                            T=Thrust,
                            D=Drag,
                            v=tas,
                            mass=mass,
                            ESF=ESF,
                            DeltaTemp=DeltaTemp,
                            reducedPower=True,
                        )
                    )
                    * 60
                )

                if ROC_val < 0:
                    ROC_val = float("Nan")

                ROC.append(ROC_val)
                tas_list.append(tas)
                ff_list.append(ff)

            TAS_CL_complet.append(conv.ms2kt(tas_list[1]))
            ROCD_CL_LO_complet.append(ROC[0])
            ROCD_CL_NOM_complet.append(ROC[1])
            ROCD_CL_HI_complet.append(ROC[2])
            FF_CL_NOM_complet.append(ff_list[1])

        CLList = [
            TAS_CL_complet,
            ROCD_CL_LO_complet,
            ROCD_CL_NOM_complet,
            ROCD_CL_HI_complet,
            FF_CL_NOM_complet,
        ]

        return CLList

    def PTF_descent(self, massList, altitudeList, DeltaTemp):
        """This function calculates the BADA3 PTF data in DESCENT

        :param massList: list of aircraft mass [kg]
        :param altitudeList: aircraft altitude list [ft]
        :param DeltaTemp: deviation from ISA temperature [K]
        :type massList: list.
        :type altitudeList: list of int.
        :type DeltaTemp: float.
        :returns: list of PTF DESCENT data [-]
        :rtype: list
        """

        TAS_DES_complet = []
        ROCD_DES_NOM_complet = []
        FF_DES_NOM_complet = []

        phase = "des"
        massNominal = massList[1]

        Vdes1 = self.AC.V1[phase]
        Vdes2 = self.AC.V2[phase]
        Mdes = self.AC.M[phase]

        Vdes1 = min(Vdes1, conv.kt2ms(250))
        crossAlt = atm.crossOver(cas=Vdes2, Mach=Mdes)

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
            [cas, speedUpdated] = self.ARPM.descentSpeed(
                theta=theta,
                delta=delta,
                h=H_m,
                mass=massNominal,
                DeltaTemp=DeltaTemp,
                speedSchedule_default=[Vdes1, Vdes2, Mdes],
                applyLimits=False,
            )
            tas_nominal = atm.cas2Tas(cas=cas, delta=delta, sigma=sigma)
            M = atm.tas2Mach(v=tas_nominal, theta=theta)
            FL = h / 100

            config = self.flightEnvelope.getConfig(
                h=H_m, phase="Descent", v=cas, mass=massNominal, DeltaTemp=DeltaTemp
            )

            CL = self.CL(tas=tas_nominal, sigma=sigma, mass=massNominal)
            CD = self.CD(CL=CL, config=config)
            Drag = self.D(tas=tas_nominal, sigma=sigma, CD=CD)

            if self.AC.engineType == "PISTON" or self.AC.engineType == "ELECTRIC":
                # PISTON  and ELECTRIC uses LIDL throughout the whole descent phase
                Thrust_nominal = self.Thrust(
                    rating="LIDL",
                    v=tas_nominal,
                    h=H_m,
                    config="CR",
                    DeltaTemp=DeltaTemp,
                )
                ff_nominal = (
                    self.ff(
                        flightPhase="Descent",
                        v=tas_nominal,
                        h=H_m,
                        T=Thrust_nominal,
                        config="CR",
                        adapted=False,
                    )
                    * 60
                )
            else:
                Thrust_nominal = self.Thrust(
                    rating="LIDL",
                    v=tas_nominal,
                    h=H_m,
                    config=config,
                    DeltaTemp=DeltaTemp,
                )
                ff_nominal = (
                    self.ff(
                        flightPhase="Descent",
                        v=tas_nominal,
                        h=H_m,
                        T=Thrust_nominal,
                        config=config,
                        adapted=False,
                    )
                    * 60
                )

            if H_m < crossAlt:
                ESF = self.esf(
                    h=H_m, flightEvolution="constCAS", M=M, DeltaTemp=DeltaTemp
                )
            else:
                ESF = self.esf(
                    h=H_m, flightEvolution="constM", M=M, DeltaTemp=DeltaTemp
                )

            ROCD = -1 * (
                conv.m2ft(
                    self.ROCD(
                        h=H_m,
                        T=Thrust_nominal,
                        D=Drag,
                        v=tas_nominal,
                        mass=massNominal,
                        ESF=ESF,
                        DeltaTemp=DeltaTemp,
                    )
                )
                * 60
            )

            TAS_DES_complet.append(conv.ms2kt(tas_nominal))
            ROCD_DES_NOM_complet.append(ROCD)
            FF_DES_NOM_complet.append(ff_nominal)

        DESList = [TAS_DES_complet, ROCD_DES_NOM_complet, FF_DES_NOM_complet]

        return DESList


class Bada3Aircraft(BADA3):
    """This class implements the BADA3 performance model following the BADA3 manual.

    :param filePath: path to the BADA3 ascii formatted file.
    :param acName: ICAO aircraft designation
    :type filePath: str.
    :type acName: str
    """

    def __init__(self, badaVersion, acName, filePath=None, allData=None):
        super().__init__(self)

        self.APFavailable = False
        self.OPFavailable = False
        self.ACModelAvailable = False
        self.ACinSynonymFile = False

        self.BADAFamilyName = "BADA3"
        self.BADAFamily = BadaFamily(BADA3=True)
        self.BADAVersion = badaVersion

        if filePath == None:
            self.filePath = configuration.getAircraftPath()
        else:
            self.filePath = filePath

        # check if the aircraft is in the allData dataframe data
        if allData is not None and acName in allData["acName"].values:
            filtered_df = allData[allData["acName"] == acName]

            self.acName = Parser.safe_get(filtered_df, "acName", None)
            self.xmlFiles = Parser.safe_get(filtered_df, "xmlFiles", None)

            self.modificationDateOPF = Parser.safe_get(
                filtered_df, "modificationDateOPF", None
            )
            self.modificationDateAPF = Parser.safe_get(
                filtered_df, "modificationDateAPF", None
            )

            self.ICAO = Parser.safe_get(filtered_df, "ICAO", None)
            self.numberOfEngines = Parser.safe_get(filtered_df, "numberOfEngines", None)
            self.engineType = Parser.safe_get(filtered_df, "engineType", None)
            self.engines = Parser.safe_get(filtered_df, "engines", None)
            self.WTC = Parser.safe_get(filtered_df, "WTC", None)
            self.mass = Parser.safe_get(filtered_df, "mass", None)

            self.MTOW = Parser.safe_get(filtered_df, "MTOW", None)
            self.OEW = Parser.safe_get(filtered_df, "OEW", None)
            self.MPL = Parser.safe_get(filtered_df, "MPL", None)
            self.MREF = Parser.safe_get(filtered_df, "MREF", None)
            self.VMO = Parser.safe_get(filtered_df, "VMO", None)
            self.MMO = Parser.safe_get(filtered_df, "MMO", None)
            self.hmo = Parser.safe_get(filtered_df, "hmo", None)
            self.Hmax = Parser.safe_get(filtered_df, "Hmax", None)
            self.tempGrad = Parser.safe_get(filtered_df, "tempGrad", None)

            self.S = Parser.safe_get(filtered_df, "S", None)
            self.Clbo = Parser.safe_get(filtered_df, "Clbo", None)
            self.k = Parser.safe_get(filtered_df, "k", None)
            self.Vstall = Parser.safe_get(filtered_df, "Vstall", None)
            self.CD0 = Parser.safe_get(filtered_df, "CD0", None)
            self.CD2 = Parser.safe_get(filtered_df, "CD2", None)
            self.HLids = Parser.safe_get(filtered_df, "HLids", None)
            self.Ct = Parser.safe_get(filtered_df, "Ct", None)
            self.CTdeslow = Parser.safe_get(filtered_df, "CTdeslow", None)
            self.CTdeshigh = Parser.safe_get(filtered_df, "CTdeshigh", None)
            self.CTdesapp = Parser.safe_get(filtered_df, "CTdesapp", None)
            self.CTdesld = Parser.safe_get(filtered_df, "CTdesld", None)
            self.HpDes = Parser.safe_get(filtered_df, "HpDes", None)
            self.Cf = Parser.safe_get(filtered_df, "Cf", None)
            self.CfDes = Parser.safe_get(filtered_df, "CfDes", None)
            self.CfCrz = Parser.safe_get(filtered_df, "CfCrz", None)
            self.TOL = Parser.safe_get(filtered_df, "TOL", None)
            self.LDL = Parser.safe_get(filtered_df, "LDL", None)
            self.span = Parser.safe_get(filtered_df, "span", None)
            self.length = Parser.safe_get(filtered_df, "length", None)

            self.V1 = Parser.safe_get(filtered_df, "V1", None)
            self.V2 = Parser.safe_get(filtered_df, "V2", None)
            self.M = Parser.safe_get(filtered_df, "M", None)

            self.GPFdata = Parser.safe_get(filtered_df, "GPFdata", None)

            self.drone = Parser.safe_get(filtered_df, "drone", None)

            self.DeltaCD = Parser.safe_get(filtered_df, "DeltaCD", None)
            self.speedSchedule = Parser.safe_get(filtered_df, "speedSchedule", None)
            self.aeroConfig = Parser.safe_get(filtered_df, "aeroConfig", None)

            self.flightEnvelope = FlightEnvelope(self)
            self.ARPM = ARPM(self)
            self.PTD = PTD(self)
            self.PTF = PTF(self)

        else:
            # read BADA3 GPF file
            GPFDataFrame = Parser.parseGPF(self.filePath, badaVersion)

            # check if SYNONYM file exist
            synonymFile = os.path.join(
                self.filePath, "BADA3", badaVersion, "SYNONYM.NEW"
            )
            synonymFileXML = os.path.join(
                self.filePath, "BADA3", badaVersion, "SYNONYM.xml"
            )
            if os.path.isfile(synonymFile) or os.path.isfile(synonymFileXML):
                self.synonymFileAvailable = True

                self.SearchedACName = Parser.parseSynonym(
                    self.filePath, badaVersion, acName
                )

                if self.SearchedACName == None:
                    # look for file name directly, which consists of added "_" at the end of file
                    fileName = acName
                    while len(fileName) < 6:
                        fileName += "_"
                    self.SearchedACName = fileName
                else:
                    self.ACinSynonymFile = True
            else:
                # if doesn't exist - look for full name based on acName (may not be ICAO designator)
                self.SearchedACName = acName

            # look for either found synonym or original full BADA3 model name designator
            if self.SearchedACName is not None:
                # check for existence of OPF and APF files
                OPFfile = (
                    os.path.join(
                        self.filePath, "BADA3", badaVersion, self.SearchedACName
                    )
                    + ".OPF"
                )
                APFfile = (
                    os.path.join(
                        self.filePath, "BADA3", badaVersion, self.SearchedACName
                    )
                    + ".APF"
                )
                if os.path.isfile(OPFfile):
                    self.OPFavailable = True
                if os.path.isfile(APFfile):
                    self.APFavailable = True

                if self.OPFavailable and self.APFavailable:
                    self.ACModelAvailable = True

                    OPFDataFrame = Parser.parseOPF(
                        self.filePath, badaVersion, self.SearchedACName
                    )
                    APFDataFrame = Parser.parseAPF(
                        self.filePath, badaVersion, self.SearchedACName
                    )

                    OPF_APF_combined_df = Parser.combineOPF_APF(
                        OPFDataFrame, APFDataFrame
                    )
                    combined_df = Parser.combineACDATA_GPF(
                        OPF_APF_combined_df, GPFDataFrame
                    )

                    self.acName = Parser.safe_get(combined_df, "acName", None)
                    self.xmlFiles = Parser.safe_get(combined_df, "xmlFiles", None)

                    self.modificationDateOPF = Parser.safe_get(
                        combined_df, "modificationDateOPF", None
                    )
                    self.modificationDateAPF = Parser.safe_get(
                        combined_df, "modificationDateAPF", None
                    )

                    self.ICAO = Parser.safe_get(combined_df, "ICAO", None)
                    self.numberOfEngines = Parser.safe_get(
                        combined_df, "numberOfEngines", None
                    )
                    self.engineType = Parser.safe_get(combined_df, "engineType", None)
                    self.engines = Parser.safe_get(combined_df, "engines", None)
                    self.WTC = Parser.safe_get(combined_df, "WTC", None)
                    self.mass = Parser.safe_get(combined_df, "mass", None)

                    self.MTOW = Parser.safe_get(combined_df, "MTOW", None)
                    self.OEW = Parser.safe_get(combined_df, "OEW", None)
                    self.MPL = Parser.safe_get(combined_df, "MPL", None)
                    self.MREF = Parser.safe_get(combined_df, "MREF", None)
                    self.VMO = Parser.safe_get(combined_df, "VMO", None)
                    self.MMO = Parser.safe_get(combined_df, "MMO", None)
                    self.hmo = Parser.safe_get(combined_df, "hmo", None)
                    self.Hmax = Parser.safe_get(combined_df, "Hmax", None)
                    self.tempGrad = Parser.safe_get(combined_df, "tempGrad", None)

                    self.S = Parser.safe_get(combined_df, "S", None)
                    self.Clbo = Parser.safe_get(combined_df, "Clbo", None)
                    self.k = Parser.safe_get(combined_df, "k", None)
                    self.Vstall = Parser.safe_get(combined_df, "Vstall", None)
                    self.CD0 = Parser.safe_get(combined_df, "CD0", None)
                    self.CD2 = Parser.safe_get(combined_df, "CD2", None)
                    self.HLids = Parser.safe_get(combined_df, "HLids", None)
                    self.Ct = Parser.safe_get(combined_df, "Ct", None)
                    self.CTdeslow = Parser.safe_get(combined_df, "CTdeslow", None)
                    self.CTdeshigh = Parser.safe_get(combined_df, "CTdeshigh", None)
                    self.CTdesapp = Parser.safe_get(combined_df, "CTdesapp", None)
                    self.CTdesld = Parser.safe_get(combined_df, "CTdesld", None)
                    self.HpDes = Parser.safe_get(combined_df, "HpDes", None)
                    self.Cf = Parser.safe_get(combined_df, "Cf", None)
                    self.CfDes = Parser.safe_get(combined_df, "CfDes", None)
                    self.CfCrz = Parser.safe_get(combined_df, "CfCrz", None)
                    self.TOL = Parser.safe_get(combined_df, "TOL", None)
                    self.LDL = Parser.safe_get(combined_df, "LDL", None)
                    self.span = Parser.safe_get(combined_df, "span", None)
                    self.length = Parser.safe_get(combined_df, "length", None)

                    self.V1 = Parser.safe_get(combined_df, "V1", None)
                    self.V2 = Parser.safe_get(combined_df, "V2", None)
                    self.M = Parser.safe_get(combined_df, "M", None)

                    self.GPFdata = Parser.safe_get(combined_df, "GPFdata", None)

                    self.drone = Parser.safe_get(combined_df, "drone", None)

                    self.DeltaCD = Parser.safe_get(combined_df, "DeltaCD", None)
                    self.speedSchedule = Parser.safe_get(
                        combined_df, "speedSchedule", None
                    )
                    self.aeroConfig = Parser.safe_get(combined_df, "aeroConfig", None)

                    self.flightEnvelope = FlightEnvelope(self)
                    self.ARPM = ARPM(self)
                    self.PTD = PTD(self)
                    self.PTF = PTF(self)

                elif not self.OPFavailable and not self.APFavailable:
                    # search for xml files

                    XMLDataFrame = Parser.parseXML(
                        self.filePath, badaVersion, self.SearchedACName
                    )

                    combined_df = Parser.combineACDATA_GPF(XMLDataFrame, GPFDataFrame)

                    self.acName = Parser.safe_get(combined_df, "acName", None)
                    self.xmlFiles = Parser.safe_get(combined_df, "xmlFiles", None)

                    self.modificationDateOPF = Parser.safe_get(
                        combined_df, "modificationDateOPF", None
                    )
                    self.modificationDateAPF = Parser.safe_get(
                        combined_df, "modificationDateAPF", None
                    )

                    self.ICAO = Parser.safe_get(combined_df, "ICAO", None)
                    self.numberOfEngines = Parser.safe_get(
                        combined_df, "numberOfEngines", None
                    )
                    self.engineType = Parser.safe_get(combined_df, "engineType", None)
                    self.engines = Parser.safe_get(combined_df, "engines", None)
                    self.WTC = Parser.safe_get(combined_df, "WTC", None)
                    self.mass = Parser.safe_get(combined_df, "mass", None)

                    self.MTOW = Parser.safe_get(combined_df, "MTOW", None)
                    self.OEW = Parser.safe_get(combined_df, "OEW", None)
                    self.MPL = Parser.safe_get(combined_df, "MPL", None)
                    self.MREF = Parser.safe_get(combined_df, "MREF", None)
                    self.VMO = Parser.safe_get(combined_df, "VMO", None)
                    self.MMO = Parser.safe_get(combined_df, "MMO", None)
                    self.hmo = Parser.safe_get(combined_df, "hmo", None)
                    self.Hmax = Parser.safe_get(combined_df, "Hmax", None)
                    self.tempGrad = Parser.safe_get(combined_df, "tempGrad", None)

                    self.S = Parser.safe_get(combined_df, "S", None)
                    self.Clbo = Parser.safe_get(combined_df, "Clbo", None)
                    self.k = Parser.safe_get(combined_df, "k", None)
                    self.Vstall = Parser.safe_get(combined_df, "Vstall", None)
                    self.CD0 = Parser.safe_get(combined_df, "CD0", None)
                    self.CD2 = Parser.safe_get(combined_df, "CD2", None)
                    self.HLids = Parser.safe_get(combined_df, "HLids", None)
                    self.Ct = Parser.safe_get(combined_df, "Ct", None)
                    self.CTdeslow = Parser.safe_get(combined_df, "CTdeslow", None)
                    self.CTdeshigh = Parser.safe_get(combined_df, "CTdeshigh", None)
                    self.CTdesapp = Parser.safe_get(combined_df, "CTdesapp", None)
                    self.CTdesld = Parser.safe_get(combined_df, "CTdesld", None)
                    self.HpDes = Parser.safe_get(combined_df, "HpDes", None)
                    self.Cf = Parser.safe_get(combined_df, "Cf", None)
                    self.CfDes = Parser.safe_get(combined_df, "CfDes", None)
                    self.CfCrz = Parser.safe_get(combined_df, "CfCrz", None)
                    self.TOL = Parser.safe_get(combined_df, "TOL", None)
                    self.LDL = Parser.safe_get(combined_df, "LDL", None)
                    self.span = Parser.safe_get(combined_df, "span", None)
                    self.length = Parser.safe_get(combined_df, "length", None)

                    self.V1 = Parser.safe_get(combined_df, "V1", None)
                    self.V2 = Parser.safe_get(combined_df, "V2", None)
                    self.M = Parser.safe_get(combined_df, "M", None)

                    self.GPFdata = Parser.safe_get(combined_df, "GPFdata", None)

                    self.drone = Parser.safe_get(combined_df, "drone", None)

                    self.DeltaCD = Parser.safe_get(combined_df, "DeltaCD", None)
                    self.speedSchedule = Parser.safe_get(
                        combined_df, "speedSchedule", None
                    )
                    self.aeroConfig = Parser.safe_get(combined_df, "aeroConfig", None)

                    self.flightEnvelope = FlightEnvelope(self)
                    self.ARPM = ARPM(self)
                    self.PTD = PTD(self)
                    self.PTF = PTF(self)

                else:
                    # AC name cannot be found
                    raise ValueError(acName + " Cannot be found")

    def __str__(self):
        return f"(BADA3, AC_name: {self.acName}, searched_AC_name: {self.SearchedACName}, model_ICAO: {self.ICAO}, ID: {id(self.AC)})"
