"""Generic BADA3 aircraft performance module."""

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
from pyBADA.aircraft import Airplane, BadaFamily, Bada


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


class Parser:
    """This class implements the BADA3 parsing mechanism to parse APF, OPF and
    GPF BADA3 files."""

    def __init__(self):
        pass

    @staticmethod
    def parseXML(filePath, acName):
        """Parses a BADA3 XML formatted file for aircraft performance data.

        :param filePath: Path to the XML file containing BADA data.
        :param acName: Name of the aircraft for which data is being parsed
            from the XML file.
        :type filePath: str
        :type acName: str
        :raises IOError: If the file cannot be found or read.
        :raises ValueError: If the BADA version is unsupported or if parsing
            fails.
        :returns: A pandas DataFrame containing the parsed aircraft
            performance data.
        :rtype: pd.DataFrame
        """

        filename = os.path.join(filePath, acName, acName) + ".xml"

        try:
            tree = ET.parse(filename)
            root = tree.getroot()
        except Exception:
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

            if (
                LGDN is not None
            ):  # Landing gear NOT allowed in clean configuration
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
            "MMO": [MMO],
            "VMO": [VMO],
            "MTOW": [MTOW],
            "OEW": [OEW],
            "MPL": [MPL],
            "TOL": [TOL],
            "LDL": [LDL],
            "span": [span],
            "length": [length],
            "aeroConfig": [aeroConfig],
            "speedSchedule": [speedSchedule],
            "V1": [V1],
            "V2": [V2],
            "M": [M],
            "xmlFiles": [xmlFiles],
        }
        df_single = pd.DataFrame(data)

        return df_single

    @staticmethod
    def findData(f):
        """Searches for specific data lines in an open file stream.

        :param f: An open file object from which lines are read.
        :type f: file object
        :returns: A tuple containing the file object and a parsed line split
            into a list, or None if no relevant line is found.
        :rtype: tuple(file object, list of str or None) This function reads
            the file line by line until it finds a line that starts with "CD".
            Once found, the line is stripped of extra spaces, split into a
            list, and returned. If no such line is found, it returns None for
            the line.
        """

        line = f.readline()
        while line is not None and not line.startswith("CD"):
            line = f.readline()

        if line is None:
            return f, None

        line = " ".join(line.split())
        line = line.strip().split(" ")
        return f, line

    @staticmethod
    def parseOPF(filePath, acName):
        """Parses a BADA3 OPF (Operational Performance File) ASCII formatted
        file for aircraft performance data.

        :param filePath: Path to the BADA3 OPF ASCII formatted file.
        :param acName: ICAO aircraft designation (e.g., 'A320').
        :type filePath: str
        :type acName: str
        :raises IOError: If the file cannot be opened or read.
        :returns: A pandas DataFrame containing the parsed aircraft
            performance data.
        :rtype: pd.DataFrame
        """

        filename = (
            os.path.join(
                filePath,
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
                        engines = (
                            line.split("with")[1].split("engines")[0].strip()
                        )
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
            "CTdesapp": [CTdesapp],
            "CTdesld": [CTdesld],
            "HpDes": [HpDes],
            "CfDes": [CfDes],
            "CfCrz": [CfCrz],
            "Cf": [Cf],
            "hmo": [hmo],
            "Hmax": [Hmax],
            "tempGrad": [tempGrad],
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
        }
        df_single = pd.DataFrame(data)

        return df_single

    @staticmethod
    def parseAPF(filePath, acName):
        """Parses a BADA3 APF ASCII formatted file for aircraft performance
        data.

        :param filePath: Path to the BADA3 APF ASCII formatted file.
        :param acName: ICAO aircraft designation (e.g., 'A320').
        :type filePath: str
        :type acName: str
        :raises IOError: If the file cannot be opened or read.
        :returns: A pandas DataFrame containing the parsed aircraft
            performance data.
        :rtype: pd.DataFrame
        """

        filename = os.path.join(filePath, acName) + ".APF"

        dataLines = list()
        with open(filename, "r", encoding="latin-1") as f:
            while True:
                line = f.readline()

                if line.startswith("CC"):
                    if "Modification_date" in line:
                        data = line.split(":")[1].strip().split(" ")
                        modificationDateAPF = " ".join(
                            [data[0], data[1], data[2]]
                        )
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
        dataLines.pop(
            0
        )  # remove first line that does not contain usefull data

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
        """Combines data from OPF and APF DataFrames.

        :param OPFDataFrame: DataFrame containing parsed data from the OPF
            file.
        :param APFDataFrame: DataFrame containing parsed data from the APF
            file.
        :type OPFDataFrame: pd.DataFrame
        :type APFDataFrame: pd.DataFrame
        :returns: A single DataFrame combining both OPF and APF data.
        :rtype: pd.DataFrame
        """

        # Combine data with GPF data (temporary solution)
        combined_df = pd.concat(
            [
                OPFDataFrame.reset_index(drop=True),
                APFDataFrame.reset_index(drop=True),
            ],
            axis=1,
        )

        return combined_df

    @staticmethod
    def readSynonym(filePath):
        """Reads a BADA3 SYNONYM.NEW ASCII file and returns a dictionary of
        model-synonym pairs.

        :param filePath: Path to the directory containing BADA3 files.
        :type filePath: str
        :returns: A dictionary where the keys are aircraft models and the
            values are the corresponding file names.
        :rtype: dict
        """

        filename = os.path.join(filePath, "SYNONYM.NEW")

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
    def readSynonymXML(filePath):
        """Reads a BADA3 SYNONYM.xml file and returns a dictionary of model-
        synonym pairs.

        :param filePath: Path to the directory containing BADA3 files.
        :type filePath: str
        :returns: A dictionary where the keys are aircraft models (codes) and
            the values are the corresponding file names.
        :rtype: dict
        :raises IOError: If the XML file is not found or cannot be read. This
            function parses the 'SYNONYM.xml' file to extract aircraft model
            codes and their associated file names. If the XML file is not
            found or is improperly formatted, an IOError is raised.
        """

        filename = os.path.join(filePath, "SYNONYM.xml")

        # synonym - file name pair dictionary
        synonym_fileName = {}

        if os.path.isfile(filename):
            try:
                tree = ET.parse(filename)
                root = tree.getroot()
            except Exception:
                raise IOError(filename + " not found or in correct format")

            for child in root.iter("SYN"):
                code = child.find("code").text
                manufacturer = child.find("manu").text
                file = child.find("file").text
                ICAO = child.find("ICAO").text

                synonym_fileName[code] = file

        return synonym_fileName

    @staticmethod
    def parseSynonym(filePath, acName):
        """Parses either the ASCII or XML synonym file and returns the file
        name corresponding to the aircraft.

        :param filePath: Path to the directory containing BADA3 files.
        :param acName: ICAO aircraft designation for which the file name is
            needed.
        :type filePath: str
        :type acName: str
        :returns: The file name corresponding to the aircraft, or None if not
            found.
        :rtype: str or None This function first attempts to read the aircraft
            synonym from the ASCII file ('SYNONYM.NEW'). If the synonym is not
            found, it then tries to read the XML version ('SYNONYM.xml'). It
            returns the associated file name or None if the aircraft synonym
            is not found.
        """

        synonym_fileName = Parser.readSynonym(filePath)

        # if ASCI synonym does not exist, try XML synonym file
        if not synonym_fileName:
            synonym_fileName = Parser.readSynonymXML(filePath)

        if synonym_fileName and acName in synonym_fileName:
            fileName = synonym_fileName[acName]
            return fileName
        else:
            return None

    @staticmethod
    def readGPF(filePath):
        """Parses a BADA3 GPF ASCII formatted file.

        :param filePath: Path to the directory containing BADA3 files.
        :type filePath: str
        :raises IOError: If the GPF file cannot be opened or read.
        :returns: A list of dictionaries, each containing GPF parameters like
            engine type, flight phase, and parameter values.
        :rtype: list of dict
        """

        filename = os.path.join(filePath, "BADA.GPF")

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
    def readGPFXML(filePath):
        """Parses a BADA3 GPF XML formatted file.

        :param filePath: Path to the directory containing BADA3 files.
        :type filePath: str
        :raises IOError: If the XML file is not found or cannot be read.
        :returns: A list of dictionaries, each containing GPF parameters such
            as engine type, flight phase, and performance values.
        :rtype: list of dict This function reads the 'GPF.xml' file and
            extracts general performance parameters for the aircraft,
            including maximum acceleration, bank angles, thrust coefficients,
            speed limits, and more. It parses the XML structure and returns a
            list of dictionaries representing these parameters.
        """

        filename = os.path.join(filePath, "GPF.xml")

        GPFparamList = list()

        if os.path.isfile(filename):
            try:
                tree = ET.parse(filename)
                root = tree.getroot()
            except Exception:
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
                    "value": float(
                        AngBank.find("Nom").find("Civ").find("ToLd").text
                    ),
                    "engine": allEngines,
                    "phase": ["to", "ld"],
                    "flight": ["civ"],
                }
            )
            GPFparamList.append(
                {
                    "name": "ang_bank_nom",
                    "value": float(
                        AngBank.find("Nom").find("Civ").find("Others").text
                    ),
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
                    "value": float(
                        AngBank.find("Max").find("Civ").find("ToLd").text
                    ),
                    "engine": allEngines,
                    "phase": ["to", "ld"],
                    "flight": ["civ"],
                }
            )
            GPFparamList.append(
                {
                    "name": "ang_bank_max",
                    "value": float(
                        AngBank.find("Max").find("Civ").find("Hold").text
                    ),
                    "engine": allEngines,
                    "phase": ["hold"],
                    "flight": ["civ"],
                }
            )
            GPFparamList.append(
                {
                    "name": "ang_bank_max",
                    "value": float(
                        AngBank.find("Max").find("Civ").find("Others").text
                    ),
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
                HmaxList[phase.find("Phase").text] = float(
                    phase.find("Hmax").text
                )

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
    def parseGPF(filePath):
        """Parses a BADA3 (GPF) from either ASCII or XML format.

        :param filePath: Path to the directory containing BADA3 files.
        :type filePath: str
        :returns: A pandas DataFrame containing GPF data.
        :rtype: pd.DataFrame
        """

        GPFdata = Parser.readGPF(filePath)

        # if ASCI GPF does not exist, try XML GPF file
        if not GPFdata:
            GPFdata = Parser.readGPFXML(filePath)

        # Single row dataframe
        data = {"GPFdata": [GPFdata]}
        df_single = pd.DataFrame(data)

        return df_single

    @staticmethod
    def getGPFValue(GPFdata, name, engine="JET", phase="cr", flight="civ"):
        """Retrieves the value of a specified GPF parameter based on engine
        type, flight phase, and flight type.

        :param GPFdata: List of dictionaries containing GPF parameters.
        :param name: Name of the GPF parameter to retrieve.
        :param engine: Engine type to filter by (e.g., 'JET', 'TURBOPROP',
            'PISTON', 'ELECTRIC'). Default is 'JET'.
        :param phase: Flight phase to filter by (e.g., 'cr', 'cl', 'des').
            Default is 'cr'.
        :param flight: Flight type to filter by ('civ' or 'mil'). Default is
            'civ'.
        :type GPFdata: list
        :type name: str
        :type engine: str
        :type phase: str
        :type flight: str
        :returns: The value of the specified GPF parameter or None if not
            found.
        :rtype: float or None
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
        """
        Combines two DataFrames: one containing aircraft-specific data (ACData) and another containing (GPF) data.

        :param ACDataFrame: DataFrame containing parsed aircraft data.
        :param GPFDataframe: DataFrame containing parsed GPF data.
        :type ACDataFrame: pd.DataFrame
        :type GPFDataframe: pd.DataFrame
        :returns: A combined DataFrame containing both ACData and GPF data.
        :rtype: pd.DataFrame
        """

        # Combine data with GPF data (temporary solution)
        combined_df = pd.concat(
            [
                ACDataFrame.reset_index(drop=True),
                GPFDataframe.reset_index(drop=True),
            ],
            axis=1,
        )

        return combined_df

    @staticmethod
    def parseAll(badaVersion, filePath=None):
        """Parses all BADA3 formatted files and combines them into a final
        DataFrame.

        :param badaVersion: BADA version being used.
        :param filePath: Path to the BADA3 formatted files. If not provided,
            the default path is used.
        :type badaVersion: str
        :type filePath: str, optional
        :returns: A pandas DataFrame containing all parsed BADA3 data.
        :rtype: pd.DataFrame
        :raises IOError: If any of the required files cannot be opened or
            read.
        """

        if filePath is None:
            filePath = configuration.getBadaVersionPath(
                badaFamily="BADA3", badaVersion=badaVersion
            )
        else:
            filePath = filePath

        # parsing GPF file
        GPFparsedDataframe = Parser.parseGPF(filePath)

        # try to get subfolders, if they exist
        # get names of all the folders in the main BADA model folder to search for XML files
        subfolders = configuration.list_subfolders(filePath)

        if not subfolders:
            # use APF and OPF files
            merged_df = pd.DataFrame()

            # get synonym-filename pairs
            synonym_fileName = Parser.readSynonym(filePath)

            for synonym in synonym_fileName:
                file = synonym_fileName[synonym]

                # parse the original data of a model
                OPFDataFrame = Parser.parseOPF(filePath, file)
                APFDataFrame = Parser.parseAPF(filePath, file)

                df = Parser.combineOPF_APF(OPFDataFrame, APFDataFrame)

                # rename acName in the dateaframe to match the synonym model name
                df.at[0, "acName"] = synonym

                # Combine data with GPF data (temporary solution)
                combined_df = Parser.combineACDATA_GPF(df, GPFparsedDataframe)

                # Merge DataFrames
                merged_df = pd.concat(
                    [merged_df, combined_df], ignore_index=True
                )

            return merged_df

        else:
            # use xml files inside those subfolders
            merged_df = pd.DataFrame()

            # get synonym-filename pairs
            synonym_fileName = Parser.readSynonymXML(filePath)

            for synonym in synonym_fileName:
                file = synonym_fileName[synonym]

                if file in subfolders:
                    # parse the original XML of a model
                    df = Parser.parseXML(filePath, file)

                    # rename acName in the dateaframe to match the synonym model name
                    df.at[0, "acName"] = synonym

                    # Combine data with GPF data (temporary solution)
                    combined_df = Parser.combineACDATA_GPF(
                        df, GPFparsedDataframe
                    )

                    # Merge DataFrames
                    merged_df = pd.concat(
                        [merged_df, combined_df], ignore_index=True
                    )

            return merged_df


class BADA3(Airplane, Bada):
    """This class implements the part of BADA3 performance model that will be
    used in other classes following the BADA3 manual.

    :param AC: Aircraft object {BADA3}.
    :type AC: bada3Aircraft.
    """

    def __init__(self, AC):
        super().__init__()
        self.AC = AC

    def CL(self, sigma, mass, tas, nz=1.0):
        """Computes the lift coefficient for the aircraft.

        :param sigma: Normalized air density [-].
        :param mass: Aircraft mass in kilograms [kg].
        :param tas: True airspeed in meters per second [m/s].
        :param nz: Load factor [-], default is 1.0 (straight and level
            flight).
        :type sigma: float
        :type mass: float
        :type tas: float
        :type nz: float
        :returns: Lift coefficient [-].
        :rtype: float
        """

        return (
            2
            * mass
            * const.g
            * nz
            / (sigma * const.rho_0 * tas * tas * self.AC.S)
        )

    def CD(
        self,
        CL,
        config,
        expedite=False,
        speedBrakes={"deployed": False, "value": 0.03},
    ):
        """Computes the drag coefficient based on the lift coefficient and
        aircraft configuration.

        :param CL: Lift coefficient [-].
        :param config: Aircraft aerodynamic configuration (e.g., 'CR', 'IC',
            'TO', 'AP', 'LD').
        :param expedite: Flag indicating if expedite descent is used (default
            is False).
        :param speedBrakes: Dictionary indicating if speed brakes are deployed
            and their effect.
        :type CL: float
        :type config: str
        :type expedite: bool
        :type speedBrakes: dict
        :returns: Drag coefficient [-].
        :rtype: float
        :raises: ValueError if an invalid configuration is provided.
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
                CD = self.AC.CD0[HLid_CR][LG_CR] + self.AC.CD2[HLid_CR][
                    LG_CR
                ] * (CL * CL)
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
            C_des_exp = Parser.getGPFValue(
                self.AC.GPFdata, "C_des_exp", phase="des"
            )
            CD = CD * C_des_exp

        return CD

    def D(self, sigma, tas, CD):
        """Computes the aerodynamic drag force.

        :param sigma: Normalized air density [-].
        :param tas: True airspeed in meters per second [m/s].
        :param CD: Drag coefficient [-].
        :type sigma: float
        :type tas: float
        :type CD: float
        :returns: Aerodynamic drag in Newtons [N].
        :rtype: float
        """

        return 0.5 * sigma * const.rho_0 * tas * tas * self.AC.S * CD

    def L(self, sigma, tas, CL):
        """Computes the aerodynamic lift force.

        :param sigma: Normalized air density [-].
        :param tas: True airspeed in meters per second [m/s].
        :param CL: Lift coefficient [-].
        :type sigma: float
        :type tas: float
        :type CL: float
        :returns: Aerodynamic lift in Newtons [N].
        :rtype: float
        """

        return 0.5 * sigma * const.rho_0 * tas * tas * self.AC.S * CL

    def Thrust(self, h, DeltaTemp, rating, v, config, **kwargs):
        """Computes the aircraft thrust based on engine rating and flight
        conditions.

        :param rating: Engine rating ('MCMB', 'MCRZ', 'MTKF', 'LIDL',
            'ADAPTED').
        :param h: Altitude in meters [m].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param v: True airspeed (TAS) in meters per second [m/s].
        :param config: Aircraft aerodynamic configuration (e.g., 'CR', 'IC',
            'TO', 'AP', 'LD').
        :type rating: str
        :type h: float
        :type DeltaTemp: float
        :type v: float
        :type config: str
        :returns: Thrust in Newtons [N].
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
                h=h,
                DeltaTemp=DeltaTemp,
                ROCD=ROCD,
                mass=mass,
                v=v,
                acc=acc,
                Drag=Drag,
            )

        else:
            T = float("Nan")

        return T

    def TAdapted(self, h, DeltaTemp, ROCD, mass, v, acc, Drag):
        """Computes adapted thrust for non-standard flight conditions (e.g.,
        climb, acceleration).

        :param h: Altitude in meters [m].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param ROCD: Rate of climb or descent in meters per second [m/s].
        :param mass: Aircraft mass in kilograms [kg].
        :param v: True airspeed (TAS) in meters per second [m/s].
        :param acc: Acceleration in meters per second squared [m/s²].
        :param Drag: Aerodynamic drag in Newtons [N].
        :type h: float
        :type DeltaTemp: float
        :type ROCD: float
        :type mass: float
        :type v: float
        :type acc: float
        :type Drag: float
        :returns: Adapted thrust in Newtons [N].
        :rtype: float
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        tau_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)
        Tadapted = ROCD * mass * const.g * tau_const / v + mass * acc + Drag

        return Tadapted

    def TMax(self, h, DeltaTemp, rating, v):
        """Computes the maximum thrust based on engine type, altitude, and
        temperature deviation.

        :param h: Altitude in meters [m].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param rating: Engine rating ('MCMB', 'MCRZ', 'MTKF').
        :param v: True airspeed (TAS) in meters per second [m/s].
        :type h: float
        :type DeltaTemp: float
        :type rating: str
        :type v: float
        :returns: Maximum thrust in Newtons [N].
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
            return TMax * Parser.getGPFValue(
                self.AC.GPFdata, "C_th_cr", phase="cr"
            )

    def TDes(self, h, DeltaTemp, v, config):
        """Computes descent thrust based on altitude, temperature deviation,
        and configuration.

        :param h: Altitude in meters [m].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param v: True airspeed (TAS) in meters per second [m/s].
        :param config: Aircraft aerodynamic configuration (e.g., 'CR', 'IC',
            'TO', 'AP', 'LD').
        :type h: float
        :type DeltaTemp: float
        :type v: float
        :type config: str
        :returns: Descent thrust in Newtons [N].
        :rtype: float
        """

        H_max_app = Parser.getGPFValue(
            self.AC.GPFdata, "H_max_app", phase="app"
        )

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
        """Computes the nominal fuel flow based on airspeed and thrust.

        :param v: True airspeed (TAS) in meters per second [m/s].
        :param T: Thrust in Newtons [N].
        :type v: float
        :type T: float
        :returns: Nominal fuel flow in kilograms per second [kg/s].
        :rtype: float
        """

        if self.AC.engineType == "JET":
            eta = (
                self.AC.Cf[0]
                * (1 + conv.ms2kt(v) / self.AC.Cf[1])
                / (1000 * 60)
            )
            ffnom = eta * T

        elif self.AC.engineType == "TURBOPROP":
            eta = (
                self.AC.Cf[0]
                * (1 - conv.ms2kt(v) / self.AC.Cf[1])
                * (conv.ms2kt(v) / 1000)
                / (1000 * 60)
            )
            ffnom = eta * T

        elif (
            self.AC.engineType == "PISTON" or self.AC.engineType == "ELECTRIC"
        ):
            ffnom = self.AC.Cf[0] / 60

        return ffnom

    def ffMin(self, h):
        """Computes the minimum fuel flow based on altitude.

        :param h: Altitude in meters [m].
        :type h: float
        :returns: Minimum fuel flow in kilograms per second [kg/s].
        :rtype: float
        """

        if self.AC.engineType == "JET" or self.AC.engineType == "TURBOPROP":
            ffmin = (
                self.AC.CfDes[0] * (1 - (conv.m2ft(h)) / self.AC.CfDes[1]) / 60
            )
        elif (
            self.AC.engineType == "PISTON" or self.AC.engineType == "ELECTRIC"
        ):
            ffmin = self.AC.CfDes[0] / 60  # Cf3 param

        return ffmin

    def ff(self, h, v, T, config=None, flightPhase=None, adapted=False):
        """Computes the fuel flow based on flight phase and current flight
        conditions.

        :param h: Altitude in meters [m].
        :param v: True airspeed (TAS) in meters per second [m/s].
        :param T: Thrust in Newtons [N].
        :param config: Aircraft aerodynamic configuration (e.g., 'CR', 'AP',
            'LD'). Optional.
        :param flightPhase: Flight phase (e.g., 'Climb', 'Cruise', 'Descent').
            Optional.
        :param adapted: If True, computes fuel flow for adapted thrust.
            Default is False.
        :type h: float
        :type v: float
        :type T: float
        :type config: str, optional
        :type flightPhase: str, optional
        :type adapted: bool, optional
        :returns: Fuel flow in kilograms per second [kg/s].
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
        """Computes the reduced climb power coefficient based on altitude,
        mass, and temperature deviation.

        :param h: Altitude in meters [m].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param mass: Aircraft mass in kilograms [kg].
        :type h: float
        :type DeltaTemp: float
        :type mass: float
        :returns: Reduced climb power coefficient [-].
        :rtype: float
        """

        hMax = self.flightEnvelope.maxAltitude(mass=mass, DeltaTemp=DeltaTemp)
        mMax = self.AC.mass["maximum"]
        mMin = self.AC.mass["minimum"]

        ep = 1e-6  # floating point precision
        if (h + ep) < 0.8 * hMax:
            if self.AC.engineType == "JET":
                CRed = Parser.getGPFValue(
                    self.AC.GPFdata,
                    "C_red_jet",
                    engine=self.AC.engineType,
                    phase="cl",
                )
            elif self.AC.engineType == "TURBOPROP":
                CRed = Parser.getGPFValue(
                    self.AC.GPFdata,
                    "C_red_turbo",
                    engine="TURBOPROP",
                    phase="cl",
                )
            elif self.AC.engineType == "PISTON":
                CRed = Parser.getGPFValue(
                    self.AC.GPFdata,
                    "C_red_piston",
                    engine="PISTON",
                    phase="cl",
                )
            elif self.AC.engineType == "ELECTRIC":
                CRed = Parser.getGPFValue(
                    self.AC.GPFdata,
                    "C_red_elec",
                    engine="ELECTRIC",
                    phase="cl",
                )
        else:
            CRed = 0

        CPowRed = 1 - CRed * (mMax - mass) / (mMax - mMin)
        return CPowRed

    def ROCD(self, T, D, v, mass, ESF, h, DeltaTemp, reducedPower=False):
        """Computes the rate of climb or descent (ROCD) based on thrust, drag,
        airspeed, and other flight parameters.

        :param T: Aircraft thrust in Newtons [N].
        :param D: Aircraft drag in Newtons [N].
        :param v: True airspeed (TAS) in meters per second [m/s].
        :param mass: Aircraft mass in kilograms [kg].
        :param ESF: Energy share factor [-].
        :param h: Altitude in meters [m].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param reducedPower: Whether to account for reduced power in the
            calculation. Default is False.
        :type T: float
        :type D: float
        :type v: float
        :type mass: float
        :type ESF: float
        :type h: float
        :type DeltaTemp: float
        :type reducedPower: bool, optional
        :returns: Rate of climb or descent in meters per second [m/s].
        :rtype: float
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        temp = theta * const.temp_0

        CPowRed = 1.0
        if reducedPower:
            CPowRed = self.reducedPower(h=h, mass=mass, DeltaTemp=DeltaTemp)

        ROCD = (
            ((temp - DeltaTemp) / temp)
            * (T - D)
            * v
            * ESF
            * CPowRed
            / (mass * const.g)
        )

        return ROCD


class FlightEnvelope(BADA3):
    """This class is a BADA3 aircraft subclass and implements the flight
    envelope caclulations following the BADA3 manual.

    :param AC: Aircraft object {BADA3}.
    :type AC: bada3Aircraft.
    """

    def __init__(self, AC):
        super().__init__(AC)

    def maxAltitude(self, mass, DeltaTemp):
        """Computes the maximum altitude for a given aircraft mass and
        deviation from ISA.

        :param mass: Actual aircraft mass in kilograms [kg].
        :param DeltaTemp: Deviation from International Standard Atmosphere
            (ISA) temperature in Kelvin [K].
        :type mass: float
        :type DeltaTemp: float
        :returns: Maximum altitude in meters [m].
        :rtype: float This function calculates the maximum possible altitude
            based on the aircraft's mass, ISA temperature deviation, and
            engine-specific parameters such as temperature and mass gradients.
            It considers the maximum operational altitude and adjusts for the
            given conditions.
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
            hMax = min(
                self.AC.hmo, self.AC.Hmax + Gt * var + Gw * (mMax - mass)
            )

        return conv.ft2m(hMax)

    def VStall(self, mass, config):
        """Computes the stall speed based on the aircraft configuration and
        mass.

        :param config: Aircraft configuration (e.g., 'CR', 'TO', 'AP', 'LD').
        :param mass: Aircraft mass in kilograms [kg].
        :type config: str
        :type mass: float
        :returns: Stall speed in meters per second [m/s].
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

    def VMin(
        self, h, mass, config, DeltaTemp, nz=1.2, envelopeType="OPERATIONAL"
    ):
        """Computes the minimum speed for a given configuration and
        conditions.

        :param h: Altitude in meters [m].
        :param mass: Aircraft mass in kilograms [kg].
        :param config: Aircraft configuration (e.g., 'CR', 'TO', 'LD').
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param nz: Load factor, default is 1.2.
        :param envelopeType: Type of flight envelope ('OPERATIONAL' or
            'CERTIFIED').
        :type h: float
        :type mass: float
        :type config: str
        :type DeltaTemp: float
        :type nz: float, optional
        :type envelopeType: str, optional
        :returns: Minimum speed in meters per second [m/s].
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
            return VminStall

        if self.AC.Clbo == 0.0 and self.AC.k == 0.0:
            Vmin = VminStall
        else:
            if h < conv.ft2m(15000):
                Vmin = VminStall
            elif h >= conv.ft2m(15000):
                # low speed buffeting limit applies only for JET and TURBOPROP
                if (
                    self.AC.engineType == "JET"
                    or self.AC.engineType == "TURBOPROP"
                ):
                    [theta, delta, sigma] = atm.atmosphereProperties(
                        h=h, DeltaTemp=DeltaTemp
                    )
                    buffetLimit = self.lowSpeedBuffetLimit(
                        h=h, mass=mass, DeltaTemp=DeltaTemp, nz=nz
                    )
                    if isnan(buffetLimit):
                        Vmin = VminStall
                    else:
                        Vmin = max(
                            VminStall,
                            atm.mach2Cas(
                                buffetLimit,
                                theta=theta,
                                delta=delta,
                                sigma=sigma,
                            ),
                        )
                elif (
                    self.AC.engineType == "PISTON"
                    or self.AC.engineType == "ELECTRIC"
                ):
                    Vmin = VminStall

        return Vmin

    def Vmax_thrustLimited(self, h, mass, DeltaTemp, rating, config):
        """Computes the maximum CAS speed considering thrust limitations
        within the certified flight envelope.

        :param h: Altitude in meters [m].
        :param mass: Aircraft mass in kilograms [kg].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param rating: Aircraft engine rating (e.g., 'MTKF', 'MCMB', 'MCRZ').
        :param config: Aircraft configuration (e.g., 'TO', 'CR').
        :type h: float
        :type mass: float
        :type DeltaTemp: float
        :type rating: str
        :type config: str
        :returns: Maximum thrust-limited speed in meters per second [m/s].
        :rtype: float
        """

        [theta, delta, sigma] = atm.atmosphereProperties(
            h=h, DeltaTemp=DeltaTemp
        )

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
        for CAS in np.linspace(
            VminCertified, VmaxCertified, num=200, endpoint=True
        ):
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
        """Computes the best angle of climb (Vx) speed.

        :param h: Altitude in meters [m].
        :param mass: Aircraft mass in kilograms [kg].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param rating: Aircraft engine rating (e.g., 'MTKF', 'MCMB', 'MCRZ').
        :param config: Aircraft configuration (e.g., 'TO', 'CR').
        :type h: float
        :type mass: float
        :type DeltaTemp: float
        :type rating: str
        :type config: str
        :returns: Best angle of climb speed (Vx) in meters per second [m/s].
        :rtype: float
        """

        [theta, delta, sigma] = atm.atmosphereProperties(
            h=h, DeltaTemp=DeltaTemp
        )

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

        for CAS in np.linspace(
            VminCertified, VmaxCertified, num=200, endpoint=True
        ):
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
        """Computes the maximum speed based on altitude and temperature
        deviation.

        :param h: Altitude in meters [m].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :type h: float
        :type DeltaTemp: float
        :returns: Maximum speed in meters per second [m/s].
        :rtype: float
        """

        crossoverAlt = atm.crossOver(
            cas=conv.kt2ms(self.AC.VMO), Mach=self.AC.MMO
        )

        if h >= crossoverAlt:
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=h, DeltaTemp=DeltaTemp
            )
            VMax = atm.mach2Cas(
                Mach=self.AC.MMO, theta=theta, delta=delta, sigma=sigma
            )
        else:
            VMax = conv.kt2ms(self.AC.VMO)

        return VMax

    def lowSpeedBuffetLimit(self, h, mass, DeltaTemp, nz=1.2):
        """Computes the low-speed buffet limit using numerical methods.

        :param h: Altitude in meters [m].
        :param mass: Aircraft mass in kilograms [kg].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param nz: Load factor, default is 1.2.
        :type h: float
        :type mass: float
        :type DeltaTemp: float
        :type nz: float, optional
        :returns: Low-speed buffet limit as Mach number [-].
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
        """Returns the aerodynamic configuration based on altitude, speed, and
        phase of flight.

        :param phase: Phase of flight (e.g., 'Climb', 'Cruise', 'Descent').
        :param h: Altitude in meters [m].
        :param v: Calibrated airspeed in meters per second [m/s].
        :param mass: Aircraft mass in kilograms [kg].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param hRWY: Runway elevation above mean sea level in meters [m],
            default is 0.
        :param nz: Load factor, default is 1.2.
        :type phase: str
        :type h: float
        :type v: float
        :type mass: float
        :type DeltaTemp: float
        :type hRWY: float, optional
        :type nz: float, optional
        :returns: Aerodynamic configuration (e.g., 'TO', 'IC', 'CR', 'AP',
            'LD').
        :rtype: str
        """

        config = None

        # aircraft AGL altitude assuming being close to the RWY [m]
        h_AGL = h - hRWY

        HmaxTO_AGL = (
            conv.ft2m(
                Parser.getGPFValue(self.AC.GPFdata, "H_max_to", phase="to")
            )
            - hRWY
        )
        HmaxIC_AGL = (
            conv.ft2m(
                Parser.getGPFValue(self.AC.GPFdata, "H_max_ic", phase="ic")
            )
            - hRWY
        )
        HmaxAPP_AGL = (
            conv.ft2m(
                Parser.getGPFValue(self.AC.GPFdata, "H_max_app", phase="app")
            )
            - hRWY
        )
        HmaxLD_AGL = (
            conv.ft2m(
                Parser.getGPFValue(self.AC.GPFdata, "H_max_ld", phase="lnd")
            )
            - hRWY
        )

        if phase == "Climb" and h_AGL <= HmaxTO_AGL:
            config = "TO"
        elif phase == "Climb" and (h_AGL > HmaxTO_AGL and h_AGL <= HmaxIC_AGL):
            config = "IC"
        else:
            vMinCR = self.VMin(
                h=h, mass=mass, config="CR", DeltaTemp=DeltaTemp, nz=nz
            )
            vMinAPP = self.VMin(
                h=h, mass=mass, config="AP", DeltaTemp=DeltaTemp, nz=nz
            )
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
        """Returns the aerodynamic configuration ID for a given configuration.

        :param config: Aircraft configuration (e.g., 'CR', 'IC', 'TO', 'AP',
            'LD').
        :type config: str
        :returns: A list containing the HLid and LG for the given
            configuration.
        :rtype: [int, str]
        """

        HLid = self.AC.aeroConfig[config]["HLid"]
        LG = self.AC.aeroConfig[config]["LG"]

        return [HLid, LG]

    def getHoldSpeed(self, h, theta, delta, sigma, DeltaTemp):
        """Computes the aircraft's holding speed (CAS) based on the current
        altitude.

        :param h: Altitude in meters [m].
        :param theta: Normalized temperature [-].
        :param delta: Normalized pressure [-].
        :param sigma: Normalized air density [-].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :type h: float
        :type theta: float
        :type delta: float
        :type sigma: float
        :type DeltaTemp: float
        :returns: Holding calibrated airspeed (CAS) in meters per second
            [m/s].
        :rtype: float
        """

        if h <= conv.ft2m(14000):
            vHold = Parser.getGPFValue(
                self.AC.GPFdata, "V_hold_1", phase="hold"
            )
        elif h > conv.ft2m(14000) and h <= conv.ft2m(20000):
            vHold = Parser.getGPFValue(
                self.AC.GPFdata, "V_hold_2", phase="hold"
            )
        elif h > conv.ft2m(20000) and h <= conv.ft2m(34000):
            vHold = Parser.getGPFValue(
                self.AC.GPFdata, "V_hold_3", phase="hold"
            )
        elif h > conv.ft2m(34000):
            MHold = Parser.getGPFValue(
                self.AC.GPFdata, "V_hold_4", phase="hold"
            )
            vHold = atm.mach2Cas(Mach=M, theta=theta, delta=delta, sigma=sigma)

        return conv.kt2ms(vHold)

    def getGroundMovementSpeed(self, pos):
        """Returns the ground movement speed based on the aircraft's position
        on the ground.

        :param pos: Aircraft position on the airport ground (e.g.,
            'backtrack', 'taxi', 'apron', 'gate').
        :type pos: str
        :returns: Ground movement speed in meters per second [m/s].
        :rtype: float
        """

        if pos == "backtrack":
            vGround = Parser.getGPFValue(
                self.AC.GPFdata, "V_backtrack", phase="gnd"
            )
        elif pos == "taxi":
            vGround = Parser.getGPFValue(
                self.AC.GPFdata, "V_taxi", phase="gnd"
            )
        elif pos == "apron":
            vGround = Parser.getGPFValue(
                self.AC.GPFdata, "V_apron", phase="gnd"
            )
        elif pos == "gate":
            vGround = Parser.getGPFValue(
                self.AC.GPFdata, "V_gate", phase="gnd"
            )

        return conv.kt2ms(vGround)

    def getBankAngle(self, phase, flightUnit, value):
        """Returns the nominal or maximum bank angle for the given flight
        phase and unit type.

        :param phase: Phase of flight (e.g., 'to', 'ic', 'cl', 'cr', etc.).
        :param flightUnit: Flight unit (e.g., 'civ' for civilian, 'mil' for
            military).
        :param value: Desired value, either 'nom' for nominal or 'max' for
            maximum bank angle.
        :type phase: str
        :type flightUnit: str
        :type value: str
        :returns: Bank angle in degrees [deg].
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
        """Checks whether the acceleration between two time steps is within
        allowable limits.

        :param v1: Airspeed (or vertical speed for 'norm') at the previous
            time step [m/s].
        :param v2: Airspeed (or vertical speed for 'norm') at the current time
            step [m/s].
        :param type: Type of acceleration to check ('long' for longitudinal,
            'norm' for normal).
        :param flightUnit: Flight unit type ('civ' for civilian, 'mil' for
            military).
        :param deltaTime: Time difference between the two time steps in
            seconds [s].
        :type v1: float
        :type v2: float
        :type type: str
        :type flightUnit: str
        :type deltaTime: float
        :returns: True if the acceleration is within limits, False otherwise.
        :rtype: bool
        """

        OK = False

        if flightUnit == "civ":
            if type == "long":
                if (
                    abs(v2 - v1)
                    <= conv.ft2m(
                        Parser.getGPFValue(self.AC.GPFdata, "acc_long_max")
                    )
                    * deltaTime
                ):
                    OK = True

        elif type == "norm":
            if (
                abs(v2 - v1)
                <= conv.ft2m(
                    Parser.getGPFValue(self.AC.GPFdata, "acc_norm_max")
                )
                * deltaTime
            ):
                OK = True

        # currently undefined for BADA3
        elif flightUnit == "mil":
            OK = True

        return OK

    def getSpeedSchedule(self, phase):
        """Returns the speed schedule for a given phase of flight.

        :param phase: Flight phase ('Climb', 'Cruise', 'Descent').
        :type phase: str
        :returns: A list containing CAS1, CAS2, and Mach number for the
            specified phase [m/s, m/s, -].
        :rtype: list[float, float, float]
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

    def checkConfigurationContinuity(
        self, phase, previousConfig, currentConfig
    ):
        """Ensures the continuity of aerodynamic configuration changes based
        on the phase of flight.

        :param phase: Current flight phase ('Climb', 'Cruise', 'Descent').
        :param previousConfig: The previous aerodynamic configuration.
        :param currentConfig: The current aerodynamic configuration.
        :type phase: str
        :type previousConfig: str
        :type currentConfig: str
        :returns: Updated aerodynamic configuration.
        :rtype: str This function ensures that the aerodynamic configuration
            transitions logically based on the phase of flight. For example,
            during descent, the configuration should not revert to a clean
            configuration after deploying flaps for approach or landing.
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
    """This class is a BADA3 aircraft subclass and implements the Airline
    Procedure Model (ARPM) following the BADA3 user manual.

    :param AC: Aircraft object {BADA3}.
    :type AC: bada3Aircraft.
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
        """Computes the climb speed schedule (CAS) for the given altitude
        based on various procedures and aircraft parameters.

        :param theta: Normalized air temperature [-].
        :param delta: Normalized air pressure [-].
        :param mass: Aircraft mass in kilograms [kg].
        :param h: Altitude in meters [m].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param speedSchedule_default: Optional, a default speed schedule that overrides the BADA schedule. It should be in the form [Vcl1, Vcl2, Mcl].
        :param applyLimits: Boolean flag indicating whether to apply the minimum and maximum speed limits based on the flight envelope.
        :param config: Optional, current aircraft aerodynamic configuration (TO/IC/CR/AP/LD).
        :param procedure: Climb procedure to be followed, e.g., 'BADA', 'NADP1', 'NADP2'. Default is 'BADA'.
        :param NADP1_ALT: Altitude in feet for NADP1 procedure. Default is 3000 feet.
        :param NADP2_ALT: Altitude range in feet for NADP2 procedure. Default is [1000, 3000].
        :type theta: float
        :type delta: float
        :type mass: float
        :type h: float
        :type DeltaTemp: float
        :type speedSchedule_default: list[float, float, float], optional
        :type applyLimits: bool
        :type config: str, optional
        :type procedure: str
        :type NADP1_ALT: float
        :type NADP2_ALT: list[float, float]
        :returns: A tuple containing the climb calibrated airspeed (CAS) in meters per second [m/s] and a status flag indicating whether the calculated CAS is constrained ('C'), unconstrained ('V' or 'v'), or not altered ('').
        :rtype: tuple[float, str]

        This function computes the climb speed schedule for different phases of flight and aircraft types.
        It supports BADA, NADP1, and NADP2 procedures for both jet and turboprop/piston/electric aircraft.

        The climb schedule uses specific speed profiles depending on altitude and aircraft model. For jet engines, the speed is constrained
        below 250 knots below 10,000 feet, and then it follows a defined speed schedule, either from BADA or NADP procedures.

        Additionally, the function applies speed limits based on the aircraft's flight envelope, adjusting the calculated climb speed if necessary.

        - For `procedure='BADA'`, it uses the BADA climb speed schedule.
        - For `procedure='NADP1'`, it implements the Noise Abatement Departure Procedure 1.
        - For `procedure='NADP2'`, it implements the Noise Abatement Departure Procedure 2.

        The function also ensures that the calculated CAS remains within the bounds of the aircraft's minimum and maximum speeds.
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
                        Parser.getGPFValue(
                            self.AC.GPFdata, "V_cl_5", phase=phase
                        )
                    )
                )
                speed.append(
                    Cvmin * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(
                            self.AC.GPFdata, "V_cl_4", phase=phase
                        )
                    )
                )
                speed.append(
                    Cvmin * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(
                            self.AC.GPFdata, "V_cl_3", phase=phase
                        )
                    )
                )
                speed.append(
                    Cvmin * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(
                            self.AC.GPFdata, "V_cl_2", phase=phase
                        )
                    )
                )
                speed.append(
                    Cvmin * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(
                            self.AC.GPFdata, "V_cl_1", phase=phase
                        )
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
                    cas = atm.mach2Cas(
                        Mach=Mcl, theta=theta, delta=delta, sigma=sigma
                    )

            elif (
                acModel == "TURBOPROP"
                or acModel == "PISTON"
                or acModel == "ELECTRIC"
            ):
                speed = list()
                speed.append(min(Vcl1, conv.kt2ms(250)))
                speed.append(
                    Cvmin * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(
                            self.AC.GPFdata,
                            "V_cl_8",
                            engine="TURBOPROP",
                            phase=phase,
                        )
                    )
                )
                speed.append(
                    Cvmin * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(
                            self.AC.GPFdata,
                            "V_cl_7",
                            engine="TURBOPROP",
                            phase=phase,
                        )
                    )
                )
                speed.append(
                    Cvmin * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(
                            self.AC.GPFdata,
                            "V_cl_6",
                            engine="TURBOPROP",
                            phase=phase,
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
                    cas = atm.mach2Cas(
                        Mach=Mcl, theta=theta, delta=delta, sigma=sigma
                    )

        elif procedure == "NADP1":
            if acModel == "JET":
                speed = list()
                speed.append(min(Vcl1, conv.kt2ms(250)))
                speed.append(
                    CvminTO * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(
                            self.AC.GPFdata, "V_cl_2", phase=phase
                        )
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
                    cas = atm.mach2Cas(
                        Mach=Mcl, theta=theta, delta=delta, sigma=sigma
                    )

            elif acModel == "TURBOPROP" or acModel == "PISTON":
                speed = list()
                speed.append(min(Vcl1, conv.kt2ms(250)))
                speed.append(
                    CvminTO * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(
                            self.AC.GPFdata, "V_cl_1", phase=phase
                        )
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
                    cas = atm.mach2Cas(
                        Mach=Mcl, theta=theta, delta=delta, sigma=sigma
                    )

        elif procedure == "NADP2":
            if acModel == "JET":
                speed = list()
                speed.append(min(Vcl1, conv.kt2ms(250)))
                speed.append(
                    Cvmin * VstallCR
                    + conv.kt2ms(
                        Parser.getGPFValue(
                            self.AC.GPFdata, "V_cl_2", phase=phase
                        )
                    )
                )
                speed.append(
                    CvminTO * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(
                            self.AC.GPFdata, "V_cl_2", phase=phase
                        )
                    )
                )

                n = 1
                while n < len(speed):
                    if speed[n] > speed[n - 1]:
                        speed[n] = speed[n - 1]
                    n = n + 1

                if h < conv.ft2m(NADP2_ALT[0]):
                    cas = speed[2]
                elif h >= conv.ft2m(NADP2_ALT[0]) and h < conv.ft2m(
                    NADP2_ALT[1]
                ):
                    cas = speed[1]
                elif h >= conv.ft2m(NADP2_ALT[1]) and h < conv.ft2m(10000):
                    cas = speed[0]
                elif h >= conv.ft2m(10000) and h < crossOverAlt:
                    cas = Vcl2
                elif h >= crossOverAlt:
                    sigma = atm.sigma(theta=theta, delta=delta)
                    cas = atm.mach2Cas(
                        Mach=Mcl, theta=theta, delta=delta, sigma=sigma
                    )

            elif acModel == "TURBOPROP" or acModel == "PISTON":
                speed = list()
                speed.append(min(Vcl1, conv.kt2ms(250)))
                speed.append(
                    Cvmin * VstallCR
                    + conv.kt2ms(
                        Parser.getGPFValue(
                            self.AC.GPFdata, "V_cl_2", phase=phase
                        )
                    )
                )
                speed.append(
                    CvminTO * VstallTO
                    + conv.kt2ms(
                        Parser.getGPFValue(
                            self.AC.GPFdata, "V_cl_1", phase=phase
                        )
                    )
                )

                n = 1
                while n < len(speed):
                    if speed[n] > speed[n - 1]:
                        speed[n] = speed[n - 1]
                    n = n + 1

                if h < conv.ft2m(NADP2_ALT[0]):
                    cas = speed[2]
                elif h >= conv.ft2m(NADP2_ALT[0]) and h < conv.ft2m(
                    NADP2_ALT[1]
                ):
                    cas = speed[1]
                elif h >= conv.ft2m(NADP2_ALT[1]) and h < conv.ft2m(10000):
                    cas = speed[0]
                elif h >= conv.ft2m(10000) and h < crossOverAlt:
                    cas = Vcl2
                elif h >= crossOverAlt:
                    sigma = atm.sigma(theta=theta, delta=delta)
                    cas = atm.mach2Cas(
                        Mach=Mcl, theta=theta, delta=delta, sigma=sigma
                    )

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
        """Computes the cruise speed schedule (CAS) for a given altitude based
        on aircraft parameters and procedures.

        :param h: Altitude in meters [m].
        :param mass: Aircraft mass in kilograms [kg].
        :param theta: Normalized air temperature [-].
        :param delta: Normalized air pressure [-].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param speedSchedule_default: Optional, a default speed schedule that overrides the BADA schedule. It should be in the form [Vcr1, Vcr2, Mcr].
        :param applyLimits: Boolean flag indicating whether to apply the minimum and maximum speed limits based on the flight envelope.
        :param config: Optional, current aircraft aerodynamic configuration (TO/IC/CR/AP/LD).
        :type h: float
        :type mass: float
        :type theta: float
        :type delta: float
        :type DeltaTemp: float
        :type speedSchedule_default: list[float, float, float], optional
        :type applyLimits: bool
        :type config: str, optional
        :returns: A tuple containing the cruise calibrated airspeed (CAS) in meters per second [m/s] and a status flag indicating whether the calculated CAS is constrained ('C'), unconstrained ('V' or 'v'), or not altered ('').
        :rtype: tuple[float, str]

        This function computes the cruise speed schedule for various phases of flight and aircraft models.
        It supports both jet and turboprop/piston/electric aircraft models by using the BADA (Base of Aircraft Data) speed schedules.

        - If a `speedSchedule_default` is provided, it overwrites the BADA speed schedule.
        - For jet engines, the speed is constrained based on altitude, starting with 170 knots below 3000 feet, 220 knots below 6000 feet, and then follows the standard speed schedule.
        - For other aircraft types (TURBOPROP, PISTON, ELECTRIC), the speed limits are lower, starting with 150 knots below 3000 feet.

        The function also applies limits based on the aircraft's flight envelope, ensuring the calculated speed does not exceed the minimum or maximum allowable speeds.
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
                cas = atm.mach2Cas(
                    Mach=Mcr, theta=theta, delta=delta, sigma=sigma
                )

        elif (
            acModel == "TURBOPROP"
            or acModel == "PISTON"
            or acModel == "ELECTRIC"
        ):
            if h < conv.ft2m(3000):
                cas = min(Vcr1, conv.kt2ms(150))
            elif h >= conv.ft2m(3000) and h < conv.ft2m(6000):
                cas = min(Vcr1, conv.kt2ms(180))
            elif h >= conv.ft2m(6000) and h < conv.ft2m(10000):
                cas = min(Vcr1, conv.kt2ms(250))
            elif h >= conv.ft2m(10000) and h < crossOverAlt:
                cas = Vcr2
            elif h >= crossOverAlt:
                cas = atm.mach2Cas(
                    Mach=Mcr, theta=theta, delta=delta, sigma=sigma
                )

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
        """Computes the descent speed schedule (CAS) for a given altitude
        based on aircraft parameters and procedures.

        :param h: Altitude in meters [m].
        :param mass: Aircraft mass in kilograms [kg].
        :param theta: Normalized air temperature [-].
        :param delta: Normalized air pressure [-].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param speedSchedule_default: Optional, a default speed schedule that overrides the BADA schedule. It should be in the form [Vdes1, Vdes2, Mdes].
        :param applyLimits: Boolean flag indicating whether to apply the minimum and maximum speed limits based on the flight envelope.
        :param config: Optional, current aircraft aerodynamic configuration (TO/IC/CR/AP/LD).
        :type h: float
        :type mass: float
        :type theta: float
        :type delta: float
        :type DeltaTemp: float
        :type speedSchedule_default: list[float, float, float], optional
        :type applyLimits: bool
        :type config: str, optional
        :returns: A tuple containing the descent calibrated airspeed (CAS) in meters per second [m/s] and a status flag indicating whether the calculated CAS is constrained ('C'), unconstrained ('V' or 'v'), or not altered ('').
        :rtype: tuple[float, str]

        This function computes the descent speed schedule for various phases of flight and aircraft models.
        It supports both jet and turboprop/piston/electric aircraft models using the BADA (Base of Aircraft Data) speed schedules.

        - If a `speedSchedule_default` is provided, it overwrites the BADA speed schedule.
        - For jet and turboprop engines, the speed schedule is constrained based on altitude, starting from 220 knots below 3000 feet and then following the standard speed schedule.
        - For piston and electric engines, lower speed limits are applied based on stall speeds.

        The function also applies limits based on the aircraft's flight envelope, ensuring that the calculated speed does not exceed the minimum or maximum allowable speeds.
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
                cas = atm.mach2Cas(
                    Mach=Mdes, theta=theta, delta=delta, sigma=sigma
                )

        elif acModel == "PISTON" or acModel == "ELECTRIC":
            speed = list()
            speed.append(Vdes1)
            speed.append(
                Cvmin * VstallDES
                + conv.kt2ms(
                    Parser.getGPFValue(
                        self.AC.GPFdata,
                        "V_des_7",
                        engine="PISTON",
                        phase=phase,
                    )
                )
            )
            speed.append(
                Cvmin * VstallDES
                + conv.kt2ms(
                    Parser.getGPFValue(
                        self.AC.GPFdata,
                        "V_des_6",
                        engine="PISTON",
                        phase=phase,
                    )
                )
            )
            speed.append(
                Cvmin * VstallDES
                + conv.kt2ms(
                    Parser.getGPFValue(
                        self.AC.GPFdata,
                        "V_des_5",
                        engine="PISTON",
                        phase=phase,
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
                cas = atm.mach2Cas(
                    Mach=Mdes, theta=theta, delta=delta, sigma=sigma
                )

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
    """This class implements the PTD file creator for BADA3 aircraft following
    BADA3 manual.

    :param AC: Aircraft object {BADA3}.
    :type AC: bada3Aircraft.
    """

    def __init__(self, AC):
        super().__init__(AC)

        self.flightEnvelope = FlightEnvelope(AC)
        self.ARPM = ARPM(AC)

    def create(self, DeltaTemp, saveToPath):
        """Creates a BADA3 PTD file based on specified temperature deviation
        from ISA and saves it to the provided directory path. It generates
        performance data for different aircraft mass levels (low, medium,
        high) in both climb and descent phases.

        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param saveToPath: Path to the directory where the PTD file will be
            stored.
        :type DeltaTemp: float.
        :type saveToPath: str.
        :returns: None
        :rtype: None
        """

        # 3 different mass levels [kg]
        if 1.2 * self.AC.mass["minimum"] > self.AC.mass["reference"]:
            massLow = self.AC.mass["minimum"]
        else:
            massLow = 1.2 * self.AC.mass["minimum"]

        massList = [
            massLow,
            self.AC.mass["reference"],
            self.AC.mass["maximum"],
        ]
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
        self,
        saveToPath,
        CLList_low,
        CLList_med,
        CLList_high,
        DESList_med,
        DeltaTemp,
    ):
        """Saves BADA3 (PTD) to a file. It stores performance data for low,
        medium, and high aircraft masses during the climb phase, and medium
        aircraft mass during the descent phase. The file is saved in a
        predefined format.

        :param saveToPath: Path to the directory where the PTD file should be
            saved.
        :param CLList_low: List containing PTD data for CLIMB at low aircraft
            mass.
        :param CLList_med: List containing PTD data for CLIMB at medium
            aircraft mass.
        :param CLList_high: List containing PTD data for CLIMB at high
            aircraft mass.
        :param DESList_med: List containing PTD data for DESCENT at medium
            aircraft mass.
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :type saveToPath: str.
        :type CLList_low: list.
        :type CLList_med: list.
        :type CLList_high: list.
        :type DESList_med: list.
        :type DeltaTemp: float.
        :returns: None
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
        file.write(
            "=============================\n=============================\n\n"
        )
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
        """Calculates the BADA3 PTD data in climb phase.

        :param mass: Aircraft mass [kg]
        :param altitudeList: List of altitude levels for calculation (in feet)
        :param DeltaTemp: Deviation from International Standard Atmosphere
            (ISA) temperature [K]
        :type mass: float
        :type altitudeList: list of int
        :type DeltaTemp: float
        :returns: A list of calculated PTD data for the climb phase
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
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=H_m, DeltaTemp=DeltaTemp
            )
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
        """Calculates the BADA3 PTD data in descent phase.

        This function generates a detailed list of descent performance metrics
        for different altitudes and mass configurations based on BADA3
        performance models.

        :param mass: Aircraft mass [kg].
        :param altitudeList: List of aircraft altitudes in feet [ft].
        :param DeltaTemp: Deviation from ISA temperature [K].
        :type mass: float.
        :type altitudeList: list of int.
        :type DeltaTemp: float.
        :returns: List of descent performance data.
        :rtype: list.
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
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=H_m, DeltaTemp=DeltaTemp
            )
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

            if (
                self.AC.engineType == "PISTON"
                or self.AC.engineType == "ELECTRIC"
            ):
                # PISTON  and ELECTRIC uses LIDL throughout the whole descent phase
                Thrust = self.Thrust(
                    rating="LIDL",
                    v=tas,
                    h=H_m,
                    config="CR",
                    DeltaTemp=DeltaTemp,
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
                    rating="LIDL",
                    v=tas,
                    h=H_m,
                    config=config,
                    DeltaTemp=DeltaTemp,
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

            tau_const = (theta * const.temp_0) / (
                theta * const.temp_0 - DeltaTemp
            )
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
    """This class implements the PTF file creator for BADA3 aircraft following
    BADA3 manual.

    :param AC: Aircraft object {BADA3}.
    :type AC: bada3Aircraft.
    """

    def __init__(self, AC):
        super().__init__(AC)

        self.flightEnvelope = FlightEnvelope(AC)
        self.ARPM = ARPM(AC)

    def create(self, DeltaTemp, saveToPath):
        """Creates a BADA3 PTF file based on specified temperature deviation
        from ISA and saves it to the provided directory path. It generates
        performance data for different aircraft mass levels (low, medium,
        high) in both climb and descent phases.

        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param saveToPath: Path to the directory where the PTF file will be
            stored.
        :type DeltaTemp: float.
        :type saveToPath: str.
        :returns: None
        :rtype: None
        """

        # 3 different mass levels [kg]
        if 1.2 * self.AC.mass["minimum"] > self.AC.mass["reference"]:
            massLow = self.AC.mass["minimum"]
        else:
            massLow = 1.2 * self.AC.mass["minimum"]

        massList = [
            massLow,
            self.AC.mass["reference"],
            self.AC.mass["maximum"],
        ]
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
        self,
        saveToPath,
        CRList,
        CLList,
        DESList,
        altitudeList,
        massList,
        DeltaTemp,
    ):
        """Saves performance data to a PTF file.

        :param saveToPath: Directory path where the PTF file will be stored.
        :param CRList: List of cruise phase data.
        :param CLList: List of climb phase data.
        :param DESList: List of descent phase data.
        :param altitudeList: List of aircraft altitudes [ft].
        :param massList: List of aircraft masses [kg].
        :param DeltaTemp: Deviation from ISA temperature [K].
        :type saveToPath: string.
        :type CRList: list.
        :type CLList: list.
        :type DESList: list.
        :type altitudeList: list of int.
        :type massList: list of int.
        :type DeltaTemp: float.
        :returns: None
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
            "BADA PERFORMANCE FILE                                        %s\n\n"
            % (d3)
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
        """Calculates BADA3 PTF data for the cruise phase.

        :param massList: List of aircraft masses [kg] (low, nominal, and
            high).
        :param altitudeList: List of aircraft altitudes [ft].
        :param DeltaTemp: Deviation from the International Standard Atmosphere
            (ISA) temperature [K].
        :type massList: list of float.
        :type altitudeList: list of int.
        :type DeltaTemp: float.
        :returns: List containing cruise phase TAS and fuel flow data.
        :rtype: list.
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
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=H_m, DeltaTemp=DeltaTemp
            )
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
                ff.append(
                    self.ff(flightPhase="Cruise", v=tas, h=H_m, T=Thrust) * 60
                )

            TAS_CR_complet.append(conv.ms2kt(tas_nominal))
            FF_CR_LO_complet.append(ff[0])
            FF_CR_NOM_complet.append(ff[1])
            FF_CR_HI_complet.append(ff[2])

        CRList = [
            TAS_CR_complet,
            FF_CR_LO_complet,
            FF_CR_NOM_complet,
            FF_CR_HI_complet,
        ]

        return CRList

    def PTF_climb(self, massList, altitudeList, DeltaTemp):
        """Calculates BADA3 PTF data for the climb phase.

        :param massList: List of aircraft masses [kg] (low, nominal, high).
        :param altitudeList: List of aircraft altitudes [ft].
        :param DeltaTemp: Deviation from the International Standard Atmosphere
            (ISA) temperature [K].
        :type massList: list of float.
        :type altitudeList: list of int.
        :type DeltaTemp: float.
        :returns: List containing climb phase TAS, ROCD, and fuel flow data.
        :rtype: list.
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
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=H_m, DeltaTemp=DeltaTemp
            )
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
                    h=H_m,
                    phase="Climb",
                    v=cas,
                    mass=massNominal,
                    DeltaTemp=DeltaTemp,
                )
                CD = self.CD(CL=CL, config=config)
                Drag = self.D(tas=tas, sigma=sigma, CD=CD)
                Thrust = self.Thrust(
                    rating="MCMB",
                    v=tas,
                    h=H_m,
                    config=config,
                    DeltaTemp=DeltaTemp,
                )
                ff = self.ff(flightPhase="Climb", v=tas, h=H_m, T=Thrust) * 60

                if H_m < crossAlt:
                    ESF = self.esf(
                        h=H_m,
                        flightEvolution="constCAS",
                        M=M,
                        DeltaTemp=DeltaTemp,
                    )
                else:
                    ESF = self.esf(
                        h=H_m,
                        flightEvolution="constM",
                        M=M,
                        DeltaTemp=DeltaTemp,
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
        """Calculates BADA3 PTF data for the descent phase.

        :param massList: List of aircraft masses [kg] (low, nominal, high).
        :param altitudeList: List of aircraft altitudes [ft].
        :param DeltaTemp: Deviation from the International Standard Atmosphere
            (ISA) temperature [K].
        :type massList: list of float.
        :type altitudeList: list of int.
        :type DeltaTemp: float.
        :returns: List containing descent phase TAS, ROCD, and fuel flow data.
        :rtype: list.
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
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=H_m, DeltaTemp=DeltaTemp
            )
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
                h=H_m,
                phase="Descent",
                v=cas,
                mass=massNominal,
                DeltaTemp=DeltaTemp,
            )

            CL = self.CL(tas=tas_nominal, sigma=sigma, mass=massNominal)
            CD = self.CD(CL=CL, config=config)
            Drag = self.D(tas=tas_nominal, sigma=sigma, CD=CD)

            if (
                self.AC.engineType == "PISTON"
                or self.AC.engineType == "ELECTRIC"
            ):
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
    """Implements the BADA3 performance model for an aircraft following the
    BADA3 manual.

    This class handles the loading of aircraft-specific data from either a
    predefined dataset or a set of BADA3 performance model files (e.g., OPF
    and APF files). It initializes various parameters such as mass, speed
    schedules, and engine type necessary for simulating the aircraft's
    performance.

    :param badaVersion: The BADA version being used.
    :param acName: The ICAO aircraft designation (e.g., "A320").
    :param filePath: Optional path to the BADA3 formatted file. If not
        provided, the default aircraft directory is used.
    :param allData: Optional DataFrame containing all aircraft data. If
        provided, the class will try to load the aircraft data from this
        DataFrame.
    :type badaVersion: str.
    :type acName: str.
    :type filePath: str, optional.
    :type allData: pd.DataFrame, optional.
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

        if filePath is None:
            self.filePath = configuration.getBadaVersionPath(
                badaFamily="BADA3", badaVersion=badaVersion
            )
        else:
            self.filePath = filePath

        # check if the aircraft is in the allData dataframe data
        if allData is not None and acName in allData["acName"].values:
            filtered_df = allData[allData["acName"] == acName]

            self.acName = configuration.safe_get(filtered_df, "acName", None)
            self.xmlFiles = configuration.safe_get(
                filtered_df, "xmlFiles", None
            )

            self.modificationDateOPF = configuration.safe_get(
                filtered_df, "modificationDateOPF", None
            )
            self.modificationDateAPF = configuration.safe_get(
                filtered_df, "modificationDateAPF", None
            )

            self.ICAO = configuration.safe_get(filtered_df, "ICAO", None)
            self.numberOfEngines = configuration.safe_get(
                filtered_df, "numberOfEngines", None
            )
            self.engineType = configuration.safe_get(
                filtered_df, "engineType", None
            )
            self.engines = configuration.safe_get(filtered_df, "engines", None)
            self.WTC = configuration.safe_get(filtered_df, "WTC", None)
            self.mass = configuration.safe_get(filtered_df, "mass", None)

            self.MTOW = configuration.safe_get(filtered_df, "MTOW", None)
            self.OEW = configuration.safe_get(filtered_df, "OEW", None)
            self.MPL = configuration.safe_get(filtered_df, "MPL", None)
            self.MREF = configuration.safe_get(filtered_df, "MREF", None)
            self.VMO = configuration.safe_get(filtered_df, "VMO", None)
            self.MMO = configuration.safe_get(filtered_df, "MMO", None)
            self.hmo = configuration.safe_get(filtered_df, "hmo", None)
            self.Hmax = configuration.safe_get(filtered_df, "Hmax", None)
            self.tempGrad = configuration.safe_get(
                filtered_df, "tempGrad", None
            )

            self.S = configuration.safe_get(filtered_df, "S", None)
            self.Clbo = configuration.safe_get(filtered_df, "Clbo", None)
            self.k = configuration.safe_get(filtered_df, "k", None)
            self.Vstall = configuration.safe_get(filtered_df, "Vstall", None)
            self.CD0 = configuration.safe_get(filtered_df, "CD0", None)
            self.CD2 = configuration.safe_get(filtered_df, "CD2", None)
            self.HLids = configuration.safe_get(filtered_df, "HLids", None)
            self.Ct = configuration.safe_get(filtered_df, "Ct", None)
            self.CTdeslow = configuration.safe_get(
                filtered_df, "CTdeslow", None
            )
            self.CTdeshigh = configuration.safe_get(
                filtered_df, "CTdeshigh", None
            )
            self.CTdesapp = configuration.safe_get(
                filtered_df, "CTdesapp", None
            )
            self.CTdesld = configuration.safe_get(filtered_df, "CTdesld", None)
            self.HpDes = configuration.safe_get(filtered_df, "HpDes", None)
            self.Cf = configuration.safe_get(filtered_df, "Cf", None)
            self.CfDes = configuration.safe_get(filtered_df, "CfDes", None)
            self.CfCrz = configuration.safe_get(filtered_df, "CfCrz", None)
            self.TOL = configuration.safe_get(filtered_df, "TOL", None)
            self.LDL = configuration.safe_get(filtered_df, "LDL", None)
            self.span = configuration.safe_get(filtered_df, "span", None)
            self.length = configuration.safe_get(filtered_df, "length", None)

            self.V1 = configuration.safe_get(filtered_df, "V1", None)
            self.V2 = configuration.safe_get(filtered_df, "V2", None)
            self.M = configuration.safe_get(filtered_df, "M", None)

            self.GPFdata = configuration.safe_get(filtered_df, "GPFdata", None)

            self.drone = configuration.safe_get(filtered_df, "drone", None)

            self.DeltaCD = configuration.safe_get(filtered_df, "DeltaCD", None)
            self.speedSchedule = configuration.safe_get(
                filtered_df, "speedSchedule", None
            )
            self.aeroConfig = configuration.safe_get(
                filtered_df, "aeroConfig", None
            )

            self.flightEnvelope = FlightEnvelope(self)
            self.ARPM = ARPM(self)
            self.PTD = PTD(self)
            self.PTF = PTF(self)

        else:
            # read BADA3 GPF file
            GPFDataFrame = Parser.parseGPF(self.filePath)

            # check if SYNONYM file exist
            synonymFile = os.path.join(self.filePath, "SYNONYM.NEW")
            synonymFileXML = os.path.join(self.filePath, "SYNONYM.xml")

            if os.path.isfile(synonymFile) or os.path.isfile(synonymFileXML):
                self.synonymFileAvailable = True

                self.SearchedACName = Parser.parseSynonym(
                    self.filePath, acName
                )

                if self.SearchedACName is None:
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
                        self.filePath,
                        self.SearchedACName,
                    )
                    + ".OPF"
                )
                APFfile = (
                    os.path.join(
                        self.filePath,
                        self.SearchedACName,
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
                        self.filePath, self.SearchedACName
                    )
                    APFDataFrame = Parser.parseAPF(
                        self.filePath, self.SearchedACName
                    )

                    OPF_APF_combined_df = Parser.combineOPF_APF(
                        OPFDataFrame, APFDataFrame
                    )
                    combined_df = Parser.combineACDATA_GPF(
                        OPF_APF_combined_df, GPFDataFrame
                    )

                    self.acName = configuration.safe_get(
                        combined_df, "acName", None
                    )
                    self.xmlFiles = configuration.safe_get(
                        combined_df, "xmlFiles", None
                    )

                    self.modificationDateOPF = configuration.safe_get(
                        combined_df, "modificationDateOPF", None
                    )
                    self.modificationDateAPF = configuration.safe_get(
                        combined_df, "modificationDateAPF", None
                    )

                    self.ICAO = configuration.safe_get(
                        combined_df, "ICAO", None
                    )
                    self.numberOfEngines = configuration.safe_get(
                        combined_df, "numberOfEngines", None
                    )
                    self.engineType = configuration.safe_get(
                        combined_df, "engineType", None
                    )
                    self.engines = configuration.safe_get(
                        combined_df, "engines", None
                    )
                    self.WTC = configuration.safe_get(combined_df, "WTC", None)
                    self.mass = configuration.safe_get(
                        combined_df, "mass", None
                    )

                    self.MTOW = configuration.safe_get(
                        combined_df, "MTOW", None
                    )
                    self.OEW = configuration.safe_get(combined_df, "OEW", None)
                    self.MPL = configuration.safe_get(combined_df, "MPL", None)
                    self.MREF = configuration.safe_get(
                        combined_df, "MREF", None
                    )
                    self.VMO = configuration.safe_get(combined_df, "VMO", None)
                    self.MMO = configuration.safe_get(combined_df, "MMO", None)
                    self.hmo = configuration.safe_get(combined_df, "hmo", None)
                    self.Hmax = configuration.safe_get(
                        combined_df, "Hmax", None
                    )
                    self.tempGrad = configuration.safe_get(
                        combined_df, "tempGrad", None
                    )

                    self.S = configuration.safe_get(combined_df, "S", None)
                    self.Clbo = configuration.safe_get(
                        combined_df, "Clbo", None
                    )
                    self.k = configuration.safe_get(combined_df, "k", None)
                    self.Vstall = configuration.safe_get(
                        combined_df, "Vstall", None
                    )
                    self.CD0 = configuration.safe_get(combined_df, "CD0", None)
                    self.CD2 = configuration.safe_get(combined_df, "CD2", None)
                    self.HLids = configuration.safe_get(
                        combined_df, "HLids", None
                    )
                    self.Ct = configuration.safe_get(combined_df, "Ct", None)
                    self.CTdeslow = configuration.safe_get(
                        combined_df, "CTdeslow", None
                    )
                    self.CTdeshigh = configuration.safe_get(
                        combined_df, "CTdeshigh", None
                    )
                    self.CTdesapp = configuration.safe_get(
                        combined_df, "CTdesapp", None
                    )
                    self.CTdesld = configuration.safe_get(
                        combined_df, "CTdesld", None
                    )
                    self.HpDes = configuration.safe_get(
                        combined_df, "HpDes", None
                    )
                    self.Cf = configuration.safe_get(combined_df, "Cf", None)
                    self.CfDes = configuration.safe_get(
                        combined_df, "CfDes", None
                    )
                    self.CfCrz = configuration.safe_get(
                        combined_df, "CfCrz", None
                    )
                    self.TOL = configuration.safe_get(combined_df, "TOL", None)
                    self.LDL = configuration.safe_get(combined_df, "LDL", None)
                    self.span = configuration.safe_get(
                        combined_df, "span", None
                    )
                    self.length = configuration.safe_get(
                        combined_df, "length", None
                    )

                    self.V1 = configuration.safe_get(combined_df, "V1", None)
                    self.V2 = configuration.safe_get(combined_df, "V2", None)
                    self.M = configuration.safe_get(combined_df, "M", None)

                    self.GPFdata = configuration.safe_get(
                        combined_df, "GPFdata", None
                    )

                    self.drone = configuration.safe_get(
                        combined_df, "drone", None
                    )

                    self.DeltaCD = configuration.safe_get(
                        combined_df, "DeltaCD", None
                    )
                    self.speedSchedule = configuration.safe_get(
                        combined_df, "speedSchedule", None
                    )
                    self.aeroConfig = configuration.safe_get(
                        combined_df, "aeroConfig", None
                    )

                    self.flightEnvelope = FlightEnvelope(self)
                    self.ARPM = ARPM(self)
                    self.PTD = PTD(self)
                    self.PTF = PTF(self)

                elif not self.OPFavailable and not self.APFavailable:
                    # search for xml files

                    XMLDataFrame = Parser.parseXML(
                        self.filePath, self.SearchedACName
                    )

                    combined_df = Parser.combineACDATA_GPF(
                        XMLDataFrame, GPFDataFrame
                    )

                    self.acName = configuration.safe_get(
                        combined_df, "acName", None
                    )
                    self.xmlFiles = configuration.safe_get(
                        combined_df, "xmlFiles", None
                    )

                    self.modificationDateOPF = configuration.safe_get(
                        combined_df, "modificationDateOPF", None
                    )
                    self.modificationDateAPF = configuration.safe_get(
                        combined_df, "modificationDateAPF", None
                    )

                    self.ICAO = configuration.safe_get(
                        combined_df, "ICAO", None
                    )
                    self.numberOfEngines = configuration.safe_get(
                        combined_df, "numberOfEngines", None
                    )
                    self.engineType = configuration.safe_get(
                        combined_df, "engineType", None
                    )
                    self.engines = configuration.safe_get(
                        combined_df, "engines", None
                    )
                    self.WTC = configuration.safe_get(combined_df, "WTC", None)
                    self.mass = configuration.safe_get(
                        combined_df, "mass", None
                    )

                    self.MTOW = configuration.safe_get(
                        combined_df, "MTOW", None
                    )
                    self.OEW = configuration.safe_get(combined_df, "OEW", None)
                    self.MPL = configuration.safe_get(combined_df, "MPL", None)
                    self.MREF = configuration.safe_get(
                        combined_df, "MREF", None
                    )
                    self.VMO = configuration.safe_get(combined_df, "VMO", None)
                    self.MMO = configuration.safe_get(combined_df, "MMO", None)
                    self.hmo = configuration.safe_get(combined_df, "hmo", None)
                    self.Hmax = configuration.safe_get(
                        combined_df, "Hmax", None
                    )
                    self.tempGrad = configuration.safe_get(
                        combined_df, "tempGrad", None
                    )

                    self.S = configuration.safe_get(combined_df, "S", None)
                    self.Clbo = configuration.safe_get(
                        combined_df, "Clbo", None
                    )
                    self.k = configuration.safe_get(combined_df, "k", None)
                    self.Vstall = configuration.safe_get(
                        combined_df, "Vstall", None
                    )
                    self.CD0 = configuration.safe_get(combined_df, "CD0", None)
                    self.CD2 = configuration.safe_get(combined_df, "CD2", None)
                    self.HLids = configuration.safe_get(
                        combined_df, "HLids", None
                    )
                    self.Ct = configuration.safe_get(combined_df, "Ct", None)
                    self.CTdeslow = configuration.safe_get(
                        combined_df, "CTdeslow", None
                    )
                    self.CTdeshigh = configuration.safe_get(
                        combined_df, "CTdeshigh", None
                    )
                    self.CTdesapp = configuration.safe_get(
                        combined_df, "CTdesapp", None
                    )
                    self.CTdesld = configuration.safe_get(
                        combined_df, "CTdesld", None
                    )
                    self.HpDes = configuration.safe_get(
                        combined_df, "HpDes", None
                    )
                    self.Cf = configuration.safe_get(combined_df, "Cf", None)
                    self.CfDes = configuration.safe_get(
                        combined_df, "CfDes", None
                    )
                    self.CfCrz = configuration.safe_get(
                        combined_df, "CfCrz", None
                    )
                    self.TOL = configuration.safe_get(combined_df, "TOL", None)
                    self.LDL = configuration.safe_get(combined_df, "LDL", None)
                    self.span = configuration.safe_get(
                        combined_df, "span", None
                    )
                    self.length = configuration.safe_get(
                        combined_df, "length", None
                    )

                    self.V1 = configuration.safe_get(combined_df, "V1", None)
                    self.V2 = configuration.safe_get(combined_df, "V2", None)
                    self.M = configuration.safe_get(combined_df, "M", None)

                    self.GPFdata = configuration.safe_get(
                        combined_df, "GPFdata", None
                    )

                    self.drone = configuration.safe_get(
                        combined_df, "drone", None
                    )

                    self.DeltaCD = configuration.safe_get(
                        combined_df, "DeltaCD", None
                    )
                    self.speedSchedule = configuration.safe_get(
                        combined_df, "speedSchedule", None
                    )
                    self.aeroConfig = configuration.safe_get(
                        combined_df, "aeroConfig", None
                    )

                    self.flightEnvelope = FlightEnvelope(self)
                    self.ARPM = ARPM(self)
                    self.PTD = PTD(self)
                    self.PTF = PTF(self)

                else:
                    # AC name cannot be found
                    raise ValueError(acName + " Cannot be found")

    def __str__(self):
        return f"(BADA3, AC_name: {self.acName}, searched_AC_name: {self.SearchedACName}, model_ICAO: {self.ICAO}, ID: {id(self.AC)})"
