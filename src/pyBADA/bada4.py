# -*- coding: utf-8 -*-
"""
pyBADA
Generic BADA4 aircraft performance module
Developed @EUROCONTROL (EIH)
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


from math import sqrt, pow, pi, isnan, sin, asin

import numpy as np
import os
from datetime import date
import xml.etree.ElementTree as ET

from scipy.optimize import fmin, fminbound
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


class Parser:
    """This class implements the BADA4 parsing mechanism to parse xml and GPF(xml) BADA4 files."""

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
    def readMappingFile(filePath, badaVersion):
        """This function parses BADA4 mapping xml formatted file and stores
        a dictionary of code names and files to be used for that specific code name

        :param filePath: path to the BADA4 mapping xml formatted file.
        :type filePath: str.
        :raises: IOError
        """

        filename = os.path.join(
            filePath, "BADA4", badaVersion, "aircraft_model_default.xml"
        )

        code_fileName = {}

        if os.path.isfile(filename):
            try:
                tree = ET.parse(filename)
                root = tree.getroot()
            except:
                raise IOError(filename + " not found or in correct format")

            # synonym - file name pair dictionary
            code_fileName = {}

            for child in root.iter("MAP"):
                code = child.find("code").text
                file = child.find("file").text

                code_fileName[code] = file

        return code_fileName

    @staticmethod
    def parseMappingFile(filePath, badaVersion, acName):

        code_fileName = Parser.readMappingFile(filePath, badaVersion)
        if acName in code_fileName:
            fileName = code_fileName[acName]
            return fileName
        else:
            return None

    @staticmethod
    def parseXML(filePath, badaVersion, acName):
        """This function parses BADA4 xml formatted file

        :param filePath: path to the folder with BADA4 xml formatted file.
        :param acName: name of Aircraft for BADA4 xml formatted file.
        :type filePath: str.
        :type acName: str.
        :raises: IOError
        """

        acXmlFile = (
            os.path.join(filePath, "BADA4", badaVersion, acName, acName) + ".xml"
        )

        try:
            tree = ET.parse(acXmlFile)
            root = tree.getroot()
        except:
            raise IOError(acXmlFile + " not found or in correct format")

        # Parse general aircraft data
        model = root.find("model").text  # aircraft model
        engineType = root.find("type").text  # aircraft type
        engines = root.find("engine").text  # engine type

        ICAO_desig = {}  # ICAO designator and WTC
        ICAO = root.find("ICAO").find("designator").text
        WTC = root.find("ICAO").find("WTC").text

        # Parse engine data
        PFM = root.find("PFM")  # get PFM

        MREF = float(PFM.find("MREF").text)  # reference mass
        WREF = MREF * const.g
        LHV = float(PFM.find("LHV").text)
        n_eng = int(PFM.find("n_eng").text)  # number of engines

        rho = float(PFM.find("rho").text)

        TFA = None
        if PFM.find("TFA") is not None:
            TFA = float(PFM.find("TFA").text)

        # parameters introduced with BADA 4.3
        p_delta = None
        if PFM.find("p_delta") is not None:
            p_delta = float(PFM.find("p_delta").text)

        p_theta = None
        if PFM.find("p_theta") is not None:
            p_theta = float(PFM.find("p_theta").text)

        TFM = PFM.find("TFM")  # get TFM

        # set all the parameters to NONE as a default
        a = None
        f = None
        b = None
        c = None
        ti = None
        fi = None
        throttle = None
        prop_dia = None
        max_eff = None
        p = None
        Hd_turbo = None
        CPSFC = None
        P = None

        kink = {}
        max_power = {}

        if engineType == "JET":
            CT = TFM.find("CT")  # Thrust coefficients

            a = []
            for i in CT.findall("a"):
                a.append(float(i.text))  # C_T polynomial coefficients

            CF = TFM.find("CF")  # Fuel flow coefficients

            f = []
            for i in CF.findall("f"):
                f.append(float(i.text))  # FF polynomial coefficients

            b = {}
            c = {}

            for rating in ["MCMB", "MCRZ", "MTKF"]:
                ENG = TFM.find(rating)
                if ENG is not None:
                    flat_rating = ENG.find("flat_rating")
                    temp_rating = ENG.find("temp_rating")

                    kink[rating] = float(
                        ENG.find("kink").text
                    )  # kink point for Max Climb

                    b[rating] = []
                    for i in flat_rating.findall("b"):
                        b[rating].append(float(i.text))

                    c[rating] = []
                    for i in temp_rating.findall("c"):
                        c[rating].append(float(i.text))

            # Idle data
            LIDL = TFM.find("LIDL")
            CT = LIDL.find("CT")

            ti = []
            for i in CT.findall("ti"):
                ti.append(float(i.text))  # idle thrust coefficients

            CF = LIDL.find("CF")

            fi = []
            for i in CF.findall("fi"):
                fi.append(float(i.text))  # idle fuel flow coefficients

            throttle = {}
            throttle["low"] = float(TFM.find("throttle").find("low").text)
            throttle["high"] = float(TFM.find("throttle").find("high").text)

        elif engineType == "TURBOPROP":
            TPM = PFM.find("TPM")  # get TFM

            prop_dia = float(TPM.find("prop_dia").text)
            max_eff = float(TPM.find("max_eff").text)
            p = {}

            CP = TPM.find("CP")  # Thrust coefficients

            a = []
            for i in CP.findall("a"):
                a.append(float(i.text))  # C_P polynomial coefficients

            CF = TPM.find("CF")  # Fuel flow coefficients

            f = []
            for i in CF.findall("f"):
                f.append(float(i.text))  # FF polynomial coefficients

            for rating in ["MCMB", "MCRZ"]:
                ENG = TPM.find(rating)
                if ENG is not None:
                    r = ENG.find("rating")

                    max_power[rating] = float(ENG.find("max_power").text)

                    p[rating] = []
                    for i in r.findall("p"):
                        p[rating].append(float(i.text))

            # Idle data
            LIDL = TPM.find("LIDL")
            CT = LIDL.find("CT")

            ti = []
            for i in CT.findall("ti"):
                ti.append(float(i.text))  # idle thrust coefficients

            CF = LIDL.find("CF")

            fi = []
            for i in CF.findall("fi"):
                fi.append(float(i.text))  # idle fuel flow coefficients

            throttle = {}
            throttle["low"] = float(TPM.find("throttle").find("low").text)
            throttle["high"] = float(TPM.find("throttle").find("high").text)

        elif engineType == "PISTON":
            PEM = PFM.find("PEM")

            prop_dia = float(PEM.find("prop_dia").text)
            max_eff = float(PEM.find("max_eff").text)
            Hd_turbo = float(PEM.find("Hd_turbo").text)
            CPSFC = float(PEM.find("CPSFC").text)
            P = float(PEM.find("P").text)

        # Parse aerodynamic data
        AFCM = root.find("AFCM")  # get AFCM

        S = float(AFCM.find("S").text)

        HLPosition = {}
        configName = {}
        VFE = {}
        d = {}
        CL_max = {}
        bf = []
        HLids = []

        Mmin = None
        Mmax = None
        CL_Mach0 = None

        CL_clean = None

        for conf in AFCM.findall("Configuration"):
            HLid = float(conf.get("HLid"))
            HLids.append(str(HLid))
            HLPosition[HLid] = float(conf.find("HLPosition").text)
            configName[HLid] = conf.find("name").text
            VFE[HLid] = float(conf.find("vfe").text)
            d[HLid] = {}

            LGUP = conf.find("LGUP")
            LGDN = conf.find("LGDN")

            if LGUP is not None:
                DPM = LGUP.find("DPM_nonclean")

                if DPM is not None:  # DPM is not clean
                    CD = DPM.find("CD_nonclean")
                    if CD is not None:
                        if CD.find("d") is not None:
                            d[HLid]["LGUP"] = []
                            for i in CD.findall("d"):
                                d[HLid]["LGUP"].append(float(i.text))

                else:  # DPM is clean
                    DPM = LGUP.find("DPM_clean")
                    if DPM.find("M_max") is not None:
                        M_max = float(DPM.find("M_max").text)
                    if DPM.find("scalar") is not None:
                        scalar = float(DPM.find("scalar").text)

                    CD = DPM.find("CD_clean")
                    if CD is not None:
                        if CD.find("d") is not None:
                            d[HLid]["LGUP"] = []
                            for i in CD.findall("d"):
                                d[HLid]["LGUP"].append(float(i.text))

                BLM = LGUP.find("BLM")

                if BLM is not None:  # BLM is not clean
                    if BLM.find("CL_max") is not None:
                        if HLid not in CL_max:
                            CL_max[HLid] = {}

                        CL_max[HLid]["LGUP"] = float(BLM.find("CL_max").text)

                else:  # BLM is clean
                    BLM = LGUP.find("BLM_clean")

                    if BLM.find("Mmin") is not None:
                        Mmin = float(BLM.find("Mmin").text)
                    if BLM.find("Mmax") is not None:
                        Mmax = float(BLM.find("Mmax").text)
                    if BLM.find("CL_Mach0") is not None:
                        CL_Mach0 = float(BLM.find("CL_Mach0").text)

                    CL_clean = BLM.find("CL_clean")
                    CL_clean = BLM.find("CL_clean")

                    if CL_clean is not None:
                        for i in CL_clean.findall("bf"):
                            bf.append(float(i.text))

            if LGDN is not None:  # Landing gear NOT allowed in clean configuration
                DPM = LGDN.find("DPM_nonclean")

                if DPM is not None:  # DPM is not clean
                    CD = DPM.find("CD_nonclean")
                    if CD.find("d") is not None:
                        d[HLid]["LGDN"] = []
                        for i in CD.findall("d"):
                            d[HLid]["LGDN"].append(float(i.text))

                BLM = LGDN.find("BLM")

                if BLM is not None:  # BLM is not clean
                    if HLid not in CL_max:
                        CL_max[HLid] = {}
                    CL_max[HLid]["LGDN"] = float(BLM.find("CL_max").text)

        ALM = root.find("ALM")  # get ALM
        DLM = ALM.find("DLM")  # get DLM

        MTOW = float(DLM.find("MTOW").text)
        OEW = float(DLM.find("OEW").text)
        MFL = float(DLM.find("MFL").text)

        MTW = None
        MZFW = None
        MPL = None
        MLW = None

        if DLM.find("MTW") is not None:
            MTW = float(DLM.find("MTW").text)
        if DLM.find("MZFW") is not None:
            MZFW = float(DLM.find("MZFW").text)
        if DLM.find("MPL") is not None:
            MPL = float(DLM.find("MPL").text)
        if DLM.find("MLW") is not None:
            MLW = float(DLM.find("MLW").text)

        GLM = ALM.find("GLM")  # get GLM
        hmo = float(GLM.find("hmo").text)

        if GLM.find("mfa") is not None:
            mfa = float(GLM.find("mfa").text)
        else:
            mfa = None

        KLM = ALM.find("KLM")  # get GLM

        MMO = None
        MLE = None
        VLE = None
        VMO = None

        if KLM.find("mmo") is not None:
            MMO = float(KLM.find("mmo").text)

        if KLM.find("mle") is not None:
            MLE = float(KLM.find("mle").text)

        if KLM.find("vmo") is not None:
            VMO = float(KLM.find("vmo").text)

        if KLM.find("vle") is not None:
            VLE = float(KLM.find("vle").text)

        ground = root.find("Ground")
        dimensions = ground.find("Dimensions")

        span = float(dimensions.find("span").text)
        length = float(dimensions.find("length").text)

        # ARPM model
        ARPM = root.find("ARPM")
        aeroConfSchedule = ARPM.find("AeroConfSchedule")

        # all aerodynamic configurations
        aeroConfig = {}

        for conf in aeroConfSchedule.findall("AeroPhase"):
            name = conf.find("name").text
            HLid = float(conf.find("HLid").text)
            LG = "LG" + conf.find("LG").text
            aeroConfig[name] = {"HLid": HLid, "LG": LG}

        speedScheduleList = ARPM.find("SpeedScheduleList")
        SpeedSchedule = speedScheduleList.find("SpeedSchedule")

        # all phases of flight
        speedSchedule = {}
        for phaseOfFlight in SpeedSchedule.findall("SpeedPhase"):
            name = phaseOfFlight.find("name").text
            CAS1 = float(phaseOfFlight.find("CAS1").text)
            CAS2 = float(phaseOfFlight.find("CAS2").text)
            M = float(phaseOfFlight.find("M").text)
            speedSchedule[name] = {"CAS1": CAS1, "CAS2": CAS2, "M": M}

        # Single row dataframe
        data = {
            "acName": [acName],
            "model": [model],
            "engineType": [engineType],
            "engines": [engines],
            "ICAO": [ICAO],
            "WTC": [WTC],
            "MREF": [MREF],
            "WREF": [WREF],
            "LHV": [LHV],
            "n_eng": [n_eng],
            "rho": [rho],
            "TFA": [TFA],
            "p_delta": [p_delta],
            "p_theta": [p_theta],
            "a": [a],
            "f": [f],
            "b": [b],
            "c": [c],
            "ti": [ti],
            "fi": [fi],
            "throttle": [throttle],
            "prop_dia": [prop_dia],
            "max_eff": [max_eff],
            "p": [p],
            "Hd_turbo": [Hd_turbo],
            "CPSFC": [CPSFC],
            "P": [P],
            "kink": [kink],
            "max_power": [max_power],
            "S": [S],
            "HLPosition": [HLPosition],
            "configName": [configName],
            "VFE": [VFE],
            "d": [d],
            "CL_max": [CL_max],
            "bf": [bf],
            "HLids": [HLids],
            "Mmin": [Mmin],
            "Mmax": [Mmax],
            "CL_Mach0": [CL_Mach0],
            "CL_clean": [CL_clean],
            "M_max": [M_max],
            "scalar": [scalar],
            "MTOW": [MTOW],
            "OEW": [OEW],
            "MFL": [MFL],
            "MTW": [MTW],
            "MZFW": [MZFW],
            "MPL": [MPL],
            "MLW": [MLW],
            "hmo": [hmo],
            "mfa": [mfa],
            "MMO": [MMO],
            "MLE": [MLE],
            "VLE": [VLE],
            "VMO": [VMO],
            "span": [span],
            "length": [length],
            "aeroConfig": [aeroConfig],
            "speedSchedule": [speedSchedule],
        }
        df_single = pd.DataFrame(data)

        return df_single

    @staticmethod
    def parseGPF(filePath, badaVersion):
        """This function parses BADA4 GPF xml formatted file

        :param filePath: path to the BADA4 xml formatted file.
        :param acName: name of Aircraft for BADA4 xml formatted file.
        :type filePath: str.
        :type acName: str.
        :raises: IOError
        """

        filename = os.path.join(filePath, "BADA4", badaVersion, "GPF.xml")

        tree = ET.parse(filename)

        try:
            tree = ET.parse(filename)
            root = tree.getroot()
        except:
            raise IOError(filename + " not found or in correct format")

        CVminTO = float(root.find("CVminTO").text)  # CVminTO
        CVmin = float(root.find("CVmin").text)  # CVmin

        HmaxList = root.find("HmaxList")
        HmaxPhase = {}  # phase of flight
        for Hmax_Phase in HmaxList.findall("HmaxPhase"):
            phaseName = Hmax_Phase.find("Phase").text
            Hmax = float(Hmax_Phase.find("Hmax").text)
            HmaxPhase[phaseName] = Hmax

        V_des = {}  # V_des
        V_cl = {}  # V_cl
        VdList = root.find("VdList")
        for VdPhase in VdList.findall("VdPhase"):
            Phase = VdPhase.find("Phase")
            name = Phase.find("name").text
            index = int(Phase.find("index").text)
            Vd = float(VdPhase.find("Vd").text)

            if name == "CL":
                V_cl[index] = Vd
            elif name == "DES":
                V_des[index] = Vd

        # Single row dataframe
        data = {
            "CVminTO": [CVminTO],
            "CVmin": [CVmin],
            "HmaxPhase": [HmaxPhase],
            "V_des": [V_des],
            "V_cl": [V_cl],
        }
        df_single = pd.DataFrame(data)

        return df_single

    @staticmethod
    def combineXML_GPF(XMLDataFrame, GPFDataframe):
        """This function combines 2 dataframes, the parsed aircraft XML file
        and parsed GPF file
        """

        # Combine data with GPF data (temporary solution)
        combined_df = pd.concat(
            [XMLDataFrame.reset_index(drop=True), GPFDataframe.reset_index(drop=True)],
            axis=1,
        )

        return combined_df

    @staticmethod
    def parseAll(badaVersion, filePath=None):
        """This function parses all BADA4 xml formatted file and stores
        all data in the final dataframe containing all the BADA data.

        :param filePath: path to the BADAH Synonym xml formatted file.
        :type filePath: str.
        :raises: IOError
        """

        if filePath == None:
            filePath = configuration.getAircraftPath()
        else:
            filePath = filePath

        # parsing GPF file
        GPFparsedDataframe = Parser.parseGPF(filePath, badaVersion)

        # retrieving mapping data
        code_fileName = Parser.readMappingFile(filePath, badaVersion)

        # get names of all the folders in the main BADA model folder to search for XML files
        subfolders = Parser.list_subfolders(
            os.path.join(filePath, "BADA4", badaVersion)
        )

        # Initialize an empty list to collect DataFrames
        mapping_dfs = []

        if code_fileName:
            for code in code_fileName:
                file = code_fileName[code]

                if file in subfolders:
                    # Parse the original XML of a model
                    df = Parser.parseXML(filePath, badaVersion, file)

                    # Rename 'acName' in the DataFrame to match the code model name
                    df.at[0, "acName"] = code

                    # Combine data with GPF data (temporary solution)
                    combined_df = Parser.combineXML_GPF(df, GPFparsedDataframe)

                    # Drop columns that are all NaN
                    combined_df = combined_df.dropna(axis=1, how="all")

                    # Check if combined_df is not empty
                    if not combined_df.empty:
                        mapping_dfs.append(combined_df)

        # Concatenate all collected DataFrames
        if mapping_dfs:
            merged_mapping_df = pd.concat(mapping_dfs, ignore_index=True)
        else:
            merged_mapping_df = pd.DataFrame()

        # Initialize an empty list to collect DataFrames
        original_dfs = []

        for file in subfolders:
            # Parse the original XML of a model
            df = Parser.parseXML(filePath, badaVersion, file)

            # Combine data with GPF data (temporary solution)
            combined_df = pd.concat(
                [df.reset_index(drop=True), GPFparsedDataframe.reset_index(drop=True)],
                axis=1,
            )

            # Drop columns that are all NaN
            combined_df = combined_df.dropna(axis=1, how="all")

            # Check if combined_df is not empty
            if not combined_df.empty:
                original_dfs.append(combined_df)

        # Concatenate all collected DataFrames
        if original_dfs:
            merged_original_df = pd.concat(original_dfs, ignore_index=True)
        else:
            merged_original_df = pd.DataFrame()

        # Merge mapping and original aircraft models
        merged_final_df = pd.concat(
            [merged_original_df, merged_mapping_df], ignore_index=True
        )

        return merged_final_df

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


class BADA4(Airplane):
    """This class implements the part of BADA4 performance model that will be used in other classes following the BADA4 manual.

    :param AC: parsed aircraft.
    :type AC: bada4.Parse.
    """

    def __init__(self, AC):
        super().__init__()
        self.AC = AC

    def CL(self, delta, mass, M, nz=1.0):
        """This function computes the lift coefficient

        :param M: mach airspeed [-].
        :param delta: normalised air pressure [-].
        :param mass: aircraft mass [kg].
        :param nz: load factor [-].
        :type M: float.
        :type delta: float.
        :type mass: float.
        :type nz: float.
        :returns: Lift coefficient [-].
        :rtype: float.
        """

        return (
            2
            * mass
            * const.g
            * nz
            / (delta * const.p_0 * const.Agamma * pow(M, 2) * self.AC.S)
        )

    def CLPoly(self, M):
        """This function computes the lift coefficient polynom

        :param M: mach airspeed [-].
        :type M: float.
        :returns: Lift coefficient polynom [-].
        :rtype: float.
        """

        CLpoly = 0.0
        for i in range(5):
            CLpoly += self.AC.bf[i] * pow(M, i)

        return CLpoly

    def CLmax(self, M, HLid, LG):
        """This function computes the maximum lift coefficient in clean and non-clean configuration

        :param M: Mach number [-].
        :param HLid: high lift devices position [-]
        :param LG: landing gear position, [LGUP/LGDN] [-]
        :type M: float.
        :type HLid: float.
        :type LG: string.
        :returns: maximum lift coefficient [-].
        :rtype: float.
        """

        CLmax = 0.0

        # CLmax available - non clean configuration
        if HLid in self.AC.CL_max and LG in self.AC.CL_max[HLid]:
            CLmax = self.AC.CL_max[HLid][LG]

        # CLmax unavailable - clean configuration
        elif HLid == 0 and LG == "LGUP":
            if self.AC.CL_clean is not None:
                if M < self.AC.Mmin:
                    CLmax = self.AC.CL_Mach0 + (M / self.AC.Mmin) * (
                        self.AC.CLPoly(self.AC.Mmin) - self.AC.CL_Mach0
                    )
                elif M > self.AC.Mmax:
                    CLder = (
                        self.AC.bf[1]
                        + 2 * self.AC.bf[2] * self.AC.Mmax
                        + 3 * self.AC.bf[3] * self.AC.Mmax * self.AC.Mmax
                        + 4 * self.AC.bf[4] * self.AC.Mmax * self.AC.Mmax * self.AC.Mmax
                    )
                    CLmax = (
                        self.CLPoly(self.AC.Mmax)
                        + self.CLPoly(M - self.AC.Mmax) * CLder
                    )
                else:
                    CLmax = self.CLPoly(M)
            else:
                CLmax = self.AC.CL_max[0]["LGUP"]

        return CLmax

    def CF_idle(self, delta, theta, M):
        """This function computes the idle rating model for JET and TURBOPROP engine

        :param delta: Normalised pressure [-].
        :param theta: Normalised temperature [-].
        :param M: Mach speed [-].
        :type delta: float.
        :type theta: float.
        :type M: float.
        :returns: idle Fuel flow coefficient [-].
        :rtype: float
        """
        if self.AC.engineType == "JET":
            CF_idle = 0.0
            for i in range(0, 3):
                for j in range(0, 3):
                    CF_idle += self.AC.fi[i * 3 + j] * (delta**j) * (M**i)

            if self.AC.BADAVersion == "4.2" or self.AC.BADAVersion == "DUMMY":
                CF_idle = CF_idle * pow(delta, -1) * pow(theta, -0.5)
            elif self.AC.BADAVersion == "4.3":
                CF_idle = CF_idle * pow(delta, -1)

        elif self.AC.engineType == "TURBOPROP":
            CF_idle = 0.0
            for i in range(0, 3):
                for j in range(0, 3):
                    CF_idle += self.AC.fi[i * 3 + j] * (delta**j) * (M**i)
            CF_idle += (
                self.AC.fi[9] * theta
                + self.AC.fi[10] * (theta**2)
                + self.AC.fi[11] * M * theta
                + self.AC.fi[12] * M * delta * sqrt(theta)
                + self.AC.fi[13] * M * delta * theta
            )

            if self.AC.BADAVersion == "4.2" or self.AC.BADAVersion == "DUMMY":
                CF_idle = CF_idle * pow(delta, -1) * pow(theta, -0.5)
            elif self.AC.BADAVersion == "4.3":
                CF_idle = CF_idle * pow(delta, -1)

        return CF_idle

    def CF(self, delta, theta, DeltaTemp, **kwargs):
        """This function computes the fuel flow coefficient

        :param rating: Throttle setting {MCMB,MCRZ,LIDL}
        :param deltaT: direct throttle parameter [-]
        :param delta: Normalised pressure [-].
        :param theta: Normalised temperature [-].
        :param M: Mach number [-].
        :param DeltaTemp: Temperature deviation with respect to ISA [K].
        :type rating: string
        :type deltaT: float.
        :type delta: float.
        :type theta: float.
        :type M: float.
        :type DeltaTemp: float.
        :returns: Fuel flow coefficient [-].
        :rtype: float
        """

        if self.AC.engineType == "JET":
            M = checkArgument("M", **kwargs)

            # when idle rating is used
            CF_idle = self.CF_idle(delta=delta, theta=theta, M=M)

            # for adaptive thrust calcualation if CT is an input
            if "CT" in kwargs:
                CT = kwargs.get("CT")

                CF_gen_rating = 0.0
                for i in range(0, 5):
                    for j in range(0, 5):
                        CF_gen_rating += self.AC.f[i * 5 + j] * (CT**j) * (M**i)

                CF = max(CF_gen_rating, CF_idle)

            # rating as input parameter
            elif "rating" in kwargs:
                rating = checkArgument("rating", **kwargs)

                # in case MCRZ rating is not defined, we switch to MCMB
                if rating == "MCRZ" and rating not in self.AC.kink.keys():
                    rating = "MCMB"

                if rating not in self.AC.kink.keys() and rating != "LIDL":
                    raise ValueError("Unknown engine rating " + rating)

                if rating == "LIDL":
                    CF = CF_idle

                elif rating in self.AC.kink.keys():
                    # when non-idle rating is used
                    CF_gen_rating = 0.0
                    CT_rating = self.CT(
                        rating=rating,
                        theta=theta,
                        delta=delta,
                        M=M,
                        DeltaTemp=DeltaTemp,
                    )
                    for i in range(0, 5):
                        for j in range(0, 5):
                            CF_gen_rating += (
                                self.AC.f[i * 5 + j] * (CT_rating**j) * (M**i)
                            )

                    CF = max(CF_gen_rating, CF_idle)

            # deltaT - direct throttle as input parameter
            elif "deltaT" in kwargs:
                # when no rating is used
                CF_gen_deltaT = 0.0
                CT_deltaT = self.CT(
                    deltaT=deltaT, theta=theta, delta=delta, M=M, DeltaTemp=DeltaTemp
                )
                for i in range(0, 5):
                    for j in range(0, 5):
                        CF_gen_deltaT += self.AC.f[i * 5 + j] * (CT_deltaT**j) * (M**i)

                CF = max(CF_gen_deltaT, CF_idle)

        elif self.AC.engineType == "TURBOPROP":
            M = checkArgument("M", **kwargs)

            # when idle rating is used
            CF_idle = self.CF_idle(delta=delta, theta=theta, M=M)

            # for adaptive thrust calcualation if CT is an input
            if "CT" in kwargs:
                CT = kwargs.get("CT")
                CP = self.CP(CT=CT, M=M)

                CF_gen_rating = 0.0
                for i in range(0, 5):
                    for j in range(0, 5):
                        CF_gen_rating += self.AC.f[i * 5 + j] * (CP**j) * (M**i)

                CF = max(CF_gen_rating, CF_idle)

            # rating as input parameter
            elif "rating" in kwargs:
                rating = checkArgument("rating", **kwargs)

                # in case MCRZ rating is not defined, we switch to MCMB
                if rating == "MCRZ" and rating not in self.AC.max_power.keys():
                    rating = "MCMB"

                if rating not in self.AC.max_power.keys() and rating != "LIDL":
                    raise ValueError("Unknown engine rating " + rating)

                if rating == "LIDL":
                    CF = CF_idle

                elif rating in self.AC.max_power.keys():
                    # when non-idle rating is used
                    CF_gen_rating = 0.0
                    CP_rating = self.CP(rating=rating, theta=theta, delta=delta, M=M)
                    for i in range(0, 5):
                        for j in range(0, 5):
                            CF_gen_rating += (
                                self.AC.f[i * 5 + j] * (CP_rating**j) * (M**i)
                            )

                    CF = max(CF_gen_rating, CF_idle)

            # deltaT - direct throttle as input parameter
            elif "deltaT" in kwargs:
                # when non-idle rating is used
                CF_gen_deltaT = 0.0
                CP_deltaT = self.CP(deltaT=deltaT, M=M)

                CT_deltaT = self.CT(
                    deltaT=deltaT, theta=theta, delta=delta, M=M, DeltaTemp=DeltaTemp
                )
                for i in range(0, 5):
                    for j in range(0, 5):
                        CF_gen_deltaT += self.AC.f[i * 5 + j] * (CP_deltaT**j) * (M**i)

                CF = max(CF_gen_deltaT, CF_idle)

        elif self.AC.engineType == "PISTON":
            sigma = atm.sigma(theta=theta, delta=delta)

            # for adaptive thrust calcualation if CT is an input
            if "CT" in kwargs:
                CT = kwargs.get("CT")
                M = checkArgument("M", **kwargs)

                deltaT_vec = np.arange(0.01, 1.01, 0.01)

                CT_diff = []
                for k in range(len(deltaT_vec)):
                    CT_k = self.CT(
                        theta=theta,
                        delta=delta,
                        deltaT=deltaT_vec[k],
                        M=M,
                        DeltaTemp=DeltaTemp,
                    )
                    CT_diff.append(abs(CT_k - CT))

                CT_min_idx = CT_diff.index(min(CT_diff))
                deltaT = deltaT_vec[CT_min_idx]

                CP = self.CP(
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                    deltaT=deltaT,
                    DeltaTemp=DeltaTemp,
                )

            # rating as input parameter
            elif "rating" in kwargs:
                rating = checkArgument("rating", **kwargs)

                CP = self.CP(
                    rating=rating,
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                    DeltaTemp=DeltaTemp,
                )

            CF = self.AC.CPSFC * CP / (delta * sqrt(theta))

        return CF

    def CT(self, delta, **kwargs):
        """This function computes the thrust coefficient

        :param rating: Throttle setting {MCMB,MTKF,MCRZ,LIDL}
        :param deltaT: direct throttle parameter [-]
        :param delta: Normalised pressure [-].
        :param theta: Normalised temperature [-].
        :param M: Mach number [-].
        :param DeltaTemp: Temperature deviation with respect to ISA [K].
        :type rating: string
        :type deltaT: float.
        :type delta: float.
        :type theta: float.
        :type M: float.
        :type DeltaTemp: float.
        :returns: thrust coefficient [-].
        :rtype: float
        :raise: ValueError
        """

        if "Thrust" in kwargs:
            Thrust = kwargs.get("Thrust")
            CT = Thrust / (delta * self.AC.WREF)

            return CT

        else:
            theta = checkArgument("theta", **kwargs)
            DeltaTemp = checkArgument("DeltaTemp", **kwargs)
            M = checkArgument("M", **kwargs)

        if self.AC.engineType == "JET":
            # rating as input parameter
            if "rating" in kwargs:
                rating = checkArgument("rating", **kwargs)

                # in case MCRZ rating is not defined, we switch to MCMB
                if rating == "MCRZ" and rating not in self.AC.kink.keys():
                    rating = "MCMB"

                if rating == "MTKF" and rating not in self.AC.kink.keys():
                    rating = "MCMB"

                if rating not in self.AC.kink.keys() and rating != "LIDL":
                    raise ValueError("Unknown engine rating " + rating)

                if rating == "LIDL":
                    CT = self.CT_LIDL(delta=delta, M=M)

                elif rating in self.AC.kink.keys():
                    CT = self.CT_nonLIDL(
                        rating=rating,
                        theta=theta,
                        delta=delta,
                        M=M,
                        DeltaTemp=DeltaTemp,
                    )

            # deltaT - direct throttle as input parameter
            elif "deltaT" in kwargs:
                deltaT = checkArgument("deltaT", **kwargs)

                CT = 0.0
                for i in range(0, 6):
                    for j in range(0, 6):
                        CT += self.AC.a[i * 6 + j] * (M**j) * (deltaT**i)

                # limit CT with CT_LIDL and CT_MCMB
                if CT > self.CT_nonLIDL(
                    rating="MCMB", theta=theta, delta=delta, M=M, DeltaTemp=DeltaTemp
                ):
                    raise ValueError(
                        "Throttle parameter value result in CT > CT_MCMB" + deltaT
                    )
                elif CT < self.CT_LIDL(delta=delta, M=M):
                    raise ValueError(
                        "Throttle parameter value result in CT < CT_LIDL" + deltaT
                    )

        elif self.AC.engineType == "TURBOPROP":
            # rating as input parameter
            if "rating" in kwargs:
                rating = checkArgument("rating", **kwargs)

                # in case MCRZ rating is not defined, we switch to MCMB
                if rating == "MCRZ" and rating not in self.AC.max_power.keys():
                    rating = "MCMB"

                if rating == "MTKF" and rating not in self.AC.kink.keys():
                    rating = "MCMB"

                if rating not in self.AC.max_power.keys() and rating != "LIDL":
                    raise ValueError("Unknown engine rating " + rating)

                if rating == "LIDL":
                    CT = self.CT_LIDL(theta=theta, delta=delta, M=M)

                elif rating in self.AC.max_power.keys():
                    CT = self.CT_nonLIDL(rating=rating, theta=theta, delta=delta, M=M)

            # deltaT - direct throttle as input parameter
            elif "deltaT" in kwargs:
                deltaT = checkArgument("deltaT", **kwargs)

                CP = self.CP(deltaT=deltaT, M=M)
                CT = CP / M

                if CT > self.CT_nonLIDL(rating="MCMB", theta=theta, delta=delta, M=M):
                    raise ValueError(
                        "Throttle parameter value result in CT > CT_MCMB" + deltaT
                    )
                elif CT < self.CT_LIDL(theta=theta, delta=delta, M=M):
                    raise ValueError(
                        "Throttle parameter value result in CT < CT_LIDL" + deltaT
                    )

        elif self.AC.engineType == "PISTON":
            sigma = atm.sigma(theta=theta, delta=delta)

            # rating as input parameter
            if "rating" in kwargs:
                rating = checkArgument("rating", **kwargs)
                if rating not in ["MCMB", "MCRZ"] and rating != "LIDL":
                    raise ValueError("Unknown engine rating " + rating)

                CP = self.CP(
                    rating=rating,
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                    DeltaTemp=DeltaTemp,
                )

            # deltaT - direct throttle as input parameter
            elif "deltaT" in kwargs:
                deltaT = checkArgument("deltaT", **kwargs)

                CP = self.CP(
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                    deltaT=deltaT,
                    DeltaTemp=DeltaTemp,
                )

            CT = self.CT_nonLIDL(theta=theta, delta=delta, M=M, CP=CP)

        return CT

    def CT_LIDL(self, **kwargs):
        """This function computes the thrust coefficient for LIDL rating

        :param delta: Normalised pressure [-].
        :param theta: Normalised temperature [-].
        :param M: Mach number [-].
        :param CP: power coefficient [-].
        :type delta: float.
        :type theta: float.
        :type M: float.
        :type CP: float.
        :returns: LIDL thrust coefficient [-].
        :rtype: float
        """

        if self.AC.engineType == "JET":
            delta = checkArgument("delta", **kwargs)
            M = checkArgument("M", **kwargs)

            CT = 0.0
            for i in range(0, 3):
                for j in range(0, 4):
                    CT += self.AC.ti[i * 4 + j] * pow(delta, j - 1) * (M**i)

        elif self.AC.engineType == "TURBOPROP":
            theta = checkArgument("theta", **kwargs)
            delta = checkArgument("delta", **kwargs)
            M = checkArgument("M", **kwargs)

            CT = 0.0
            for i in range(0, 3):
                for j in range(0, 4):
                    CT += self.AC.ti[i * 4 + j] * pow(delta, j - 1) * (M**i)

            CT += (
                self.AC.ti[12] * sqrt(theta)
                + self.AC.ti[13] * theta
                + self.AC.ti[14] / sqrt(theta)
                + self.AC.ti[15] * theta**2
            )
            CT += (
                self.AC.ti[16] / delta
                + self.AC.ti[17] * delta
                + self.AC.ti[18] * delta**2
                + self.AC.ti[19] * M
                + self.AC.ti[20] * M**2
            ) / sqrt(theta)
            CT += (
                self.AC.ti[21] / M + self.AC.ti[22] * delta / M + self.AC.ti[23] * M**3
            )
            CT += (
                self.AC.ti[24] * M
                + self.AC.ti[25] * M**2
                + self.AC.ti[26]
                + self.AC.ti[27] * M / delta
            ) / theta
            CT += (
                self.AC.ti[28] * M / (delta * theta**2)
                + self.AC.ti[29] * M**2 / (delta * theta**2)
                + self.AC.ti[30] * M**2 / (delta * sqrt(theta))
                + self.AC.ti[31] * delta / theta
            )

        elif self.AC.engineType == "PISTON":
            theta = checkArgument("theta", **kwargs)
            delta = checkArgument("delta", **kwargs)
            CP = checkArgument("CP", **kwargs)
            M = checkArgument("M", **kwargs)

            CT = self.CT_nonLIDL(theta=theta, delta=delta, M=M, CP=CP)

        return CT

    def CT_nonLIDL(self, theta, delta, M, **kwargs):
        """This function computes the thrust coefficient for non-LIDL rating {MCMB,MCRZ}

        :param delta: Normalised pressure [-].
        :param theta: Normalised temperature [-].
        :param M: Mach number [-].
        :param CP: power coefficient [-].
        :type delta: float.
        :type theta: float.
        :type M: float.
        :type CP: float.
        :returns: LIDL thrust coefficient [-].
        :rtype: float
        """

        if self.AC.engineType == "JET":
            rating = checkArgument("rating", **kwargs)
            DeltaTemp = checkArgument("DeltaTemp", **kwargs)

            # deltaT below kink point -> flat-rated area
            deltaTFlat = 0.0
            for i in range(0, 6):
                for j in range(0, 6):
                    deltaTFlat += self.AC.b[rating][i * 6 + j] * (M**j) * (delta**i)

            # deltaT above kink point -> temperature-rated area
            thetaT = theta * (1 + (M**2) * (const.Agamma - 1.0) / 2.0)
            deltaTTemp = 0.0
            for i in range(0, 5):
                for j in range(0, 5):
                    deltaTTemp += self.AC.c[rating][i * 5 + j] * (M**j) * (thetaT**i)

            for i in range(5, 9):
                for j in range(0, 5):
                    deltaTTemp += (
                        self.AC.c[rating][i * 5 + j] * (M**j) * (delta ** (i - 4))
                    )

            # compute deltaT according to DeltaTemp with respect to kink point
            if DeltaTemp <= self.AC.kink[rating]:
                deltaT = deltaTFlat
            else:
                deltaT = deltaTTemp

            CT = 0.0
            for i in range(0, 6):
                for j in range(0, 6):
                    CT += self.AC.a[i * 6 + j] * (M**j) * (deltaT**i)

        elif self.AC.engineType == "TURBOPROP":
            rating = checkArgument("rating", **kwargs)

            CP = self.CP(rating=rating, theta=theta, delta=delta, M=M)
            CT = CP / M

        elif self.AC.engineType == "PISTON":
            CP = checkArgument("CP", **kwargs)

            sigma = atm.sigma(theta=theta, delta=delta)
            Wp = self.AC.WREF * const.a_0 * CP
            TAS = atm.mach2Tas(Mach=M, theta=theta)
            propEff = self.propEfficiency(Wp=Wp * self.AC.n_eng, sigma=sigma, tas=TAS)
            CT = propEff * const.a_0 * CP / (delta * TAS)

        return CT

    def CPmax(self, rating, delta, theta, M):
        """This function computes the maximum engine power coefficient for Turboprop

        :param rating: Throttle setting {MCMB,MCRZ}
        :param delta: Normalised pressure [-].
        :param theta: Normalised temperature [-].
        :param M: Mach number [-].
        :type rating: str
        :type delta: float.
        :type theta: float.
        :type M: float.
        :returns: Maximum engine power coefficient [-].
        :rtype: float
        :raise: ValueError
        """

        if self.AC.engineType == "TURBOPROP":
            if rating not in self.AC.max_power.keys():
                raise ValueError("Unknown engine rating " + rating)

            Wpmax = self.AC.max_power[rating]
            aSound = atm.aSound(theta=theta)
            tas = atm.mach2Tas(Mach=M, theta=theta)
            sigma = atm.sigma(theta=theta, delta=delta)
            propEff = self.propEfficiency(Wp=Wpmax, sigma=sigma, tas=tas)
            if propEff is None:
                return None
            CPmax = Wpmax * propEff / (delta * self.AC.WREF * aSound)
        else:
            raise ValueError("CPmax implemented only for turboprop")

        return CPmax

    def CP(self, **kwargs):
        """This function computes the power coefficient

        :param rating: Throttle setting {MCMB,MCRZ,LIDL}
        :param deltaT: direct throttle parameter [-]
        :param delta: Normalised pressure [-].
        :param theta: Normalised temperature [-].
        :param sigma: Normalised density [-].
        :param M: Mach number [-].
        :type rating: string
        :type deltaT: float.
        :type delta: float.
        :type theta: float.
        :type sigma: float.
        :type M: float.
        :returns: thrust coefficient [-].
        :rtype: float
        :raise: ValueError
        """

        if self.AC.engineType == "TURBOPROP":
            M = checkArgument("M", **kwargs)

            # CT as input parameter
            # computes the power coefficient from thrust coefficient assuming efficiency of 1
            if "CT" in kwargs:
                CT = kwargs.get("CT")
                CP = CT * M
                return CP

            # rating as input parameter
            elif "rating" in kwargs:
                rating = checkArgument("rating", **kwargs)

                # in case MCRZ rating is not defined, we switch to MCMB
                if rating == "MCRZ" and rating not in self.AC.max_power.keys():
                    rating = "MCMB"

                delta = checkArgument("delta", **kwargs)
                theta = checkArgument("theta", **kwargs)

                deltaT = 0.0
                for i in range(0, 6):
                    for j in range(0, 6):
                        deltaT += self.AC.p[rating][i * 6 + j] * (M**j) * (theta**i)

                CP = 0.0
                for i in range(0, 6):
                    for j in range(0, 6):
                        CP += self.AC.a[i * 6 + j] * (M**j) * (deltaT**i)

                CPmax = self.CPmax(rating=rating, theta=theta, delta=delta, M=M)
                CP = min(CP, CPmax)

            # deltaT - direct throttle as input parameter
            elif "deltaT" in kwargs:
                deltaT = checkArgument("deltaT", **kwargs)

                CP = 0.0
                for i in range(0, 6):
                    for j in range(0, 6):
                        CP += self.AC.a[i * 6 + j] * (M**j) * (deltaT**i)

        elif self.AC.engineType == "PISTON":
            delta = checkArgument("delta", **kwargs)
            theta = checkArgument("theta", **kwargs)
            sigma = checkArgument("sigma", **kwargs)
            DeltaTemp = checkArgument("DeltaTemp", **kwargs)

            if self.AC.Hd_turbo <= 0:
                theta_turbo = atm.theta(
                    conv.ft2m(self.AC.Hd_turbo), DeltaTemp=DeltaTemp
                )
                delta_turbo = atm.delta(
                    conv.ft2m(self.AC.Hd_turbo), DeltaTemp=DeltaTemp
                )

                # ensure that all real sigmas are smaller than sigma_turbo
                sigma_turbo = float("Inf")

            else:
                sigma_turbo = atm.sigma(theta=theta_turbo, delta=delta_turbo)

            CPmaxStdMSL = (
                conv.hp2W(self.AC.P) * self.AC.n_eng / (self.AC.WREF * const.a_0)
            )

            # deltaT - direct throttle as input parameter
            if "deltaT" in kwargs:
                deltaT = checkArgument("deltaT", **kwargs)

            # rating as input parameter
            elif "rating" in kwargs:
                rating = checkArgument("rating", **kwargs)

                if rating == "LIDL":
                    if self.AC.BADAVersion == "4.2" or self.AC.BADAVersion == "DUMMY":
                        deltaT = 0.0
                    elif self.AC.BADAVersion == "4.3":
                        deltaT = 0.1
                elif rating == "MCMB" or rating == "MCRZ":
                    deltaT = 1.0

            CPstdMSL = CPmaxStdMSL * deltaT

            if sigma >= sigma_turbo:
                CP = CPstdMSL
            else:
                CP = min(
                    CPstdMSL,
                    CPmaxStdMSL
                    * delta
                    * sqrt(theta_turbo)
                    / (delta_turbo * sqrt(theta)),
                )

        return CP

    def CDClean(self, CL, M):
        """This function computes the drag coefficient in clean configuration

        :param M: Mach number [-].
        :param CL: Lift coefficient [-].
        :type M: float.
        :type CL: float.
        :returns: Drag coefficient in clean configuration [-].
        :rtype: float
        """

        param = 1 - M * M

        d = self.AC.d[0]["LGUP"]
        C0 = (
            d[0]
            + d[1] / sqrt(param)
            + d[2] / (param)
            + d[3] / pow(param, 3.0 / 2.0)
            + d[4] / pow(param, 2.0)
        )
        C2 = (
            d[5]
            + d[6] / pow(param, 3.0 / 2.0)
            + d[7] / pow(param, 3.0)
            + d[8] / pow(param, 9.0 / 2.0)
            + d[9] / pow(param, 6.0)
        )
        C6 = (
            d[10]
            + d[11] / pow(param, 7.0)
            + d[12] / pow(param, 15.0 / 2.0)
            + d[13] / pow(param, 8.0)
            + d[14] / pow(param, 17.0 / 2.0)
        )

        CD_clean = self.AC.scalar * (C0 + C2 * CL * CL + C6 * pow(CL, 6))
        return CD_clean

    def CD(
        self, HLid, LG, CL, M, speedBrakes={"deployed": False, "value": 0.03}, **kwargs
    ):
        """This function computes the drag coefficient

        :param M: Mach number [-].
        :param CL: Lift coefficient [-].
        :param config: aircraft aerodynamic configuration [TO/IC/CR/AP/LD][-]
        :param HLid: high lift devices position [-]
        :param LG: landing gear position, [LGUP/LGDN] [-]
        :param speedBrakes: speed brakes used or not [-].
        :type M: float.
        :type CL: float.
        :type speedBrakes: boolean.
        :type config: string.
        :type HLid: float.
        :type LG: string.
        :returns: Drag coefficient [-].
        :rtype: float
        """

        # clean configuration
        if HLid == 0 and LG == "LGUP":
            # below Mmax
            if M <= self.AC.M_max:
                CD = self.AC.CDClean(M=M, CL=CL)
            # above Mmax (accounting for air compresibility)
            else:
                CD = self.CDClean(M=self.AC.M_max - 0.01, CL=CL) + pow(
                    (M - (self.AC.M_max - 0.01)) / 0.01, 3 / 2
                ) * (
                    self.CDClean(M=self.AC.M_max, CL=CL)
                    - self.CDClean(M=self.AC.M_max - 0.01, CL=CL)
                )
        # non-clean configuration
        else:
            CD = (
                self.AC.d[HLid][LG][0]
                + self.AC.d[HLid][LG][1] * CL
                + self.AC.d[HLid][LG][2] * CL * CL
            )

        # implementation of a simple speed brakes model
        if speedBrakes["deployed"]:
            if speedBrakes["value"] is not None:
                CD = CD + speedBrakes["value"]
            else:
                CD = CD + 0.03
            return CD

        # calculation of drag coefficient in transition for HLid assuming LG is not changing
        if "HLid_init" in kwargs and "HLid_final" in kwargs:
            HLid_init = checkArgument("HLid_init", **kwargs)
            HLid_final = checkArgument("HLid_final", **kwargs)
            LG_init = LG
            LG_final = LG

            if HLid_init == 0 and LG_init == "LGUP":
                CD_init = self.AC.CDClean(M=M, CL=CL)
            else:
                CD_init = (
                    self.AC.d[HLid_init][LG_init][0]
                    + self.AC.d[HLid_init][LG_init][1] * CL
                    + self.AC.d[HLid_init][LG_init][2] * CL
                )

            if HLid_final == 0 and LG_final == "LGUP":
                CD_final = self.CDClean(M=M, CL=CL)
            else:
                CD_final = (
                    self.AC.d[HLid_final][LG_final][0]
                    + self.AC.d[HLid_final][LG_final][1] * CL
                    + self.AC.d[HLid_final][LG_final][2] * CL
                )

            # linear interpolation
            xp = [HLid_init, HLid_final]
            fp = [CD_init, CD_final]
            CD = np.interp(HLid, xp, fp)

        # calculation of drag coefficient in transition for LG assuming HLid is not changing
        if "LG_init" in kwargs and "LG_final" in kwargs:
            LG_init = checkArgument("LG_init", **kwargs)
            LG_final = checkArgument("LG_final", **kwargs)
            HLid_init = HLid
            HLid_final = HLid

            if HLid_init == 0 and LG_init == "LGUP":
                CD_init = self.CDClean(M=M, CL=CL)
            else:
                CD_init = (
                    self.AC.d[HLid_init][LG_init][0]
                    + self.AC.d[HLid_init][LG_init][1] * CL
                    + self.AC.d[HLid_init][LG_init][2] * CL
                )

            if HLid_final == 0 and LG_final == "LGUP":
                CD_final = self.CDClean(M=M, CL=CL)
            else:
                CD_final = (
                    self.AC.d[HLid_final][LG_final][0]
                    + self.AC.d[HLid_final][LG_final][1] * CL
                    + self.AC.d[HLid_final][LG_final][2] * CL
                )

            # linear interpolation
            xp = [HLid_init, HLid_final]
            fp = [CD_init, CD_final]
            CD = np.interp(HLid, xp, fp)

        return CD

    def L(self, delta, M, CL):
        """This function computes the aerodynamic lift

        :param M: Mach airspeed [-].
        :param delta: normalised air pressure [-].
        :param CL: Lift coefficient [-].
        :type M: float.
        :type delta: float.
        :type CL: float.
        :returns: Aerodynamic lift [N].
        :rtype: float.
        """

        return 0.5 * delta * const.p_0 * const.Agamma * M * M * self.AC.S * CL

    def D(self, delta, M, CD):
        """This function computes the aerodynamic drag

        :param M: Mach airspeed [-].
        :param delta: normalised air pressure [-].
        :param CD: Drag coefficient [-].
        :type M: float.
        :type delta: float.
        :type CD: float.
        :returns: Aerodynamic drag [N].
        :rtype: float.
        """

        return 0.5 * delta * const.p_0 * const.Agamma * M * M * self.AC.S * CD

    def Thrust(self, delta, **kwargs):
        """This function computes the maximum thrust coefficient

        :param rating: Throttle setting {MCMB,MCRZ,LIDL}
        :param deltaT: direct throttle parameter [-]
        :param delta: Normalised pressure [-].
        :param theta: Normalised temperature [-].
        :param M: Mach number [-].
        :param DeltaTemp: Temperature deviation with respect to ISA [K].
        :type rating: string
        :type deltaT: float.
        :type delta: float.
        :type theta: float.
        :type M: float.
        :type DeltaTemp: float.
        :returns: Thrust [N].
        :rtype: float
        """

        CT = self.CT(delta=delta, **kwargs)

        return delta * self.AC.WREF * CT

    def ff(self, delta, theta, DeltaTemp, **kwargs):
        """This function computes the fuel flow

        :param rating: Throttle setting {MCMB,MCRZ,LIDL,TAXI}
        :param deltaT: direct throttle parameter [-]
        :param delta: Normalised pressure [-].
        :param theta: Normalised temperature [-].
        :param M: Mach number [-].
        :param DeltaTemp: Temperature deviation with respect to ISA [K].
        :type rating: string
        :type deltaT: float.
        :type delta: float.
        :type theta: float.
        :type M: float.
        :type DeltaTemp: float
        :returns: Fuel flow [kg s^-1].
        :rtype: float
        """

        if "rating" in kwargs:
            rating = checkArgument("rating", **kwargs)
            if rating == "TAXI":
                if self.AC.TFA is not None:
                    return self.AC.TFA / 60
                else:
                    return None

        CF = self.CF(delta=delta, theta=theta, DeltaTemp=DeltaTemp, **kwargs)

        if self.AC.BADAVersion == "4.2" or self.AC.BADAVersion == "DUMMY":
            return delta * pow(theta, 0.5) * self.AC.WREF * const.a_0 * CF / self.AC.LHV

        elif self.AC.BADAVersion == "4.3":
            return (
                pow(delta, self.AC.p_delta)
                * pow(theta, self.AC.p_theta)
                * self.AC.WREF
                * const.a_0
                * CF
                / self.AC.LHV
            )

    def ROCD(self, T, D, v, mass, ESF, h, DeltaTemp):
        """This function computes the Rate of Climb or Descent

        :param h: altitude [m].
        :param T: aircraft thrust [N].
        :param D: aircraft drag [N].
        :param v: aircraft true airspeed [TAS] [m s^-1].
        :param mass: actual aircraft mass  [kg].
        :param ESF: energy share factor [-].
        :param DeltaTemp: deviation with respect to ISA [K]
        :type h: float.
        :type T: float.
        :type D: float.
        :type v: float.
        :type mass: float.
        :type ESF: float.
        :type theta: float
        :type DeltaTemp: float.
        :returns: rate of climb/descend [m/s].
        :rtype: float
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        temp = theta * const.temp_0

        ROCD = ((temp - DeltaTemp) / temp) * (T - D) * v * ESF / (mass * const.g)

        return ROCD

    def controlLawThrust(self, ROCD, D, v, mass, ESF, h, DeltaTemp):
        """This function computes the Thrust based on the TEM control Law

        :param h: altitude [m].
        :param ROCD: rate of climb/descend [m/s].
        :param D: aircraft drag [N].
        :param v: aircraft true airspeed [TAS] [m s^-1].
        :param mass: actual aircraft mass  [kg].
        :param ESF: energy share factor [-].
        :param DeltaTemp: deviation with respect to ISA [K]
        :type h: float.
        :type ROCD: float.
        :type D: float.
        :type v: float.
        :type mass: float.
        :type ESF: float.
        :type theta: float
        :type DeltaTemp: float.
        :returns: thrust [N].
        :rtype: float
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        temp = theta * const.temp_0

        if ROCD == 0.0 or ESF == 0.0:
            thrust = (temp / (temp - DeltaTemp)) * (ROCD * mass * const.g) / v + D
        else:
            thrust = (temp / (temp - DeltaTemp)) * (ROCD * mass * const.g) / (
                ESF * v
            ) + D

        return thrust

    def propEfficiency(self, Wp, sigma, tas):
        """This function computes the propeller efficiency of a piston or turboprop engine through the momentum theory

        :param Wp: all-engine power [W]
        :param sigma: Normalised density [-].
        :param tas: true airspeed TAS [m/s].
        :type Wp: float
        :type sigma: float.
        :type tas: float.
        :returns: Propeller efficiency [-].
        :rtype: float
        """

        if self.AC.engineType == "TURBOPROP" or self.AC.engineType == "PISTON":
            a1 = (
                2
                * (Wp / self.AC.n_eng)
                / (
                    sigma
                    * const.rho_0
                    * self.AC.prop_dia
                    * self.AC.prop_dia
                    * pi
                    * tas
                    * tas
                    * tas
                    * self.AC.max_eff
                )
            )
            a2 = 0.0
            a3 = 1.0
            a4 = -1 * (self.AC.max_eff)

            coef = np.array([a1, a2, a3, a4])
            roots = np.roots(coef)

            for root in roots:
                if not np.iscomplex(root):
                    eff = float(np.real(root))
            return eff


class FlightEnvelope(BADA4):
    """This class is a BADA4 aircraft subclass and implements the flight envelope caclulations
    following the BADA4 manual.

    :param AC: parsed aircraft.
    :type AC: bada4.Parse.
    """

    def __init__(self, AC):
        super().__init__(AC)

    def maxM(self, LG):
        """This function computes the maximum M speed based on kinematic limitation (KLM) model

        :param LG: landing gear position, [LGUP/LGDN] [-]
        :type LG: string.
        :returns: maximum M speed [-].
        :rtype: float.
        """

        if LG == "LGUP":
            Mmax = self.AC.MMO
        else:
            if self.AC.MLE is not None:
                Mmax = self.AC.MLE
            else:
                Mmax = self.AC.MMO

        return Mmax

    def maxCAS(self, HLid, LG):
        """This function computes the maximum CAS speed based on kinematic limitation (KLM) model

        :param HLid: high lift devices position [-]
        :param LG: landing gear position, [LGUP/LGDN] [-]
        :type HLid: float.
        :type LG: string.
        :returns: maximum CAS speed [m s^-1].
        :rtype: float.
        """

        if HLid == 0 and LG == "LGUP":
            CASmax = conv.kt2ms(self.AC.VMO)

        elif HLid > 0 and LG == "LGUP":
            if self.AC.VFE[HLid] is not None:
                CASmax = conv.kt2ms(self.AC.VFE[HLid])
            else:
                CASmax = conv.kt2ms(self.AC.VMO)

        elif HLid == 0 and LG == "LGDN":
            if self.AC.VLE is not None:
                CASmax = conv.kt2ms(self.AC.VLE)
            else:
                CASmax = self.AC.VMO

        elif HLid > 0 and LG == "LGDN":
            if self.AC.VLE is not None:
                CASmax = conv.kt2ms(min(self.AC.VFE[HLid], self.AC.VLE))
            else:
                CASmax = self.AC.VMO

        return CASmax

    def VMax(self, h, HLid, LG, delta, theta, mass, nz=1.0):
        """This function computes the maximum speed

        :param HLid: high lift devices position [-]
        :param LG: landing gear position, [LGUP/LGDN] [-]
        :param theta: normalised air temperature [-].
        :param delta: normalised air pressure [-].
        :param mass: aircraft mass [kg].
        :param nz: load factor [-].
        :type HLid: float.
        :type LG: string.
        :type theta: float.
        :type delta: float.
        :type mass: float.
        :type nz: float.
        :returns: maximum CAS speed [m s^-1].
        :rtype: float.
        """

        if self.AC.MMO is not None:
            crossoverAlt = atm.crossOver(
                cas=self.maxCAS(HLid=HLid, LG=LG), Mach=self.maxM(LG=LG)
            )

            if h >= crossoverAlt:
                M = self.maxM(LG=LG)
                M_buffet = self.maxMbuffet(
                    HLid=HLid, LG=LG, delta=delta, mass=mass, nz=nz
                )

                if M_buffet is None:
                    return None

                sigma = atm.sigma(theta=theta, delta=delta)
                VMax = atm.mach2Cas(
                    Mach=min(M, M_buffet), theta=theta, delta=delta, sigma=sigma
                )

                # if M_buffet == float('-inf'):
                # VMax = float('-inf')
            else:
                VMax = self.maxCAS(HLid=HLid, LG=LG)
        else:
            VMax = self.maxCAS(HLid=HLid, LG=LG)

        return VMax

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

        [HLid, LG] = self.getAeroConfig(config=config)

        [theta, delta, sigma] = atm.atmosphereProperties(h=h, DeltaTemp=DeltaTemp)

        VmaxCertified = self.VMax(
            h=h, HLid=HLid, LG=LG, delta=delta, theta=theta, mass=mass, nz=1.0
        )
        VminCertified = self.VStall(
            theta=theta, delta=delta, mass=mass, HLid=HLid, LG=LG, nz=1.0
        )

        maxCASList = []
        for CAS in np.linspace(VminCertified, VmaxCertified, num=200, endpoint=True):
            [M, CAS, TAS] = atm.convertSpeed(
                v=conv.ms2kt(CAS),
                speedType="CAS",
                theta=theta,
                delta=delta,
                sigma=sigma,
            )

            maxThrust = self.Thrust(
                delta=delta, theta=theta, M=M, rating=rating, DeltaTemp=DeltaTemp
            )
            CL = self.CL(delta=delta, mass=mass, M=M, nz=1.0)
            CD = self.CD(HLid=HLid, LG=LG, CL=CL, M=M)
            Drag = self.D(delta=delta, M=M, CD=CD)

            if maxThrust >= Drag:
                maxCASList.append(CAS)

        if not maxCASList:
            return None
        else:
            return max(maxCASList)

    def maxMbuffet(self, HLid, LG, delta, mass, nz=1.0):
        """This function computes the maximum M speed applying buffet limitation

        :param HLid: high lift devices position [-]
        :param LG: landing gear position, [LGUP/LGDN] [-]
        :param delta: normalised air pressure [-].
        :param mass: aircraft mass [kg].
        :param nz: load factor [-].
        :type HLid: float.
        :type LG: string.
        :type delta: float.
        :type mass: float.
        :type nz: float.
        :returns: maximum M speed [-].
        :rtype: float.
        """

        # if CLMax model exist for aircraft, additional limitation apply
        if self.AC.CL_clean is not None:
            M_list = np.arange(0.01, self.AC.MMO + 0.001, 0.001)

            # start from maximum value, since we are looking for max M
            idx = -1
            M = M_list[idx]
            while True:
                CL = self.CL(delta=delta, mass=mass, M=M, nz=nz)
                CL_max = self.CLmax(M=M, HLid=HLid, LG=LG)

                if CL_max - CL > 0:
                    return M
                else:
                    # didn't find any M satisfying the CL < CLmax condition
                    if abs(idx) == len(M_list):
                        return None
                    else:
                        idx -= 1
                        M = M_list[idx]
        else:
            return self.AC.MMO

    def minMbuffet(self, HLid, LG, theta, delta, mass, nz=1.0):
        """This function computes the minimum M speed applying buffet limitation

        :param HLid: high lift devices position [-]
        :param LG: landing gear position, [LGUP/LGDN] [-]
        :param delta: normalised air pressure [-].
        :param mass: aircraft mass [kg].
        :param nz: load factor [-].
        :type HLid: float.
        :type LG: string.
        :type delta: float.
        :type mass: float.
        :type nz: float.
        :returns: maximum M speed [-].
        :rtype: float.
        """

        if HLid in self.AC.CL_max and LG in self.AC.CL_max[HLid]:
            CLmax = self.AC.CL_max[HLid][LG]
            # estimation of min M where CLmax = CL
            Mmin = sqrt(
                2
                * mass
                * const.g
                / (delta * const.p_0 * const.Agamma * CLmax * self.AC.S)
            )
            return Mmin

        else:
            if self.AC.MMO is not None:
                MMO_max = self.AC.MMO
            else:
                sigma = atm.sigma(theta=theta, delta=delta)
                MMO_max = atm.cas2Mach(
                    cas=conv.kt2ms(self.AC.VMO), theta=theta, delta=delta, sigma=sigma
                )

            M_list = np.arange(0.1, MMO_max + 0.001, 0.001)

            idx = 0
            M = M_list[idx]
            while True:
                CL = self.CL(delta=delta, mass=mass, M=M, nz=nz)
                CL_max = self.CLmax(M=M, HLid=HLid, LG=LG)

                if CL_max - CL > 0:
                    return M
                else:
                    # didn't find any M satisfying the CL < CLmax condition
                    if idx == len(M_list) - 1:
                        return None
                    else:
                        idx += 1
                        M = M_list[idx]

    def VMin(self, config, theta, delta, mass):
        """This function computes the minimum speed

        :param config: aircraft configuration [CR/IC/TO/AP/LD][-]
        :param delta: normalised air pressure [-].
        :param mass: aircraft mass [kg].
        :param mass: aircraft operating mass [kg]
        :type config: string.
        :type delta: float.
        :type theta: string.
        :type mass: float
        :returns: minimum CAS speed [m s^-1].
        :rtype: float
        """

        aeroConf = self.getAeroConfig(config=config)
        HLid = aeroConf[0]
        LG = aeroConf[1]

        if (HLid == 0 and LG == "LGUP") and self.AC.CL_clean is not None:
            Vmin = self.VStall(
                theta=theta, delta=delta, mass=mass, HLid=HLid, LG=LG, nz=1.2
            )
        else:
            if config == "TO":
                Vmin = self.AC.CVminTO * self.VStall(
                    theta=theta, delta=delta, mass=mass, HLid=HLid, LG=LG, nz=1.0
                )
            else:
                Vmin = self.AC.CVmin * self.VStall(
                    theta=theta, delta=delta, mass=mass, HLid=HLid, LG=LG, nz=1.0
                )

        return Vmin

    def VStall(self, mass, HLid, LG, nz=1.0, **kwargs):
        """This function calculates the stall speed based on the aerodynamic configuration

        :param HLid: high lift devices position [-]
        :param LG: landing gear position, [LGUP/LGDN] [-]
        :param delta: normalised air pressure [-].
        :param theta: normalised air temperature [-].
        :param mass: aircraft operating mass [kg]
        :param nz: load factor [-].
        :param h: altitude AMSL [m].
        :param DeltaTemp: deviation with respect to ISA [K]
        :type HLid: float.
        :type LG: string.
        :type delta: float.
        :type theta: float.
        :type mass: float.
        :type nz: float.
        :type h: float.
        :type DeltaTemp: float.
        :returns: stall CAS speed [m s^-1].
        :rtype: float
        """

        if "h" in kwargs:
            h = kwargs.get("h")
            DeltaTemp = checkArgument("DeltaTemp", **kwargs)
            delta = atm.delta(h, DeltaTemp)
            theta = atm.theta(h, DeltaTemp)
        else:
            theta = checkArgument("theta", **kwargs)
            delta = checkArgument("delta", **kwargs)

        sigma = atm.sigma(theta=theta, delta=delta)
        minM = self.minMbuffet(
            theta=theta, delta=delta, mass=mass, HLid=HLid, LG=LG, nz=nz
        )

        if minM is None:
            return None

        minCAS = atm.mach2Cas(Mach=minM, theta=theta, delta=delta, sigma=sigma)

        return minCAS

    def maxAltitude(self, HLid, LG, M, DeltaTemp, mass, nz=1.0):
        """This function computes the maximum altitude taking into account impact of buffet

        :param M: mach airspeed [-].
        :param mass: aircraft mass [kg].
        :param nz: load factor [-].
        :param HLid: high lift devices position [-]
        :param LG: landing gear position, [LGUP/LGDN] [-]
        :type M: float.
        :type mass: float.
        :type nz: float.
        :type HLid: float.
        :type LG: string.
        :returns: maximum altitude [m].
        :rtype: float
        """

        if HLid > 0:
            if self.AC.mfa is not None:
                hMax = self.AC.mfa
            else:
                hMax = self.AC.hmo
        else:
            hMax = self.AC.hmo

        if self.AC.CL_clean is not None:

            def f(H):
                delta = atm.delta(h=H[0], DeltaTemp=DeltaTemp)
                CL = self.CL(delta=delta, mass=mass, M=M, nz=nz)
                CL_max = self.CLmax(M=M, HLid=HLid, LG=LG)
                return -CL - CL_max

            hMax = float(
                fminbound(
                    f, x1=np.array([0]), x2=np.array([conv.ft2m(hMax)]), disp=False
                )
            )

        return hMax

    def getConfig(self, phase, h, mass, v, DeltaTemp=0.0, hRWY=0.0):
        """This function returns the aircraft aerodynamic configuration
        based on the aircraft altitude and speed and phase of flight

        :param hRWY: runway elevation AMSL [m].
        :param phase: aircraft phase of flight [Climb,Cruise,Descent][-].
        :param h: altitude [m].
        :param v: calibrated airspeed (CAS) [m s^-1].
        :param mass: aircraft mass [kg]
        :type hRWY: float.
        :type phase: string.
        :type h: float.
        :type v: float.
        :type mass: float.
        :returns: aircraft aerodynamic configuration [TO/IC/CR/AP/LD][-].
        :rtype: string
        """

        config = None

        [theta, delta, sigma] = atm.atmosphereProperties(h=h, DeltaTemp=DeltaTemp)

        # aircraft AGL altitude assuming being close to the RWY [m]
        h_AGL = h - hRWY

        HmaxTO_AGL = conv.ft2m(self.AC.HmaxPhase["TO"]) - hRWY
        HmaxIC_AGL = conv.ft2m(self.AC.HmaxPhase["IC"]) - hRWY
        HmaxAPP_AGL = conv.ft2m(self.AC.HmaxPhase["AP"]) - hRWY
        HmaxLD_AGL = conv.ft2m(self.AC.HmaxPhase["LD"]) - hRWY

        if phase == "Climb" and h_AGL <= HmaxTO_AGL:
            config = "TO"
            return config
        elif phase == "Climb" and (h_AGL > HmaxTO_AGL and h_AGL <= HmaxIC_AGL):
            config = "IC"
            return config
        elif (
            phase == "Cruise"
            or (phase == "Climb" and h_AGL >= HmaxIC_AGL)
            or (phase == "Descent" and h_AGL >= HmaxAPP_AGL)
        ):
            config = "CR"
            return config

        else:
            vMinCR = self.VMin(config="CR", mass=mass, theta=theta, delta=delta)
            vMinAPP = self.VMin(config="AP", mass=mass, theta=theta, delta=delta)

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
                and v < (vMinCR + conv.kt2ms(10))
            ) or (
                phase == "Descent"
                and (h_AGL + ep) < HmaxLD_AGL
                and (
                    (v + ep) < (vMinCR + conv.kt2ms(10))
                    and v >= (vMinAPP + conv.kt2ms(10))
                )
            ):
                config = "AP"

            elif (
                phase == "Descent"
                and (h_AGL + ep) < HmaxAPP_AGL
                and v >= (vMinCR + conv.kt2ms(10))
            ):
                config = "CR"

        if config is None:
            raise TypeError("Unable to determine aircraft configuration")

        return config

    def getAeroConfig(self, config):
        """This function returns the aircraft aerodynamic configuration
        based on the aerodynamic configuration ID, or None if not existent

        :param config: aircraft configuration [CR/IC/TO/AP/LD][-]
        :type config: string
        :returns: aircraft aerodynamic configuration combination of HLID and LG [-].
        :rtype: [float, string]
        """

        configDict = self.AC.aeroConfig.get(config)

        return [configDict["HLid"], configDict["LG"]]

    def getSpeedSchedule(self, phase):
        """This function returns the speed schedule
        based on the phase of flight {Climb, Cruise, Descent}

        :param phase: aircraft phase of flight {Climb, Cruise, Descent}
        :type phase: string
        :returns: speed schedule combination of CAS1, CAS2 and M [m s^-1, m s^-1, -].
        :rtype: [float, float, float]
        """

        speedScheduleDict = self.AC.speedSchedule[phase]

        return [
            conv.kt2ms(speedScheduleDict["CAS1"]),
            conv.kt2ms(speedScheduleDict["CAS2"]),
            speedScheduleDict["M"],
        ]

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


class ARPM(BADA4):
    """This class is a BADA4 aircraft subclass and implements the Airline Procedure Model (ARPM)
    following the BADA4 user manual.

    :param AC: parsed aircraft.
    :type AC: bada4.Parse.
    """

    def __init__(self, AC):
        super().__init__(AC)

        self.flightEnvelope = FlightEnvelope(AC)

    def climbSpeed(
        self,
        theta,
        delta,
        mass,
        h,
        hRWY=0.0,
        speedSchedule_default=None,
        procedure="BADA",
        config=None,
        NADP1_ALT=3000,
        NADP2_ALT=[1000, 3000],
        DeltaTemp=0.0,
    ):
        """This function computes the climb speed schedule CAS speed for any given altitude

        :param h: altitude [m].
        :param mass: aircraft mass [kg].
        :param theta: normalised air temperature [-].
        :param delta: normalised air pressure [-].
        :param speedSchedule_default: default speed schedule that will overwrite the BADA schedule [Vcl1, Vcl2, Mcl].
        :param hRWY: runway elevation AMSL [m].
        :param DeltaTemp: deviation with respect to ISA [K]
        :type h: float.
        :type mass: float.
        :type theta: float.
        :type delta: float.
        :type speedSchedule_default: [float, float, float].
        :type hRWY: float.
        :type DeltaTemp: float.
        :returns: climb calibrated airspeed (CAS) [m s^-1] & updatedSpeed if speed has been updated.
        :rtype: [float boolean]
        """

        # aircraft AGL altitude assuming being close to the RWY [m]
        h_AGL = h - hRWY

        phase = "Climb"
        acModel = self.AC.engineType

        [HLidTO, LGTO] = self.flightEnvelope.getAeroConfig(config="TO")
        VstallTO = self.flightEnvelope.VStall(
            h=h_AGL, mass=mass, HLid=HLidTO, LG=LGTO, nz=1.0, DeltaTemp=DeltaTemp
        )

        [HLidCR, LGCR] = self.flightEnvelope.getAeroConfig(config="CR")
        VstallCR = self.flightEnvelope.VStall(
            h=h_AGL, mass=mass, HLid=HLidCR, LG=LGCR, nz=1.0, DeltaTemp=DeltaTemp
        )
        [Vcl1, Vcl2, Mcl] = self.flightEnvelope.getSpeedSchedule(phase=phase)

        if speedSchedule_default is not None:
            Vcl1 = speedSchedule_default[0]
            Vcl2 = speedSchedule_default[1]
            Mcl = speedSchedule_default[2]

        crossOverAlt = atm.crossOver(cas=Vcl2, Mach=Mcl)

        if procedure == "BADA":
            if acModel == "JET":
                speed = list()
                speed.append(min(Vcl1, conv.kt2ms(250)))
                speed.append(self.AC.CVmin * VstallTO + conv.kt2ms(self.AC.V_cl[5]))
                speed.append(self.AC.CVmin * VstallTO + conv.kt2ms(self.AC.V_cl[4]))
                speed.append(self.AC.CVmin * VstallTO + conv.kt2ms(self.AC.V_cl[3]))
                speed.append(self.AC.CVmin * VstallTO + conv.kt2ms(self.AC.V_cl[2]))
                speed.append(self.AC.CVmin * VstallTO + conv.kt2ms(self.AC.V_cl[1]))

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
                    sigma = atm.sigma(theta=theta, delta=delta)
                    cas = atm.mach2Cas(Mach=Mcl, theta=theta, delta=delta, sigma=sigma)

            elif acModel == "TURBOPROP" or acModel == "PISTON":
                speed = list()
                speed.append(min(Vcl1, conv.kt2ms(250)))
                speed.append(self.AC.CVmin * VstallTO + conv.kt2ms(self.AC.V_cl[8]))
                speed.append(self.AC.CVmin * VstallTO + conv.kt2ms(self.AC.V_cl[7]))
                speed.append(self.AC.CVmin * VstallTO + conv.kt2ms(self.AC.V_cl[6]))

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
                    sigma = atm.sigma(theta=theta, delta=delta)
                    cas = atm.mach2Cas(Mach=Mcl, theta=theta, delta=delta, sigma=sigma)

        elif procedure == "NADP1":
            if acModel == "JET":
                speed = list()
                speed.append(min(Vcl1, conv.kt2ms(250)))
                speed.append(self.AC.CVminTO * VstallTO + conv.kt2ms(self.AC.V_cl[2]))

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
                speed.append(self.AC.CVminTO * VstallTO + conv.kt2ms(self.AC.V_cl[1]))

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
                speed.append(self.AC.CVmin * VstallCR + conv.kt2ms(self.AC.V_cl[2]))
                speed.append(self.AC.CVminTO * VstallTO + conv.kt2ms(self.AC.V_cl[2]))

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
                speed.append(self.AC.CVmin * VstallCR + conv.kt2ms(self.AC.V_cl[2]))
                speed.append(self.AC.CVminTO * VstallTO + conv.kt2ms(self.AC.V_cl[1]))

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

        # check if the speed is within the limits of minimum and maximum speed from the flight envelope, if not, overwrite calculated speed with flight envelope min/max speed
        if config is None:
            config = self.flightEnvelope.getConfig(
                h=h, phase=phase, v=cas, mass=mass, DeltaTemp=DeltaTemp, hRWY=hRWY
            )
        minSpeed = self.flightEnvelope.VMin(
            config=config, mass=mass, theta=theta, delta=delta
        )
        [HLid, LG] = self.flightEnvelope.getAeroConfig(config=config)
        maxSpeed = self.flightEnvelope.VMax(
            h=h, HLid=HLid, LG=LG, theta=theta, delta=delta, mass=mass, nz=1.2
        )

        eps = 1e-6  # float calculation precision
        # empty envelope - keep the original calculated CAS speed

        if minSpeed is None or maxSpeed is None:
            return [cas, "vV"]

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
        hRWY=0.0,
        speedSchedule_default=None,
        config=None,
        DeltaTemp=0.0,
    ):
        """This function computes the cruise speed schedule CAS speed for any given altitude

        :param h: altitude [m].
        :param mass: aircraft mass [kg].
        :param theta: normalised air temperature [-].
        :param delta: normalised air pressure [-].
        :param speedSchedule_default: default speed schedule that will overwrite the BADA schedule [Vcr1, Vcr2, Mcr].
        :type h: float.
        :type mass: float.
        :type theta: float.
        :type delta: float.
        :type speedSchedule_default: [float, float, float].
        :returns: climb calibrated airspeed (CAS) [m s^-1] & updatedSpeed if speed has been updated.
        :rtype: [float boolean]
        """

        phase = "Cruise"
        acModel = self.AC.engineType
        [Vcr1, Vcr2, Mcr] = self.flightEnvelope.getSpeedSchedule(phase=phase)

        if speedSchedule_default is not None:
            Vcr1 = speedSchedule_default[0]
            Vcr2 = speedSchedule_default[1]
            Mcr = speedSchedule_default[2]

        crossOverAlt = atm.crossOver(cas=Vcr2, Mach=Mcr)

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
                sigma = atm.sigma(theta=theta, delta=delta)
                cas = atm.mach2Cas(Mach=Mcr, theta=theta, delta=delta, sigma=sigma)

        elif acModel == "TURBOPROP" or acModel == "PISTON":
            if h < conv.ft2m(3000):
                cas = min(Vcr1, conv.kt2ms(150))
            elif h >= conv.ft2m(3000) and h < conv.ft2m(6000):
                cas = min(Vcr1, conv.kt2ms(180))
            elif h >= conv.ft2m(6000) and h < conv.ft2m(10000):
                cas = min(Vcr1, conv.kt2ms(250))
            elif h >= conv.ft2m(10000) and h < crossOverAlt:
                cas = Vcr2
            elif h >= crossOverAlt:
                sigma = atm.sigma(theta=theta, delta=delta)
                cas = atm.mach2Cas(Mach=Mcr, theta=theta, delta=delta, sigma=sigma)

        # check if the speed is within the limits of minimum and maximum speed from the flight envelope, if not, overwrite calculated speed with flight envelope min/max speed
        if config is None:
            config = self.flightEnvelope.getConfig(
                h=h, phase=phase, v=cas, mass=mass, DeltaTemp=DeltaTemp, hRWY=hRWY
            )

        minSpeed = self.flightEnvelope.VMin(
            config=config, mass=mass, theta=theta, delta=delta
        )
        [HLid, LG] = self.flightEnvelope.getAeroConfig(config=config)
        maxSpeed = self.flightEnvelope.VMax(
            h=h, HLid=HLid, LG=LG, theta=theta, delta=delta, mass=mass, nz=1.2
        )

        eps = 1e-6  # float calculation precision
        # empty envelope - keep the original calculated CAS speed

        if minSpeed is None or maxSpeed is None:
            return [cas, "vV"]

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
        hRWY=0.0,
        speedSchedule_default=None,
        config=None,
        DeltaTemp=0.0,
    ):
        """This function computes the descent speed schedule CAS speed for any given altitude

        :param h: altitude [m].
        :param mass: aircraft mass [kg].
        :param theta: normalised air temperature [-].
        :param delta: normalised air pressure [-].
        :param speedSchedule_default: default speed schedule that will overwrite the BADA schedule [Vdes1, Vdes2, Mdes].
        :param hRWY: runway elevation AMSL [m].
        :param DeltaTemp: deviation with respect to ISA [K]
        :type h: float.
        :type mass: float.
        :type theta: float.
        :type delta: float.
        :type speedSchedule_default: [float, float, float].
        :type hRWY: float.
        :type DeltaTemp: float.
        :returns: climb calibrated airspeed (CAS) [m s^-1] & updatedSpeed if speed has been updated.
        :rtype: [float, boolean]
        """

        # aircraft AGL altitude assuming being close to the RWY [m]
        h_AGL = h - hRWY

        phase = "Descent"
        acModel = self.AC.engineType

        [HLid, LG] = self.flightEnvelope.getAeroConfig(config="LD")
        VstallDES = self.flightEnvelope.VStall(
            h=h_AGL, mass=mass, HLid=HLid, LG=LG, nz=1.0, DeltaTemp=DeltaTemp
        )
        [Vdes1, Vdes2, Mdes] = self.flightEnvelope.getSpeedSchedule(phase=phase)

        if speedSchedule_default is not None:
            Vdes1 = speedSchedule_default[0]
            Vdes2 = speedSchedule_default[1]
            Mdes = speedSchedule_default[2]

        crossOverAlt = atm.crossOver(cas=Vdes2, Mach=Mdes)

        if acModel == "JET" or acModel == "TURBOPROP":
            speed = []
            speed.append(min(Vdes1, conv.kt2ms(220)))
            speed.append(self.AC.CVmin * VstallDES + conv.kt2ms(self.AC.V_des[4]))
            speed.append(self.AC.CVmin * VstallDES + conv.kt2ms(self.AC.V_des[3]))
            speed.append(self.AC.CVmin * VstallDES + conv.kt2ms(self.AC.V_des[2]))
            speed.append(self.AC.CVmin * VstallDES + conv.kt2ms(self.AC.V_des[1]))

            n = 1
            while n < len(speed):
                if speed[n] > speed[n - 1]:
                    speed[n] = speed[n - 1]
                n = n + 1

            epsilon = 1e-6

            if h < conv.ft2m(1000):
                cas = speed[4]
            elif h >= conv.ft2m(1000) and h < conv.ft2m(1500):
                cas = speed[3]
            elif h >= conv.ft2m(1500) and h < conv.ft2m(2000):
                cas = speed[2]
            elif h >= conv.ft2m(2000) and h < conv.ft2m(3000):
                cas = speed[1]
            elif h >= conv.ft2m(3000) and h < conv.ft2m(6000):
                cas = min(Vdes1, conv.kt2ms(220))
            elif h >= conv.ft2m(6000) and h < conv.ft2m(10000):
                cas = min(Vdes1, conv.kt2ms(250))
            elif h >= conv.ft2m(10000) and h < crossOverAlt:
                cas = Vdes2
            elif h >= crossOverAlt:
                sigma = atm.sigma(theta=theta, delta=delta)
                cas = atm.mach2Cas(Mach=Mdes, theta=theta, delta=delta, sigma=sigma)

        elif acModel == "PISTON":
            speed = list()
            speed.append(Vdes1)
            speed.append(self.AC.CVmin * VstallDES + conv.kt2ms(self.AC.V_des[7]))
            speed.append(self.AC.CVmin * VstallDES + conv.kt2ms(self.AC.V_des[6]))
            speed.append(self.AC.CVmin * VstallDES + conv.kt2ms(self.AC.V_des[5]))

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
                sigma = atm.sigma(theta=theta, delta=delta)
                cas = atm.mach2Cas(Mach=Mdes, theta=theta, delta=delta, sigma=sigma)

        # check if the speed is within the limits of minimum and maximum speed from the flight envelope, if not, overwrite calculated speed with flight envelope min/max speed
        if config is None:
            config = self.flightEnvelope.getConfig(
                h=h, phase=phase, v=cas, mass=mass, DeltaTemp=DeltaTemp, hRWY=hRWY
            )

        minSpeed = self.flightEnvelope.VMin(
            config=config, mass=mass, theta=theta, delta=delta
        )
        [HLid, LG] = self.flightEnvelope.getAeroConfig(config=config)
        maxSpeed = self.flightEnvelope.VMax(
            h=h, HLid=HLid, LG=LG, theta=theta, delta=delta, mass=mass, nz=1.2
        )

        eps = 1e-6  # float calculation precision
        # empty envelope - keep the original calculated CAS speed

        if minSpeed is None or maxSpeed is None:
            return [cas, "vV"]

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


class Optimization(BADA4):
    """This class implements the BADA4 optimization following the BADA4 manual.

    :param AC: parsed aircraft.
    :type AC: bada4.Parse.
    """

    def __init__(self, AC):
        super().__init__(AC)
        self.flightEnvelope = FlightEnvelope(AC)

    def CCI(self, theta, delta, cI):
        """This function computes the cost index coefficient for given flight conditions

        :param cI: cost index [kg min^-1].
        :param delta: normalised pressure [-].
        :param theta: normalised temperature [-].
        :type cI: float.
        :type delta: float.
        :type theta: float.
        :returns: cost index coefficient [-]
        :rtype: float.
        """

        if self.AC.BADAVersion == "4.2" or self.AC.BADAVersion == "DUMMY":
            return (
                (cI / 60.0)
                * self.AC.LHV
                / (self.AC.MTOW * delta * const.g * const.a_0 * sqrt(theta))
            )
        elif self.AC.BADAVersion == "4.3":
            return (
                (cI / 60.0)
                * self.AC.LHV
                / (
                    self.AC.MTOW
                    * pow(delta, self.AC.p_delta)
                    * const.g
                    * const.a_0
                    * pow(theta, self.AC.p_theta)
                )
            )

    def CW(self, mass, delta):
        """This function computes the weight coefficient at a given mass and pressure

        :param mass: aircraft mass [kg].
        :param delta: normalised pressure [-].
        :type mass: float.
        :type delta: float.
        :returns: weight coefficient
        """

        return mass * const.g / (self.AC.MTOW * delta * const.g)

    def SR(self, M, CF):
        """This function computes the specific range (SR) for given flight conditions

        :param M: mach ground speed [-].
        :param CF: fuel coefficient [-].
        :type M: float.
        :type CF: float.
        :returns: specific range [NM kg^-1]
        :rtype: float.
        """

        return M / CF

    def econMach(self, theta, delta, mass, DeltaTemp, cI, wS):
        """This function computes the economic mach
        for a given flight condition and cost index

        :param delta: normalised pressure [-].
        :param theta: normalised temperature [-].
        :param mass: aircraft weight [kg].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param cI: cost index [kg min^-1].
        :param wS: longitudinal wind speed [m s^-1].
        :type delta: float.
        :type theta: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type cI: float.
        :type wS: float.
        :returns: Maximum Range Cruise (MRC) in M [-]
        :rtype: float.
        """

        # clean configuration during CR
        HLid = 0
        LG = "LGUP"

        ccI = self.CCI(cI=cI, delta=delta, theta=theta)
        Mws = atm.tas2Mach(v=wS, theta=theta)

        # min/max M speed limitation
        Mmin = self.flightEnvelope.minMbuffet(
            theta=theta, delta=delta, mass=mass, HLid=HLid, LG=LG
        )
        Mmax = self.flightEnvelope.maxMbuffet(delta=delta, mass=mass, HLid=HLid, LG=LG)

        epsilon = 0.001
        M_list = np.arange(Mmin, Mmax + epsilon, epsilon)

        M_econ = []
        cost_econ = []
        for M in M_list:
            CL = self.CL(M=M, delta=delta, mass=mass)
            CD = self.CD(M=M, CL=CL, HLid=HLid, LG=LG)
            Drag = self.D(M=M, delta=delta, CD=CD)
            Thrust = Drag
            ThrustMax = self.Thrust(
                rating="MCRZ", delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
            )

            # max Thrust limitation
            if Thrust > ThrustMax:
                continue

            CT = self.CT(Thrust=Thrust, delta=delta)
            CF = self.CF(CT=CT, delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp)

            # maximize the cost function
            cost = self.SR(M=M + Mws, CF=ccI + CF)

            M_econ.append(M)
            cost_econ.append(cost)

        if not cost_econ:
            return float("Nan")

        econM = M_econ[cost_econ.index(max(cost_econ))]

        return proper_round(econM, 10)

        # def f(M):
        #     CL = self.CL(M=M[0], delta=delta, mass=mass)
        #     CD = self.CD(M=M[0], CL=CL, HLid=HLid, LG=LG)
        #     Drag = self.D(M=M[0], delta=delta, CD=CD)

        #     CT = self.CT(Thrust=Drag, delta=delta)
        #     CF = self.CF(CT=CT, delta=delta, theta=theta, M=M[0], DeltaTemp=DeltaTemp)

        # maximize the cost function -> to minimize, change the sign to -1 what was originally a maximization
        #     cost = - (self.SR(M=M[0]+Mws, CF=ccI+CF))
        #     return cost

        # bnds = Bounds([Mmin],[Mmax])
        # ThrustMax - Thrust >= 0
        # cons = ({'type': 'ineq','fun': lambda M: self.Thrust(rating='MCRZ', delta=delta, theta=theta, M=M[0], DeltaTemp=DeltaTemp) - self.D(M=M[0], delta=delta, CD=self.CD(M=M[0], CL=self.CL(M=M[0], delta=delta, mass=mass), HLid=HLid, LG=LG))})

        # econ = minimize(f, np.array([Mmin]), method='SLSQP', bounds=bnds, constraints=cons)
        # return float(econ.x)

    def MRC(self, theta, delta, mass, DeltaTemp, wS):
        """This function computes the Mach reperesenting Maximum Range Cruise (MRC) for given flight conditions

        :param delta: normalised pressure [-].
        :param theta: normalised temperature [-].
        :param mass: aircraft weight [kg].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param wS: longitudinal wind speed [m s^-1].
        :type delta: float.
        :type theta: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type wS: float.
        :returns: Maximum Range Cruise (MRC) in M [-]
        :rtype: float.
        """

        mrcM = self.econMach(
            cI=0.0, theta=theta, delta=delta, mass=mass, DeltaTemp=DeltaTemp, wS=wS
        )

        if isnan(mrcM):
            return float("Nan")

        return mrcM

    def LRC(self, theta, delta, mass, DeltaTemp, wS):
        """This function computes the Mach reperesenting Long Range Cruise (LRC) for given flight conditions

        :param delta: normalised pressure [-].
        :param theta: normalised temperature [-].
        :param mass: aircraft weight [kg].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param wS: longitudinal wind speed [m s^-1].
        :type delta: float.
        :type theta: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type wS: float.
        :returns: Long Range Cruise (LRC) in M [-]
        :rtype: float.
        """

        Mws = atm.tas2Mach(v=wS, theta=theta)

        MRC = self.MRC(theta=theta, delta=delta, mass=mass, DeltaTemp=DeltaTemp, wS=wS)

        if isnan(MRC):
            return float("Nan")

        # clean configuration during CR
        HLid = 0
        LG = "LGUP"

        CL = self.CL(M=MRC, delta=delta, mass=mass)
        CD = self.CD(M=MRC, CL=CL, HLid=HLid, LG=LG)
        Drag = self.D(M=MRC, delta=delta, CD=CD)
        CT = self.CT(Thrust=Drag, delta=delta)
        CF = self.CF(CT=CT, delta=delta, theta=theta, M=MRC, DeltaTemp=DeltaTemp)
        SR_LRC = 0.99 * self.SR(M=MRC + Mws, CF=CF)

        # min/max M speed limitation
        Mmax = self.flightEnvelope.maxMbuffet(delta=delta, mass=mass, HLid=HLid, LG=LG)

        # LRC > MRC
        epsilon = 0.001
        M_list = np.arange(MRC, Mmax + epsilon, epsilon)

        M_LRC = []
        cost_LRC = []

        for M in M_list:
            CL = self.CL(M=M, delta=delta, mass=mass)
            CL_max = self.CLmax(M=M, HLid=HLid, LG=LG)
            CD = self.CD(M=M, CL=CL, HLid=HLid, LG=LG)
            Drag = self.D(M=M, delta=delta, CD=CD)
            Thrust = Drag
            ThrustMax = self.Thrust(
                rating="MCRZ", delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
            )

            # max Thrust limitation
            if Thrust > ThrustMax:
                continue

            if CL > CL_max:
                continue

            CT = self.CT(Thrust=Thrust, delta=delta)
            CF = self.CF(CT=CT, delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp)

            # specific range for LRC (definition)
            SR = self.SR(M=M + Mws, CF=CF)
            # minimize the cost function
            cost_LRC.append(sqrt((SR - SR_LRC) ** 2))
            M_LRC.append(M)

        lrcM = M_LRC[cost_LRC.index(min(cost_LRC))]

        return lrcM

        # def f(M):
        #     CL = self.CL(delta=delta, mass=mass, M=M[0])
        #     CD = self.CD(M=M[0], CL=CL, HLid=HLid, LG=LG)
        #     Drag = self.D(M=M[0], delta=delta, CD=CD)

        #     CT = self.CT(Thrust=Drag, delta=delta)
        #     CF = self.CF(CT=CT, delta=delta, theta=theta, M=M[0], DeltaTemp=DeltaTemp)
        #     SR = self.SR(M=M[0]+Mws, CF=CF)

        #     return sqrt((SR - SR_LRC)**2)

        # bnds = Bounds([MRC],[Mmax])
        # ThrustMax - Thrust >= 0
        # cons = ({'type': 'ineq','fun': lambda M: self.Thrust(rating='MCRZ', delta=delta, theta=theta, M=M[0], DeltaTemp=DeltaTemp) - self.D(M=M[0], delta=delta, CD=self.CD(M=M[0], CL=self.CL(M=M[0], delta=delta, mass=mass), HLid=HLid, LG=LG))})

        # lrc = minimize(f, np.array([0.1]), method='SLSQP', bounds=bnds, constraints=cons)
        # return float(lrc.x)

        # minimum = float(fmin(f, x0=np.array([MRC]), disp=False))
        # return minimum

    def MEC(self, theta, delta, mass, DeltaTemp, wS):
        """This function computes the Mach reperesenting Maximum Endurance Cruise (MEC) for given flight conditions

        :param delta: normalised pressure [-].
        :param theta: normalised temperature [-].
        :param mass: aircraft weight [kg].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param wS: longitudinal wind speed [m s^-1].
        :type delta: float.
        :type theta: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type wS: float.
        :returns: Maximum Endurance Cruise (MEC) in M [-]
        :rtype: float.
        """

        # clean configuration during CR
        HLid = 0
        LG = "LGUP"

        Mws = atm.tas2Mach(v=wS, theta=theta)

        # min/max M speed limitation
        Mmin = self.flightEnvelope.minMbuffet(
            theta=theta, delta=delta, mass=mass, HLid=HLid, LG=LG
        )
        Mmax = self.flightEnvelope.maxMbuffet(delta=delta, mass=mass, HLid=HLid, LG=LG)

        epsilon = 0.001
        M_list = np.arange(Mmin, Mmax + epsilon, epsilon)

        M_mec = []
        CF_mec = []
        for M in M_list:
            CL = self.CL(M=M, delta=delta, mass=mass)
            CD = self.CD(M=M, CL=CL, HLid=HLid, LG=LG)
            Drag = self.D(M=M, delta=delta, CD=CD)
            Thrust = Drag
            ThrustMax = self.Thrust(
                rating="MCRZ", delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
            )

            # max Thrust limitation
            if Thrust > ThrustMax:
                continue

            CT = self.CT(Thrust=Thrust, delta=delta)
            CF = self.CF(CT=CT, delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp)

            # minimize the cost function
            CF_mec.append(CF)
            M_mec.append(M)

        if not CF_mec:
            return float("Nan")

        mecM = M_mec[CF_mec.index(min(CF_mec))]

        return proper_round(mecM, 10)

        # def f(M):
        #     CL = self.CL(M=M[0], delta=delta, mass=mass)
        #     CD = self.CD(M=M[0], CL=CL, HLid=HLid, LG=LG)
        #     Drag = self.D(M=M[0], delta=delta, CD=CD)

        #     CT = self.CT(Thrust=Drag, delta=delta)
        #     CF = self.CF(CT=CT, delta=delta, theta=theta, M=M[0], DeltaTemp=DeltaTemp)
        #     return CF

        # bnds = Bounds([Mmin],[Mmax + 1e-8])
        # ThrustMax - Thrust >= 0
        # cons = ({'type': 'ineq','fun': lambda M: self.Thrust(rating='MCRZ', delta=delta, theta=theta, M=M[0], DeltaTemp=DeltaTemp) - self.D(M=M[0], delta=delta, CD=self.CD(M=M[0], CL=self.CL(M=M[0], delta=delta, mass=mass), HLid=HLid, LG=LG))})

        # mecM = minimize(f, np.array([Mmin]), method='SLSQP', bounds=bnds, constraints=cons)
        # return float(mecM.x)

    def optAltitude(self, M, mass, DeltaTemp):
        """This function computes optimum altitude for a given flight condition

        :param delta: normalised pressure [-].
        :param theta: normalised temperature [-].
        :param mass: aircraft weight [kg].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param cI: cost index [kg min^-1].
        :param wS: longitudinal wind speed [m s^-1].
        :type delta: float.
        :type theta: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type cI: float.
        :type wS: float.
        :returns: optium altitude [ft]
        :rtype: float.
        """

        # clean configuration during CR
        HLid = 0
        LG = "LGUP"

        Hmax = self.flightEnvelope.maxAltitude(
            HLid=HLid, LG=LG, M=M, DeltaTemp=DeltaTemp, mass=mass, nz=1.0
        )

        epsilon = 100
        H_list = np.arange(0, Hmax, epsilon)

        H_opt = []
        cost_opt = []
        for H in H_list:
            theta = atm.theta(H, DeltaTemp=DeltaTemp)
            delta = atm.delta(H, DeltaTemp=DeltaTemp)

            # min/max M speed limitation
            Mmin = self.flightEnvelope.minMbuffet(
                theta=theta, delta=delta, mass=mass, HLid=HLid, LG=LG
            )
            Mmax = self.flightEnvelope.maxMbuffet(
                delta=delta, mass=mass, HLid=HLid, LG=LG
            )

            if M < Mmin or M > Mmax:
                continue

            CL = self.CL(M=M, delta=delta, mass=mass)
            CD = self.CD(M=M, CL=CL, HLid=HLid, LG=LG)

            Drag = self.D(M=M, delta=delta, CD=CD)
            Thrust = Drag
            ThrustMax = self.Thrust(
                rating="MCRZ", delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
            )

            # max Thrust limitation
            if Thrust > ThrustMax:
                continue

            CT = self.CT(Thrust=Thrust, delta=delta)
            CF = self.CF(CT=CT, delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp)
            ff = self.ff(CT=CT, delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp)
            a = atm.aSound(theta=theta)

            # maximize the cost function
            cost = (M * a) / ff
            # cost = CL/CD

            H_opt.append(H)
            cost_opt.append(cost)

        if not cost_opt:
            return float("Nan")

        optH = conv.m2ft(H_opt[cost_opt.index(max(cost_opt))])
        # bound the optimum altitude at 2000ft
        if optH < 2000.0:
            return 2000.0

        return proper_round(optH, 10)

        # def f(H):
        #     theta = atm.theta(h=H[0],DeltaTemp=DeltaTemp)
        #     delta = atm.delta(h=H[0],DeltaTemp=DeltaTemp)

        # min/max M speed limitation
        #     Mmin = self.flightEnvelope.minMbuffet(theta=theta, delta=delta, mass=mass, HLid=HLid, LG=LG)
        #     Mmax = self.flightEnvelope.maxMbuffet(delta=delta, mass=mass, HLid=HLid, LG=LG)

        # if M < Mmin or M > Mmax:
        #     return float('Inf')

        #     CL = self.CL(M=M, delta=delta, mass=mass)
        #     CD = self.CD(M=M, CL=CL, HLid=HLid, LG=LG)

        #     Drag = self.D(M=M, delta=delta, CD=CD)
        # ThrustMax = self.Thrust(rating='MCRZ', delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp)

        # max Thrust limitation
        # if Thrust > ThrustMax:
        #     return float('Inf')

        #     CT = self.CT(Thrust=Drag, delta=delta)
        #     ff = self.ff(CT=CT, delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp)
        #     a = atm.aSound(theta=theta)

        # maximize the cost function
        #     return -((M*a) / ff)

        # optAlt = conv.m2ft(float(fminbound(f, x1=np.array([0]), x2=np.array([Hmax]), disp=False)))

        # bound the optimum altitude at 2000ft
        # print(optAlt)
        # if optAlt < 1e-3:
        #     return float('Nan')
        # if optAlt < 2000:
        #     return 2000

        # return optAlt

    def getOPTParam(self, optParam, var_1, var_2=None):
        """This function returns value of the OPT parameter based on the input value from OPT file
        like LRC, MEC, MRC, ECON speed or OPTALT altitude

        .. note::
                array used in this fuction is expected to be sorted (design of OPT files)

        :param optParam: name of optimization file {LRC,MEC,MRC,ECON,OPTALT}.
        :param var_1: value of the first optimizing factor.
        :param var_2: (optional) value of the second optimizing factor.
        :type optParam: string.
        :type var_1: float.
        :type var_2: float.
        :returns: None
        """

        filename = os.path.join(
            self.AC.filePath,
            "BADA4",
            self.AC.BADAVersion,
            self.AC.acName,
            optParam + ".OPT",
        )

        def findNearest(value, array):
            """This function returns indices of the nearest value in the array.
            if the value is lower/higher than lowest/highest value in array, only one idx is returned
            otherwise if the value is somewhere in between, 2 closest (left and right) idx are returned

            .. note::
                    array used in this fuction is expected to be sorted (design of OPT files)

            :param value: value to which the array value will be comapred
            :param array: list of values

            :type value: float.
            :type array: array of float.
            """

            nearestIdx = list()

            idx = np.searchsorted(array, value, side="left")

            if idx == len(array):
                nearestIdx = idx - 1

            elif idx == 0 or value == array[idx]:
                nearestIdx = idx

            elif value < array[idx] or value > array[idx]:
                nearestIdx = [idx - 1, idx]

            return nearestIdx

        def parseOPT(filename):
            """This function parses BADA4 OPT ascii formatted files

            :param filename: path to the ___.OPT ascii formatted file.
            :type filename: str.
            """

            file = open(filename, "r")
            lines = file.readlines()

            self.tableTypes = lines[8].split(":")[1].strip()
            self.tableDimension = lines[10].split(":")[1].strip()

            self.var_1 = list()
            self.var_2 = list()
            self.var_3 = list()

            if self.tableTypes == "2D":
                self.tableDimensionColumns = int(self.tableDimension.split("x")[0])
                self.tableDimensionRows = int(self.tableDimension.split("x")[1])

                self.var_2 = [
                    float(i)
                    for i in list(
                        filter(None, lines[13].split("|")[1].strip().split(" "))
                    )
                ]

                for j in range(16, 16 + self.tableDimensionRows, 1):
                    self.var_1.append(float(lines[j].split("|")[0].strip()))
                    self.var_3.extend(
                        [
                            float(i) if i != "---" else float("nan")
                            for i in list(
                                filter(None, lines[j].split("|")[1].strip().split(" "))
                            )
                        ]
                    )

            if self.tableTypes == "1D":
                self.tableDimensionColumns = int(self.tableDimension.split("x")[1])
                self.tableDimensionRows = int(self.tableDimension.split("x")[0])

                for j in range(15, 15 + self.tableDimensionRows, 1):
                    self.var_1.append(float(lines[j].split("|")[0].strip()))
                    self.var_2.append(float(lines[j].split("|")[1].strip()))

        parseOPT(filename=filename)

        if self.tableTypes == "2D":
            if var_2 is None:
                return float("NaN")

            nearestIdx_1 = np.array(findNearest(var_1, self.var_1))
            nearestIdx_2 = np.array(findNearest(var_2, self.var_2))

            # if nearestIdx_1 & nearestIdx_2 [1] [1]
            if (nearestIdx_1.size == 1) & (nearestIdx_2.size == 1):
                return self.var_3[
                    nearestIdx_1 * (self.tableDimensionColumns) + nearestIdx_2
                ]

            # if nearestIdx_1 & nearestIdx_2 [1] [1,2]
            if (nearestIdx_1.size == 1) & (nearestIdx_2.size == 2):
                varTemp_1 = self.var_3[
                    nearestIdx_1 * (self.tableDimensionColumns) + nearestIdx_2[0]
                ]
                varTemp_2 = self.var_3[
                    nearestIdx_1 * (self.tableDimensionColumns) + nearestIdx_2[1]
                ]

                # interpolation between the 2 found points
                interpVar = np.interp(
                    var_2,
                    [self.var_2[nearestIdx_2[0]], self.var_2[nearestIdx_2[1]]],
                    [varTemp_1, varTemp_2],
                )
                return interpVar

            # if nearestIdx_1 & nearestIdx_2 [1,2] [1]
            if (nearestIdx_1.size == 2) & (nearestIdx_2.size == 1):
                varTemp_1 = self.var_3[
                    nearestIdx_1[0] * (self.tableDimensionColumns) + nearestIdx_2
                ]
                varTemp_2 = self.var_3[
                    nearestIdx_1[1] * (self.tableDimensionColumns) + nearestIdx_2
                ]

                # interpolation between the 2 found points
                interpVar = np.interp(
                    var_1,
                    [self.var_1[nearestIdx_1[0]], self.var_1[nearestIdx_1[1]]],
                    [varTemp_1, varTemp_2],
                )
                return interpVar

            # if nearestIdx_1 & nearestIdx_2 [1,2] [1,2]
            if (nearestIdx_1.size == 2) & (nearestIdx_2.size == 2):
                varTemp_1 = self.var_3[
                    nearestIdx_1[0] * (self.tableDimensionColumns) + nearestIdx_2[0]
                ]
                varTemp_2 = self.var_3[
                    nearestIdx_1[0] * (self.tableDimensionColumns) + nearestIdx_2[1]
                ]

                varTemp_3 = self.var_3[
                    nearestIdx_1[1] * (self.tableDimensionColumns) + nearestIdx_2[0]
                ]
                varTemp_4 = self.var_3[
                    nearestIdx_1[1] * (self.tableDimensionColumns) + nearestIdx_2[1]
                ]

                # interpolation between the 4 found points
                interpVar_1 = np.interp(
                    var_2,
                    [self.var_2[nearestIdx_2[0]], self.var_2[nearestIdx_2[1]]],
                    [varTemp_1, varTemp_2],
                )
                interpVar_2 = np.interp(
                    var_2,
                    [self.var_2[nearestIdx_2[0]], self.var_2[nearestIdx_2[1]]],
                    [varTemp_3, varTemp_4],
                )
                interpVar_3 = np.interp(
                    var_1,
                    [self.var_1[nearestIdx_1[0]], self.var_1[nearestIdx_1[1]]],
                    [interpVar_1, interpVar_2],
                )

                return interpVar_3

        if self.tableTypes == "1D":
            nearestIdx_1 = np.array(findNearest(var_1, self.var_1))
            # if nearestIdx_1 & nearestIdx_2 [1] [1]
            if nearestIdx_1.size == 1:
                return self.var_2[nearestIdx_1]

            if nearestIdx_1.size == 2:
                varTemp_1 = self.var_2[nearestIdx_1[0]]
                varTemp_2 = self.var_2[nearestIdx_1[1]]

                interpVar = np.interp(
                    var_1,
                    [self.var_1[nearestIdx_1[0]], self.var_1[nearestIdx_1[1]]],
                    [varTemp_1, varTemp_2],
                )
                return interpVar


class PTD(BADA4):
    """This class implements the PTD file creator for BADA4 aircraft following BADA4 manual.

    :param AC: parsed aircraft.
    :type AC: bada4.Parse.
    """

    def __init__(self, AC):
        super().__init__(AC)

        self.flightEnvelope = FlightEnvelope(AC)
        self.ARPM = ARPM(AC)

    def create(self, DeltaTemp, saveToPath):
        """This function creates the BADA4 PTD file

        :param saveToPath: path to directory where PTF should be stored [-]
        :param DeltaTemp: deviation from ISA temperature [K]
        :type saveToPath: string.
        :type DeltaTemp: float.
        :returns: NONE
        """

        # 3 different mass levels [kg]
        massList = [
            1.2 * self.AC.OEW,
            self.AC.OEW + 0.7 * (self.AC.MTOW - self.AC.OEW),
            self.AC.MTOW,
        ]
        max_alt_ft = self.AC.hmo

        if max_alt_ft > 51000:
            max_alt_ft = 51000

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

        CRList = []
        CLList = []
        DESList = []

        for mass in massList:
            CLList.append(
                self.PTD_climb(
                    mass=mass, altitudeList=altitudeList, DeltaTemp=DeltaTemp
                )
            )
            CRList.append(
                self.PTD_cruise(
                    mass=mass, altitudeList=altitudeList, DeltaTemp=DeltaTemp
                )
            )
            DESList.append(
                self.PTD_descent(
                    mass=mass, altitudeList=altitudeList, DeltaTemp=DeltaTemp
                )
            )

        self.save2PTD(
            saveToPath=saveToPath,
            CLList=CLList,
            CRList=CRList,
            DESList=DESList,
            DeltaTemp=DeltaTemp,
        )

    def save2PTD(self, saveToPath, CLList, CRList, DESList, DeltaTemp):
        """This function saves data to PTD file

        :param saveToPath: path to directory where PTD should be stored [-]
        :param CLList: list of PTD data in CLIMB [-].
        :param CRList: list of PTD data in CRUISE [-].
        :param DESList_med: list of PTD data in DESCENT [-].
        :param DeltaTemp: deviation from ISA temperature [K]
        :type saveToPath: string.
        :type CLList: list.
        :type CRList: list.
        :type DESList: list.
        :type DeltaTemp: float.
        :returns: NONE
        """

        def Nan2Zero(list):
            # replace NAN values by 0 for printing purposes
            for n in range(len(list)):
                for k in range(len(list[n])):
                    for m in range(len(list[n][k])):
                        if isinstance(list[n][k][m], float):
                            if isnan(list[n][k][m]):
                                list[n][k][m] = 0
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

        filename = saveToPath + self.AC.acName + "_ISA" + ISA + ".PTD"

        file = open(filename, "w")
        file.write("BADA PERFORMANCE FILE RESULTS\n")
        file = open(filename, "a")
        file.write("=============================\n=============================\n\n")
        file.write("Low mass CLIMB\n")
        file.write("==============\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass   Thrust    Drag     Fuel    ESF    ROCD   gamma Conf  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [N]     [N]     [kgm]    [-]   [fpm]   [deg]  [-]     \n"
        )

        # replace NAN values by 0 for printing purposes
        CLList = Nan2Zero(CLList)
        CRList = Nan2Zero(CRList)
        DESList = Nan2Zero(DESList)

        # low mass
        list_mass = CLList[0]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %7.0f %7.0f %9.2f %6.3f %6.0f %7.2f  %s   %s\n"
                % (
                    list_mass[0][k],
                    list_mass[1][k],
                    list_mass[2][k],
                    list_mass[3][k],
                    list_mass[4][k],
                    list_mass[5][k],
                    list_mass[6][k],
                    list_mass[7][k],
                    list_mass[8][k],
                    list_mass[9][k],
                    list_mass[10][k],
                    list_mass[11][k],
                    list_mass[12][k],
                    list_mass[13][k],
                    list_mass[14][k],
                    list_mass[15][k],
                    list_mass[16][k],
                )
            )

        file.write("\n\nMedium mass CLIMB\n")
        file.write("=================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass   Thrust    Drag     Fuel    ESF    ROCD   gamma Conf  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [N]     [N]     [kgm]    [-]   [fpm]   [deg]  [-]     \n"
        )

        # medium mass
        list_mass = CLList[1]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %7.0f %7.0f %9.2f %6.3f %6.0f %7.2f  %s   %s\n"
                % (
                    list_mass[0][k],
                    list_mass[1][k],
                    list_mass[2][k],
                    list_mass[3][k],
                    list_mass[4][k],
                    list_mass[5][k],
                    list_mass[6][k],
                    list_mass[7][k],
                    list_mass[8][k],
                    list_mass[9][k],
                    list_mass[10][k],
                    list_mass[11][k],
                    list_mass[12][k],
                    list_mass[13][k],
                    list_mass[14][k],
                    list_mass[15][k],
                    list_mass[16][k],
                )
            )

        file.write("\n\nHigh mass CLIMB\n")
        file.write("===============\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass   Thrust    Drag     Fuel    ESF    ROCD   gamma Conf  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [N]     [N]     [kgm]    [-]   [fpm]   [deg]  [-]     \n"
        )

        # high mass
        list_mass = CLList[2]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %7.0f %7.0f %9.2f %6.3f %6.0f %7.2f  %s   %s\n"
                % (
                    list_mass[0][k],
                    list_mass[1][k],
                    list_mass[2][k],
                    list_mass[3][k],
                    list_mass[4][k],
                    list_mass[5][k],
                    list_mass[6][k],
                    list_mass[7][k],
                    list_mass[8][k],
                    list_mass[9][k],
                    list_mass[10][k],
                    list_mass[11][k],
                    list_mass[12][k],
                    list_mass[13][k],
                    list_mass[14][k],
                    list_mass[15][k],
                    list_mass[16][k],
                )
            )

        file.write("\n\nLow mass DESCENT\n")
        file.write("================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass   Thrust    Drag     Fuel    ESF    ROCD   gamma Conf  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [N]     [N]     [kgm]    [-]   [fpm]   [deg]  [-]     \n"
        )

        # low mass
        list_mass = DESList[0]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %7.0f %7.0f %9.2f %6.3f %6.0f %7.2f  %s   %s\n"
                % (
                    list_mass[0][k],
                    list_mass[1][k],
                    list_mass[2][k],
                    list_mass[3][k],
                    list_mass[4][k],
                    list_mass[5][k],
                    list_mass[6][k],
                    list_mass[7][k],
                    list_mass[8][k],
                    list_mass[9][k],
                    list_mass[10][k],
                    list_mass[11][k],
                    list_mass[12][k],
                    list_mass[13][k],
                    list_mass[14][k],
                    list_mass[15][k],
                    list_mass[16][k],
                )
            )

        file.write("\n\nMedium mass DESCENT\n")
        file.write("===================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass   Thrust    Drag     Fuel    ESF    ROCD   gamma Conf  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [N]     [N]     [kgm]    [-]   [fpm]   [deg]  [-]     \n"
        )

        # medium mass
        list_mass = DESList[1]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %7.0f %7.0f %9.2f %6.3f %6.0f %7.2f  %s   %s\n"
                % (
                    list_mass[0][k],
                    list_mass[1][k],
                    list_mass[2][k],
                    list_mass[3][k],
                    list_mass[4][k],
                    list_mass[5][k],
                    list_mass[6][k],
                    list_mass[7][k],
                    list_mass[8][k],
                    list_mass[9][k],
                    list_mass[10][k],
                    list_mass[11][k],
                    list_mass[12][k],
                    list_mass[13][k],
                    list_mass[14][k],
                    list_mass[15][k],
                    list_mass[16][k],
                )
            )

        file.write("\n\nHigh mass DESCENT\n")
        file.write("=================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass   Thrust    Drag     Fuel    ESF    ROCD   gamma Conf  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [N]     [N]     [kgm]    [-]   [fpm]   [deg]  [-]     \n"
        )

        # high mass
        list_mass = DESList[2]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %7.0f %7.0f %9.2f %6.3f %6.0f %7.2f  %s   %s\n"
                % (
                    list_mass[0][k],
                    list_mass[1][k],
                    list_mass[2][k],
                    list_mass[3][k],
                    list_mass[4][k],
                    list_mass[5][k],
                    list_mass[6][k],
                    list_mass[7][k],
                    list_mass[8][k],
                    list_mass[9][k],
                    list_mass[10][k],
                    list_mass[11][k],
                    list_mass[12][k],
                    list_mass[13][k],
                    list_mass[14][k],
                    list_mass[15][k],
                    list_mass[16][k],
                )
            )

        file.write("\n\nLow mass CRUISE\n")
        file.write("===============\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass   Thrust    Drag     Fuel    ESF    ROCD   gamma Conf  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [N]     [N]     [kgm]    [-]   [fpm]   [deg]  [-]     \n"
        )

        # low mass
        list_mass = CRList[0]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %7.0f %7.0f %9.2f %6.3f %6.0f %7.2f  %s   %s\n"
                % (
                    list_mass[0][k],
                    list_mass[1][k],
                    list_mass[2][k],
                    list_mass[3][k],
                    list_mass[4][k],
                    list_mass[5][k],
                    list_mass[6][k],
                    list_mass[7][k],
                    list_mass[8][k],
                    list_mass[9][k],
                    list_mass[10][k],
                    list_mass[11][k],
                    list_mass[12][k],
                    list_mass[13][k],
                    list_mass[14][k],
                    list_mass[15][k],
                    list_mass[16][k],
                )
            )

        file.write("\n\nMedium mass CRUISE\n")
        file.write("==================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass   Thrust    Drag     Fuel    ESF    ROCD   gamma Conf  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [N]     [N]     [kgm]    [-]   [fpm]   [deg]  [-]     \n"
        )

        # medium mass
        list_mass = CRList[1]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %7.0f %7.0f %9.2f %6.3f %6.0f %7.2f  %s   %s\n"
                % (
                    list_mass[0][k],
                    list_mass[1][k],
                    list_mass[2][k],
                    list_mass[3][k],
                    list_mass[4][k],
                    list_mass[5][k],
                    list_mass[6][k],
                    list_mass[7][k],
                    list_mass[8][k],
                    list_mass[9][k],
                    list_mass[10][k],
                    list_mass[11][k],
                    list_mass[12][k],
                    list_mass[13][k],
                    list_mass[14][k],
                    list_mass[15][k],
                    list_mass[16][k],
                )
            )

        file.write("\n\nHigh mass CRUISE\n")
        file.write("================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass   Thrust    Drag     Fuel    ESF    ROCD   gamma Conf  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [N]     [N]     [kgm]    [-]   [fpm]   [deg]  [-]     \n"
        )

        # high mass
        list_mass = CRList[2]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %7.0f %7.0f %9.2f %6.3f %6.0f %7.2f  %s   %s\n"
                % (
                    list_mass[0][k],
                    list_mass[1][k],
                    list_mass[2][k],
                    list_mass[3][k],
                    list_mass[4][k],
                    list_mass[5][k],
                    list_mass[6][k],
                    list_mass[7][k],
                    list_mass[8][k],
                    list_mass[9][k],
                    list_mass[10][k],
                    list_mass[11][k],
                    list_mass[12][k],
                    list_mass[13][k],
                    list_mass[14][k],
                    list_mass[15][k],
                    list_mass[16][k],
                )
            )

    def PTD_climb(self, mass, altitudeList, DeltaTemp):
        """This function calculates the BADA4 PTD data in CLIMB

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
        gamma_complet = []
        conf_complet = []
        Lim_complet = []

        phase = "Climb"

        [Vcl1, Vcl2, Mcl] = self.flightEnvelope.getSpeedSchedule(phase=phase)
        Vcl1 = min(Vcl1, conv.kt2ms(250))
        crossAlt = atm.crossOver(cas=Vcl2, Mach=Mcl)

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
            [cas, speedUpdated] = self.ARPM.climbSpeed(
                h=H_m,
                mass=mass,
                theta=theta,
                delta=delta,
                DeltaTemp=DeltaTemp,
                speedSchedule_default=[Vcl1, Vcl2, Mcl],
            )
            tas = atm.cas2Tas(cas=cas, delta=delta, sigma=sigma)
            M = atm.tas2Mach(v=tas, theta=theta)
            a = atm.aSound(theta=theta)
            FL = h / 100

            # add limitation that has been applied (if some has been applied)
            limitation = speedUpdated

            ff = (
                self.ff(
                    rating="MCMB", delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
                )
                * 60
            )
            Thrust = self.Thrust(
                rating="MCMB", delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
            )
            config = self.flightEnvelope.getConfig(
                h=H_m,
                phase=phase,
                v=cas,
                mass=mass,
                DeltaTemp=DeltaTemp,
            )

            # ensure the continuity of aerodyamic configuration during Climb phase of flight
            if conf_complet:
                prevConf = conf_complet[-1]

                if config == "TO" and (prevConf == "IC" or prevConf == "CR"):
                    config = prevConf
                elif config == "IC" and prevConf == "CR":
                    config = prevConf

            [HLid, LG] = self.flightEnvelope.getAeroConfig(config=config)
            CL = self.CL(M=M, delta=delta, mass=mass)
            CD = self.CD(M=M, CL=CL, HLid=HLid, LG=LG)
            Drag = self.D(M=M, delta=delta, CD=CD)

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

            temp_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)
            dhdt = (conv.ft2m(ROCD / 60)) * temp_const
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
            ROCD_complet.append(ROCD)
            gamma_complet.append(gamma)
            conf_complet.append(config)
            Lim_complet.append(limitation)

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
            gamma_complet,
            conf_complet,
            Lim_complet,
        ]

        return CLList

    def PTD_descent(self, mass, altitudeList, DeltaTemp):
        """This function calculates the BADA4 PTD data in DESCENT

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
        ff_complet = []
        ESF_complet = []
        ROCD_complet = []
        gamma_complet = []
        conf_complet = []
        Lim_complet = []

        phase = "Descent"

        [Vdes1, Vdes2, Mdes] = self.flightEnvelope.getSpeedSchedule(phase=phase)
        Vdes1 = min(Vdes1, conv.kt2ms(250))
        crossAlt = atm.crossOver(cas=Vdes2, Mach=Mdes)

        for h in reversed(altitudeList):
            H_m = conv.ft2m(h)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
            [cas, speedUpdated] = self.ARPM.descentSpeed(
                h=H_m,
                mass=mass,
                theta=theta,
                delta=delta,
                DeltaTemp=DeltaTemp,
                speedSchedule_default=[Vdes1, Vdes2, Mdes],
            )
            tas = atm.cas2Tas(cas=cas, delta=delta, sigma=sigma)
            M = atm.tas2Mach(v=tas, theta=theta)
            a = atm.aSound(theta=theta)
            FL = h / 100

            # add limitation that has been applied (if some has been applied)
            limitation = speedUpdated

            Thrust = self.Thrust(
                rating="LIDL", delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
            )
            config = self.flightEnvelope.getConfig(
                h=H_m, phase="Descent", v=cas, mass=mass, DeltaTemp=DeltaTemp
            )

            # ensure the continuity of aerodyamic configuration during Descent phase of flight
            if conf_complet:
                prevConf = conf_complet[0]

                if config == "CR" and (prevConf == "AP" or prevConf == "LD"):
                    config = prevConf
                elif config == "AP" and prevConf == "LD":
                    config = prevConf

            [HLid, LG] = self.flightEnvelope.getAeroConfig(config=config)
            CL = self.CL(M=M, delta=delta, mass=mass)
            CD = self.CD(M=M, CL=CL, HLid=HLid, LG=LG)
            Drag = self.D(M=M, delta=delta, CD=CD)

            ff = (
                self.ff(
                    rating="LIDL", delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
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
            temp_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)
            dhdt = (conv.ft2m(ROCD / 60)) * temp_const
            gamma = conv.rad2deg(asin(dhdt / tas))

            minSpeed = self.flightEnvelope.VMin(
                config=config, mass=mass, theta=theta, delta=delta
            )
            [HLid, LG] = self.flightEnvelope.getAeroConfig(config=config)
            maxSpeed = self.flightEnvelope.VMax(
                h=h, HLid=HLid, LG=LG, theta=theta, delta=delta, mass=mass
            )

            # in case of AP & LD thrust is computed to fly a 3deg slope
            if config == "AP" or config == "LD":
                gamma = -3.0
                temp_const = (theta * const.temp_0) / (theta * const.temp_0 - DeltaTemp)

                ROCD_gamma = sin(conv.deg2rad(gamma)) * tas * (1 / temp_const)
                ROCD = conv.m2ft(ROCD_gamma) * 60  # [ft/min]

                n = 1.0  # aircraft.loadFactor(gamma) - use this in case of L = W * (cos(gamma))
                CL = self.CL(M=M, delta=delta, mass=mass, nz=n)
                CD = self.CD(M=M, CL=CL, HLid=HLid, LG=LG)
                Drag = self.D(M=M, delta=delta, CD=CD)
                Thrust = (ROCD_gamma * mass * const.g) * temp_const / (ESF * tas) + Drag
                CT = self.CT(Thrust=Thrust, delta=delta)
                ff = (
                    self.ff(CT=CT, delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp)
                    * 60
                )

            FL_complet = [proper_round(FL)] + FL_complet
            T_complet = [theta * const.temp_0] + T_complet
            p_complet = [delta * const.p_0] + p_complet
            rho_complet = [sigma * const.rho_0] + rho_complet
            a_complet = [a] + a_complet
            TAS_complet = [conv.ms2kt(tas)] + TAS_complet
            CAS_complet = [conv.ms2kt(cas)] + CAS_complet
            M_complet = [M] + M_complet
            mass_complet = [mass] + mass_complet
            Thrust_complet = [Thrust] + Thrust_complet
            Drag_complet = [Drag] + Drag_complet
            ff_complet = [ff] + ff_complet
            ESF_complet = [ESF] + ESF_complet
            ROCD_complet = [-1 * ROCD] + ROCD_complet
            gamma_complet = [gamma] + gamma_complet
            conf_complet = [config] + conf_complet
            Lim_complet = [limitation] + Lim_complet

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
            ff_complet,
            ESF_complet,
            ROCD_complet,
            gamma_complet,
            conf_complet,
            Lim_complet,
        ]

        return DESList

    def PTD_cruise(self, mass, altitudeList, DeltaTemp):
        """This function calculates the BADA4 PTD data in CRUISE

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
        gamma_complet = []
        conf_complet = []
        Lim_complet = []

        phase = "Cruise"

        [Vcr1, Vcr2, Mcr] = self.flightEnvelope.getSpeedSchedule(phase=phase)
        Vcr1 = min(Vcr1, conv.kt2ms(250))

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
            [cas, speedUpdated] = self.ARPM.cruiseSpeed(
                h=H_m,
                mass=mass,
                theta=theta,
                delta=delta,
                speedSchedule_default=[Vcr1, Vcr2, Mcr],
            )
            tas = atm.cas2Tas(cas=cas, delta=delta, sigma=sigma)
            M = atm.tas2Mach(v=tas, theta=theta)
            a = atm.aSound(theta=theta)
            FL = h / 100

            # add limitation that has been applied (if some has been applied)
            limitation = ""

            config = "CR"
            HLid = 0
            LG = "LGUP"
            CL = self.CL(M=M, delta=delta, mass=mass)
            CD = self.CD(M=M, CL=CL, HLid=HLid, LG=LG)
            Drag = self.D(M=M, delta=delta, CD=CD)
            Thrust = Drag
            ThrustMax = self.Thrust(
                rating="MCRZ", delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
            )
            CT = self.CT(Thrust=Thrust, delta=delta)
            ff = self.ff(CT=CT, delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp) * 60

            if Thrust > ThrustMax:
                # "(T)" - as thrust limited
                limitation += "T"

            limitation += speedUpdated

            ESF = 0.0
            ROCD = 0.0
            gamma = 0.0

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
            gamma_complet.append(gamma)
            conf_complet.append(config)
            Lim_complet.append(limitation)

        CRList = [
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
            gamma_complet,
            conf_complet,
            Lim_complet,
        ]

        return CRList


class PTF(BADA4):
    """This class implements the PTF file creator for BADA4 aircraft following BADA4 manual.

    :param AC: parsed aircraft.
    :type AC: bada4.Parse.
    """

    def __init__(self, AC):
        super().__init__(AC)

        self.flightEnvelope = FlightEnvelope(AC)
        self.ARPM = ARPM(AC)

    def create(self, DeltaTemp, saveToPath):
        """This function creates the BADA4 PTF file

        :param saveToPath: path to directory where PTF should be stored [-]
        :param DeltaTemp: deviation from ISA temperature [K]
        :type saveToPath: string.
        :type DeltaTemp: float.
        :returns: NONE
        """

        # 3 different mass levels [kg]
        massList = [
            1.2 * self.AC.OEW,
            self.AC.OEW + 0.7 * (self.AC.MTOW - self.AC.OEW),
            self.AC.MTOW,
        ]
        max_alt_ft = self.AC.hmo

        if max_alt_ft > 51000:
            max_alt_ft = 51000

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
        self, saveToPath, CLList, CRList, DESList, DeltaTemp, massList, altitudeList
    ):
        """This function saves data to PTF file

        :param saveToPath: path to directory where PTF should be stored [-]
        :param CRList: list of PTF data in CRUISE [-].
        :param CLList: list of PTF data in CLIMB [-].
        :param DESList: list of PTF data in DESCENT [-].
        :param DeltaTemp: deviation from ISA temperature [K]
        :param massList: list of aircraft mass [kg]
        :type saveToPath: string.
        :type CRList: list.
        :type CLList: list.
        :type DESList: list.
        :type DeltaTemp: float.
        :type massList: list(float).
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

        filename = saveToPath + self.AC.acName + "_ISA" + ISA + ".PTF"

        [Vcl1, Vcl2, Mcl] = self.flightEnvelope.getSpeedSchedule(phase="Climb")
        [Vcr1, Vcr2, Mcr] = self.flightEnvelope.getSpeedSchedule(phase="Cruise")
        [Vdes1, Vdes2, Mdes] = self.flightEnvelope.getSpeedSchedule(phase="Descent")

        V1cl = min(250, conv.ms2kt(Vcl1))
        V2cl = conv.ms2kt(Vcl2)
        V1des = min(250, conv.ms2kt(Vdes1))
        V2des = conv.ms2kt(Vdes2)
        V1cr = min(250, conv.ms2kt(Vcr1))
        V2cr = conv.ms2kt(Vcr2)

        today = date.today()
        d3 = today.strftime("%b %d %Y")

        acModel = self.AC.model

        file = open(filename, "w")
        file.write(
            "BADA PERFORMANCE FILE                                        %s\n\n" % (d3)
        )
        file = open(filename, "a")
        file.write("AC/Type: %s\n\n" % (acModel))
        file.write(
            " Speeds:   CAS(LO/HI)  Mach   Mass Levels [kg]         Temperature: ISA%s\n"
            % (ISA)
        )
        file.write(
            " climb   - %3d/%3d     %4.3f  low     -  %.0f\n"
            % (V1cl, V2cl, Mcl, massList[0])
        )
        file.write(
            " cruise  - %3d/%3d     %4.3f  nominal -  %-5.0f        Max Alt. [ft]:%7d\n"
            % (V1cr, V2cr, Mcr, massList[1], altitudeList[-1])
        )
        file.write(
            " descent - %3d/%3d     %4.3f  high    -  %0.f\n"
            % (V1des, V2des, Mdes, massList[2])
        )
        file.write(
            "======================================================================================================\n"
        )
        file.write(
            " FL |          CRUISE           |               CLIMB               |             DESCENT             \n"
        )
        file.write(
            "    |  TAS          fuel        |  TAS          ROCD         fuel   |  TAS          ROCD        fuel  \n"
        )
        file.write(
            "    | [kts]       [kg/min]      | [kts]        [fpm]       [kg/min] | [kts]        [fpm]      [kg/min]\n"
        )
        file.write(
            "    |          lo   nom    hi   |         lo    nom    hi    nom    |         lo    nom    hi    nom  \n"
        )
        file.write(
            "======================================================================================================\n"
        )

        # replace NAN values by 0 for printing purposes
        CLList = Nan2Zero(CLList)
        DESList = Nan2Zero(DESList)

        for k in range(0, len(altitudeList)):
            FL = proper_round(altitudeList[k] / 100)
            file.write(
                "%3.0f |  %s   %s %s %s  |  %3.0f   %5.0f %5.0f %5.0f   %5.1f  |  %3.0f   %5.0f %5.0f %5.0f   %5.1f\n"
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
                    DESList[3][k],
                    DESList[4][k],
                )
            )
            file.write(
                "    |                           |                                   | \n"
            )

        file.write(
            "======================================================================================================\n"
        )

    def PTF_cruise(self, massList, altitudeList, DeltaTemp):
        """This function calculates the BADA4 PTF data in CRUISE

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

        massNominal = massList[1]

        [Vcr1, Vcr2, Mcr] = self.flightEnvelope.getSpeedSchedule(phase="Cruise")
        Vcr1 = min(Vcr1, conv.kt2ms(250))

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
            [cas, speedUpdated] = self.ARPM.cruiseSpeed(
                h=H_m,
                mass=massNominal,
                theta=theta,
                delta=delta,
                speedSchedule_default=[Vcr1, Vcr2, Mcr],
                DeltaTemp=DeltaTemp,
            )
            tas_nominal = atm.cas2Tas(cas=cas, delta=delta, sigma=sigma)
            FL = h / 100
            ff = []

            for mass in massList:
                [cas, speedUpdated] = self.ARPM.cruiseSpeed(
                    h=H_m,
                    mass=mass,
                    theta=theta,
                    delta=delta,
                    speedSchedule_default=[Vcr1, Vcr2, Mcr],
                    DeltaTemp=DeltaTemp,
                )
                tas = atm.cas2Tas(cas=cas, delta=delta, sigma=sigma)
                M = atm.tas2Mach(v=tas, theta=theta)

                config = "CR"
                HLid = 0
                LG = "LGUP"
                CL = self.CL(M=M, delta=delta, mass=mass)
                CD = self.CD(M=M, CL=CL, HLid=HLid, LG=LG)
                Drag = self.D(M=M, delta=delta, CD=CD)
                Thrust = Drag
                ThrustMax = self.Thrust(
                    rating="MCRZ", delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
                )

                CL = self.flightEnvelope.CL(delta=delta, mass=mass, M=M, nz=1.2)
                CL_max = self.flightEnvelope.CLmax(M=M, HLid=HLid, LG=LG)

                epsilon = 0.01
                if Thrust > ThrustMax:
                    # "(T)" - as thrust limited
                    ff.append("(T)")

                elif CL > (CL_max + epsilon):
                    # "(B)" - as buffet limited
                    ff.append("(B)")
                else:
                    CT = self.CT(Thrust=Thrust, delta=delta)
                    ff.append(
                        self.ff(
                            CT=CT, delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
                        )
                        * 60
                    )

            TAS_CR_complet.append(f"{conv.ms2kt(tas_nominal):3.0f}")
            if isinstance(ff[0], str):
                FF_CR_LO_complet.append(" " + ff[0] + " ")
            else:
                FF_CR_LO_complet.append(f"{ff[0]:5.1f}")
            if isinstance(ff[1], str):
                FF_CR_NOM_complet.append(" " + ff[1] + " ")
            else:
                FF_CR_NOM_complet.append(f"{ff[1]:5.1f}")
            if isinstance(ff[2], str):
                FF_CR_HI_complet.append(" " + ff[2] + " ")
            else:
                FF_CR_HI_complet.append(f"{ff[2]:5.1f}")

        CRList = [TAS_CR_complet, FF_CR_LO_complet, FF_CR_NOM_complet, FF_CR_HI_complet]

        return CRList

    def PTF_climb(self, massList, altitudeList, DeltaTemp):
        """This function calculates the BADA4 PTF data in CLIMB

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
        conf_LO_complet = []
        conf_NOM_complet = []
        conf_HI_complet = []
        conf_complet = {}

        for mass in massList:
            conf_complet[str(mass)] = []

        massNominal = massList[1]

        [Vcl1, Vcl2, Mcl] = self.flightEnvelope.getSpeedSchedule(phase="Climb")
        Vcl1 = min(Vcl1, conv.kt2ms(250))
        crossAlt = atm.crossOver(cas=Vcl2, Mach=Mcl)

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
            [cas, speedUpdated] = self.ARPM.climbSpeed(
                h=H_m,
                mass=massNominal,
                theta=theta,
                delta=delta,
                DeltaTemp=DeltaTemp,
                speedSchedule_default=[Vcl1, Vcl2, Mcl],
            )
            tas_nominal = atm.cas2Tas(cas=cas, delta=delta, sigma=sigma)
            FL = h / 100

            M_nominal = atm.tas2Mach(v=tas_nominal, theta=theta)
            ff_nominal = (
                self.ff(
                    rating="MCMB",
                    delta=delta,
                    theta=theta,
                    M=M_nominal,
                    DeltaTemp=DeltaTemp,
                )
                * 60
            )

            ROC = []
            for mass in massList:
                [cas, speedUpdated] = self.ARPM.climbSpeed(
                    h=H_m,
                    mass=mass,
                    theta=theta,
                    delta=delta,
                    DeltaTemp=DeltaTemp,
                    speedSchedule_default=[Vcl1, Vcl2, Mcl],
                )
                tas = atm.cas2Tas(cas=cas, delta=delta, sigma=sigma)
                M = atm.tas2Mach(v=tas, theta=theta)

                Thrust = self.Thrust(
                    rating="MCMB", delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
                )
                config = self.flightEnvelope.getConfig(
                    h=H_m,
                    phase="Climb",
                    v=cas,
                    mass=mass,
                    DeltaTemp=DeltaTemp,
                )

                # ensure the continuity of aerodyamic configuration during CLimb phase of flight
                if conf_complet[str(mass)]:
                    prevConf = conf_complet[str(mass)][-1]

                    if config == "TO" and (prevConf == "IC" or prevConf == "CR"):
                        config = prevConf
                    elif config == "IC" and prevConf == "CR":
                        config = prevConf

                [HLid, LG] = self.flightEnvelope.getAeroConfig(config=config)
                CL = self.CL(M=M, delta=delta, mass=mass)
                CD = self.CD(M=M, CL=CL, HLid=HLid, LG=LG)
                Drag = self.D(M=M, delta=delta, CD=CD)

                if H_m < crossAlt:
                    ESF = self.esf(
                        h=H_m, flightEvolution="constCAS", M=M, DeltaTemp=DeltaTemp
                    )
                else:
                    ESF = self.esf(
                        h=H_m, flightEvolution="constM", M=M, DeltaTemp=DeltaTemp
                    )

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
                        )
                    )
                    * 60
                )

                if ROC_val < 0:
                    ROC_val = float("Nan")

                ROC.append(ROC_val)
                conf_complet[str(mass)].append(config)

            TAS_CL_complet.append(conv.ms2kt(tas_nominal))
            ROCD_CL_LO_complet.append(ROC[0])
            ROCD_CL_NOM_complet.append(ROC[1])
            ROCD_CL_HI_complet.append(ROC[2])
            FF_CL_NOM_complet.append(ff_nominal)

        CLList = [
            TAS_CL_complet,
            ROCD_CL_LO_complet,
            ROCD_CL_NOM_complet,
            ROCD_CL_HI_complet,
            FF_CL_NOM_complet,
        ]

        return CLList

    def PTF_descent(self, massList, altitudeList, DeltaTemp):
        """This function calculates the BADA4 PTF data in DESCENT

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
        ROCD_DES_LO_complet = []
        ROCD_DES_NOM_complet = []
        ROCD_DES_HI_complet = []
        FF_DES_NOM_complet = []
        conf_complet = {}

        for mass in massList:
            conf_complet[str(mass)] = []

        massNominal = massList[1]

        [Vdes1, Vdes2, Mdes] = self.flightEnvelope.getSpeedSchedule(phase="Descent")
        Vdes1 = min(Vdes1, conv.kt2ms(250))
        crossAlt = atm.crossOver(cas=Vdes2, Mach=Mdes)

        # for h in altitudeList:
        for h in reversed(altitudeList):
            H_m = conv.ft2m(h)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)
            [cas, speedUpdated] = self.ARPM.descentSpeed(
                h=H_m,
                mass=massNominal,
                theta=theta,
                delta=delta,
                DeltaTemp=DeltaTemp,
                speedSchedule_default=[Vdes1, Vdes2, Mdes],
            )
            tas_nominal = atm.cas2Tas(cas=cas, delta=delta, sigma=sigma)
            FL = h / 100

            M_nominal = atm.tas2Mach(v=tas_nominal, theta=theta)
            ff_nominal = (
                self.ff(
                    rating="LIDL",
                    delta=delta,
                    theta=theta,
                    M=M_nominal,
                    DeltaTemp=DeltaTemp,
                )
                * 60
            )

            ROD = []
            ff_gamma_list = []
            for mass in massList:
                [cas, speedUpdated] = self.ARPM.descentSpeed(
                    h=H_m,
                    mass=mass,
                    theta=theta,
                    delta=delta,
                    DeltaTemp=DeltaTemp,
                    speedSchedule_default=[Vdes1, Vdes2, Mdes],
                )
                tas = atm.cas2Tas(cas=cas, delta=delta, sigma=sigma)
                M = atm.tas2Mach(v=tas, theta=theta)

                Thrust = self.Thrust(
                    rating="LIDL", delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
                )
                config = self.flightEnvelope.getConfig(
                    h=H_m,
                    phase="Descent",
                    v=cas,
                    mass=mass,
                    DeltaTemp=DeltaTemp,
                )

                # ensure the continuity of aerodyamic configuration during Descent phase of flight
                if conf_complet[str(mass)]:
                    prevConf = conf_complet[str(mass)][0]

                    if config == "CR" and (prevConf == "AP" or prevConf == "LD"):
                        config = prevConf
                    elif config == "AP" and prevConf == "LD":
                        config = prevConf

                [HLid, LG] = self.flightEnvelope.getAeroConfig(config=config)
                CL = self.CL(M=M, delta=delta, mass=mass)
                CD = self.CD(M=M, CL=CL, HLid=HLid, LG=LG)
                Drag = self.D(M=M, delta=delta, CD=CD)

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

                # in case of AP & LD thrust is computed to fly a 3? slope
                if config == "AP" or config == "LD":
                    gamma = -3.0
                    temp_const = (theta * const.temp_0) / (
                        theta * const.temp_0 - DeltaTemp
                    )
                    ROCD_gamma = sin(conv.deg2rad(gamma)) * tas * (1 / temp_const)
                    ROCD = conv.m2ft(ROCD_gamma) * 60  # [ft/min]

                    n = 1.0  # aircraft.loadFactor(gamma) - use this in case of L = W * (cos(gamma))
                    CL = self.CL(M=M, delta=delta, mass=mass, nz=n)
                    CD = self.CD(M=M, CL=CL, HLid=HLid, LG=LG)
                    Drag = self.D(M=M, delta=delta, CD=CD)
                    Thrust = (ROCD_gamma * mass * const.g) * temp_const / (
                        ESF * tas
                    ) + Drag
                    CT = self.CT(Thrust=Thrust, delta=delta)
                    ff_gamma = (
                        self.ff(
                            CT=CT, delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
                        )
                        * 60
                    )
                    ff_gamma_list.append(ff_gamma)

                else:
                    ff_gamma_list.append(float("Nan"))

                ROD.append(ROCD)
                conf_complet[str(mass)] = [config] + conf_complet[str(mass)]

            if not isnan(ff_gamma_list[1]):
                ff_nominal = ff_gamma_list[1]

            TAS_DES_complet = [conv.ms2kt(tas_nominal)] + TAS_DES_complet
            ROCD_DES_LO_complet = [-1 * ROD[0]] + ROCD_DES_LO_complet
            ROCD_DES_NOM_complet = [-1 * ROD[1]] + ROCD_DES_NOM_complet
            ROCD_DES_HI_complet = [-1 * ROD[2]] + ROCD_DES_HI_complet
            FF_DES_NOM_complet = [ff_nominal] + FF_DES_NOM_complet

        DESList = [
            TAS_DES_complet,
            ROCD_DES_LO_complet,
            ROCD_DES_NOM_complet,
            ROCD_DES_HI_complet,
            FF_DES_NOM_complet,
        ]

        return DESList


class Bada4Aircraft(BADA4):
    """This class groups the BADA4 performance model classes.

    :param filePath: path to the BADA4 xml formatted file.
    :param acName: ICAO aircraft designation
    :type filePath: str.
    :type acName: str
    """

    def __init__(self, badaVersion, acName, filePath=None, allData=None):
        super().__init__(self)

        self.BADAFamily = BadaFamily(BADA4=True)
        self.BADAFamilyName = "BADA4"
        self.BADAVersion = badaVersion
        self.acName = acName

        if filePath == None:
            self.filePath = configuration.getAircraftPath()
        else:
            self.filePath = filePath

        # check if the aircraft is in the allData dataframe data
        if allData is not None and acName in allData["acName"].values:
            filtered_df = allData[allData["acName"] == acName]

            self.model = Parser.safe_get(filtered_df, "model", None)
            self.engineType = Parser.safe_get(filtered_df, "engineType", None)
            self.engines = Parser.safe_get(filtered_df, "engines", None)
            self.WTC = Parser.safe_get(filtered_df, "WTC", None)
            self.ICAO = Parser.safe_get(filtered_df, "ICAO", None)
            self.MREF = Parser.safe_get(filtered_df, "MREF", None)
            self.WREF = Parser.safe_get(filtered_df, "WREF", None)
            self.LHV = Parser.safe_get(filtered_df, "LHV", None)
            self.n_eng = Parser.safe_get(filtered_df, "n_eng", None)
            self.rho = Parser.safe_get(filtered_df, "rho", None)
            self.TFA = Parser.safe_get(filtered_df, "TFA", None)
            self.p_delta = Parser.safe_get(filtered_df, "p_delta", None)
            self.p_theta = Parser.safe_get(filtered_df, "p_theta", None)
            self.kink = Parser.safe_get(filtered_df, "kink", None)
            self.b = Parser.safe_get(filtered_df, "b", None)
            self.c = Parser.safe_get(filtered_df, "c", None)
            self.max_power = Parser.safe_get(filtered_df, "max_power", None)
            self.p = Parser.safe_get(filtered_df, "p", None)
            self.a = Parser.safe_get(filtered_df, "a", None)
            self.f = Parser.safe_get(filtered_df, "f", None)
            self.ti = Parser.safe_get(filtered_df, "ti", None)
            self.fi = Parser.safe_get(filtered_df, "fi", None)
            self.throttle = Parser.safe_get(filtered_df, "throttle", None)
            self.prop_dia = Parser.safe_get(filtered_df, "prop_dia", None)
            self.max_eff = Parser.safe_get(filtered_df, "max_eff", None)
            self.Hd_turbo = Parser.safe_get(filtered_df, "Hd_turbo", None)
            self.CPSFC = Parser.safe_get(filtered_df, "CPSFC", None)
            self.P = Parser.safe_get(filtered_df, "P", None)
            self.S = Parser.safe_get(filtered_df, "S", None)
            self.HLPosition = Parser.safe_get(filtered_df, "HLPosition", None)
            self.configName = Parser.safe_get(filtered_df, "configName", None)
            self.VFE = Parser.safe_get(filtered_df, "VFE", None)
            self.d = Parser.safe_get(filtered_df, "d", None)
            self.CL_max = Parser.safe_get(filtered_df, "CL_max", None)
            self.bf = Parser.safe_get(filtered_df, "bf", None)
            self.HLids = Parser.safe_get(filtered_df, "HLids", None)
            self.CL_clean = Parser.safe_get(filtered_df, "CL_clean", None)
            self.M_max = Parser.safe_get(filtered_df, "M_max", None)
            self.scalar = Parser.safe_get(filtered_df, "scalar", None)
            self.Mmin = Parser.safe_get(filtered_df, "Mmin", None)
            self.Mmax = Parser.safe_get(filtered_df, "Mmax", None)
            self.CL_Mach0 = Parser.safe_get(filtered_df, "CL_Mach0", None)
            self.MTOW = Parser.safe_get(filtered_df, "MTOW", None)
            self.OEW = Parser.safe_get(filtered_df, "OEW", None)
            self.MFL = Parser.safe_get(filtered_df, "MFL", None)
            self.MTW = Parser.safe_get(filtered_df, "MTW", None)
            self.MZFW = Parser.safe_get(filtered_df, "MZFW", None)
            self.MPL = Parser.safe_get(filtered_df, "MPL", None)
            self.MLW = Parser.safe_get(filtered_df, "MLW", None)
            self.hmo = Parser.safe_get(filtered_df, "hmo", None)
            self.mfa = Parser.safe_get(filtered_df, "mfa", None)
            self.MMO = Parser.safe_get(filtered_df, "MMO", None)
            self.MLE = Parser.safe_get(filtered_df, "MLE", None)
            self.VLE = Parser.safe_get(filtered_df, "VLE", None)
            self.VMO = Parser.safe_get(filtered_df, "VMO", None)
            self.span = Parser.safe_get(filtered_df, "span", None)
            self.length = Parser.safe_get(filtered_df, "length", None)
            self.aeroConfig = Parser.safe_get(filtered_df, "aeroConfig", None)
            self.speedSchedule = Parser.safe_get(filtered_df, "speedSchedule", None)

            # GPF data (temporary)
            self.CVminTO = Parser.safe_get(filtered_df, "CVminTO", None)
            self.CVmin = Parser.safe_get(filtered_df, "CVmin", None)
            self.HmaxPhase = Parser.safe_get(filtered_df, "HmaxPhase", None)
            self.V_des = Parser.safe_get(filtered_df, "V_des", None)
            self.V_cl = Parser.safe_get(filtered_df, "V_cl", None)

            self.flightEnvelope = FlightEnvelope(self)
            self.ARPM = ARPM(self)
            self.OPT = Optimization(self)
            self.PTD = PTD(self)
            self.PTF = PTF(self)

        else:
            self.ACModelAvailable = False
            self.synonymFileAvailable = False
            self.ACinSynonymFile = False

            # check if SYNONYM file exist - since for BADA4 this is not a standard procedure (yet)
            synonymFile = os.path.join(
                self.filePath, "BADA4", badaVersion, "aircraft_model_default.xml"
            )

            if os.path.isfile(synonymFile):
                self.synonymFileAvailable = True

                # if SYNONYM exist - look for synonym based on defined acName
                self.SearchedACName = Parser.parseMappingFile(
                    filePath=self.filePath, badaVersion=badaVersion, acName=acName
                )

                # if cannot find - look for full name (in sub folder names) based on acName (may not be ICAO designator)
                if self.SearchedACName == None:
                    self.SearchedACName = acName
                else:
                    self.ACinSynonymFile = True

            else:
                # if it doesn't exist - look for full name (in sub folder names) based on acName (may not be ICAO designator)
                self.SearchedACName = acName

            acXmlFile = (
                os.path.join(
                    self.filePath,
                    "BADA4",
                    badaVersion,
                    self.SearchedACName,
                    self.SearchedACName,
                )
                + ".xml"
            )

            # look for either found synonym or original full BADA4 model name designator
            if self.SearchedACName is not None:
                if os.path.isfile(acXmlFile):

                    self.ACModelAvailable = True

                    XMLDataFrame = Parser.parseXML(
                        filePath=self.filePath,
                        badaVersion=badaVersion,
                        acName=self.SearchedACName,
                    )
                    GPFDataframe = Parser.parseGPF(
                        filePath=self.filePath, badaVersion=badaVersion
                    )

                    combined_df = Parser.combineXML_GPF(XMLDataFrame, GPFDataframe)

                    self.acName = self.SearchedACName

                    self.model = Parser.safe_get(combined_df, "model", None)
                    self.engineType = Parser.safe_get(combined_df, "engineType", None)
                    self.engines = Parser.safe_get(combined_df, "engines", None)
                    self.WTC = Parser.safe_get(combined_df, "WTC", None)
                    self.ICAO = Parser.safe_get(combined_df, "ICAO", None)
                    self.MREF = Parser.safe_get(combined_df, "MREF", None)
                    self.WREF = Parser.safe_get(combined_df, "WREF", None)
                    self.LHV = Parser.safe_get(combined_df, "LHV", None)
                    self.n_eng = Parser.safe_get(combined_df, "n_eng", None)
                    self.rho = Parser.safe_get(combined_df, "rho", None)
                    self.TFA = Parser.safe_get(combined_df, "TFA", None)
                    self.p_delta = Parser.safe_get(combined_df, "p_delta", None)
                    self.p_theta = Parser.safe_get(combined_df, "p_theta", None)
                    self.kink = Parser.safe_get(combined_df, "kink", None)
                    self.b = Parser.safe_get(combined_df, "b", None)
                    self.c = Parser.safe_get(combined_df, "c", None)
                    self.max_power = Parser.safe_get(combined_df, "max_power", None)
                    self.p = Parser.safe_get(combined_df, "p", None)
                    self.a = Parser.safe_get(combined_df, "a", None)
                    self.f = Parser.safe_get(combined_df, "f", None)
                    self.ti = Parser.safe_get(combined_df, "ti", None)
                    self.fi = Parser.safe_get(combined_df, "fi", None)
                    self.throttle = Parser.safe_get(combined_df, "throttle", None)
                    self.prop_dia = Parser.safe_get(combined_df, "prop_dia", None)
                    self.max_eff = Parser.safe_get(combined_df, "max_eff", None)
                    self.Hd_turbo = Parser.safe_get(combined_df, "Hd_turbo", None)
                    self.CPSFC = Parser.safe_get(combined_df, "CPSFC", None)
                    self.P = Parser.safe_get(combined_df, "P", None)
                    self.S = Parser.safe_get(combined_df, "S", None)
                    self.HLPosition = Parser.safe_get(combined_df, "HLPosition", None)
                    self.configName = Parser.safe_get(combined_df, "configName", None)
                    self.VFE = Parser.safe_get(combined_df, "VFE", None)
                    self.d = Parser.safe_get(combined_df, "d", None)
                    self.CL_max = Parser.safe_get(combined_df, "CL_max", None)
                    self.bf = Parser.safe_get(combined_df, "bf", None)
                    self.HLids = Parser.safe_get(combined_df, "HLids", None)
                    self.CL_clean = Parser.safe_get(combined_df, "CL_clean", None)
                    self.M_max = Parser.safe_get(combined_df, "M_max", None)
                    self.scalar = Parser.safe_get(combined_df, "scalar", None)
                    self.Mmin = Parser.safe_get(combined_df, "Mmin", None)
                    self.Mmax = Parser.safe_get(combined_df, "Mmax", None)
                    self.CL_Mach0 = Parser.safe_get(combined_df, "CL_Mach0", None)
                    self.MTOW = Parser.safe_get(combined_df, "MTOW", None)
                    self.OEW = Parser.safe_get(combined_df, "OEW", None)
                    self.MFL = Parser.safe_get(combined_df, "MFL", None)
                    self.MTW = Parser.safe_get(combined_df, "MTW", None)
                    self.MZFW = Parser.safe_get(combined_df, "MZFW", None)
                    self.MPL = Parser.safe_get(combined_df, "MPL", None)
                    self.MLW = Parser.safe_get(combined_df, "MLW", None)
                    self.hmo = Parser.safe_get(combined_df, "hmo", None)
                    self.mfa = Parser.safe_get(combined_df, "mfa", None)
                    self.MMO = Parser.safe_get(combined_df, "MMO", None)
                    self.MLE = Parser.safe_get(combined_df, "MLE", None)
                    self.VLE = Parser.safe_get(combined_df, "VLE", None)
                    self.VMO = Parser.safe_get(combined_df, "VMO", None)
                    self.span = Parser.safe_get(combined_df, "span", None)
                    self.length = Parser.safe_get(combined_df, "length", None)
                    self.aeroConfig = Parser.safe_get(combined_df, "aeroConfig", None)
                    self.speedSchedule = Parser.safe_get(
                        combined_df, "speedSchedule", None
                    )

                    # GPF data (temporary)
                    self.CVminTO = Parser.safe_get(combined_df, "CVminTO", None)
                    self.CVmin = Parser.safe_get(combined_df, "CVmin", None)
                    self.HmaxPhase = Parser.safe_get(combined_df, "HmaxPhase", None)
                    self.V_des = Parser.safe_get(combined_df, "V_des", None)
                    self.V_cl = Parser.safe_get(combined_df, "V_cl", None)

                    # BADA4.__init__(self, AC_parsed)
                    self.flightEnvelope = FlightEnvelope(self)
                    self.ARPM = ARPM(self)
                    self.OPT = Optimization(self)
                    self.PTD = PTD(self)
                    self.PTF = PTF(self)

                else:
                    # AC name cannot be found
                    raise ValueError(
                        acName + " Cannot be found at path " + self.filePath
                    )

    def __str__(self):
        return f"(BADA4, AC_name: {self.acName}, searched_AC_name: {self.SearchedACName}, model_ICAO: {self.ICAO}, ID: {id(self.AC)})"
