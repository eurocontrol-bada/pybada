"""
pyBADA
Generic BADA4 aircraft performance module
Developed @EUROCONTROL (EIH)
2024
"""

from math import sqrt, pow, pi, isnan, sin, asin

import numpy as np
import os
from datetime import date
import xml.etree.ElementTree as ET

from scipy.optimize import fminbound
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
    """This class implements the BADA4 parsing mechanism to parse xml and
    GPF(xml) BADA4 files."""

    def __init__(self):
        pass

    @staticmethod
    def readMappingFile(filePath):
        """Parses the BADA4 mapping XML file and stores a dictionary of
        aircraft code names and their corresponding XML file paths.

        This function processes the BADA4 aircraft model mapping XML file to
        create a dictionary that maps the aircraft code names to the XML file
        paths for their corresponding BADA models. The mapping file contains
        information about available models for the specified BADA version.

        :param filePath: The path to the directory containing the BADA4
            mapping XML file.
        :type filePath: str
        :raises IOError: If the XML file cannot be found or parsed.
        :return: A dictionary with aircraft code names as keys and
            corresponding file names as values.
        :rtype: dict
        """

        filename = os.path.join(filePath, "aircraft_model_default.xml")

        code_fileName = {}

        if os.path.isfile(filename):
            try:
                tree = ET.parse(filename)
                root = tree.getroot()
            except Exception:
                raise IOError(filename + " not found or in correct format")

            for child in root.iter("MAP"):
                code = child.find("code").text
                file = child.find("file").text

                code_fileName[code] = file

        return code_fileName

    @staticmethod
    def parseMappingFile(filePath, acName):
        """Retrieves the file name for a given aircraft name from the parsed
        BADA4 mapping file.

        This function uses the readMappingFile method to parse the BADA4 XML
        mapping file and returns the file name associated with the given
        aircraft name (acName).

        :param filePath: The path to the directory containing the BADA4
            mapping XML file.
        :param acName: The aircraft code name for which the corresponding file
            is being requested.
        :type filePath: str
        :type acName: str
        :return: The file name corresponding to the aircraft code, or None if
            the aircraft code is not found.
        :rtype: str or None
        """

        code_fileName = Parser.readMappingFile(filePath)
        if acName in code_fileName:
            fileName = code_fileName[acName]
            return fileName
        else:
            return None

    @staticmethod
    def parseXML(filePath, acName):
        """Parses the BADA4 XML file for a specific aircraft model and
        extracts various parameters.

        This function parses the BADA4 aircraft XML file for a given aircraft
        model (acName). It retrieves general information about the aircraft,
        engine type, aerodynamic configurations, performance parameters, and
        more.

        :param filePath: The path to the folder containing the BADA4 XML file.
        :param acName: The aircraft code name for which the XML file is being
            parsed.
        :type filePath: str
        :type acName: str
        :raises IOError: If the XML file cannot be found or parsed.
        :return: A pandas DataFrame containing the parsed data for the
            specified aircraft.
        :rtype: pd.DataFrame
        """

        acXmlFile = os.path.join(filePath, acName, acName) + ".xml"

        try:
            tree = ET.parse(acXmlFile)
            root = tree.getroot()
        except Exception:
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

            if (
                LGDN is not None
            ):  # Landing gear NOT allowed in clean configuration
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
    def parseGPF(filePath):
        """Parses the BADA4 GPF XML file and extracts key performance factors.

        This function processes the BADA4 GPF XML file to extract data related
        to various flight performance factors, such as minimum climb/descent
        speeds, maximum altitude limits for different phases of flight, and
        speed schedules for climb and descent.

        :param filePath: The path to the directory containing the BADA4 GPF
            XML file.
        :type filePath: str
        :raises IOError: If the GPF XML file cannot be found or parsed.
        :return: A pandas DataFrame containing the parsed GPF data.
        :rtype: pd.DataFrame
        """

        filename = os.path.join(filePath, "GPF.xml")

        tree = ET.parse(filename)

        try:
            tree = ET.parse(filename)
            root = tree.getroot()
        except Exception:
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
        """Combines the parsed aircraft XML DataFrame with the parsed GPF
        DataFrame.

        This function merges two DataFrames, one containing the parsed
        aircraft-specific data (from the XML file) and the other containing
        the parsed GPF data. This combination provides a unified set of
        aircraft performance data along with general performance factors.

        :param XMLDataFrame: A DataFrame containing the parsed aircraft XML
            data.
        :param GPFDataframe: A DataFrame containing the parsed GPF data.
        :type XMLDataFrame: pd.DataFrame
        :type GPFDataframe: pd.DataFrame
        :return: A combined DataFrame with both aircraft and GPF data.
        :rtype: pd.DataFrame
        """

        # Combine data with GPF data (temporary solution)
        combined_df = pd.concat(
            [
                XMLDataFrame.reset_index(drop=True),
                GPFDataframe.reset_index(drop=True),
            ],
            axis=1,
        )

        return combined_df

    @staticmethod
    def parseAll(badaVersion, filePath=None):
        """Parses all BADA4 XML files and combines the data into a single
        DataFrame.

        This function parses both the BADA4 aircraft XML files and the GPF
        (General Performance Factors) file. It combines the data from both
        sources and creates a unified DataFrame that contains all aircraft and
        performance-related data for the specified BADA version.

        :param badaVersion: The version of BADA (e.g., "4.3") to be parsed.
        :param filePath: The path to the folder containing the BADA4 XML files
            and GPF file. If None, the default path from the configuration
            will be used.
        :type badaVersion: str
        :type filePath: str, optional
        :raises IOError: If any of the necessary XML files cannot be found or
            parsed.
        :return: A DataFrame containing the combined data from all parsed
            aircraft and GPF files.
        :rtype: pd.DataFrame
        """

        if filePath is None:
            filePath = configuration.getBadaVersionPath(
                badaFamily="BADA4", badaVersion=badaVersion
            )
        else:
            filePath = filePath

        # parsing GPF file
        GPFparsedDataframe = Parser.parseGPF(filePath)

        # retrieving mapping data
        code_fileName = Parser.readMappingFile(filePath)

        # get names of all the folders in the main BADA model folder to search for XML files
        subfolders = configuration.list_subfolders(filePath)

        # Initialize an empty list to collect DataFrames
        mapping_dfs = []

        if code_fileName:
            for code in code_fileName:
                file = code_fileName[code]

                if file in subfolders:
                    # Parse the original XML of a model
                    df = Parser.parseXML(filePath, file)

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
            df = Parser.parseXML(filePath, file)

            # Combine data with GPF data (temporary solution)
            combined_df = pd.concat(
                [
                    df.reset_index(drop=True),
                    GPFparsedDataframe.reset_index(drop=True),
                ],
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


class BADA4(Airplane, Bada):
    """This class implements the part of BADA4 performance model that will be
    used in other classes following the BADA4 manual.

    :param AC: Aircraft object {BADA4}.
    :type AC: bada4Aircraft.
    """

    def __init__(self, AC):
        """Initializes the BADA4 class by inheriting from the parent Airplane
        class and assigns the parsed aircraft data to the class instance.

        :param AC: Aircraft object {BADA4}.
        :type AC: bada4Aircraft.
        """

        super().__init__()
        self.AC = AC

    def CL(self, delta, mass, M, nz=1.0):
        """Computes the aircraft's lift coefficient based on the current Mach
        number (M), normalized air pressure (delta), aircraft mass, and load
        factor (nz).

        :param M: Mach number [-].
        :param delta: Normalized air pressure [-].
        :param mass: Aircraft mass [kg].
        :param nz: Load factor [-]. Default is 1.0 (straight and level
            flight).
        :type M: float
        :type delta: float
        :type mass: float
        :type nz: float
        :return: Lift coefficient (CL) [-].
        :rtype: float
        """

        return (
            2
            * mass
            * const.g
            * nz
            / (delta * const.p_0 * const.Agamma * pow(M, 2) * self.AC.S)
        )

    def CLPoly(self, M):
        """Computes the lift coefficient polynomial for the given Mach number
        (M).

        This method uses a 5th-degree polynomial defined by the coefficients
        in the parsed aircraft data (self.AC.bf).

        :param M: Mach number [-].
        :type M: float
        :return: Lift coefficient (polynomial approximation) [-].
        :rtype: float
        """

        CLpoly = 0.0
        for i in range(5):
            CLpoly += self.AC.bf[i] * pow(M, i)

        return CLpoly

    def CLmax(self, M, HLid, LG):
        """Computes the maximum lift coefficient (CLmax) for the given Mach
        number (M), high-lift device (HLid) position, and landing gear (LG)
        configuration.

        If the aircraft is in a clean configuration (HLid == 0 and LG ==
        "LGUP"), the method interpolates or extrapolates the lift coefficient
        based on Mach number.

        :param M: Mach number [-].
        :param HLid: High-lift device position [-].
        :param LG: Landing gear position [LGUP/LGDN] [-].
        :type M: float
        :type HLid: float
        :type LG: str
        :return: Maximum lift coefficient (CLmax) [-].
        :rtype: float
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
                        + 4
                        * self.AC.bf[4]
                        * self.AC.Mmax
                        * self.AC.Mmax
                        * self.AC.Mmax
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
        """Computes the fuel flow coefficient at idle throttle for JET and
        TURBOPROP engines.

        :param delta: Normalized pressure [-].
        :param theta: Normalized temperature [-].
        :param M: Mach speed [-].
        :type delta: float
        :type theta: float
        :type M: float
        :return: Idle fuel flow coefficient [-].
        :rtype: float
        """

        if self.AC.engineType == "JET":
            CF_idle = 0.0
            for i in range(0, 3):
                for j in range(0, 3):
                    CF_idle += self.AC.fi[i * 3 + j] * (delta**j) * (M**i)

            if self.AC.BADAVersion == "4.2":
                CF_idle = CF_idle * pow(delta, -1) * pow(theta, -0.5)
            elif (
                self.AC.BADAVersion == "4.3" or self.AC.BADAVersion == "DUMMY"
            ):
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

            if self.AC.BADAVersion == "4.2":
                CF_idle = CF_idle * pow(delta, -1) * pow(theta, -0.5)
            elif (
                self.AC.BADAVersion == "4.3" or self.AC.BADAVersion == "DUMMY"
            ):
                CF_idle = CF_idle * pow(delta, -1)

        return CF_idle

    def CF(self, delta, theta, DeltaTemp, **kwargs):
        """Computes the fuel flow coefficient (CF) for JET, TURBOPROP, and
        PISTON engines.

        :param delta: Normalized pressure [-].
        :param theta: Normalized temperature [-].
        :param DeltaTemp: Temperature deviation from ISA [K].
        :param kwargs: Optional parameters including 'rating', 'deltaT', 'CT',
            or 'M'.
        :type delta: float
        :type theta: float
        :type DeltaTemp: float
        :return: Fuel flow coefficient [-].
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
                        CF_gen_rating += (
                            self.AC.f[i * 5 + j] * (CT**j) * (M**i)
                        )

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
                    deltaT=deltaT,
                    theta=theta,
                    delta=delta,
                    M=M,
                    DeltaTemp=DeltaTemp,
                )
                for i in range(0, 5):
                    for j in range(0, 5):
                        CF_gen_deltaT += (
                            self.AC.f[i * 5 + j] * (CT_deltaT**j) * (M**i)
                        )

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
                        CF_gen_rating += (
                            self.AC.f[i * 5 + j] * (CP**j) * (M**i)
                        )

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
                    CP_rating = self.CP(
                        rating=rating, theta=theta, delta=delta, M=M
                    )
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
                    deltaT=deltaT,
                    theta=theta,
                    delta=delta,
                    M=M,
                    DeltaTemp=DeltaTemp,
                )
                for i in range(0, 5):
                    for j in range(0, 5):
                        CF_gen_deltaT += (
                            self.AC.f[i * 5 + j] * (CP_deltaT**j) * (M**i)
                        )

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
        """Computes the thrust coefficient (CT) based on engine type, throttle
        setting, or engine rating. The thrust coefficient is calculated
        differently for JET, TURBOPROP, and PISTON engines based on normalized
        pressure, temperature, Mach number, and other inputs.

        :param delta: Normalized pressure [-].
        :type delta: float
        :param kwargs: Optional parameters:
            - 'rating': Engine rating {MCMB, MCRZ, MTKF, LIDL}.
            - 'deltaT': Direct throttle parameter [-].
            - 'Thrust': Direct thrust value [N].
            - 'M': Mach number [-].
            - 'theta': Normalized temperature [-].
            - 'DeltaTemp': Deviation from the standard ISA temperature [K].
        :return: Thrust coefficient (CT) [-].
        :rtype: float
        :raises: ValueError: If an invalid rating is provided.
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
                    rating="MCMB",
                    theta=theta,
                    delta=delta,
                    M=M,
                    DeltaTemp=DeltaTemp,
                ):
                    raise ValueError(
                        "Throttle parameter value result in CT > CT_MCMB"
                        + deltaT
                    )
                elif CT < self.CT_LIDL(delta=delta, M=M):
                    raise ValueError(
                        "Throttle parameter value result in CT < CT_LIDL"
                        + deltaT
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
                    CT = self.CT_nonLIDL(
                        rating=rating, theta=theta, delta=delta, M=M
                    )

            # deltaT - direct throttle as input parameter
            elif "deltaT" in kwargs:
                deltaT = checkArgument("deltaT", **kwargs)

                CP = self.CP(deltaT=deltaT, M=M)
                CT = CP / M

                if CT > self.CT_nonLIDL(
                    rating="MCMB", theta=theta, delta=delta, M=M
                ):
                    raise ValueError(
                        "Throttle parameter value result in CT > CT_MCMB"
                        + deltaT
                    )
                elif CT < self.CT_LIDL(theta=theta, delta=delta, M=M):
                    raise ValueError(
                        "Throttle parameter value result in CT < CT_LIDL"
                        + deltaT
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
        """Computes the thrust coefficient (CT) for the LIDL (idle) rating.

        :param delta: Normalized pressure [-].
        :param theta: Normalized temperature [-].
        :param M: Mach number [-].
        :param CP: Power coefficient (for PISTON engines) [-].
        :type delta: float.
        :type theta: float.
        :type M: float.
        :type CP: float.
        :return: Idle thrust coefficient (CT) [-].
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
                self.AC.ti[21] / M
                + self.AC.ti[22] * delta / M
                + self.AC.ti[23] * M**3
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
        """Computes the thrust coefficient (CT) for non-LIDL ratings {MCMB,
        MCRZ}.

        :param delta: Normalized pressure [-].
        :param theta: Normalized temperature [-].
        :param M: Mach number [-].
        :param CP: Power coefficient (for PISTON engines) [-].
        :type delta: float.
        :type theta: float.
        :type M: float.
        :type CP: float.
        :return: Thrust coefficient (CT) [-].
        :rtype: float
        """

        if self.AC.engineType == "JET":
            rating = checkArgument("rating", **kwargs)
            DeltaTemp = checkArgument("DeltaTemp", **kwargs)

            # deltaT below kink point -> flat-rated area
            deltaTFlat = 0.0
            for i in range(0, 6):
                for j in range(0, 6):
                    deltaTFlat += (
                        self.AC.b[rating][i * 6 + j] * (M**j) * (delta**i)
                    )

            # deltaT above kink point -> temperature-rated area
            thetaT = theta * (1 + (M**2) * (const.Agamma - 1.0) / 2.0)
            deltaTTemp = 0.0
            for i in range(0, 5):
                for j in range(0, 5):
                    deltaTTemp += (
                        self.AC.c[rating][i * 5 + j] * (M**j) * (thetaT**i)
                    )

            for i in range(5, 9):
                for j in range(0, 5):
                    deltaTTemp += (
                        self.AC.c[rating][i * 5 + j]
                        * (M**j)
                        * (delta ** (i - 4))
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
            propEff = self.propEfficiency(
                Wp=Wp * self.AC.n_eng, sigma=sigma, tas=TAS
            )
            CT = propEff * const.a_0 * CP / (delta * TAS)

        return CT

    def CPmax(self, rating, delta, theta, M):
        """Computes the maximum engine power coefficient (CPmax) for TURBOPROP
        engines.

        :param rating: Throttle setting {MCMB, MCRZ}.
        :param delta: Normalized pressure [-].
        :param theta: Normalized temperature [-].
        :param M: Mach number [-].
        :type rating: str.
        :type delta: float.
        :type theta: float.
        :type M: float.
        :return: Maximum engine power coefficient (CPmax) [-].
        :rtype: float.
        :raises: ValueError: If the rating is unknown.
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
        """Computes the power coefficient (CP) for TURBOPROP and PISTON
        engines.

        :param rating: Throttle setting {MCMB, MCRZ, LIDL}.
        :param deltaT: Direct throttle parameter [-].
        :param delta: Normalized pressure [-].
        :param theta: Normalized temperature [-].
        :param sigma: Normalized density [-] (for piston engines).
        :param M: Mach number [-].
        :type rating: str.
        :type deltaT: float.
        :type delta: float.
        :type theta: float.
        :type sigma: float.
        :type M: float.
        :return: Power coefficient (CP) [-].
        :rtype: float.
        :raises: ValueError if an unknown rating is provided.
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
                        deltaT += (
                            self.AC.p[rating][i * 6 + j] * (M**j) * (theta**i)
                        )

                CP = 0.0
                for i in range(0, 6):
                    for j in range(0, 6):
                        CP += self.AC.a[i * 6 + j] * (M**j) * (deltaT**i)

                CPmax = self.CPmax(
                    rating=rating, theta=theta, delta=delta, M=M
                )
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
                conv.hp2W(self.AC.P)
                * self.AC.n_eng
                / (self.AC.WREF * const.a_0)
            )

            # deltaT - direct throttle as input parameter
            if "deltaT" in kwargs:
                deltaT = checkArgument("deltaT", **kwargs)

            # rating as input parameter
            elif "rating" in kwargs:
                rating = checkArgument("rating", **kwargs)

                if rating == "LIDL":
                    if self.AC.BADAVersion == "4.2":
                        deltaT = 0.0
                    elif (
                        self.AC.BADAVersion == "4.3"
                        or self.AC.BADAVersion == "DUMMY"
                    ):
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
        """Computes the drag coefficient (CD) in a clean configuration based
        on the Mach number (M) and the lift coefficient (CL).

        :param M: Mach number [-].
        :param CL: Lift coefficient [-].
        :type M: float.
        :type CL: float.
        :return: Drag coefficient (CD) in clean configuration [-].
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
        self,
        HLid,
        LG,
        CL,
        M,
        expedite=False,
        speedBrakes={"deployed": False, "value": 0.03},
        **kwargs,
    ):
        """Computes the drag coefficient (CD) based on the Mach number (M),
        lift coefficient (CL), high lift devices (HLid), landing gear position
        (LG), and speed brakes status. The drag coefficient is calculated for
        both clean and non-clean configurations.

        :param M: Mach number [-].
        :param CL: Lift coefficient [-].
        :param HLid: High lift devices position [-].
        :param LG: Landing gear position, [LGUP/LGDN] [-].
        :param speedBrakes: Dictionary indicating if speed brakes are deployed
            and their value. Default: {"deployed": False, "value": 0.03}.
        :param expedite: Flag indicating if expedite descent is used (default is False).
        :type M: float.
        :type CL: float.
        :type HLid: float.
        :type LG: string.
        :type speedBrakes: dict.
        :type expedite: bool
        :returns: Drag coefficient [-].
        :rtype: float.
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

        # implementation of expedite descent
        if expedite:
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
        """Computes the aerodynamic lift based on the Mach number (M),
        normalized air pressure (delta), and lift coefficient (CL).

        :param M: Mach number [-].
        :param delta: Normalized air pressure [-].
        :param CL: Lift coefficient [-].
        :type M: float.
        :type delta: float.
        :type CL: float.
        :returns: Aerodynamic lift [N].
        :rtype: float.
        """

        return 0.5 * delta * const.p_0 * const.Agamma * M * M * self.AC.S * CL

    def D(self, delta, M, CD):
        """Computes the thrust based on throttle settings, normalized air
        pressure (delta), and other flight parameters.

        :param rating: Throttle setting {MCMB, MCRZ, LIDL}.
        :param deltaT: Direct throttle parameter [-].
        :param delta: Normalized air pressure [-].
        :param theta: Normalized temperature [-].
        :param M: Mach number [-].
        :param DeltaTemp: Temperature deviation with respect to ISA [K].
        :type rating: string.
        :type deltaT: float.
        :type delta: float.
        :type theta: float.
        :type M: float.
        :type DeltaTemp: float.
        :returns: Thrust [N].
        :rtype: float.
        """

        return 0.5 * delta * const.p_0 * const.Agamma * M * M * self.AC.S * CD

    def Thrust(self, delta, **kwargs):
        """Computes the maximum thrust produced by the aircraft based on
        throttle settings, normalized air pressure (delta), and other flight
        parameters like Mach number and temperature deviation.

        :param rating: Throttle setting {MCMB (Max Climb), MCRZ (Max Cruise),
            LIDL (Idle)}.
        :param deltaT: Direct throttle parameter for intermediate throttle
            levels [-].
        :param delta: Normalized air pressure relative to sea-level pressure
            (ISA) [-].
        :param theta: Normalized temperature relative to sea-level temperature
            (ISA) [-].
        :param M: Mach number, the ratio of the aircraft's speed to the speed
            of sound [-].
        :param DeltaTemp: Temperature deviation from the International
            Standard Atmosphere (ISA) [K].
        :type rating: string (optional).
        :type deltaT: float (optional).
        :type delta: float.
        :type theta: float (optional).
        :type M: float (optional).
        :type DeltaTemp: float (optional).
        :returns: Thrust force produced by the engine [N].
        :rtype: float.
        """

        CT = self.CT(delta=delta, **kwargs)

        return delta * self.AC.WREF * CT

    def ff(self, delta, theta, DeltaTemp, **kwargs):
        """Computes the fuel flow (FF) of the aircraft based on engine
        throttle settings, normalized air pressure (delta), normalized
        temperature (theta), and other relevant parameters.

        :param rating: Throttle setting {MCMB (Max Climb), MCRZ (Max Cruise),
            LIDL (Idle), TAXI}. If 'TAXI' is selected, the taxi fuel allowance
            is used for ground operations.
        :param deltaT: Direct throttle parameter for intermediate throttle
            levels [-].
        :param delta: Normalized air pressure relative to sea-level pressure
            (ISA) [-].
        :param theta: Normalized temperature relative to sea-level temperature
            (ISA) [-].
        :param M: Mach number, the ratio of the aircraft's speed to the speed
            of sound [-].
        :param DeltaTemp: Temperature deviation from the International
            Standard Atmosphere (ISA) [K].
        :type rating: string (optional).
        :type deltaT: float (optional).
        :type delta: float.
        :type theta: float.
        :type M: float.
        :type DeltaTemp: float.
        :returns: Fuel flow rate in kilograms per second [kg s^-1].
        :rtype: float.
        """

        if "rating" in kwargs:
            rating = checkArgument("rating", **kwargs)
            if rating == "TAXI":
                if self.AC.TFA is not None:
                    return self.AC.TFA / 60
                else:
                    return None

        CF = self.CF(delta=delta, theta=theta, DeltaTemp=DeltaTemp, **kwargs)

        if self.AC.BADAVersion == "4.2":
            return (
                delta
                * pow(theta, 0.5)
                * self.AC.WREF
                * const.a_0
                * CF
                / self.AC.LHV
            )

        elif self.AC.BADAVersion == "4.3" or self.AC.BADAVersion == "DUMMY":
            return (
                pow(delta, self.AC.p_delta)
                * pow(theta, self.AC.p_theta)
                * self.AC.WREF
                * const.a_0
                * CF
                / self.AC.LHV
            )

    def ROCD(self, T, D, v, mass, ESF, h, DeltaTemp):
        """Computes the Rate of Climb or Descent (ROCD) of the aircraft based
        on the provided thrust, drag, true airspeed, mass, and Energy Share
        Factor (ESF).

        :param T: Thrust produced by the engines [N].
        :param D: Drag acting on the aircraft [N].
        :param v: True airspeed (TAS) of the aircraft [m/s].
        :param mass: Current aircraft mass [kg].
        :param ESF: Energy Share Factor, which controls the allocation of
            excess thrust between climb/descent and acceleration [-].
        :param h: Altitude of the aircraft above sea level [m].
        :param DeltaTemp: Temperature deviation from the International
            Standard Atmosphere (ISA) [K].
        :type T: float.
        :type D: float.
        :type v: float.
        :type mass: float.
        :type ESF: float.
        :type h: float.
        :type DeltaTemp: float.
        :returns: Rate of climb or descent (ROCD) in meters per second [m/s].
        :rtype: float.
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        temp = theta * const.temp_0

        ROCD = (
            ((temp - DeltaTemp) / temp) * (T - D) * v * ESF / (mass * const.g)
        )

        return ROCD

    def controlLawThrust(self, ROCD, D, v, mass, ESF, h, DeltaTemp):
        """Computes the required thrust based on the TEM (Thrust-Equilibrium
        Method) control law. This calculation takes into account the
        aircraft's rate of climb/descent (ROCD), drag (D), true airspeed (v),
        and energy share factor (ESF), along with altitude and temperature
        deviations.

        :param h: Altitude of the aircraft above sea level [m].
        :param ROCD: Rate of climb or descent of the aircraft [m/s]. Positive
            value for climb, negative for descent.
        :param D: Drag force acting on the aircraft [N].
        :param v: True airspeed (TAS) of the aircraft [m/s].
        :param mass: Current aircraft mass [kg].
        :param ESF: Energy Share Factor, determining the fraction of excess
            thrust used for acceleration or climb/descent [-].
        :param DeltaTemp: Temperature deviation from the International
            Standard Atmosphere (ISA) [K].
        :type h: float.
        :type ROCD: float.
        :type D: float.
        :type v: float.
        :type mass: float.
        :type ESF: float.
        :type DeltaTemp: float.
        :returns: Thrust required to maintain the specified rate of climb or
            descent [N].
        :rtype: float.
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        temp = theta * const.temp_0

        if ROCD == 0.0 or ESF == 0.0:
            thrust = (temp / (temp - DeltaTemp)) * (
                ROCD * mass * const.g
            ) / v + D
        else:
            thrust = (temp / (temp - DeltaTemp)) * (ROCD * mass * const.g) / (
                ESF * v
            ) + D

        return thrust

    def propEfficiency(self, Wp, sigma, tas):
        """Computes the propeller efficiency of a piston or turboprop engine
        using momentum theory. The calculation estimates the efficiency based
        on power, air density, and true airspeed.

        :param Wp: Total engine power output for all engines [W].
        :param sigma: Normalized air density relative to sea-level density
            [-].
        :param tas: True airspeed (TAS) of the aircraft [m/s].
        :type Wp: float.
        :type sigma: float.
        :type tas: float.
        :returns: Propeller efficiency as a dimensionless ratio [-].
        :rtype: float.
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
    """This class is a BADA4 aircraft subclass and implements the flight
    envelope caclulations following the BADA4 manual.

    :param AC: Aircraft object {BADA4}.
    :type AC: bada4Aircraft.
    """

    def __init__(self, AC):
        super().__init__(AC)

    def maxM(self, LG):
        """Computes the maximum allowable Mach speed (Mmax) based on the
        kinematic limitations of the aircraft and the position of the landing
        gear (LG).

        :param LG: Landing gear position, either "LGUP" (gear up) or "LGDN"
            (gear down).
        :type LG: str
        :return: The maximum allowable Mach number (Mmax) based on the
            aircraft's limitations.
        :rtype: float
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
        """Computes the maximum allowable Calibrated Airspeed (CASmax) based
        on the position of the high-lift devices (HLid) and landing gear (LG),
        following the kinematic limitations of the aircraft.

        :param HLid: Position of high-lift devices (0 for clean configuration,
            >0 for extended).
        :type HLid: float
        :param LG: Landing gear position, either "LGUP" (gear up) or "LGDN"
            (gear down).
        :type LG: str
        :return: The maximum allowable Calibrated Airspeed (CASmax) in meters
            per second [m/s].
        :rtype: float
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
        """Computes the maximum speed (CAS)

        :param h: Altitude in meters [m].
        :param HLid: High-lift devices position [-].
        :param LG: Landing gear position, either "LGUP" (gear up) or "LGDN"
            (gear down) [-].
        :param theta: Normalized temperature [-].
        :param delta: Normalized pressure [-].
        :param mass: Aircraft mass in kilograms [kg].
        :param nz: Load factor [-], default is 1.0.
        :type h: float.
        :type HLid: float.
        :type LG: str.
        :type theta: float.
        :type delta: float.
        :type mass: float.
        :type nz: float.
        :return: Maximum allowable calibrated airspeed (CAS) in meters per
            second [m/s].
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
                    Mach=min(M, M_buffet),
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
                )

                # if M_buffet == float('-inf'):
                # VMax = float('-inf')
            else:
                VMax = self.maxCAS(HLid=HLid, LG=LG)
        else:
            VMax = self.maxCAS(HLid=HLid, LG=LG)

        return VMax

    def Vmax_thrustLimited(self, h, mass, DeltaTemp, rating, config):
        """Computes the maximum speed (CAS) within the certified flight
        envelope while accounting for thrust limitations at a given altitude,
        temperature deviation, and configuration.

        :param h: Altitude in meters [m].
        :param mass: Aircraft mass in kilograms [kg].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param rating: Engine rating, options include "MTKF", "MCMB", "MCRZ"
            [-].
        :param config: Aircraft aerodynamic configuration, such as "TO", "IC",
            "CR" [-].
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type rating: str.
        :type config: str.
        :return: Maximum thrust-limited calibrated airspeed (CAS) in meters
            per second [m/s].
        :rtype: float.
        """

        [HLid, LG] = self.getAeroConfig(config=config)

        [theta, delta, sigma] = atm.atmosphereProperties(
            h=h, DeltaTemp=DeltaTemp
        )

        VmaxCertified = self.VMax(
            h=h, HLid=HLid, LG=LG, delta=delta, theta=theta, mass=mass, nz=1.0
        )
        VminCertified = self.VStall(
            theta=theta, delta=delta, mass=mass, HLid=HLid, LG=LG, nz=1.0
        )

        if VminCertified is None or VmaxCertified is None:
            return None

        maxCASList = []
        for CAS in np.linspace(
            VminCertified, VmaxCertified, num=200, endpoint=True
        ):
            [M, CAS, TAS] = atm.convertSpeed(
                v=conv.ms2kt(CAS),
                speedType="CAS",
                theta=theta,
                delta=delta,
                sigma=sigma,
            )

            maxThrust = self.Thrust(
                delta=delta,
                theta=theta,
                M=M,
                rating=rating,
                DeltaTemp=DeltaTemp,
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
        """Computes the maximum allowable Mach number (M) under buffet
        limitations, where the lift coefficient (CL) cannot exceed the maximum
        lift coefficient (CL_max) for the given Mach number and configuration.

        :param HLid: High-lift devices position [-].
        :param LG: Landing gear position, either "LGUP" (gear up) or "LGDN"
            (gear down) [-].
        :param delta: Normalized pressure [-].
        :param mass: Aircraft mass in kilograms [kg].
        :param nz: Load factor [-], default is 1.0.
        :type HLid: float.
        :type LG: str.
        :type delta: float.
        :type mass: float.
        :type nz: float.
        :return: Maximum allowable Mach number (M) limited by buffet
            conditions.
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
        """Computes the minimum Mach number (M) applying buffet limitations,
        where the lift coefficient (CL) must not exceed the maximum lift
        coefficient (CL_max) for the given Mach number and aerodynamic
        configuration.

        :param HLid: High-lift devices position [-].
        :param LG: Landing gear position, either "LGUP" (gear up) or "LGDN"
            (gear down) [-].
        :param delta: Normalized pressure [-].
        :param theta: Normalized temperature [-].
        :param mass: Aircraft mass in kilograms [kg].
        :param nz: Load factor [-], default is 1.0.
        :type HLid: float.
        :type LG: str.
        :type delta: float.
        :type theta: float.
        :type mass: float.
        :type nz: float.
        :returns: Minimum Mach number (M) limited by buffet conditions.
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
                    cas=conv.kt2ms(self.AC.VMO),
                    theta=theta,
                    delta=delta,
                    sigma=sigma,
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
        """Computes the minimum speed (CAS) for the given aerodynamic
        configuration, accounting for stall speed and other configuration-
        based limitations.

        :param config: Aircraft configuration, options include "CR", "IC",
            "TO", "AP", "LD" [-].
        :param delta: Normalized pressure [-].
        :param theta: Normalized temperature [-].
        :param mass: Aircraft mass in kilograms [kg].
        :type config: str.
        :type delta: float.
        :type theta: float.
        :type mass: float.
        :returns: Minimum calibrated airspeed (CAS) in meters per second
            [m/s].
        :rtype: float.
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
                    theta=theta,
                    delta=delta,
                    mass=mass,
                    HLid=HLid,
                    LG=LG,
                    nz=1.0,
                )
            else:
                Vmin = self.AC.CVmin * self.VStall(
                    theta=theta,
                    delta=delta,
                    mass=mass,
                    HLid=HLid,
                    LG=LG,
                    nz=1.0,
                )

        return Vmin

    def VStall(self, mass, HLid, LG, nz=1.0, **kwargs):
        """Calculates the stall speed (CAS) for the given aerodynamic
        configuration and load factor, taking into account altitude and
        temperature deviations from ISA.

        :param HLid: High-lift devices position [-].
        :param LG: Landing gear position, either "LGUP" (gear up) or "LGDN"
            (gear down) [-].
        :param delta: Normalized pressure [-].
        :param theta: Normalized temperature [-].
        :param mass: Aircraft mass in kilograms [kg].
        :param nz: Load factor [-], default is 1.0.
        :param h: Altitude above mean sea level (AMSL) in meters [m].
        :param DeltaTemp: Temperature deviation from ISA in Kelvin [K].
        :type HLid: float.
        :type LG: str.
        :type delta: float.
        :type theta: float.
        :type mass: float.
        :type nz: float.
        :type h: float, optional.
        :type DeltaTemp: float, optional.
        :returns: Stall speed in calibrated airspeed (CAS) [m/s].
        :rtype: float.
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

    def Vx(self, h, mass, DeltaTemp, rating, HLid, LG):
        """Computes the best angle of climb (Vx) speed.

        :param h: Altitude in meters [m].
        :param mass: Aircraft mass in kilograms [kg].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param rating: Aircraft engine rating (e.g., 'MTKF', 'MCMB', 'MCRZ').
        :param HLid: High-lift device configuration identifier.
        :param LG: Landing gear configuration identifier.
        :type h: float
        :type mass: float
        :type DeltaTemp: float
        :type rating: str
        :type HLid: str
        :type LG: str
        :returns: Best angle of climb speed (Vx) in meters per second [m/s].
        :rtype: float
        """

        [theta, delta, sigma] = atm.atmosphereProperties(
            h=h, DeltaTemp=DeltaTemp
        )

        VmaxCertified = self.VMax(
            h=h, delta=delta, theta=theta, HLid=HLid, LG=LG, mass=mass, nz=1.0
        )
        VminCertified = self.VStall(
            delta=delta, theta=theta, mass=mass, HLid=HLid, LG=LG, nz=1.0
        )

        excessThrustList = []
        VxList = []

        for CAS in np.linspace(
            VminCertified, VmaxCertified, num=200, endpoint=True
        ):
            M = atm.cas2Mach(cas=CAS, theta=theta, delta=delta, sigma=sigma)

            maxThrust = self.Thrust(
                rating=rating,
                delta=delta,
                theta=theta,
                M=M,
                DeltaTemp=DeltaTemp,
            )
            CL = self.CL(M=M, delta=delta, mass=mass, nz=1.0)
            CD = self.CD(M=M, CL=CL, HLid=HLid, LG=LG)
            Drag = self.D(M=M, delta=delta, CD=CD)

            excessThrustList.append(maxThrust - Drag)
            VxList.append(CAS)

        idx = excessThrustList.index(max(excessThrustList))

        return VxList[idx]

    def maxAltitude(self, HLid, LG, M, DeltaTemp, mass, nz=1.0):
        """Computes the maximum altitude taking into account buffet
        limitations. The altitude is calculated based on the aerodynamic
        configuration and the available buffet boundary conditions.

        :param HLid: High-lift devices position [-].
        :param LG: Landing gear position, either "LGUP" (gear up) or "LGDN"
            (gear down) [-].
        :param M: Mach airspeed [-].
        :param mass: Aircraft mass [kg].
        :param nz: Load factor [-], default is 1.0.
        :param DeltaTemp: Temperature deviation from ISA [K].
        :type HLid: float.
        :type LG: str.
        :type M: float.
        :type mass: float.
        :type nz: float.
        :type DeltaTemp: float.
        :returns: Maximum altitude [m].
        :rtype: float.
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
                    f,
                    x1=np.array([0]),
                    x2=np.array([conv.ft2m(hMax)]),
                    disp=False,
                )
            )

        return hMax

    def getConfig(self, phase, h, mass, v, DeltaTemp=0.0, hRWY=0.0):
        """Returns the aircraft's aerodynamic configuration based on altitude,
        speed, and phase of flight.

        :param phase: Phase of flight [Climb, Cruise, Descent] [-].
        :param h: Altitude above mean sea level (AMSL) [m].
        :param mass: Aircraft mass [kg].
        :param v: Calibrated airspeed (CAS) [m/s].
        :param DeltaTemp: Deviation from ISA temperature [K], default is
            0.0.
        :param hRWY: Runway elevation above mean sea level (AMSL) [m],
            default is 0.0.
        :type phase: str.
        :type h: float.
        :type mass: float.
        :type v: float.
        :type DeltaTemp: float.
        :type hRWY: float.
        :returns: Aircraft aerodynamic configuration [TO/IC/CR/AP/LD]
            [-].
        :rtype: str.
        :raises: TypeError if unable to determine the configuration.
        """

        config = None

        [theta, delta, sigma] = atm.atmosphereProperties(
            h=h, DeltaTemp=DeltaTemp
        )

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
            vMinCR = self.VMin(
                config="CR", mass=mass, theta=theta, delta=delta
            )
            vMinAPP = self.VMin(
                config="AP", mass=mass, theta=theta, delta=delta
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
        """Returns the aircraft aerodynamic configuration based on the
        provided configuration ID. This includes the high-lift device (HLID)
        position and landing gear (LG) position. If the configuration is not
        found, it returns None.

        :param config: Aircraft configuration (TO/IC/CR/AP/LD).
        :type config: str
        :returns: Aerodynamic configuration as a combination of HLID and LG.
        :rtype: [float, str]
        """

        configDict = self.AC.aeroConfig.get(config)

        return [configDict["HLid"], configDict["LG"]]

    def getSpeedSchedule(self, phase):
        """Returns the speed schedule based on the phase of flight (Climb,
        Cruise, Descent). The schedule includes two CAS values (CAS1, CAS2)
        and a Mach number (M).

        :param phase: Aircraft phase of flight (Climb, Cruise, Descent).
        :type phase: str
        :returns: Speed schedule as a combination of CAS1, CAS2 (in m/s) and
            Mach number (M).
        :rtype: [float, float, float]
        """

        speedScheduleDict = self.AC.speedSchedule[phase]

        return [
            conv.kt2ms(speedScheduleDict["CAS1"]),
            conv.kt2ms(speedScheduleDict["CAS2"]),
            speedScheduleDict["M"],
        ]

    def checkConfigurationContinuity(
        self, phase, previousConfig, currentConfig
    ):
        """Ensures the continuity of the aerodynamic configuration changes
        based on the phase of flight. It prevents sudden or improper
        configuration transitions, ensuring the aerodynamic configuration does
        not change in the wrong direction during Climb, Cruise, or Descent.

        :param phase: Aircraft phase of flight (Climb, Cruise, Descent).
        :param previousConfig: Previous aerodynamic configuration.
        :param currentConfig: Current aerodynamic configuration.
        :type phase: str
        :type previousConfig: str
        :type currentConfig: str
        :returns: The appropriate configuration for the current flight phase.
        :rtype: str
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
    """This class is a BADA4 aircraft subclass and implements the Airline
    Procedure Model (ARPM) following the BADA4 user manual.

    :param AC: Aircraft object {BADA4}.
    :type AC: bada4Aircraft.
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
        """Computes the climb speed schedule (CAS) for the given altitude
        based on various procedures and aircraft parameters.

        :param theta: Normalized air temperature [-].
        :param delta: Normalized air pressure [-].
        :param mass: Aircraft mass in kilograms [kg].
        :param h: Altitude in meters [m].
        :param hRWY: Runway elevation AMSL in meters [m].
        :param speedSchedule_default: Optional, a default speed schedule that overrides the BADA schedule. It should be in the form [Vcl1, Vcl2, Mcl].
        :param procedure: Climb procedure to be followed, e.g., 'BADA', 'NADP1', 'NADP2'. Default is 'BADA'.
        :param config: Optional, current aircraft aerodynamic configuration (TO/IC/CR/AP/LD).
        :param NADP1_ALT: Altitude in feet for NADP1 procedure. Default is 3000 feet.
        :param NADP2_ALT: Altitude range in feet for NADP2 procedure. Default is [1000, 3000].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :type theta: float
        :type delta: float
        :type mass: float
        :type h: float
        :type hRWY: float, optional
        :type speedSchedule_default: list[float], optional
        :type procedure: str
        :type config: str, optional
        :type NADP1_ALT: float, optional
        :type NADP2_ALT: list[float], optional
        :type DeltaTemp: float, optional
        :returns: A tuple containing the climb calibrated airspeed (CAS) in meters per second [m/s] and a status flag indicating whether the calculated CAS is constrained ('C'), unconstrained ('V' or 'v'), or not altered ('').
        :rtype: tuple[float, str]

        This function computes the climb speed schedule for different phases of flight and aircraft types.
        It supports BADA, NADP1, and NADP2 procedures for both jet and turboprop/piston/electric aircraft.

        The climb schedule uses specific speed profiles depending on altitude and aircraft model. For jet engines, the speed is constrained
        below 250 knots below 10,000 feet, and then it follows a defined speed schedule, either from BADA or NADP procedures.

        Additionally, the function applies speed limits based on the aircraft's flight envelope, adjusting the calculated climb speed if necessary.

        - For procedure='BADA', it uses the BADA climb speed schedule.
        - For procedure='NADP1', it implements the Noise Abatement Departure Procedure 1.
        - For procedure='NADP2', it implements the Noise Abatement Departure Procedure 2.

        The function also ensures that the calculated CAS remains within the bounds of the aircraft's minimum and maximum speeds.
        """

        # aircraft AGL altitude assuming being close to the RWY [m]
        h_AGL = h - hRWY

        phase = "Climb"
        acModel = self.AC.engineType

        [HLidTO, LGTO] = self.flightEnvelope.getAeroConfig(config="TO")
        VstallTO = self.flightEnvelope.VStall(
            h=h_AGL,
            mass=mass,
            HLid=HLidTO,
            LG=LGTO,
            nz=1.0,
            DeltaTemp=DeltaTemp,
        )

        [HLidCR, LGCR] = self.flightEnvelope.getAeroConfig(config="CR")
        VstallCR = self.flightEnvelope.VStall(
            h=h_AGL,
            mass=mass,
            HLid=HLidCR,
            LG=LGCR,
            nz=1.0,
            DeltaTemp=DeltaTemp,
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
                speed.append(
                    self.AC.CVmin * VstallTO + conv.kt2ms(self.AC.V_cl[5])
                )
                speed.append(
                    self.AC.CVmin * VstallTO + conv.kt2ms(self.AC.V_cl[4])
                )
                speed.append(
                    self.AC.CVmin * VstallTO + conv.kt2ms(self.AC.V_cl[3])
                )
                speed.append(
                    self.AC.CVmin * VstallTO + conv.kt2ms(self.AC.V_cl[2])
                )
                speed.append(
                    self.AC.CVmin * VstallTO + conv.kt2ms(self.AC.V_cl[1])
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
                    sigma = atm.sigma(theta=theta, delta=delta)
                    cas = atm.mach2Cas(
                        Mach=Mcl, theta=theta, delta=delta, sigma=sigma
                    )

            elif acModel == "TURBOPROP" or acModel == "PISTON":
                speed = list()
                speed.append(min(Vcl1, conv.kt2ms(250)))
                speed.append(
                    self.AC.CVmin * VstallTO + conv.kt2ms(self.AC.V_cl[8])
                )
                speed.append(
                    self.AC.CVmin * VstallTO + conv.kt2ms(self.AC.V_cl[7])
                )
                speed.append(
                    self.AC.CVmin * VstallTO + conv.kt2ms(self.AC.V_cl[6])
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
                    sigma = atm.sigma(theta=theta, delta=delta)
                    cas = atm.mach2Cas(
                        Mach=Mcl, theta=theta, delta=delta, sigma=sigma
                    )

        elif procedure == "NADP1":
            if acModel == "JET":
                speed = list()
                speed.append(min(Vcl1, conv.kt2ms(250)))
                speed.append(
                    self.AC.CVminTO * VstallTO + conv.kt2ms(self.AC.V_cl[2])
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
                    self.AC.CVminTO * VstallTO + conv.kt2ms(self.AC.V_cl[1])
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
                    self.AC.CVmin * VstallCR + conv.kt2ms(self.AC.V_cl[2])
                )
                speed.append(
                    self.AC.CVminTO * VstallTO + conv.kt2ms(self.AC.V_cl[2])
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
                    self.AC.CVmin * VstallCR + conv.kt2ms(self.AC.V_cl[2])
                )
                speed.append(
                    self.AC.CVminTO * VstallTO + conv.kt2ms(self.AC.V_cl[1])
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

        # check if the speed is within the limits of minimum and maximum speed from the flight envelope, if not, overwrite calculated speed with flight envelope min/max speed
        if config is None:
            config = self.flightEnvelope.getConfig(
                h=h,
                phase=phase,
                v=cas,
                mass=mass,
                DeltaTemp=DeltaTemp,
                hRWY=hRWY,
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
        """Computes the cruise speed schedule (CAS) for the given altitude
        based on various aircraft parameters.

        :param theta: Normalized air temperature [-].
        :param delta: Normalized air pressure [-].
        :param mass: Aircraft mass in kilograms [kg].
        :param h: Altitude in meters [m].
        :param hRWY: Runway elevation AMSL in meters [m].
        :param speedSchedule_default: Optional, a default speed schedule
            that overrides the BADA schedule. It should be in the form
            [Vcr1, Vcr2, Mcr].
        :param config: Optional, current aircraft aerodynamic
            configuration (TO/IC/CR/AP/LD).
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :type theta: float
        :type delta: float
        :type mass: float
        :type h: float
        :type hRWY: float, optional
        :type speedSchedule_default: list[float], optional
        :type config: str, optional
        :type DeltaTemp: float, optional
        :returns: A tuple containing the cruise calibrated airspeed
            (CAS) in meters per second [m/s] and a status flag
            indicating whether the calculated CAS is constrained ('C'),
            unconstrained ('V' or 'v'), or not altered ('').
        :rtype: tuple[float, str] This function computes the cruise
            speed schedule for different phases of flight and aircraft
            types. It uses either the default speed schedule or the BADA
            speed schedule based on the aircraft model and altitude.
            The cruise speed schedule varies depending on the altitude
            and type of engine (JET, TURBOPROP, or PISTON). The
            function also applies speed limits based on the aircraft's
            flight envelope, ensuring the calculated cruise speed
            remains within the aircraft's minimum and maximum allowable
            speeds. The function ensures the calculated CAS remains
            within the aircraft's operational speed limits, adjusting
            the speed if necessary.
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
                cas = atm.mach2Cas(
                    Mach=Mcr, theta=theta, delta=delta, sigma=sigma
                )

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
                cas = atm.mach2Cas(
                    Mach=Mcr, theta=theta, delta=delta, sigma=sigma
                )

        # check if the speed is within the limits of minimum and maximum speed from the flight envelope, if not, overwrite calculated speed with flight envelope min/max speed
        if config is None:
            config = self.flightEnvelope.getConfig(
                h=h,
                phase=phase,
                v=cas,
                mass=mass,
                DeltaTemp=DeltaTemp,
                hRWY=hRWY,
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
        """Computes the descent speed schedule (CAS) for the given altitude
        based on various aircraft parameters.

        :param theta: Normalized air temperature [-].
        :param delta: Normalized air pressure [-].
        :param mass: Aircraft mass in kilograms [kg].
        :param h: Altitude in meters [m].
        :param hRWY: Runway elevation AMSL in meters [m].
        :param speedSchedule_default: Optional, a default speed schedule
            that overrides the BADA schedule. It should be in the form
            [Vdes1, Vdes2, Mdes].
        :param config: Optional, current aircraft aerodynamic
            configuration (TO/IC/CR/AP/LD).
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :type theta: float
        :type delta: float
        :type mass: float
        :type h: float
        :type hRWY: float, optional
        :type speedSchedule_default: list[float], optional
        :type config: str, optional
        :type DeltaTemp: float, optional
        :returns: A tuple containing the descent calibrated airspeed
            (CAS) in meters per second [m/s] and a status flag
            indicating whether the calculated CAS is constrained ('C'),
            unconstrained ('V' or 'v'), or not altered ('').
        :rtype: tuple[float, str] This function computes the descent
            speed schedule for different phases of flight and aircraft
            types. It uses either the default speed schedule or the BADA
            speed schedule based on the aircraft model and altitude.
            The descent speed schedule varies depending on the altitude
            and type of engine (JET, TURBOPROP, or PISTON). The
            function ensures the calculated CAS remains within the
            aircraft's operational speed limits, adjusting the speed if
            necessary.
        """

        # aircraft AGL altitude assuming being close to the RWY [m]
        h_AGL = h - hRWY

        phase = "Descent"
        acModel = self.AC.engineType

        [HLid, LG] = self.flightEnvelope.getAeroConfig(config="LD")
        VstallDES = self.flightEnvelope.VStall(
            h=h_AGL, mass=mass, HLid=HLid, LG=LG, nz=1.0, DeltaTemp=DeltaTemp
        )
        [Vdes1, Vdes2, Mdes] = self.flightEnvelope.getSpeedSchedule(
            phase=phase
        )

        if speedSchedule_default is not None:
            Vdes1 = speedSchedule_default[0]
            Vdes2 = speedSchedule_default[1]
            Mdes = speedSchedule_default[2]

        crossOverAlt = atm.crossOver(cas=Vdes2, Mach=Mdes)

        if acModel == "JET" or acModel == "TURBOPROP":
            speed = []
            speed.append(min(Vdes1, conv.kt2ms(220)))
            speed.append(
                self.AC.CVmin * VstallDES + conv.kt2ms(self.AC.V_des[4])
            )
            speed.append(
                self.AC.CVmin * VstallDES + conv.kt2ms(self.AC.V_des[3])
            )
            speed.append(
                self.AC.CVmin * VstallDES + conv.kt2ms(self.AC.V_des[2])
            )
            speed.append(
                self.AC.CVmin * VstallDES + conv.kt2ms(self.AC.V_des[1])
            )

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
                cas = atm.mach2Cas(
                    Mach=Mdes, theta=theta, delta=delta, sigma=sigma
                )

        elif acModel == "PISTON":
            speed = list()
            speed.append(Vdes1)
            speed.append(
                self.AC.CVmin * VstallDES + conv.kt2ms(self.AC.V_des[7])
            )
            speed.append(
                self.AC.CVmin * VstallDES + conv.kt2ms(self.AC.V_des[6])
            )
            speed.append(
                self.AC.CVmin * VstallDES + conv.kt2ms(self.AC.V_des[5])
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
                sigma = atm.sigma(theta=theta, delta=delta)
                cas = atm.mach2Cas(
                    Mach=Mdes, theta=theta, delta=delta, sigma=sigma
                )

        # check if the speed is within the limits of minimum and maximum speed from the flight envelope, if not, overwrite calculated speed with flight envelope min/max speed
        if config is None:
            config = self.flightEnvelope.getConfig(
                h=h,
                phase=phase,
                v=cas,
                mass=mass,
                DeltaTemp=DeltaTemp,
                hRWY=hRWY,
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
    """This class implements the BADA4 optimization following the BADA4
    manual.

    :param AC: Aircraft object {BADA4}.
    :type AC: bada4Aircraft.
    """

    def __init__(self, AC):
        super().__init__(AC)
        self.flightEnvelope = FlightEnvelope(AC)

    def CCI(self, theta, delta, cI):
        """Computes the cost index coefficient for given flight conditions.

        :param cI: Cost index in kilograms per minute [kg min^-1].
        :param delta: Normalized pressure [-].
        :param theta: Normalized temperature [-].
        :type cI: float
        :type delta: float
        :type theta: float
        :returns: Cost index coefficient [-].
        :rtype: float
        """

        if self.AC.BADAVersion == "4.2":
            return (
                (cI / 60.0)
                * self.AC.LHV
                / (self.AC.MTOW * delta * const.g * const.a_0 * sqrt(theta))
            )
        elif self.AC.BADAVersion == "4.3" or self.AC.BADAVersion == "DUMMY":
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
        """Computes the weight coefficient at a given mass and pressure.

        :param mass: Aircraft mass in kilograms [kg].
        :param delta: Normalized pressure [-].
        :type mass: float
        :type delta: float
        :returns: Weight coefficient [-].
        :rtype: float The weight coefficient is used to represent the
            aircraft's weight relative to its maximum takeoff weight (MTOW)
            under given atmospheric conditions.
        """

        return mass * const.g / (self.AC.MTOW * delta * const.g)

    def SR(self, M, CF):
        """Computes the specific range (SR) for given flight conditions.

        :param M: Mach ground speed [-].
        :param CF: Fuel coefficient [-].
        :type M: float
        :type CF: float
        :returns: Specific range in nautical miles per kilogram [NM kg^-1].
        :rtype: float Specific range is a measure of the distance that can be
            flown per unit of fuel mass. It is calculated as the ratio of Mach
            number to fuel flow coefficient.
        """

        return M / CF

    def econMach(self, theta, delta, mass, DeltaTemp, cI, wS):
        """Computes the economic Mach number for a given flight condition and
        cost index.

        :param delta: Normalized pressure [-].
        :param theta: Normalized temperature [-].
        :param mass: Aircraft mass in kilograms [kg].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param cI: Cost index in kilograms per minute [kg min^-1].
        :param wS: Longitudinal wind speed in meters per second [m/s].
        :type delta: float
        :type theta: float
        :type mass: float
        :type DeltaTemp: float
        :type cI: float
        :type wS: float
        :returns: Maximum Range Cruise (MRC) Mach number [-].
        :rtype: float

        This function calculates the economic Mach number by balancing the fuel consumption and time costs,
        taking into account the aircraft’s flight conditions and the provided cost index (cI).
        It iterates through possible Mach numbers within the aircraft's flight envelope (bounded by Mmin and Mmax),
        computing the drag, thrust, and fuel flow, then selects the Mach that maximizes the cost efficiency.

        - **Mmin**: Lower bound of Mach speed, limited by buffet constraints.
        - **Mmax**: Upper bound of Mach speed, limited by buffet constraints.
        - **SR**: Specific Range, used to calculate the most efficient cruise speed in terms of fuel consumption.
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
        Mmax = self.flightEnvelope.maxMbuffet(
            delta=delta, mass=mass, HLid=HLid, LG=LG
        )

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
                rating="MCRZ",
                delta=delta,
                theta=theta,
                M=M,
                DeltaTemp=DeltaTemp,
            )

            # max Thrust limitation
            if Thrust > ThrustMax:
                continue

            CT = self.CT(Thrust=Thrust, delta=delta)
            CF = self.CF(
                CT=CT, delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
            )

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
        """Computes the Mach number representing Maximum Range Cruise (MRC)
        for the given flight conditions.

        :param theta: Normalized air temperature [-].
        :param delta: Normalized air pressure [-].
        :param mass: Aircraft mass in kilograms [kg].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param wS: Longitudinal wind speed in meters per second [m/s].
        :type theta: float
        :type delta: float
        :type mass: float
        :type DeltaTemp: float
        :type wS: float
        :returns: Maximum Range Cruise (MRC) Mach number [-].
        :rtype: float

        This function calculates the Mach number that corresponds to the Maximum Range Cruise (MRC),
        which maximizes the specific range. The calculation assumes clean configuration during the cruise phase.
        It uses the `econMach` function with a cost index (cI) of zero to find the MRC.

        If the calculated MRC is invalid, it returns NaN.
        """

        mrcM = self.econMach(
            cI=0.0,
            theta=theta,
            delta=delta,
            mass=mass,
            DeltaTemp=DeltaTemp,
            wS=wS,
        )

        if isnan(mrcM):
            return float("Nan")

        return mrcM

    def LRC(self, theta, delta, mass, DeltaTemp, wS):
        """Computes the Mach number representing Long Range Cruise (LRC) for
        the given flight conditions.

        :param theta: Normalized air temperature [-].
        :param delta: Normalized air pressure [-].
        :param mass: Aircraft mass in kilograms [kg].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param wS: Longitudinal wind speed in meters per second [m/s].
        :type theta: float
        :type delta: float
        :type mass: float
        :type DeltaTemp: float
        :type wS: float
        :returns: Long Range Cruise (LRC) Mach number [-].
        :rtype: float The Long Range Cruise (LRC) is defined as the speed
            where fuel efficiency is 99% of the Maximum Range Cruise (MRC).
            This function calculates the LRC based on the MRC and iterates
            through possible Mach numbers to find the one that minimizes the
            difference between the specific range (SR) at LRC and 99% of the
            SR at MRC. If no valid LRC is found, it returns NaN.
        """

        Mws = atm.tas2Mach(v=wS, theta=theta)

        MRC = self.MRC(
            theta=theta, delta=delta, mass=mass, DeltaTemp=DeltaTemp, wS=wS
        )

        if isnan(MRC):
            return float("Nan")

        # clean configuration during CR
        HLid = 0
        LG = "LGUP"

        CL = self.CL(M=MRC, delta=delta, mass=mass)
        CD = self.CD(M=MRC, CL=CL, HLid=HLid, LG=LG)
        Drag = self.D(M=MRC, delta=delta, CD=CD)
        CT = self.CT(Thrust=Drag, delta=delta)
        CF = self.CF(
            CT=CT, delta=delta, theta=theta, M=MRC, DeltaTemp=DeltaTemp
        )
        SR_LRC = 0.99 * self.SR(M=MRC + Mws, CF=CF)

        # min/max M speed limitation
        Mmax = self.flightEnvelope.maxMbuffet(
            delta=delta, mass=mass, HLid=HLid, LG=LG
        )

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
                rating="MCRZ",
                delta=delta,
                theta=theta,
                M=M,
                DeltaTemp=DeltaTemp,
            )

            # max Thrust limitation
            if Thrust > ThrustMax:
                continue

            if CL > CL_max:
                continue

            CT = self.CT(Thrust=Thrust, delta=delta)
            CF = self.CF(
                CT=CT, delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
            )

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
        """Computes the Mach number representing Maximum Endurance Cruise
        (MEC) for given flight conditions.

        :param delta: Normalized pressure [-].
        :param theta: Normalized temperature [-].
        :param mass: Aircraft mass in kilograms [kg].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param wS: Longitudinal wind speed in meters per second [m/s].
        :type delta: float
        :type theta: float
        :type mass: float
        :type DeltaTemp: float
        :type wS: float
        :returns: Maximum Endurance Cruise (MEC) in Mach [-].
        :rtype: float The Maximum Endurance Cruise (MEC) Mach is the Mach
            number that minimizes the fuel consumption rate, maximizing the
            endurance of the flight. This function iterates over a range of
            possible Mach numbers within the flight envelope and returns the
            Mach number with the lowest fuel coefficient (CF). The calculation
            is subject to thrust limitations, and the function ensures that
            the resulting Mach does not exceed the maximum thrust available.
        """

        # clean configuration during CR
        HLid = 0
        LG = "LGUP"

        Mws = atm.tas2Mach(v=wS, theta=theta)

        # min/max M speed limitation
        Mmin = self.flightEnvelope.minMbuffet(
            theta=theta, delta=delta, mass=mass, HLid=HLid, LG=LG
        )
        Mmax = self.flightEnvelope.maxMbuffet(
            delta=delta, mass=mass, HLid=HLid, LG=LG
        )

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
                rating="MCRZ",
                delta=delta,
                theta=theta,
                M=M,
                DeltaTemp=DeltaTemp,
            )

            # max Thrust limitation
            if Thrust > ThrustMax:
                continue

            CT = self.CT(Thrust=Thrust, delta=delta)
            CF = self.CF(
                CT=CT, delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
            )

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
        """Computes the optimum altitude for a given flight condition and Mach
        number.

        :param M: Mach number [-].
        :param mass: Aircraft mass in kilograms [kg].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :type M: float
        :type mass: float
        :type DeltaTemp: float
        :returns: Optimum altitude in feet [ft].
        :rtype: float The optimum altitude is the altitude where the aircraft
            achieves the maximum efficiency for a given Mach number and mass,
            subject to thrust and buffet limitations. The function iterates
            over a range of altitudes and returns the one with the best fuel
            efficiency. The function also ensures that the optimum altitude is
            bounded by a minimum of 2000 feet.
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
                rating="MCRZ",
                delta=delta,
                theta=theta,
                M=M,
                DeltaTemp=DeltaTemp,
            )

            # max Thrust limitation
            if Thrust > ThrustMax:
                continue

            CT = self.CT(Thrust=Thrust, delta=delta)
            CF = self.CF(
                CT=CT, delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
            )
            ff = self.ff(
                CT=CT, delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
            )
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
        """Returns the value of an optimization (OPT) parameter from a BADA4
        OPT file for various flight conditions. The OPT file contains values
        for parameters like Long Range Cruise (LRC), Maximum Endurance Cruise
        (MEC), Maximum Range Cruise (MRC), ECON speed, or OPTALT altitude, and
        this function interpolates the appropriate value based on input
        variables.

        Note:
            The array used in this function is expected to be sorted, as per the design of OPT files.

        :param optParam: Name of the optimization parameter file {LRC, MEC, MRC, ECON, OPTALT}.
        :param var_1: First optimizing factor, typically a flight condition like mass or altitude.
        :param var_2: (Optional) Second optimizing factor if the parameter depends on more than one factor.
        :type optParam: str
        :type var_1: float
        :type var_2: float, optional
        :returns: Interpolated optimization parameter value or NaN if not found.
        :rtype: float
        """

        filename = os.path.join(
            self.AC.filePath,
            self.AC.acName,
            optParam + ".OPT",
        )

        def findNearest(value, array):
            """Returns the indices of the nearest value(s) in the array. If
            the value is lower/higher than the lowest/highest value in the
            array, only one index is returned. If the value is between two
            values, two closest indices (left and right) are returned.

            Note:
                The array is expected to be sorted.

            :param value: Value to be compared to the array elements.
            :param array: Sorted array of values to search.
            :type value: float
            :type array: list[float]
            :returns: A list containing one or two nearest indices.
            :rtype: list[int]
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
            """Parses a BADA4 OPT file and populates variables for later
            interpolation.

            :param filename: Path to the OPT file.
            :type filename: str
            """

            file = open(filename, "r")
            lines = file.readlines()

            self.tableTypes = lines[8].split(":")[1].strip()
            self.tableDimension = lines[10].split(":")[1].strip()

            self.var_1 = list()
            self.var_2 = list()
            self.var_3 = list()

            if self.tableTypes == "2D":
                self.tableDimensionColumns = int(
                    self.tableDimension.split("x")[0]
                )
                self.tableDimensionRows = int(
                    self.tableDimension.split("x")[1]
                )

                self.var_2 = [
                    float(i)
                    for i in list(
                        filter(
                            None, lines[13].split("|")[1].strip().split(" ")
                        )
                    )
                ]

                for j in range(16, 16 + self.tableDimensionRows, 1):
                    self.var_1.append(float(lines[j].split("|")[0].strip()))
                    self.var_3.extend(
                        [
                            float(i) if i != "---" else float("nan")
                            for i in list(
                                filter(
                                    None,
                                    lines[j].split("|")[1].strip().split(" "),
                                )
                            )
                        ]
                    )

            if self.tableTypes == "1D":
                self.tableDimensionColumns = int(
                    self.tableDimension.split("x")[1]
                )
                self.tableDimensionRows = int(
                    self.tableDimension.split("x")[0]
                )

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
                    nearestIdx_1 * (self.tableDimensionColumns)
                    + nearestIdx_2[0]
                ]
                varTemp_2 = self.var_3[
                    nearestIdx_1 * (self.tableDimensionColumns)
                    + nearestIdx_2[1]
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
                    nearestIdx_1[0] * (self.tableDimensionColumns)
                    + nearestIdx_2
                ]
                varTemp_2 = self.var_3[
                    nearestIdx_1[1] * (self.tableDimensionColumns)
                    + nearestIdx_2
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
                    nearestIdx_1[0] * (self.tableDimensionColumns)
                    + nearestIdx_2[0]
                ]
                varTemp_2 = self.var_3[
                    nearestIdx_1[0] * (self.tableDimensionColumns)
                    + nearestIdx_2[1]
                ]

                varTemp_3 = self.var_3[
                    nearestIdx_1[1] * (self.tableDimensionColumns)
                    + nearestIdx_2[0]
                ]
                varTemp_4 = self.var_3[
                    nearestIdx_1[1] * (self.tableDimensionColumns)
                    + nearestIdx_2[1]
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
    """This class implements the PTD file creator for BADA4 aircraft following
    BADA4 manual.

    :param AC: Aircraft object {BADA4}.
    :type AC: bada4Aircraft.
    """

    def __init__(self, AC):
        super().__init__(AC)

        self.flightEnvelope = FlightEnvelope(AC)
        self.ARPM = ARPM(AC)

    def create(self, DeltaTemp, saveToPath):
        """Creates the BADA4 Performance Table Data (PTD) file, calculating
        climb, cruise, and descent profiles for different mass levels and
        altitudes, and saves the output to the specified path.

        :param saveToPath: Directory path where the PTD file should be stored.
        :param DeltaTemp: Deviation from the ISA (International Standard Atmosphere) temperature in Kelvin [K].
        :type saveToPath: str
        :type DeltaTemp: float
        :returns: None
        :rtype: None

        The function generates PTD files by computing the performance for three mass levels:
        - 120% of the Operating Empty Weight (OEW),
        - OEW plus 70% of the difference between Maximum Take-Off Weight (MTOW) and OEW,
        - and Maximum Take-Off Weight (MTOW).

        For each mass level, climb, cruise, and descent performance is calculated at various altitudes
        (up to the maximum operating altitude or a limit of 51,000 feet). The data is then saved to a PTD file.
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
        """Saves climb, cruise, and descent performance data into a PTD
        (Performance Table Data) file.

        :param saveToPath: Directory path where the PTD file will be stored.
        :param CLList: List of climb performance data.
        :param CRList: List of cruise performance data.
        :param DESList: List of descent performance data.
        :param DeltaTemp: Deviation from ISA (International Standard
            Atmosphere) temperature in Kelvin [K].
        :type saveToPath: str
        :type CLList: list
        :type CRList: list
        :type DESList: list
        :type DeltaTemp: float
        :returns: None
        :rtype: None
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
        file.write(
            "=============================\n=============================\n\n"
        )
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
        """Calculates the BADA4 PTD (Performance Table Data) for the CLIMB
        phase of flight.

        :param mass: Aircraft mass in kilograms [kg].
        :param altitudeList: List of aircraft altitudes in feet [ft].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :type mass: float
        :type altitudeList: list of int
        :type DeltaTemp: float
        :returns: List of PTD CLIMB data.
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
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=H_m, DeltaTemp=DeltaTemp
            )
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
                    rating="MCMB",
                    delta=delta,
                    theta=theta,
                    M=M,
                    DeltaTemp=DeltaTemp,
                )
                * 60
            )
            Thrust = self.Thrust(
                rating="MCMB",
                delta=delta,
                theta=theta,
                M=M,
                DeltaTemp=DeltaTemp,
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

            temp_const = (theta * const.temp_0) / (
                theta * const.temp_0 - DeltaTemp
            )
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
        """Calculates the BADA4 PTD (Performance Table Data) for the DESCENT
        phase of flight.

        :param mass: Aircraft mass in kilograms [kg].
        :param altitudeList: List of aircraft altitudes in feet [ft].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :type mass: float
        :type altitudeList: list of int
        :type DeltaTemp: float
        :returns: List of PTD DESCENT data.
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

        [Vdes1, Vdes2, Mdes] = self.flightEnvelope.getSpeedSchedule(
            phase=phase
        )
        Vdes1 = min(Vdes1, conv.kt2ms(250))
        crossAlt = atm.crossOver(cas=Vdes2, Mach=Mdes)

        for h in reversed(altitudeList):
            H_m = conv.ft2m(h)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=H_m, DeltaTemp=DeltaTemp
            )
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
                rating="LIDL",
                delta=delta,
                theta=theta,
                M=M,
                DeltaTemp=DeltaTemp,
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
                    rating="LIDL",
                    delta=delta,
                    theta=theta,
                    M=M,
                    DeltaTemp=DeltaTemp,
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
            temp_const = (theta * const.temp_0) / (
                theta * const.temp_0 - DeltaTemp
            )
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
                ff = (
                    self.ff(
                        CT=CT,
                        delta=delta,
                        theta=theta,
                        M=M,
                        DeltaTemp=DeltaTemp,
                    )
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
        """Calculates the BADA4 PTD (Performance Table Data) for the CRUISE
        phase of flight.

        :param mass: Aircraft mass in kilograms [kg].
        :param altitudeList: List of aircraft altitudes in feet [ft].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :type mass: float
        :type altitudeList: list of int
        :type DeltaTemp: float
        :returns: List of PTD CRUISE data.
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
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=H_m, DeltaTemp=DeltaTemp
            )
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
                rating="MCRZ",
                delta=delta,
                theta=theta,
                M=M,
                DeltaTemp=DeltaTemp,
            )
            CT = self.CT(Thrust=Thrust, delta=delta)
            ff = (
                self.ff(
                    CT=CT, delta=delta, theta=theta, M=M, DeltaTemp=DeltaTemp
                )
                * 60
            )

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
    """This class implements the PTF file creator for BADA4 aircraft following
    BADA4 manual.

    :param AC: Aircraft object {BADA4}.
    :type AC: bada4Aircraft.
    """

    def __init__(self, AC):
        super().__init__(AC)

        self.flightEnvelope = FlightEnvelope(AC)
        self.ARPM = ARPM(AC)

    def create(self, DeltaTemp, saveToPath):
        """Creates the BADA4 PTF and saves it to the specified directory.

        :param saveToPath: Path to the directory where the PTF should be
            stored.
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :type saveToPath: str
        :type DeltaTemp: float
        :returns: None
        :rtype: None
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
        self,
        saveToPath,
        CLList,
        CRList,
        DESList,
        DeltaTemp,
        massList,
        altitudeList,
    ):
        """Saves the BADA4 performance data to a PTF format.

        :param saveToPath: Path to the directory where the PTF should be
            stored.
        :param CLList: List of PTD data for CLIMB.
        :param CRList: List of PTD data for CRUISE.
        :param DESList: List of PTD data for DESCENT.
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param massList: List of aircraft mass levels [kg].
        :param altitudeList: List of altitudes [ft].
        :type saveToPath: str
        :type CLList: list
        :type CRList: list
        :type DESList: list
        :type DeltaTemp: float
        :type massList: list(float)
        :returns: None
        :rtype: None This function formats and writes the climb, cruise, and
            descent data for different mass levels and altitudes into a .PTF
            file, adhering to the BADA4 performance file format.
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
        [Vcr1, Vcr2, Mcr] = self.flightEnvelope.getSpeedSchedule(
            phase="Cruise"
        )
        [Vdes1, Vdes2, Mdes] = self.flightEnvelope.getSpeedSchedule(
            phase="Descent"
        )

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
            "BADA PERFORMANCE FILE                                        %s\n\n"
            % (d3)
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
        """Calculates the BADA4 PTF for the CRUISE phase of flight.

        :param massList: List of aircraft mass levels in kilograms [kg].
        :param altitudeList: List of aircraft altitudes in feet [ft].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :type massList: list
        :type altitudeList: list of int
        :type DeltaTemp: float
        :returns: List of PTF CRUISE data.
        :rtype: list
        """

        TAS_CR_complet = []
        FF_CR_LO_complet = []
        FF_CR_NOM_complet = []
        FF_CR_HI_complet = []

        massNominal = massList[1]

        [Vcr1, Vcr2, Mcr] = self.flightEnvelope.getSpeedSchedule(
            phase="Cruise"
        )
        Vcr1 = min(Vcr1, conv.kt2ms(250))

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=H_m, DeltaTemp=DeltaTemp
            )
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
                    rating="MCRZ",
                    delta=delta,
                    theta=theta,
                    M=M,
                    DeltaTemp=DeltaTemp,
                )

                CL = self.flightEnvelope.CL(
                    delta=delta, mass=mass, M=M, nz=1.2
                )
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
                            CT=CT,
                            delta=delta,
                            theta=theta,
                            M=M,
                            DeltaTemp=DeltaTemp,
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

        CRList = [
            TAS_CR_complet,
            FF_CR_LO_complet,
            FF_CR_NOM_complet,
            FF_CR_HI_complet,
        ]

        return CRList

    def PTF_climb(self, massList, altitudeList, DeltaTemp):
        """Calculates the BADA4 PTF for the CLIMB phase of flight.

        :param massList: List of aircraft mass levels in kilograms [kg].
        :param altitudeList: List of aircraft altitudes in feet [ft].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :type massList: list
        :type altitudeList: list of int
        :type DeltaTemp: float
        :returns: List of PTF CLIMB data.
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
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=H_m, DeltaTemp=DeltaTemp
            )
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
                    rating="MCMB",
                    delta=delta,
                    theta=theta,
                    M=M,
                    DeltaTemp=DeltaTemp,
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

                    if config == "TO" and (
                        prevConf == "IC" or prevConf == "CR"
                    ):
                        config = prevConf
                    elif config == "IC" and prevConf == "CR":
                        config = prevConf

                [HLid, LG] = self.flightEnvelope.getAeroConfig(config=config)
                CL = self.CL(M=M, delta=delta, mass=mass)
                CD = self.CD(M=M, CL=CL, HLid=HLid, LG=LG)
                Drag = self.D(M=M, delta=delta, CD=CD)

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
        """Calculates the BADA4 PTF for the DESCENT phase of flight.

        :param massList: List of aircraft mass levels in kilograms [kg].
        :param altitudeList: List of aircraft altitudes in feet [ft].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :type massList: list
        :type altitudeList: list of int
        :type DeltaTemp: float
        :returns: List of PTF DESCENT data.
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

        [Vdes1, Vdes2, Mdes] = self.flightEnvelope.getSpeedSchedule(
            phase="Descent"
        )
        Vdes1 = min(Vdes1, conv.kt2ms(250))
        crossAlt = atm.crossOver(cas=Vdes2, Mach=Mdes)

        # for h in altitudeList:
        for h in reversed(altitudeList):
            H_m = conv.ft2m(h)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(
                h=H_m, DeltaTemp=DeltaTemp
            )
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
                    rating="LIDL",
                    delta=delta,
                    theta=theta,
                    M=M,
                    DeltaTemp=DeltaTemp,
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

                    if config == "CR" and (
                        prevConf == "AP" or prevConf == "LD"
                    ):
                        config = prevConf
                    elif config == "AP" and prevConf == "LD":
                        config = prevConf

                [HLid, LG] = self.flightEnvelope.getAeroConfig(config=config)
                CL = self.CL(M=M, delta=delta, mass=mass)
                CD = self.CD(M=M, CL=CL, HLid=HLid, LG=LG)
                Drag = self.D(M=M, delta=delta, CD=CD)

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
                    ROCD_gamma = (
                        sin(conv.deg2rad(gamma)) * tas * (1 / temp_const)
                    )
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
                            CT=CT,
                            delta=delta,
                            theta=theta,
                            M=M,
                            DeltaTemp=DeltaTemp,
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
    """This class encapsulates the BADA4 performance model for an aircraft,
    extending the BADA4 base class.

    :param badaVersion: The version of the BADA4 model being used.
    :param acName: The ICAO designation or name of the aircraft.
    :param filePath: (Optional) Path to the BADA4 XML file. If not provided, a
        default path is used.
    :param allData: (Optional) Dataframe containing pre-loaded aircraft data,
        typically used to initialize the aircraft parameters without needing
        to parse XML files.
    :type badaVersion: str
    :type acName: str
    :type filePath: str, optional
    :type allData: pd.DataFrame, optional This class initializes the
        aircraft's performance model using data from a dataframe or by reading
        from XML files in the BADA4 format.
    """

    def __init__(self, badaVersion, acName, filePath=None, allData=None):
        """Initializes the BADA4Aircraft class by loading aircraft-specific
        data.

        - If `allData` is provided and contains the aircraft's information, it will be used to
          initialize various parameters such as engine type, mass, thrust settings, and performance
          data.
        - If the aircraft is not found in `allData`, the class will search for the corresponding
          BADA4 XML file or synonym file (if applicable) in the specified or default file path.
        - Once the aircraft data is found, the class initializes various performance modules such
          as the flight envelope, aerodynamic model, and performance optimizations.

        :param badaVersion: Version of the BADA4 model (e.g., "4.2", "4.3").
        :param acName: ICAO aircraft designation or model name.
        :param filePath: (Optional) Custom file path to load the aircraft data. If not provided,
                         a default directory is used.
        :param allData: (Optional) Dataframe containing pre-loaded aircraft data for initialization.
        """

        super().__init__(self)

        self.BADAFamily = BadaFamily(BADA4=True)
        self.BADAFamilyName = "BADA4"
        self.BADAVersion = badaVersion
        self.acName = acName

        if filePath is None:
            self.filePath = configuration.getBadaVersionPath(
                badaFamily="BADA4", badaVersion=badaVersion
            )
        else:
            self.filePath = filePath

        # check if the aircraft is in the allData dataframe data
        if allData is not None and acName in allData["acName"].values:
            filtered_df = allData[allData["acName"] == acName]

            self.model = configuration.safe_get(filtered_df, "model", None)
            self.engineType = configuration.safe_get(
                filtered_df, "engineType", None
            )
            self.engines = configuration.safe_get(filtered_df, "engines", None)
            self.WTC = configuration.safe_get(filtered_df, "WTC", None)
            self.ICAO = configuration.safe_get(filtered_df, "ICAO", None)
            self.MREF = configuration.safe_get(filtered_df, "MREF", None)
            self.WREF = configuration.safe_get(filtered_df, "WREF", None)
            self.LHV = configuration.safe_get(filtered_df, "LHV", None)
            self.n_eng = configuration.safe_get(filtered_df, "n_eng", None)
            self.rho = configuration.safe_get(filtered_df, "rho", None)
            self.TFA = configuration.safe_get(filtered_df, "TFA", None)
            self.p_delta = configuration.safe_get(filtered_df, "p_delta", None)
            self.p_theta = configuration.safe_get(filtered_df, "p_theta", None)
            self.kink = configuration.safe_get(filtered_df, "kink", None)
            self.b = configuration.safe_get(filtered_df, "b", None)
            self.c = configuration.safe_get(filtered_df, "c", None)
            self.max_power = configuration.safe_get(
                filtered_df, "max_power", None
            )
            self.p = configuration.safe_get(filtered_df, "p", None)
            self.a = configuration.safe_get(filtered_df, "a", None)
            self.f = configuration.safe_get(filtered_df, "f", None)
            self.ti = configuration.safe_get(filtered_df, "ti", None)
            self.fi = configuration.safe_get(filtered_df, "fi", None)
            self.throttle = configuration.safe_get(
                filtered_df, "throttle", None
            )
            self.prop_dia = configuration.safe_get(
                filtered_df, "prop_dia", None
            )
            self.max_eff = configuration.safe_get(filtered_df, "max_eff", None)
            self.Hd_turbo = configuration.safe_get(
                filtered_df, "Hd_turbo", None
            )
            self.CPSFC = configuration.safe_get(filtered_df, "CPSFC", None)
            self.P = configuration.safe_get(filtered_df, "P", None)
            self.S = configuration.safe_get(filtered_df, "S", None)
            self.HLPosition = configuration.safe_get(
                filtered_df, "HLPosition", None
            )
            self.configName = configuration.safe_get(
                filtered_df, "configName", None
            )
            self.VFE = configuration.safe_get(filtered_df, "VFE", None)
            self.d = configuration.safe_get(filtered_df, "d", None)
            self.CL_max = configuration.safe_get(filtered_df, "CL_max", None)
            self.bf = configuration.safe_get(filtered_df, "bf", None)
            self.HLids = configuration.safe_get(filtered_df, "HLids", None)
            self.CL_clean = configuration.safe_get(
                filtered_df, "CL_clean", None
            )
            self.M_max = configuration.safe_get(filtered_df, "M_max", None)
            self.scalar = configuration.safe_get(filtered_df, "scalar", None)
            self.Mmin = configuration.safe_get(filtered_df, "Mmin", None)
            self.Mmax = configuration.safe_get(filtered_df, "Mmax", None)
            self.CL_Mach0 = configuration.safe_get(
                filtered_df, "CL_Mach0", None
            )
            self.MTOW = configuration.safe_get(filtered_df, "MTOW", None)
            self.OEW = configuration.safe_get(filtered_df, "OEW", None)
            self.MFL = configuration.safe_get(filtered_df, "MFL", None)
            self.MTW = configuration.safe_get(filtered_df, "MTW", None)
            self.MZFW = configuration.safe_get(filtered_df, "MZFW", None)
            self.MPL = configuration.safe_get(filtered_df, "MPL", None)
            self.MLW = configuration.safe_get(filtered_df, "MLW", None)
            self.hmo = configuration.safe_get(filtered_df, "hmo", None)
            self.mfa = configuration.safe_get(filtered_df, "mfa", None)
            self.MMO = configuration.safe_get(filtered_df, "MMO", None)
            self.MLE = configuration.safe_get(filtered_df, "MLE", None)
            self.VLE = configuration.safe_get(filtered_df, "VLE", None)
            self.VMO = configuration.safe_get(filtered_df, "VMO", None)
            self.span = configuration.safe_get(filtered_df, "span", None)
            self.length = configuration.safe_get(filtered_df, "length", None)
            self.aeroConfig = configuration.safe_get(
                filtered_df, "aeroConfig", None
            )
            self.speedSchedule = configuration.safe_get(
                filtered_df, "speedSchedule", None
            )

            # GPF data (temporary)
            self.CVminTO = configuration.safe_get(filtered_df, "CVminTO", None)
            self.CVmin = configuration.safe_get(filtered_df, "CVmin", None)
            self.HmaxPhase = configuration.safe_get(
                filtered_df, "HmaxPhase", None
            )
            self.V_des = configuration.safe_get(filtered_df, "V_des", None)
            self.V_cl = configuration.safe_get(filtered_df, "V_cl", None)

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
                self.filePath,
                "aircraft_model_default.xml",
            )

            if os.path.isfile(synonymFile):
                self.synonymFileAvailable = True

                # if SYNONYM exist - look for synonym based on defined acName
                self.SearchedACName = Parser.parseMappingFile(
                    filePath=self.filePath,
                    acName=acName,
                )

                # if cannot find - look for full name (in sub folder names) based on acName (may not be ICAO designator)
                if self.SearchedACName is None:
                    self.SearchedACName = acName
                else:
                    self.ACinSynonymFile = True

            else:
                # if it doesn't exist - look for full name (in sub folder names) based on acName (may not be ICAO designator)
                self.SearchedACName = acName

            acXmlFile = (
                os.path.join(
                    self.filePath,
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
                        acName=self.SearchedACName,
                    )
                    GPFDataframe = Parser.parseGPF(filePath=self.filePath)

                    combined_df = Parser.combineXML_GPF(
                        XMLDataFrame, GPFDataframe
                    )

                    self.acName = self.SearchedACName

                    self.model = configuration.safe_get(
                        combined_df, "model", None
                    )
                    self.engineType = configuration.safe_get(
                        combined_df, "engineType", None
                    )
                    self.engines = configuration.safe_get(
                        combined_df, "engines", None
                    )
                    self.WTC = configuration.safe_get(combined_df, "WTC", None)
                    self.ICAO = configuration.safe_get(
                        combined_df, "ICAO", None
                    )
                    self.MREF = configuration.safe_get(
                        combined_df, "MREF", None
                    )
                    self.WREF = configuration.safe_get(
                        combined_df, "WREF", None
                    )
                    self.LHV = configuration.safe_get(combined_df, "LHV", None)
                    self.n_eng = configuration.safe_get(
                        combined_df, "n_eng", None
                    )
                    self.rho = configuration.safe_get(combined_df, "rho", None)
                    self.TFA = configuration.safe_get(combined_df, "TFA", None)
                    self.p_delta = configuration.safe_get(
                        combined_df, "p_delta", None
                    )
                    self.p_theta = configuration.safe_get(
                        combined_df, "p_theta", None
                    )
                    self.kink = configuration.safe_get(
                        combined_df, "kink", None
                    )
                    self.b = configuration.safe_get(combined_df, "b", None)
                    self.c = configuration.safe_get(combined_df, "c", None)
                    self.max_power = configuration.safe_get(
                        combined_df, "max_power", None
                    )
                    self.p = configuration.safe_get(combined_df, "p", None)
                    self.a = configuration.safe_get(combined_df, "a", None)
                    self.f = configuration.safe_get(combined_df, "f", None)
                    self.ti = configuration.safe_get(combined_df, "ti", None)
                    self.fi = configuration.safe_get(combined_df, "fi", None)
                    self.throttle = configuration.safe_get(
                        combined_df, "throttle", None
                    )
                    self.prop_dia = configuration.safe_get(
                        combined_df, "prop_dia", None
                    )
                    self.max_eff = configuration.safe_get(
                        combined_df, "max_eff", None
                    )
                    self.Hd_turbo = configuration.safe_get(
                        combined_df, "Hd_turbo", None
                    )
                    self.CPSFC = configuration.safe_get(
                        combined_df, "CPSFC", None
                    )
                    self.P = configuration.safe_get(combined_df, "P", None)
                    self.S = configuration.safe_get(combined_df, "S", None)
                    self.HLPosition = configuration.safe_get(
                        combined_df, "HLPosition", None
                    )
                    self.configName = configuration.safe_get(
                        combined_df, "configName", None
                    )
                    self.VFE = configuration.safe_get(combined_df, "VFE", None)
                    self.d = configuration.safe_get(combined_df, "d", None)
                    self.CL_max = configuration.safe_get(
                        combined_df, "CL_max", None
                    )
                    self.bf = configuration.safe_get(combined_df, "bf", None)
                    self.HLids = configuration.safe_get(
                        combined_df, "HLids", None
                    )
                    self.CL_clean = configuration.safe_get(
                        combined_df, "CL_clean", None
                    )
                    self.M_max = configuration.safe_get(
                        combined_df, "M_max", None
                    )
                    self.scalar = configuration.safe_get(
                        combined_df, "scalar", None
                    )
                    self.Mmin = configuration.safe_get(
                        combined_df, "Mmin", None
                    )
                    self.Mmax = configuration.safe_get(
                        combined_df, "Mmax", None
                    )
                    self.CL_Mach0 = configuration.safe_get(
                        combined_df, "CL_Mach0", None
                    )
                    self.MTOW = configuration.safe_get(
                        combined_df, "MTOW", None
                    )
                    self.OEW = configuration.safe_get(combined_df, "OEW", None)
                    self.MFL = configuration.safe_get(combined_df, "MFL", None)
                    self.MTW = configuration.safe_get(combined_df, "MTW", None)
                    self.MZFW = configuration.safe_get(
                        combined_df, "MZFW", None
                    )
                    self.MPL = configuration.safe_get(combined_df, "MPL", None)
                    self.MLW = configuration.safe_get(combined_df, "MLW", None)
                    self.hmo = configuration.safe_get(combined_df, "hmo", None)
                    self.mfa = configuration.safe_get(combined_df, "mfa", None)
                    self.MMO = configuration.safe_get(combined_df, "MMO", None)
                    self.MLE = configuration.safe_get(combined_df, "MLE", None)
                    self.VLE = configuration.safe_get(combined_df, "VLE", None)
                    self.VMO = configuration.safe_get(combined_df, "VMO", None)
                    self.span = configuration.safe_get(
                        combined_df, "span", None
                    )
                    self.length = configuration.safe_get(
                        combined_df, "length", None
                    )
                    self.aeroConfig = configuration.safe_get(
                        combined_df, "aeroConfig", None
                    )
                    self.speedSchedule = configuration.safe_get(
                        combined_df, "speedSchedule", None
                    )

                    # GPF data (temporary)
                    self.CVminTO = configuration.safe_get(
                        combined_df, "CVminTO", None
                    )
                    self.CVmin = configuration.safe_get(
                        combined_df, "CVmin", None
                    )
                    self.HmaxPhase = configuration.safe_get(
                        combined_df, "HmaxPhase", None
                    )
                    self.V_des = configuration.safe_get(
                        combined_df, "V_des", None
                    )
                    self.V_cl = configuration.safe_get(
                        combined_df, "V_cl", None
                    )

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
