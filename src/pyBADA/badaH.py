# -*- coding: utf-8 -*-
"""
pyBADA
Generic BADAH aircraft performance module
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


import xml.etree.ElementTree as ET
from datetime import date
import os
import numpy as np
from math import sqrt, pow, pi, cos, asin, radians, isnan
import pandas as pd

from scipy.optimize import minimize, Bounds

from pyBADA import constants as const
from pyBADA import conversions as conv
from pyBADA import atmosphere as atm
from pyBADA import configuration as configuration
from pyBADA.aircraft import Helicopter, BadaFamily


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
    """This class implements the BADAH parsing mechanism to parse xml BADAH files.

    :param filePath: path to the folder with BADAH xml formatted files.
    :param acName: ICAO aircraft designation
    :type filePath: str.
    :type acName: str
    """

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
        """This function parses BADAH xml formatted file

        :param filePath: path to the BADAH xml formatted file.
        :type filePath: str.
        :raises: IOError
        """

        acXmlFile = (
            os.path.join(filePath, "BADAH", badaVersion, acName, acName) + ".xml"
        )

        try:
            tree = ET.parse(acXmlFile)
            root = tree.getroot()
        except:
            raise IOError(acXmlFile + " not found or in correct format")

        # Parse general aircraft data
        model = root.find("model").text  # aircraft model
        engineType = root.find("type").text  # aircraft engine type
        engines = root.find("engine").text  # engine

        ICAO_desig = {}  # ICAO designator and WTC
        ICAO = root.find("ICAO").find("designator").text
        WTC = root.find("ICAO").find("WTC").text

        # Parse Aerodynamic Forces Model
        AFM = root.find("AFM")  # get AFM

        MR_radius = float(AFM.find("MR_radius").text)  # Main rotor radius
        MR_Speed = float(AFM.find("MR_speed").text)  # omega_m

        CPreq = AFM.find("CPreq")
        cpr = []
        for i in CPreq.findall("cpr"):
            cpr.append(float(i.text))

        # Parse engine data
        PFM = root.find("PFM")  # get PFM

        n_eng = int(PFM.find("n_eng").text)  # number of engines

        TPM = PFM.find("TPM")  # get TPM
        P0 = float(TPM.find("P0").text)

        CF = TPM.find("CF")
        cf = []
        for i in CF.findall("cf"):
            cf.append(float(i.text))

        Pmax_ = {}
        cpa = {}
        # Maximum take-off (MTKF)
        MTKF = TPM.find("MTKF")
        Pmax_["MTKF"] = float(MTKF.find("Pmax").text)
        CPav = MTKF.find("CPav")
        cpa["MTKF"] = []
        for i in CPav.findall("cpa"):
            cpa["MTKF"].append(float(i.text))

        # Maximum continuous (MCNT)
        MCNT = TPM.find("MCNT")
        Pmax_["MCNT"] = float(MCNT.find("Pmax").text)
        CPav = MCNT.find("CPav")
        cpa["MCNT"] = []
        for i in CPav.findall("cpa"):
            cpa["MCNT"].append(float(i.text))

        # Parse Aircraft Limitation Model (ALM)
        ALM = root.find("ALM")  # get ALM
        hmo = float(ALM.find("GLM").find("hmo").text)
        vne = float(ALM.find("KLM").find("vne").text)
        MTOW = float(ALM.find("DLM").find("MTOW").text)
        OEW = float(ALM.find("DLM").find("OEW").text)
        MFL = float(ALM.find("DLM").find("MFL").text)

        MREF = 2 * (MTOW - OEW) / 3 + OEW
        MPL = None  # maximum payload weight

        # Single row dataframe
        data = {
            "acName": [acName],
            "model": [model],
            "engineType": [engineType],
            "engines": [engines],
            "ICAO": [ICAO],
            "WTC": [WTC],
            "MR_radius": [MR_radius],
            "MR_Speed": [MR_Speed],
            "cpr": [cpr],
            "n_eng": [n_eng],
            "P0": [P0],
            "cf": [cf],
            "Pmax_": [Pmax_],
            "cpa": [cpa],
            "hmo": [hmo],
            "vne": [vne],
            "MTOW": [MTOW],
            "OEW": [OEW],
            "MFL": [MFL],
            "MREF": [MREF],
            "MPL": [MPL],
        }
        df_single = pd.DataFrame(data)

        return df_single

    @staticmethod
    def readSynonym(filePath, badaVersion):
        """This function parses BADAH Synonym xml formatted file and stores
        a dictionary of code names and files to be used for that specific code name

        :param filePath: path to the BADAH Synonym xml formatted file.
        :type filePath: str.
        :raises: IOError
        """

        filename = os.path.join(filePath, "BADAH", badaVersion, "SYNONYM.xml")

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

        if acName in synonym_fileName:
            fileName = synonym_fileName[acName]
            return fileName
        else:
            return None

    @staticmethod
    def parseAll(badaVersion, filePath=None):
        """This function parses all BADAH xml formatted file and stores
        all data in the final dataframe containing all the BADA data.

        :param filePath: path to the BADAH Synonym xml formatted file.
        :type filePath: str.
        :raises: IOError
        """

        if filePath == None:
            filePath = configuration.getAircraftPath()
        else:
            filePath = filePath

        synonym_fileName = Parser.readSynonym(filePath, badaVersion)

        # get names of all the folders in the main BADA model folder to search for XML files
        folderPath = os.path.join(filePath, "BADAH", badaVersion)
        subfolders = Parser.list_subfolders(folderPath)

        merged_df = pd.DataFrame()

        if synonym_fileName:
            for synonym in synonym_fileName:
                file = synonym_fileName[synonym]

                if file in subfolders:
                    # parse the original XML of a model
                    df = Parser.parseXML(filePath, badaVersion, file)

                    # rename acName in the data frame to match the synonym model name
                    df.at[0, "acName"] = synonym

                    # Merge DataFrames
                    merged_df = pd.concat([merged_df, df], ignore_index=True)

        else:
            for file in subfolders:
                # Parse the original XML of a model
                df = Parser.parseXML(filePath, badaVersion, file)

                # Merge DataFrames
                merged_df = pd.concat([merged_df, df], ignore_index=True)

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


class BADAH(Helicopter):
    """This class implements the part of BADAH performance model that will be used in other classes following the BADAH manual.

    :param AC: parsed aircraft.
    :type AC: badaH.Parse.
    """

    def __init__(self, AC):
        super().__init__()
        self.AC = AC

    def mu(self, tas):
        """This function computes the advance ratio

        :param tas: true airspeed (TAS) [m/s].
        :param gamma: flight path angle [rad].
        :type tas: float.
        :type gamma: float.
        :return: mu: advance ratio [-].
        :rtype: float.
        """

        # mu = (tas * math.cos(gamma))/(self.AC.MR_Speed*self.AC.MR_radius) #TODO: apply gamma modification
        mu = tas / (self.AC.MR_Speed * self.AC.MR_radius)

        return mu

    def CT(self, mass, rho, phi):
        """This function computes the thrust coefficient

        :param m: aircraft mass [kg].
        :param rho: air density [kg/m³].
        :param phi: bank angle [deg].
        :type m: float.
        :type rho: float.
        :type phi: float.
        :return: thrust coefficient [-].
        :rtype: float.
        """

        CT = (mass * const.g) / (
            rho
            * pi
            * pow(self.AC.MR_radius, 2)
            * pow(self.AC.MR_Speed * self.AC.MR_radius, 2)
            * cos(radians(phi))
        )

        return CT

    def CPreq(self, mu, CT):
        """This function computes the power required coefficient

        :param mu: advance ratio [-].
        :param CT: thrust coefficient [-].
        :type mu: float
        :type CT: float
        :return: power required coefficient [-]
        :rtype: float
        """

        CPreq = (
            self.AC.cpr[0]
            + self.AC.cpr[1] * pow(mu, 2)
            + self.AC.cpr[2] * CT * sqrt(sqrt(pow(mu, 4) + pow(CT, 2)) - pow(mu, 2))
            + self.AC.cpr[3] * pow(mu, 3)
            + self.AC.cpr[4] * pow(CT, 2) * pow(mu, 3)
        )

        return CPreq

    def Preq(self, sigma, tas, mass, phi=0.0):
        """This function computes the power required

        :param sigma: Normalised density [-].
        :param tas: true airspeed (TAS) [m/s]
        :param gamma: flight path angle [rad]
        :param mass: aircraft mass [kg].
        :param phi: bank angle [deg].
        :type sigma: float.
        :type v: float.
        :type gamma: float.
        :type mass: float
        :type phi: float
        :returns: power required [W].
        :rtype: float
        """

        # gamma = checkArgument('gamma', **kwargs)

        rho = sigma * const.rho_0

        # mu = self.mu(tas=tas,gamma=gamma)
        mu = self.mu(tas=tas)
        CT = self.CT(mass=mass, rho=rho, phi=phi)
        CPreq = self.CPreq(mu=mu, CT=CT)
        Preq = (
            rho
            * pi
            * pow(self.AC.MR_radius, 2)
            * pow(self.AC.MR_Speed * self.AC.MR_radius, 3)
            * CPreq
        )

        return Preq

    def Peng_target(self, ROCD, mass, Preq, ESF, temp, DeltaTemp):
        """This function computes the engine power

        :param temp: atmoshpere temperature [K].
        :param DeltaTemp: ISA temperature deviation [K].
        :param ROCD: rate of climb [m s^-1]
        :param mass: aircraft mass [kg].
        :param Preq: required power [W].
        :param ESF: energy share factor [-].
        :type temp: float.
        :type DeltaTemp: float.
        :type ROCD: float.
        :type mass: float.
        :type Preq: float.
        :type ESF: float.
        :returns: Peng [W].
        :rtype: float
        """

        temp_const = temp / (temp - DeltaTemp)
        Peng_target = (ROCD / ESF) * mass * const.g * temp_const + Preq

        return Peng_target

    def CPav(self, rating, delta, theta):
        """This function computes the power available coefficient

        :param rating: throttle setting {MTKF,MCNT}.
        :param delta: normalised pressure [-].
        :param theta: normalised temperature [-].
        :param SOC: State of charge [%]
        :type rating: str.
        :type delta: float.
        :type theta: float.
        :type SOC: float.
        :return: power available coefficient.
        :rtype: [-].
        :raise: ValueError.
        """

        sigma = atm.sigma(delta=delta, theta=theta)

        if self.AC.engineType == "TURBOPROP":
            if rating not in self.AC.Pmax_.keys():
                raise ValueError("Unknown engine rating " + rating)

            CPav = (
                self.AC.cpa[rating][0]
                + self.AC.cpa[rating][1] * pow(delta, 0.5)
                + self.AC.cpa[rating][2] * delta
                + self.AC.cpa[rating][3] * pow(delta, 2)
                + self.AC.cpa[rating][4] * pow(delta, 3)
                + self.AC.cpa[rating][5] * pow(theta, 0.5)
                + self.AC.cpa[rating][6] * theta
                + self.AC.cpa[rating][7] * pow(theta, 2)
                + self.AC.cpa[rating][8] * pow(theta, 3)
                + self.AC.cpa[rating][9] * pow(sigma, -1)
                + self.AC.cpa[rating][10] * pow(sigma, -0.5)
                + self.AC.cpa[rating][11] * pow(sigma, 0.5)
                + self.AC.cpa[rating][12] * sigma
            )

        elif self.AC.engineType == "PISTON":
            # currently identical to TURBOPROP, but this is subject to change in future versions
            if rating not in self.AC.Pmax_.keys():
                raise ValueError("Unknown engine rating " + rating)

            CPav = (
                self.AC.cpa[rating][0]
                + self.AC.cpa[rating][1] * pow(delta, 0.5)
                + self.AC.cpa[rating][2] * delta
                + self.AC.cpa[rating][3] * pow(delta, 2)
                + self.AC.cpa[rating][4] * pow(delta, 3)
                + self.AC.cpa[rating][5] * pow(theta, 0.5)
                + self.AC.cpa[rating][6] * theta
                + self.AC.cpa[rating][7] * pow(theta, 2)
                + self.AC.cpa[rating][8] * pow(theta, 3)
                + self.AC.cpa[rating][9] * pow(sigma, -1)
                + self.AC.cpa[rating][10] * pow(sigma, -0.5)
                + self.AC.cpa[rating][11] * pow(sigma, 0.5)
                + self.AC.cpa[rating][12] * sigma
            )

        else:
            raise ValueError("Unknown engine type")

        return CPav

    def Pmax(self, rating):
        """This function computes the maximum all-engine power [W]

        :param rating: throttle setting {MTKF,MCNT}.
        :type rating: str.
        :return: maximum all-engine power.
        :rtype: float.
        :raise: ValueError.
        """

        if rating not in self.AC.Pmax_.keys():
            raise ValueError("Unknown engine rating " + rating)
        return self.AC.Pmax_[rating]

    def Pav(self, rating, delta, theta):
        """This function computes the power available

        :param rating: throttle setting {MTKF,MCNT}.
        :param delta: normalised pressure [-].
        :param theta: normalised temperature [-].
        :param SOC:
        :type rating: str.
        :type delta: float.
        :type theta: float.
        :type SOC: float.
        :return: power available.
        :rtype: [W].
        :raise: ValueError.
        """

        Pmax = self.Pmax(rating=rating)

        CPav = self.CPav(rating=rating, delta=delta, theta=theta)

        Pav = min(
            Pmax,
            const.rho_0
            * pi
            * pow(self.AC.MR_radius, 2)
            * pow(self.AC.MR_Speed * self.AC.MR_radius, 3)
            * CPav,
        )

        return Pav

    def Q(self, Peng):
        """This function computes the torque value (expressed in percentage of of a reference torque)

        :param Peng: all-engine power [W].
        :type Peng: float.
        :return: torque value [%].
        :rtype: float.
        """

        Q = Peng / self.AC.P0

        return Q

    def CP(self, Peng):
        """This function computes the engine power coefficient

        :param Peng: all-engine power [W].
        :type Peng: float.
        :return: engine power coefficent [-].
        :rtype: float.
        """

        CP = Peng / (
            const.rho_0
            * pi
            * pow(self.AC.MR_radius, 2)
            * pow(self.AC.MR_Speed * self.AC.MR_radius, 3)
        )

        return CP

    def ff(self, delta, CP):
        """This function computes the fuel flow

        :param delta: Normalised pressure [-].
        :param CP: power coefficient [-].
        :type delta: float.
        :type CP: float.
        :return: fuel flow [kg/s]
        :rtype: float
        """

        if self.AC.engineType == "TURBOPROP":
            ff = (
                self.AC.cf[0]
                + self.AC.cf[1] * delta
                + self.AC.cf[2] * CP
                + self.AC.cf[3] * delta * CP
            )

        elif self.AC.engineType == "PISTON":
            # currently identical to TURBOPROP, but this is subject to change in future versions
            ff = (
                self.AC.cf[0]
                + self.AC.cf[1] * delta
                + self.AC.cf[2] * CP
                + self.AC.cf[3] * delta * CP
            )

        elif self.AC.engineType == "ELECTRIC":
            ff = 0.0

        else:
            raise ValueError("Unknown engine type")

        return ff / 3600

    def ROCD(self, Peng, Preq, mass, ESF, theta, DeltaTemp):
        """This function computes the Rate of Climb or Descent

        :param theta: normalised temperature [-].
        :param Peng: Peng: all-engine power [W].
        :param Preq: power required [W].
        :param mass: actual aircraft mass  [kg].
        :param ESF: energy share factor [-].
        :param DeltaTemp: deviation with respect to ISA [K]
        :type theta: float.
        :type Peng: float.
        :type Preq: float.
        :type mass: float.
        :type ESF: float.
        :type DeltaTemp: float.
        :returns: rate of climb/descend [m/s].
        :rtype: float
        """

        temp = theta * const.temp_0
        ROCD = ((temp - DeltaTemp) / temp) * (Peng - Preq) * ESF / (mass * const.g)

        return ROCD


class FlightEnvelope(BADAH):
    """This class is a BADAH aircraft subclass and implements the flight envelope caclulations
    following the BADAH manual.

    :param AC: parsed aircraft.
    :type AC: badaH.Parse.
    """

    def __init__(self, AC):
        super().__init__(AC)

    def maxAltitude(self):
        """This function computes the maximum altitude

        :returns: maximum altitude [m].
        :rtype: float
        """

        hMax = conv.ft2m(self.AC.hmo)
        return hMax

    def VMax(self):
        """This function computes the maximum speed

        :returns: maximum CAS speed [m s^-1].
        :rtype: float.
        """

        Vmax = conv.kt2ms(self.AC.vne)
        return Vmax

    def speedEnvelope_powerLimited(
        self, h, mass, DeltaTemp, rating="MCNT", rateOfTurn=0
    ):
        """This function computes the maximum CAS speed within the certified flight envelope taking into account the trust limitation.

        :param h: altitude [m].
        :param mass: aircraft operating mass [kg]
        :param DeltaTemp: deviation with respect to ISA [K]
        :param rating: aircraft engine rating [MTKF/MCNT][-]
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
        Pmax = self.Pav(rating=rating, theta=theta, delta=delta)
        Pmin = 0.1 * self.AC.P0  # No minimum power model: assume 10% torque

        VminCertified = 0
        VmaxCertified = self.VMax()

        CASlist = []
        for CAS in np.linspace(VminCertified, VmaxCertified, num=200, endpoint=True):
            [M, CAS, TAS] = atm.convertSpeed(
                v=conv.ms2kt(CAS),
                speedType="CAS",
                theta=theta,
                delta=delta,
                sigma=sigma,
            )

            bankAngle = self.bankAngle(rateOfTurn=rateOfTurn, v=TAS)

            Preq = self.Preq(sigma=sigma, tas=TAS, mass=mass, phi=bankAngle)

            if Pmax >= Preq:
                CASlist.append(CAS)

        if not CASlist:
            return (None, None)
        else:
            minCAS = min(CASlist)
            maxCAS = max(CASlist)

            return (minCAS, maxCAS)


class Optimization(BADAH):
    """This class implements the BADAH optimization following the BADAH manual.

    :param AC: parsed aircraft.
    :type AC: badaH.Parse.
    """

    def __init__(self, AC):
        super().__init__(AC)

        self.flightEnvelope = FlightEnvelope(AC)

    def MRC(self, h, mass, DeltaTemp, wS):
        """This function computes the TAS reperesenting Maximum Range Cruise (MRC) for given flight conditions

        :param h: altitude [m].
        :param mass: aircraft weight [kg].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param wS: longitudinal wind speed (TAS) [m s^-1].
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type wS: float.
        :returns: Maximum Range Cruise (MRC) in TAS [m s^-1]
        :rtype: float.
        """
        # NOTE: check for precision of algorithm needed. Possible local minima, instead of global minima

        [theta, delta, sigma] = atm.atmosphereProperties(
            h=h, DeltaTemp=DeltaTemp
        )  # atmosphere properties

        # max TAS speed limitation
        Vmax = atm.cas2Tas(cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma)

        epsilon = 0.1
        TAS_list = np.arange(0, Vmax + epsilon, epsilon)

        Pav = self.Pav(rating="MCNT", theta=theta, delta=delta)

        TAS_MRC = []
        cost_MRC = []

        # def f(TAS):
        #     Preq = self.Preq(sigma=sigma, tas=TAS[0], mass=mass)
        #     CP = self.CP(Peng=Preq)
        #     ff = self.ff(delta=delta, CP=CP)

        #     cost = -((TAS[0]+wS) / ff)

        # minimize cost function
        #     return cost

        # epsilon = 0.01
        # bnds = Bounds([0],[Vmax+epsilon])
        # Pav limitation -> Preq > Pav
        # cons = ({'type': 'ineq','fun': lambda TAS: Pav - self.Preq(sigma=sigma, tas=TAS[0], mass=mass)})
        # mrc = minimize(f, np.array([0.1]), method='SLSQP', bounds=bnds, constraints=cons).x

        for TAS in TAS_list:
            Preq = self.Preq(sigma=sigma, tas=TAS, mass=mass)

            CP = self.CP(Peng=Preq)
            ff = self.ff(delta=delta, CP=CP)

            # Pav limitation
            if Preq > Pav:
                continue

            # maximize the cost function
            cost = (TAS + wS) / ff

            TAS_MRC.append(TAS)
            cost_MRC.append(cost)

        if not cost_MRC:
            return float("Nan")

        mrc = TAS_MRC[cost_MRC.index(max(cost_MRC))]

        return mrc

    def LRC(self, h, mass, DeltaTemp, wS):
        """This function computes the TAS reperesenting Long Range Cruise (LRC) for given flight conditions

        :param h: altitude [m].
        :param mass: aircraft weight [kg].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param wS: longitudinal wind speed (TAS) [m s^-1].
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type wS: float.
        :returns: Long Range Cruise (LRC) in M [-]
        :rtype: float.
        """
        # NOTE: check for precision of algorithm needed. Possible local minima, instead of global minima

        MRC = self.MRC(mass=mass, h=h, DeltaTemp=DeltaTemp, wS=wS)

        if isnan(MRC):
            return float("Nan")

        [theta, delta, sigma] = atm.atmosphereProperties(
            h=h, DeltaTemp=DeltaTemp
        )  # atmosphere properties

        Preq = self.Preq(sigma=sigma, tas=MRC, mass=mass)
        CP = self.CP(Peng=Preq)
        ff = self.ff(delta=delta, CP=CP)
        SR = (MRC + wS) / ff
        SR_LRC = 0.99 * SR

        # max TAS speed limitation
        Vmax = atm.cas2Tas(cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma)
        Pav = self.Pav(rating="MCNT", theta=theta, delta=delta)

        # LRC > MRC
        epsilon = 0.001
        TAS_list = np.arange(MRC, Vmax + epsilon, epsilon)

        TAS_LRC = []
        cost_LRC = []

        for TAS in TAS_list:
            Preq = self.Preq(sigma=sigma, tas=TAS, mass=mass)

            CP = self.CP(Peng=Preq)
            ff = self.ff(delta=delta, CP=CP)

            # Pav limitation
            if Preq > Pav:
                continue

            SR = (TAS + wS) / ff

            # minimize the cost function
            cost_LRC.append(sqrt((SR - SR_LRC) ** 2))
            TAS_LRC.append(TAS)

        lrc = TAS_LRC[cost_LRC.index(min(cost_LRC))]

        # def f(TAS):
        #     Preq = self.Preq(sigma=sigma, tas=TAS[0], mass=mass)
        #     CP = self.CP(Peng=Preq)
        #     ff = self.ff(delta=delta, CP=CP)

        #     SR = (TAS[0]+wS) / ff
        #     cost_LRC = sqrt((SR - SR_LRC)**2)

        # minimize cost function
        #     return cost_LRC

        # epsilon = 0.01
        # LRC > MRC
        # bnds = Bounds([MRC],[Vmax+epsilon])
        # Pav limitation -> Preq > Pav
        # cons = ({'type': 'ineq','fun': lambda TAS: Pav - self.Preq(sigma=sigma, tas=TAS[0], mass=mass)})
        # lrc = minimize(f, np.array([300]), method='SLSQP', bounds=bnds, constraints=cons)

        # lrc = fmin(f, x0=np.array([MRC]), disp=False)

        return lrc

    def MEC(self, h, mass, DeltaTemp, wS):
        """This function computes the TAS reperesenting Maximum Endurance Cruise (MEC) for given flight conditions

        :param h: altitude [m].
        :param mass: aircraft weight [kg].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param wS: longitudinal wind speed (TAS) [m s^-1].
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type wS: float.
        :returns: Maximum Endurance Cruise (MEC) in TAS [m s^-1]
        :rtype: float.
        """

        [theta, delta, sigma] = atm.atmosphereProperties(
            h=h, DeltaTemp=DeltaTemp
        )  # atmosphere properties

        # max TAS speed limitation
        Vmax = atm.cas2Tas(cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma)

        # def f(TAS):
        # Preq = self.Preq(sigma=sigma, tas=TAS[0], mass=mass)
        # CP = self.CP(Peng=Preq)
        # ff = self.ff(delta=delta, CP=CP)

        # minimize ff -> const function
        # return ff

        # epsilon = 0.01
        # bnds = Bounds([0], [Vmax + epsilon])
        # mec = minimize(f, np.array([epsilon]), method="SLSQP", bounds=bnds).x

        epsilon = 0.01
        TAS_list = np.arange(0, Vmax + epsilon, epsilon)

        ff_mec = []
        TAS_mec = []
        for TAS in TAS_list:
            Preq = self.Preq(sigma=sigma, tas=TAS, mass=mass)
            CP = self.CP(Peng=Preq)
            ff = self.ff(delta=delta, CP=CP)

            # minimize the cost function
            ff_mec.append(ff)
            TAS_mec.append(TAS)

        if not ff_mec:
            return float("Nan")

        mecTAS = TAS_mec[ff_mec.index(min(ff_mec))]

        return mecTAS

    def parseOPT(self, filename):
        """This function parses BADA4 OPT ascii formatted files

        :param filename: path to the ___.OPT ascii formatted file.
        :type filename: str.
        :returns: dictionary of Delta temperature and data from the OPT file for that delta temperature [kt]
        :rtype: dict.
        """

        file = open(filename, "r")
        lines = file.readlines()

        DeltaTempPos = {}

        # create a dictionary for list of DeltaTemp available in OPT file mapped to the line number in the file
        for k in range(len(lines)):
            line = lines[k]
            if "DeltaT:" in line:
                DeltaTempPos[int(line.split(":")[1].strip())] = k

        self.tableTypes = lines[7].split(":")[1].strip()
        self.tableDimension = lines[9].split(":")[1].strip()

        DeltaTempDict = {}

        if self.tableTypes == "3D":
            self.tableDimensionColumns = int(self.tableDimension.split("x")[2])
            self.tableDimensionRows = int(self.tableDimension.split("x")[1])
            self.DeltaTempNum = int(self.tableDimension.split("x")[0])

            for DeltaTemp in DeltaTempPos:
                var_1 = []
                var_2 = []
                var_3 = []

                startIdx = DeltaTempPos[DeltaTemp] + 1
                var_2 = [
                    float(i)
                    for i in list(
                        filter(None, lines[startIdx].split("|")[1].strip().split(" "))
                    )
                ]

                for j in range(startIdx + 3, startIdx + 3 + self.tableDimensionRows, 1):
                    var_1.append(float(lines[j].split("|")[0].strip()))

                    str_list = list(
                        filter(None, lines[j].split("|")[1].strip().split(" "))
                    )
                    for k in range(len(str_list)):
                        if str_list[k] == "-":
                            str_list[k] = float("Nan")

                    var_3.extend([float(i) for i in str_list])

                DeltaTempDict[DeltaTemp] = [var_1, var_2, var_3]

        return DeltaTempDict

    def findNearestIdx(self, value, array):
        """This function returns indices of the nearest value in the array.
        if the value is lower/higher than lowest/highest value in array, only one idx is returned
        otherwise if the value is somewhere in between, 2 closest (left and right) idx are returned

        .. note::
                array used in this fuction is expected to be sorted (design of OPT files)

        :param value: value to which the array value will be comapred
        :param array: list of values
        :type value: float.
        :type array: array of float.
        :returns: nearest indices
        :rtype: list[float].
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

    def calculateOPTparam(self, var_1, var_2, detaTauList):
        """This function calculate the OPT value by either selecting the existing value, or interpolating between closest 2 values

        .. note::array used in this fuction is expected to be sorted (design of OPT files)

        :param var_1: value of the first optimizing factor.
        :param var_2: value of the second optimizing factor.
        :param detaTauList: list of values belonging to specified DeltaTemp from OPT file.
        :type var_1: float.
        :type var_2: float.
        :type detaTauList: list[float].
        """

        var_1_list = detaTauList[0]
        var_2_list = detaTauList[1]
        var_3_list = detaTauList[2]

        nearestIdx_1 = np.array(self.findNearestIdx(var_1, var_1_list))
        nearestIdx_2 = np.array(self.findNearestIdx(var_2, var_2_list))

        # if nearestIdx_1 & nearestIdx_2 [1] [1]
        if (nearestIdx_1.size == 1) & (nearestIdx_2.size == 1):
            return var_3_list[
                nearestIdx_1 * (self.tableDimensionColumns) + nearestIdx_2
            ]

        # if nearestIdx_1 & nearestIdx_2 [1] [1,2]
        if (nearestIdx_1.size == 1) & (nearestIdx_2.size == 2):
            varTemp_1 = var_3_list[
                nearestIdx_1 * (self.tableDimensionColumns) + nearestIdx_2[0]
            ]
            varTemp_2 = var_3_list[
                nearestIdx_1 * (self.tableDimensionColumns) + nearestIdx_2[1]
            ]

            # interpolation between the 2 found points
            interpVar = np.interp(
                var_2,
                [var_2_list[nearestIdx_2[0]], var_2_list[nearestIdx_2[1]]],
                [varTemp_1, varTemp_2],
            )
            return interpVar

        # if nearestIdx_1 & nearestIdx_2 [1,2] [1]
        if (nearestIdx_1.size == 2) & (nearestIdx_2.size == 1):
            varTemp_1 = var_3_list[
                nearestIdx_1[0] * (self.tableDimensionColumns) + nearestIdx_2
            ]
            varTemp_2 = var_3_list[
                nearestIdx_1[1] * (self.tableDimensionColumns) + nearestIdx_2
            ]

            # interpolation between the 2 found points
            interpVar = np.interp(
                var_1,
                [var_1_list[nearestIdx_1[0]], var_1_list[nearestIdx_1[1]]],
                [varTemp_1, varTemp_2],
            )
            return interpVar

        # if nearestIdx_1 & nearestIdx_2 [1,2] [1,2]
        if (nearestIdx_1.size == 2) & (nearestIdx_2.size == 2):
            varTemp_1 = var_3_list[
                nearestIdx_1[0] * (self.tableDimensionColumns) + nearestIdx_2[0]
            ]
            varTemp_2 = var_3_list[
                nearestIdx_1[0] * (self.tableDimensionColumns) + nearestIdx_2[1]
            ]

            varTemp_3 = var_3_list[
                nearestIdx_1[1] * (self.tableDimensionColumns) + nearestIdx_2[0]
            ]
            varTemp_4 = var_3_list[
                nearestIdx_1[1] * (self.tableDimensionColumns) + nearestIdx_2[1]
            ]

            # interpolation between the 4 found points
            interpVar_1 = np.interp(
                var_2,
                [var_2_list[nearestIdx_2[0]], var_2_list[nearestIdx_2[1]]],
                [varTemp_1, varTemp_2],
            )
            interpVar_2 = np.interp(
                var_2,
                [var_2_list[nearestIdx_2[0]], var_2_list[nearestIdx_2[1]]],
                [varTemp_3, varTemp_4],
            )
            interpVar_3 = np.interp(
                var_1,
                [var_1_list[nearestIdx_1[0]], var_1_list[nearestIdx_1[1]]],
                [interpVar_1, interpVar_2],
            )

            return interpVar_3

    def getOPTParam(self, optParam, var_1, var_2, DeltaTemp):
        """This function returns value of the OPT parameter based on the input value from OPT file
        like LRC, MEC, MRC

        .. note::
                array used in this fuction is expected to be sorted (design of BADA OPT files)

        :param optParam: name of optimization file {LRC,MEC,MRC}.
        :param var_1: value of the first optimizing factor.
        :param var_2: value of the second optimizing factor.
        :type optParam: string.
        :type var_1: float.
        :type var_2: float.
        """

        filename = os.path.join(
            self.AC.filePath,
            self.AC.BADAFamilyName,
            self.AC.BADAVersion,
            self.AC.acName,
            optParam + ".OPT",
        )

        detaTauDict = self.parseOPT(filename=filename)

        if DeltaTemp in detaTauDict:
            # value of DeltaTemp exist in the OPT file
            optVal = self.calculateOPTparam(var_1, var_2, detaTauDict[DeltaTemp])
        else:
            # value of DeltaTemp does not exist in OPT file - will be interpolated. But only within the range of <-20;20>
            nearestIdx = np.array(self.findNearestIdx(DeltaTemp, list(detaTauDict)))

            if nearestIdx.size == 1:
                # DeltaTemp value is either outside of the <-20;20> DeltaTemp range
                DeltaTemp_new = list(detaTauDict)[nearestIdx]
                optVal = self.calculateOPTparam(
                    var_1, var_2, detaTauDict[DeltaTemp_new]
                )
            else:
                # DeltaTemp value is within the <-20;20> DeltaTemp range
                # calculate the interpolation between 2 closest DeltaTemp values from the OPT file
                DeltaTemp_new_1 = list(detaTauDict)[nearestIdx[0]]
                DeltaTemp_new_2 = list(detaTauDict)[nearestIdx[1]]

                optVal_1 = self.calculateOPTparam(
                    var_1, var_2, detaTauDict[DeltaTemp_new_1]
                )
                optVal_2 = self.calculateOPTparam(
                    var_1, var_2, detaTauDict[DeltaTemp_new_2]
                )

                optVal = np.interp(
                    DeltaTemp, [DeltaTemp_new_1, DeltaTemp_new_2], [optVal_1, optVal_2]
                )

        return optVal


class ARPM(BADAH):
    """This class is a BADAH aircraft subclass and implements the Airline Procedure Model (ARPM)
    following the BADAH user manual.

    :param AC: parsed aircraft.
    :type AC: badaH.Parse.
    """

    def __init__(self, AC):
        super().__init__(AC)

        self.flightEnvelope = FlightEnvelope(AC)
        self.OPT = Optimization(AC)

    def takeoff(
        self, h, mass, DeltaTemp, rating="ARPM", speedLimit=None, ROCDDefault=None
    ):
        """This function computes parameters for the takeoff ARPM

        :param h: altitude [m].
        :param mass: aircraft weight [kg].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param rating: engine rating {MTKF,MCNT,ARPM} [-].
        :param speedLimit: decision to apply or not the speed limit {applyLimit,''} [-].
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type rating: string.
        :type speedLimit: string.
        :returns: [Pav, Peng, Preq, tas, ROCD, ESF, limitation] [W, W, W, m/s, m/s, -, -]
        :rtype: float.
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        delta = atm.delta(h=h, DeltaTemp=DeltaTemp)
        sigma = atm.sigma(theta=theta, delta=delta)

        temp = theta * const.temp_0

        # control parameters
        tas = 0
        if ROCDDefault is None:
            ROCD = conv.ft2m(100) / 60  # [m/s]
        else:
            ROCD = ROCDDefault

        # check for speed envelope limitations
        eps = 1e-6  # float calculation precision
        maxSpeed = atm.cas2Tas(cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma)
        minSpeed = 0
        limitation = ""

        # empty envelope - keep the original calculated TAS speed
        if maxSpeed < minSpeed:
            if (tas - eps) > maxSpeed and (tas - eps) > minSpeed:
                limitation = "V"
            elif (tas + eps) < minSpeed and (tas + eps) < maxSpeed:
                limitation = "v"
            else:
                limitation = "vV"

        elif minSpeed > (tas + eps):
            if speedLimit == "applyLimit":
                tas = minSpeed
                limitation = "C"
            else:
                limitation = "v"

        elif maxSpeed < (tas - eps):
            if speedLimit == "applyLimit":
                tas = maxSpeed
                limitation = "C"
            else:
                limitation = "V"

        ESF = self.esf(flightEvolution="constTAS")

        Preq = self.Preq(sigma=sigma, tas=tas, mass=mass)

        if rating == "ARPM":
            Peng_target = self.Peng_target(
                temp=temp, DeltaTemp=DeltaTemp, ROCD=ROCD, mass=mass, Preq=Preq, ESF=ESF
            )
            Pav = self.Pav(rating="MTKF", theta=theta, delta=delta)
            Peng = min(Peng_target, Pav)

            ROCD_TEM = self.ROCD(
                Peng=Peng,
                Preq=Preq,
                mass=mass,
                ESF=ESF,
                theta=theta,
                DeltaTemp=DeltaTemp,
            )

            if ROCD_TEM < ROCD:
                ROCD = ROCD_TEM

        elif rating == "MTKF":
            Pav = self.Pav(rating="MTKF", theta=theta, delta=delta)
            Peng = Pav
            ROCD = self.ROCD(
                Peng=Peng,
                Preq=Preq,
                mass=mass,
                ESF=ESF,
                theta=theta,
                DeltaTemp=DeltaTemp,
            )

        elif rating == "MCNT":
            Pav = self.Pav(rating="MCNT", theta=theta, delta=delta)
            Peng = Pav
            ROCD = self.ROCD(
                Peng=Peng,
                Preq=Preq,
                mass=mass,
                ESF=ESF,
                theta=theta,
                DeltaTemp=DeltaTemp,
            )

        if Pav < Peng:
            limitation += "P"

        return [Pav, Peng, Preq, tas, ROCD, ESF, limitation]

    def accelerationToClimb(self):
        pass

    def climb(
        self,
        h,
        mass,
        DeltaTemp,
        rating="ARPM",
        speedLimit=None,
        ROCDDefault=None,
        tasDefault=None,
    ):
        """This function computes parameters for the climb ARPM

        :param h: altitude [m].
        :param mass: aircraft weight [kg].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param wS: longitudinal wind speed (TAS) [m s^-1].
        :param rating: engine rating {MTKF,MCNT,ARPM} [-].
        :param speedLimit: decision to apply or not the speed limit {applyLimit,''} [-].
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type wS: float.
        :type rating: string.
        :type speedLimit: string.
        :returns: [Pav, Peng, Preq, tas, ROCD, ESF, limitation] [W, W, W, m/s, m/s, -, -]
        :rtype: float.
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        delta = atm.delta(h=h, DeltaTemp=DeltaTemp)
        sigma = atm.sigma(theta=theta, delta=delta)

        temp = theta * const.temp_0

        # MEC = self.OPT.MEC(mass=mass, h=h, DeltaTemp=DeltaTemp, wS=0)
        MEC = conv.kt2ms(self.OPT.getOPTParam("MEC", conv.m2ft(h), mass, DeltaTemp))

        # control parameters
        if tasDefault is None:
            tas = MEC
        else:
            tas = tasDefault

        if ROCDDefault is None:
            ROCD = conv.ft2m(1000) / 60  # [m/s]
        else:
            ROCD = ROCDDefault

        # check for speed envelope limitations
        eps = 1e-6  # float calculation precision
        maxSpeed = atm.cas2Tas(cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma)
        minSpeed = 0
        limitation = ""

        # empty envelope - keep the original calculated TAS speed
        if maxSpeed < minSpeed:
            if (tas - eps) > maxSpeed and (tas - eps) > minSpeed:
                limitation = "V"
            elif (tas + eps) < minSpeed and (tas + eps) < maxSpeed:
                limitation = "v"
            else:
                limitation = "vV"

        elif minSpeed > (tas + eps):
            if speedLimit == "applyLimit":
                tas = minSpeed
                limitation = "C"
            else:
                limitation = "v"

        elif maxSpeed < (tas - eps):
            if speedLimit == "applyLimit":
                tas = maxSpeed
                limitation = "C"
            else:
                limitation = "V"

        ESF = self.esf(flightEvolution="constTAS")
        Preq = self.Preq(sigma=sigma, tas=tas, mass=mass)

        if rating == "ARPM":
            Peng_target = self.Peng_target(
                temp=temp, DeltaTemp=DeltaTemp, ROCD=ROCD, mass=mass, Preq=Preq, ESF=ESF
            )
            Pav = self.Pav(rating="MTKF", theta=theta, delta=delta)
            Peng = min(Peng_target, Pav)

            ROCD_TEM = self.ROCD(
                Peng=Peng,
                Preq=Preq,
                mass=mass,
                ESF=ESF,
                theta=theta,
                DeltaTemp=DeltaTemp,
            )

            if ROCD_TEM < ROCD:
                ROCD = ROCD_TEM

        elif rating == "MTKF":
            Pav = self.Pav(rating="MTKF", theta=theta, delta=delta)
            Peng = Pav
            ROCD = self.ROCD(
                Peng=Peng,
                Preq=Preq,
                mass=mass,
                ESF=ESF,
                theta=theta,
                DeltaTemp=DeltaTemp,
            )

        elif rating == "MCNT":
            Pav = self.Pav(rating="MCNT", theta=theta, delta=delta)
            Peng = Pav
            ROCD = self.ROCD(
                Peng=Peng,
                Preq=Preq,
                mass=mass,
                ESF=ESF,
                theta=theta,
                DeltaTemp=DeltaTemp,
            )

        if Pav < Peng:
            limitation += "P"

        return [Pav, Peng, Preq, tas, ROCD, ESF, limitation]

    def accelerationToCruise(self):
        pass

    def cruise(self, h, mass, DeltaTemp, speedLimit=None, tasDefault=None):
        """This function computes parameters for the cruise ARPM

        :param h: altitude [m].
        :param mass: aircraft weight [kg].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param wS: longitudinal wind speed (TAS) [m s^-1].
        :param speedLimit: decision to apply or not the speed limit {applyLimit,''} [-].
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type wS: float.
        :type speedLimit: string.
        :returns: [Pav, Peng, Preq, tas, ROCD, ESF, limitation] [W, W, W, m/s, m/s, -, -]
        :rtype: float.
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        delta = atm.delta(h=h, DeltaTemp=DeltaTemp)
        sigma = atm.sigma(theta=theta, delta=delta)

        # LRC = self.OPT.LRC(mass=mass, h=h, DeltaTemp=DeltaTemp, wS=0)
        LRC = conv.kt2ms(self.OPT.getOPTParam("LRC", conv.m2ft(h), mass, DeltaTemp))

        # control parameters
        if tasDefault is None:
            if isnan(LRC):
                MEC = conv.kt2ms(
                    self.OPT.getOPTParam("MEC", conv.m2ft(h), mass, DeltaTemp)
                )
                tas = MEC
            else:
                tas = LRC
        else:
            tas = tasDefault

        ROCD = 0  # [m/s]

        # check for speed envelope limitations
        eps = 1e-6  # float calculation precision
        maxSpeed = atm.cas2Tas(cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma)
        minSpeed = 0
        limitation = ""

        # empty envelope - keep the original calculated TAS speed
        if maxSpeed < minSpeed:
            if (tas - eps) > maxSpeed and (tas - eps) > minSpeed:
                limitation = "V"
            elif (tas + eps) < minSpeed and (tas + eps) < maxSpeed:
                limitation = "v"
            else:
                limitation = "vV"

        elif minSpeed > (tas + eps):
            if speedLimit == "applyLimit":
                tas = minSpeed
                limitation = "C"
            else:
                limitation = "v"

        elif maxSpeed < (tas - eps):
            if speedLimit == "applyLimit":
                tas = maxSpeed
                limitation = "C"
            else:
                limitation = "V"

        # ESF is N/A for cruise
        ESF = 0

        Preq = self.Preq(sigma=sigma, tas=tas, mass=mass)
        Pav = self.Pav(rating="MCNT", theta=theta, delta=delta)
        Peng = min(Preq, Pav)

        if Pav < Peng:
            limitation += "P"

        return [Pav, Peng, Preq, tas, ROCD, ESF, limitation]

    def descent(
        self, h, mass, DeltaTemp, speedLimit=None, ROCDDefault=None, tasDefault=None
    ):
        """This function computes parameters for the descent ARPM

        :param h: altitude [m].
        :param mass: aircraft weight [kg].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param wS: longitudinal wind speed (TAS) [m s^-1].
        :param speedLimit: decision to apply or not the speed limit {applyLimit,''} [-].
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type wS: float.
        :type speedLimit: string.
        :returns: [Pav, Peng, Preq, tas, ROCD, ESF, limitation] [W, W, W, m/s, m/s, -, -]
        :rtype: float.
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        delta = atm.delta(h=h, DeltaTemp=DeltaTemp)
        sigma = atm.sigma(theta=theta, delta=delta)

        temp = theta * const.temp_0

        # LRC = self.OPT.LRC(mass=mass, h=h, DeltaTemp=DeltaTemp, wS=0)
        LRC = conv.kt2ms(self.OPT.getOPTParam("LRC", conv.m2ft(h), mass, DeltaTemp))

        # control parameters
        if tasDefault is None:
            if isnan(LRC):
                MEC = conv.kt2ms(
                    self.OPT.getOPTParam("MEC", conv.m2ft(h), mass, DeltaTemp)
                )
                tas = MEC
            else:
                tas = LRC
        else:
            tas = tasDefault

        if ROCDDefault is None:
            ROCD = conv.ft2m(-500) / 60  # [m/s]
        else:
            ROCD = ROCDDefault

        # check for speed envelope limitations
        eps = 1e-6  # float calculation precision
        maxSpeed = atm.cas2Tas(cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma)
        minSpeed = 0
        limitation = ""

        # empty envelope - keep the original calculated TAS speed
        if maxSpeed < minSpeed:
            if (tas - eps) > maxSpeed and (tas - eps) > minSpeed:
                limitation = "V"
            elif (tas + eps) < minSpeed and (tas + eps) < maxSpeed:
                limitation = "v"
            else:
                limitation = "vV"

        elif minSpeed > (tas + eps):
            if speedLimit == "applyLimit":
                tas = minSpeed
                limitation = "C"
            else:
                limitation = "v"

        elif maxSpeed < (tas - eps):
            if speedLimit == "applyLimit":
                tas = maxSpeed
                limitation = "C"
            else:
                limitation = "V"

        ESF = self.esf(flightEvolution="constTAS")

        Pav = self.Pav(
            rating="MTKF", theta=theta, delta=delta
        )  # verify if Pav is calualted based on MTKF rating
        Preq = self.Preq(sigma=sigma, tas=tas, mass=mass)
        Peng_target = self.Peng_target(
            temp=temp, DeltaTemp=DeltaTemp, ROCD=ROCD, mass=mass, Preq=Preq, ESF=ESF
        )
        Peng = Peng_target

        if Pav < Peng:
            limitation += "P"

        return [Pav, Peng, Preq, tas, ROCD, ESF, limitation]

    def decelerationToApproach(self):
        pass

    def approach(
        self, h, mass, DeltaTemp, speedLimit=None, ROCDDefault=None, tasDefault=None
    ):
        """This function computes parameters for the approach ARPM

        :param h: altitude [m].
        :param mass: aircraft weight [kg].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param wS: longitudinal wind speed (TAS) [m s^-1].
        :param speedLimit: decision to apply or not the speed limit {applyLimit,''} [-].
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type wS: float.
        :type speedLimit: string.
        :returns: [Pav, Peng, Preq, tas, ROCD, ESF, limitation] [W, W, W, m/s, m/s, -, -]
        :rtype: float.
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        delta = atm.delta(h=h, DeltaTemp=DeltaTemp)
        sigma = atm.sigma(theta=theta, delta=delta)

        temp = theta * const.temp_0

        # MEC = self.OPT.MEC(mass=mass, h=h, DeltaTemp=DeltaTemp, wS=0)
        MEC = conv.kt2ms(self.OPT.getOPTParam("MEC", conv.m2ft(h), mass, DeltaTemp))

        # control parameters
        if tasDefault is None:
            tas = MEC
        else:
            tas = tasDefault

        if ROCDDefault is None:
            ROCD = conv.ft2m(-300) / 60  # [m/s]
        else:
            ROCD = ROCDDefault

        # check for speed envelope limitations
        eps = 1e-6  # float calculation precision
        maxSpeed = atm.cas2Tas(cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma)
        minSpeed = 0
        limitation = ""

        # empty envelope - keep the original calculated TAS speed
        if maxSpeed < minSpeed:
            if (tas - eps) > maxSpeed and (tas - eps) > minSpeed:
                limitation = "V"
            elif (tas + eps) < minSpeed and (tas + eps) < maxSpeed:
                limitation = "v"
            else:
                limitation = "vV"

        elif minSpeed > (tas + eps):
            if speedLimit == "applyLimit":
                tas = minSpeed
                limitation = "C"
            else:
                limitation = "v"

        elif maxSpeed < (tas - eps):
            if speedLimit == "applyLimit":
                tas = maxSpeed
                limitation = "C"
            else:
                limitation = "V"

        ESF = self.esf(flightEvolution="constTAS")

        Pav = Pav = self.Pav(
            rating="MTKF", theta=theta, delta=delta
        )  # verify if Pav is calualted based on MTKF rating
        Preq = self.Preq(sigma=sigma, tas=tas, mass=mass)
        Peng_target = self.Peng_target(
            temp=temp, DeltaTemp=DeltaTemp, ROCD=ROCD, mass=mass, Preq=Preq, ESF=ESF
        )
        Peng = Peng_target

        if Pav < Peng:
            limitation += "P"

        return [Pav, Peng, Preq, tas, ROCD, ESF, limitation]

    def decelerationToFinalApproach(self):
        pass

    def finalApproach(
        self, h, mass, DeltaTemp, speedLimit=None, ROCDDefault=None, tasDefault=None
    ):
        """This function computes parameters for the final approach ARPM

        :param h: altitude [m].
        :param mass: aircraft weight [kg].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param wS: longitudinal wind speed (TAS) [m s^-1].
        :param speedLimit: decision to apply or not the speed limit {applyLimit,''} [-].
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type wS: float.
        :type speedLimit: string.
        :returns: [Pav, Peng, Preq, tas, ROCD, ESF, limitation] [W, W, W, m/s, m/s, -, -]
        :rtype: float.
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        delta = atm.delta(h=h, DeltaTemp=DeltaTemp)
        sigma = atm.sigma(theta=theta, delta=delta)

        temp = theta * const.temp_0

        # control parameters
        if tasDefault is None:
            tas = conv.kt2ms(30)
        else:
            tas = tasDefault

        if ROCDDefault is None:
            ROCD = conv.ft2m(-200) / 60  # [m/s]
        else:
            ROCD = ROCDDefault

        # check for speed envelope limitations
        eps = 1e-6  # float calculation precision
        maxSpeed = atm.cas2Tas(cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma)
        minSpeed = 0
        limitation = ""

        # empty envelope - keep the original calculated TAS speed
        if maxSpeed < minSpeed:
            if (tas - eps) > maxSpeed and (tas - eps) > minSpeed:
                limitation = "V"
            elif (tas + eps) < minSpeed and (tas + eps) < maxSpeed:
                limitation = "v"
            else:
                limitation = "vV"

        elif minSpeed > (tas + eps):
            if speedLimit == "applyLimit":
                tas = minSpeed
                limitation = "C"
            else:
                limitation = "v"

        elif maxSpeed < (tas - eps):
            if speedLimit == "applyLimit":
                tas = maxSpeed
                limitation = "C"
            else:
                limitation = "V"

        ESF = self.esf(flightEvolution="constTAS")

        Pav = Pav = self.Pav(
            rating="MTKF", theta=theta, delta=delta
        )  # verify if Pav is calualted based on MTKF rating
        Preq = self.Preq(sigma=sigma, tas=tas, mass=mass)
        Peng_target = self.Peng_target(
            temp=temp, DeltaTemp=DeltaTemp, ROCD=ROCD, mass=mass, Preq=Preq, ESF=ESF
        )
        Peng = Peng_target

        if Pav < Peng:
            limitation += "P"

        return [Pav, Peng, Preq, tas, ROCD, ESF, limitation]

    def decelerationToLanding(self):
        pass

    def landing(self, h, mass, DeltaTemp, ROCDDefault=None):
        """This function computes parameters for the landing ARPM

        :param h: altitude [m].
        :param mass: aircraft weight [kg].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param wS: longitudinal wind speed (TAS) [m s^-1].
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type wS: float.
        :returns: [Pav, Peng, Preq, tas, ROCD, ESF, limitation] [W, W, W, m/s, m/s, -, -]
        :rtype: float.
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        delta = atm.delta(h=h, DeltaTemp=DeltaTemp)
        sigma = atm.sigma(theta=theta, delta=delta)

        temp = theta * const.temp_0

        # control parameters
        if ROCDDefault is None:
            ROCD = conv.ft2m(-100) / 60  # [m/s]
        else:
            ROCD = ROCDDefault

        tas = 0

        limitation = ""

        ESF = self.esf(flightEvolution="constTAS")

        Pav = self.Pav(
            rating="MTKF", theta=theta, delta=delta
        )  # verify if Pav is calualted based on MTKF rating
        Preq = self.Preq(sigma=sigma, tas=tas, mass=mass)
        Peng_target = self.Peng_target(
            temp=temp, DeltaTemp=DeltaTemp, ROCD=ROCD, mass=mass, Preq=Preq, ESF=ESF
        )
        Peng = Peng_target

        if Pav < Peng:
            limitation += "P"

        return [Pav, Peng, Preq, tas, ROCD, ESF, limitation]

    def hover(self, h, mass, DeltaTemp):
        """This function computes parameters for the hover ARPM

        :param h: altitude [m].
        :param mass: aircraft weight [kg].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param wS: longitudinal wind speed (TAS) [m s^-1].
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type wS: float.
        :returns: [Pav, Peng, Preq, tas, ROCD, ESF, limitation] [W, W, W, m/s, m/s, -, -]
        :rtype: float.
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        delta = atm.delta(h=h, DeltaTemp=DeltaTemp)
        sigma = atm.sigma(theta=theta, delta=delta)

        # control parameters
        tas = 0
        ROCD = 0  # [m/s]

        limitation = ""

        # ESF is N/A for cruise
        ESF = 0

        Pav = self.Pav(rating="MTKF", theta=theta, delta=delta)
        Preq = self.Preq(sigma=sigma, tas=tas, mass=mass)
        Peng = Preq

        if Pav < Peng:
            limitation += "P"

        return [Pav, Peng, Preq, tas, ROCD, ESF, limitation]

    def ARPMProcedure(
        self,
        h,
        mass,
        phase,
        DeltaTemp,
        rating="ARPM",
        speedLimit=None,
        ROCDDefault=None,
        tasDefault=None,
    ):
        """This function computes parameters for the ARPM

        :param h: altitude [m].
        :param mass: aircraft weight [kg].
        :param DeltaTemp: deviation with respect to ISA [K]
        :param wS: longitudinal wind speed (TAS) [m s^-1].
        :param rating: engine rating {MTKF,MCNT,ARPM} [-].
        :param speedLimit: decision to apply or not the speed limit {applyLimit,''} [-].
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type wS: float.
        :type rating: string.
        :type speedLimit: string.
        :returns: [Pav, Peng, Preq, tas, ROCD, ESF, limitation] [W, W, W, m/s, m/s, -, -]
        :rtype: float.
        """

        if phase == "Climb":
            if h <= conv.ft2m(5):
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.takeoff(
                    h=h,
                    mass=mass,
                    DeltaTemp=DeltaTemp,
                    rating=rating,
                    speedLimit=speedLimit,
                    ROCDDefault=ROCDDefault,
                )
            elif h > conv.ft2m(5):
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.climb(
                    h=h,
                    mass=mass,
                    DeltaTemp=DeltaTemp,
                    rating=rating,
                    speedLimit=speedLimit,
                    ROCDDefault=ROCDDefault,
                    tasDefault=tasDefault,
                )

        elif phase == "Cruise":
            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.cruise(
                h=h,
                mass=mass,
                DeltaTemp=DeltaTemp,
                speedLimit=speedLimit,
                tasDefault=tasDefault,
            )

        elif phase == "Descent":
            if h >= conv.ft2m(500):
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.descent(
                    h=h,
                    mass=mass,
                    DeltaTemp=DeltaTemp,
                    speedLimit=speedLimit,
                    ROCDDefault=ROCDDefault,
                    tasDefault=tasDefault,
                )
            elif h < conv.ft2m(500) and h >= conv.ft2m(150):
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.approach(
                    h=h,
                    mass=mass,
                    DeltaTemp=DeltaTemp,
                    speedLimit=speedLimit,
                    ROCDDefault=ROCDDefault,
                    tasDefault=tasDefault,
                )
            elif h < conv.ft2m(150) and h >= conv.ft2m(5):
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.finalApproach(
                    h=h,
                    mass=mass,
                    DeltaTemp=DeltaTemp,
                    speedLimit=speedLimit,
                    ROCDDefault=ROCDDefault,
                    tasDefault=tasDefault,
                )
            elif h < conv.ft2m(5):
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.landing(
                    h=h, mass=mass, DeltaTemp=DeltaTemp, ROCDDefault=ROCDDefault
                )

        elif phase == "Hover":
            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.hover(
                h=h, mass=mass, DeltaTemp=DeltaTemp
            )

        return [Pav, Peng, Preq, tas, ROCD, ESF, limitation]


class PTD(BADAH):
    """This class implements the PTD file creator for BADAH aircraft following BADAH manual.

    :param AC: parsed aircraft.
    :type AC: badaH.Parse.
    """

    def __init__(self, AC):
        super().__init__(AC)

        self.flightEnvelope = FlightEnvelope(AC)
        self.ARPM = ARPM(AC)

    def create(self, saveToPath, DeltaTemp):
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

        # original PTD altitude list
        altitudeList = list(range(0, 500, 100))
        altitudeList.extend(range(500, 3000, 500))
        altitudeList.extend(range(3000, int(max_alt_ft), 1000))
        altitudeList.append(max_alt_ft)

        CLList_ARPM = []
        CLList_MTKF = []
        CLList_MCNT = []
        CLList = []
        DESList = []
        CRList = []
        HOVERList = []

        for mass in massList:
            CLList_ARPM.append(
                self.PTD_climb(
                    mass=mass,
                    altitudeList=altitudeList,
                    DeltaTemp=DeltaTemp,
                    rating="ARPM",
                )
            )
            CLList_MTKF.append(
                self.PTD_climb(
                    mass=mass,
                    altitudeList=altitudeList,
                    DeltaTemp=DeltaTemp,
                    rating="MTKF",
                )
            )
            CLList_MCNT.append(
                self.PTD_climb(
                    mass=mass,
                    altitudeList=altitudeList,
                    DeltaTemp=DeltaTemp,
                    rating="MCNT",
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
            HOVERList.append(
                self.PTD_hover(
                    mass=mass, altitudeList=altitudeList, DeltaTemp=DeltaTemp
                )
            )

        self.save2PTD(
            saveToPath=saveToPath,
            CLList_ARPM=CLList_ARPM,
            CLList_MTKF=CLList_MTKF,
            CLList_MCNT=CLList_MCNT,
            CRList=CRList,
            DESList=DESList,
            HOVERList=HOVERList,
            DeltaTemp=DeltaTemp,
        )

    def save2PTD(
        self,
        saveToPath,
        CLList_ARPM,
        CLList_MTKF,
        CLList_MCNT,
        CRList,
        DESList,
        HOVERList,
        DeltaTemp,
    ):
        """This function saves data to PTD file

        :param saveToPath: path to directory where PTD should be stored [-]
        :param CLList_ARPM: list of PTD data in CLIMB for BADA ARPM[-].
        :param CLList_MTKF: list of PTD data in CLIMB for BADA MTKF rating[-].
        :param CLList_MCNT: list of PTD data in CLIMB for BADA MCNT rating[-].
        :param CRList: list of PTD data in CRUISE [-].
        :param DESList: list of PTD data in DESCENT [-].
        :param HOVERList: list of PTD data in HOVER [-].
        :param DeltaTemp: deviation from ISA temperature [K]
        :type saveToPath: string.
        :type CLList_ARPM: list.
        :type CLList_MTKF: list.
        :type CLList_MCNT: list.
        :type CRList: list.
        :type DESList: list.
        :type HOVERList: list.
        :type DeltaTemp: float.
        :returns: NONE
        """

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
        file.write("Low mass CLIMB (MTKF)\n")
        file.write("=====================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass    Peng     Preq      Fuel   ESF    ROCD   gamma  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [W]      [W]     [kgm]   [-]   [fpm]   [deg]     \n"
        )

        # low mass
        list_mass = CLList_MTKF[0]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %8.0f %8.0f %7.2f %6.3f %6.0f %7.2f  %s\n"
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
                )
            )

        file.write("\n\nMedium mass CLIMB (MTKF)\n")
        file.write("========================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass    Peng     Preq      Fuel   ESF    ROCD   gamma  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [W]      [W]     [kgm]   [-]   [fpm]   [deg]     \n"
        )

        # medium mass
        list_mass = CLList_MTKF[1]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %8.0f %8.0f %7.2f %6.3f %6.0f %7.2f  %s\n"
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
                )
            )

        file.write("\n\nHigh mass CLIMB (MTKF)\n")
        file.write("======================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass    Peng     Preq      Fuel   ESF    ROCD   gamma  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [W]      [W]     [kgm]   [-]   [fpm]   [deg]     \n"
        )

        # high mass
        list_mass = CLList_MTKF[2]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %8.0f %8.0f %7.2f %6.3f %6.0f %7.2f  %s\n"
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
                )
            )

        file.write("\n\nLow mass CLIMB (MCNT)\n")
        file.write("=====================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass    Peng     Preq      Fuel   ESF    ROCD   gamma  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [W]      [W]     [kgm]   [-]   [fpm]   [deg]     \n"
        )

        # low mass
        list_mass = CLList_MCNT[0]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %8.0f %8.0f %7.2f %6.3f %6.0f %7.2f  %s\n"
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
                )
            )

        file.write("\n\nMedium mass CLIMB (MCNT)\n")
        file.write("========================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass    Peng     Preq      Fuel   ESF    ROCD   gamma  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [W]      [W]     [kgm]   [-]   [fpm]   [deg]     \n"
        )

        # medium mass
        list_mass = CLList_MCNT[1]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %8.0f %8.0f %7.2f %6.3f %6.0f %7.2f  %s\n"
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
                )
            )

        file.write("\n\nHigh mass CLIMB (MCNT)\n")
        file.write("======================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass    Peng     Preq      Fuel   ESF    ROCD   gamma  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [W]      [W]     [kgm]   [-]   [fpm]   [deg]     \n"
        )

        # high mass
        list_mass = CLList_MCNT[2]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %8.0f %8.0f %7.2f %6.3f %6.0f %7.2f  %s\n"
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
                )
            )

        file.write("\n\nLow mass CLIMB (ARPM)\n")
        file.write("=====================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass    Peng     Preq      Fuel   ESF    ROCD   gamma  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [W]      [W]     [kgm]   [-]   [fpm]   [deg]     \n"
        )

        # low mass
        list_mass = CLList_ARPM[0]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %8.0f %8.0f %7.2f %6.3f %6.0f %7.2f  %s\n"
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
                )
            )

        file.write("\n\nMedium mass CLIMB (ARPM)\n")
        file.write("========================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass    Peng     Preq      Fuel   ESF    ROCD   gamma  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [W]      [W]     [kgm]   [-]   [fpm]   [deg]     \n"
        )

        # medium mass
        list_mass = CLList_ARPM[1]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %8.0f %8.0f %7.2f %6.3f %6.0f %7.2f  %s\n"
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
                )
            )

        file.write("\n\nHigh mass CLIMB (ARPM)\n")
        file.write("======================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass    Peng     Preq      Fuel   ESF    ROCD   gamma  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [W]      [W]     [kgm]   [-]   [fpm]   [deg]     \n"
        )

        # high mass
        list_mass = CLList_ARPM[2]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %8.0f %8.0f %7.2f %6.3f %6.0f %7.2f  %s\n"
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
                )
            )

        file.write("\n\nLow mass DESCENT\n")
        file.write("================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass    Peng     Preq      Fuel   ESF    ROCD   gamma  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [W]      [W]     [kgm]   [-]   [fpm]   [deg]     \n"
        )

        # low mass
        list_mass = DESList[0]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %8.0f %8.0f %7.2f %6.3f %6.0f %7.2f  %s\n"
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
                )
            )

        file.write("\n\nMedium mass DESCENT\n")
        file.write("===================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass    Peng     Preq      Fuel   ESF    ROCD   gamma  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [W]      [W]     [kgm]   [-]   [fpm]   [deg]     \n"
        )

        # medium mass
        list_mass = DESList[1]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %8.0f %8.0f %7.2f %6.3f %6.0f %7.2f  %s\n"
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
                )
            )

        file.write("\n\nHigh mass DESCENT\n")
        file.write("=================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass    Peng     Preq      Fuel   ESF    ROCD   gamma  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [W]      [W]     [kgm]   [-]   [fpm]   [deg]     \n"
        )

        # high mass
        list_mass = DESList[2]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %8.0f %8.0f %7.2f %6.3f %6.0f %7.2f  %s\n"
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
                )
            )

        file.write("\n\nLow mass CRUISE\n")
        file.write("===============\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass    Peng     Preq      Fuel   ESF    ROCD   gamma  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [W]      [W]     [kgm]   [-]   [fpm]   [deg]     \n"
        )

        # low mass
        list_mass = CRList[0]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %8.0f %8.0f %7.2f %6.3f %6.0f %7.2f  %s\n"
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
                )
            )

        file.write("\n\nMedium mass CRUISE\n")
        file.write("==================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass    Peng     Preq      Fuel   ESF    ROCD   gamma  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [W]      [W]     [kgm]   [-]   [fpm]   [deg]     \n"
        )

        # medium mass
        list_mass = CRList[1]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %8.0f %8.0f %7.2f %6.3f %6.0f %7.2f  %s\n"
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
                )
            )

        file.write("\n\nHigh mass CRUISE\n")
        file.write("================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass    Peng     Preq      Fuel   ESF    ROCD   gamma  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [W]      [W]     [kgm]   [-]   [fpm]   [deg]     \n"
        )

        # high mass
        list_mass = CRList[2]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %8.0f %8.0f %7.2f %6.3f %6.0f %7.2f  %s\n"
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
                )
            )

        file.write("\n\nLow mass HOVER\n")
        file.write("==============\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass    Peng     Preq      Fuel   ESF    ROCD   gamma  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [W]      [W]     [kgm]   [-]   [fpm]   [deg]     \n"
        )

        # low mass
        list_mass = HOVERList[0]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %8.0f %8.0f %7.2f %6.3f %6.0f %7.2f  %s\n"
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
                )
            )

        file.write("\n\nMedium mass HOVER\n")
        file.write("=================\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass    Peng     Preq      Fuel   ESF    ROCD   gamma  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [W]      [W]     [kgm]   [-]   [fpm]   [deg]     \n"
        )

        # medium mass
        list_mass = HOVERList[1]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %8.0f %8.0f %7.2f %6.3f %6.0f %7.2f  %s\n"
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
                )
            )

        file.write("\n\nHigh mass HOVER\n")
        file.write("===============\n\n")
        file.write(
            " FL    T       p      rho     a      TAS     CAS     M     mass    Peng     Preq      Fuel   ESF    ROCD   gamma  Lim\n"
        )
        file.write(
            "[-]   [K]     [Pa]  [kg/m3] [m/s]   [kt]    [kt]    [-]    [kg]     [W]      [W]     [kgm]   [-]   [fpm]   [deg]     \n"
        )

        # high mass
        list_mass = HOVERList[2]
        for k in range(0, len(list_mass[0])):
            file.write(
                "%3d %7.2f %7.0f %6.3f %6.1f %7.2f %7.2f %6.3f %7.0f %8.0f %8.0f %7.2f %6.3f %6.0f %7.2f  %s\n"
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
                )
            )

    def PTD_climb(self, mass, altitudeList, DeltaTemp, rating):
        """This function calculates the BADAH PTD data in CLIMB

        :param mass: aircraft mass [kg]
        :param altitudeList: aircraft altitude list [ft]
        :param DeltaTemp: deviation from ISA temperature [K]
        :param rating: engine rating {MTKF,MCNT,ARPM}[-]
        :type mass: float.
        :type altitudeList: list of int.
        :type DeltaTemp: float.
        :type rating: string.
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
        Peng_complet = []
        Preq_complet = []
        ff_complet = []
        ESF_complet = []
        ROCD_complet = []
        gamma_complet = []
        Lim_complet = []

        phase = "Climb"

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            theta = atm.theta(H_m, DeltaTemp)
            delta = atm.delta(H_m, DeltaTemp)
            sigma = atm.sigma(theta=theta, delta=delta)

            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.ARPM.ARPMProcedure(
                phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=mass, rating=rating
            )

            cas = atm.tas2Cas(tas=tas, delta=delta, sigma=sigma)
            M = atm.tas2Mach(v=tas, theta=theta)
            a = atm.aSound(theta=theta)
            FL = h / 100

            CP = self.CP(Peng=Peng)
            ff = self.ff(delta=delta, CP=CP) * 60  # [kg/min]

            temp = theta * const.temp_0
            temp_const = (temp) / (temp - DeltaTemp)
            dhdt = ROCD * temp_const
            if tas == 0:
                if ROCD >= 0:
                    gamma = 90
                else:
                    gamma = -90
            else:
                gamma = conv.rad2deg(asin(dhdt / tas))

            FL_complet.append(proper_round(FL))
            T_complet.append(temp)
            p_complet.append(delta * const.p_0)
            rho_complet.append(sigma * const.rho_0)
            a_complet.append(a)
            TAS_complet.append(conv.ms2kt(tas))
            CAS_complet.append(conv.ms2kt(cas))
            M_complet.append(M)
            mass_complet.append(proper_round(mass))
            Peng_complet.append(Peng)
            Preq_complet.append(Preq)
            ff_complet.append(ff)
            ESF_complet.append(ESF)
            ROCD_complet.append(conv.m2ft(ROCD) * 60)
            gamma_complet.append(gamma)
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
            Peng_complet,
            Preq_complet,
            ff_complet,
            ESF_complet,
            ROCD_complet,
            gamma_complet,
            Lim_complet,
        ]

        return CLList

    def PTD_descent(self, mass, altitudeList, DeltaTemp):
        """This function calculates the BADAH PTD data in DESCENT

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
        Peng_complet = []
        Preq_complet = []
        ff_comlet = []
        ESF_complet = []
        ROCD_complet = []
        gamma_complet = []
        Lim_complet = []

        phase = "Descent"

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            theta = atm.theta(H_m, DeltaTemp)
            delta = atm.delta(H_m, DeltaTemp)
            sigma = atm.sigma(theta=theta, delta=delta)

            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.ARPM.ARPMProcedure(
                phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=mass
            )

            cas = atm.tas2Cas(tas=tas, delta=delta, sigma=sigma)
            M = atm.tas2Mach(v=tas, theta=theta)
            a = atm.aSound(theta=theta)
            FL = h / 100

            CP = self.CP(Peng=Peng)
            ff = self.ff(delta=delta, CP=CP) * 60  # [kg/min]

            temp = theta * const.temp_0
            temp_const = (temp) / (temp - DeltaTemp)
            dhdt = ROCD * temp_const
            if tas == 0:
                gamma = -90
            else:
                gamma = conv.rad2deg(asin(dhdt / tas))

            FL_complet.append(proper_round(FL))
            T_complet.append(temp)
            p_complet.append(delta * const.p_0)
            rho_complet.append(sigma * const.rho_0)
            a_complet.append(a)
            TAS_complet.append(conv.ms2kt(tas))
            CAS_complet.append(conv.ms2kt(cas))
            M_complet.append(M)
            mass_complet.append(proper_round(mass))
            Peng_complet.append(Peng)
            Preq_complet.append(Preq)
            ff_comlet.append(ff)
            ESF_complet.append(ESF)
            ROCD_complet.append((-1) * conv.m2ft(ROCD) * 60)
            gamma_complet.append(gamma)
            Lim_complet.append(limitation)

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
            Peng_complet,
            Preq_complet,
            ff_comlet,
            ESF_complet,
            ROCD_complet,
            gamma_complet,
            Lim_complet,
        ]

        return DESList

    def PTD_cruise(self, mass, altitudeList, DeltaTemp):
        """This function calculates the BADAH PTD data in Cruise

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
        Peng_complet = []
        Preq_complet = []
        ff_complet = []
        ESF_complet = []
        ROCD_complet = []
        gamma_complet = []
        Lim_complet = []

        phase = "Cruise"

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            theta = atm.theta(H_m, DeltaTemp)
            delta = atm.delta(H_m, DeltaTemp)
            sigma = atm.sigma(theta=theta, delta=delta)

            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.ARPM.ARPMProcedure(
                phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=mass
            )

            cas = atm.tas2Cas(tas=tas, delta=delta, sigma=sigma)
            M = atm.tas2Mach(v=tas, theta=theta)
            a = atm.aSound(theta=theta)
            FL = h / 100

            CP = self.CP(Peng=Peng)
            ff = self.ff(delta=delta, CP=CP) * 60  # [kg/min]

            temp = theta * const.temp_0
            gamma = 0

            FL_complet.append(proper_round(FL))
            T_complet.append(temp)
            p_complet.append(delta * const.p_0)
            rho_complet.append(sigma * const.rho_0)
            a_complet.append(a)
            TAS_complet.append(conv.ms2kt(tas))
            CAS_complet.append(conv.ms2kt(cas))
            M_complet.append(M)
            mass_complet.append(proper_round(mass))
            Peng_complet.append(Peng)
            Preq_complet.append(Preq)
            ff_complet.append(ff)
            ESF_complet.append(ESF)
            ROCD_complet.append(conv.m2ft(ROCD) * 60)
            gamma_complet.append(gamma)
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
            Peng_complet,
            Preq_complet,
            ff_complet,
            ESF_complet,
            ROCD_complet,
            gamma_complet,
            Lim_complet,
        ]

        return CRList

    def PTD_hover(self, mass, altitudeList, DeltaTemp):
        """This function calculates the BADAH PTD data in Cruise

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
        Peng_complet = []
        Preq_complet = []
        ff_comlet = []
        ESF_complet = []
        ROCD_complet = []
        gamma_complet = []
        Lim_complet = []

        phase = "Hover"

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            theta = atm.theta(H_m, DeltaTemp)
            delta = atm.delta(H_m, DeltaTemp)
            sigma = atm.sigma(theta=theta, delta=delta)

            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.ARPM.ARPMProcedure(
                phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=mass
            )

            cas = atm.tas2Cas(tas=tas, delta=delta, sigma=sigma)
            M = atm.tas2Mach(v=tas, theta=theta)
            a = atm.aSound(theta=theta)
            FL = h / 100

            CP = self.CP(Peng=Peng)
            ff = self.ff(delta=delta, CP=CP) * 60  # [kg/min]

            temp = theta * const.temp_0
            gamma = 0

            FL_complet.append(proper_round(FL))
            T_complet.append(temp)
            p_complet.append(delta * const.p_0)
            rho_complet.append(sigma * const.rho_0)
            a_complet.append(a)
            TAS_complet.append(conv.ms2kt(tas))
            CAS_complet.append(conv.ms2kt(cas))
            M_complet.append(M)
            mass_complet.append(proper_round(mass))
            Peng_complet.append(Peng)
            Preq_complet.append(Preq)
            ff_comlet.append(ff)
            ESF_complet.append(ESF)
            ROCD_complet.append(conv.m2ft(ROCD) * 60)
            gamma_complet.append(gamma)
            Lim_complet.append(limitation)

        HOVERList = [
            FL_complet,
            T_complet,
            p_complet,
            rho_complet,
            a_complet,
            TAS_complet,
            CAS_complet,
            M_complet,
            mass_complet,
            Peng_complet,
            Preq_complet,
            ff_comlet,
            ESF_complet,
            ROCD_complet,
            gamma_complet,
            Lim_complet,
        ]

        return HOVERList


class PTF(BADAH):
    """This class implements the PTF file creator for BADAH aircraft following BADAH manual.

    :param AC: parsed aircraft.
    :type AC: badaH.Parse.
    """

    def __init__(self, AC):
        super().__init__(AC)

        self.flightEnvelope = FlightEnvelope(AC)
        self.ARPM = ARPM(AC)

    def create(self, saveToPath, DeltaTemp):
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

        # original PTF altitude list
        altitudeList = list(range(0, 500, 100))
        altitudeList.extend(range(500, 3000, 500))
        altitudeList.extend(range(3000, int(max_alt_ft), 1000))
        altitudeList.append(max_alt_ft)

        CRList = self.PTF_cruise(
            massList=massList, altitudeList=altitudeList, DeltaTemp=DeltaTemp
        )
        CLList = self.PTF_climb(
            massList=massList,
            altitudeList=altitudeList,
            DeltaTemp=DeltaTemp,
            rating="ARPM",
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
        self, saveToPath, altitudeList, CLList, CRList, DESList, DeltaTemp, massList
    ):
        """This function saves data to PTF file

        :param saveToPath: path to directory where PTF should be stored [-]
        :param CRList: list of PTF data in CRUISE [-].
        :param CLList: list of PTF data in CLIMB [-].
        :param DESList: list of PTF data in DESCENT [-].
        :param DeltaTemp: deviation from ISA temperature [K]
        :type saveToPath: string.
        :type CRList: list.
        :type CLList: list.
        :type DESList: list.
        :type DeltaTemp: float.
        :returns: NONE
        """

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

        today = date.today()
        d3 = today.strftime("%b %d %Y")

        acICAO = self.AC.ICAO

        file = open(filename, "w")
        file.write(
            "BADA PERFORMANCE FILE                                        %s\n\n" % (d3)
        )
        file = open(filename, "a")
        file.write("AC/Type: %s\n\n" % (acICAO))
        file.write(
            " Speeds:                      Masses [kg]:             Temperature: ISA%s\n"
            % (ISA)
        )
        file.write(
            " climb   - MEC                low     -    %.0f\n"
            % (proper_round(massList[0]))
        )
        file.write(
            " cruise  - LRC                nominal -    %-4.0f        Max Alt. [ft]:%7d\n"
            % (proper_round(massList[1]), altitudeList[-1])
        )
        file.write(
            " descent - LRC                high    -    %0.f\n"
            % (proper_round(massList[2]))
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
            "    |  nom     lo   nom    hi   |  nom     lo   nom    hi    nom    |  nom     lo   nom    hi    nom  \n"
        )
        file.write(
            "======================================================================================================\n"
        )

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
        """This function calculates the BADAH PTF data in CRUISE

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

        phase = "Cruise"
        massNominal = massList[1]

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            delta = atm.delta(H_m, DeltaTemp)

            [
                Pav,
                Peng,
                Preq,
                tas_nominal,
                ROCD,
                ESF,
                limitation,
            ] = self.ARPM.ARPMProcedure(
                phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=massNominal
            )

            ff = []
            for mass in massList:
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.ARPM.ARPMProcedure(
                    phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=mass
                )

                if isnan(tas):
                    ff.append("(P)")

                else:
                    CP = self.CP(Peng=Peng)
                    ff.append(self.ff(delta=delta, CP=CP) * 60)  # [kg/min]

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

    def PTF_climb(self, massList, altitudeList, DeltaTemp, rating):
        """This function calculates the BADAH PTF data in CLIMB

        :param massList: list of aircraft mass [kg]
        :param altitudeList: aircraft altitude list [ft]
        :param DeltaTemp: deviation from ISA temperature [K]
        :param rating: engine rating {MTKF,MCNT,ARPM}[-]
        :type massList: list.
        :type altitudeList: list of int.
        :type DeltaTemp: float.
        :type rating: string.
        :returns: list of PTF CLIMB data [-]
        :rtype: list
        """

        TAS_CL_complet = []
        ROCD_CL_LO_complet = []
        ROCD_CL_NOM_complet = []
        ROCD_CL_HI_complet = []
        FF_CL_NOM_complet = []

        phase = "Climb"
        massNominal = massList[1]

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            delta = atm.delta(H_m, DeltaTemp)

            [
                Pav,
                Peng,
                Preq,
                tas_nominal,
                ROCD,
                ESF,
                limitation,
            ] = self.ARPM.ARPMProcedure(
                phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=massNominal, rating=rating
            )

            CP = self.CP(Peng=Peng)
            ff_nominal = self.ff(delta=delta, CP=CP) * 60  # [kg/min]

            ROC = []
            for mass in massList:
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.ARPM.ARPMProcedure(
                    phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=mass, rating=rating
                )

                ROC.append(conv.m2ft(ROCD) * 60)

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
        """This function calculates the BADAH PTF data in DESCENT

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

        phase = "Descent"
        massNominal = massList[1]

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            delta = atm.delta(H_m, DeltaTemp)

            [
                Pav,
                Peng,
                Preq,
                tas_nominal,
                ROCD,
                ESF,
                limitation,
            ] = self.ARPM.ARPMProcedure(
                phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=massNominal
            )

            CP = self.CP(Peng=Peng)
            ff_nominal = self.ff(delta=delta, CP=CP) * 60  # [kg/min]

            ROD = []
            ff_gamma_list = []
            for mass in massList:
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.ARPM.ARPMProcedure(
                    phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=mass
                )

                ROD.append(-conv.m2ft(ROCD) * 60)

            TAS_DES_complet.append(conv.ms2kt(tas_nominal))
            ROCD_DES_LO_complet.append(ROD[0])
            ROCD_DES_NOM_complet.append(ROD[1])
            ROCD_DES_HI_complet.append(ROD[2])
            FF_DES_NOM_complet.append(ff_nominal)

        DESList = [
            TAS_DES_complet,
            ROCD_DES_LO_complet,
            ROCD_DES_NOM_complet,
            ROCD_DES_HI_complet,
            FF_DES_NOM_complet,
        ]

        return DESList


class BadaHAircraft(BADAH):
    """This class implements the BADAH performance model following the BADAH manual.
       At the moment only TurboProp Model (TPM) is implemented

    :param filePath: path to the folder with BADAH xml formatted file.
    :param acName: ICAO aircraft designation
    :type filePath: str.
    :type acName: str
    """

    def __init__(self, badaVersion, acName, filePath=None, allData=None):
        super().__init__(self)

        self.BADAFamily = BadaFamily(BADAH=True)
        self.BADAFamilyName = "BADAH"
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
            self.MR_radius = Parser.safe_get(filtered_df, "MR_radius", None)
            self.MR_Speed = Parser.safe_get(filtered_df, "MR_Speed", None)
            self.cpr = Parser.safe_get(filtered_df, "cpr", None)
            self.n_eng = Parser.safe_get(filtered_df, "n_eng", None)
            self.P0 = Parser.safe_get(filtered_df, "P0", None)
            self.cf = Parser.safe_get(filtered_df, "cf", None)
            self.Pmax_ = Parser.safe_get(filtered_df, "Pmax_", None)
            self.cpa = Parser.safe_get(filtered_df, "cpa", None)
            self.hmo = Parser.safe_get(filtered_df, "hmo", None)
            self.vne = Parser.safe_get(filtered_df, "vne", None)
            self.MTOW = Parser.safe_get(filtered_df, "MTOW", None)
            self.OEW = Parser.safe_get(filtered_df, "OEW", None)
            self.MFL = Parser.safe_get(filtered_df, "MFL", None)
            self.MREF = Parser.safe_get(filtered_df, "MREF", None)
            self.MPL = Parser.safe_get(filtered_df, "MPL", None)
            self.VMO = Parser.safe_get(filtered_df, "VMO", None)
            self.MMO = Parser.safe_get(filtered_df, "MMO", None)

            self.flightEnvelope = FlightEnvelope(self)
            self.OPT = Optimization(self)
            self.ARPM = ARPM(self)
            self.PTD = PTD(self)
            self.PTF = PTF(self)

        # search file by file and using Synonym file
        else:
            self.ACModelAvailable = False
            self.synonymFileAvailable = False
            self.ACinSynonymFile = False

            # check if SYNONYM file exist - since for BADAH this is not a standard procedure (yet)
            synonymFile = os.path.join(
                self.filePath, "BADAH", badaVersion, "SYNONYM.xml"
            )
            if os.path.isfile(synonymFile):
                self.synonymFileAvailable = True

                # if SYNONYM exist - look for synonym based on defined acName
                self.SearchedACName = Parser.parseSynonym(
                    self.filePath, badaVersion, acName
                )

                # if cannot find - look for full name (in sub folder names) based on acName (may not be ICAO designator)
                if self.SearchedACName == None:
                    self.SearchedACName = acName
                else:
                    self.ACinSynonymFile = True

            else:
                # if doesn't exist - look for full name (in sub folder names) based on acName (may not be ICAO designator)
                self.SearchedACName = acName

            if self.SearchedACName is not None:
                acXmlFile = (
                    os.path.join(
                        self.filePath,
                        "BADAH",
                        badaVersion,
                        self.SearchedACName,
                        self.SearchedACName,
                    )
                    + ".xml"
                )
                OPTFilePath = os.path.join(self.filePath, "BADAH", badaVersion, acName)

                if os.path.isfile(acXmlFile):
                    self.ACModelAvailable = True

                    ACparsed_df = Parser.parseXML(
                        self.filePath, badaVersion, self.SearchedACName
                    )

                    self.OPTFilePath = OPTFilePath

                    self.model = Parser.safe_get(ACparsed_df, "model", None)
                    self.engineType = Parser.safe_get(ACparsed_df, "engineType", None)
                    self.engines = Parser.safe_get(ACparsed_df, "engines", None)
                    self.WTC = Parser.safe_get(ACparsed_df, "WTC", None)
                    self.ICAO = Parser.safe_get(ACparsed_df, "ICAO", None)
                    self.MR_radius = Parser.safe_get(ACparsed_df, "MR_radius", None)
                    self.MR_Speed = Parser.safe_get(ACparsed_df, "MR_Speed", None)
                    self.cpr = Parser.safe_get(ACparsed_df, "cpr", None)
                    self.n_eng = Parser.safe_get(ACparsed_df, "n_eng", None)
                    self.P0 = Parser.safe_get(ACparsed_df, "P0", None)
                    self.cf = Parser.safe_get(ACparsed_df, "cf", None)
                    self.Pmax_ = Parser.safe_get(ACparsed_df, "Pmax_", None)
                    self.cpa = Parser.safe_get(ACparsed_df, "cpa", None)
                    self.hmo = Parser.safe_get(ACparsed_df, "hmo", None)
                    self.vne = Parser.safe_get(ACparsed_df, "vne", None)
                    self.MTOW = Parser.safe_get(ACparsed_df, "MTOW", None)
                    self.OEW = Parser.safe_get(ACparsed_df, "OEW", None)
                    self.MFL = Parser.safe_get(ACparsed_df, "MFL", None)
                    self.MREF = Parser.safe_get(ACparsed_df, "MREF", None)
                    self.MPL = Parser.safe_get(ACparsed_df, "MPL", None)
                    self.VMO = Parser.safe_get(ACparsed_df, "VMO", None)
                    self.MMO = Parser.safe_get(ACparsed_df, "MMO", None)

                    self.flightEnvelope = FlightEnvelope(self)
                    self.OPT = Optimization(self)
                    self.ARPM = ARPM(self)
                    self.PTD = PTD(self)
                    self.PTF = PTF(self)

                else:
                    # AC name cannot be found
                    raise ValueError(acName + " Cannot be found")

    def __str__(self):
        return f"(BADAH, AC_name: {self.acName}, searched_AC_name: {self.SearchedACName}, model_ICAO: {self.ICAO}, ID: {id(self.AC)})"
