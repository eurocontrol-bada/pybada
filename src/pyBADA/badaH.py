"""
Generic BADAH aircraft performance module
"""

import xml.etree.ElementTree as ET
from datetime import date
import os
import numpy as np
from math import sqrt, pow, pi, cos, asin, radians, isnan
import pandas as pd


from pyBADA import constants as const
from pyBADA import conversions as conv
from pyBADA import atmosphere as atm
from pyBADA import configuration as configuration
from pyBADA.aircraft import Helicopter, BadaFamily, Bada


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
    """This class implements the BADAH parsing mechanism to parse xml BADAH
    files."""

    def __init__(self):
        pass

    @staticmethod
    def parseXML(filePath, acName):
        """Parses the BADAH XML file for a specific aircraft model and
        extracts various parameters.

        This function parses the BADAH aircraft XML file for a given aircraft
        model (acName). It retrieves general information about the aircraft,
        engine type, aerodynamic configurations, performance parameters, and
        more.

        :param filePath: The path to the folder containing the BADAH XML file.
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
    def readSynonym(filePath):
        """Parses the BADAH Synonym XML file and returns a dictionary mapping
        aircraft code names to their respective model files.

        :param filePath: Path to the directory containing the BADA4 synonym
            XML file.
        :type filePath: str
        :returns: A dictionary where the keys are aircraft codes and the
            values are associated file names.
        :rtype: dict
        :raises IOError: If the XML file is missing or has an invalid format.
            This function attempts to read the synonym XML file, parse its
            contents, and store the mappings in a dictionary. The file
            contains aircraft code, manufacturer, ICAO designation, and file
            name data for each aircraft in the synonym list.
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
        """Retrieves the file name associated with a given aircraft code from
        the BADAH synonym file.

        :param filePath: Path to the directory containing the BADAH synonym XML file.
        :param acName: The ICAO aircraft code or name to search for in the synonym file.
        :type filePath: str
        :type acName: str
        :returns: The associated file name if found, otherwise None.
        :rtype: str

        This function uses the `readSynonym` function to load the synonym dictionary and looks up the
        given aircraft code (acName) to return the associated file name. If no match is found, it returns None.
        """
        synonym_fileName = Parser.readSynonym(filePath)

        if acName in synonym_fileName:
            fileName = synonym_fileName[acName]
            return fileName
        else:
            return None

    @staticmethod
    def parseAll(badaVersion, filePath=None):
        """Parses all BADAH XML-formatted files and compiles the data into a
        single DataFrame.

        This function reads the BADAH aircraft performance model data by
        parsing the XML files for each aircraft model found in the specified
        directory. If the synonym XML file is present, it maps synonyms
        (alternative names for aircraft) to their respective model files and
        includes them in the output DataFrame.

        :param badaVersion: The version of BADAH being used (e.g., '1.1').
        :param filePath: Optional path to the directory containing the BADAH
            files. If not provided, it uses the default path.
        :type badaVersion: str
        :type filePath: str, optional
        :returns: A pandas DataFrame containing all parsed BADAH model data,
            including any synonyms found.
        :rtype: pd.DataFrame
        :raises IOError: If an error occurs while reading or parsing the XML
            files. This function first checks if a synonym XML file exists to
            map synonyms to model files. Then, it parses all XML files in the
            directory and its subfolders, merges the parsed data into a final
            DataFrame, and returns it.
        """

        if filePath is None:
            filePath = configuration.getBadaVersionPath(
                badaFamily="BADAH", badaVersion=badaVersion
            )
        else:
            filePath = filePath

        synonym_fileName = Parser.readSynonym(filePath)

        # get names of all the folders in the main BADA model folder to search for XML files
        subfolders = configuration.list_subfolders(filePath)

        merged_df = pd.DataFrame()

        if synonym_fileName:
            for synonym in synonym_fileName:
                file = synonym_fileName[synonym]

                if file in subfolders:
                    # parse the original XML of a model
                    df = Parser.parseXML(filePath, file)

                    # rename acName in the data frame to match the synonym model name
                    df.at[0, "acName"] = synonym

                    # Merge DataFrames
                    merged_df = pd.concat([merged_df, df], ignore_index=True)

        else:
            for file in subfolders:
                # Parse the original XML of a model
                df = Parser.parseXML(filePath, file)

                # Merge DataFrames
                merged_df = pd.concat([merged_df, df], ignore_index=True)

        return merged_df


class BADAH(Helicopter, Bada):
    """This class implements the part of BADAH performance model that will be
    used in other classes following the BADAH manual.

    :param AC: Aircraft object {BADAH}.
    :type AC: badaHAircraft.
    """

    def __init__(self, AC):
        super().__init__()
        self.AC = AC

    def mu(self, tas):
        """Computes the advance ratio (mu) for the aircraft based on true
        airspeed (TAS) and rotor speed.

        The advance ratio is a non-dimensional parameter that relates the
        forward speed of the aircraft to the rotational speed of its main
        rotor.

        :param tas: True airspeed (TAS) in meters per second [m/s].
        :type tas: float
        :return: Advance ratio (mu) [-].
        :rtype: float
        """

        # mu = (tas * math.cos(gamma))/(self.AC.MR_Speed*self.AC.MR_radius) #TODO: apply gamma modification
        mu = tas / (self.AC.MR_Speed * self.AC.MR_radius)

        return mu

    def CT(self, mass, rho, phi):
        """Computes the thrust coefficient (CT) for the aircraft.

        The thrust coefficient is a dimensionless quantity that represents the
        thrust produced by the aircraft's rotor in relation to the air
        density, rotor radius, and rotor speed.

        :param mass: Aircraft mass in kilograms [kg].
        :param rho: Air density in kilograms per cubic meter [kg/m³].
        :param phi: Bank angle in degrees [deg].
        :type mass: float
        :type rho: float
        :type phi: float
        :return: Thrust coefficient (CT) [-].
        :rtype: float
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
        """Computes the power required coefficient (CPreq) based on the
        advance ratio (mu) and thrust coefficient (CT).

        The power required coefficient relates to the total power required to
        maintain flight, factoring in the aerodynamic performance of the rotor
        in different operating regimes.

        :param mu: Advance ratio [-].
        :param CT: Thrust coefficient [-].
        :type mu: float
        :type CT: float
        :return: Power required coefficient (CPreq) [-].
        :rtype: float
        """

        CPreq = (
            self.AC.cpr[0]
            + self.AC.cpr[1] * pow(mu, 2)
            + self.AC.cpr[2]
            * CT
            * sqrt(sqrt(pow(mu, 4) + pow(CT, 2)) - pow(mu, 2))
            + self.AC.cpr[3] * pow(mu, 3)
            + self.AC.cpr[4] * pow(CT, 2) * pow(mu, 3)
        )

        return CPreq

    def Preq(self, sigma, tas, mass, phi=0.0):
        """Computes the power required for the aircraft to maintain flight
        based on various factors such as air density, true airspeed (TAS),
        aircraft mass, and bank angle.

        :param sigma: Normalized air density [-], which is the ratio of the
            current air density to sea level air density.
        :param tas: True airspeed (TAS) in meters per second [m/s].
        :param mass: Aircraft mass in kilograms [kg].
        :param phi: Bank angle in degrees [deg], default is 0 for straight
            flight.
        :type sigma: float
        :type tas: float
        :type mass: float
        :type phi: float
        :returns: Power required for the aircraft in watts [W].
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
        """Computes the target engine power required to achieve a specific
        rate of climb or descent.

        :param ROCD: Rate of climb or descent in meters per second [m/s].
        :param mass: Aircraft mass in kilograms [kg].
        :param Preq: Power required in watts [W].
        :param ESF: Energy share factor, a dimensionless factor [-].
        :param temp: Atmospheric temperature in kelvins [K].
        :param DeltaTemp: Deviation from the International Standard Atmosphere
            (ISA) temperature in kelvins [K].
        :type ROCD: float
        :type mass: float
        :type Preq: float
        :type ESF: float
        :type temp: float
        :type DeltaTemp: float
        :returns: Target engine power in watts [W].
        :rtype: float
        """

        temp_const = temp / (temp - DeltaTemp)
        Peng_target = (ROCD / ESF) * mass * const.g * temp_const + Preq

        return Peng_target

    def CPav(self, rating, delta, theta):
        """Computes the power available coefficient (CPav) based on engine
        type, throttle rating, normalized air pressure (delta), and normalized
        temperature (theta).

        :param rating: Engine throttle setting, e.g., {MTKF, MCNT}.
        :param delta: Normalized air pressure [-].
        :param theta: Normalized air temperature [-].
        :type rating: str
        :type delta: float
        :type theta: float
        :return: Power available coefficient [-].
        :rtype: float
        :raises ValueError: If the engine rating or type is unknown.
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
        """Computes the maximum power available for all engines at a given
        throttle setting.

        :param rating: Throttle setting, e.g., {MTKF, MCNT}.
        :type rating: str
        :return: Maximum all-engine power in watts [W].
        :rtype: float
        :raises ValueError: If the specified throttle setting is not
            recognized.
        """

        if rating not in self.AC.Pmax_.keys():
            raise ValueError("Unknown engine rating " + rating)
        return self.AC.Pmax_[rating]

    def Pav(self, rating, delta, theta):
        """Computes the power available at the given throttle setting, based
        on normalized pressure and temperature.

        :param rating: Throttle setting, e.g., {MTKF, MCNT}.
        :param delta: Normalized pressure [-], ratio of actual pressure to
            standard sea level pressure.
        :param theta: Normalized temperature [-], ratio of actual temperature
            to standard sea level temperature.
        :type rating: str
        :type delta: float
        :type theta: float
        :return: Available power in watts [W].
        :rtype: float
        :raises ValueError: If the specified throttle setting is not
            recognized.
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
        """Computes the torque value as a percentage of the reference torque
        (P0).

        :param Peng: All-engine power in watts [W].
        :type Peng: float
        :return: Torque value as a percentage [%] of the reference torque.
        :rtype: float
        """

        Q = Peng / self.AC.P0

        return Q

    def CP(self, Peng):
        """Computes the engine power coefficient (CP) based on the given all-
        engine power.

        :param Peng: All-engine power in watts [W].
        :type Peng: float
        :return: Engine power coefficient [-].
        :rtype: float
        """

        CP = Peng / (
            const.rho_0
            * pi
            * pow(self.AC.MR_radius, 2)
            * pow(self.AC.MR_Speed * self.AC.MR_radius, 3)
        )

        return CP

    def ff(self, delta, CP):
        """Computes the fuel flow rate based on normalized pressure and power
        coefficient.

        :param delta: Normalized pressure [-], which is the ratio of actual
            air pressure to standard sea-level pressure.
        :param CP: Power coefficient [-], representing the power output in
            relation to the engine's maximum power.
        :type delta: float
        :type CP: float
        :return: Fuel flow rate in kilograms per second [kg/s].
        :rtype: float
        :raises ValueError: If the engine type is unknown.
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
        """Computes the Rate of Climb or Descent (ROCD) for an aircraft.

        :param Peng: All-engine power available [W].
        :param Preq: Power required for steady flight [W].
        :param mass: Aircraft's current mass [kg].
        :param ESF: Energy share factor [-], a multiplier used to adjust power
            distribution in different flight phases.
        :param theta: Normalized temperature [-], ratio of actual temperature
            to standard sea-level temperature.
        :param DeltaTemp: Deviation from the International Standard Atmosphere
            (ISA) temperature [K].
        :type Peng: float
        :type Preq: float
        :type mass: float
        :type ESF: float
        :type theta: float
        :type DeltaTemp: float
        :return: Rate of Climb or Descent (ROCD) in meters per second [m/s].
        :rtype: float
        """

        temp = theta * const.temp_0
        ROCD = (
            ((temp - DeltaTemp) / temp)
            * (Peng - Preq)
            * ESF
            / (mass * const.g)
        )

        return ROCD


class FlightEnvelope(BADAH):
    """This class is a BADAH aircraft subclass and implements the flight
    envelope caclulations following the BADAH manual.

    :param AC: Aircraft object {BADAH}.
    :type AC: badaHAircraft.
    """

    def __init__(self, AC):
        super().__init__(AC)

    def maxAltitude(self):
        """Computes the maximum operational altitude for the aircraft.

        :return: Maximum altitude in meters [m].
        :rtype: float.
        """

        hMax = conv.ft2m(self.AC.hmo)
        return hMax

    def VMax(self):
        """Computes the maximum speed in Calibrated Airspeed (CAS) as limited
        by the flight envelope.

        :return: Maximum CAS speed in meters per second [m/s].
        :rtype: float.
        """

        Vmax = conv.kt2ms(self.AC.vne)
        return Vmax

    def speedEnvelope_powerLimited(
        self, h, mass, DeltaTemp, rating="MCNT", rateOfTurn=0
    ):
        """Computes the maximum and minimum speeds (CAS) within the certified
        flight envelope, taking into account engine thrust limitations.

        :param h: Altitude in meters [m].
        :param mass: Aircraft mass in kilograms [kg].
        :param DeltaTemp: Deviation from the International Standard Atmosphere
            (ISA) temperature [K].
        :param rating: Engine rating (e.g., "MTKF", "MCNT") determining the
            power output [-].
        :param rateOfTurn: Rate of turn in degrees per second, which affects
            bank angle [°/s].
        :type h: float
        :type mass: float
        :type DeltaTemp: float
        :type rating: str
        :type rateOfTurn: float
        :return: A tuple containing the minimum and maximum thrust- limited
            CAS speeds [m/s].
        :rtype: tuple(float, float)
        """

        [theta, delta, sigma] = atm.atmosphereProperties(
            h=h, DeltaTemp=DeltaTemp
        )
        Pmax = self.Pav(rating=rating, theta=theta, delta=delta)
        Pmin = 0.1 * self.AC.P0  # No minimum power model: assume 10% torque

        VminCertified = 0
        VmaxCertified = self.VMax()

        CASlist = []
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

    def Vx(self, h, mass, DeltaTemp, rating="MTKF", rateOfTurn=0):
        """Computes the best angle climb speed (TAS) by finding the speed that
        maximizes the excess power per unit speed within the helicopter's
        performance envelope.

        :param h: Altitude in meters [m].
        :param mass: Aircraft mass in kilograms [kg].
        :param DeltaTemp: Deviation from the International Standard Atmosphere
            (ISA) temperature [K].
        :param rating: Engine rating (e.g., "MTKF", "MCNT") determining the
            power output [-].
        :param rateOfTurn: Rate of turn in degrees per second, which affects
            the bank angle [°/s].
        :type h: float
        :type mass: float
        :type DeltaTemp: float
        :type rating: str
        :type rateOfTurn: float
        :return: The true airspeed (TAS) corresponding to the best angle climb
            speed [m/s].
        :rtype: float
        """
        [theta, delta, sigma] = atm.atmosphereProperties(
            h=h, DeltaTemp=DeltaTemp
        )

        VminCertified = 0 + 5  # putting some margin to not start at 0 speed
        VmaxCertified = self.VMax()

        excessPowerList = []
        VxList = []

        for CAS in np.linspace(
            VminCertified, VmaxCertified, num=200, endpoint=True
        ):
            TAS = atm.cas2Tas(cas=CAS, delta=delta, sigma=sigma)
            bankAngle = self.bankAngle(rateOfTurn=rateOfTurn, v=TAS)
            Preq = self.Preq(sigma=sigma, tas=TAS, mass=mass, phi=bankAngle)
            Pav = self.Pav(rating=rating, theta=theta, delta=delta)

            tempConst = (theta * const.temp_0 - DeltaTemp) / (
                theta * const.temp_0
            )

            excessPowerList.append(
                (Pav - Preq) * tempConst / TAS
            )  # including speed and impact of the temperature deviation from ISA conditions
            VxList.append(CAS)

        VxCAS = VxList[excessPowerList.index(max(excessPowerList))]
        [VxM, VxCAS, VxTAS] = atm.convertSpeed(
            v=conv.ms2kt(VxCAS),
            speedType="CAS",
            theta=theta,
            delta=delta,
            sigma=sigma,
        )

        return VxTAS


class Optimization(BADAH):
    """This class implements the BADAH optimization following the BADAH
    manual.

    :param AC: Aircraft object {BADAH}.
    :type AC: badaHAircraft.
    """

    def __init__(self, AC):
        super().__init__(AC)

        self.flightEnvelope = FlightEnvelope(AC)

    def MRC(self, h, mass, DeltaTemp, wS):
        """Computes the True Airspeed (TAS) representing Maximum Range Cruise
        (MRC) for given flight conditions.

        The Maximum Range Cruise speed is the speed that maximizes the
        aircraft's range per unit of fuel, which is determined by balancing
        the fuel flow rate and airspeed. The algorithm ensures that the
        computed TAS stays within the power available limitations of the
        aircraft.

        :param h: Altitude in meters [m].
        :param mass: Aircraft mass in kilograms [kg].
        :param DeltaTemp: Deviation from International Standard Atmosphere
            (ISA) temperature in Kelvin [K].
        :param wS: Longitudinal wind speed (TAS) in meters per second [m/s].
        :type h: float
        :type mass: float
        :type DeltaTemp: float
        :type wS: float
        :return: Maximum Range Cruise (MRC) speed in True Airspeed (TAS)
            [m/s].
        :rtype: float.
        :raises ValueError: If no valid MRC speed is found, the function will
            return NaN.
        """

        # NOTE: check for precision of algorithm needed. Possible local minima, instead of global minima

        [theta, delta, sigma] = atm.atmosphereProperties(
            h=h, DeltaTemp=DeltaTemp
        )  # atmosphere properties

        # max TAS speed limitation
        Vmax = atm.cas2Tas(
            cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma
        )

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
        """Computes the True Airspeed (TAS) representing Long Range Cruise
        (LRC) for the given flight conditions.

        The Long Range Cruise speed is the speed that allows for 99% of the
        specific range (range per unit of fuel) of the Maximum Range Cruise
        (MRC) speed while offering a higher cruise speed. This function
        ensures that the computed TAS remains within the aircraft's power
        limitations.

        :param h: Altitude in meters [m].
        :param mass: Aircraft mass in kilograms [kg].
        :param DeltaTemp: Deviation from International Standard Atmosphere
            (ISA) temperature in Kelvin [K].
        :param wS: Longitudinal wind speed (TAS) in meters per second [m/s].
        :type h: float
        :type mass: float
        :type DeltaTemp: float
        :type wS: float
        :return: Long Range Cruise (LRC) speed in True Airspeed (TAS) [m/s].
        :rtype: float.
        :raises ValueError: If no valid LRC speed is found, the function will
            return NaN. The algorithm starts by computing the MRC speed. Using
            the MRC as a reference, it then calculates the LRC by finding the
            TAS that achieves 99% of the specific range of the MRC while
            staying within the aircraft’s thrust limitations.
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
        Vmax = atm.cas2Tas(
            cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma
        )
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
        """Computes the True Airspeed (TAS) representing Maximum Endurance
        Cruise (MEC) for the given flight conditions.

        The Maximum Endurance Cruise speed is the speed that maximizes the
        time an aircraft can stay in the air for a given amount of fuel,
        making it ideal for loiter operations. This function minimizes fuel
        flow (ff) to determine the most fuel-efficient speed.

        :param h: Altitude in meters [m].
        :param mass: Aircraft weight in kilograms [kg].
        :param DeltaTemp: Deviation from the International Standard Atmosphere
            (ISA) temperature in Kelvin [K].
        :param wS: Longitudinal wind speed (TAS) in meters per second [m/s].
        :type h: float
        :type mass: float
        :type DeltaTemp: float
        :type wS: float
        :return: Maximum Endurance Cruise (MEC) speed in True Airspeed (TAS)
            [m/s].
        :rtype: float
        :raises: If no valid MEC speed is found, the function returns NaN. The
            algorithm iterates over possible True Airspeeds (TAS) and computes
            the fuel flow for each, aiming to minimize fuel consumption and
            return the TAS that achieves this.
        """

        [theta, delta, sigma] = atm.atmosphereProperties(
            h=h, DeltaTemp=DeltaTemp
        )  # atmosphere properties

        # max TAS speed limitation
        Vmax = atm.cas2Tas(
            cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma
        )

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
        """Parses BADAH OPT ASCII formatted files and stores data for each
        available delta temperature in the file.

        :param filename: Path to the ___.OPT ASCII formatted file.
        :type filename: str
        :return: Dictionary of delta temperature values and corresponding data
            from the OPT file.
        :rtype: dict This function reads and processes a BADAH OPT file,
            extracting delta temperature values and the corresponding
            performance data. The data is stored in a dictionary where each
            delta temperature is mapped to its respective dataset of
            performance values.
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
                        filter(
                            None,
                            lines[startIdx].split("|")[1].strip().split(" "),
                        )
                    )
                ]

                for j in range(
                    startIdx + 3, startIdx + 3 + self.tableDimensionRows, 1
                ):
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
        """Finds the nearest index or indices for a given value in a sorted
        array.

        If the value is lower or higher than the array’s bounds, a single
        index is returned. If the value lies between two elements, two closest
        indices (left and right) are returned.

        :param value: The value to find the nearest match for.
        :param array: The sorted array of values.
        :type value: float
        :type array: list[float]
        :return: A list of nearest index or indices.
        :rtype: list[float] The function uses binary search to efficiently
            find the nearest value or values, ensuring precise interpolation
            when needed.
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
        """Calculates the interpolated value of an OPT parameter based on two
        optimizing factors.

        If the exact values of the factors exist in the data, the function
        returns the corresponding OPT value. Otherwise, it interpolates
        between the nearest two values to provide a more accurate result.

        :param var_1: The first optimizing factor.
        :param var_2: The second optimizing factor.
        :param detaTauList: List of values belonging to the specified delta
            temperature from the OPT file.
        :type var_1: float
        :type var_2: float
        :type detaTauList: list[float]
        :return: Interpolated or exact OPT value based on the input factors.
        :rtype: float This function handles both single-index and two- index
            cases for the nearest values, ensuring correct interpolation in
            the case of multiple values being found.
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
                nearestIdx_1[0] * (self.tableDimensionColumns)
                + nearestIdx_2[0]
            ]
            varTemp_2 = var_3_list[
                nearestIdx_1[0] * (self.tableDimensionColumns)
                + nearestIdx_2[1]
            ]

            varTemp_3 = var_3_list[
                nearestIdx_1[1] * (self.tableDimensionColumns)
                + nearestIdx_2[0]
            ]
            varTemp_4 = var_3_list[
                nearestIdx_1[1] * (self.tableDimensionColumns)
                + nearestIdx_2[1]
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
        """Retrieves the value of the specified optimization parameter (e.g.,
        LRC, MEC, MRC) from the BADA OPT file, either directly or through
        interpolation based on the given flight conditions.

        The function searches for the requested optimization parameter value using two optimizing factors.
        If the exact DeltaTemp exists in the OPT file, it retrieves the value. Otherwise, the function interpolates
        between the closest available DeltaTemp values.

        :param optParam: Name of the optimization parameter file to query. Possible values include {LRC, MEC, MRC}.
        :param var_1: First optimizing factor (e.g., speed, altitude) used to retrieve the value from the OPT file.
        :param var_2: Second optimizing factor used to retrieve the value from the OPT file.
        :param DeltaTemp: Deviation from the International Standard Atmosphere (ISA) temperature [K], used to retrieve or interpolate the value.
        :type optParam: str
        :type var_1: float
        :type var_2: float
        :type DeltaTemp: float
        :return: The optimization parameter value, either directly from the file or interpolated.
        :rtype: float

        .. note::
           The function assumes that the arrays in the OPT file are sorted (as per the design of BADA OPT files).
           If the exact DeltaTemp value is not present in the file, the function interpolates between the nearest
           DeltaTemp values within the range of [-20, 20] degrees.
        """

        filename = os.path.join(
            self.AC.filePath,
            self.AC.acName,
            optParam + ".OPT",
        )

        detaTauDict = self.parseOPT(filename=filename)

        if DeltaTemp in detaTauDict:
            # value of DeltaTemp exist in the OPT file
            optVal = self.calculateOPTparam(
                var_1, var_2, detaTauDict[DeltaTemp]
            )
        else:
            # value of DeltaTemp does not exist in OPT file - will be interpolated. But only within the range of <-20;20>
            nearestIdx = np.array(
                self.findNearestIdx(DeltaTemp, list(detaTauDict))
            )

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
                    DeltaTemp,
                    [DeltaTemp_new_1, DeltaTemp_new_2],
                    [optVal_1, optVal_2],
                )

        return optVal


class ARPM(BADAH):
    """This class is a BADAH aircraft subclass and implements the Airline
    Procedure Model (ARPM) following the BADAH user manual.

    :param AC: Aircraft object {BADAH}.
    :type AC: badaHAircraft.
    """

    def __init__(self, AC):
        super().__init__(AC)

        self.flightEnvelope = FlightEnvelope(AC)
        self.OPT = Optimization(AC)

    def takeoff(
        self,
        h,
        mass,
        DeltaTemp,
        rating="ARPM",
        speedLimit=None,
        ROCDDefault=None,
    ):
        """Computes various parameters for the aircraft takeoff phase using
        the ARPM model (or other specified engine ratings).

        This function calculates key takeoff parameters, including the available and required power, true airspeed, rate of climb (ROCD), and other
        performance metrics. It also checks for speed limitations based on the flight envelope and applies them as necessary.

        :param h: Altitude above sea level [m].
        :param mass: Aircraft weight [kg].
        :param DeltaTemp: Deviation from the International Standard Atmosphere (ISA) temperature [K].
        :param rating: Engine rating mode, defaults to "ARPM". Other options include {MTKF, MCNT}.
        :param speedLimit: Optional parameter to specify if speed limits should be applied {"applyLimit", None}.
        :param ROCDDefault: Default rate of climb or descent [m/s], optional.
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type rating: str.
        :type speedLimit: str.
        :type ROCDDefault: float, optional.

        :returns: A list of computed values for:
                  - Pav: Available power [W].
                  - Peng: Engine power [W].
                  - Preq: Required power [W].
                  - tas: True airspeed [m/s].
                  - ROCD: Rate of Climb or Descent [m/s].
                  - ESF: Energy share factor [-].
                  - limitation: Speed and power limitations encountered during takeoff.
        :rtype: list[float]

        The function calculates these values by:
        - Determining atmosphere conditions at the given altitude (using temperature, pressure, and density).
        - Computing power requirements and available engine power based on the engine rating (e.g., ARPM, MTKF, MCNT).
        - Applying optional speed limitations from the flight envelope and checking if engine power or speed limits constrain performance.

        If the engine rating is "ARPM", the function tries to reach a target power for the rate of climb (ROCD), and adjusts it based on the available power (Pav).
        If the rating is "MTKF" or "MCNT", it simply uses the maximum available power for that setting.

        .. note::
           The function automatically handles speed envelope limitations by applying adjustments to the true airspeed (tas) if required.

        .. warning::
           The accuracy of the output depends on the precision of the flight envelope model and other internal aircraft parameters.
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
        maxSpeed = atm.cas2Tas(
            cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma
        )
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
                temp=temp,
                DeltaTemp=DeltaTemp,
                ROCD=ROCD,
                mass=mass,
                Preq=Preq,
                ESF=ESF,
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
        """Computes various parameters for the aircraft climb phase using the
        ARPM model or other engine ratings.

        This function calculates key climb parameters, including available and required power, true airspeed (TAS),
        rate of climb (ROCD), and performance limitations. It takes into account speed envelope constraints
        and engine power limits based on the flight altitude and aircraft mass.

        :param h: Altitude above sea level [m].
        :param mass: Aircraft weight [kg].
        :param DeltaTemp: Deviation from the International Standard Atmosphere (ISA) temperature [K].
        :param rating: Engine rating mode, defaults to "ARPM". Other options include {MTKF, MCNT}.
        :param speedLimit: Optional parameter to apply speed limits. Use {"applyLimit", None}.
        :param ROCDDefault: Default rate of climb or descent [m/s], optional.
        :param tasDefault: Default true airspeed (TAS) [m/s], optional.
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type rating: str.
        :type speedLimit: str.
        :type ROCDDefault: float, optional.
        :type tasDefault: float, optional.

        :returns: A list of computed values:
                  - Pav: Available power [W].
                  - Peng: Engine power [W].
                  - Preq: Required power [W].
                  - tas: True airspeed [m/s].
                  - ROCD: Rate of climb [m/s].
                  - ESF: Energy share factor [-].
                  - limitation: Performance limitations encountered during the climb (e.g., speed, power limits).
        :rtype: list[float]

        The function calculates these values by:
        - Determining atmospheric conditions at the given altitude (using temperature, pressure, and density).
        - Calculating the maximum endurance cruise speed (MEC) for the given altitude and mass, and setting TAS accordingly.
        - Checking speed envelope limitations (minimum and maximum allowable speeds) and applying them if necessary.
        - Computing the required power (Preq) and available power (Pav) based on the engine rating (e.g., ARPM, MTKF, MCNT).
        - Adjusting the rate of climb (ROCD) if engine power is insufficient to reach the desired target climb rate.

        If the engine rating is "ARPM", the function attempts to reach the target climb rate by adjusting the power (Peng).
        If the rating is "MTKF" or "MCNT", it uses the maximum available power for that setting.

        .. note::
           The function handles speed envelope limitations automatically by applying speed adjustments if necessary.

        .. warning::
           The output accuracy depends on the precision of the atmospheric model and the flight envelope data.
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        delta = atm.delta(h=h, DeltaTemp=DeltaTemp)
        sigma = atm.sigma(theta=theta, delta=delta)

        temp = theta * const.temp_0

        # MEC = self.OPT.MEC(mass=mass, h=h, DeltaTemp=DeltaTemp, wS=0)
        MEC = conv.kt2ms(
            self.OPT.getOPTParam("MEC", conv.m2ft(h), mass, DeltaTemp)
        )

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
        maxSpeed = atm.cas2Tas(
            cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma
        )
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
                temp=temp,
                DeltaTemp=DeltaTemp,
                ROCD=ROCD,
                mass=mass,
                Preq=Preq,
                ESF=ESF,
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
        """Computes various parameters for the aircraft cruise phase using the
        ARPM model or default speed.

        This function calculates key cruise parameters, including available and required power, true airspeed (TAS),
        and potential limitations due to the flight envelope or engine power. The calculations take into account
        atmospheric conditions, altitude, and aircraft mass.

        :param h: Altitude above sea level [m].
        :param mass: Aircraft weight [kg].
        :param DeltaTemp: Deviation from the International Standard Atmosphere (ISA) temperature [K].
        :param speedLimit: Optional parameter to apply speed limits. Use {"applyLimit", None}.
        :param tasDefault: Optional true airspeed (TAS) [m/s].
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type speedLimit: str, optional.
        :type tasDefault: float, optional.

        :returns: A list of computed values:
                  - Pav: Available power [W].
                  - Peng: Engine power [W].
                  - Preq: Required power [W].
                  - tas: True airspeed [m/s].
                  - ROCD: Rate of climb or descent, set to 0 for cruise [m/s].
                  - ESF: Energy share factor [-], set to 0 for cruise.
                  - limitation: Any performance limitations encountered during the cruise (e.g., speed, power limits).
        :rtype: list[float]

        The function determines the Long Range Cruise (LRC) speed or a default speed if provided. It checks for any
        speed envelope limitations and calculates power requirements for the given conditions. If the available power
        is less than the required power, performance limitations are recorded.

        .. note::
           ESF (Energy Share Factor) is not applicable in cruise mode and is therefore set to 0.
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        delta = atm.delta(h=h, DeltaTemp=DeltaTemp)
        sigma = atm.sigma(theta=theta, delta=delta)

        # LRC = self.OPT.LRC(mass=mass, h=h, DeltaTemp=DeltaTemp, wS=0)
        LRC = conv.kt2ms(
            self.OPT.getOPTParam("LRC", conv.m2ft(h), mass, DeltaTemp)
        )

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
        maxSpeed = atm.cas2Tas(
            cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma
        )
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
        self,
        h,
        mass,
        DeltaTemp,
        speedLimit=None,
        ROCDDefault=None,
        tasDefault=None,
    ):
        """Computes various parameters for the aircraft descent phase using
        the ARPM model or default speed.

        This function calculates key descent parameters, including available and required power, true airspeed (TAS),
        rate of descent (ROD), and potential performance limitations. The calculations take into account atmospheric
        conditions, altitude, and aircraft mass.

        :param h: Altitude above sea level [m].
        :param mass: Aircraft weight [kg].
        :param DeltaTemp: Deviation from the International Standard Atmosphere (ISA) temperature [K].
        :param speedLimit: Optional parameter to apply speed limits. Use {"applyLimit", None}.
        :param ROCDDefault: Default rate of climb or descent [m/s], optional.
        :param tasDefault: Optional true airspeed (TAS) [m/s].
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type speedLimit: str, optional.
        :type ROCDDefault: float, optional.
        :type tasDefault: float, optional.

        :returns: A list of computed values:
                  - Pav: Available power [W].
                  - Peng: Engine power [W].
                  - Preq: Required power [W].
                  - tas: True airspeed [m/s].
                  - ROCD: Rate of descent [m/s].
                  - ESF: Energy share factor [-].
                  - limitation: Any performance limitations encountered during the descent (e.g., speed, power limits).
        :rtype: list[float]

        The function determines the Long Range Cruise (LRC) or Maximum Endurance Cruise (MEC) TAS for descent, applying
        the default values when necessary. It checks for speed envelope limitations and adjusts TAS accordingly.
        It calculates the power available, required power, and engine power needed for the descent.

        .. note::
           Power limitations are handled by adjusting TAS and calculating the rate of descent (ROCD).
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        delta = atm.delta(h=h, DeltaTemp=DeltaTemp)
        sigma = atm.sigma(theta=theta, delta=delta)

        temp = theta * const.temp_0

        # LRC = self.OPT.LRC(mass=mass, h=h, DeltaTemp=DeltaTemp, wS=0)
        LRC = conv.kt2ms(
            self.OPT.getOPTParam("LRC", conv.m2ft(h), mass, DeltaTemp)
        )

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
        maxSpeed = atm.cas2Tas(
            cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma
        )
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
            temp=temp,
            DeltaTemp=DeltaTemp,
            ROCD=ROCD,
            mass=mass,
            Preq=Preq,
            ESF=ESF,
        )
        Peng = Peng_target

        if Pav < Peng:
            limitation += "P"

        return [Pav, Peng, Preq, tas, ROCD, ESF, limitation]

    def decelerationToApproach(self):
        pass

    def approach(
        self,
        h,
        mass,
        DeltaTemp,
        speedLimit=None,
        ROCDDefault=None,
        tasDefault=None,
    ):
        """Computes various parameters for the aircraft approach phase using
        the ARPM model.

        This function calculates key approach parameters, including available and required power, true airspeed (TAS),
        rate of descent (ROCD), and potential performance limitations. The calculations take into account atmospheric
        conditions, altitude, and aircraft mass.

        :param h: Altitude above sea level [m].
        :param mass: Aircraft weight [kg].
        :param DeltaTemp: Deviation from the International Standard Atmosphere (ISA) temperature [K].
        :param speedLimit: Optional parameter to apply speed limits. Use {"applyLimit", None}.
        :param ROCDDefault: Default rate of climb or descent [m/s], optional.
        :param tasDefault: Optional true airspeed (TAS) [m/s].
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type speedLimit: str, optional.
        :type ROCDDefault: float, optional.
        :type tasDefault: float, optional.

        :returns: A list of computed values:
                  - Pav: Available power [W].
                  - Peng: Engine power [W].
                  - Preq: Required power [W].
                  - tas: True airspeed [m/s].
                  - ROCD: Rate of descent [m/s].
                  - ESF: Energy share factor [-].
                  - limitation: Any performance limitations encountered during the approach (e.g., speed, power limits).
        :rtype: list[float]

        The function determines the Maximum Endurance Cruise (MEC) TAS for the approach, applying the default
        values when necessary. It checks for speed envelope limitations and adjusts TAS accordingly. It calculates
        the power available, required power, and engine power needed for the approach.

        .. note::
           Power limitations are handled by adjusting TAS and calculating the rate of descent (ROCD).
        """

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        delta = atm.delta(h=h, DeltaTemp=DeltaTemp)
        sigma = atm.sigma(theta=theta, delta=delta)

        temp = theta * const.temp_0

        # MEC = self.OPT.MEC(mass=mass, h=h, DeltaTemp=DeltaTemp, wS=0)
        MEC = conv.kt2ms(
            self.OPT.getOPTParam("MEC", conv.m2ft(h), mass, DeltaTemp)
        )

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
        maxSpeed = atm.cas2Tas(
            cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma
        )
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
            temp=temp,
            DeltaTemp=DeltaTemp,
            ROCD=ROCD,
            mass=mass,
            Preq=Preq,
            ESF=ESF,
        )
        Peng = Peng_target

        if Pav < Peng:
            limitation += "P"

        return [Pav, Peng, Preq, tas, ROCD, ESF, limitation]

    def decelerationToFinalApproach(self):
        pass

    def finalApproach(
        self,
        h,
        mass,
        DeltaTemp,
        speedLimit=None,
        ROCDDefault=None,
        tasDefault=None,
    ):
        """Computes various parameters for the final approach phase using the
        ARPM model.

        This function calculates key final approach parameters, including available and required power, true airspeed (TAS),
        rate of descent (ROCD), and potential performance limitations. The calculations take into account atmospheric
        conditions, altitude, and aircraft mass.

        :param h: Altitude above sea level [m].
        :param mass: Aircraft weight [kg].
        :param DeltaTemp: Deviation from the International Standard Atmosphere (ISA) temperature [K].
        :param speedLimit: Optional parameter to apply speed limits. Use {"applyLimit", None}.
        :param ROCDDefault: Default rate of climb or descent [m/s], optional.
        :param tasDefault: Optional true airspeed (TAS) [m/s].
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type speedLimit: str, optional.
        :type ROCDDefault: float, optional.
        :type tasDefault: float, optional.

        :returns: A list of computed values:
                  - Pav: Available power [W].
                  - Peng: Engine power [W].
                  - Preq: Required power [W].
                  - tas: True airspeed [m/s].
                  - ROCD: Rate of descent [m/s].
                  - ESF: Energy share factor [-].
                  - limitation: Any performance limitations encountered during the final approach (e.g., speed, power limits).
        :rtype: list[float]

        The function sets the default true airspeed (TAS) for the final approach to 30 knots (converted to meters per second),
        or uses a specified value if provided. It also sets a default rate of descent (ROCD) to -200 feet per minute, or
        takes an optional value if available. The function checks the speed envelope for limitations, adjusts TAS if necessary,
        and calculates the required and available power. If there are power limitations, they are flagged in the output.

        .. note::
           This function uses the constant TAS evolution for the final approach and computes engine power and rate of descent
           based on available power and atmospheric conditions.
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
        maxSpeed = atm.cas2Tas(
            cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma
        )
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
            temp=temp,
            DeltaTemp=DeltaTemp,
            ROCD=ROCD,
            mass=mass,
            Preq=Preq,
            ESF=ESF,
        )
        Peng = Peng_target

        if Pav < Peng:
            limitation += "P"

        return [Pav, Peng, Preq, tas, ROCD, ESF, limitation]

    def decelerationToLanding(self):
        pass

    def landing(self, h, mass, DeltaTemp, ROCDDefault=None):
        """Computes various parameters for the landing phase using the ARPM
        model.

        This function calculates key landing parameters, including available and required power, true airspeed (TAS),
        rate of descent (ROCD), and potential performance limitations. The calculations take into account atmospheric
        conditions, altitude, and aircraft mass.

        :param h: Altitude above sea level [m].
        :param mass: Aircraft weight [kg].
        :param DeltaTemp: Deviation from the International Standard Atmosphere (ISA) temperature [K].
        :param ROCDDefault: Default rate of descent [m/s], optional.
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.
        :type ROCDDefault: float, optional.

        :returns: A list of computed values:
                  - Pav: Available power [W].
                  - Peng: Engine power [W].
                  - Preq: Required power [W].
                  - tas: True airspeed (set to 0 for landing) [m/s].
                  - ROCD: Rate of descent [m/s].
                  - ESF: Energy share factor [-].
                  - limitation: Any performance limitations encountered during the landing (e.g., power limits).
        :rtype: list[float]

        This function sets the rate of descent (ROCD) to a default value of -100 feet per minute or an optional value
        if provided. The true airspeed (TAS) is set to 0 for landing calculations. The function checks if available
        power meets the required power and flags any limitations in performance, such as power limitations.

        .. note::
           The ESF (Energy Share Factor) is calculated for constant TAS during landing, and engine power is computed
           accordingly.
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
            temp=temp,
            DeltaTemp=DeltaTemp,
            ROCD=ROCD,
            mass=mass,
            Preq=Preq,
            ESF=ESF,
        )
        Peng = Peng_target

        if Pav < Peng:
            limitation += "P"

        return [Pav, Peng, Preq, tas, ROCD, ESF, limitation]

    def hover(self, h, mass, DeltaTemp):
        """Computes various parameters for the hover phase using the ARPM
        model.

        This function calculates key hover parameters, including available and required power, true airspeed (TAS),
        and any potential performance limitations. The calculations take into account atmospheric conditions, altitude,
        and aircraft mass.

        :param h: Altitude above sea level [m].
        :param mass: Aircraft weight [kg].
        :param DeltaTemp: Deviation from the International Standard Atmosphere (ISA) temperature [K].
        :type h: float.
        :type mass: float.
        :type DeltaTemp: float.

        :returns: A list of computed values:
                  - Pav: Available power [W].
                  - Peng: Engine power [W].
                  - Preq: Required power [W].
                  - tas: True airspeed (set to 0 for hover) [m/s].
                  - ROCD: Rate of climb or descent (set to 0 for hover) [m/s].
                  - ESF: Energy share factor [-], set to 0 for hover.
                  - limitation: Any performance limitations encountered during the hover (e.g., power limits).
        :rtype: list[float]

        This function calculates the hover parameters where both true airspeed (TAS) and rate of climb or descent (ROCD)
        are set to 0. The available and required power are computed, and performance limitations such as power limitations
        are flagged if applicable.

        .. note::
           The energy share factor (ESF) is not applicable during hover, and is therefore set to 0.
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
        """Computes various parameters for different flight phases using the
        ARPM model.

        This function calculates the available power (Pav), engine power (Peng), required power (Preq),
        true airspeed (TAS), rate of climb or descent (ROCD), energy share factor (ESF), and any limitations
        encountered during the specified flight phase. The phases include climb, cruise, descent, and hover.
        The calculations take into account atmospheric conditions, altitude, and aircraft mass.

        :param h: Altitude above sea level [m].
        :param mass: Aircraft weight [kg].
        :param phase: The flight phase being calculated, one of {"Climb", "Cruise", "Descent", "Hover"}.
        :param DeltaTemp: Deviation from the International Standard Atmosphere (ISA) temperature [K].
        :param rating: Engine rating {MTKF, MCNT, ARPM}, default is ARPM [-].
        :param speedLimit: Optional parameter to apply speed limits. Use {"applyLimit", None}.
        :param ROCDDefault: Default rate of climb or descent [m/s], optional.
        :param tasDefault: Optional true airspeed (TAS) [m/s].
        :type h: float.
        :type mass: float.
        :type phase: str.
        :type DeltaTemp: float.
        :type rating: str, optional.
        :type speedLimit: str, optional.
        :type ROCDDefault: float, optional.
        :type tasDefault: float, optional.

        :returns: A list of computed values:
                  - Pav: Available power [W].
                  - Peng: Engine power [W].
                  - Preq: Required power [W].
                  - tas: True airspeed [m/s].
                  - ROCD: Rate of climb or descent [m/s].
                  - ESF: Energy share factor [-].
                  - limitation: Any performance limitations encountered (e.g., power, speed limits).
        :rtype: list[float]

        The function determines the appropriate flight phase and computes the required and available power, TAS,
        and ROCD accordingly. It handles various flight phases, including:

        - **Climb**: For altitudes ≤ 5 meters, it uses the takeoff ARPM procedure. For altitudes > 5 meters, it uses the climb procedure.
        - **Cruise**: Computes cruise parameters.
        - **Descent**: Handles descent, approach, final approach, and landing depending on the altitude.
          - For h ≥ 500 feet, descent parameters are computed.
          - For 150 feet ≤ h < 500 feet, the approach procedure is used.
          - For 5 feet ≤ h < 150 feet, the final approach is computed.
          - For h < 5 feet, landing parameters are computed.
        - **Hover**: Computes hover parameters.

        .. note::
           Power limitations, speed envelope constraints, and other performance-related limitations are flagged and returned
           as part of the limitation output.
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
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = (
                    self.finalApproach(
                        h=h,
                        mass=mass,
                        DeltaTemp=DeltaTemp,
                        speedLimit=speedLimit,
                        ROCDDefault=ROCDDefault,
                        tasDefault=tasDefault,
                    )
                )
            elif h < conv.ft2m(5):
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.landing(
                    h=h,
                    mass=mass,
                    DeltaTemp=DeltaTemp,
                    ROCDDefault=ROCDDefault,
                )

        elif phase == "Hover":
            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.hover(
                h=h, mass=mass, DeltaTemp=DeltaTemp
            )

        return [Pav, Peng, Preq, tas, ROCD, ESF, limitation]


class PTD(BADAH):
    """This class implements the PTD file creator for BADAH aircraft following
    BADAH manual.

    :param AC: Aircraft object {BADAH}.
    :type AC: badaHAircraft.
    """

    def __init__(self, AC):
        super().__init__(AC)

        self.flightEnvelope = FlightEnvelope(AC)
        self.ARPM = ARPM(AC)

    def create(self, saveToPath, DeltaTemp):
        """Creates a BADAH PTD file based on aircraft performance data at
        different mass levels, altitudes, and temperatures.

        This function calculates performance data for three different mass levels (low, medium, high), at various
        altitudes, and for different temperature deviations from the International Standard Atmosphere (ISA).
        It computes climb, cruise, descent, and hover performance data, then saves this information into a
        PTD file.

        :param saveToPath: Path to the directory where the PTD file should be saved.
        :param DeltaTemp: Deviation from the ISA temperature [K].
        :type saveToPath: str.
        :type DeltaTemp: float.
        :returns: None
        :rtype: None

        The function generates data for different flight phases (climb, cruise, descent, and hover) for the
        BADA engine ratings (ARPM, MTKF, MCNT). It stores the computed data in lists and then calls
        the `save2PTD` method to save the data into a PTD file.
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
        """Saves the computed performance data to a BADAH PTD file.

        This function saves the performance data generated during different
        flight phases (climb, cruise, descent, hover) and for different engine
        ratings (ARPM, MTKF, MCNT) into a PTD file. The file is named based on
        the aircraft name and ISA deviation.

        :param saveToPath: Path to the directory where the PTD file should be
            saved.
        :param CLList_ARPM: List of climb data for the BADA ARPM rating.
        :param CLList_MTKF: List of climb data for the BADA MTKF rating.
        :param CLList_MCNT: List of climb data for the BADA MCNT rating.
        :param CRList: List of cruise data.
        :param DESList: List of descent data.
        :param HOVERList: List of hover data.
        :param DeltaTemp: Deviation from ISA temperature [K].
        :type saveToPath: str.
        :type CLList_ARPM: list.
        :type CLList_MTKF: list.
        :type CLList_MCNT: list.
        :type CRList: list.
        :type DESList: list.
        :type HOVERList: list.
        :type DeltaTemp: float.
        :returns: None
        :rtype: None
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
        file.write(
            "=============================\n=============================\n\n"
        )
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
        """Calculates the BADAH PTD (Performance Table Data) for the climb
        phase.

        This function computes the aircraft's performance parameters during
        the climb phase for each altitude level in the given altitude list.
        Parameters such as temperature, pressure, density, true airspeed
        (TAS), and rate of climb/descent (ROCD) are calculated and returned in
        a list format that can be used for generating PTD files.

        :param mass: Aircraft mass [kg].
        :param altitudeList: List of altitude values [ft].
        :param DeltaTemp: Deviation from ISA temperature [K].
        :param rating: Engine rating, e.g., {MTKF, MCNT, ARPM}.
        :type mass: float.
        :type altitudeList: list of int.
        :type DeltaTemp: float.
        :type rating: str.
        :returns: List of PTD climb data.
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

            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = (
                self.ARPM.ARPMProcedure(
                    phase=phase,
                    h=H_m,
                    DeltaTemp=DeltaTemp,
                    mass=mass,
                    rating=rating,
                )
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
        """Calculates the BADAH PTD (Performance Table Data) for the descent
        phase.

        This function computes the aircraft's performance parameters during
        the descent phase for each altitude level in the given altitude list.
        It calculates values such as temperature, pressure, density, true
        airspeed (TAS), and rate of descent (ROD), and returns the data in a
        structured list format for PTD file generation.

        :param mass: Aircraft mass [kg].
        :param altitudeList: List of altitude values [ft].
        :param DeltaTemp: Deviation from ISA temperature [K].
        :type mass: float.
        :type altitudeList: list of int.
        :type DeltaTemp: float.
        :returns: List of PTD descent data.
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

            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = (
                self.ARPM.ARPMProcedure(
                    phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=mass
                )
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
        """Calculates the BADAH PTD (Performance Table Data) for the cruise
        phase.

        This function computes the aircraft's performance parameters during
        the cruise phase for each altitude level in the given altitude list.
        Key performance metrics like temperature, pressure, density, TAS, and
        fuel consumption are calculated and stored in a structured list for
        PTD file generation.

        :param mass: Aircraft mass [kg].
        :param altitudeList: List of altitude values [ft].
        :param DeltaTemp: Deviation from ISA temperature [K].
        :type mass: float.
        :type altitudeList: list of int.
        :type DeltaTemp: float.
        :returns: List of PTD cruise data.
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

            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = (
                self.ARPM.ARPMProcedure(
                    phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=mass
                )
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
        """Calculates the BADAH PTD (Performance Table Data) for the hover
        phase.

        This function computes the aircraft's performance parameters during
        the hover phase for each altitude level in the given altitude list. It
        calculates values like temperature, pressure, density, and fuel
        consumption during hover and returns the data in a structured list
        format for PTD generation.

        :param mass: Aircraft mass [kg].
        :param altitudeList: List of altitude values [ft].
        :param DeltaTemp: Deviation from ISA temperature [K].
        :type mass: float.
        :type altitudeList: list of int.
        :type DeltaTemp: float.
        :returns: List of PTD hover data.
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

            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = (
                self.ARPM.ARPMProcedure(
                    phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=mass
                )
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
    """This class implements the PTF file creator for BADAH aircraft following
    BADAH manual.

    :param AC: Aircraft object {BADAH}.
    :type AC: badaHAircraft.
    """

    def __init__(self, AC):
        super().__init__(AC)

        self.flightEnvelope = FlightEnvelope(AC)
        self.ARPM = ARPM(AC)

    def create(self, saveToPath, DeltaTemp):
        """Creates the BADAH PTF and saves it to the specified directory.

        :param saveToPath: Path to the directory where the PTF should be
            stored.
        :param DeltaTemp: Deviation from ISA temperature [K].
        :type saveToPath: str
        :type DeltaTemp: float
        :returns: None
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
        self,
        saveToPath,
        altitudeList,
        CLList,
        CRList,
        DESList,
        DeltaTemp,
        massList,
    ):
        """Saves the BADAH performance data to a PTF format.

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
            file, adhering to the BADAH performance file format.
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
            "BADA PERFORMANCE FILE                                        %s\n\n"
            % (d3)
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
        """Calculates the BADAH PTF for the CRUISE phase of flight.

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
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = (
                    self.ARPM.ARPMProcedure(
                        phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=mass
                    )
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

        CRList = [
            TAS_CR_complet,
            FF_CR_LO_complet,
            FF_CR_NOM_complet,
            FF_CR_HI_complet,
        ]

        return CRList

    def PTF_climb(self, massList, altitudeList, DeltaTemp, rating):
        """Calculates the BADAH PTF for the CLIMB phase of flight.

        :param massList: List of aircraft mass levels in kilograms [kg].
        :param altitudeList: List of aircraft altitudes in feet [ft].
        :param DeltaTemp: Deviation from ISA temperature in Kelvin [K].
        :param rating: Engine rating {MTKF, MCNT, ARPM} [-].
        :type massList: list
        :type altitudeList: list of int
        :type DeltaTemp: float
        :type rating: str
        :returns: List of PTF CLIMB data, including True Airspeed, Rates of
            Climb, and Fuel Flow for each mass level.
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
                phase=phase,
                h=H_m,
                DeltaTemp=DeltaTemp,
                mass=massNominal,
                rating=rating,
            )

            CP = self.CP(Peng=Peng)
            ff_nominal = self.ff(delta=delta, CP=CP) * 60  # [kg/min]

            ROC = []
            for mass in massList:
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = (
                    self.ARPM.ARPMProcedure(
                        phase=phase,
                        h=H_m,
                        DeltaTemp=DeltaTemp,
                        mass=mass,
                        rating=rating,
                    )
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
        """Calculates the BADAH PTF for the DESCENT phase of flight.

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
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = (
                    self.ARPM.ARPMProcedure(
                        phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=mass
                    )
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
    """This class encapsulates the BADAH performance model for an aircraft,
    extending the BADAH base class.

    :param badaVersion: The version of the BADAH model being used.
    :param acName: The ICAO designation or name of the aircraft.
    :param filePath: (Optional) Path to the BADAH XML file. If not provided, a
        default path is used.
    :param allData: (Optional) Dataframe containing pre-loaded aircraft data,
        typically used to initialize the aircraft parameters without needing
        to parse XML files.
    :type badaVersion: str
    :type acName: str
    :type filePath: str, optional
    :type allData: pd.DataFrame, optional This class initializes the
        aircraft's performance model using data from a dataframe or by reading
        from XML files in the BADAH format.
    """

    def __init__(self, badaVersion, acName, filePath=None, allData=None):
        """Initializes the BADAHAircraft class by loading aircraft-specific
        data.

        - If `allData` is provided and contains the aircraft's information, it will be used to
          initialize various parameters such as engine type, mass, thrust settings, and performance
          data.
        - If the aircraft is not found in `allData`, the class will search for the corresponding
          BADAH XML file or synonym file (if applicable) in the specified or default file path.
        - Once the aircraft data is found, the class initializes various performance modules such
          as the flight envelope, aerodynamic model, and performance optimizations.

        :param badaVersion: Version of the BADAH model (e.g., "1.1").
        :param acName: ICAO aircraft designation or model name.
        :param filePath: (Optional) Custom file path to load the aircraft data. If not provided,
                         a default directory is used.
        :param allData: (Optional) Dataframe containing pre-loaded aircraft data for initialization.
        """
        super().__init__(self)

        self.BADAFamily = BadaFamily(BADAH=True)
        self.BADAFamilyName = "BADAH"
        self.BADAVersion = badaVersion
        self.acName = acName

        if filePath is None:
            self.filePath = configuration.getBadaVersionPath(
                badaFamily="BADAH", badaVersion=badaVersion
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
            self.MR_radius = configuration.safe_get(
                filtered_df, "MR_radius", None
            )
            self.MR_Speed = configuration.safe_get(
                filtered_df, "MR_Speed", None
            )
            self.cpr = configuration.safe_get(filtered_df, "cpr", None)
            self.n_eng = configuration.safe_get(filtered_df, "n_eng", None)
            self.P0 = configuration.safe_get(filtered_df, "P0", None)
            self.cf = configuration.safe_get(filtered_df, "cf", None)
            self.Pmax_ = configuration.safe_get(filtered_df, "Pmax_", None)
            self.cpa = configuration.safe_get(filtered_df, "cpa", None)
            self.hmo = configuration.safe_get(filtered_df, "hmo", None)
            self.vne = configuration.safe_get(filtered_df, "vne", None)
            self.MTOW = configuration.safe_get(filtered_df, "MTOW", None)
            self.OEW = configuration.safe_get(filtered_df, "OEW", None)
            self.MFL = configuration.safe_get(filtered_df, "MFL", None)
            self.MREF = configuration.safe_get(filtered_df, "MREF", None)
            self.MPL = configuration.safe_get(filtered_df, "MPL", None)
            self.VMO = configuration.safe_get(filtered_df, "VMO", None)
            self.MMO = configuration.safe_get(filtered_df, "MMO", None)

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
            synonymFile = os.path.join(self.filePath, "SYNONYM.xml")
            if os.path.isfile(synonymFile):
                self.synonymFileAvailable = True

                # if SYNONYM exist - look for synonym based on defined acName
                self.SearchedACName = Parser.parseSynonym(
                    self.filePath, acName
                )

                # if cannot find - look for full name (in sub folder names) based on acName (may not be ICAO designator)
                if self.SearchedACName is None:
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
                        self.SearchedACName,
                        self.SearchedACName,
                    )
                    + ".xml"
                )
                OPTFilePath = os.path.join(self.filePath, acName)

                if os.path.isfile(acXmlFile):
                    self.ACModelAvailable = True

                    ACparsed_df = Parser.parseXML(
                        self.filePath, self.SearchedACName
                    )

                    self.OPTFilePath = OPTFilePath

                    self.model = configuration.safe_get(
                        ACparsed_df, "model", None
                    )
                    self.engineType = configuration.safe_get(
                        ACparsed_df, "engineType", None
                    )
                    self.engines = configuration.safe_get(
                        ACparsed_df, "engines", None
                    )
                    self.WTC = configuration.safe_get(ACparsed_df, "WTC", None)
                    self.ICAO = configuration.safe_get(
                        ACparsed_df, "ICAO", None
                    )
                    self.MR_radius = configuration.safe_get(
                        ACparsed_df, "MR_radius", None
                    )
                    self.MR_Speed = configuration.safe_get(
                        ACparsed_df, "MR_Speed", None
                    )
                    self.cpr = configuration.safe_get(ACparsed_df, "cpr", None)
                    self.n_eng = configuration.safe_get(
                        ACparsed_df, "n_eng", None
                    )
                    self.P0 = configuration.safe_get(ACparsed_df, "P0", None)
                    self.cf = configuration.safe_get(ACparsed_df, "cf", None)
                    self.Pmax_ = configuration.safe_get(
                        ACparsed_df, "Pmax_", None
                    )
                    self.cpa = configuration.safe_get(ACparsed_df, "cpa", None)
                    self.hmo = configuration.safe_get(ACparsed_df, "hmo", None)
                    self.vne = configuration.safe_get(ACparsed_df, "vne", None)
                    self.MTOW = configuration.safe_get(
                        ACparsed_df, "MTOW", None
                    )
                    self.OEW = configuration.safe_get(ACparsed_df, "OEW", None)
                    self.MFL = configuration.safe_get(ACparsed_df, "MFL", None)
                    self.MREF = configuration.safe_get(
                        ACparsed_df, "MREF", None
                    )
                    self.MPL = configuration.safe_get(ACparsed_df, "MPL", None)
                    self.VMO = configuration.safe_get(ACparsed_df, "VMO", None)
                    self.MMO = configuration.safe_get(ACparsed_df, "MMO", None)

                    self.flightEnvelope = FlightEnvelope(self)
                    self.OPT = Optimization(self)
                    self.ARPM = ARPM(self)
                    self.PTD = PTD(self)
                    self.PTF = PTF(self)

                else:
                    # AC name cannot be found
                    raise ValueError(
                        acName + " Cannot be found at: " + self.filePath
                    )

    def __str__(self):
        return f"(BADAH, AC_name: {self.acName}, searched_AC_name: {self.SearchedACName}, model_ICAO: {self.ICAO}, ID: {id(self.AC)})"
