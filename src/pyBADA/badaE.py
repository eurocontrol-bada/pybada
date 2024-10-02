# -*- coding: utf-8 -*-
"""
pyBADA
Generic BADAE aircraft performance module
Developped @EUROCONTROL (EIH) by Henrich Glaser-Opitz, PhD
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
from math import sqrt, pow, pi, cos, atan, radians, isnan

from scipy.optimize import fminbound

from pyBADA import constants as const
from pyBADA import conversions as conv
from pyBADA import atmosphere as atm
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


class Parse(object):
    """This class implements the BADAE parsing mechanism to parse xml BADAE files.

    :param filePath: path to the folder with BADAE xml formatted file.
    :param acName: ICAO aircraft designation
    :type filePath: str.
    :type acName: str
    """

    def __init__(self):
        pass

    def parse(self, filePath, badaFamily, badaVersion, acName):
        """This function parses BADAE xml formatted file

        :param filename: path to the BADAE xml formatted file.
        :type filename: str.
        :raises: IOError
        """

        self.filePath = filePath
        self.acName = acName
        self.BADAVersion = badaVersion
        selfBADAFamily = badaFamily

        acXmlFile = (
            os.path.join(filePath, badaFamily, badaVersion, acName, acName) + ".xml"
        )

        try:
            tree = ET.parse(acXmlFile)
            root = tree.getroot()
        except:
            raise IOError(acXmlFile + " not found or in correct format")

        # Parse general aircraft data
        self.model = root.find("model").text  # aircraft model
        self.engineType = root.find("type").text  # aircraft type
        self.engines = root.find("engine").text  # engine type

        self.ICAO_desig = {}  # ICAO designator and WTC
        self.ICAO_desig["designator"] = root.find("ICAO").find("designator").text
        self.ICAO_desig["WTC"] = root.find("ICAO").find("WTC").text

        # Parse Aerodynamic Forces Model
        AFM = root.find("AFM")  # get AFM

        self.MR_radius = float(AFM.find("MR_radius").text)  # Main rotor radius

        CRS = AFM.find("CRS")
        self.crs = []
        for i in CRS.findall("crs"):
            self.crs.append(float(i.text))

        CPreq = AFM.find("CPreq")
        self.cpr = []
        for i in CPreq.findall("cpr"):
            self.cpr.append(float(i.text))

        # Parse engine data
        PFM = root.find("PFM")  # get PFM

        self.n_eng = int(PFM.find("n_eng").text)  # number of engines

        EEM = PFM.find("EEM")  # get EEM
        self.P0 = float(EEM.find("P0").text)

        CVoc = EEM.find("CVoc")
        self.cVoc = []  # Matrix of open-circuit voltage coefficients [V]
        for i in CVoc.findall("cVoc"):
            self.cVoc.append(float(i.text))

        CR0 = EEM.find("CR0")
        self.cR0 = []  # matrix of first internal resistance coefficients [ohm]
        for i in CR0.findall("cR0"):
            self.cR0.append(float(i.text))

        CRi = EEM.find("CRi")
        self.cRi = []  # matrix of second internal resistance coefficients [ohm]
        for i in CRi.findall("cRi"):
            self.cRi.append(float(i.text))

        self.Imax = float(
            EEM.find("Imax").text
        )  # Maximum output current of the battery [A]
        self.Vmin = float(
            EEM.find("Vmin").text
        )  # Minimum operating voltage of the motor [V]
        self.capacity = float(EEM.find("capacity").text)  # Battery capacity [Wh]
        self.eta = float(EEM.find("eta").text)  # Efficiency of the motor [-]

        self.Pmax_ = {}
        # Maximum take-off (MTKF)
        MTKF = EEM.find("MTKF")
        self.Pmax_["MTKF"] = float(MTKF.find("Pmax").text)

        # Maximum continuous (MCNT)
        MCNT = EEM.find("MCNT")
        self.Pmax_["MCNT"] = float(MCNT.find("Pmax").text)

        # Parse Aircraft Limitation Model (ALM)
        ALM = root.find("ALM")  # get ALM
        self.hmo = float(ALM.find("GLM").find("hmo").text)
        self.vne = float(ALM.find("KLM").find("vne").text)
        self.MTOW = float(ALM.find("DLM").find("MTOW").text)
        self.OEW = float(ALM.find("DLM").find("OEW").text)
        self.MFL = float(ALM.find("DLM").find("MFL").text)


class BADAE(Helicopter):
    """This class implements the part of BADAE performance model that will be used in other classes following the BADAE manual.

    :param AC: parsed aircraft.
    :type AC: badaE.Parse.
    """

    def __init__(self, AC):
        Helicopter.__init__(object)

        self.AC = AC

    def RotorSpd_reg(self, tas):
        """This function computes rotor rotational speed regression corresponding to the TAS speed

        :param tas: true airspeed (TAS) [kt].
        :type tas: float.
        :return: rotor rotational speed regression [rad/s].
        :rtype: float.
        """

        rotorSpd_reg = (
            pow(tas, 5) * self.AC.crs[0]
            + pow(tas, 4) * self.AC.crs[1]
            + pow(tas, 3) * self.AC.crs[2]
            + pow(tas, 2) * self.AC.crs[3]
            + pow(tas, 1) * self.AC.crs[4]
            + self.AC.crs[5]
        )

        return rotorSpd_reg

    def RotorSpd(self, tas):
        """This function computes rotor rotational speed corresponding to the TAS speed

        :param tas: true airspeed (TAS) [m s^-1].
        :type tas: float.
        :return: rotor rotational speed [rad/min].
        :rtype: float.
        """

        # check if rotor speed is variable or constant
        if any(self.AC.crs[0:-1]):
            # regresion was performed using TAS in [kt] and rots speed in rad/s
            TAS = conv.ms2kt(tas)

            if TAS <= 0:
                rotorSpd = self.RotorSpd_reg(tas=0.0)
            elif TAS >= 100:
                rotorSpd = self.RotorSpd_reg(tas=100.0)
            else:
                rotorSpd = self.RotorSpd_reg(tas=TAS)

            rotorSpd = rotorSpd * 60 / (2 * pi)
        else:
            rotorSpd = self.AC.crs[-1]

        return rotorSpd

    def TipSpd(self, tas):
        """This function computes rotor blade tip speed corresponding to the TAS speed

        :param tas: true airspeed (TAS) [m s^-1].
        :type tas: float.
        :return: rotor blade tip speed [m s^-1].
        :rtype: float.
        """

        tipSpd = self.RotorSpd(tas=tas) * (2 * pi / 60) * self.AC.MR_radius

        return tipSpd

    def Const2(self, tas):
        """This function computes constant2

        :param tas: true airspeed (TAS) [m s^-1].
        :type tas: float.
        :return: const2 [N].
        :rtype: float.
        """

        const2 = (
            const.rho_0 * pi * pow(self.AC.MR_radius, 2) * pow(self.TipSpd(tas=tas), 2)
        )

        return const2

    def Const3(self, tas):
        """This function computes constant3

        :param tas: true airspeed (TAS) [m s^-1].
        :type tas: float.
        :return: const3 [W].
        :rtype: float.
        """

        const3 = self.Const2(tas=tas) * self.TipSpd(tas=tas)

        return const3

    def mu(self, tas):
        """This function computes the advance ratio

        :param tas: true airspeed (TAS) [m/s].
        :param gamma: flight path angle [rad].
        :type tas: float.
        :type gamma: float.
        :return: mu: advance ratio [-].
        :rtype: float.
        """

        mu = tas / self.TipSpd(tas=tas)

        return mu

    def CT(self, tas, sigma, mass, phi=0.0):
        """This function computes the thrust coefficient

        :param tas: true airspeed (TAS) [m/s]
        :param sigma: Normalised density [-].
        :param mass: aircraft mass [kg].
        :param phi: bank angle [deg].
        :type tas: float.
        :type sigma: float.
        :type mass: float.
        :type phi: float.
        :return: thrust coefficient [-].
        :rtype: float.
        """

        CT = (mass * const.g) / (sigma * self.Const2(tas=tas) * cos(radians(phi)))

        return CT

    def Thrust(self, tas, sigma, mass, phi=0.0):
        """This function computes the thrust

        :param tas: true airspeed (TAS) [m/s]
        :param sigma: Normalised density [-].
        :param mass: aircraft mass [kg].
        :type tas: float.
        :type sigma: float.
        :type mass: float.
        :return: thrust [N].
        :rtype: float.
        """

        CT = self.CT(sigma=sigma, mass=mass, tas=tas, phi=phi)
        thrust = sigma * self.Const2(tas) * CT

        return thrust

    def CPreq(self, mu, CT):
        """This function computes the power required coefficient

        :param mu: advance ratio [-].
        :param Ct: thrust coefficient [-].
        :type mu: float
        :type Ct: float
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

    def Preq(self, tas, sigma, mass, phi=0.0):
        """This function computes the power required

        :param sigma: Normalised density [-].
        :param tas: true airspeed (TAS) [m/s]
        :param gamma: flight path angle [rad]
        :param mass: aircraft mass [kg].
        :param phi: bank angle [deg].
        :type sigma: float.
        :type tas: float.
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
        CT = self.CT(sigma=sigma, mass=mass, tas=tas, phi=phi)
        CPreq = self.CPreq(mu=mu, CT=CT)
        Preq = sigma * self.Const3(tas=tas) * CPreq

        return Preq

    def Peng_target(self, Preq, ROCD, mass, ESF, temp, DeltaTemp):
        """This function computes the targeted required power

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

    def Vocbat(self, SOC):
        """This function computes the battery open-circuit voltage

        :param SOC: State of charge [%]
        :type SOC: float.
        :return: battery open-circuit voltage [V].
        :rtype: float.
        """

        # avoid computations with negative SOC
        if SOC <= 0:
            vocBat = float("Nan")
        else:
            vocBat = (
                self.AC.cVoc[0]
                + self.AC.cVoc[1] * pow(SOC, self.AC.cVoc[2])
                + self.AC.cVoc[3] * SOC / (SOC + 0.1)
                + self.AC.cVoc[4] / (100.1 - SOC)
            )

        return max(vocBat, 0)

    def Vbat(self, SOC, I):
        """This function computes the battery voltage

        :param I: required current [A] .
        :param SOC: State of charge [%]
        :type I: float.
        :type SOC: float.
        :return: battery voltage [V].
        :rtype: float.
        """

        vBat = self.Vocbat(SOC=SOC) - self.Rtbat(SOC=SOC) * I

        return vBat

    def Pbat_fromCurrent(self, SOC, I):
        """This function computes battery power delivered with current I

        :param I: required current [A] .
        :param SOC: State of charge [%]
        :type I: float.
        :type SOC: float.
        :return: battery power [W].
        :rtype: float.
        """

        Vbat = self.Vbat(SOC=SOC, I=I)
        Pbat = Vbat * I

        return Pbat

    def Ibat(self, SOC, P):
        """This function computes battery current required to deliver electrical power P

        :param P: electrical power [W] .
        :param SOC: State of charge [%]
        :type P: float.
        :type SOC: float.
        :return: battery current required [A].
        :rtype: float.
        """

        # avoid computations with negative SOC
        if SOC <= 0:
            Ibat = float("Nan")
        else:
            V0 = self.Vocbat(SOC=SOC)
            Rt = self.Rtbat(SOC=SOC)

            if (V0 * V0 - 4 * P * Rt) < 0:
                Ibat = float("Nan")
            else:
                Ibat = (V0 - sqrt(V0 * V0 - 4 * P * Rt)) / (2 * Rt)

        return Ibat

    def PbatLoss_fromCurrent(self, SOC, I):
        """This function computes power loss in the battery when delivering current I

        :param I: required current [A].
        :param SOC: State of charge [%]
        :type I: float.
        :type SOC: float.
        :return: battery power loss [W].
        :rtype: float.
        """

        # NB: not clear why, but the Vahana software computes the losses using only second internal resistance rather than the total one
        PLoss = self.Ribat(SOC=SOC) * I * I

        return PLoss

    def PbatLoss(self, SOC, P):
        """This function computes power loss in the battery when delivering electrical power P

        :param P: electrical power [W] .
        :param SOC: State of charge [%]
        :type P: float.
        :type SOC: float.
        :return: battery power loss [W].
        :rtype: float.
        """

        I = self.Ibat(P=P, SOC=SOC)
        PbatLoss = self.PbatLoss_fromCurrent(I=I, SOC=SOC)

        return PbatLoss

    def Ribat(self, SOC):
        """This function computes battery second internal resistance

        :param SOC: State of charge [%]
        :type SOC: float.
        :return: battery second internal resistance [ohm].
        :rtype: float.
        """

        Ri = self.AC.cRi[0] + self.AC.cRi[1] * SOC + self.AC.cRi[2] * pow(SOC, 2)

        return Ri

    def R0bat(self, SOC):
        """This function computes battery first internal resistance

        :param SOC: State of charge [%]
        :type SOC: float.
        :return: battery first internal resistance [ohm].
        :rtype: float.
        """

        R0bat = (
            self.AC.cR0[0]
            + self.AC.cR0[1] * pow(SOC, self.AC.cR0[2])
            + self.AC.cR0[3] * SOC / (SOC + 0.1)
            + self.AC.cR0[4] / (100.1 - SOC)
        )

        return R0bat

    def Rtbat(self, SOC):
        """This function computes battery total internal resistance

        :param SOC: State of charge [%]
        :type SOC: float.
        :return: battery total internal resistance [ohm].
        :rtype: float.
        """

        Rtbat = self.R0bat(SOC=SOC) + self.Ribat(SOC=SOC)

        return Rtbat

    def PavBat(self, SOC):
        """This function computes electrical power available from the battery

        :param SOC: State of charge [%]
        :type SOC: float.
        :return: electrical power available [W].
        :rtype: float.
        """

        # Compute the max power limit from Vmin
        Pmax_V = (
            self.AC.Vmin * (self.Vocbat(SOC=SOC) - self.AC.Vmin) / self.Rtbat(SOC=SOC)
        )

        # Compute the max power limit from Imax
        Pmax_I = self.AC.Imax * self.Vbat(SOC=SOC, I=self.AC.Imax)

        # Compute most limiting value
        Pmax = max(Pmax_V, Pmax_I)

        # Replace negative values by zeros
        PavBat = max(Pmax, 0)

        return PavBat

    def PavEng(self, SOC):
        """This function computes mechanical power available from the engine

        :param SOC: State of charge [%]
        :type SOC: float.
        :return: mechanical power available [W].
        :rtype: float.
        """

        Pavbat = self.PavBat(SOC=SOC)
        PavEng = self.AC.eta * Pavbat

        return PavEng

    def Pbat(self, SOC, Preq):
        """This function computes total power consumption in the battery

        :param SOC: State of charge [%]
        :param Preq: mechanical Power required [W]
        :type SOC: float.
        :type Preq: float.
        :return: power consumption [W].
        :rtype: float.
        """

        # compute electrical power
        Pelec = Preq / self.AC.eta
        # compute teh losses in the battery
        Ploss = self.PbatLoss(P=Pelec, SOC=SOC)

        Pbat = Pelec + Ploss

        return Pbat

    def SOCrate(self, SOC, Preq):
        """This function computes SOC rate (SOC = State of charge)

        :param SOC: State of charge [%]
        :param Preq: mechanical Power required [W]
        :type SOC: float.
        :type Preq: float.
        :return: discharge rate [%/h].
        :rtype: float.
        """

        Pbat = self.Pbat(Preq=Preq, SOC=SOC)
        SOCrate = 100 * Pbat / self.AC.capacity

        return SOCrate

    def CPav(self, SOC):
        """This function computes the power available coefficient

        :param SOC: State of charge [%]
        :type SOC: float.
        :return: power available coefficient.
        :rtype: [-].
        """

        Pav = self.PavEng(SOC=SOC)
        CPav = Pav / self.Const3(tas=0.0)

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

    def Pav(self, rating, SOC=100.0):
        """This function computes the power available

        :param rating: throttle setting {MTKF,MCNT}.
        :param SOC: State of charge [%]
        :type rating: str.
        :type SOC: float.
        :return: power available.
        :rtype: [W].
        :raise: ValueError.
        """

        Pmax = self.Pmax(rating=rating)
        CPav = self.CPav(SOC=SOC)
        Pav = min(Pmax, self.Const3(tas=0.0) * CPav)

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

    def ff(self):
        """This function computes the fuel flow

        :return fuel flow [kg/s]
        :rtype float
        """

        if self.AC.type == "ELECTRIC":
            ff = 0.0
        else:
            raise ValueError("Unknown engine type")

        return ff

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


class FlightEnvelope(BADAE):
    """This class is a BADAE aircraft subclass and implements the flight envelope caclulations
    following the BADAE manual.

    :param AC: parsed aircraft.
    :type AC: badaE.Parse.
    """

    def __init__(self, AC):
        BADAE.__init__(self, AC)

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


class Optimization(BADAE):
    """This class implements the BADAE optimization following the BADAE manual.

    :param AC: parsed aircraft.
    :type AC: badaE.Parse.
    """

    def __init__(self, AC):
        BADAE.__init__(self, AC)
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

        MRC = float("Nan")

        return MRC

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

        LRC = float("Nan")

        return LRC

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

        theta = atm.theta(h=h, DeltaTemp=DeltaTemp)
        delta = atm.delta(h=h, DeltaTemp=DeltaTemp)
        sigma = atm.sigma(theta=theta, delta=delta)

        # max TAS speed limitation
        Vmax = atm.cas2Tas(cas=self.flightEnvelope.VMax(), delta=delta, sigma=sigma)

        def f(TAS):
            Preq = self.Preq(sigma=sigma, tas=TAS[0], mass=mass)

            # minimize Preq -> const function
            return Preq

        epsilon = 0.01
        mec = float(
            fminbound(f, x1=np.array([0]), x2=np.array([Vmax + epsilon]), disp=False)
        )

        return mec

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
        :param DeltaTemp: deviation with respect to ISA [K]
        :type optParam: string.
        :type var_1: float.
        :type var_2: float.
        :type DeltaTemp: float.
        """

        filename = self.AC.filePath + "/" + self.AC.acName + "/" + optParam + ".OPT"
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


class ARPM(BADAE):
    """This class is a BADAE aircraft subclass and implements the Airline Procedure Model (ARPM)
    following the BADAE user manual.

    :param AC: parsed aircraft.
    :type AC: badaE.Parse.
    """

    def __init__(self, AC):
        BADAE.__init__(self, AC)
        self.flightEnvelope = FlightEnvelope(AC)
        self.OPT = Optimization(AC)

    def takeoff(self, h, mass, DeltaTemp, rating="ARPM", speedLimit=None):
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
        ROCD = conv.ft2m(100) / 60  # [m/s]

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
            Pav = self.Pav(rating="MTKF")
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
            Pav = self.Pav(rating="MTKF")
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
            Pav = self.Pav(rating="MCNT")
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

    def climb(self, h, mass, DeltaTemp, rating="ARPM", speedLimit=None):
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

        MEC = self.OPT.MEC(mass=mass, h=h, DeltaTemp=DeltaTemp, wS=0)
        # MEC = conv.kt2ms(self.OPT.getOPTParam('MEC', conv.m2ft(h), mass, DeltaTemp))

        # control parameters
        tas = MEC
        ROCD = conv.ft2m(1000) / 60  # [m/s]

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
            Pav = self.Pav(rating="MTKF")
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
            Pav = self.Pav(rating="MTKF")
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
            Pav = self.Pav(rating="MCNT")
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

    def cruise(self, h, mass, DeltaTemp, speedLimit=None):
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
        # LRC = conv.kt2ms(self.OPT.getOPTParam('LRC', conv.m2ft(h), mass, DeltaTemp))
        MEC = self.OPT.MEC(mass=mass, h=h, DeltaTemp=DeltaTemp, wS=0)
        # MEC = conv.kt2ms(self.OPT.getOPTParam('MEC', conv.m2ft(h), mass, DeltaTemp))

        # control parameters
        tas = MEC
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
        Pav = self.Pav(rating="MCNT")
        Peng = min(Preq, Pav)

        if Pav < Peng:
            limitation += "P"

        return [Pav, Peng, Preq, tas, ROCD, ESF, limitation]

    def descent(self, h, mass, DeltaTemp, speedLimit=None):
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
        # LRC = conv.kt2ms(self.OPT.getOPTParam('LRC', conv.m2ft(h), mass, DeltaTemp))
        MEC = self.OPT.MEC(mass=mass, h=h, DeltaTemp=DeltaTemp, wS=0)
        # MEC = conv.kt2ms(self.OPT.getOPTParam('MEC', conv.m2ft(h), mass, DeltaTemp))

        # control parameters
        tas = MEC
        ROCD = conv.ft2m(-500) / 60  # [m/s]

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
            rating="MTKF"
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

    def approach(self, h, mass, DeltaTemp, speedLimit=None):
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

        MEC = self.OPT.MEC(mass=mass, h=h, DeltaTemp=DeltaTemp, wS=0)
        # MEC = conv.kt2ms(self.OPT.getOPTParam('MEC', conv.m2ft(h), mass, DeltaTemp))

        # control parameters
        tas = MEC
        ROCD = conv.ft2m(-300) / 60  # [m/s]

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
            rating="MTKF"
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

    def finalApproach(self, h, mass, DeltaTemp, speedLimit=None):
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
        tas = conv.kt2ms(30)
        ROCD = conv.ft2m(-200) / 60  # [m/s]

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
            rating="MTKF"
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

    def landing(self, h, mass, DeltaTemp):
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
        tas = 0
        ROCD = conv.ft2m(-100) / 60  # [m/s]

        limitation = ""

        ESF = self.esf(flightEvolution="constTAS")

        Pav = self.Pav(rating="MTKF")  # verify if Pav is calualted based on MTKF rating
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

        Pav = self.Pav(rating="MTKF")
        Preq = self.Preq(sigma=sigma, tas=tas, mass=mass)
        Peng = Preq

        if Pav < Peng:
            limitation += "P"

        return [Pav, Peng, Preq, tas, ROCD, ESF, limitation]

    def ARPMProcedure(self, phase, h, mass, DeltaTemp, rating="ARPM", speedLimit=None):
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
                )
            elif h > conv.ft2m(5):
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.climb(
                    h=h,
                    mass=mass,
                    DeltaTemp=DeltaTemp,
                    rating=rating,
                    speedLimit=speedLimit,
                )

        elif phase == "Cruise":
            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.cruise(
                h=h, mass=mass, DeltaTemp=DeltaTemp, speedLimit=speedLimit
            )

        elif phase == "Descent":
            if h >= conv.ft2m(500):
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.descent(
                    h=h, mass=mass, DeltaTemp=DeltaTemp, speedLimit=speedLimit
                )
            elif h < conv.ft2m(500) and h >= conv.ft2m(150):
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.approach(
                    h=h, mass=mass, DeltaTemp=DeltaTemp, speedLimit=speedLimit
                )
            elif h < conv.ft2m(150) and h >= conv.ft2m(5):
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.finalApproach(
                    h=h, mass=mass, DeltaTemp=DeltaTemp, speedLimit=speedLimit
                )
            elif h < conv.ft2m(5):
                [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.landing(
                    h=h, mass=mass, DeltaTemp=DeltaTemp
                )

        elif phase == "Hover":
            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.hover(
                h=h, mass=mass, DeltaTemp=DeltaTemp
            )

        return [Pav, Peng, Preq, tas, ROCD, ESF, limitation]


class PTD(BADAE):
    """This class implements the PTD file creator for BADAE aircraft following BADAE manual.

    :param AC: parsed aircraft.
    :type AC: badaE.Parse.
    """

    def __init__(self, AC):
        BADAE.__init__(self, AC)
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
            self.AC.OEW,
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
        """This function calculates the BADAE PTD data in CLIMB

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
        ff_comlet = []
        ESF_complet = []
        ROCD_complet = []
        gamma_complet = []
        Lim_complet = []

        phase = "Climb"

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)

            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.ARPM.ARPMProcedure(
                phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=mass, rating=rating
            )

            cas = atm.tas2Cas(tas=tas, delta=delta, sigma=sigma)
            M = atm.tas2Mach(v=tas, theta=theta)
            a = atm.aSound(theta=theta)
            FL = h / 100

            ff = self.ff() / 60  # [kg/min]

            temp = theta * const.temp_0
            temp_const = (temp) / (temp - DeltaTemp)
            dhdt = ROCD * temp_const

            if tas == 0:
                if ROCD >= 0:
                    gamma = 90
                else:
                    gamma = -90
            else:
                gamma = conv.rad2deg(atan(dhdt / tas))

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
            ff_comlet,
            ESF_complet,
            ROCD_complet,
            gamma_complet,
            Lim_complet,
        ]

        return CLList

    def PTD_descent(self, mass, altitudeList, DeltaTemp):
        """This function calculates the BADAE PTD data in DESCENT

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
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)

            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.ARPM.ARPMProcedure(
                phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=mass
            )

            cas = atm.tas2Cas(tas=tas, delta=delta, sigma=sigma)
            M = atm.tas2Mach(v=tas, theta=theta)
            a = atm.aSound(theta=theta)
            FL = h / 100

            ff = self.ff() / 60  # [kg/min]

            temp = theta * const.temp_0
            temp_const = (temp) / (temp - DeltaTemp)
            dhdt = ROCD * temp_const
            if tas == 0:
                gamma = -90
            else:
                gamma = conv.rad2deg(atan(dhdt / tas))

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
        """This function calculates the BADAE PTD data in Cruise

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

        phase = "Cruise"

        for h in altitudeList:
            H_m = conv.ft2m(h)  # altitude [m]
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)

            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.ARPM.ARPMProcedure(
                phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=mass
            )

            cas = atm.tas2Cas(tas=tas, delta=delta, sigma=sigma)
            M = atm.tas2Mach(v=tas, theta=theta)
            a = atm.aSound(theta=theta)
            FL = h / 100

            ff = self.ff() / 60  # [kg/min]

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
            ff_comlet,
            ESF_complet,
            ROCD_complet,
            gamma_complet,
            Lim_complet,
        ]

        return CRList

    def PTD_hover(self, mass, altitudeList, DeltaTemp):
        """This function calculates the BADAE PTD data in Cruise

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
            [theta, delta, sigma] = atm.atmosphereProperties(h=H_m, DeltaTemp=DeltaTemp)

            [Pav, Peng, Preq, tas, ROCD, ESF, limitation] = self.ARPM.ARPMProcedure(
                phase=phase, h=H_m, DeltaTemp=DeltaTemp, mass=mass
            )

            cas = atm.tas2Cas(tas=tas, delta=delta, sigma=sigma)
            M = atm.tas2Mach(v=tas, theta=theta)
            a = atm.aSound(theta=theta)
            FL = h / 100

            ff = self.ff() / 60  # [kg/min]

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


class PTF(BADAE):
    """This class implements the PTF file creator for BADAE aircraft following BADAE manual.

    :param AC: parsed aircraft.
    :type AC: badaE.Parse.
    """

    def __init__(self, AC):
        BADAE.__init__(self, AC)
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
            self.AC.OEW,
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

        acModel = self.AC.model

        file = open(filename, "w")
        file.write(
            "BADA PERFORMANCE FILE                                        %s\n\n" % (d3)
        )
        file = open(filename, "a")
        file.write("AC/Type: %s\n\n" % (acModel))
        file.write(
            " Speeds:                      Masses [kg]:             Temperature: ISA%s\n"
            % (ISA)
        )
        file.write(
            " climb   - MEC                low     -    %.0f\n"
            % (proper_round(massList[0]))
        )
        file.write(
            " cruise  - MEC                nominal -    %-4.0f        Max Alt. [ft]:%7d\n"
            % (proper_round(massList[1]), altitudeList[-1])
        )
        file.write(
            " descent - MEC                high    -    %0.f\n"
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
        """This function calculates the BADAE PTF data in CRUISE

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
                    ff.append(self.ff() / 60)  # [kg/min]

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
        """This function calculates the BADAE PTF data in CLIMB

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

            ff_nominal = self.ff() / 60  # [kg/min]

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
        """This function calculates the BADAE PTF data in DESCENT

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

            ff_nominal = self.ff() / 60  # [kg/min]

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


class BadaEAircraft(BADAE, BadaFamily):
    """This class implements the BADAE performance model following the BADAE manual.

    :param filePath: path to the folder with BADAE xml formatted file.
    :param acName: ICAO aircraft designation
    :type filePath: str.
    :type acName: str
    """

    def __init__(self, filePath, badaVersion, acName):
        AC_parsed = Parse()
        self.ACModelAvailable = False

        self.BADAFamily = BadaFamily(BADAE=True)
        self.BADAFamilyName = "BADAE"
        self.BADAVersion = badaVersion

        acXmlFile = (
            os.path.join(filePath, "BADAE", badaVersion, acName, acName) + ".xml"
        )

        if os.path.isfile(acXmlFile):
            self.ACModelAvailable = True

            AC_parsed.parse(
                filePath=filePath,
                badaFamily="BADAE",
                badaVersion=badaVersion,
                acName=acName,
            )

            self.filePath = filePath
            self.acName = acName

            self.model = AC_parsed.model
            self.engineType = AC_parsed.engineType
            self.engines = AC_parsed.engines
            self.ICAO_desig = AC_parsed.ICAO_desig
            self.WTC = AC_parsed.ICAO_desig["WTC"]
            self.ICAO = AC_parsed.ICAO_desig["designator"]
            self.MR_radius = AC_parsed.MR_radius
            self.crs = AC_parsed.crs
            self.cpr = AC_parsed.cpr
            self.n_eng = AC_parsed.n_eng
            self.P0 = AC_parsed.P0
            self.cVoc = AC_parsed.cVoc
            self.cR0 = AC_parsed.cR0
            self.cRi = AC_parsed.cRi
            self.Imax = AC_parsed.Imax
            self.Vmin = AC_parsed.Vmin
            self.capacity = AC_parsed.capacity
            self.eta = AC_parsed.eta
            self.Pmax_ = AC_parsed.Pmax_
            self.hmo = AC_parsed.hmo
            self.vne = AC_parsed.vne
            self.MTOW = AC_parsed.MTOW
            self.OEW = AC_parsed.OEW
            self.MFL = AC_parsed.MFL

            self.VMO = None
            self.MMO = None

            BADAE.__init__(self, AC_parsed)
            self.flightEnvelope = FlightEnvelope(AC_parsed)
            self.OPT = Optimization(AC_parsed)
            self.ARPM = ARPM(AC_parsed)
            self.PTD = PTD(AC_parsed)
            self.PTF = PTF(AC_parsed)

            self.BADAFamily = BadaFamily(BADAE=True)
            self.BADAFamilyName = "BADAE"

        else:
            # AC name cannot be found
            raise ValueError(acName + " Cannot be found")

    def __str__(self):
        return f"(BADAE, AC_name: {self.acName}, searched_AC_name: {self.SearchedACName}, model_ICAO: {self.ICAO}, ID: {id(self.AC)})"
