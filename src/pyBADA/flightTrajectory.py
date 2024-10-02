# -*- coding: utf-8 -*-
"""
pyBADA
Generic flight trajectory module
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


import os
import datetime
import pandas as pd
import simplekml

from pyBADA import conversions as conv


class FlightTrajectory:
    """This class implements the flight trajectory module and handles
    all the operations on the flight trajectory

    """

    def __init__(self):
        self.flightData = {}

    def createFT(self):
        """This function creates a flight trajectory and populate it with input data

        :param AC: aircraft {BADA3/4/H/E}
        :param Hp: altitude
        :param TAS: True Air Speed (TAS)
        :param CAS: Calibrated Air Speed (CAS)
        :param M: Mach speed (M)
        :param ROCD: Rate of Climb/Descent
        :param FUEL: fuel consumption
        :param P: Power
        :param slope: trajectory slope
        :param acc: acceleration
        :param THR: thrust
        :param config: aerodynamic configuration
        :param HLid: High Lift Device - Level of deployement
        :param LG: Landing Gear level of deployment
        :param mass: aircraft mass
        :param LAT: Geographical Latitude
        :param LON: Geographical Longitude
        :param HDG: aircraft heading
        :param time: time flown
        :param dist: distance flown
        :param comment: comment describing the trajectory segment
        :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}.
        :type Hp: float
        :type TAS: float
        :type CAS: float
        :type M: float
        :type ROCD: float
        :type FUEL: float
        :type P: float
        :type slope: float
        :type acc: float
        :type THR: float
        :type config: string
        :type HLid: float
        :type LG: string
        :type mass: float
        :type LAT: float
        :type LON: float
        :type HDG: float
        :type time: float
        :type dist: float
        :type comment: string

        :returns: aircraft flight trajectory
        :rtype: dict{list[float]}.
        """

        # Define the empty DataFrame with columns
        flightTrajectory = pd.DataFrame(
            columns=[
                "Hp",
                "TAS",
                "CAS",
                "GS",
                "M",
                "ROCD",
                "ESF",
                "FUEL",
                "FUELCONSUMED",
                "Preq",
                "Peng",
                "Pav",
                "slope",
                "acc",
                "THR",
                "DRAG",
                "config",
                "HLid",
                "LG",
                "mass",
                "LAT",
                "LON",
                "HDGTrue",
                "HDGMagnetic",
                "time",
                "dist",
                "comment",
                "BankAngle",
                "ROT",
            ]
        )
        return flightTrajectory

    @staticmethod
    def createFlightTrajectoryDataframe(flight_data):
        """This function creates a pandas dataframe from the trajectory data in form of a list of data,
        and fill in the missing values with None, to ensure the same size columns

        :param flight_data: trajectory data
        :type flight_data: dict{list[float]}.

        :returns: aircraft flight trajectory
        :rtype: pandas dataframe.
        """
        # Find the maximum length of all lists in the flight data (ignore the Aircraft object)
        max_length = max(
            len(lst) if isinstance(lst, list) else 0
            for key, lst in flight_data.items()
            if key != "Aircraft"
        )

        # Function to pad lists with None to ensure all lists are of equal length
        def pad_list(lst, max_length):
            return lst + [None] * (max_length - len(lst))

        # Pad each list to the same length
        for key in flight_data:
            flight_data[key] = (
                pad_list(flight_data[key], max_length)
                if isinstance(flight_data[key], list)
                else [None] * max_length
            )

        # Convert the padded data to a DataFrame
        flightTrajectory = pd.DataFrame(flight_data)

        # Explode all columns that contain lists
        columns_to_explode = [key for key in flight_data]

        # Explode the DataFrame
        flightTrajectory_exploded = flightTrajectory.explode(columns_to_explode)

        return flightTrajectory_exploded

    def getACList(self):
        """This function return the list of aircraft in the flightTrajectory object

        :returns: list of aircraft in the current flight trajectory object
        :rtype: list[BadaAircraft].
        """

        return list(self.flightData.keys())

    def addFT(self, AC, flightTrajectory):
        """This function adds the flight trajectory based on the aircraft

        .. note::this will overwrite the stored data for the same aircraft

        :param AC: BadaAircraft {BADA3/4/H/E}
        :param flightTrajectory: aircraft full flight trajectory
        :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}.
        :type flightTrajectory: pandas dataframe.
        """

        self.flightData[AC] = flightTrajectory

    def getFT(self, AC):
        """This function returns the flight trajectory based on the aircraft

        :param AC: BadaAircraft {BADA3/4/H/E}
        :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}.

        :returns: aircraft flight trajectory
        :rtype: pandas dataframe.
        """

        return self.flightData.get(AC)

    def getAllValues(self, AC, parameter):
        """This function returns the list of values corresponding to the aircarft trajectory
        and defined parameter name

        :param AC: BadaAircraft {BADA3/4/H/E}
        :param parameter: name of the parameter to search for in the flight trajectory
        :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}.
        :type parameter: string

        :returns: value of a selected parameter for the whole trajectory
        :rtype: list[float]
        """

        values = self.getFT(AC).get(parameter)

        if values is not None:
            return values.tolist()
        else:
            return []

    def getFinalValue(self, AC, parameter):
        """This function returns last value corresponding to the aircarft trajectory
        and defined parameter name

        :param AC: BadaAircraft {BADA3/4/H/E}
        :param parameter: name of the parameter to search for in the flight trajectory
        :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}.
        :type parameter: list[string]

        :returns: final value in the list of a selected parameter for the whole trajectory
        :rtype: float or string
        """

        if isinstance(parameter, list):
            finalValueList = []
            for param in parameter:
                parameterValues = self.getAllValues(AC, param)

                if not parameterValues:
                    finalValueList.append(None)
                else:
                    finalValueList.append(parameterValues[-1])
            return finalValueList

        else:
            parameterValues = self.getAllValues(AC, parameter)
            if not parameterValues:
                return None
            else:
                return self.getAllValues(AC, parameter)[-1]

    def append(self, AC, flightTrajectoryToAppend):
        """This function will append data of 2 consecutive flight trajectories and merge them in terms of time and distance
        if the aircraft is not in the list, the trajectory will be added to the flightTrajectory object

        :param AC: BadaAircraft {BADA3/4/H/E}
        :param flightTrajectoryToAppend: second flight trajectory to combine with original one [dict]
        :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}.
        :type flightTrajectoryToAppend: dict{list[float]}.
        """

        # retrieve the original trajectory
        flightTrajectory = self.getFT(AC)

        # Drop columns with all NaN values from both DataFrames before concatenating
        if flightTrajectory is not None:
            flightTrajectory = flightTrajectory.dropna(axis=1, how="all")

        if flightTrajectoryToAppend is not None:
            flightTrajectoryToAppend = flightTrajectoryToAppend.dropna(
                axis=1, how="all"
            )

        # Make a deep copy of flightTrajectoryToAppend to avoid SettingWithCopyWarning
        flightTrajectoryToAppend = flightTrajectoryToAppend.copy()

        # Handle cumulative columns (time, distance, fuelConsumed)
        cumulative_columns = ["time", "dist", "FUELCONSUMED"]
        if flightTrajectory is not None and not flightTrajectory.empty:
            # For cumulative columns, add the last value of the original to the subsequent values in the appended trajectory
            for col in cumulative_columns:
                if (
                    col in flightTrajectory.columns
                    and col in flightTrajectoryToAppend.columns
                ):
                    last_value = flightTrajectory[col].iloc[
                        -1
                    ]  # Get last value of original trajectory

                    # Ensure both columns are cast to float64 before performing the addition
                    flightTrajectoryToAppend[col] = flightTrajectoryToAppend[
                        col
                    ].astype(float)

                    # Perform the cumulative addition using .loc[]
                    flightTrajectoryToAppend.loc[:, col] = flightTrajectoryToAppend[
                        col
                    ] + float(last_value)

        # Concatenating the two trajectory data
        flightTrajectoryCombined = pd.concat(
            [flightTrajectory, flightTrajectoryToAppend], ignore_index=True
        )

        # rewrite the original trajectory data
        self.addFT(AC, flightTrajectoryCombined)

    def cut(self, AC, parameter, threshold, direction="BELOW"):
        """This function cuts data from aircraft flight trajectory based on input field name and value.

        .. note::The value should be sorted to work as expected.

        :param AC: BadaAircraft {BADA3/4/H/E}
        :param parameter: name of the parameter to take into account
        :param threshold: value of the parameter where the cut shall be performed []
        :param direction: cut above or below set threshold value []
        :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}.
        :type parameter: string.
        :type threshold: float.
        :type direction: string {BELOW/ABOVE}.
        """

        flightTrajectory = self.getFT(AC)

        if direction == "ABOVE":
            flightTrajectoryCut = flightTrajectory[
                flightTrajectory[parameter] < threshold
            ]
        elif direction == "BELOW":
            flightTrajectoryCut = flightTrajectory[
                flightTrajectory[parameter] > threshold
            ]

        self.addFT(AC, flightTrajectoryCut)

    def save2csv(self, saveToPath, separator=","):
        """
        This function saves the trajectory into a CSV file.

        :param saveToPath: Path to directory where the file should be stored.
        :param separator: Separator to be used in the CSV file (only applicable for CSV). Default is a comma (',').
        :type saveToPath: string
        :type separator: string
        :returns: None
        """

        # Get the current time in a suitable format for filenames
        currentTime = "_".join(
            str(datetime.datetime.now()).split(".")[0].split(" ")
        ).replace(":", "-")

        # Create the full directory path
        filepath = os.path.join(saveToPath, f"export_{currentTime}")

        # Check if the directory exists, if not create it
        if not os.path.exists(filepath):
            os.makedirs(filepath)

        # Loop through the aircraft list
        for AC in self.getACList():
            # Get the aircraft ID
            AC_ID = str(id(AC))

            # Flight Trajectory data
            flightTrajectory = self.getFT(AC)

            filename = os.path.join(filepath, f"{AC.acName}_ID{AC_ID}.csv")

            # get custom header based on the BADA Family and some other calculation specificities
            if AC.BADAFamily.BADA3:
                if (
                    "LAT" in flightTrajectory.columns
                    and "LON" in flightTrajectory.columns
                ):
                    customHeader = [
                        "Hp [ft]",
                        "TAS [kt]",
                        "CAS [kt]",
                        "GS [kt]",
                        "M [-]",
                        "acc [m/s^2]",
                        "ROCD [ft/min]",
                        "ESF []",
                        "FUEL [kg/s]",
                        "FUELCONSUMED [kg]",
                        "THR [N]",
                        " DRAG [N]",
                        "t [s]",
                        "d [NM]",
                        "slope [deg]",
                        "m [kg]",
                        "config",
                        "LAT [deg]",
                        "LON [deg]",
                        "HDG True [deg]",
                        "HDG Magnetic [deg]",
                        "bankAngle [deg]",
                        " ROT [deg/s]",
                        "COMMENT",
                    ]
                else:
                    customHeader = [
                        "Hp [ft]",
                        "TAS [kt]",
                        "CAS [kt]",
                        "GS [kt]",
                        "M [-]",
                        "acc [m/s^2]",
                        "ROCD [ft/min]",
                        "ESF []",
                        "FUEL [kg/s]",
                        "FUELCONSUMED [kg]",
                        "THR [N]",
                        " DRAG [N]",
                        "t [s]",
                        "d [NM]",
                        "slope [deg]",
                        "m [kg]",
                        "config",
                        "bankAngle [deg]",
                        " ROT [deg/s]",
                        "COMMENT",
                    ]

            elif AC.BADAFamily.BADA4:
                if (
                    "LAT" in flightTrajectory.columns
                    and "LON" in flightTrajectory.columns
                ):
                    customHeader = [
                        "Hp [ft]",
                        "TAS [kt]",
                        "CAS [kt]",
                        "GS [kt]",
                        "M [-]",
                        "acc [m/s^2]",
                        "ROCD [ft/min]",
                        "ESF []",
                        "FUEL [kg/s]",
                        "FUELCONSUMED [kg]",
                        "THR [N]",
                        " DRAG [N]",
                        "t [s]",
                        "d [NM]",
                        "slope [deg]",
                        "m [kg]",
                        "config",
                        "HLid",
                        "LG",
                        "LAT [deg]",
                        "LON [deg]",
                        "HDG True [deg]",
                        "HDG Magnetic [deg]",
                        "bankAngle [deg]",
                        " ROT [deg/s]",
                        "COMMENT",
                    ]
                else:
                    customHeader = [
                        "Hp [ft]",
                        "TAS [kt]",
                        "CAS [kt]",
                        "GS [kt]",
                        "M [-]",
                        "acc [m/s^2]",
                        "ROCD [ft/min]",
                        "ESF []",
                        "FUEL [kg/s]",
                        "FUELCONSUMED [kg]",
                        "THR [N]",
                        " DRAG [N]",
                        "t [s]",
                        "d [NM]",
                        "slope [deg]",
                        "m [kg]",
                        "config",
                        "HLid",
                        "LG",
                        "bankAngle [deg]",
                        " ROT [deg/s]",
                        "COMMENT",
                    ]

            elif AC.BADAFamily.BADAH:
                if (
                    "LAT" in flightTrajectory.columns
                    and "LON" in flightTrajectory.columns
                ):
                    customHeader = [
                        "Hp [ft]",
                        "TAS [kt]",
                        "CAS [kt]",
                        "GS [kt]",
                        "M [-]",
                        "acc [m/s^2]",
                        "ROCD [ft/min]",
                        "ESF []",
                        "FUEL [kg/s]",
                        "FUELCONSUMED [kg]",
                        "Peng [W]",
                        "Preq [W]",
                        "Pav [W]",
                        "t [s]",
                        "d [NM]",
                        "slope [deg]",
                        "m [kg]",
                        "LAT [deg]",
                        "LON [deg]",
                        "HDG True [deg]",
                        "HDG Magnetic [deg]",
                        "bankAngle [deg]",
                        " ROT [deg/s]",
                        "COMMENT",
                    ]
                else:
                    customHeader = [
                        "Hp [ft]",
                        "TAS [kt]",
                        "CAS [kt]",
                        "GS [kt]",
                        "M [-]",
                        "acc [m/s^2]",
                        "ROCD [ft/min]",
                        "ESF []",
                        "FUEL [kg/s]",
                        "FUELCONSUMED [kg]",
                        "Peng [W]",
                        "Preq [W]",
                        "Pav [W]",
                        "t [s]",
                        "d [NM]",
                        "slope [deg]",
                        "m [kg]",
                        "bankAngle [deg]",
                        " ROT [deg/s]",
                        "COMMENT",
                    ]

            elif AC.BADAFamily.BADAE:
                if (
                    "LAT" in flightTrajectory.columns
                    and "LON" in flightTrajectory.columns
                ):
                    customHeader = [
                        "Hp [ft]",
                        "TAS [kt]",
                        "CAS [kt]",
                        "GS [kt]",
                        "M [-]",
                        "acc [m/s^2]",
                        "ROCD [ft/min]",
                        "ESF []",
                        "Pmec [W]",
                        "Pelc [W]",
                        "Pbat, [W]",
                        "SOCr [%/h]",
                        "SOC [%]",
                        "Ibat [A]",
                        "Vbat [V];",
                        "Vgbat [V]",
                        "t [s]",
                        "d [NM]",
                        "slope [deg]",
                        "m [kg]",
                        "LAT [deg]",
                        "LON [deg]",
                        "HDG True [deg]",
                        "HDG Magnetic [deg]",
                        "bankAngle [deg]",
                        " ROT [deg/s]",
                        "COMMENT",
                    ]
                else:
                    customHeader = [
                        "Hp [ft]",
                        "TAS [kt]",
                        "CAS [kt]",
                        "GS [kt]",
                        "M [-]",
                        "acc [m/s^2]",
                        "ROCD [ft/min]",
                        "ESF []",
                        "Pmec [W]",
                        "Pelc [W]",
                        "Pbat, [W]",
                        "SOCr [%/h]",
                        "SOC [%]",
                        "Ibat [A]",
                        "Vbat [V];",
                        "Vgbat [V]",
                        "t [s]",
                        "d [NM]",
                        "slope [deg]",
                        "m [kg]",
                        "bankAngle [deg]",
                        " ROT [deg/s]",
                        "COMMENT",
                    ]

            # Save to CSV file with custom header and separator
            flightTrajectory.to_csv(
                filename, sep=separator, index=False, header=customHeader
            )

    def save2xlsx(self, saveToPath):
        """
        This function saves the trajectory into a Excel/xlsx file.

        :param saveToPath: Path to directory where the file should be stored.
        :param separator: Separator to be used in the CSV file (only applicable for CSV). Default is a comma (',').
        :type saveToPath: string
        :type separator: string
        :returns: None
        """

        # Get the current time in a suitable format for filenames
        currentTime = "_".join(
            str(datetime.datetime.now()).split(".")[0].split(" ")
        ).replace(":", "-")

        # Create the full directory path
        filepath = os.path.join(saveToPath, f"export_{currentTime}")

        # Check if the directory exists, if not create it
        if not os.path.exists(filepath):
            os.makedirs(filepath)

        # Loop through the aircraft list
        for AC in self.getACList():
            # Get the aircraft ID
            AC_ID = str(id(AC))

            # Flight Trajectory data
            flightTrajectory = self.getFT(AC)

            filename = os.path.join(filepath, f"{AC.acName}_ID{AC_ID}.xlsx")

            # get custom header based on the BADA Family and some other calculation specificities
            if AC.BADAFamily.BADA3:
                if (
                    "LAT" in flightTrajectory.columns
                    and "LON" in flightTrajectory.columns
                ):
                    customHeader = [
                        "Hp [ft]",
                        "TAS [kt]",
                        "CAS [kt]",
                        "GS [kt]",
                        "M [-]",
                        "acc [m/s^2]",
                        "ROCD [ft/min]",
                        "ESF []",
                        "FUEL [kg/s]",
                        "FUELCONSUMED [kg]",
                        "THR [N]",
                        " DRAG [N]",
                        "t [s]",
                        "d [NM]",
                        "slope [deg]",
                        "m [kg]",
                        "config",
                        "LAT [deg]",
                        "LON [deg]",
                        "HDG True [deg]",
                        "HDG Magnetic [deg]",
                        "bankAngle [deg]",
                        " ROT [deg/s]",
                        "COMMENT",
                    ]
                else:
                    customHeader = [
                        "Hp [ft]",
                        "TAS [kt]",
                        "CAS [kt]",
                        "GS [kt]",
                        "M [-]",
                        "acc [m/s^2]",
                        "ROCD [ft/min]",
                        "ESF []",
                        "FUEL [kg/s]",
                        "FUELCONSUMED [kg]",
                        "THR [N]",
                        " DRAG [N]",
                        "t [s]",
                        "d [NM]",
                        "slope [deg]",
                        "m [kg]",
                        "config",
                        "bankAngle [deg]",
                        " ROT [deg/s]",
                        "COMMENT",
                    ]

            elif AC.BADAFamily.BADA4:
                if (
                    "LAT" in flightTrajectory.columns
                    and "LON" in flightTrajectory.columns
                ):
                    customHeader = [
                        "Hp [ft]",
                        "TAS [kt]",
                        "CAS [kt]",
                        "GS [kt]",
                        "M [-]",
                        "acc [m/s^2]",
                        "ROCD [ft/min]",
                        "ESF []",
                        "FUEL [kg/s]",
                        "FUELCONSUMED [kg]",
                        "THR [N]",
                        " DRAG [N]",
                        "t [s]",
                        "d [NM]",
                        "slope [deg]",
                        "m [kg]",
                        "config",
                        "HLid",
                        "LG",
                        "LAT [deg]",
                        "LON [deg]",
                        "HDG True [deg]",
                        "HDG Magnetic [deg]",
                        "bankAngle [deg]",
                        " ROT [deg/s]",
                        "COMMENT",
                    ]
                else:
                    customHeader = [
                        "Hp [ft]",
                        "TAS [kt]",
                        "CAS [kt]",
                        "GS [kt]",
                        "M [-]",
                        "acc [m/s^2]",
                        "ROCD [ft/min]",
                        "ESF []",
                        "FUEL [kg/s]",
                        "FUELCONSUMED [kg]",
                        "THR [N]",
                        " DRAG [N]",
                        "t [s]",
                        "d [NM]",
                        "slope [deg]",
                        "m [kg]",
                        "config",
                        "HLid",
                        "LG",
                        "bankAngle [deg]",
                        " ROT [deg/s]",
                        "COMMENT",
                    ]

            elif AC.BADAFamily.BADAH:
                if (
                    "LAT" in flightTrajectory.columns
                    and "LON" in flightTrajectory.columns
                ):
                    customHeader = [
                        "Hp [ft]",
                        "TAS [kt]",
                        "CAS [kt]",
                        "GS [kt]",
                        "M [-]",
                        "acc [m/s^2]",
                        "ROCD [ft/min]",
                        "ESF []",
                        "FUEL [kg/s]",
                        "FUELCONSUMED [kg]",
                        "Peng [W]",
                        "Preq [W]",
                        "Pav [W]",
                        "t [s]",
                        "d [NM]",
                        "slope [deg]",
                        "m [kg]",
                        "LAT [deg]",
                        "LON [deg]",
                        "HDG True [deg]",
                        "HDG Magnetic [deg]",
                        "bankAngle [deg]",
                        " ROT [deg/s]",
                        "COMMENT",
                    ]
                else:
                    customHeader = [
                        "Hp [ft]",
                        "TAS [kt]",
                        "CAS [kt]",
                        "GS [kt]",
                        "M [-]",
                        "acc [m/s^2]",
                        "ROCD [ft/min]",
                        "ESF []",
                        "FUEL [kg/s]",
                        "FUELCONSUMED [kg]",
                        "Peng [W]",
                        "Preq [W]",
                        "Pav [W]",
                        "t [s]",
                        "d [NM]",
                        "slope [deg]",
                        "m [kg]",
                        "bankAngle [deg]",
                        " ROT [deg/s]",
                        "COMMENT",
                    ]

            elif AC.BADAFamily.BADAE:
                if (
                    "LAT" in flightTrajectory.columns
                    and "LON" in flightTrajectory.columns
                ):
                    customHeader = [
                        "Hp [ft]",
                        "TAS [kt]",
                        "CAS [kt]",
                        "GS [kt]",
                        "M [-]",
                        "acc [m/s^2]",
                        "ROCD [ft/min]",
                        "ESF []",
                        "Pmec [W]",
                        "Pelc [W]",
                        "Pbat, [W]",
                        "SOCr [%/h]",
                        "SOC [%]",
                        "Ibat [A]",
                        "Vbat [V];",
                        "Vgbat [V]",
                        "t [s]",
                        "d [NM]",
                        "slope [deg]",
                        "m [kg]",
                        "LAT [deg]",
                        "LON [deg]",
                        "HDG True [deg]",
                        "HDG Magnetic [deg]",
                        "bankAngle [deg]",
                        " ROT [deg/s]",
                        "COMMENT",
                    ]
                else:
                    customHeader = [
                        "Hp [ft]",
                        "TAS [kt]",
                        "CAS [kt]",
                        "GS [kt]",
                        "M [-]",
                        "acc [m/s^2]",
                        "ROCD [ft/min]",
                        "ESF []",
                        "Pmec [W]",
                        "Pelc [W]",
                        "Pbat, [W]",
                        "SOCr [%/h]",
                        "SOC [%]",
                        "Ibat [A]",
                        "Vbat [V];",
                        "Vgbat [V]",
                        "t [s]",
                        "d [NM]",
                        "slope [deg]",
                        "m [kg]",
                        "bankAngle [deg]",
                        " ROT [deg/s]",
                        "COMMENT",
                    ]

            # Save to xlsx file, since xlsx format doesn’t use a separator
            with pd.ExcelWriter(filename, engine="xlsxwriter") as writer:
                flightTrajectory.to_excel(writer, index=False, header=customHeader)

    def save2kml(self, saveToPath):
        """
        This function saves the trajectory into a KML file.

        :param saveToPath: Path to directory where the file should be stored.
        :param separator: Separator to be used in the CSV file (only applicable for CSV). Default is a comma (',').
        :type saveToPath: string
        :type separator: string
        :returns: None
        """
        # Create a KML object
        kml = simplekml.Kml()

        # Get the current time in a suitable format for filenames
        currentTime = "_".join(
            str(datetime.datetime.now()).split(".")[0].split(" ")
        ).replace(":", "-")

        # Create the full directory path
        filepath = os.path.join(saveToPath, f"export_{currentTime}")

        # Check if the directory exists, if not create it
        if not os.path.exists(filepath):
            os.makedirs(filepath)

        # Loop through the aircraft list
        for AC in self.getACList():
            # Get the aircraft ID
            AC_ID = str(id(AC))

            # Flight Trajectory data
            flightTrajectory = self.getFT(AC)

            if not all(col in flightTrajectory.columns for col in ["LAT", "LON", "Hp"]):
                print(f"Skipping {AC_ID}: Required columns (LAT, LON, Hp) are missing.")
                continue

            filename = os.path.join(filepath, f"{AC.acName}_ID{AC_ID}.kml")

            # Create a LineString for each aircraft's trajectory
            linestring = kml.newlinestring(name=f"{AC.acName} Trajectory")
            linestring.coords = [
                (row["LON"], row["LAT"], conv.ft2m(row["Hp"]))
                for _, row in flightTrajectory.iterrows()
            ]
            linestring.altitudemode = (
                simplekml.AltitudeMode.absolute
            )  # Set altitude mode to absolute

            # Customize the line style for altitude extrusion and color (Yellow)
            linestring.style.linestyle.color = simplekml.Color.yellow  # Yellow line
            linestring.style.linestyle.width = 3  # Line width in pixels
            linestring.extrude = 1  # Enable altitude extrusion

            # Customize the fill color (extruded space) between the line and the ground
            linestring.style.polystyle.color = simplekml.Color.changealpha(
                "80", simplekml.Color.yellow
            )  # 50% transparent yellow

            # Save the KML file
            kml.save(filename)
