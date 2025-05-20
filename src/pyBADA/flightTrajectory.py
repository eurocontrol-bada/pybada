"""
Generic flight trajectory module
"""

import os
import datetime
import pandas as pd
import simplekml

from pyBADA import conversions as conv


class FlightTrajectory:
    """This class implements the flight trajectory module and handles all the
    operations on the flight trajectory."""

    def __init__(self):
        self.flightData = {}

    def createFT(self):
        """Creates and returns an empty DataFrame for storing aircraft flight
        trajectory data. This DataFrame includes various flight parameters
        such as altitude, speed, fuel consumption, acceleration, and more. The
        columns are predefined to match the data typically recorded during a
        flight.

        :param AC: Aircraft object from the BADA family (BADA3/4/H/E).
        :param Hp: Pressure altitude [ft].
        :param TAS: True Airspeed [kt].
        :param CAS: Calibrated Airspeed [kt].
        :param M: Mach number [-].
        :param ROCD: Rate of Climb/Descent [ft/min].
        :param FUEL: Fuel consumption rate [kg/s].
        :param P: Power output of the engines [W].
        :param slope: Flight path slope [degrees].
        :param acc: Aircraft acceleration [m/s^2].
        :param THR: Thrust produced by the engines [N].
        :param config: Aerodynamic configuration of the aircraft (e.g., clean, takeoff, landing).
        :param HLid: High Lift Device deployment level [-].
        :param LG: Landing gear deployment status (e.g., retracted, deployed).
        :param mass: Aircraft mass [kg].
        :param LAT: Geographical latitude [degrees].
        :param LON: Geographical longitude [degrees].
        :param HDG: Aircraft heading [degrees].
        :param time: Elapsed flight time [s].
        :param dist: Distance traveled [NM].
        :param comment: Optional comment describing the trajectory segment.

        :type AC: BadaAircraft {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}.
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
        :type config: str
        :type HLid: float
        :type LG: str
        :type mass: float
        :type LAT: float
        :type LON: float
        :type HDG: float
        :type time: float
        :type dist: float
        :type comment: str

        :returns: An empty DataFrame for flight trajectory data.
        :rtype: pd.DataFrame
        """

        # Create an empty DataFrame with the required flight parameters as columns
        flightTrajectory = pd.DataFrame(
            columns=[
                "Hp",  # Pressure altitude [ft]
                "TAS",  # True Airspeed [kt]
                "CAS",  # Calibrated Airspeed [kt]
                "GS",  # Ground Speed [kt]
                "M",  # Mach number [-]
                "ROCD",  # Rate of Climb/Descent [ft/min]
                "ESF",  # Engine specific fuel consumption [-]
                "FUEL",  # Fuel flow rate [kg/s]
                "FUELCONSUMED",  # Total fuel consumed [kg]
                "Preq",  # Required power [W]
                "Peng",  # Generated power [W]
                "Pav",  # Available power [W]
                "slope",  # Flight path slope [degrees]
                "acc",  # Acceleration [m/s^2]
                "THR",  # Thrust [N]
                "DRAG",  # Drag force [N]
                "config",  # Aircraft aerodynamic configuration (clean, takeoff, etc.)
                "HLid",  # High Lift Device deployment level [-]
                "LG",  # Landing gear deployment status (e.g., up, down)
                "mass",  # Aircraft mass [kg]
                "LAT",  # Geographical latitude [degrees]
                "LON",  # Geographical longitude [degrees]
                "HDGTrue",  # True heading [degrees]
                "HDGMagnetic",  # Magnetic heading [degrees]
                "time",  # Time flown [s]
                "dist",  # Distance traveled [NM]
                "comment",  # Optional comment about the flight segment
                "BankAngle",  # Bank angle during the turn [degrees]
                "ROT",  # Rate of turn [degrees/s]
            ]
        )
        return flightTrajectory

    @staticmethod
    def createFlightTrajectoryDataframe(flight_data):
        """Creates a pandas DataFrame from flight trajectory data, ensuring
        that all lists of data have the same length by padding shorter lists
        with None. This makes sure the resulting DataFrame has equal column
        lengths for each parameter.

        :param flight_data: Dictionary containing flight trajectory data,
            where values are lists of float values representing various
            parameters.
        :type flight_data: dict{list[float]}
        :returns: A pandas DataFrame representing the aircraft's flight
            trajectory.
        :rtype: pandas.DataFrame
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
        flightTrajectory_exploded = flightTrajectory.explode(
            columns_to_explode
        )

        return flightTrajectory_exploded

    def getACList(self):
        """Returns a list of aircraft present in the flight trajectory object.

        :returns: A list of BadaAircraft objects corresponding to the aircraft
            in the current flight trajectory.
        :rtype: list[BadaAircraft]
        """

        return list(self.flightData.keys())

    def addFT(self, AC, flightTrajectory):
        """Adds a flight trajectory for a specific aircraft to the internal
        data structure.

        .. note:: This will overwrite any previously stored flight trajectory for the same aircraft.

        :param AC: Aircraft object (BADA3/4/H/E) whose trajectory is being stored.
        :param flightTrajectory: Pandas DataFrame containing the full flight trajectory for the aircraft.
        :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}
        :type flightTrajectory: pandas.DataFrame
        """

        self.flightData[AC] = flightTrajectory

    def getFT(self, AC):
        """Returns the flight trajectory DataFrame for a specific aircraft.

        :param AC: Aircraft object (BADA3/4/H/E) whose flight trajectory is
            being retrieved.
        :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}
        :returns: A pandas DataFrame containing the flight trajectory data of
            the aircraft.
        :rtype: pandas.DataFrame
        """

        return self.flightData.get(AC)

    def getAllValues(self, AC, parameter):
        """Retrieves all values for a specific parameter from the aircraft's
        flight trajectory.

        :param AC: Aircraft object (BADA3/4/H/E) whose flight data is being
            queried.
        :param parameter: The name of the parameter to retrieve values for
            (e.g., 'altitude', 'speed').
        :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}
        :type parameter: str
        :returns: A list of values corresponding to the specified parameter
            throughout the flight.
        :rtype: list[float]
        """

        values = self.getFT(AC).get(parameter)

        if values is not None:
            return values.tolist()
        else:
            return []

    def getFinalValue(self, AC, parameter):
        """Retrieves the last value for a specific parameter or a list of
        parameters from the aircraft's flight trajectory.

        :param AC: Aircraft object (BADA3/4/H/E) whose flight data is being
            queried.
        :param parameter: The name or list of names of the parameter(s) to
            retrieve the final value(s) for.
        :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}
        :type parameter: list[str] or str
        :returns: The last value (or list of last values) for the specified
            parameter(s).
        :rtype: float or list[float]
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

    def append(self, AC, flightTrajectoryToAppend, overwriteLastRow=False):
        """Appends two consecutive flight trajectories and merges them,
        adjusting cumulative fields such as time, distance, and fuel consumed.
        If the aircraft is not already present, the new trajectory will be
        added.

        If overwriteLastRow is True, the last point of the existing trajectory
        is removed so that the first row of the appended trajectory replaces
        it.

        :param AC: Aircraft object (BADA3/4/H/E) whose trajectory is being
            appended.
        :param flightTrajectoryToAppend: The second flight trajectory to
            append, in the form of a DataFrame.
        :param overwriteLastRow: Flag to indicate whether the last point of
            the existing trajectory should be overwritten.
        :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}
        :type flightTrajectoryToAppend: pd.DataFrame
        :type overwriteLastRow: bool
        """

        # Retrieve the original trajectory
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
            # Determine offset values from the last row of flightTrajectory
            offset_values = {}
            for col in cumulative_columns:
                if (
                    col in flightTrajectory.columns
                    and col in flightTrajectoryToAppend.columns
                ):
                    offset_values[col] = flightTrajectory[col].iloc[-1]

            # If overwriteLastRow is True, remove the last row so it gets replaced by the new data.
            if overwriteLastRow:
                flightTrajectory = flightTrajectory.iloc[:-1]

            # Apply the cumulative addition using the offset values
            for col, offset in offset_values.items():
                flightTrajectoryToAppend[col] = flightTrajectoryToAppend[
                    col
                ].astype(float)
                flightTrajectoryToAppend.loc[:, col] = (
                    flightTrajectoryToAppend[col] + float(offset)
                )

        # Concatenate the original (or modified) trajectory with the new trajectory data
        flightTrajectoryCombined = pd.concat(
            [flightTrajectory, flightTrajectoryToAppend], ignore_index=True
        )

        # Rewrite the original trajectory data
        self.addFT(AC, flightTrajectoryCombined)

    def cut(self, AC, parameter, threshold, direction="BELOW"):
        """Cuts the aircraft's flight trajectory based on a specified
        parameter and threshold value, keeping either the data above or below
        the threshold, depending on the direction.

        .. note:: The data must be sorted by the parameter for the cut to work as expected.

        :param AC: Aircraft object (BADA3/4/H/E) whose flight trajectory is being modified.
        :param parameter: The name of the parameter (e.g., altitude, speed) used for filtering the data.
        :param threshold: The value of the parameter that defines the cutting point.
        :param direction: The direction of the cut. 'ABOVE' removes values above the threshold, 'BELOW' removes values below it.
        :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}
        :type parameter: str
        :type threshold: float
        :type direction: str {'ABOVE', 'BELOW'}
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

    def removeLines(self, AC, numberOfLines=1):
        """Removes from the aircraft's flight trajectory list X number of
        lines,

        :param AC: Aircraft object (BADA3/4/H/E) whose flight trajectory is
            being modified.
        :param numberOfLines: How many lines should be removed from teh end of
            the trajectory
        :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}
        :type numberOfLines: int
        """

        flightTrajectory = self.getFT(AC)

        if numberOfLines <= 0:
            return

        if numberOfLines < len(flightTrajectory):
            flightTrajectoryCut = flightTrajectory.iloc[:-numberOfLines]
        else:
            flightTrajectoryCut = flightTrajectory.iloc[0:0]

        self.addFT(AC, flightTrajectoryCut)

    def overwriteLastValue(self, AC, parameter, new_value):
        """Overwrites the last value of a specified parameter (column) in the
        aircraft's flight trajectory with a new value.

        :param AC: Aircraft object (BADA3/4/H/E) whose flight trajectory is
            being modified.
        :param parameter: The name of the parameter (column) whose last value
            is to be overwritten.
        :param new_value: The new value to assign to the last entry of the
            specified parameter.
        :type AC: {Bada3Aircraft, Bada4Aircraft, BadaEAircraft, BadaHAircraft}
        :type parameter: str
        :type new_value: Depends on the data type of the column (e.g., int,
            float, str)
        """

        # Retrieve the current flight trajectory
        flightTrajectory = self.getFT(AC)

        # Check if the flight trajectory exists and is not empty
        if flightTrajectory is None or flightTrajectory.empty:
            return

        # Ensure the parameter exists in the DataFrame
        if parameter not in flightTrajectory.columns:
            raise ValueError(
                f"The parameter '{parameter}' is not in the flight trajectory columns."
            )

        # Overwrite the last value in the specified column
        flightTrajectory.loc[flightTrajectory.index[-1], parameter] = new_value

        # Rewrite the modified flight trajectory
        self.addFT(AC, flightTrajectory)

    def save2csv(self, saveToPath, separator=","):
        """Saves the aircraft flight trajectory data into a CSV file with a
        custom header depending on the BADA family. The CSV file will be saved
        with a timestamp in the filename.

        :param saveToPath: Path to the directory where the file should be
            stored.
        :param separator: Separator to be used in the CSV file. Default is a
            comma (',').
        :type saveToPath: str
        :type separator: str
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
        """Saves the aircraft flight trajectory data into an Excel (.xlsx)
        file with a custom header depending on the BADA family. The Excel file
        will be saved with a timestamp in the filename.

        :param saveToPath: Path to the directory where the file should be
            stored.
        :type saveToPath: str
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
                flightTrajectory.to_excel(
                    writer, index=False, header=customHeader
                )

    def save2kml(self, saveToPath):
        """Saves the aircraft flight trajectory data into a KML (Keyhole
        Markup Language) file for visualization in tools like Google Earth.
        The KML file is generated with a timestamp in the filename and
        includes aircraft trajectory details with altitude extrusion.

        :param saveToPath: Path to the directory where the file should be
            stored.
        :type saveToPath: str
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

            if not all(
                col in flightTrajectory.columns for col in ["LAT", "LON", "Hp"]
            ):
                print(
                    f"Skipping {AC_ID}: Required columns (LAT, LON, Hp) are missing."
                )
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
            linestring.style.linestyle.color = (
                simplekml.Color.yellow
            )  # Yellow line
            linestring.style.linestyle.width = 3  # Line width in pixels
            linestring.extrude = 1  # Enable altitude extrusion

            # Customize the fill color (extruded space) between the line and the ground
            linestring.style.polystyle.color = simplekml.Color.changealpha(
                "80", simplekml.Color.yellow
            )  # 50% transparent yellow

            # Save the KML file
            kml.save(filename)
