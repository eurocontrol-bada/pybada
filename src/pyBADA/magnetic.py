"""
Magnetic declination module
"""

import json
import bisect
from pyBADA import configuration


class Grid:
    """This class provides methods to calculate the magnetic declination at a given latitude and longitude
    using pre-loaded grid data.

    The grid data is loaded from a JSON file, and the closest grid point is used to determine the
    magnetic declination.
    """

    def __init__(self, inputJSON=None):
        """Initializes the grid with magnetic declination data.

        :param inputJSON: Path to the JSON file containing grid data. Defaults to a pre-configured path.
        :type inputJSON: str, optional
        """

        if inputJSON is None:
            inputJSON = (
                configuration.getDataPath()
                + "/magneticDeclinationGridData.json"
            )

        f = open(inputJSON)
        grid = json.load(f)

        self.gridData = {}
        latitudeList = []
        longitudeList = []
        magneticDeclinationList = []

        for result in grid["result"]:
            latitudeList.append(result["latitude"])
            longitudeList.append(result["longitude"])
            magneticDeclinationList.append(result["declination"])

        self.gridData["LAT"] = latitudeList
        self.gridData["LON"] = longitudeList
        self.gridData["declination"] = magneticDeclinationList

    def getClosestLatitude(self, LAT_target):
        """Finds the closest latitude in the grid to the target latitude.

        :param LAT_target: Target latitude to search for.
        :type LAT_target: float
        :return: The closest latitude from the grid or None if the target is
            out of bounds.
        :rtype: float or None
        """

        latitudeList = sorted(self.gridData["LAT"])

        if LAT_target < latitudeList[0] or LAT_target > latitudeList[-1]:
            return None

        index = bisect.bisect_left(latitudeList, LAT_target)
        if index == 0:
            closest = latitudeList[0]
        elif index == len(latitudeList):
            closest = latitudeList[-1]
        else:
            before = latitudeList[index - 1]
            after = latitudeList[index]
            closest = (
                before if after - LAT_target > LAT_target - before else after
            )

        return closest

    def getClosestLongitude(self, LON_target):
        """Finds the closest longitude in the grid to the target longitude.

        :param LON_target: Target longitude to search for.
        :type LON_target: float
        :return: The closest longitude from the grid or None if the target is
            out of bounds.
        :rtype: float or None
        """

        longitudeList = sorted(self.gridData["LON"])

        if LON_target < longitudeList[0] or LON_target > longitudeList[-1]:
            return None

        index = bisect.bisect_left(longitudeList, LON_target)
        if index == 0:
            closest = longitudeList[0]
        elif index == len(longitudeList):
            closest = longitudeList[-1]
        else:
            before = longitudeList[index - 1]
            after = longitudeList[index]
            closest = (
                before if after - LON_target > LON_target - before else after
            )

        return closest

    def getClosestIdx(self, LAT_target, LON_target):
        """Finds the index of the closest grid point for a given latitude and
        longitude.

        :param LAT_target: Target latitude.
        :param LON_target: Target longitude.
        :type LAT_target: float
        :type LON_target: float
        :return: Index of the closest point in the grid or None if no point is
            found.
        :rtype: int or None
        """

        closestLAT = self.getClosestLatitude(LAT_target=LAT_target)
        closestLON = self.getClosestLongitude(LON_target=LON_target)

        indicesLAT = [
            i
            for i in range(len(self.gridData["LAT"]))
            if self.gridData["LAT"][i] == closestLAT
        ]
        indicesLON = [
            i
            for i in range(len(self.gridData["LON"]))
            if self.gridData["LON"][i] == closestLON
        ]

        for idx in indicesLAT:
            if idx in indicesLON:
                return idx

        return None

    def getMagneticDeclination(self, LAT_target, LON_target):
        """Returns the magnetic declination for the closest grid point to the
        target coordinates.

        :param LAT_target: Target latitude.
        :param LON_target: Target longitude.
        :type LAT_target: float
        :type LON_target: float
        :return: Magnetic declination at the closest grid point or None if no
            point is found.
        :rtype: float or None
        """

        idx = self.getClosestIdx(LAT_target=LAT_target, LON_target=LON_target)

        if idx is None:
            return None
        else:
            magneticDeclination = self.gridData["declination"][idx]

            return magneticDeclination
