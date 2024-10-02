# -*- coding: utf-8 -*-
"""
pyBADA
Magnetic declination module
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


import json
import bisect
from pyBADA import configuration


class Grid:
    """This class implements the calculations on the magnetic declination using the grid data

    .. note::

    """

    def __init__(self, inputJSON=None):
        if inputJSON is None:
            inputJSON = (
                configuration.getDataPath() + "/magneticDeclinationGridData.json"
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
            closest = before if after - LAT_target > LAT_target - before else after

        return closest

    def getClosestLongitude(self, LON_target):
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
            closest = before if after - LON_target > LON_target - before else after

        return closest

    def getClosestIdx(self, LAT_target, LON_target):
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
        idx = self.getClosestIdx(LAT_target=LAT_target, LON_target=LON_target)

        if idx is None:
            return None
        else:
            magneticDeclination = self.gridData["declination"][idx]

            return magneticDeclination
