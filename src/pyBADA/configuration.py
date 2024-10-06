# -*- coding: utf-8 -*-
"""
pyBADA
Common configuration module
Developped @EUROCONTROL (EIH)
2024
"""

import os
import importlib.resources


def getVersionsList(badaFamily):
    # list file and directories
    path = getBadaFamilyPath(badaFamily)
    items = os.listdir(path)

    versionsList = []
    for item in items:
        if os.path.isdir(os.path.join(path, item)):
            versionsList.append(item)

    return versionsList


def getAircraftList(badaFamily, badaVersion):
    path = getBadaVersionPath(badaFamily, badaVersion)

    if not os.path.exists(path):
        return []
    else:
        items = os.listdir(path)

    # check if I have BADA3 xml or standard ACSII files
    xml = True
    for item in items:
        if "OPF" in item:
            xml = False
            break

    if (
        badaFamily == "BADA4"
        or badaFamily == "BADAH"
        or badaFamily == "BADAE"
        or xml
    ):
        aircraftList = []
        for item in items:
            if os.path.isdir(os.path.join(path, item)):
                aircraftList.append(item)

    elif badaFamily == "BADA3":
        aircraftList = []
        for item in items:
            if len(item.split(".")) == 2:
                if item.split(".")[0].rstrip("_") not in aircraftList and (
                    item.split(".")[1] == "PTD"
                    or item.split(".")[1] == "PTF"
                    or item.split(".")[1] == "OPF"
                    or item.split(".")[1] == "APF"
                ):
                    aircraftList.append(item.split(".")[0].rstrip("_"))

    return aircraftList


def getBadaFamilyPath(badaFamily):
    path = os.path.join(getAircraftPath(), badaFamily)
    return path


def getBadaVersionPath(badaFamily, badaVersion):
    path = os.path.join(getAircraftPath(), badaFamily, badaVersion)
    return path


def getAircraftPath():
    package_name = "pyBADA"
    resource_name = "aircraft"

    # Get the path to the 'aircraft' resource directory
    with importlib.resources.as_file(
        importlib.resources.files(package_name) / resource_name
    ) as resource_path:
        return str(resource_path)


def getDataPath():
    package_name = "pyBADA"
    resource_name = "data"

    # Get the path to the 'data' resource directory
    with importlib.resources.as_file(
        importlib.resources.files(package_name) / resource_name
    ) as resource_path:
        return str(resource_path)
