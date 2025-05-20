"""
Common configuration module
"""

import os
import importlib.resources
import pandas as pd


@staticmethod
def list_subfolders(folderPath):
    """
    Lists all subfolders within a specified directory.

    :param folderPath: Path to the directory where subfolders are to be listed.
    :type folderPath: str
    :returns: A list of subfolder names within the specified directory.
    :rtype: list of str

    This function retrieves all entries in the given directory and filters out
    the ones that are not directories. Only the names of the subfolders are returned.
    """

    # List all entries in the directory
    entries = os.listdir(folderPath)

    # Filter out entries that are directories
    subfolders = [
        entry
        for entry in entries
        if os.path.isdir(os.path.join(folderPath, entry))
    ]

    return subfolders


@staticmethod
def safe_get(df, column_name, default_value=None):
    """Safely retrieves a column's value from a DataFrame, returning a default
    value if the column does not exist or if the value is NaN.

    :param df: DataFrame to retrieve the value from.
    :param column_name: Name of the column to retrieve.
    :param default_value: Value to return if the column does not exist.
        Default is None.
    :type df: pd.DataFrame
    :type column_name: str
    :type default_value: any
    :returns: The value from the specified column or the default value if the
        column is missing.
    :rtype: any
    """

    if column_name in df.columns:
        value = df[column_name].iloc[0]

        if isinstance(value, list):
            if pd.isna(value).all():
                return default_value
            else:
                return value
        else:
            if pd.isna(value):
                return default_value
            else:
                return value
    else:
        return default_value


def getVersionsList(badaFamily):
    """Retrieve a list of available BADA versions for a given BADA family.

    This function scans the directory corresponding to the specified BADA
    family and returns a list of all subdirectories (which represent different
    versions of BADA).

    :param badaFamily: The BADA family (e.g., BADA3, BADA4) for which versions
        are being retrieved.
    :type badaFamily: str.
    :returns: List of available BADA versions.
    :rtype: list of str.
    """

    # list file and directories
    path = getBadaFamilyPath(badaFamily)
    items = os.listdir(path)

    versionsList = []
    for item in items:
        if os.path.isdir(os.path.join(path, item)):
            versionsList.append(item)

    return versionsList


def getAircraftList(badaFamily, badaVersion):
    """Retrieve a list of available aircraft for a given BADA family and
    version.

    This function checks if the specified BADA family and version
    directory exists, and if so, determines whether the aircraft data is
    stored in XML format or as standard ASCII files (like OPF, APF, PTD,
    or PTF). It then returns a list of available aircraft.

    :param badaFamily: The BADA family (e.g., BADA3, BADA4) for which
        aircraft are being retrieved.
    :param badaVersion: The specific version of the BADA family (e.g.,
        3.10, 4.2).
    :type badaFamily: str.
    :type badaVersion: str.
    :returns: List of available aircraft names.
    :rtype: list of str.
    """

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
    """Get the full path to the specified BADA family directory.

    :param badaFamily: The BADA family (e.g., BADA3, BADA4) for which the path
        is required.
    :type badaFamily: str.
    :returns: The path to the BADA family directory.
    :rtype: str.
    """

    path = os.path.join(getAircraftPath(), badaFamily)
    return path


def getBadaVersionPath(badaFamily, badaVersion):
    """Get the full path to the specified BADA version directory.

    :param badaFamily: The BADA family (e.g., BADA3, BADA4) for which the path
        is required.
    :param badaVersion: The specific version of the BADA family.
    :type badaFamily: str.
    :type badaVersion: str.
    :returns: The path to the BADA version directory.
    :rtype: str.
    """

    path = os.path.join(getAircraftPath(), badaFamily, badaVersion)
    return path


def getAircraftPath():
    """Get the path to the 'aircraft' resource directory.

    This function locates the 'aircraft' directory within the pyBADA package
    and returns its absolute path.

    :returns: The absolute path to the 'aircraft' resource directory.
    :rtype: str.
    """

    package_name = "pyBADA"
    resource_name = "aircraft"

    # Get the path to the 'aircraft' resource directory
    with importlib.resources.as_file(
        importlib.resources.files(package_name) / resource_name
    ) as resource_path:
        return str(resource_path)


def getDataPath():
    """Get the path to the 'data' resource directory.

    This function locates the 'data' directory within the pyBADA package and
    returns its absolute path.

    :returns: The absolute path to the 'data' resource directory.
    :rtype: str.
    """

    package_name = "pyBADA"
    resource_name = "data"

    # Get the path to the 'data' resource directory
    with importlib.resources.as_file(
        importlib.resources.files(package_name) / resource_name
    ) as resource_path:
        return str(resource_path)
