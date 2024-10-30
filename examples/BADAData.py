"""
BADA Data Retrieval
===================

Example of BADA parametes retrieval for specific aircraft
"""

from pyBADA.bada3 import Parser as Bada3Parser
from pyBADA.aircraft import Bada


# loading all the BADA data into a dataframe
allData = Bada3Parser.parseAll(badaVersion="DUMMY")

# retrieve specific data from the whole database, including synonyms
params = Bada.getBADAParameters(
    df=allData,
    acName=["A1", "P38", "AT45", "DA42", "B789", "J2H"],
    parameters=["VMO", "MMO", "MTOW", "engineType"],
)
print(params)
print("\n")

params = Bada.getBADAParameters(
    df=allData,
    acName=["B737"],
    parameters=["VMO", "MMO", "MTOW", "engineType"],
)
print(params)
print("\n")

params = Bada.getBADAParameters(
    df=allData, acName="DA42", parameters=["VMO", "MMO", "MTOW", "engineType"]
)
print(params)
print("\n")

params = Bada.getBADAParameters(df=allData, acName="DA42", parameters="VMO")
print(params)
print("\n")
