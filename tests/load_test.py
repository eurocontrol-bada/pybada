from pyBADA.bada3 import Parser as Bada3Parser
from pyBADA.aircraft import Bada


def test_data_load():
    # loading all the BADA data into a dataframe
    allData = Bada3Parser.parseAll(badaVersion="DUMMY")

    # retrieve specific data from the whole database, including synonyms
    params = Bada.getBADAParameters(
        df=allData,
        acName=["A1", "P38", "AT45", "DA42", "B789", "J2H"],
        parameters=["VMO", "MMO", "MTOW", "engineType"],
    )

    assert params.iloc[0]["VMO"] == 250.0
