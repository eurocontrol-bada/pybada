"""Generic BADA Family Wrapper."""

from typing import Union

from pyBADA.bada3 import Bada3Aircraft
from pyBADA.bada4 import Bada4Aircraft
from pyBADA.badaH import BadaHAircraft


class BadaAircraft:
    """
    A wrapper to handle the instantiation of different BADA family models.
    """

    _FAMILY_MAP = {
        "BADA3": Bada3Aircraft,
        "BADA4": Bada4Aircraft,
        "BADAH": BadaHAircraft,
    }

    def __new__(
        cls, badaFamily, badaVersion, acName, **kwargs
    ) -> Union[Bada3Aircraft, Bada4Aircraft, BadaHAircraft]:
        target_class = cls._FAMILY_MAP.get(badaFamily)

        if not target_class:
            valid_keys = list(cls._FAMILY_MAP.keys())
            raise ValueError(
                f"Unknown BADA family: '{badaFamily}'. Valid options: {valid_keys}"
            )

        return target_class(badaVersion=badaVersion, acName=acName, **kwargs)

    def __init__(self, *args, **kwargs):
        pass
