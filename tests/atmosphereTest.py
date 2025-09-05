"""
Atmosphere Module Test
======================

Test script for atmosphere module functions: theta, delta, sigma, and aSound
"""

import numpy as np
import pandas as pd
import xarray as xr

from pyBADA import atmosphere as atm
from pyBADA import constants as const


def run_theta_tests():
    """Test atm.theta with various input types"""
    tests = {
        "float": (1000.0, 2.0),
        "numpy": (
            np.array([0.0, 5000.0, const.h_11]),
            np.array([0.0, -5.0, 0.0]),
        ),
        "pandas.Series": (
            pd.Series([0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]),
            pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
        ),
        "pandas.DataFrame": (
            pd.DataFrame(
                {
                    "h": [0.0, 5000.0, const.h_11],
                    "DeltaTemp": [0.0, -5.0, 0.0],
                },
                index=["SL", "5km", "tropo"],
            )["h"],
            pd.DataFrame(
                {
                    "h": [0.0, 5000.0, const.h_11],
                    "DeltaTemp": [0.0, -5.0, 0.0],
                },
                index=["SL", "5km", "tropo"],
            )["DeltaTemp"],
        ),
        "xarray.DataArray": (
            xr.DataArray(
                [0.0, 5000.0, const.h_11],
                dims=("z",),
                coords={"z": ["SL", "5km", "tropo"]},
            ),
            xr.DataArray(
                [0.0, -5.0, 0.0],
                dims=("z",),
                coords={"z": ["SL", "5km", "tropo"]},
            ),
        ),
    }
    for name, (h_in, dT_in) in tests.items():
        th = atm.theta(h_in, dT_in)
        print(f"\n=== Theta Test: {name} ===")
        print(f"Input h ({type(h_in)}): {h_in}")
        print(f"Input ΔT ({type(dT_in)}): {dT_in}")
        print(f"Output θ ({type(th)}): {th}")


def run_delta_tests():
    """Test atm.delta with various input types"""
    tests = {
        "float": (1000.0, 2.0),
        "numpy": (
            np.array([0.0, 5000.0, const.h_11]),
            np.array([0.0, -5.0, 0.0]),
        ),
        "pandas.Series": (
            pd.Series([0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]),
            pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
        ),
        "pandas.DataFrame": (
            pd.DataFrame(
                {
                    "h": [0.0, 5000.0, const.h_11],
                    "DeltaTemp": [0.0, -5.0, 0.0],
                },
                index=["SL", "5km", "tropo"],
            )["h"],
            pd.DataFrame(
                {
                    "h": [0.0, 5000.0, const.h_11],
                    "DeltaTemp": [0.0, -5.0, 0.0],
                },
                index=["SL", "5km", "tropo"],
            )["DeltaTemp"],
        ),
        "xarray.DataArray": (
            xr.DataArray(
                [0.0, 5000.0, const.h_11],
                dims=("z",),
                coords={"z": ["SL", "5km", "tropo"]},
            ),
            xr.DataArray(
                [0.0, -5.0, 0.0],
                dims=("z",),
                coords={"z": ["SL", "5km", "tropo"]},
            ),
        ),
    }
    for name, (h_in, dT_in) in tests.items():
        dp = atm.delta(h_in, dT_in)
        print(f"\n=== Delta Test: {name} ===")
        print(f"Input h ({type(h_in)}): {h_in}")
        print(f"Input ΔT ({type(dT_in)}): {dT_in}")
        print(f"Output δ ({type(dp)}): {dp}")


def run_sigma_tests():
    """Test atm.sigma with two input modes"""
    tests = {
        "float": (1000.0, 2.0),
        "numpy": (
            np.array([0.0, 5000.0, const.h_11]),
            np.array([0.0, -5.0, 0.0]),
        ),
        "pandas.Series": (
            pd.Series([0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]),
            pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
        ),
        "pandas.DataFrame": (
            pd.DataFrame(
                {
                    "h": [0.0, 5000.0, const.h_11],
                    "DeltaTemp": [0.0, -5.0, 0.0],
                },
                index=["SL", "5km", "tropo"],
            )["h"],
            pd.DataFrame(
                {
                    "h": [0.0, 5000.0, const.h_11],
                    "DeltaTemp": [0.0, -5.0, 0.0],
                },
                index=["SL", "5km", "tropo"],
            )["DeltaTemp"],
        ),
        "xarray.DataArray": (
            xr.DataArray(
                [0.0, 5000.0, const.h_11],
                dims=("z",),
                coords={"z": ["SL", "5km", "tropo"]},
            ),
            xr.DataArray(
                [0.0, -5.0, 0.0],
                dims=("z",),
                coords={"z": ["SL", "5km", "tropo"]},
            ),
        ),
    }
    # Case A: h and DeltaTemp
    print("\n*** Sigma Tests: inputs (h, DeltaTemp) ***")
    for name, (h_in, dT_in) in tests.items():
        sigma_val = atm.sigma(h=h_in, DeltaTemp=dT_in)
        print(f"\n--- {name} ---")
        print(f"Input h ({type(h_in)}): {h_in}")
        print(f"Input ΔT ({type(dT_in)}): {dT_in}")
        print(f"Output σ ({type(sigma_val)}): {sigma_val}")

    # Case B: theta and delta
    theta_delta = {
        name: (atm.theta(h, dT), atm.delta(h, dT))
        for name, (h, dT) in tests.items()
    }
    print("\n*** Sigma Tests: inputs (theta, delta) ***")
    for name, (th, dp) in theta_delta.items():
        sigma_val = atm.sigma(theta=th, delta=dp)
        print(f"\n--- {name} ---")
        print(f"Input θ ({type(th)}): {th}")
        print(f"Input δ ({type(dp)}): {dp}")
        print(f"Output σ ({type(sigma_val)}): {sigma_val}")


def run_aSound_tests():
    """Test atm.aSound with theta inputs from previous tests"""
    tests = {
        "float": atm.theta(1000.0, 2.0),
        "numpy": atm.theta(
            np.array([0.0, 5000.0, const.h_11]), np.array([0.0, -5.0, 0.0])
        ),
        "pandas.Series": atm.theta(
            pd.Series([0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]),
            pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
        ),
        "pandas.DataFrame": atm.theta(
            pd.DataFrame(
                {
                    "h": [0.0, 5000.0, const.h_11],
                    "DeltaTemp": [0.0, -5.0, 0.0],
                },
                index=["SL", "5km", "tropo"],
            )["h"],
            pd.DataFrame(
                {
                    "h": [0.0, 5000.0, const.h_11],
                    "DeltaTemp": [0.0, -5.0, 0.0],
                },
                index=["SL", "5km", "tropo"],
            )["DeltaTemp"],
        ),
        "xarray.DataArray": atm.theta(
            xr.DataArray(
                [0.0, 5000.0, const.h_11],
                dims=("z",),
                coords={"z": ["SL", "5km", "tropo"]},
            ),
            xr.DataArray(
                [0.0, -5.0, 0.0],
                dims=("z",),
                coords={"z": ["SL", "5km", "tropo"]},
            ),
        ),
    }
    print("\n*** aSound Tests: input theta only ***")
    for name, th in tests.items():
        a_val = atm.aSound(th)
        print(f"\n--- {name} ---")
        print(f"Input θ ({type(th)}): {th}")
        print(f"Output aSound ({type(a_val)}): {a_val}")


def run_mach2Tas_tests():
    """Test atm.mach2Tas with Mach and theta inputs"""
    tests = {
        "float": (0.8, atm.theta(1000.0, 2.0)),
        "numpy": (
            np.array([0.5, 0.8, 1.0]),
            atm.theta(
                np.array([0.0, 5000.0, const.h_11]), np.array([0.0, -5.0, 0.0])
            ),
        ),
        "pandas.Series": (
            pd.Series([0.5, 0.8, 1.0], index=["low", "med", "high"]),
            atm.theta(
                pd.Series(
                    [0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]
                ),
                pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
            ),
        ),
        "pandas.DataFrame": (
            pd.DataFrame(
                {"Mach": [0.5, 0.8, 1.0]}, index=["low", "med", "high"]
            )["Mach"],
            atm.theta(
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["h"],
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["DeltaTemp"],
            ),
        ),
        "xarray.DataArray": (
            xr.DataArray(
                [0.5, 0.8, 1.0],
                dims=("m",),
                coords={"m": ["low", "med", "high"]},
            ),
            atm.theta(
                xr.DataArray(
                    [0.0, 5000.0, const.h_11],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
                xr.DataArray(
                    [0.0, -5.0, 0.0],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
            ),
        ),
    }
    print("\n*** mach2Tas Tests: inputs (Mach, theta) ***")
    for name, (M_in, th_in) in tests.items():
        tas = atm.mach2Tas(M_in, th_in)
        print(f"\n--- {name} ---")
        print(f"Input Mach ({type(M_in)}): {M_in}")
        print(f"Input θ ({type(th_in)}): {th_in}")
        print(f"Output TAS ({type(tas)}): {tas}")


def run_tas2Mach_tests():
    """Test atm.tas2Mach with TAS and theta inputs"""
    tests = {
        "float": (
            atm.mach2Tas(0.8, atm.theta(1000.0, 2.0)),
            atm.theta(1000.0, 2.0),
        ),
        "numpy": (
            atm.mach2Tas(
                np.array([0.5, 0.8, 1.0]),
                atm.theta(
                    np.array([0.0, 5000.0, const.h_11]),
                    np.array([0.0, -5.0, 0.0]),
                ),
            ),
            atm.theta(
                np.array([0.0, 5000.0, const.h_11]), np.array([0.0, -5.0, 0.0])
            ),
        ),
        "pandas.Series": (
            atm.mach2Tas(
                pd.Series([0.5, 0.8, 1.0], index=["low", "med", "high"]),
                atm.theta(
                    pd.Series(
                        [0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]
                    ),
                    pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
                ),
            ),
            atm.theta(
                pd.Series(
                    [0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]
                ),
                pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
            ),
        ),
        "pandas.DataFrame": (
            atm.mach2Tas(
                pd.DataFrame(
                    {"Mach": [0.5, 0.8, 1.0]}, index=["low", "med", "high"]
                )["Mach"],
                atm.theta(
                    pd.DataFrame(
                        {
                            "h": [0.0, 5000.0, const.h_11],
                            "DeltaTemp": [0.0, -5.0, 0.0],
                        },
                        index=["SL", "5km", "tropo"],
                    )["h"],
                    pd.DataFrame(
                        {
                            "h": [0.0, 5000.0, const.h_11],
                            "DeltaTemp": [0.0, -5.0, 0.0],
                        },
                        index=["SL", "5km", "tropo"],
                    )["DeltaTemp"],
                ),
            ),
            atm.theta(
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["h"],
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["DeltaTemp"],
            ),
        ),
        "xarray.DataArray": (
            atm.mach2Tas(
                xr.DataArray(
                    [0.5, 0.8, 1.0],
                    dims=("m",),
                    coords={"m": ["low", "med", "high"]},
                ),
                atm.theta(
                    xr.DataArray(
                        [0.0, 5000.0, const.h_11],
                        dims=("z",),
                        coords={"z": ["SL", "5km", "tropo"]},
                    ),
                    xr.DataArray(
                        [0.0, -5.0, 0.0],
                        dims=("z",),
                        coords={"z": ["SL", "5km", "tropo"]},
                    ),
                ),
            ),
            atm.theta(
                xr.DataArray(
                    [0.0, 5000.0, const.h_11],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
                xr.DataArray(
                    [0.0, -5.0, 0.0],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
            ),
        ),
    }
    print("\n*** tas2Mach Tests: inputs (TAS, theta) ***")
    for name, (tas_in, th_in) in tests.items():
        M_out = atm.tas2Mach(tas_in, th_in)
        print(f"\n--- {name} ---")
        print(f"Input TAS ({type(tas_in)}): {tas_in}")
        print(f"Input θ ({type(th_in)}): {th_in}")
        print(f"Output Mach ({type(M_out)}): {M_out}")


def run_tas2Cas_tests():
    """Test atm.tas2Cas with TAS, delta, and sigma inputs"""
    tests = {
        "float": (300.0, atm.delta(1000.0, 2.0), atm.sigma(1000.0, 2.0)),
        "numpy": (
            np.array([200.0, 300.0, 400.0]),
            atm.delta(
                np.array([0.0, 5000.0, const.h_11]), np.array([0.0, -5.0, 0.0])
            ),
            atm.sigma(
                np.array([0.0, 5000.0, const.h_11]), np.array([0.0, -5.0, 0.0])
            ),
        ),
        "pandas.Series": (
            pd.Series([200, 300, 400], index=["low", "med", "high"]),
            atm.delta(
                pd.Series(
                    [0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]
                ),
                pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
            ),
            atm.sigma(
                pd.Series(
                    [0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]
                ),
                pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
            ),
        ),
        "pandas.DataFrame": (
            pd.DataFrame(
                {"TAS": [200, 300, 400]}, index=["low", "med", "high"]
            )["TAS"],
            atm.delta(
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["h"],
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["DeltaTemp"],
            ),
            atm.sigma(
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["h"],
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["DeltaTemp"],
            ),
        ),
        "xarray.DataArray": (
            xr.DataArray(
                [200.0, 300.0, 400.0],
                dims=("t",),
                coords={"t": ["low", "med", "high"]},
            ),
            atm.delta(
                xr.DataArray(
                    [0.0, 5000.0, const.h_11],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
                xr.DataArray(
                    [0.0, -5.0, 0.0],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
            ),
            atm.sigma(
                xr.DataArray(
                    [0.0, 5000.0, const.h_11],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
                xr.DataArray(
                    [0.0, -5.0, 0.0],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
            ),
        ),
    }
    print("\n*** tas2Cas Tests: inputs (TAS, delta, sigma) ***")
    for name, (tas, dp, sg) in tests.items():
        cas = atm.tas2Cas(tas, dp, sg)
        print(f"\n--- {name} ---")
        print(f"Input TAS ({type(tas)}): {tas}")
        print(f"Input δ ({type(dp)}): {dp}")
        print(f"Input σ ({type(sg)}): {sg}")
        print(f"Output CAS ({type(cas)}): {cas}")


def run_cas2Tas_tests():
    """Test atm.cas2Tas with CAS, delta, and sigma inputs"""
    tests = {
        "float": (
            atm.tas2Cas(300.0, atm.delta(1000.0, 2.0), atm.sigma(1000.0, 2.0)),
            atm.delta(1000.0, 2.0),
            atm.sigma(1000.0, 2.0),
        ),
        "numpy": (
            atm.tas2Cas(
                np.array([200.0, 300.0, 400.0]),
                atm.delta(
                    np.array([0.0, 5000.0, const.h_11]),
                    np.array([0.0, -5.0, 0.0]),
                ),
                atm.sigma(
                    np.array([0.0, 5000.0, const.h_11]),
                    np.array([0.0, -5.0, 0.0]),
                ),
            ),
            atm.delta(
                np.array([0.0, 5000.0, const.h_11]), np.array([0.0, -5.0, 0.0])
            ),
            atm.sigma(
                np.array([0.0, 5000.0, const.h_11]), np.array([0.0, -5.0, 0.0])
            ),
        ),
        "pandas.Series": (
            atm.tas2Cas(
                pd.Series([200, 300, 400], index=["low", "med", "high"]),
                atm.delta(
                    pd.Series(
                        [0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]
                    ),
                    pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
                ),
                atm.sigma(
                    pd.Series(
                        [0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]
                    ),
                    pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
                ),
            ),
            atm.delta(
                pd.Series(
                    [0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]
                ),
                pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
            ),
            atm.sigma(
                pd.Series(
                    [0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]
                ),
                pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
            ),
        ),
        "pandas.DataFrame": (
            atm.tas2Cas(
                pd.DataFrame(
                    {"CAS": [200, 300, 400]}, index=["low", "med", "high"]
                )["CAS"],
                atm.delta(
                    pd.DataFrame(
                        {
                            "h": [0.0, 5000.0, const.h_11],
                            "DeltaTemp": [0.0, -5.0, 0.0],
                        },
                        index=["SL", "5km", "tropo"],
                    )["h"],
                    pd.DataFrame(
                        {
                            "h": [0.0, 5000.0, const.h_11],
                            "DeltaTemp": [0.0, -5.0, 0.0],
                        },
                        index=["SL", "5km", "tropo"],
                    )["DeltaTemp"],
                ),
                atm.sigma(
                    pd.DataFrame(
                        {
                            "h": [0.0, 5000.0, const.h_11],
                            "DeltaTemp": [0.0, -5.0, 0.0],
                        },
                        index=["SL", "5km", "tropo"],
                    )["h"],
                    pd.DataFrame(
                        {
                            "h": [0.0, 5000.0, const.h_11],
                            "DeltaTemp": [0.0, -5.0, 0.0],
                        },
                        index=["SL", "5km", "tropo"],
                    )["DeltaTemp"],
                ),
            ),
            atm.delta(
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["h"],
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["DeltaTemp"],
            ),
            atm.sigma(
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["h"],
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["DeltaTemp"],
            ),
        ),
        "xarray.DataArray": (
            atm.tas2Cas(
                xr.DataArray(
                    [200.0, 300.0, 400.0],
                    dims=("t",),
                    coords={"t": ["low", "med", "high"]},
                ),
                atm.delta(
                    xr.DataArray(
                        [0.0, 5000.0, const.h_11],
                        dims=("z",),
                        coords={"z": ["SL", "5km", "tropo"]},
                    ),
                    xr.DataArray(
                        [0.0, -5.0, 0.0],
                        dims=("z",),
                        coords={"z": ["SL", "5km", "tropo"]},
                    ),
                ),
                atm.sigma(
                    xr.DataArray(
                        [0.0, 5000.0, const.h_11],
                        dims=("z",),
                        coords={"z": ["SL", "5km", "tropo"]},
                    ),
                    xr.DataArray(
                        [0.0, -5.0, 0.0],
                        dims=("z",),
                        coords={"z": ["SL", "5km", "tropo"]},
                    ),
                ),
            ),
            atm.delta(
                xr.DataArray(
                    [0.0, 5000.0, const.h_11],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
                xr.DataArray(
                    [0.0, -5.0, 0.0],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
            ),
            atm.sigma(
                xr.DataArray(
                    [0.0, 5000.0, const.h_11],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
                xr.DataArray(
                    [0.0, -5.0, 0.0],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
            ),
        ),
    }
    print("\n*** cas2Tas Tests: inputs (CAS, delta, sigma) ***")
    for name, (cas, dp, sg) in tests.items():
        tas = atm.cas2Tas(cas, dp, sg)
        print(f"\n--- {name} ---")
        print(f"Input CAS ({type(cas)}): {cas}")
        print(f"Input δ ({type(dp)}): {dp}")
        print(f"Input σ ({type(sg)}): {sg}")
        print(f"Output TAS ({type(tas)}): {tas}")


def run_mach2Cas_tests():
    """Test atm.mach2Cas with Mach, theta, delta, sigma inputs"""
    tests = {
        "float": (
            0.8,
            atm.theta(5000.0, -5.0),
            atm.delta(5000.0, -5.0),
            atm.sigma(5000.0, -5.0),
        ),
        "numpy": (
            np.array([0.5, 0.8, 1.2]),
            atm.theta(
                np.array([0.0, 5000.0, const.h_11]), np.array([0.0, -5.0, 0.0])
            ),
            atm.delta(
                np.array([0.0, 5000.0, const.h_11]), np.array([0.0, -5.0, 0.0])
            ),
            atm.sigma(
                np.array([0.0, 5000.0, const.h_11]), np.array([0.0, -5.0, 0.0])
            ),
        ),
        "pandas.Series": (
            pd.Series([0.5, 0.8, 1.2], index=["M0.5", "M0.8", "M1.2"]),
            atm.theta(
                pd.Series(
                    [0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]
                ),
                pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
            ),
            atm.delta(
                pd.Series(
                    [0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]
                ),
                pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
            ),
            atm.sigma(
                pd.Series(
                    [0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]
                ),
                pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
            ),
        ),
        "pandas.DataFrame": (
            pd.DataFrame(
                {"Mach": [0.5, 0.8, 1.2]}, index=["M0.5", "M0.8", "M1.2"]
            )["Mach"],
            atm.theta(
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["h"],
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["DeltaTemp"],
            ),
            atm.delta(
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["h"],
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["DeltaTemp"],
            ),
            atm.sigma(
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["h"],
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["DeltaTemp"],
            ),
        ),
        "xarray.DataArray": (
            xr.DataArray(
                [0.5, 0.8, 1.2],
                dims=("m",),
                coords={"m": ["M0.5", "M0.8", "M1.2"]},
            ),
            atm.theta(
                xr.DataArray(
                    [0.0, 5000.0, const.h_11],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
                xr.DataArray(
                    [0.0, -5.0, 0.0],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
            ),
            atm.delta(
                xr.DataArray(
                    [0.0, 5000.0, const.h_11],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
                xr.DataArray(
                    [0.0, -5.0, 0.0],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
            ),
            atm.sigma(
                xr.DataArray(
                    [0.0, 5000.0, const.h_11],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
                xr.DataArray(
                    [0.0, -5.0, 0.0],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
            ),
        ),
    }
    print("\n*** mach2Cas Tests: inputs (Mach, theta, delta, sigma) ***")
    for name, (m_in, th_in, dp_in, sg_in) in tests.items():
        cas_out = atm.mach2Cas(m_in, th_in, dp_in, sg_in)
        print(f"\n--- {name} ---")
        print(f"Input Mach ({type(m_in)}): {m_in}")
        print(f"Input θ ({type(th_in)}): {th_in}")
        print(f"Input δ ({type(dp_in)}): {dp_in}")
        print(f"Input σ ({type(sg_in)}): {sg_in}")
        print(f"Output CAS ({type(cas_out)}): {cas_out}")


def run_cas2Mach_tests():
    """Test atm.cas2Mach with CAS, theta, delta, sigma inputs"""
    tests = {
        "float": (
            300.0,
            atm.theta(5000.0, -5.0),
            atm.delta(5000.0, -5.0),
            atm.sigma(5000.0, -5.0),
        ),
        "numpy": (
            np.array([200.0, 240.0, 280.0]),
            atm.theta(
                np.array([0.0, 5000.0, const.h_11]), np.array([0.0, -5.0, 0.0])
            ),
            atm.delta(
                np.array([0.0, 5000.0, const.h_11]), np.array([0.0, -5.0, 0.0])
            ),
            atm.sigma(
                np.array([0.0, 5000.0, const.h_11]), np.array([0.0, -5.0, 0.0])
            ),
        ),
        "pandas.Series": (
            pd.Series(
                [200.0, 240.0, 280.0], index=["CAS200", "CAS240", "CAS280"]
            ),
            atm.theta(
                pd.Series(
                    [0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]
                ),
                pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
            ),
            atm.delta(
                pd.Series(
                    [0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]
                ),
                pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
            ),
            atm.sigma(
                pd.Series(
                    [0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]
                ),
                pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
            ),
        ),
        "pandas.DataFrame": (
            pd.DataFrame(
                {"CAS": [200.0, 240.0, 280.0]},
                index=["CAS200", "CAS240", "CAS280"],
            )["CAS"],
            atm.theta(
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["h"],
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["DeltaTemp"],
            ),
            atm.delta(
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["h"],
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["DeltaTemp"],
            ),
            atm.sigma(
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["h"],
                pd.DataFrame(
                    {
                        "h": [0.0, 5000.0, const.h_11],
                        "DeltaTemp": [0.0, -5.0, 0.0],
                    },
                    index=["SL", "5km", "tropo"],
                )["DeltaTemp"],
            ),
        ),
        "xarray.DataArray": (
            xr.DataArray(
                [200.0, 240.0, 280.0],
                dims=("cas",),
                coords={"cas": ["CAS200", "CAS240", "CAS280"]},
            ),
            atm.theta(
                xr.DataArray(
                    [0.0, 5000.0, const.h_11],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
                xr.DataArray(
                    [0.0, -5.0, 0.0],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
            ),
            atm.delta(
                xr.DataArray(
                    [0.0, 5000.0, const.h_11],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
                xr.DataArray(
                    [0.0, -5.0, 0.0],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
            ),
            atm.sigma(
                xr.DataArray(
                    [0.0, 5000.0, const.h_11],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
                xr.DataArray(
                    [0.0, -5.0, 0.0],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
            ),
        ),
    }
    print("\n*** cas2Mach Tests: inputs (CAS, theta, delta, sigma) ***")
    for name, (cas_in, th_in, dp_in, sg_in) in tests.items():
        m_out = atm.cas2Mach(cas_in, th_in, dp_in, sg_in)
        print(f"\n--- {name} ---")
        print(f"Input CAS ({type(cas_in)}): {cas_in}")
        print(f"Input θ ({type(th_in)}): {th_in}")
        print(f"Input δ ({type(dp_in)}): {dp_in}")
        print(f"Input σ ({type(sg_in)}): {sg_in}")
        print(f"Output Mach ({type(m_out)}): {m_out}")


def run_pressureAltitude_tests():
    """Test atm.pressureAltitude with pressure and default QNH (101325.0 Pa)"""
    tests = {
        "float": 101325.0,
        "numpy": np.array([90000.0, 101325.0, 110000.0]),
        "pandas.Series": pd.Series(
            [90000.0, 101325.0, 110000.0], index=["P90k", "P101k", "P110k"]
        ),
        "xarray.DataArray": xr.DataArray(
            [90000.0, 101325.0, 110000.0],
            dims=("p",),
            coords={"p": ["P90k", "P101k", "P110k"]},
        ),
    }
    print("\n*** pressureAltitude Tests: inputs (pressure, QNH=101325 Pa) ***")
    for name, p_in in tests.items():
        pa_out = atm.pressureAltitude(p_in, QNH=101325.0)
        print(f"\n--- {name} ---")
        print(f"Input pressure ({type(p_in)}): {p_in}")
        print(f"Output pressure altitude ({type(pa_out)}): {pa_out}")


def run_ISATemperatureDeviation_tests():
    """Test atm.ISATemperatureDeviation with temperature and pressure altitude inputs"""
    tests = {
        "float": (298.15, 0.0),  # Temperature (K), altitude (m)
        "numpy": (
            np.array([280.0, 290.0, 300.0]),
            np.array([0.0, 5000.0, 11000.0]),
        ),
        "pandas.Series": (
            pd.Series([280.0, 290.0, 300.0], index=["T280", "T290", "T300"]),
            pd.Series([0.0, 5000.0, 11000.0], index=["SL", "5km", "11km"]),
        ),
        "xarray.DataArray": (
            xr.DataArray(
                [280.0, 290.0, 300.0],
                dims=("t",),
                coords={"t": ["T280", "T290", "T300"]},
            ),
            xr.DataArray(
                [0.0, 5000.0, 11000.0],
                dims=("h",),
                coords={"h": ["SL", "5km", "11km"]},
            ),
        ),
    }
    print(
        "\n*** ISATemperatureDeviation Tests: inputs (temperature, altitude) ***"
    )
    for name, (T_in, h_in) in tests.items():
        dev_out = atm.ISATemperatureDeviation(T_in, h_in)
        print(f"\n--- {name} ---")
        print(f"Input temperature ({type(T_in)}): {T_in}")
        print(f"Input altitude ({type(h_in)}): {h_in}")
        print(f"Output deviation ({type(dev_out)}): {dev_out}")


def run_crossOver_tests():
    """Test atm.crossOver with CAS and Mach inputs"""
    tests = {
        "float": (150.0, 0.8),
        "numpy": (np.array([200.0, 250.0, 300.0]), np.array([0.5, 0.8, 1.0])),
        "pandas.Series": (
            pd.Series(
                [200.0, 250.0, 300.0], index=["CAS200", "CAS250", "CAS300"]
            ),
            pd.Series([0.5, 0.8, 1.0], index=["M0.5", "M0.8", "M1.0"]),
        ),
        "xarray.DataArray": (
            xr.DataArray(
                [200.0, 250.0, 300.0],
                dims=("cas",),
                coords={"cas": ["C200", "C250", "C300"]},
            ),
            xr.DataArray(
                [0.5, 0.8, 1.0],
                dims=("m",),
                coords={"m": ["M0.5", "M0.8", "M1.0"]},
            ),
        ),
    }
    print("\n*** crossOver Tests: inputs (CAS, Mach) ***")
    for name, (cas_in, m_in) in tests.items():
        cross_out = atm.crossOver(cas_in, m_in)
        print(f"\n--- {name} ---")
        print(f"Input CAS ({type(cas_in)}): {cas_in}")
        print(f"Input Mach ({type(m_in)}): {m_in}")
        print(f"Output crossover ({type(cross_out)}): {cross_out}")


def run_atmosphereProperties_tests():
    """Test atm.atmosphereProperties with altitude and DeltaTemp inputs"""
    tests = {
        "float": (1000.0, 10.0),
        "numpy": (
            np.array([0.0, 5000.0, 11000.0]),
            np.array([0.0, -5.0, 10.0]),
        ),
        "pandas.Series": (
            pd.Series([0.0, 5000.0, 11000.0], index=["SL", "5km", "11km"]),
            pd.Series([0.0, -5.0, 10.0], index=["SL", "5km", "11km"]),
        ),
        "xarray.DataArray": (
            xr.DataArray(
                [0.0, 5000.0, 11000.0],
                dims=("h",),
                coords={"h": ["SL", "5km", "11km"]},
            ),
            xr.DataArray(
                [0.0, -5.0, 10.0],
                dims=("dT",),
                coords={"dT": ["SL", "5km", "11km"]},
            ),
        ),
    }
    print("\n*** atmosphereProperties Tests: inputs (h, DeltaTemp) ***")
    for name, (h_in, dT_in) in tests.items():
        props = atm.atmosphereProperties(h_in, dT_in)
        print(f"\n--- {name} ---")
        print(f"Input altitude ({type(h_in)}): {h_in}")
        print(f"Input DeltaTemp ({type(dT_in)}): {dT_in}")
        print(f"Output properties ({type(props)}): {props}")


def run_convertSpeed_tests():
    """Test atm.convertSpeed with speed, speedType, theta, delta, sigma inputs"""
    # base theta, delta, sigma for scalar tests
    base_th = atm.theta(5000.0, -5.0)
    base_dp = atm.delta(5000.0, -5.0)
    base_sg = atm.sigma(5000.0, -5.0)

    tests = {
        "float": (250.0, "TAS", base_th, base_dp, base_sg),
        "numpy": (
            np.array([200.0, 250.0, 300.0]),
            "CAS",
            atm.theta(
                np.array([0.0, 5000.0, const.h_11]), np.array([0.0, -5.0, 0.0])
            ),
            atm.delta(
                np.array([0.0, 5000.0, const.h_11]), np.array([0.0, -5.0, 0.0])
            ),
            atm.sigma(
                np.array([0.0, 5000.0, const.h_11]), np.array([0.0, -5.0, 0.0])
            ),
        ),
        "pandas.Series": (
            pd.Series([200.0, 250.0, 300.0], index=["V200", "V250", "V300"]),
            "MACH",
            atm.theta(
                pd.Series(
                    [0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]
                ),
                pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
            ),
            atm.delta(
                pd.Series(
                    [0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]
                ),
                pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
            ),
            atm.sigma(
                pd.Series(
                    [0.0, 5000.0, const.h_11], index=["SL", "5km", "tropo"]
                ),
                pd.Series([0.0, -5.0, 0.0], index=["SL", "5km", "tropo"]),
            ),
        ),
        "xarray.DataArray": (
            xr.DataArray(
                [200.0, 250.0, 300.0],
                dims=("v",),
                coords={"v": ["V200", "V250", "V300"]},
            ),
            "MACH",
            atm.theta(
                xr.DataArray(
                    [0.0, 5000.0, const.h_11],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
                xr.DataArray(
                    [0.0, -5.0, 0.0],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
            ),
            atm.delta(
                xr.DataArray(
                    [0.0, 5000.0, const.h_11],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
                xr.DataArray(
                    [0.0, -5.0, 0.0],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
            ),
            atm.sigma(
                xr.DataArray(
                    [0.0, 5000.0, const.h_11],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
                xr.DataArray(
                    [0.0, -5.0, 0.0],
                    dims=("z",),
                    coords={"z": ["SL", "5km", "tropo"]},
                ),
            ),
        ),
    }
    print(
        "\n*** convertSpeed Tests: inputs (v, speedType, theta, delta, sigma) ***"
    )
    for name, (v_in, stype, th_in, dp_in, sg_in) in tests.items():
        conv_out = atm.convertSpeed(v_in, stype, th_in, dp_in, sg_in)
        print(f"\n--- {name} ---")
        print(f"Input speed ({type(v_in)}): {v_in}")
        print(f"Speed type: {stype}")
        print(f"Input θ ({type(th_in)}): {th_in}")
        print(f"Input δ ({type(dp_in)}): {dp_in}")
        print(f"Input σ ({type(sg_in)}): {sg_in}")
        print(f"Output speed ({type(conv_out)}): {conv_out}")


if __name__ == "__main__":
    run_theta_tests()
    run_delta_tests()
    run_sigma_tests()
    run_aSound_tests()
    run_mach2Tas_tests()
    run_tas2Mach_tests()
    run_tas2Cas_tests()
    run_cas2Tas_tests()
    run_mach2Cas_tests()
    run_cas2Mach_tests()
    run_pressureAltitude_tests()
    run_ISATemperatureDeviation_tests()
    run_crossOver_tests()
    run_atmosphereProperties_tests()
    run_convertSpeed_tests()
