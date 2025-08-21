"""
pytest for field/teleconnection/index_PNA.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import TEST_DATA_PATH

z500_data = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_z500_mon.nc")))["hgt"]


def test_calc_index_PNA_modified_pointwise():
    result_data = ecl.field.teleconnection.calc_index_PNA_modified_pointwise(
        z500_data
    ).data[:20]
    refer_data = np.array(
        [
            3.2386305,
            0.5221999,
            2.2075872,
            -1.0701959,
            0.97871405,
            -0.11170971,
            -0.354215,
            0.74102503,
            0.5632406,
            -0.76577324,
            1.5145679,
            0.42980206,
            -2.4111109,
            -1.3430984,
            -1.9059343,
            -1.4022393,
            -0.11519514,
            0.32205728,
            0.17806418,
            -0.13098219,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data, atol=0.01).all()


def test_calc_index_PNA_Wallace_Gutzler_1981():
    result_data = ecl.field.teleconnection.calc_index_PNA_Wallace_Gutzler_1981(
        z500_data
    ).data[:20]
    refer_data = np.array(
        [
            2.976821,
            0.5509343,
            2.1246343,
            -1.1577514,
            1.2946671,
            -0.01380592,
            -0.16604027,
            0.9383815,
            0.3553485,
            -0.90358675,
            1.2770851,
            0.18722364,
            -2.1787858,
            -1.6771985,
            -2.0632799,
            -1.564345,
            0.03422692,
            0.41718328,
            0.07791588,
            -0.19194871,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data, atol=0.01).all()


def test_calc_index_PNA_NH_REOF():
    result_data = ecl.field.teleconnection.calc_index_PNA_NH_REOF(
        z500_data, solver="randomized", random_state=1
    ).data[:20]
    refer_data = np.array(
        [
            -0.16801542,
            -0.02420407,
            -0.09339072,
            -0.00363626,
            -0.01715254,
            -0.0204162,
            0.01833936,
            -0.02896684,
            -0.00893314,
            0.01331341,
            -0.09232175,
            -0.02148104,
            0.04944776,
            0.07689335,
            0.0658079,
            0.02772096,
            -0.00564644,
            -0.02970626,
            -0.03875269,
            -0.00865372,
        ]
    )
    assert np.isclose(result_data, refer_data, atol=0.01).all()
