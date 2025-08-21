"""
pytest for field/teleconnection/index_SRP.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import TEST_DATA_PATH

v200_data = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_v200_mon.nc")))["vwnd"]


def test_calc_index_SRP_EOF1_Yasui_Watanabe_2010():
    result_data = ecl.field.teleconnection.calc_index_SRP_EOF1_Yasui_Watanabe_2010(
        v200_data
    ).data[:20]
    refer_data = np.array(
        [
            -0.02578742,
            0.03841788,
            -0.03999095,
            -0.06775311,
            -0.03113178,
            0.01515116,
            0.03510514,
            0.02475167,
            0.00957336,
            0.07978619,
            -0.06067549,
            0.06182475,
            0.03466507,
            -0.06414709,
            -0.02530923,
            0.02098187,
            -0.00306072,
            -0.03611952,
            -0.00469981,
            -0.02526983,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data, atol=0.01).all()


def test_calc_index_SRP_EOF1_Kosaka_2009():
    result_data = ecl.field.teleconnection.calc_index_SRP_EOF1_Kosaka_2009(
        v200_data
    ).data[:20]
    refer_data = np.array(
        [
            0.03537321,
            0.04315557,
            -0.03528476,
            -0.08501737,
            -0.0489338,
            -0.00131242,
            0.0199855,
            0.05098543,
            0.01787614,
            0.07411663,
            0.01407047,
            0.09179481,
            0.04242809,
            -0.09060911,
            -0.03505788,
            0.03431226,
            0.00939577,
            -0.01901042,
            0.02938179,
            -0.02199294,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data, atol=0.01).all()


def test_calc_index_SRP_EOF1_Chen_Huang_2012():
    result_data = ecl.field.teleconnection.calc_index_SRP_EOF1_Chen_Huang_2012(
        v200_data
    ).data[:20]
    refer_data = np.array(
        [
            0.05603457,
            0.02980342,
            -0.03444552,
            -0.08626485,
            -0.06271177,
            -0.00129038,
            0.02879648,
            0.06670321,
            0.0231079,
            0.0784103,
            0.0303948,
            0.09838655,
            0.0040923,
            -0.07906909,
            -0.02068619,
            0.03899505,
            0.00151308,
            -0.01865891,
            0.02655686,
            -0.02087915,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data, atol=0.01).all()


def test_calc_index_SRP_EOF1_Sato_Takahashi_2006():
    result_data = ecl.field.teleconnection.calc_index_SRP_EOF1_Sato_Takahashi_2006(
        v200_data
    ).data[:20]
    refer_data = np.array(
        [
            0.02072678,
            0.00800326,
            0.01399657,
            0.03224734,
            -0.02559619,
            -0.03504753,
            0.00114697,
            -0.09218121,
            -0.01913933,
            0.01030175,
            -0.00055497,
            -0.02262819,
            -0.05380599,
            0.02889996,
            0.06257547,
            0.00733362,
            -0.03627269,
            -0.06349816,
            0.01302304,
            0.03247051,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data, atol=0.01).all()


def test_calc_index_SRP_1point_Lu_2002():
    result_data = ecl.field.teleconnection.calc_index_SRP_1point_Lu_2002(
        v200_data
    ).data[:20]
    refer_data = np.array(
        [
            -0.08426504,
            -0.55714005,
            -0.32565936,
            0.1075617,
            -2.3177822,
            -0.6324829,
            1.123698,
            0.71030414,
            0.39257553,
            -0.15486863,
            0.28079462,
            0.5026882,
            -1.4582796,
            1.0875138,
            1.089992,
            0.1700169,
            -0.65231,
            -0.21611482,
            -0.43768218,
            1.0225801,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data, atol=0.01).all()
