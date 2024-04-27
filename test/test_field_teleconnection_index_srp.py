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
            -0.57720051,
            0.85990862,
            -0.89511854,
            -1.51651977,
            -0.69682355,
            0.33912893,
            0.7857593,
            0.55401739,
            0.21428069,
            1.78585668,
            -1.35810133,
            1.38382513,
            0.77590937,
            -1.43580617,
            -0.56649715,
            0.46963778,
            -0.06850823,
            -0.80846422,
            -0.10519598,
            -0.56561528,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_index_SRP_EOF1_Kosaka_2009():
    result_data = ecl.field.teleconnection.calc_index_SRP_EOF1_Kosaka_2009(
        v200_data
    ).data[:20]
    refer_data = np.array(
        [
            0.79175967,
            0.96595249,
            -0.78977978,
            -1.90294635,
            -1.0952867,
            -0.02937602,
            0.44733614,
            1.14120837,
            0.4001222,
            1.6589547,
            0.31493984,
            2.05464599,
            0.94966914,
            -2.02810631,
            -0.78470165,
            0.76801241,
            0.21030586,
            -0.42551081,
            0.65765345,
            -0.49226852,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_index_SRP_EOF1_Chen_Huang_2012():
    result_data = ecl.field.teleconnection.calc_index_SRP_EOF1_Chen_Huang_2012(
        v200_data
    ).data[:20]
    refer_data = np.array(
        [
            1.25422347,
            0.6670908,
            -0.77099506,
            -1.93086864,
            -1.40367947,
            -0.02888256,
            0.64455236,
            1.49301984,
            0.51722475,
            1.75506001,
            0.68032771,
            2.20218895,
            0.09159797,
            -1.76980581,
            -0.46301961,
            0.87282726,
            0.03386721,
            -0.41764294,
            0.59442292,
            -0.46733867,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_index_SRP_EOF1_Sato_Takahashi_2006():
    result_data = ecl.field.teleconnection.calc_index_SRP_EOF1_Sato_Takahashi_2006(
        v200_data
    ).data[:20]
    refer_data = np.array(
        [
            0.46392815,
            0.17913725,
            0.31328571,
            0.7217932,
            -0.57292016,
            -0.78446998,
            0.02567265,
            -2.06329465,
            -0.42839608,
            0.23058431,
            -0.01242201,
            -0.50648755,
            -1.20434099,
            0.6468687,
            1.40062858,
            0.1641486,
            -0.81189267,
            -1.42128125,
            0.29149499,
            0.72678846,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()


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
    assert np.isclose(result_data, refer_data).all()
