"""
pytest for field.air_sea_interaction.index_enso.py
"""

import pytest

import easyclimate as ecl
import xarray as xr
import numpy as np
from pathlib import Path
from .const_define import TEST_DATA_PATH

data_sst = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_oisst_data.nc")))["sst"]


def test_calc_index_nino1and2_1():
    result_data = ecl.field.air_sea_interaction.calc_index_nino1and2(data_sst).data[:20]
    refer_data = np.array(
        [
            0.02753249,
            -0.18468605,
            -0.51840377,
            -1.190344,
            -0.86654085,
            -0.2813302,
            0.05576498,
            0.7211977,
            0.996652,
            1.4282483,
            2.0311863,
            2.8239908,
            3.1566236,
            2.8307683,
            1.9838692,
            2.077838,
            2.9220955,
            3.745942,
            4.3321285,
            4.054379,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_index_nino1and2_2():
    result_data = ecl.field.air_sea_interaction.calc_index_nino1and2(
        data_sst, normalized=True
    ).data[:20]
    refer_data = np.array(
        [
            0.02533327,
            -0.16993383,
            -0.47699508,
            -1.0952625,
            -0.79732394,
            -0.25885832,
            0.05131062,
            0.66359043,
            0.91704214,
            1.3141637,
            1.8689407,
            2.5984182,
            2.9044812,
            2.6046543,
            1.8254031,
            1.911866,
            2.6886866,
            3.4467266,
            3.98609,
            3.7305262,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_index_nino3_1():
    result_data = ecl.field.air_sea_interaction.calc_index_nino3(data_sst).data[:20]
    refer_data = np.array(
        [
            -0.08004957,
            0.29165897,
            0.01038647,
            -0.17457151,
            0.1112155,
            0.59582627,
            0.858595,
            0.67968374,
            1.0168155,
            1.508404,
            1.9924703,
            2.39,
            2.9911504,
            3.0670993,
            2.3803065,
            1.7990283,
            1.4186555,
            1.6728662,
            1.5577952,
            1.0040839,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_index_nino3_2():
    result_data = ecl.field.air_sea_interaction.calc_index_nino3(
        data_sst, normalized=True
    ).data[:20]
    refer_data = np.array(
        [
            -0.09108563,
            0.3318686,
            0.01181841,
            -0.19863886,
            0.12654825,
            0.67797,
            0.97696537,
            0.7733885,
            1.156999,
            1.7163604,
            2.2671626,
            2.719498,
            3.4035258,
            3.4899457,
            2.708468,
            2.0470517,
            1.6142387,
            1.9034963,
            1.772561,
            1.1425121,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_index_nino34_1():
    result_data = ecl.field.air_sea_interaction.calc_index_nino34(data_sst).data[:20]
    refer_data = np.array(
        [
            np.nan,
            np.nan,
            0.02695771,
            0.20392498,
            0.35752553,
            0.5115632,
            0.71459955,
            0.92648214,
            1.1543809,
            1.3398094,
            1.66371,
            1.9834204,
            2.170162,
            2.1392753,
            1.9423673,
            1.646401,
            1.2343637,
            0.7817538,
            0.41924733,
            0.191187,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data, equal_nan=True).all()


def test_calc_index_nino34_2():
    result_data = ecl.field.air_sea_interaction.calc_index_nino34(
        data_sst, normalized=True
    ).data[:20]
    refer_data = np.array(
        [
            np.nan,
            np.nan,
            0.03280235,
            0.2481375,
            0.43503982,
            0.622474,
            0.86953026,
            1.1273507,
            1.4046596,
            1.6302905,
            2.0244153,
            2.4134412,
            2.6406698,
            2.6030867,
            2.3634875,
            2.0033536,
            1.5019833,
            0.95124406,
            0.5101434,
            0.23263781,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data, equal_nan=True).all()


def test_calc_index_OMI_1():
    result_data = ecl.field.air_sea_interaction.calc_index_OMI(data_sst).data[:20]
    refer_data = np.array(
        [
            np.nan,
            -0.01099452,
            0.06220957,
            0.05326815,
            0.26699948,
            0.60919774,
            0.7966813,
            0.9133341,
            1.0015292,
            1.3383623,
            1.651214,
            2.0253942,
            2.2950668,
            2.38545,
            2.1309714,
            1.582792,
            1.1524597,
            0.7868804,
            0.4510941,
            0.08571508,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data, equal_nan=True).all()


def test_calc_index_OMI_2():
    result_data = ecl.field.air_sea_interaction.calc_index_OMI(
        data_sst, normalized=True
    ).data[:20]
    refer_data = np.array(
        [
            np.nan,
            -0.01296635,
            0.07336665,
            0.06282161,
            0.31488496,
            0.7184553,
            0.93956345,
            1.0771375,
            1.1811502,
            1.5783932,
            1.9473537,
            2.3886418,
            2.7066793,
            2.8132725,
            2.513154,
            1.8666606,
            1.3591496,
            0.92800474,
            0.5319963,
            0.10108779,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data, equal_nan=True).all()


def test_calc_index_nino4_1():
    result_data = ecl.field.air_sea_interaction.calc_index_nino4(data_sst).data[:20]
    refer_data = np.array(
        [
            -0.12145369,
            0.04968018,
            0.0811892,
            0.1305759,
            0.31786987,
            0.6067503,
            0.80823463,
            0.5229595,
            0.28205344,
            0.37909153,
            0.60953766,
            0.49988833,
            0.5470226,
            0.6256802,
            0.5576654,
            0.4084806,
            0.17920312,
            0.17779803,
            -0.04043187,
            -0.19818328,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_index_nino4_2():
    result_data = ecl.field.air_sea_interaction.calc_index_nino4(
        data_sst, normalized=True
    ).data[:20]
    refer_data = np.array(
        [
            -0.18653181,
            0.07630015,
            0.12469253,
            0.20054193,
            0.488193,
            0.93186325,
            1.2413082,
            0.80317503,
            0.43318516,
            0.5822188,
            0.9361441,
            0.7677418,
            0.8401318,
            0.9609363,
            0.8564774,
            0.6273554,
            0.27522492,
            0.27306694,
            -0.06209634,
            -0.30437514,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()
