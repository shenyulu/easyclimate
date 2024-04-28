"""
pytest for field.air_sea_interaction.index_iod.py
"""

import pytest

import easyclimate as ecl
import xarray as xr
import numpy as np
from pathlib import Path
from .const_define import TEST_DATA_PATH

data_sst = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_oisst_data.nc")))["sst"]


def test_calc_index_IOBM_1point_1():
    result_data = ecl.field.air_sea_interaction.calc_index_IOBM_1point(data_sst).data[
        :20
    ]
    refer_data = np.array(
        [
            -0.3269011,
            -0.37761196,
            -0.33824706,
            -0.4004367,
            -0.28636456,
            -0.09277426,
            -0.13330695,
            -0.148778,
            -0.04841796,
            -0.04041302,
            0.04719116,
            0.13628912,
            0.20324852,
            0.09735067,
            0.23684023,
            0.16799226,
            0.07154068,
            -0.00504362,
            0.23734145,
            0.1825187,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_index_IOBM_1point_2():
    result_data = ecl.field.air_sea_interaction.calc_index_IOBM_1point(
        data_sst, normalized=True
    ).data[:20]
    refer_data = np.array(
        [
            -1.2821093,
            -1.4809978,
            -1.3266082,
            -1.5705166,
            -1.1231246,
            -0.36386153,
            -0.52283114,
            -0.5835088,
            -0.1898957,
            -0.15850025,
            0.18508418,
            0.53452724,
            0.7971426,
            0.38181028,
            0.9288897,
            0.6588673,
            0.28058326,
            -0.01978113,
            0.93085545,
            0.7158401,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_index_IOBM_EOF1_1():
    result_data = ecl.field.air_sea_interaction.calc_index_IOBM_EOF1(data_sst).data[:20]
    refer_data = np.array(
        [
            -1.23310228,
            -1.40070626,
            -1.22917518,
            -1.42958303,
            -0.94668675,
            -0.20262406,
            -0.35992287,
            -0.45552472,
            -0.10674089,
            -0.05047559,
            0.29850964,
            0.60768597,
            0.78026958,
            0.23512296,
            0.84072217,
            0.50297767,
            0.15172536,
            -0.06115037,
            0.92514055,
            0.76630208,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_index_IOBM_EOF1_2():
    result_data = ecl.field.air_sea_interaction.calc_index_IOBM_EOF1(
        data_sst, normalized=False
    ).data[:20]
    refer_data = np.array(
        [
            -0.05547989,
            -0.06302075,
            -0.05530321,
            -0.06431998,
            -0.04259345,
            -0.00911649,
            -0.0161937,
            -0.02049503,
            -0.0048025,
            -0.002271,
            0.01343058,
            0.02734108,
            0.03510599,
            0.01057868,
            0.03782588,
            0.02263004,
            0.00682645,
            -0.00275129,
            0.04162404,
            0.03447756,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()
