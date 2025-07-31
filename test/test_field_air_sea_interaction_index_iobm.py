"""
pytest for field.air_sea_interaction.index_iod.py
"""

import pytest

import easyclimate as ecl
import xarray as xr
import numpy as np
from pathlib import Path
from .const_define import DOCS_DATA_PATH

data_sst = xr.open_dataset(str(Path(DOCS_DATA_PATH, "test_input_oisst_data.nc")))["sst"]


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
    assert np.isclose(result_data, refer_data, atol=0.01).all()


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
    assert np.isclose(result_data, refer_data, atol=0.01).all()


def test_calc_index_IOBM_EOF1_1():
    result_data = ecl.field.air_sea_interaction.calc_index_IOBM_EOF1(
        data_sst, normalized=False
    ).data[:20]
    refer_data = np.array(
        [
            -6.35371888,
            -7.217317,
            -6.33348235,
            -7.36610971,
            -4.87792433,
            -1.04404388,
            -1.85454885,
            -2.34715087,
            -0.54999718,
            -0.26008055,
            1.53811037,
            3.13117964,
            4.02043694,
            1.21150212,
            4.33192829,
            2.59165628,
            0.78178441,
            -0.3150826,
            4.76690465,
            3.94846789,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data, atol=0.01).all()


def test_calc_index_IOBM_EOF1_2():
    result_data = ecl.field.air_sea_interaction.calc_index_IOBM_EOF1(
        data_sst, normalized=True
    ).data[:20]
    refer_data = np.array(
        [
            -0.05547991,
            -0.06302075,
            -0.05530321,
            -0.06431999,
            -0.04259345,
            -0.00911647,
            -0.0161937,
            -0.02049504,
            -0.00480251,
            -0.00227099,
            0.0134306,
            0.02734108,
            0.03510597,
            0.01057869,
            0.03782588,
            0.02263003,
            0.00682645,
            -0.00275126,
            0.04162404,
            0.03447755,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data, atol=0.01).all()
