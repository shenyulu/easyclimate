"""
pytest for field.air_sea_interaction.index_atlantic_nino.py
"""

import pytest

import easyclimate as ecl
import xarray as xr
import numpy as np
from pathlib import Path
from .const_define import TEST_DATA_PATH

data_sst = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_oisst_data.nc")))["sst"]


def test_calc_index_ATL3_1():
    result_data = ecl.field.air_sea_interaction.calc_index_ATL3(data_sst).data[:20]
    refer_data = np.array(
        [
            0.3006445,
            0.17614102,
            0.23786576,
            0.02363809,
            -0.38653436,
            -0.4954776,
            -0.8741972,
            -1.0537297,
            -1.0030487,
            -0.6035976,
            -0.4867375,
            -0.66027457,
            -0.6431055,
            -0.45385924,
            -0.21088402,
            -0.62552875,
            -0.69945097,
            -1.0854776,
            -1.2366971,
            -0.7416465,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_index_ATL3_2():
    result_data = ecl.field.air_sea_interaction.calc_index_ATL3(
        data_sst, normalized=True
    ).data[:20]
    refer_data = np.array(
        [
            0.66056675,
            0.3870116,
            0.5226312,
            0.05193688,
            -0.8492813,
            -1.088648,
            -1.9207588,
            -2.315222,
            -2.2038674,
            -1.3262058,
            -1.0694445,
            -1.4507347,
            -1.4130114,
            -0.99720544,
            -0.46334782,
            -1.3743923,
            -1.536812,
            -2.3849776,
            -2.7172325,
            -1.6295227,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()
