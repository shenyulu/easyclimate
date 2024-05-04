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


def test_calc_index_IOD_Saji_1999_1():
    result_data = ecl.field.air_sea_interaction.calc_index_IOD_Saji_1999(data_sst).data[
        :20
    ]
    refer_data = np.array(
        [
            -0.21140042,
            0.1964276,
            0.26546946,
            0.3319999,
            0.28186905,
            0.51109624,
            0.7092762,
            0.5106905,
            0.21804449,
            0.23889592,
            0.5342278,
            0.45708993,
            0.06397119,
            -0.2634442,
            -0.4175273,
            -0.71197975,
            -0.72233343,
            -0.07891706,
            0.5505295,
            0.69957227,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_index_IOD_Saji_1999_2():
    result_data = ecl.field.air_sea_interaction.calc_index_IOD_Saji_1999(
        data_sst, normalized=True
    ).data[:20]
    refer_data = np.array(
        [
            -0.45415312,
            0.42198688,
            0.57031006,
            0.7132379,
            0.6055414,
            1.0979918,
            1.5237434,
            1.0971203,
            0.46842661,
            0.5132219,
            1.1476855,
            0.9819697,
            0.13742979,
            -0.5659591,
            -0.896977,
            -1.5295514,
            -1.5517943,
            -0.16953811,
            1.1827066,
            1.5028963,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()
