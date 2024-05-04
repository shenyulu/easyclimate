"""
pytest for field/teleconnection/index_WP.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import TEST_DATA_PATH

z500_data = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_z500_mon.nc")))["hgt"]


def test_calc_index_WP_Wallace_Gutzler_1981():
    result_data = ecl.field.teleconnection.calc_index_WP_Wallace_Gutzler_1981(
        z500_data
    ).data[:20]
    refer_data = np.array(
        [
            -0.15485036,
            0.781371,
            -0.69214886,
            -0.75852716,
            0.04201345,
            1.0340837,
            0.05557111,
            0.56930345,
            0.22329152,
            0.5570654,
            -0.710228,
            -0.96328807,
            -0.59580755,
            1.9051006,
            -0.47878247,
            0.6639153,
            -0.27092388,
            -1.5263127,
            -0.21469297,
            -0.22726436,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()
