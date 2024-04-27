"""
pytest for field/teleconnection/index_EU.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import TEST_DATA_PATH

z500_data = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_z500_mon.nc")))["hgt"]


def test_calc_index_EU_Wallace_Gutzler_1981():
    result_data = ecl.field.teleconnection.calc_index_EU_Wallace_Gutzler_1981(
        z500_data
    ).data[:20]
    refer_data = np.array(
        [
            1.9321834,
            0.5821855,
            -0.08236959,
            -0.01967423,
            0.08009898,
            1.0334079,
            -0.24247321,
            1.2930799,
            0.34333292,
            2.167064,
            2.2324498,
            2.876514,
            -0.16640146,
            -1.9561028,
            -1.3748313,
            1.3194227,
            -0.12642841,
            1.2132866,
            0.04400028,
            -0.10597672,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()
