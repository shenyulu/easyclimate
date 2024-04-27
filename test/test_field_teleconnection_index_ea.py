"""
pytest for field/teleconnection/index_EA.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import TEST_DATA_PATH

z500_data = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_z500_mon.nc")))["hgt"]


def test_calc_index_EA_Wallace_Gutzler_1981():
    result_data = ecl.field.teleconnection.calc_index_EA_Wallace_Gutzler_1981(
        z500_data
    ).data[:20]
    refer_data = np.array(
        [
            1.9084573,
            -0.8049945,
            -0.92284137,
            1.8036525,
            -1.7997904,
            0.92130464,
            0.49461925,
            1.674643,
            -1.5498532,
            -0.6422282,
            1.0929708,
            -1.1923141,
            0.8744651,
            -0.78989977,
            -0.20583953,
            1.5470413,
            -1.2639258,
            -0.53534126,
            1.0229365,
            -0.40843603,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()
