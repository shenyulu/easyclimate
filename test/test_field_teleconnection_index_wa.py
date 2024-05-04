"""
pytest for field/teleconnection/index_WA.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import TEST_DATA_PATH

z500_data = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_z500_mon.nc")))["hgt"]


def test_calc_index_WA_Wallace_Gutzler_1981():
    result_data = ecl.field.teleconnection.calc_index_WA_Wallace_Gutzler_1981(
        z500_data
    ).data[:20]
    refer_data = np.array(
        [
            0.8997387,
            2.7402146,
            1.7097857,
            -0.27913448,
            0.1245844,
            0.23812951,
            -0.1879792,
            0.00441704,
            0.697674,
            0.67117834,
            0.10534014,
            1.0711093,
            0.3301192,
            -1.46412,
            -1.3281851,
            -0.78094214,
            0.76201576,
            0.06181871,
            -0.14729209,
            -1.0941349,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()
