"""
pytest for field/teleconnection/index_CGT.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import TEST_DATA_PATH

z200_data = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_z200_mon.nc")))["hgt"]


def test_calc_index_CGT_1point_Ding_Wang_2005():
    result_data = ecl.field.teleconnection.calc_index_CGT_1point_Ding_Wang_2005(
        z200_data
    ).data[:20]
    refer_data = np.array(
        [
            3.7474614e-01,
            -1.7882865e-02,
            2.4786426e-02,
            5.6887108e-01,
            6.2733793e-01,
            -1.6665318e00,
            8.7259931e-04,
            7.5344339e-02,
            -2.5329903e-01,
            -1.2525458e00,
            1.1981032e-01,
            -3.7184653e-01,
            -5.5903262e-01,
            -1.9215151e00,
            -8.5765535e-01,
            1.2265400e00,
            4.5251453e-01,
            -3.0124459e-01,
            -1.1771039e00,
            1.8773080e-01,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_index_CGT_NH_EOF2_Ding_Wang_2005():
    result_data = ecl.field.teleconnection.calc_index_CGT_NH_EOF2_Ding_Wang_2005(
        z200_data, solver="randomized", random_state=1
    ).data[:20]
    refer_data = np.array(
        [
            -2.53768398,
            -0.65121434,
            -1.38997029,
            -0.13367014,
            -0.85138727,
            -0.50259467,
            0.24486781,
            0.04657358,
            -0.40475156,
            0.37799115,
            -0.70745496,
            -0.22954206,
            0.03327428,
            -0.49410987,
            -0.07576824,
            -0.24865819,
            -0.49270629,
            -0.71086944,
            -1.34930355,
            -0.98151028,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()
