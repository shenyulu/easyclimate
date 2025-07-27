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
    assert np.isclose(result_data, refer_data, atol=0.01).all()


def test_calc_index_CGT_NH_EOF2_Ding_Wang_2005():
    result_data = ecl.field.teleconnection.calc_index_CGT_NH_EOF2_Ding_Wang_2005(
        z200_data, solver="randomized", random_state=1
    ).data[:20]
    refer_data = np.array(
        [
            -0.11139196,
            -0.02858514,
            -0.06101292,
            -0.00586747,
            -0.03737175,
            -0.02206146,
            0.0107485,
            0.00204435,
            -0.01776662,
            0.01659197,
            -0.03105383,
            -0.01007578,
            0.00146058,
            -0.02168902,
            -0.00332586,
            -0.01091488,
            -0.02162741,
            -0.03120371,
            -0.05922785,
            -0.04308352,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data, atol=0.01).all()
