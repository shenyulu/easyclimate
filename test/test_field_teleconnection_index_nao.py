"""
pytest for field/teleconnection/index_NAO.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import TEST_DATA_PATH

z500_data = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_z500_mon.nc")))["hgt"]


def test_calc_index_NAO_NH_REOF():
    result_data = ecl.field.teleconnection.calc_index_NAO_NH_REOF(
        z500_data, solver="randomized", random_state=1
    ).data[:20]
    refer_data = np.array(
        [
            -0.5077775140625181,
            0.11957886922929109,
            1.2380946548775624,
            -0.5615255050006085,
            -0.27775159494020857,
            0.028812082567457313,
            -0.47275613504495234,
            0.04246679813223278,
            1.0316122702407973,
            1.425237663384147,
            -0.16562373262685118,
            1.4018558862417667,
            0.6147274649398006,
            -0.5479282336758031,
            -0.6066990256071071,
            -1.252849931045222,
            -0.027362574454467985,
            0.4674920578596379,
            -0.5018266288676613,
            -0.5066262180848115,
        ]
    )
    assert np.isclose(result_data, refer_data).all()
