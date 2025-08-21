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
            -0.02228896,
            0.00524893,
            0.05434632,
            -0.02464823,
            -0.01219194,
            0.00126471,
            -0.02075169,
            0.00186409,
            0.04528275,
            0.06256099,
            -0.00727007,
            0.06153464,
            0.02698354,
            -0.02405138,
            -0.02663113,
            -0.05499401,
            -0.00120108,
            0.02052062,
            -0.02202774,
            -0.02223842,
        ]
    )
    assert np.isclose(result_data, refer_data, atol=0.01).all()
