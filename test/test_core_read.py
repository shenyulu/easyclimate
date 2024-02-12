"""
pytest for variability.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
import os
from .const_define import TEST_DATA_PATH
from pathlib import Path


def test_open_muliti_dataset():
    result_data = ecl.open_muliti_dataset(
        str(Path(TEST_DATA_PATH, "test_input_core_read*.nc")), dim="time"
    ).sst.data.compute()
    refer_data = np.array(
        [29.347528, 29.097271, 29.414158, 29.765087, 29.57674, 29.308327]
    )
    assert np.isclose(result_data, refer_data).all()
