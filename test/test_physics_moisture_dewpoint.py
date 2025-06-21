"""
pytest for physics.moisture.dewpoint.py
"""

import pytest

import numpy as np
import xarray as xr
from easyclimate.physics import calc_dewpoint


def test_calc_dewpoint():
    result_data = calc_dewpoint(
        vapor_pressure_data=xr.DataArray([22, 23]), vapor_pressure_data_units="hPa"
    ).data
    refer_data = np.array([19.02910179, 19.74308483])
    assert np.isclose(result_data, refer_data).all()
