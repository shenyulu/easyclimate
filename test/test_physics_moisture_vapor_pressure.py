"""
pytest for physics.moisture.vapor_pressure.py
"""

import pytest

# import easyclimate as ecl
import numpy as np
import xarray as xr
from easyclimate.physics import calc_vapor_pressure, calc_saturation_vapor_pressure


def test_calc_vapor_pressure():
    result_data = calc_vapor_pressure(
        pressure_data=xr.DataArray([988, 991]),
        mixing_ratio_data=xr.DataArray([0.018, 0.0159]),
        pressure_data_units="hPa",
    ).data
    refer_data = np.array([27.789371, 24.70287576])
    assert np.isclose(result_data, refer_data).all()


def test_calc_saturation_vapor_pressure():
    result_data = calc_saturation_vapor_pressure(
        temperature_data=xr.DataArray([24, 36, 40]), temperature_data_units="celsius"
    ).data
    refer_data = np.array([29.83254305, 59.5118433, 73.94900581])
    assert np.isclose(result_data, refer_data).all()
