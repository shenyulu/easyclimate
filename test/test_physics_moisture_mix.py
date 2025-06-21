"""
pytest for physics.moisture.mix.py
"""

import pytest

import numpy as np
import xarray as xr
from easyclimate.physics import calc_mixing_ratio, calc_saturation_mixing_ratio
from .util import round_sf_np


def test_calc_mixing_ratio():
    result_data = calc_mixing_ratio(
        partial_pressure_data=xr.DataArray([25, 30]),
        total_pressure_data=xr.DataArray([1000, 1000]),
    ).data
    refer_data = np.array([0.01594761, 0.01923578])
    assert np.isclose(result_data, refer_data).all()


def test_calc_saturation_mixing_ratio():
    result_data = calc_saturation_mixing_ratio(
        total_pressure_data=xr.DataArray([983, 991]),
        temperature_data=xr.DataArray([25, 28]),
        temperature_data_units="celsius",
        total_pressure_data_units="hPa",
    ).data
    refer_data = np.array([0.02070799, 0.02467106])
    assert np.isclose(result_data, refer_data).all()
