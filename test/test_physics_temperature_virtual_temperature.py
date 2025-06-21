"""
pytest for physics.temperature.virtual_temperature.py
"""

import pytest

import numpy as np
from easyclimate.physics import (
    calc_virtual_temperature,
    calc_virtual_temperature_Hobbs2006,
)
from .const_data import q_data_Tv, t_data_Tv


def test_calc_virtual_temperature(q_data_Tv, t_data_Tv):
    result_data = calc_virtual_temperature(
        temper_data=t_data_Tv,
        specific_humidity_data=q_data_Tv,
        specific_humidity_data_units="g/g",
    ).data
    refer_data = np.array([303.1976896, 298.08551685])
    assert np.isclose(result_data, refer_data).all()


def test_calc_virtual_temperature_Hobbs2006(q_data_Tv, t_data_Tv):
    result_data = calc_virtual_temperature_Hobbs2006(
        temper_data=t_data_Tv,
        specific_humidity_data=q_data_Tv,
        specific_humidity_data_units="g/g",
    ).data
    refer_data = np.array([303.1345747, 298.04520656])
    assert np.isclose(result_data, refer_data).all()
