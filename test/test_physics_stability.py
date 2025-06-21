"""
pytest for physics.stability
"""

import pytest

import numpy as np
from easyclimate.physics import calc_brunt_vaisala_frequency_atm, calc_static_stability
from .const_data import pv_data, z_data, t_data
from .util import round_sf_np


def test_calc_brunt_vaisala_frequency_atm(z_data, pv_data):
    N2_data = calc_brunt_vaisala_frequency_atm(
        potential_temperature_data=pv_data, z_data=z_data, vertical_dim="level"
    )

    result_data = round_sf_np(N2_data).data
    refer_data = np.array(
        [
            0.00945,
            0.01092,
            0.0122,
            0.01232,
            0.01239,
            0.01233,
            0.01125,
            0.009368,
            0.007528,
            0.007636,
            0.0123,
            0.02002,
            0.02498,
            0.02534,
            0.02418,
            0.0241,
            0.02189,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


# calc_static_stability
def test_calc_static_stability(t_data):
    x = calc_static_stability(t_data, vertical_dim="level", vertical_dim_units="hPa")

    result_data = round_sf_np(x).data
    refer_data = np.array(
        [
            0.0002514,
            0.0003402,
            0.0004607,
            0.0005096,
            0.0005877,
            0.0006617,
            0.0006387,
            0.0004954,
            0.0003722,
            0.0004301,
            0.00128,
            0.004139,
            0.009856,
            0.01633,
            0.02436,
            0.0444,
            0.06889,
        ]
    )
    assert np.isclose(result_data, refer_data).all()
