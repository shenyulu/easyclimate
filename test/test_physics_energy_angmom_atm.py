"""
pytest for physics/energy/angmom_atm.py
"""

import pytest
import numpy as np
import xarray as xr
import easyclimate as ecl
from easyclimate.physics.energy.angmom_atm import calc_relative_angular_momentum


def test_calc_relative_angular_momentum():
    udata = ecl.tutorial.open_tutorial_dataset("uwnd_202201_mon_mean").uwnd.isel(time=0)
    result_data = ecl.physics.energy.calc_relative_angular_momentum(
        udata, vertical_dim="level", vertical_dim_units="hPa"
    ).data
    refer_data = np.array([-1.60258e26])
    assert np.isclose(result_data, refer_data, atol=0.01).all()
