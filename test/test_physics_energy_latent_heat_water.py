"""
pytest for physics/latent_heat_water.py
"""

import pytest
import numpy as np
import xarray as xr
from easyclimate.physics.energy.latent_heat_water import calc_latent_heat_water
from easyclimate.core.utility import transfer_data_multiple_units


# https://www.ncl.ucar.edu/Document/Functions/Contributed/latent_heat_water.shtml
def test_calc_latent_heat_water1():
    lh = calc_latent_heat_water(xr.DataArray([0]), "degC", "melting_freezing")
    result_data = transfer_data_multiple_units(lh, "J/kg", "J/g").data
    refer_data = np.array([333.14900488])
    assert np.isclose(result_data, refer_data).all()


def test_calc_latent_heat_water2():
    lh = calc_latent_heat_water(xr.DataArray([0]), "degC", "evaporation_condensation")
    result_data = transfer_data_multiple_units(lh, "J/kg", "J/g").data
    refer_data = np.array([2500.72402552])
    assert np.isclose(result_data, refer_data).all()


def test_calc_latent_heat_water3():
    lh = calc_latent_heat_water(xr.DataArray([0]), "degC", "sublimation_deposition")
    result_data = transfer_data_multiple_units(lh, "J/kg", "J/g").data
    refer_data = np.array([2833.80703553])
    assert np.isclose(result_data, refer_data).all()
