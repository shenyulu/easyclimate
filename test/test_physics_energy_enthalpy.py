"""
pytest for physics/enthalpy.py
"""

import pytest
import numpy as np
import xarray as xr
from easyclimate.physics.energy.enthalpy import calc_enthalpy


def test_calc_enthalpy():
    temp = xr.DataArray(np.array([20, 32.2222]), dims=["time"], name="temperature")
    mix_ratio = xr.DataArray(
        np.array([0.005875, 0.021717]), dims=["time"], name="mixing_ratio"
    )

    # Calculate enthalpy (input in °C, kg/kg; output in J/kg)
    result_data = calc_enthalpy(temp, mix_ratio, "degC", "g/g").data
    refer_data = np.array([35.109575, 88.15948639])
    assert np.isclose(result_data, refer_data).all()
