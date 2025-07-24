"""
pytest for physics/energy/enthalpy.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr

from easyclimate.physics.energy.enthalpy import calc_enthalpy


def test_calc_enthalpy1():
    # Example data
    temp = xr.DataArray(np.array([20, 32.2222]), dims=["time"], name="temperature")
    mix_ratio = xr.DataArray(
        np.array([0.005875, 0.021717]), dims=["time"], name="mixing_ratio"
    )
    result_data = calc_enthalpy(temp, mix_ratio, "degC", "kg/kg").data
    refer_data = np.array([35.109575, 88.15948639])
    assert np.isclose(result_data, refer_data, atol=0.01).all()


def test_calc_enthalpy2():
    # Example data
    temp = xr.DataArray(np.linspace(0, 30, 10), dims=["time"], name="temperature")
    mix_ratio = xr.DataArray(
        np.linspace(0.005, 0.02, 10), dims=["time"], name="mixing_ratio"
    )
    result_data = calc_enthalpy(temp, mix_ratio, "degC", "kg/kg").data
    refer_data = np.array(
        [
            12.5,
            20.07533333,
            27.67166667,
            35.289,
            42.92733333,
            50.58666667,
            58.267,
            65.96833333,
            73.69066667,
            81.434,
        ]
    )
    assert np.isclose(result_data, refer_data, atol=0.01).all()
