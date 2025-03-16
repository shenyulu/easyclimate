"""
pytest for field.detection.aerobulk.py
"""

import pytest
import easyclimate as ecl
import xarray as xr
import numpy as np


def test_calc_turbulent_fluxes_without_skin_correction():
    tmp = ecl.field.boundary_layer.calc_turbulent_fluxes_without_skin_correction(
        xr.DataArray(np.array([273.15 + 22.0, 273.15 + 22.0])),
        "degC",
        xr.DataArray(np.array([273.15 + 20.0, 273.15 + 20.0])),
        "degC",
        xr.DataArray(np.array([0.012, 0.012])),
        "g/kg",
        xr.DataArray(np.array([4, 4])),
        "m/s",
        xr.DataArray(np.array([9, 9])),
        "m/s",
        xr.DataArray(np.array([101000.0, 101000.0])),
        "Pa",
        algorithm="ncar",
    )
    result_data = tmp["ql"].data[0]
    refer_data = np.float64([4082.766776880013])
    assert np.isclose(result_data, refer_data).all()


def test_calc_turbulent_fluxes_skin_correction():
    tmp = ecl.field.boundary_layer.calc_turbulent_fluxes_skin_correction(
        xr.DataArray(
            np.array([273.15 + 22.0, 273.15 + 22.0]),
        ),
        "degC",
        xr.DataArray(np.array([273.15 + 20.0, 273.15 + 20.0])),
        "degC",
        xr.DataArray(np.array([0.012, 0.012])),
        "g/kg",
        xr.DataArray(np.array([4, 4])),
        "m/s",
        xr.DataArray(np.array([9, 9])),
        "m/s",
        xr.DataArray(np.array([101000.0, 101000.0])),
        "Pa",
        xr.DataArray(np.array([0, 0])),
        "W/m^2",
        xr.DataArray(np.array([350.0, 350.0])),
        "W/m^2",
        algorithm="coare3p6",
    )
    result_data = tmp["t_s"].data[0]
    refer_data = np.float64([566.6892344670376])
    assert np.isclose(result_data, refer_data).all()
