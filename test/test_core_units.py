"""
pytest for core/units.py
"""

import pytest
import warnings

import numpy as np
import xarray as xr
import easyclimate as ecl

ds3 = xr.DataArray(
    np.array([[3e-5, 5.4e-5], [5.2, -75.5]]),
    dims=("lon", "lat"),
    coords={"lon": np.array([160, 70]), "lat": np.array([87.5, -87.5])},
)


def test_transfer_units_coeff():
    from easyclimate.core.units import transfer_units_coeff

    result_data1 = transfer_units_coeff("m/s", "km/h")
    result_data2 = transfer_units_coeff("hPa", "mbar")
    result_data3 = transfer_units_coeff("mm/day", "m/month")
    refer_data1 = 3.5999999999999996
    refer_data2 = 1.0000000000000002
    refer_data3 = 0.0304375
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()


def test_transfer_data_multiple_units():
    from easyclimate.core.units import transfer_data_multiple_units

    result_data = transfer_data_multiple_units(ds3, "mm/day", "m/day").data.flatten()
    refer_data = np.array([3.00e-08, 5.40e-08, 5.20e-03, -7.55e-02])
    assert np.isclose(result_data, refer_data).all()


def test_transfer_data_difference_units():
    from easyclimate.core.units import transfer_data_difference_units

    result_data = transfer_data_difference_units(
        xr.DataArray(15), "celsius", "kelvin"
    ).data
    refer_data = 288.15
    assert np.isclose(result_data, refer_data).all()


def test_transfer_data_units():
    from easyclimate.core.units import transfer_data_units

    result_data = transfer_data_units(
        xr.DataArray([104, 100, 92, 92, 86, 80, 80, 60, 30]), "degC", "degF"
    ).data
    refer_data = np.array(
        [219.2, 212.0, 197.6, 197.6, 186.8, 176.0, 176.0, 140.0, 86.0]
    )
    assert np.isclose(result_data, refer_data).all()
