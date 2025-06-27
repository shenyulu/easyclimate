"""
pytest for physics.lcl.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from easyclimate.physics.condensation.lcl import (
    calc_lifting_condensation_level_bolton1980,
    calc_lifting_condensation_level_Bohren_Albrecht2023,
)


# Test data for Bolton 1980 function
@pytest.fixture
def bolton_test_data():
    temperature = xr.DataArray(
        data=[20, 25, 30],  # Celsius
        dims=["time"],
        coords={"time": [1, 2, 3]},
        attrs={"units": "celsius"},
    )
    rh = xr.DataArray(
        data=[50, 70, 90],  # percentage
        dims=["time"],
        coords={"time": [1, 2, 3]},
        attrs={"units": "%"},
    )
    return temperature, rh


# Test data for Bohren & Albrecht 2023 function
@pytest.fixture
def bohren_test_data():
    pressure = xr.DataArray(
        data=[1000, 950, 900],  # hPa
        dims=["time"],
        coords={"time": [1, 2, 3]},
        attrs={"units": "hPa"},
    )
    temperature = xr.DataArray(
        data=[20, 25, 30],  # Celsius
        dims=["time"],
        coords={"time": [1, 2, 3]},
        attrs={"units": "celsius"},
    )
    dewpoint = xr.DataArray(
        data=[10, 15, 20],  # Celsius
        dims=["time"],
        coords={"time": [1, 2, 3]},
        attrs={"units": "celsius"},
    )
    return pressure, temperature, dewpoint


def test_bolton1980_basic_calculation(bolton_test_data):
    """Test basic calculation with valid inputs"""
    temp, rh = bolton_test_data
    result = calc_lifting_condensation_level_bolton1980(temp, rh, "celsius", "%")

    # Check output type and attributes
    assert isinstance(result, xr.DataArray)
    assert result.attrs["standard_name"] == "atmosphere_lifting_condensation_level"
    assert result.attrs["units"] == "K"


def test_bolton1980_unit_conversions(bolton_test_data):
    """Test different unit inputs"""
    temp, rh = bolton_test_data

    # Test Kelvin input
    temp_k = temp + 273.15
    result_k = calc_lifting_condensation_level_bolton1980(temp_k, rh, "kelvin", "%")

    # Test Fahrenheit input
    temp_f = temp * 9 / 5 + 32
    result_f = calc_lifting_condensation_level_bolton1980(temp_f, rh, "fahrenheit", "%")

    # Test dimensionless RH input (0-1 scale)
    rh_dimless = rh / 100
    result_dimless = calc_lifting_condensation_level_bolton1980(
        temp, rh_dimless, "celsius", "dimensionless"
    )

    # All results should be approximately equal
    xr.testing.assert_allclose(result_k, result_f, rtol=1e-3)
    xr.testing.assert_allclose(result_k, result_dimless, rtol=1e-3)


def test_bolton1980_edge_cases():
    """Test edge cases like extreme values"""
    # Very high RH (near 100%)
    temp = xr.DataArray([20], attrs={"units": "celsius"})
    rh = xr.DataArray([99.9], attrs={"units": "%"})
    result = calc_lifting_condensation_level_bolton1980(temp, rh, "celsius", "%")
    assert result > 54  # Should be above freezing

    # Very low RH
    rh = xr.DataArray([1], attrs={"units": "%"})
    result = calc_lifting_condensation_level_bolton1980(temp, rh, "celsius", "%")
    assert result < 273.15  # Should be below freezing


def test_bohren_albrecht_basic_calculation(bohren_test_data):
    """Test basic calculation with valid inputs"""
    p, t, td = bohren_test_data
    result = calc_lifting_condensation_level_Bohren_Albrecht2023(
        p, t, td, "hPa", "celsius", "celsius"
    )

    # Check output type and structure
    assert isinstance(result, xr.Dataset)
    assert "p_lcl" in result
    assert "t_lcl" in result

    # Check units
    assert result["p_lcl"].attrs["units"] == "hPa"
    assert result["t_lcl"].attrs["units"] == "K"

    # Check values are reasonable
    assert all(result["p_lcl"] > 100)  # Pressure should be > 100 hPa
    assert all(result["p_lcl"] < 1500)  # and < 1500 hPa
    assert all(result["t_lcl"] > 200)  # Temperature should be > 200K
    assert all(result["t_lcl"] < 400)  # and < 400K


def test_bohren_albrecht_unit_conversions(bohren_test_data):
    """Test different unit inputs"""
    p, t, td = bohren_test_data

    # Test Pa pressure input
    p_pa = p * 100
    result_pa = calc_lifting_condensation_level_Bohren_Albrecht2023(
        p_pa, t, td, "Pa", "celsius", "celsius"
    )

    # Test Kelvin temperature inputs
    t_k = t + 273.15
    td_k = td + 273.15
    result_k = calc_lifting_condensation_level_Bohren_Albrecht2023(
        p, t_k, td_k, "hPa", "kelvin", "kelvin"
    )

    # Test Fahrenheit temperature inputs
    t_f = t * 9 / 5 + 32
    td_f = td * 9 / 5 + 32
    result_f = calc_lifting_condensation_level_Bohren_Albrecht2023(
        p, t_f, td_f, "hPa", "fahrenheit", "fahrenheit"
    )

    # All results should be approximately equal
    xr.testing.assert_allclose(result_pa["t_lcl"], result_k["t_lcl"], rtol=1e-3)
    xr.testing.assert_allclose(result_k["t_lcl"], result_f["t_lcl"], rtol=1e-3)
    xr.testing.assert_allclose(result_pa["p_lcl"], result_k["p_lcl"], rtol=1e-3)


def test_bohren_albrecht_saturation_case():
    """Test case where T = Td (already saturated)"""
    p = xr.DataArray([1000], attrs={"units": "hPa"})
    t = xr.DataArray([20], attrs={"units": "celsius"})
    td = t.copy()  # Dewpoint equals temperature

    result = calc_lifting_condensation_level_Bohren_Albrecht2023(
        p, t, td, "hPa", "celsius", "celsius"
    )

    # LCL should be at surface level
    xr.testing.assert_allclose(result["p_lcl"], p)
    xr.testing.assert_allclose(result["t_lcl"], t + 273.15)
