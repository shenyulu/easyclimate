"""
pytest for physics.temperature.equivalent_potential_temperature.py
"""

import pytest

import numpy as np
import xarray as xr
from easyclimate.physics import (
    calc_equivalent_potential_temperature,
    calc_equivalent_potential_temperature_davies_jones2009,
)
from pint.errors import UndefinedUnitError
from .util import round_sf_np


class TestCalcEquivalentPotentialTemperatureBolton1980:
    # Fixture for creating test DataArrays
    @pytest.fixture
    def sample_data_calc_equivalent_potential_temperature(self):
        pressure = xr.DataArray(
            np.array([1000, 950, 900, 850]),
            dims=["level"],
            coords={"level": [1000, 950, 900, 850]},
            name="pressure",
        )
        temperature = xr.DataArray(
            np.array([25, 20, 15, 10]),
            dims=["level"],
            coords={"level": [1000, 950, 900, 850]},
            name="temperature",
        )
        dewpoint = xr.DataArray(
            np.array([20, 15, 10, 5]),
            dims=["level"],
            coords={"level": [1000, 950, 900, 850]},
            name="dewpoint",
        )
        return pressure, temperature, dewpoint

    def test_basic_calculation_calc_equivalent_potential_temperature(
        self,
        sample_data_calc_equivalent_potential_temperature,
    ):
        """Test basic calculation with hPa and Celsius units"""
        pressure, temperature, dewpoint = (
            sample_data_calc_equivalent_potential_temperature
        )

        result = calc_equivalent_potential_temperature(
            pressure_data=pressure,
            temperature_data=temperature,
            dewpoint_data=dewpoint,
            pressure_data_units="hPa",
            temperature_data_units="celsius",
            dewpoint_data_units="celsius",
        )

        assert isinstance(result, xr.DataArray)
        assert result.name == "theta_e"
        assert result.attrs["units"] == "K"
        assert result.shape == temperature.shape

    def test_different_units_calc_equivalent_potential_temperature(
        self,
        sample_data_calc_equivalent_potential_temperature,
    ):
        """Test with different input units"""
        pressure, temperature, dewpoint = (
            sample_data_calc_equivalent_potential_temperature
        )

        # Test with Pa and Kelvin
        result_pa_k = calc_equivalent_potential_temperature(
            pressure_data=pressure * 100,  # convert to Pa
            temperature_data=temperature + 273.15,  # convert to Kelvin
            dewpoint_data=dewpoint + 273.15,  # convert to Kelvin
            pressure_data_units="Pa",
            temperature_data_units="kelvin",
            dewpoint_data_units="kelvin",
        )

        # Test with Fahrenheit
        result_f = calc_equivalent_potential_temperature(
            pressure_data=pressure,
            temperature_data=(temperature * 9 / 5) + 32,  # convert to Fahrenheit
            dewpoint_data=(dewpoint * 9 / 5) + 32,  # convert to Fahrenheit
            pressure_data_units="hPa",
            temperature_data_units="fahrenheit",
            dewpoint_data_units="fahrenheit",
        )

        # All results should be approximately equal
        xr.testing.assert_allclose(result_pa_k, result_f, rtol=1e-3)

    def test_known_value_calc_equivalent_potential_temperature(self):
        """Test against a known value from literature or manual calculation"""
        # Using example values that should produce a known result
        pressure = xr.DataArray([1000.0], dims=["level"])
        temperature = xr.DataArray([25.0], dims=["level"])
        dewpoint = xr.DataArray([20.0], dims=["level"])

        result = calc_equivalent_potential_temperature(
            pressure_data=pressure,
            temperature_data=temperature,
            dewpoint_data=dewpoint,
            pressure_data_units="hPa",
            temperature_data_units="celsius",
            dewpoint_data_units="celsius",
        )

        # This expected value should be replaced with a known correct value
        expected = xr.DataArray([340.0], dims=["level"])  # Approximate value
        xr.testing.assert_allclose(result, expected, rtol=0.1)

    def test_edge_cases_calc_equivalent_potential_temperature(self):
        """Test edge cases like very low/high temperatures"""
        pressure = xr.DataArray([1000.0, 1000.0], dims=["level"])

        # Very cold case
        temperature_cold = xr.DataArray([-40.0, -50.0], dims=["level"])
        dewpoint_cold = xr.DataArray([-45.0, -55.0], dims=["level"])

        result_cold = calc_equivalent_potential_temperature(
            pressure_data=pressure,
            temperature_data=temperature_cold,
            dewpoint_data=dewpoint_cold,
            pressure_data_units="hPa",
            temperature_data_units="celsius",
            dewpoint_data_units="celsius",
        )

        # Very hot case
        temperature_hot = xr.DataArray([40.0, 45.0], dims=["level"])
        dewpoint_hot = xr.DataArray([35.0, 40.0], dims=["level"])

        result_hot = calc_equivalent_potential_temperature(
            pressure_data=pressure,
            temperature_data=temperature_hot,
            dewpoint_data=dewpoint_hot,
            pressure_data_units="hPa",
            temperature_data_units="celsius",
            dewpoint_data_units="celsius",
        )

        # Just check that calculations complete without errors
        assert isinstance(result_cold, xr.DataArray)
        assert isinstance(result_hot, xr.DataArray)

    def test_dimensional_consistency_calc_equivalent_potential_temperature(self):
        """Test that input arrays must have same dimensions"""
        pressure = xr.DataArray([1000.0, 950.0], dims=["level"])
        temperature = xr.DataArray([25.0, 20.0, 15.0], dims=["level"])  # different size
        dewpoint = xr.DataArray([20.0, 15.0], dims=["level"])

        with pytest.raises(ValueError):
            calc_equivalent_potential_temperature(
                pressure_data=pressure,
                temperature_data=temperature,
                dewpoint_data=dewpoint,
                pressure_data_units="hPa",
                temperature_data_units="celsius",
                dewpoint_data_units="celsius",
            )


class TestEquivalentPotentialTemperatureCalculator:
    @pytest.fixture
    def sample_data(self):
        """Create sample xarray DataArrays for testing"""
        pressure = xr.DataArray(
            [1000, 950, 900], dims=["level"], attrs={"units": "hPa"}
        )
        temperature = xr.DataArray(
            [25, 20, 15], dims=["level"], attrs={"units": "celsius"}
        )
        dewpoint = xr.DataArray(
            [20, 15, 10], dims=["level"], attrs={"units": "celsius"}
        )
        return pressure, temperature, dewpoint

    def test_calc_equivalent_potential_temperature_basic(self, sample_data):
        """Test basic functionality with valid inputs"""
        pressure, temperature, dewpoint = sample_data

        result = calc_equivalent_potential_temperature_davies_jones2009(
            pressure_data=pressure,
            temperature_data=temperature,
            dewpoint_data=dewpoint,
            pressure_data_units="hPa",
            temperature_data_units="celsius",
            dewpoint_data_units="celsius",
        )

        assert isinstance(result, xr.DataArray)
        assert result.name == "theta_e"
        assert result.attrs["units"] == "K"
        assert len(result) == len(pressure)
        assert all(result > 0)  # All values should be positive

    def test_calc_equivalent_potential_temperature_units(self, sample_data):
        """Test with different unit inputs"""
        pressure, temperature, dewpoint = sample_data

        # Test with Pa instead of hPa
        result_pa = calc_equivalent_potential_temperature_davies_jones2009(
            pressure_data=pressure * 100,
            temperature_data=temperature,
            dewpoint_data=dewpoint,
            pressure_data_units="Pa",
            temperature_data_units="celsius",
            dewpoint_data_units="celsius",
        )

        # Test with Kelvin inputs
        result_kelvin = calc_equivalent_potential_temperature_davies_jones2009(
            pressure_data=pressure,
            temperature_data=temperature + 273.15,
            dewpoint_data=dewpoint + 273.15,
            pressure_data_units="hPa",
            temperature_data_units="kelvin",
            dewpoint_data_units="kelvin",
        )

        # Test with Fahrenheit inputs
        result_fahrenheit = calc_equivalent_potential_temperature_davies_jones2009(
            pressure_data=pressure,
            temperature_data=temperature * 9 / 5 + 32,
            dewpoint_data=dewpoint * 9 / 5 + 32,
            pressure_data_units="hPa",
            temperature_data_units="fahrenheit",
            dewpoint_data_units="fahrenheit",
        )

        # All results should be similar (within some tolerance due to floating point operations)
        xr.testing.assert_allclose(result_pa, result_kelvin, rtol=1e-3)
        xr.testing.assert_allclose(result_kelvin, result_fahrenheit, rtol=1e-3)

    def test_calc_equivalent_potential_temperature_edge_cases(self):
        """Test edge cases like very low temperatures"""

        # Very cold case
        pressure = xr.DataArray([1000], attrs={"units": "hPa"})
        temperature = xr.DataArray([-40], attrs={"units": "celsius"})
        dewpoint = xr.DataArray([-45], attrs={"units": "celsius"})

        result = calc_equivalent_potential_temperature_davies_jones2009(
            pressure_data=pressure,
            temperature_data=temperature,
            dewpoint_data=dewpoint,
            pressure_data_units="hPa",
            temperature_data_units="celsius",
            dewpoint_data_units="celsius",
        )

        assert result[0] > 0  # Should still be positive

    def test_calc_equivalent_potential_temperature_invalid_inputs(self, sample_data):
        """Test handling of invalid inputs"""
        pressure, temperature, dewpoint = sample_data

        # Test invalid pressure units
        with pytest.raises(UndefinedUnitError):
            calc_equivalent_potential_temperature_davies_jones2009(
                pressure_data=pressure,
                temperature_data=temperature,
                dewpoint_data=dewpoint,
                pressure_data_units="invalid",
                temperature_data_units="celsius",
                dewpoint_data_units="celsius",
            )

        # Test invalid temperature units
        with pytest.raises(ValueError):
            calc_equivalent_potential_temperature_davies_jones2009(
                pressure_data=pressure,
                temperature_data=temperature,
                dewpoint_data=dewpoint,
                pressure_data_units="hPa",
                temperature_data_units="invalid",
                dewpoint_data_units="celsius",
            )

        # Test mismatched array dimensions
        with pytest.raises(ValueError):
            calc_equivalent_potential_temperature_davies_jones2009(
                pressure_data=pressure,
                temperature_data=temperature[:-1],  # Different length
                dewpoint_data=dewpoint,
                pressure_data_units="hPa",
                temperature_data_units="celsius",
                dewpoint_data_units="celsius",
            )

    def test_calc_equivalent_potential_temperature_known_value(self):
        """Test against a known value (if available)"""
        # Example test case with known output (adjust values as needed)
        pressure = xr.DataArray([1000], attrs={"units": "hPa"})
        temperature = xr.DataArray([25], attrs={"units": "celsius"})
        dewpoint = xr.DataArray([20], attrs={"units": "celsius"})

        result = calc_equivalent_potential_temperature_davies_jones2009(
            pressure_data=pressure,
            temperature_data=temperature,
            dewpoint_data=dewpoint,
            pressure_data_units="hPa",
            temperature_data_units="celsius",
            dewpoint_data_units="celsius",
        )

        # Check if result is within expected range for this input
        assert pytest.approx(result[0].item(), rel=0.1) == 340.0
