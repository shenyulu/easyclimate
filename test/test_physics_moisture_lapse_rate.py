"""
pytest for physics/moisture/lapse_rate.py
"""

import pytest
import xarray as xr
from easyclimate.physics.moisture.lapse_rate import calc_moist_adiabatic_lapse_rate
from pint.errors import UndefinedUnitError


class TestCalcMoistAdiabaticLapseRate:
    """Test suite for calc_moist_adiabatic_lapse_rate function."""

    @pytest.fixture
    def sample_data(self):
        """Create sample xarray DataArrays for testing."""
        pressure_data = xr.DataArray([1000, 850, 700], dims=["level"])
        temperature_data = xr.DataArray([20, 10, 0], dims=["level"])
        return pressure_data, temperature_data

    def test_returns_dataarray(self, sample_data):
        """Test that function returns an xarray DataArray."""
        pressure, temp = sample_data
        result = calc_moist_adiabatic_lapse_rate(pressure, temp, "hPa", "celsius")
        assert isinstance(result, xr.DataArray)

    def test_output_units(self, sample_data):
        """Test that output has correct units attribute."""
        pressure, temp = sample_data
        result = calc_moist_adiabatic_lapse_rate(pressure, temp, "hPa", "celsius")
        assert result.attrs["units"] == "K/hPa"

    def test_pressure_unit_conversions(self, sample_data):
        """Test function works with different pressure units."""
        pressure, temp = sample_data
        for unit in ["hPa", "Pa", "mbar"]:
            result = calc_moist_adiabatic_lapse_rate(pressure, temp, unit, "celsius")
            assert isinstance(result, xr.DataArray)

    def test_temperature_unit_conversions(self, sample_data):
        """Test function works with different temperature units."""
        pressure, temp = sample_data
        for unit in ["celsius", "kelvin", "fahrenheit"]:
            result = calc_moist_adiabatic_lapse_rate(pressure, temp, "hPa", unit)
            assert isinstance(result, xr.DataArray)

    def test_output_shape_matches_input(self, sample_data):
        """Test output shape matches input shape."""
        pressure, temp = sample_data
        result = calc_moist_adiabatic_lapse_rate(pressure, temp, "hPa", "celsius")
        assert result.shape == pressure.shape

    def test_invalid_pressure_unit_raises_error(self, sample_data):
        """Test invalid pressure unit raises error."""
        pressure, temp = sample_data
        with pytest.raises(UndefinedUnitError):
            calc_moist_adiabatic_lapse_rate(pressure, temp, "invalid_unit", "celsius")

    def test_invalid_temperature_unit_raises_error(self, sample_data):
        """Test invalid temperature unit raises error."""
        pressure, temp = sample_data
        with pytest.raises(ValueError):
            calc_moist_adiabatic_lapse_rate(pressure, temp, "hPa", "invalid_unit")

    def test_output_values_are_positive(self, sample_data):
        """Test all output values are positive (as lapse rates should be)."""
        pressure, temp = sample_data
        result = calc_moist_adiabatic_lapse_rate(pressure, temp, "hPa", "celsius")
        assert (result > 0).all()

    def test_output_decreases_with_pressure(self):
        """Test that lapse rate increases with decreasing pressure."""
        pressure = xr.DataArray([1000, 900, 800, 700], dims=["level"])
        temp = xr.DataArray([15, 10, 5, 0], dims=["level"])
        result = calc_moist_adiabatic_lapse_rate(pressure, temp, "hPa", "celsius")
        # Lapse rate should increase with height (decreasing pressure)
        assert (result.diff("level") > 0).all()
