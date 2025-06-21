"""
pytest for physics/moisture/wet_bulb.py
"""

import pytest
import numpy as np
import xarray as xr
from easyclimate.physics.moisture.wet_bulb import (
    calc_wet_bulb_potential_temperature_iteration,
    calc_wet_bulb_potential_temperature_davies_jones2008,
    calc_wet_bulb_temperature_stull2011,
    calc_wet_bulb_temperature_sadeghi2013,
)


class TestCalcWetBulbPotentialTemperatureIteration:
    """Test cases for calc_wet_bulb_potential_temperature_iteration function."""

    def test_basic_calculation(self):
        """Test basic calculation with simple inputs."""
        temp = xr.DataArray(np.array([20, 25, 30]), dims=["point"])
        rh = xr.DataArray(np.array([50, 60, 70]), dims=["point"])
        pressure = xr.DataArray(np.array([1000, 950, 900]), dims=["point"])

        result = calc_wet_bulb_potential_temperature_iteration(
            temperature_data=temp,
            relative_humidity_data=rh,
            pressure_data=pressure,
            temperature_data_units="celsius",
            relative_humidity_data_units="%",
            pressure_data_units="hPa",
        )

        assert isinstance(result, xr.DataArray)
        assert result.dims == ("point",)
        assert result.attrs["units"] == "degC"
        assert result.attrs["standard_name"] == "wet_bulb_temperature"

    def test_2d_input(self):
        """Test with 2D input data."""
        temp_2d = xr.DataArray(np.random.rand(2, 2) * 30, dims=["lat", "lon"])
        rh_2d = xr.DataArray(np.random.rand(2, 2) * 100, dims=["lat", "lon"])
        pres_2d = xr.DataArray(np.random.rand(2, 2) * 200 + 800, dims=["lat", "lon"])

        result = calc_wet_bulb_potential_temperature_iteration(
            temp_2d, rh_2d, pres_2d, "celsius", "%", "hPa"
        )

        assert isinstance(result, xr.DataArray)
        assert result.dims == ("lat", "lon")

    def test_unit_conversions(self):
        """Test that unit conversions work properly."""
        temp = xr.DataArray(
            np.array([293.15, 298.15]), dims=["point"]
        )  # 20C, 25C in Kelvin
        rh = xr.DataArray(
            np.array([0.5, 0.6]), dims=["point"]
        )  # 50%, 60% as dimensionless
        pressure = xr.DataArray(
            np.array([100000, 95000]), dims=["point"]
        )  # 1000hPa, 950hPa in Pa

        result = calc_wet_bulb_potential_temperature_iteration(
            temperature_data=temp,
            relative_humidity_data=rh,
            pressure_data=pressure,
            temperature_data_units="kelvin",
            relative_humidity_data_units="dimensionless",
            pressure_data_units="Pa",
        )

        assert isinstance(result, xr.DataArray)
        assert result.attrs["units"] == "degC"


class TestCalcWetBulbPotentialTemperatureDaviesJones2008:
    """Test cases for calc_wet_bulb_potential_temperature_davies_jones2008 function."""

    def test_basic_calculation(self):
        """Test basic calculation with simple inputs."""
        pressure = xr.DataArray(np.array([1000, 950]), dims=["point"])
        temp = xr.DataArray(np.array([293.15, 298.15]), dims=["point"])  # 20C, 25C
        dewpoint = xr.DataArray(np.array([283.15, 288.15]), dims=["point"])  # 10C, 15C

        result = calc_wet_bulb_potential_temperature_davies_jones2008(
            pressure_data=pressure,
            temperature_data=temp,
            dewpoint_data=dewpoint,
            pressure_data_units="hPa",
            temperature_data_units="kelvin",
            dewpoint_data_units="kelvin",
        )

        assert isinstance(result, xr.DataArray)
        assert result.dims == ("point",)
        assert result.attrs["units"] == "K"
        assert result.attrs["standard_name"] == "wet_bulb_temperature"

    def test_unit_conversions(self):
        """Test that unit conversions work properly."""
        pressure = xr.DataArray(
            np.array([100000, 95000]), dims=["point"]
        )  # 1000hPa, 950hPa in Pa
        temp = xr.DataArray(np.array([20, 25]), dims=["point"])  # in Celsius
        dewpoint = xr.DataArray(np.array([10, 15]), dims=["point"])  # in Celsius

        result = calc_wet_bulb_potential_temperature_davies_jones2008(
            pressure_data=pressure,
            temperature_data=temp,
            dewpoint_data=dewpoint,
            pressure_data_units="Pa",
            temperature_data_units="celsius",
            dewpoint_data_units="celsius",
        )

        assert isinstance(result, xr.DataArray)
        assert result.attrs["units"] == "K"


class TestCalcWetBulbTemperatureStull2011:
    """Test cases for calc_wet_bulb_temperature_stull2011 function."""

    def test_basic_calculation(self):
        """Test basic calculation with simple inputs."""
        temp = xr.DataArray(np.array([20, 25]), dims=["point"])
        rh = xr.DataArray(np.array([50, 60]), dims=["point"])

        result = calc_wet_bulb_temperature_stull2011(
            temperature_data=temp,
            relative_humidity_data=rh,
            temperature_data_units="celsius",
            relative_humidity_data_units="%",
        )

        assert isinstance(result, xr.DataArray)
        assert result.dims == ("point",)
        assert result.attrs["units"] == "K"
        assert result.attrs["standard_name"] == "wet_bulb_temperature"

    def test_unit_conversions(self):
        """Test that unit conversions work properly."""
        temp = xr.DataArray(
            np.array([293.15, 298.15]), dims=["point"]
        )  # 20C, 25C in Kelvin
        rh = xr.DataArray(
            np.array([0.5, 0.6]), dims=["point"]
        )  # 50%, 60% as dimensionless

        result = calc_wet_bulb_temperature_stull2011(
            temperature_data=temp,
            relative_humidity_data=rh,
            temperature_data_units="kelvin",
            relative_humidity_data_units="dimensionless",
        )

        assert isinstance(result, xr.DataArray)
        assert result.attrs["units"] == "K"


class TestCalcWetBulbTemperatureSadeghi2013:
    """Test cases for calc_wet_bulb_temperature_sadeghi2013 function."""

    def test_basic_calculation(self):
        """Test basic calculation with simple inputs."""
        temp = xr.DataArray(np.array([20, 25]), dims=["point"])
        height = xr.DataArray(np.array([100, 200]), dims=["point"])
        rh = xr.DataArray(np.array([50, 60]), dims=["point"])

        result = calc_wet_bulb_temperature_sadeghi2013(
            temperature_data=temp,
            height_data=height,
            relative_humidity_data=rh,
            temperature_data_units="celsius",
            height_data_units="m",
            relative_humidity_data_units="%",
        )

        assert isinstance(result, xr.DataArray)
        assert result.dims == ("point",)
        assert result.attrs["units"] == "degC"
        assert result.attrs["standard_name"] == "wet_bulb_temperature"

    def test_unit_conversions(self):
        """Test that unit conversions work properly."""
        temp = xr.DataArray(
            np.array([293.15, 298.15]), dims=["point"]
        )  # 20C, 25C in Kelvin
        height = xr.DataArray(np.array([0.1, 0.2]), dims=["point"])  # 100m, 200m in km
        rh = xr.DataArray(
            np.array([0.5, 0.6]), dims=["point"]
        )  # 50%, 60% as dimensionless

        result = calc_wet_bulb_temperature_sadeghi2013(
            temperature_data=temp,
            height_data=height,
            relative_humidity_data=rh,
            temperature_data_units="kelvin",
            height_data_units="km",
            relative_humidity_data_units="dimensionless",
        )

        assert isinstance(result, xr.DataArray)
        assert result.attrs["units"] == "degC"
