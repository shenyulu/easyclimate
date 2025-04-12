"""
pytest for field.detection.aerobulk.py
"""

import pytest
import easyclimate as ecl
import xarray as xr
import numpy as np
import warnings
from easyclimate.field.boundary_layer.aerobulk import check_value_range


def test_calc_turbulent_fluxes_without_skin_correction():
    tmp = ecl.field.boundary_layer.calc_turbulent_fluxes_without_skin_correction(
        sst_data=xr.DataArray(np.array([273.15 + 22.0, 273.15 + 22.0])),
        sst_data_units="degC",
        absolute_temperature_data=xr.DataArray(
            np.array([273.15 + 20.0, 273.15 + 20.0])
        ),
        absolute_temperature_data_units="degC",
        specific_humidity_data=xr.DataArray(np.array([0.012, 0.012])),
        specific_humidity_data_units="g/kg",
        zonal_wind_speed_data=xr.DataArray(np.array([4, 4])),
        zonal_wind_speed_data_units="m/s",
        meridional_wind_speed_data=xr.DataArray(np.array([9, 9])),
        meridional_wind_speed_data_units="m/s",
        mean_sea_level_pressure_data=xr.DataArray(np.array([101000.0, 101000.0])),
        mean_sea_level_pressure_data_units="Pa",
        algorithm="ncar",
    )
    result_data = tmp["ql"].data[0]
    refer_data = np.float64([4082.766776880013])
    assert np.isclose(result_data, refer_data).all()


def test_calc_turbulent_fluxes_skin_correction():
    tmp = ecl.field.boundary_layer.calc_turbulent_fluxes_skin_correction(
        sst_data=xr.DataArray(
            np.array([273.15 + 22.0, 273.15 + 22.0]),
        ),
        sst_data_units="degC",
        absolute_temperature_data=xr.DataArray(
            np.array([273.15 + 20.0, 273.15 + 20.0])
        ),
        absolute_temperature_data_units="degC",
        specific_humidity_data=xr.DataArray(np.array([0.012, 0.012])),
        specific_humidity_data_units="g/kg",
        zonal_wind_speed_data=xr.DataArray(np.array([4, 4])),
        zonal_wind_speed_data_units="m/s",
        meridional_wind_speed_data=xr.DataArray(np.array([9, 9])),
        meridional_wind_speed_data_units="m/s",
        mean_sea_level_pressure_data=xr.DataArray(np.array([101000.0, 101000.0])),
        mean_sea_level_pressure_data_units="Pa",
        downwelling_shortwave_radiation=xr.DataArray(np.array([0, 0])),
        downwelling_shortwave_radiation_units="W/m^2",
        downwelling_longwave_radiation=xr.DataArray(np.array([350.0, 350.0])),
        downwelling_longwave_radiation_units="W/m^2",
        algorithm="coare3p6",
    )
    result_data = tmp["t_s"].data[0]
    refer_data = np.float64([566.6892344670376])
    assert np.isclose(result_data, refer_data).all()


# Fixture for the valid ranges (if needed in multiple tests)
@pytest.fixture
def valid_ranges():
    return {
        "sst": np.array([270, 320]),
        "t_zt": np.array([180, 330]),
        "hum_zt": np.array([0, 0.08]),
        "u_zu": np.array([-50, 50]),
        "v_zu": np.array([-50, 50]),
        "slp": np.array([80000, 110000]),
        "rad_sw": np.array([0, 1500]),
        "rad_lw": np.array([0, 750]),
    }


def test_valid_values():
    """Test values within valid ranges"""
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # Treat warnings as errors
        # Test boundary and middle values
        check_value_range("sst", xr.DataArray(270))  # Lower bound
        check_value_range("sst", xr.DataArray(320))  # Upper bound
        check_value_range("sst", xr.DataArray(300))  # Middle value
        check_value_range("hum_zt", xr.DataArray(0.04))
        check_value_range("u_zu", xr.DataArray(0))


def test_invalid_values():
    """Test values outside valid ranges"""
    # Test that warnings are raised for invalid values
    with pytest.warns(UserWarning):
        check_value_range("sst", xr.DataArray(250))  # Below minimum
        check_value_range("sst", xr.DataArray(330))  # Above maximum
        check_value_range("hum_zt", xr.DataArray(-0.1))
        check_value_range("u_zu", xr.DataArray(60))


def test_unknown_variable():
    """Test with undefined variable name"""
    with pytest.warns(UserWarning, match="don't have valid range"):
        check_value_range("unknown_var", xr.DataArray(100))


def test_lazy_evaluation():
    """Test with lazy evaluated DataArray"""
    lazy_array = xr.DataArray(np.array([300])).chunk()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        check_value_range("sst", lazy_array)


def test_array_input():
    """Test with array inputs"""
    # Test with valid array
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        check_value_range("sst", xr.DataArray([280]))

    # Test with array containing invalid values
    with pytest.warns(UserWarning):
        check_value_range("sst", xr.DataArray([350]))


def test_edge_cases():
    """Test edge cases"""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        check_value_range("rad_sw", xr.DataArray(0))
        check_value_range("rad_sw", xr.DataArray(1500))
        check_value_range("slp", xr.DataArray(80000))
        check_value_range("slp", xr.DataArray(110000))


# Parametrized test for more comprehensive coverage
@pytest.mark.parametrize(
    "var_name,value,should_warn",
    [
        ("sst", 270, False),  # exact lower bound
        ("sst", 320, False),  # exact upper bound
        ("sst", 250, True),  # below lower bound
        ("sst", 330, True),  # above upper bound
        ("hum_zt", 0.04, False),  # valid middle value
        ("hum_zt", -0.1, True),  # invalid
        ("u_zu", 0, False),  # valid
        ("unknown", 100, True),  # unknown variable
    ],
)
def test_parametrized(var_name, value, should_warn):
    """Parametrized test covering multiple cases"""
    if should_warn:
        with pytest.warns(UserWarning):
            check_value_range(var_name, xr.DataArray(value))
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            check_value_range(var_name, xr.DataArray(value))
