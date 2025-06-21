"""
pytest for physics.temperature.potential_temperature.py
"""

import pytest

import numpy as np
import xarray as xr
from easyclimate.physics import (
    calc_potential_temperature,
    calc_potential_temperature_vertical,
)
from .const_data import (
    t_data,
    sample_temperature_data_calc_potential_temperature,
    sample_pressure_data_hPa_calc_potential_temperature,
    sample_pressure_data_Pa_calc_potential_temperature,
)
from .util import round_sf_np


# --------------------------------------------
# calc_potential_temperature
# --------------------------------------------
# START ↓
# --------------------------------------------
def test_calc_potential_temperature_hPa(
    sample_temperature_data_calc_potential_temperature,
    sample_pressure_data_hPa_calc_potential_temperature,
):
    """Test potential temperature calculation with hPa units"""
    result = calc_potential_temperature(
        temper_data=sample_temperature_data_calc_potential_temperature,
        pressure_data=sample_pressure_data_hPa_calc_potential_temperature,
        pressure_data_units="hPa",
    )

    # Check the result is a DataArray
    assert isinstance(result, xr.DataArray)

    # Check dimensions are preserved
    assert result.dims == sample_temperature_data_calc_potential_temperature.dims

    # Check coordinates are preserved
    assert all(result.level == sample_temperature_data_calc_potential_temperature.level)

    # Check values (reference values calculated manually)
    expected_values = np.array([283.15, 307.06608908, 369.45667504])
    np.testing.assert_allclose(result.values, expected_values, rtol=1e-3)


def test_calc_potential_temperature_Pa(
    sample_temperature_data_calc_potential_temperature,
    sample_pressure_data_Pa_calc_potential_temperature,
):
    """Test potential temperature calculation with Pa units"""
    result = calc_potential_temperature(
        temper_data=sample_temperature_data_calc_potential_temperature,
        pressure_data=sample_pressure_data_Pa_calc_potential_temperature,
        pressure_data_units="Pa",
    )

    # The result should be identical to the hPa case
    expected_result = calc_potential_temperature(
        temper_data=sample_temperature_data_calc_potential_temperature,
        pressure_data=sample_pressure_data_Pa_calc_potential_temperature
        / 100,  # Convert to hPa
        pressure_data_units="hPa",
    )

    xr.testing.assert_allclose(result, expected_result)


def test_calc_potential_temperature_custom_kappa(
    sample_temperature_data_calc_potential_temperature,
    sample_pressure_data_hPa_calc_potential_temperature,
):
    """Test potential temperature calculation with custom kappa value"""
    custom_kappa = 0.285  # Different kappa value
    result = calc_potential_temperature(
        temper_data=sample_temperature_data_calc_potential_temperature,
        pressure_data=sample_pressure_data_hPa_calc_potential_temperature,
        pressure_data_units="hPa",
        kappa=custom_kappa,
    )

    # Calculate expected values with custom kappa
    expected_values = (
        sample_temperature_data_calc_potential_temperature.values
        * (1000.0 / sample_pressure_data_hPa_calc_potential_temperature.values)
        ** custom_kappa
    )
    np.testing.assert_allclose(result.values, expected_values)


def test_calc_potential_temperature_edge_cases():
    """Test edge cases like very low/high pressures"""
    # Very low pressure case
    temp_data = xr.DataArray([300.0])
    press_data = xr.DataArray([1.0])  # 1 hPa

    result = calc_potential_temperature(
        temper_data=temp_data, pressure_data=press_data, pressure_data_units="hPa"
    )

    expected_value = 300.0 * (1000.0 / 1.0) ** (287 / 1005.7)
    np.testing.assert_allclose(result.values[0], expected_value)

    # Reference pressure case (should return input temperature)
    press_data = xr.DataArray([1000.0])
    result = calc_potential_temperature(
        temper_data=temp_data, pressure_data=press_data, pressure_data_units="hPa"
    )
    np.testing.assert_allclose(result.values[0], 300.0)


# --------------------------------------------
# calc_potential_temperature
# --------------------------------------------
# END ↑
# --------------------------------------------


def test_calc_potential_temperature_vertical(t_data):
    pv = calc_potential_temperature_vertical(
        t_data, vertical_dim="level", vertical_dim_units="hPa"
    )

    result_data = round_sf_np(pv).data
    refer_data = np.array(
        [
            299.8,
            302.0,
            305.0,
            313.0,
            319.1,
            326.6,
            335.2,
            343.2,
            345.6,
            348.7,
            352.4,
            371.2,
            418.6,
            477.3,
            577.5,
            673.4,
            863.8,
        ]
    )
    assert np.isclose(result_data, refer_data).all()
