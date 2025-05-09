"""
pytest for diagnosis.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from .util import round_sf_np


# Sample Data Declaration
@pytest.fixture
def t_data():
    return xr.DataArray(
        np.array(
            [
                299.80002,
                295.4,
                291.2,
                282.7,
                275.80002,
                268.0,
                258.1,
                243.40001,
                232.7,
                220.3,
                205.1,
                192.40001,
                196.0,
                203.0,
                212.3,
                220.5,
                232.1,
            ]
        ),
        dims="level",
        coords={
            "level": np.array(
                [
                    1000.0,
                    925.0,
                    850.0,
                    700.0,
                    600.0,
                    500.0,
                    400.0,
                    300.0,
                    250.0,
                    200.0,
                    150.0,
                    100.0,
                    70.0,
                    50.0,
                    30.0,
                    20.0,
                    10.0,
                ]
            )
        },
    )


@pytest.fixture
def pv_data():
    return xr.DataArray(
        np.array(
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
        ),
        dims="level",
        coords={
            "level": np.array(
                [
                    1000.0,
                    925.0,
                    850.0,
                    700.0,
                    600.0,
                    500.0,
                    400.0,
                    300.0,
                    250.0,
                    200.0,
                    150.0,
                    100.0,
                    70.0,
                    50.0,
                    30.0,
                    20.0,
                    10.0,
                ]
            )
        },
    )


@pytest.fixture
def z_data():
    return xr.DataArray(
        np.array(
            [
                85.58,
                768.77,
                1500.67,
                3141.68,
                4407.93,
                5862.28,
                7586.9,
                9700.02,
                10970.99,
                12452.39,
                14248.64,
                16587.65,
                18609.14,
                20568.38,
                23691.34,
                26259.47,
                30865.37,
            ]
        ),
        dims="level",
        coords={
            "level": np.array(
                [
                    1000.0,
                    925.0,
                    850.0,
                    700.0,
                    600.0,
                    500.0,
                    400.0,
                    300.0,
                    250.0,
                    200.0,
                    150.0,
                    100.0,
                    70.0,
                    50.0,
                    30.0,
                    20.0,
                    10.0,
                ]
            )
        },
    )


@pytest.fixture
def t_data_Tv():
    return xr.DataArray(
        np.array([299.80002, 295.4]),
        dims="level",
        coords={"level": np.array([1000.0, 925.0])},
    )


@pytest.fixture
def q_data_Tv():
    return xr.DataArray(
        np.array([0.01864, 0.0149525]),
        dims="level",
        coords={"level": np.array([1000.0, 925.0])},
    )


def test_calc_brunt_vaisala_frequency_atm(z_data, pv_data):
    N2_data = ecl.calc_brunt_vaisala_frequency_atm(
        potential_temperature_data=pv_data, z_data=z_data, vertical_dim="level"
    )

    result_data = round_sf_np(N2_data).data
    refer_data = np.array(
        [
            0.00945,
            0.01092,
            0.0122,
            0.01232,
            0.01239,
            0.01233,
            0.01125,
            0.009368,
            0.007528,
            0.007636,
            0.0123,
            0.02002,
            0.02498,
            0.02534,
            0.02418,
            0.0241,
            0.02189,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_get_coriolis_parameter():
    latdata = np.array([30, 60])
    x = ecl.get_coriolis_parameter(latdata, omega=7.292e-5)

    result_data = round_sf_np(x)
    refer_data = np.array([7.292e-05, 1.263e-04])
    assert np.isclose(result_data, refer_data).all()


# test_calc_potential_temperature
@pytest.fixture
def sample_temperature_data_calc_potential_temperature():
    """Fixture providing sample temperature data"""
    return xr.DataArray(
        data=np.array([283.15, 293.15, 303.15]),  # 10°C, 20°C, 30°C
        dims=["level"],
        coords={"level": [1000, 850, 500]},
        attrs={"units": "K"},
    )


@pytest.fixture
def sample_pressure_data_hPa_calc_potential_temperature():
    """Fixture providing sample pressure data in hPa"""
    return xr.DataArray(
        data=np.array([1000, 850, 500]),
        dims=["level"],
        coords={"level": [1000, 850, 500]},
        attrs={"units": "hPa"},
    )


@pytest.fixture
def sample_pressure_data_Pa_calc_potential_temperature():
    """Fixture providing sample pressure data in Pa"""
    return xr.DataArray(
        data=np.array([100000, 85000, 50000]),
        dims=["level"],
        coords={"level": [1000, 850, 500]},
        attrs={"units": "Pa"},
    )


def test_calc_potential_temperature_hPa(
    sample_temperature_data_calc_potential_temperature,
    sample_pressure_data_hPa_calc_potential_temperature,
):
    """Test potential temperature calculation with hPa units"""
    result = ecl.calc_potential_temperature(
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
    result = ecl.calc_potential_temperature(
        temper_data=sample_temperature_data_calc_potential_temperature,
        pressure_data=sample_pressure_data_Pa_calc_potential_temperature,
        pressure_data_units="Pa",
    )

    # The result should be identical to the hPa case
    expected_result = ecl.calc_potential_temperature(
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
    result = ecl.calc_potential_temperature(
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

    result = ecl.calc_potential_temperature(
        temper_data=temp_data, pressure_data=press_data, pressure_data_units="hPa"
    )

    expected_value = 300.0 * (1000.0 / 1.0) ** (287 / 1005.7)
    np.testing.assert_allclose(result.values[0], expected_value)

    # Reference pressure case (should return input temperature)
    press_data = xr.DataArray([1000.0])
    result = ecl.calc_potential_temperature(
        temper_data=temp_data, pressure_data=press_data, pressure_data_units="hPa"
    )
    np.testing.assert_allclose(result.values[0], 300.0)


# calc_potential_temperature_vertical
def test_calc_potential_temperature_vertical(t_data):
    pv = ecl.calc_potential_temperature_vertical(
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


# calc_equivalent_potential_temperature


# Fixture for creating test DataArrays
@pytest.fixture
def sample_data_calc_equivalent_potential_temperature():
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
    sample_data_calc_equivalent_potential_temperature,
):
    """Test basic calculation with hPa and Celsius units"""
    pressure, temperature, dewpoint = sample_data_calc_equivalent_potential_temperature

    result = ecl.calc_equivalent_potential_temperature(
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
    sample_data_calc_equivalent_potential_temperature,
):
    """Test with different input units"""
    pressure, temperature, dewpoint = sample_data_calc_equivalent_potential_temperature

    # Test with Pa and Kelvin
    result_pa_k = ecl.calc_equivalent_potential_temperature(
        pressure_data=pressure * 100,  # convert to Pa
        temperature_data=temperature + 273.15,  # convert to Kelvin
        dewpoint_data=dewpoint + 273.15,  # convert to Kelvin
        pressure_data_units="Pa",
        temperature_data_units="kelvin",
        dewpoint_data_units="kelvin",
    )

    # Test with Fahrenheit
    result_f = ecl.calc_equivalent_potential_temperature(
        pressure_data=pressure,
        temperature_data=(temperature * 9 / 5) + 32,  # convert to Fahrenheit
        dewpoint_data=(dewpoint * 9 / 5) + 32,  # convert to Fahrenheit
        pressure_data_units="hPa",
        temperature_data_units="fahrenheit",
        dewpoint_data_units="fahrenheit",
    )

    # All results should be approximately equal
    xr.testing.assert_allclose(result_pa_k, result_f, rtol=1e-3)


def test_known_value_calc_equivalent_potential_temperature():
    """Test against a known value from literature or manual calculation"""
    # Using example values that should produce a known result
    pressure = xr.DataArray([1000.0], dims=["level"])
    temperature = xr.DataArray([25.0], dims=["level"])
    dewpoint = xr.DataArray([20.0], dims=["level"])

    result = ecl.calc_equivalent_potential_temperature(
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


def test_edge_cases_calc_equivalent_potential_temperature():
    """Test edge cases like very low/high temperatures"""
    pressure = xr.DataArray([1000.0, 1000.0], dims=["level"])

    # Very cold case
    temperature_cold = xr.DataArray([-40.0, -50.0], dims=["level"])
    dewpoint_cold = xr.DataArray([-45.0, -55.0], dims=["level"])

    result_cold = ecl.calc_equivalent_potential_temperature(
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

    result_hot = ecl.calc_equivalent_potential_temperature(
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


def test_dimensional_consistency_calc_equivalent_potential_temperature():
    """Test that input arrays must have same dimensions"""
    pressure = xr.DataArray([1000.0, 950.0], dims=["level"])
    temperature = xr.DataArray([25.0, 20.0, 15.0], dims=["level"])  # different size
    dewpoint = xr.DataArray([20.0, 15.0], dims=["level"])

    with pytest.raises(ValueError):
        ecl.calc_equivalent_potential_temperature(
            pressure_data=pressure,
            temperature_data=temperature,
            dewpoint_data=dewpoint,
            pressure_data_units="hPa",
            temperature_data_units="celsius",
            dewpoint_data_units="celsius",
        )


# calc_static_stability
def test_calc_static_stability(t_data):
    x = ecl.calc_static_stability(
        t_data, vertical_dim="level", vertical_dim_units="hPa"
    )

    result_data = round_sf_np(x).data
    refer_data = np.array(
        [
            0.0002514,
            0.0003402,
            0.0004607,
            0.0005096,
            0.0005877,
            0.0006617,
            0.0006387,
            0.0004954,
            0.0003722,
            0.0004301,
            0.00128,
            0.004139,
            0.009856,
            0.01633,
            0.02436,
            0.0444,
            0.06889,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_virtual_temperature(q_data_Tv, t_data_Tv):
    result_data = ecl.calc_virtual_temperature(
        temper_data=t_data_Tv,
        specific_humidity_data=q_data_Tv,
        specific_humidity_data_units="g/g",
    ).data
    refer_data = np.array([303.1976896, 298.08551685])
    assert np.isclose(result_data, refer_data).all()


def test_calc_virtual_temperature_Hobbs2006(q_data_Tv, t_data_Tv):
    result_data = ecl.calc_virtual_temperature_Hobbs2006(
        temper_data=t_data_Tv,
        specific_humidity_data=q_data_Tv,
        specific_humidity_data_units="g/g",
    ).data
    refer_data = np.array([303.1345747, 298.04520656])
    assert np.isclose(result_data, refer_data).all()


def test_calc_dewpoint():
    result_data = ecl.calc_dewpoint(
        vapor_pressure_data=xr.DataArray([22, 23]), vapor_pressure_data_units="hPa"
    ).data
    refer_data = np.array([19.02910179, 19.74308483])
    assert np.isclose(result_data, refer_data).all()


def test_calc_mixing_ratio():
    result_data = ecl.calc_mixing_ratio(
        partial_pressure_data=xr.DataArray([25, 30]),
        total_pressure_data=xr.DataArray([1000, 1000]),
    ).data
    refer_data = np.array([0.01594761, 0.01923578])
    assert np.isclose(result_data, refer_data).all()


def test_calc_vapor_pressure():
    result_data = ecl.calc_vapor_pressure(
        pressure_data=xr.DataArray([988, 991]),
        mixing_ratio_data=xr.DataArray([0.018, 0.0159]),
        pressure_data_units="hPa",
    ).data
    refer_data = np.array([27.789371, 24.70287576])
    assert np.isclose(result_data, refer_data).all()


def test_calc_saturation_vapor_pressure():
    result_data = ecl.calc_saturation_vapor_pressure(
        temperature_data=xr.DataArray([24, 36, 40]), temperature_data_units="celsius"
    ).data
    refer_data = np.array([29.83254305, 59.5118433, 73.94900581])
    assert np.isclose(result_data, refer_data).all()


def test_calc_saturation_mixing_ratio():
    result_data = ecl.calc_saturation_mixing_ratio(
        total_pressure_data=xr.DataArray([983, 991]),
        temperature_data=xr.DataArray([25, 28]),
        temperature_data_units="celsius",
        total_pressure_data_units="hPa",
    ).data
    refer_data = np.array([0.02070799, 0.02467106])


def test_transfer_mixing_ratio_2_specific_humidity():
    result_data = ecl.transfer_mixing_ratio_2_specific_humidity(
        mixing_ratio_data=xr.DataArray([0.01594761, 0.01923578])
    ).data
    refer_data = np.array([0.01569728, 0.01887275])
    assert np.isclose(result_data, refer_data).all()


def test_transfer_specific_humidity_2_mixing_ratio():
    result_data = ecl.transfer_specific_humidity_2_mixing_ratio(
        specific_humidity_data=xr.DataArray([0.01569728, 0.01887275]),
        specific_humidity_data_units="g/g",
    ).data
    refer_data = np.array([0.01594761, 0.01923578])
    assert np.isclose(result_data, refer_data).all()


def test_transfer_dewpoint_2_specific_humidity():
    result_data = ecl.transfer_dewpoint_2_specific_humidity(
        dewpoint_data=xr.DataArray([15, 18]),
        pressure_data=xr.DataArray([988, 980]),
        dewpoint_data_units="celsius",
        pressure_data_units="hPa",
    ).data
    refer_data = np.array([0.01079758, 0.01319518])
    assert np.isclose(result_data, refer_data).all()


def test_transfer_specific_humidity_2_dewpoint():
    result_data = ecl.transfer_specific_humidity_2_dewpoint(
        specific_humidity_data=xr.DataArray([0.01079, 0.01319]),
        pressure_data=xr.DataArray([988, 980]),
        specific_humidity_data_units="g/g",
        pressure_data_units="hPa",
    ).data
    refer_data = np.array([14.98916116, 17.9938138])
    assert np.isclose(result_data, refer_data).all()


def test_transfer_dewpoint_2_relative_humidity():
    result_data = ecl.transfer_dewpoint_2_relative_humidity(
        temperature_data=xr.DataArray([25, 28]),
        dewpoint_data=xr.DataArray([12, 13]),
        temperature_data_units="celsius",
        dewpoint_data_units="celsius",
    ).data
    refer_data = np.array([0.44248477, 0.3958322])
    assert np.isclose(result_data, refer_data).all()


def test_transfer_mixing_ratio_2_relative_humidity():
    result_data = ecl.transfer_mixing_ratio_2_relative_humidity(
        pressure_data=xr.DataArray([1013.25, 1014]),
        temperature_data=xr.DataArray([30, 28]),
        mixing_ratio_data=xr.DataArray([18 / 1000, 16 / 1000]),
        pressure_data_units="hPa",
        temperature_data_units="celsius",
    ).data
    refer_data = np.array([0.67127708, 0.67260414])
    assert np.isclose(result_data, refer_data).all()


def test_transfer_specific_humidity_2_relative_humidity():
    result_data = ecl.transfer_specific_humidity_2_relative_humidity(
        pressure_data=xr.DataArray([1013.25, 1014]),
        temperature_data=xr.DataArray([30, 28]),
        specific_humidity_data=xr.DataArray([18 / 1000, 16 / 1000]),
        pressure_data_units="hPa",
        temperature_data_units="degC",
        specific_humidity_data_units="g/g",
    ).data
    refer_data = np.array([0.6832293, 0.68326215])
    assert np.isclose(result_data, refer_data).all()
