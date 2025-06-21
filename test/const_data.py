"""
Sample Data Declaration
"""

import pytest
import numpy as np
import xarray as xr


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
