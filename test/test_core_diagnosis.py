"""
pytest for diagnosis.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from .util import round_sf_np

# Sample Data Declaration
t_data = xr.DataArray(
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

pv_data = xr.DataArray(
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

z_data = xr.DataArray(
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

t_data_Tv = xr.DataArray(
    np.array([299.80002, 295.4]),
    dims="level",
    coords={"level": np.array([1000.0, 925.0])},
)

q_data_Tv = xr.DataArray(
    np.array([0.01864, 0.0149525]),
    dims="level",
    coords={"level": np.array([1000.0, 925.0])},
)


def test_calc_brunt_vaisala_frequency_atm():
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


def test_calc_potential_temperature():
    pv = ecl.calc_potential_temperature(
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


def test_calc_static_stability():
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


def test_calc_virtual_temperature():
    result_data = ecl.calc_virtual_temperature(
        temper_data=t_data_Tv,
        specific_humidity_data=q_data_Tv,
        specific_humidity_units="g/g",
    ).data
    refer_data = np.array([303.1976896, 298.08551685])
    assert np.isclose(result_data, refer_data).all()


def test_calc_virtual_temperature_Hobbs2006():
    result_data = ecl.calc_virtual_temperature_Hobbs2006(
        temper_data=t_data_Tv,
        specific_humidity_data=q_data_Tv,
        specific_humidity_units="g/g",
    ).data
    refer_data = np.array([303.1345747, 298.04520656])
    assert np.isclose(result_data, refer_data).all()
