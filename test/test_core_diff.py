"""
pytest for diff.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
import pandas as pd

# Sample Data Declaration
t_data = xr.DataArray(
    np.array(
        [[268.03, 268.13, 268.19], [268.28, 268.33, 268.28], [268.46, 268.42, 268.25]]
    ),
    dims=("lat", "lon"),
    coords={"lat": np.array([5.0, 2.5, 0.0]), "lon": np.array([0.0, 2.5, 5.0])},
)

t_level_data = xr.DataArray(
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

t_time_data = xr.DataArray(
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
    dims="time",
    coords={"time": pd.date_range("2023-01-01", freq="D", periods=17)},
)

t_data_delta_pressure = xr.DataArray(
    np.array(
        [
            [299.9145, 299.78305],
            [295.8137, 295.64435],
            [291.14514, 291.22095],
            [282.74918, 282.86853],
            [275.30966, 275.2653],
            [268.46048, 268.41855],
            [257.65402, 257.7661],
            [242.13304, 242.2766],
            [232.0645, 232.07822],
            [220.43225, 220.11288],
            [206.3774, 206.2516],
            [193.91046, 194.45079],
            [196.07178, 195.95001],
            [204.80402, 204.64435],
            [209.65886, 209.90562],
            [220.95563, 221.15643],
            [226.5129, 226.2258],
        ]
    ),
    dims=("level", "lon"),
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
        ),
        "lon": np.array([0.0, 2.5]),
    },
)

msl_data_delta_pressure = xr.DataArray(
    np.array(np.array([1012.017, 1011.9201])),
    dims=("lon"),
    coords={"lon": np.array([0.0, 2.5])},
)

u_data_500hpa = xr.DataArray(
    np.array([[-3.8, -3.99, -3.77], [-6.51, -6.6, -6.33], [-7.78, -7.67, -7.41]]),
    dims=("lat", "lon"),
    coords={"lat": np.array([5.0, 2.5, 0.0]), "lon": np.array([0.0, 2.5, 5.0])},
)

v_data_500hpa = xr.DataArray(
    np.array([[0.13, -0.1, -0.18], [-0.67, -0.79, -0.62], [-0.76, -0.99, -0.87]]),
    dims=("lat", "lon"),
    coords={"lat": np.array([5.0, 2.5, 0.0]), "lon": np.array([0.0, 2.5, 5.0])},
)

z_data_500hpa = xr.DataArray(
    np.array(
        [
            [5868.7583, 5869.121, 5869.4756],
            [5868.121, 5868.113, 5868.2095],
            [5867.411, 5867.4355, 5867.5806],
        ]
    ),
    dims=("lat", "lon"),
    coords={"lat": np.array([5.0, 2.5, 0.0]), "lon": np.array([0.0, 2.5, 5.0])},
)

q_data_500hpa = xr.DataArray(
    np.array(
        [
            [1.24, 1.2658951, 1.5115161],
            [1.326887, 1.5464114, 1.8371367],
            [1.5063628, 1.7480158, 2.0217905],
        ]
    ),
    dims=("lat", "lon"),
    coords={"lat": np.array([5.0, 2.5, 0.0]), "lon": np.array([0.0, 2.5, 5.0])},
)

omega_data_sample1 = xr.DataArray(
    np.array(
        [
            0.02437097,
            0.03512904,
            0.04098387,
            0.02112902,
            0.01658062,
            0.02814513,
            0.03577418,
            0.04062096,
            0.03781452,
            0.02902419,
            0.00769354,
            -0.01137904,
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
            ]
        )
    },
)

q_data_sample1 = xr.DataArray(
    np.array(
        [
            17.74024,
            13.399196,
            9.72387,
            5.151291,
            3.0942986,
            1.7616936,
            1.09875,
            0.32160482,
        ]
    ),
    dims="level",
    coords={
        "level": np.array([1000.0, 925.0, 850.0, 700.0, 600.0, 500.0, 400.0, 300.0])
    },
)

t_data_sample1 = xr.DataArray(
    np.array(
        [
            299.68384,
            294.31772,
            290.3024,
            283.88385,
            277.28384,
            269.40726,
            259.0016,
            245.04675,
            234.3129,
            221.59273,
            206.91853,
            191.03226,
            194.09918,
            203.64514,
            211.79112,
            221.36288,
            229.66693,
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

# -----------------------------------------------------------
# TEST


def test_calc_gradient():
    result_data = ecl.calc_gradient(t_data, dim="lon").data
    refer_data = np.array(
        [[0.12, 0.08, 0.04], [0.1, 0.0, -0.1], [0.025, -0.105, -0.235]]
    )
    assert np.isclose(result_data.data, refer_data).all()


def test_calc_lon_gradient():
    result_data = ecl.calc_lon_gradient(t_data).data
    refer_data = np.array(
        [
            [4.33391322e-07, 2.88927548e-07, 1.44463774e-07],
            [3.60127877e-07, 0.00000000e00, -3.60127877e-07],
            [8.99462787e-08, -3.77774370e-07, -8.45495020e-07],
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_lat_gradient():
    result_data = ecl.calc_lat_gradient(t_data).data
    refer_data = np.array(
        [
            [-1.02538758e-06, -9.17452042e-07, -5.39677672e-07],
            [-7.73537997e-07, -5.21688416e-07, -1.07935534e-07],
            [-5.21688416e-07, -1.25924790e-07, 3.23806603e-07],
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_lon_laplacian():
    result_data = ecl.calc_lon_laplacian(t_data).data
    refer_data = np.array(
        [
            [-5.21744551e-13, -5.21744551e-13, -5.21744551e-13],
            [-1.29692088e-12, -1.29692088e-12, -1.29692088e-12],
            [-1.68278927e-12, -1.68278927e-12, -1.68278927e-12],
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_lat_laplacian():
    result_data = ecl.calc_lat_laplacian(t_data).data
    refer_data = np.array(
        [
            [-9.06117301e-13, -1.42389862e-12, -1.55334394e-12],
            [-9.06117301e-13, -1.42389862e-12, -1.55334394e-12],
            [-9.06117301e-13, -1.42389862e-12, -1.55334394e-12],
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_lon_lat_mixed_derivatives():
    result_data = ecl.calc_lon_lat_mixed_derivatives(t_data).data
    refer_data = np.array(
        [
            [-9.74548417e-14, 8.77093575e-13, 1.85164199e-12],
            [6.15451085e-13, 1.19851001e-12, 1.78156893e-12],
            [1.32681462e-12, 1.52098261e-12, 1.71515061e-12],
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_p_gradient():
    result_data = ecl.calc_p_gradient(
        t_level_data, vertical_dim="level", vertical_dim_units="hPa"
    ).data
    refer_data = np.array(
        [
            0.0006,
            0.00057333,
            0.00056444,
            0.000616,
            0.000735,
            0.000885,
            0.00123,
            0.00169333,
            0.00231,
            0.00276,
            0.00279,
            0.0011375,
            -0.00212,
            -0.004075,
            -0.00583333,
            -0.0099,
            -0.0133,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_time_gradient():
    result_data = ecl.calc_time_gradient(t_time_data, time_units="day").data
    refer_data = np.array(
        [
            -5.20836806e-05,
            -4.97686343e-05,
            -7.34953704e-05,
            -8.91202546e-05,
            -8.50694444e-05,
            -1.02430671e-04,
            -1.42361053e-04,
            -1.46990741e-04,
            -1.33680613e-04,
            -1.59722222e-04,
            -1.61458275e-04,
            -5.26620370e-05,
            6.13425347e-05,
            9.43287037e-05,
            1.01273148e-04,
            1.14583333e-04,
            1.53935185e-04,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_delta_pressure():
    result_data = ecl.calc_delta_pressure(
        data_input=t_data_delta_pressure,
        surface_pressure_data=msl_data_delta_pressure,
        vertical_dim="level",
        vertical_dim_units="hPa",
        surface_pressure_data_units="Pa",
    ).data.flatten()

    refer_data = np.array(
        [
            -95237.983,
            -95238.0799,
            7500.0,
            7500.0,
            11250.0,
            11250.0,
            12500.0,
            12500.0,
            10000.0,
            10000.0,
            10000.0,
            10000.0,
            10000.0,
            10000.0,
            7500.0,
            7500.0,
            5000.0,
            5000.0,
            5000.0,
            5000.0,
            5000.0,
            5000.0,
            4000.0,
            4000.0,
            2500.0,
            2500.0,
            2000.0,
            2000.0,
            1500.0,
            1500.0,
            1000.0,
            1000.0,
            500.0,
            500.0,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_p_integral():
    result_data = ecl.calc_p_integral(t_level_data, vertical_dim="level").data
    refer_data = 258.63737768
    assert np.isclose(result_data, refer_data).all()


def test_calc_top2surface_integral():
    result_data = ecl.calc_top2surface_integral(
        data_input=t_data_delta_pressure,
        surface_pressure_data=msl_data_delta_pressure,
        vertical_dim="level",
        vertical_dim_units="hPa",
        surface_pressure_data_units="Pa",
    ).data.flatten()
    refer_data = np.array([224.38077318, 224.09617167])
    assert np.isclose(result_data, refer_data).all()


def test_calc_laplacian():
    result_data = ecl.calc_laplacian(t_data).data.flatten()
    refer_data = np.array([-3.23613322e-12])
    assert np.isclose(result_data, refer_data).all()


def test_calc_divergence():
    result_data = ecl.calc_divergence(u_data_500hpa, v_data_500hpa).data.flatten()
    refer_data = np.array(
        [
            2.72715282e-06,
            3.41953819e-06,
            3.46225017e-06,
            6.33290775e-07,
            1.93057363e-06,
            2.86608366e-06,
            -8.27505764e-07,
            5.03699161e-07,
            1.76294706e-06,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_vorticity():
    result_data = ecl.calc_vorticity(u_data_500hpa, v_data_500hpa).data.flatten()
    refer_data = np.array(
        [
            -1.34943568e-05,
            -1.27753346e-05,
            -1.19427458e-05,
            -8.15868318e-06,
            -6.57525154e-06,
            -5.45707305e-06,
            -3.43594785e-06,
            -1.27723716e-06,
            -1.61903302e-07,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_geostrophic_wind():
    result_data = ecl.calc_geostrophic_wind(z_data_500hpa)
    result_data1 = result_data["ug"].data.flatten()[:3]
    result_data2 = result_data["vg"].data.flatten()[:3]
    refer_data1 = np.array([-1.66699748, -3.25452166, -4.39585806])
    refer_data2 = np.array([1.0212275, 0.99867278, 0.97611806])
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()


def test_calc_horizontal_water_flux():
    result_data = ecl.calc_horizontal_water_flux(
        specific_humidity_data=q_data_500hpa, u_data=u_data_500hpa, v_data=v_data_500hpa
    )
    result_data1 = result_data["qu"].data
    result_data2 = result_data["qv"].data
    refer_data1 = np.array(
        [
            [-0.48081633, -0.51540015, -0.58147099],
            [-0.88143208, -1.04146074, -1.18664034],
            [-1.19586761, -1.36808992, -1.52872118],
        ]
    )
    refer_data2 = np.array(
        [
            [0.01644898, -0.0129173, -0.02776254],
            [-0.09071574, -0.12465969, -0.11622702],
            [-0.11681997, -0.17658527, -0.17948548],
        ]
    )
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()


def test_calc_vertical_water_flux():
    result_data = ecl.calc_vertical_water_flux(q_data_sample1, omega_data_sample1).data
    refer_data = np.array(
        [
            -0.04411703,
            -0.0480307,
            -0.04066549,
            -0.0111063,
            -0.00523524,
            -0.0050595,
            -0.00401091,
            -0.00133305,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_water_flux_top2surface_integral():
    specific_humidity_data_sample = (
        ecl.open_tutorial_dataset("shum_202201_mon_mean")
        .shum.isel(time=0)
        .sel(lat=slice(5, 0), lon=slice(0, 5))
    )
    u_data_sample = (
        ecl.open_tutorial_dataset("uwnd_202201_mon_mean")
        .uwnd.isel(time=0)
        .sel(lat=slice(5, 0), lon=slice(0, 5))
    )
    v_data_sample = (
        ecl.open_tutorial_dataset("vwnd_202201_mon_mean")
        .vwnd.isel(time=0)
        .sel(lat=slice(5, 0), lon=slice(0, 5))
    )
    surface_pressure_data_sample = (
        ecl.open_tutorial_dataset("pressfc_202201_mon_mean")
        .pres.isel(time=0)
        .sel(lat=slice(5, 0), lon=slice(0, 5))
    )

    result_data = ecl.calc_water_flux_top2surface_integral(
        specific_humidity_data_sample,
        u_data_sample,
        v_data_sample,
        surface_pressure_data_sample,
        surface_pressure_data_units="hPa",
        vertical_dim="level",
        vertical_dim_units="hPa",
    )

    result_data1 = result_data["qu"].data.flatten()
    result_data2 = result_data["qv"].data.flatten()
    refer_data1 = np.array(
        [
            -61719.65,
            -53384.867,
            -46633.44,
            -94774.73,
            -84097.77,
            -68822.93,
            -127009.766,
            -119459.65,
            -109226.625,
        ]
    )
    refer_data2 = np.array(
        [
            1478.0461,
            4159.2134,
            -2113.73,
            19404.598,
            21102.201,
            16411.918,
            31748.395,
            31434.625,
            23815.54,
        ]
    )

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()


def test_calc_divergence_watervaporflux():
    result_data = ecl.calc_divergence_watervaporflux(
        q_data_500hpa,
        u_data_500hpa,
        v_data_500hpa,
        specific_humidity_units="g/kg",
    ).data.flatten()
    refer_data = np.array(
        [
            4.63115041e-10,
            3.28054031e-10,
            6.85267982e-11,
            -3.62682915e-10,
            -2.54289074e-10,
            -2.22359484e-10,
            -6.92384309e-10,
            -5.19564327e-10,
            -3.74825155e-10,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_divergence_watervaporflux_top2surface_integral():
    specific_humidity_data_sample = (
        ecl.open_tutorial_dataset("shum_202201_mon_mean")
        .shum.isel(time=0)
        .sel(lat=slice(5, 0), lon=slice(0, 5))
    )
    u_data_sample = (
        ecl.open_tutorial_dataset("uwnd_202201_mon_mean")
        .uwnd.isel(time=0)
        .sel(lat=slice(5, 0), lon=slice(0, 5))
    )
    v_data_sample = (
        ecl.open_tutorial_dataset("vwnd_202201_mon_mean")
        .vwnd.isel(time=0)
        .sel(lat=slice(5, 0), lon=slice(0, 5))
    )
    surface_pressure_data_sample = (
        ecl.open_tutorial_dataset("pressfc_202201_mon_mean")
        .pres.isel(time=0)
        .sel(lat=slice(5, 0), lon=slice(0, 5))
    )

    result_data = ecl.calc_divergence_watervaporflux_top2surface_integral(
        specific_humidity_data_sample,
        u_data_sample,
        v_data_sample,
        surface_pressure_data_sample,
        vertical_dim="level",
        specific_humidity_units="g/kg",
        surface_pressure_data_units="hPa",
        vertical_dim_units="hPa",
    ).data.flatten()
    refer_data = np.array(
        [
            -3.16235016e-08,
            -5.25733773e-08,
            -6.11471233e-08,
            -1.84145484e-08,
            3.67623974e-10,
            1.46561626e-08,
            -1.93079642e-08,
            1.67279636e-08,
            3.28154620e-08,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_u_advection():
    result_data = ecl.calc_u_advection(u_data_500hpa, t_data).data.flatten()
    refer_data = np.array(
        [
            1.64688703e-06,
            1.15282092e-06,
            5.44628429e-07,
            2.34443248e-06,
            0.00000000e00,
            -2.27960946e-06,
            6.99782048e-07,
            -2.89752942e-06,
            -6.26511809e-06,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_v_advection():
    result_data = ecl.calc_v_advection(v_data_500hpa, t_data).data.flatten()
    refer_data = np.array(
        [
            1.33300385e-07,
            -9.17452042e-08,
            -9.71419810e-08,
            -5.18270458e-07,
            -4.12133849e-07,
            -6.69200313e-08,
            -3.96483196e-07,
            -1.24665542e-07,
            2.81711745e-07,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_p_advection():
    result_data = ecl.calc_p_advection(
        omega_data_sample1,
        t_data_sample1,
        vertical_dim="level",
        vertical_dim_units="hPa",
    )
    refer_data = np.array(
        [
            -1.96316937e-05,
            -2.19707321e-05,
            -1.90053499e-05,
            -1.10027766e-05,
            -1.20015419e-05,
            -2.57278011e-05,
            -4.35738635e-05,
            -6.68585797e-05,
            -8.86902508e-05,
            -7.95099400e-05,
            -2.35118198e-05,
            1.82339871e-05,
        ]
    )
    assert np.isclose(result_data, refer_data).all()
