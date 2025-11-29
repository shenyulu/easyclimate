"""
pytest for eddy.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from easyclimate.core.eddy import (
    calc_Plumb_wave_activity_horizontal_flux,
    calc_TN_wave_activity_horizontal_flux,
)

time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)

uwnd_daily = (
    ecl.open_tutorial_dataset("uwnd_2022_day5")
    .uwnd.sortby("lat")
    .sel(lon=slice(100, 110), lat=slice(20, 30))
)
z_daily = (
    ecl.open_tutorial_dataset("hgt_2022_day5")
    .hgt.sortby("lat")
    .sel(lon=slice(100, 110), lat=slice(20, 30))
)
temp_daily = (
    ecl.open_tutorial_dataset("air_2022_day5")
    .air.sortby("lat")
    .sel(lon=slice(100, 110), lat=slice(20, 30))
)
vwnd_daily = (
    ecl.open_tutorial_dataset("vwnd_2022_day5")
    .vwnd.sortby("lat")
    .sel(lon=slice(100, 110), lat=slice(20, 30))
)
omega_daily = (
    ecl.open_tutorial_dataset("omega_2022_day5")
    .omega.sortby("lat")
    .sel(lon=slice(100, 110), lat=slice(20, 30))
)
q_daily = (
    ecl.open_tutorial_dataset("shum_2022_day5")
    .shum.sortby("lat")
    .sel(lon=slice(100, 110), lat=slice(20, 30))
)

z_climate_data = (
    ecl.open_tutorial_dataset(
        "hgt_day_ltm_1991_2020_0to6day.nc", decode_times=time_coder
    )
    .hgt.sortby("lat")
    .sel(lon=slice(100, 110), lat=slice(20, 30))
)
u_climate_data = (
    ecl.open_tutorial_dataset(
        "uwnd_day_ltm_1991_2020_0to6day.nc", decode_times=time_coder
    )
    .uwnd.sortby("lat")
    .sel(lon=slice(100, 110), lat=slice(20, 30))
)
v_climate_data = (
    ecl.open_tutorial_dataset(
        "vwnd_day_ltm_1991_2020_0to6day.nc", decode_times=time_coder
    )
    .vwnd.sortby("lat")
    .sel(lon=slice(100, 110), lat=slice(20, 30))
)

z_prime_data = ecl.calc_daily_climatological_anomaly(
    data_daily=z_daily, data_climatology_daily_data=z_climate_data
)
u_prime_data = ecl.calc_daily_climatological_anomaly(
    data_daily=uwnd_daily, data_climatology_daily_data=u_climate_data
)
v_prime_data = ecl.calc_daily_climatological_anomaly(
    data_daily=vwnd_daily, data_climatology_daily_data=v_climate_data
)

u_climate_data_need = ecl.populate_daymean2everyday(
    data_daily=z_prime_data, data_climatology_daily_data=u_climate_data
)
v_climate_data_need = ecl.populate_daymean2everyday(
    data_daily=z_prime_data, data_climatology_daily_data=v_climate_data
)


def test_calc_eady_growth_rate():
    result_data = ecl.calc_eady_growth_rate(
        u_daily_data=uwnd_daily.isel(time=0),
        z_daily_data=z_daily.isel(time=0),
        temper_daily_data=temp_daily.isel(time=0),
        vertical_dim="level",
        vertical_dim_units="hPa",
    )
    result_data1 = result_data["eady_growth_rate"].sel(level=500).data.flatten()
    result_data2 = result_data["dudz"].sel(level=500).data.flatten()
    result_data3 = result_data["brunt_vaisala_frequency"].sel(level=500).data.flatten()

    refer_data1 = np.array(
        [
            4.28113958e-06,
            5.32133720e-06,
            6.06614392e-06,
            5.79035187e-06,
            4.77603961e-06,
            4.66012543e-06,
            4.90554702e-06,
            5.63481462e-06,
            5.95356535e-06,
            5.52850071e-06,
            8.36675438e-06,
            7.70591812e-06,
            7.90070712e-06,
            8.41119332e-06,
            8.44577501e-06,
            1.21285166e-05,
            1.17882178e-05,
            1.22591802e-05,
            1.31184010e-05,
            1.33873526e-05,
            1.26371017e-05,
            1.31981298e-05,
            1.43505631e-05,
            1.55429447e-05,
            1.60624608e-05,
        ]
    )
    refer_data2 = np.array(
        [
            0.00243453,
            0.00320051,
            0.00387504,
            0.00396458,
            0.00349661,
            0.00223617,
            0.0025231,
            0.00307384,
            0.00345408,
            0.00341682,
            0.00410493,
            0.00393663,
            0.004093,
            0.00441058,
            0.00456375,
            0.00639221,
            0.00624336,
            0.00636579,
            0.00669221,
            0.00690693,
            0.00669102,
            0.00683194,
            0.0071831,
            0.00772351,
            0.00815703,
        ]
    )
    refer_data3 = np.array(
        [
            0.00878751,
            0.00929412,
            0.00987128,
            0.0105804,
            0.01131331,
            0.00829668,
            0.00889294,
            0.00943189,
            0.0100312,
            0.01068593,
            0.00936818,
            0.00975453,
            0.00989196,
            0.01001254,
            0.01031783,
            0.01099529,
            0.01104928,
            0.01083314,
            0.01064271,
            0.01076351,
            0.01196115,
            0.0116939,
            0.01130761,
            0.0112256,
            0.01147223,
        ]
    )

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()


def test_calc_apparent_heat_source():
    result_data = ecl.calc_apparent_heat_source(
        u_data=uwnd_daily,
        v_data=vwnd_daily,
        omega_data=omega_daily,
        temper_data=temp_daily,
        vertical_dim="level",
        vertical_dim_units="hPa",
        time_units="day",
    )
    result_data1 = result_data.sel(level=500).isel(time=4).data.flatten()
    refer_data1 = np.array(
        [
            -0.05698965,
            -0.04312353,
            -0.02415958,
            -0.01748332,
            -0.03120244,
            -0.02978676,
            -0.01563223,
            0.00150462,
            0.00404528,
            -0.0112371,
            0.00908426,
            0.03226015,
            0.04679398,
            0.03276806,
            -0.00866827,
            0.04798555,
            0.06669663,
            0.06485339,
            0.04368869,
            0.00017938,
            0.03583697,
            0.04942774,
            0.04090438,
            0.01674637,
            -0.01563623,
        ]
    )
    assert np.isclose(result_data1, refer_data1, atol=0.1, equal_nan=True).all()


def test_calc_total_diabatic_heating():
    result_data = ecl.calc_total_diabatic_heating(
        u_data=uwnd_daily,
        v_data=vwnd_daily,
        omega_data=omega_daily,
        temper_data=temp_daily,
        vertical_dim="level",
        vertical_dim_units="hPa",
        time_units="day",
    )
    result_data1 = result_data.sel(level=500).isel(time=4).data.flatten()
    refer_data1 = np.array(
        [
            -0.05698965,
            -0.04312353,
            -0.02415958,
            -0.01748332,
            -0.03120244,
            -0.02978676,
            -0.01563223,
            0.00150462,
            0.00404528,
            -0.0112371,
            0.00908426,
            0.03226015,
            0.04679398,
            0.03276806,
            -0.00866827,
            0.04798555,
            0.06669663,
            0.06485339,
            0.04368869,
            0.00017938,
            0.03583697,
            0.04942774,
            0.04090438,
            0.01674637,
            -0.01563623,
        ]
    )
    assert np.isclose(result_data1, refer_data1, atol=0.1, equal_nan=True).all()


def test_calc_apparent_moisture_sink():
    result_data = ecl.calc_apparent_moisture_sink(
        u_data=uwnd_daily,
        v_data=vwnd_daily,
        omega_data=omega_daily,
        specific_humidity_data=q_daily,
        vertical_dim="level",
        vertical_dim_units="hPa",
        time_units="day",
        specific_humidity_data_units="kg/kg",
    )
    result_data1 = result_data.sel(level=500).isel(time=3).data.flatten()
    refer_data1 = np.array(
        [
            -0.07268905,
            -0.050381,
            -0.01393382,
            -0.00042548,
            0.00794625,
            0.0075059,
            -0.04543636,
            -0.05116534,
            -0.02060133,
            -0.00412779,
            0.01671084,
            -0.01663374,
            -0.04205875,
            -0.02959669,
            0.00021359,
            -0.08286668,
            -0.02533198,
            -0.02680379,
            -0.03902021,
            -0.05392436,
            -0.20653985,
            -0.05570121,
            -0.02365952,
            -0.04623829,
            -0.09568803,
        ]
    )
    assert np.isclose(result_data1, refer_data1, atol=0.1, equal_nan=True).all()


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_calc_TN_wave_activity_horizontal_flux():
    z500_prime_data = (
        ecl.open_tutorial_dataset("era5_daily_z500_prime_201411_N15.nc").z / 9.8
    )
    u500_climatology_data = ecl.open_tutorial_dataset(
        "era5_ymean_monthly_u500_199101_202012_N15.nc"
    ).u
    v500_climatology_data = ecl.open_tutorial_dataset(
        "era5_ymean_monthly_v500_199101_202012_N15.nc"
    ).v

    tn01_result = calc_TN_wave_activity_horizontal_flux(
        z_prime_data=z500_prime_data,
        u_climatology_data=u500_climatology_data,
        v_climatology_data=v500_climatology_data,
        vertical_dim="level",
        vertical_dim_units="hPa",
    )
    tn01_mean_result = tn01_result.mean(dim="time").sortby("lat")

    draw_shaded = tn01_mean_result["psi_p"].sel(lat=slice(5, None))
    draw_quiver = tn01_mean_result[["fx", "fy"]].sel(lat=slice(5, None))

    fig, ax = plt.subplots()
    draw_shaded.plot.contourf(
        ax=ax, levels=np.linspace(-1e7, 1e7, 21), add_colorbar=False, zorder=1
    )
    draw_quiver.thin(lon=3, lat=1).sel(lat=slice(5, None)).plot.quiver(
        ax=ax,
        x="lon",
        y="lat",
        u="fx",
        v="fy",
        scale=800,
        add_guide=False,
    )
    ax.set_title("")
    return fig


def test_calc_EP_horizontal_flux():
    result_data = ecl.calc_EP_horizontal_flux(
        u_prime_data=u_prime_data.isel(time=slice(None, 3)),
        v_prime_data=v_prime_data.isel(time=slice(None, 3)),
    )
    result_data1 = result_data.sel(level=500).fx.data.flatten()
    result_data2 = result_data.sel(level=500).fy.data.flatten()
    refer_data1 = np.array(
        [
            1.91168248e-01,
            4.28245410e-01,
            6.43983948e-01,
            3.46764909e-02,
            -1.31919014e00,
            5.76269722e-01,
            -8.08320470e-01,
            -1.82824446e00,
            -1.74764765e00,
            -2.15184385e00,
            1.12821221e00,
            -1.86050541e-01,
            -1.18780357e00,
            -9.11238963e-01,
            -1.74250649e-01,
            3.68609400e-04,
            -2.44389936e-01,
            -1.08952057e00,
            -1.91617293e00,
            -1.73710680e00,
            4.28209128e-01,
            1.24807639e-01,
            -1.11164536e00,
            -2.10441846e00,
            -2.12381239e00,
        ]
    )
    refer_data2 = np.array(
        [
            -0.41150835,
            0.23702698,
            0.26110698,
            1.137125,
            2.86586619,
            0.71956657,
            1.68935054,
            2.94322918,
            4.60303644,
            5.285263,
            1.29498139,
            1.71268423,
            2.12802845,
            3.43608455,
            3.62735392,
            1.59939143,
            1.17824546,
            1.08769549,
            2.96109134,
            4.47672289,
            2.41702292,
            1.64206503,
            3.25435964,
            6.7395262,
            8.74013486,
        ]
    )
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_calc_Plumb_wave_activity_horizontal_flux():
    z500_prime_data = (
        ecl.open_tutorial_dataset("era5_daily_z500_prime_201411_N15.nc").z / 9.8
    )

    pwaf_result = calc_Plumb_wave_activity_horizontal_flux(
        z_prime_data=z500_prime_data,
        vertical_dim="level",
        vertical_dim_units="hPa",
    )
    pwaf_mean_result = pwaf_result.mean(dim="time").sortby("lat")

    draw_shaded = pwaf_mean_result["psi_p"].sel(lat=slice(5, None))
    draw_quiver = pwaf_mean_result[["fx", "fy"]].sel(lat=slice(5, None))

    fig, ax = plt.subplots()
    draw_shaded.plot.contourf(
        ax=ax, levels=np.linspace(-1e7, 1e7, 21), add_colorbar=False, zorder=1
    )
    draw_quiver.thin(lon=3, lat=1).sel(lat=slice(5, None)).plot.quiver(
        ax=ax,
        x="lon",
        y="lat",
        u="fx",
        v="fy",
        scale=800,
        add_guide=False,
    )
    ax.set_title("")
    return fig
