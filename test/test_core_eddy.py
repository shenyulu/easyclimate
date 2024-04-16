"""
pytest for eddy.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
import pandas as pd

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
    ecl.open_tutorial_dataset("hgt_day_ltm_1991_2020_0to6day.nc")
    .hgt.sortby("lat")
    .sel(lon=slice(100, 110), lat=slice(20, 30))
)
u_climate_data = (
    ecl.open_tutorial_dataset("uwnd_day_ltm_1991_2020_0to6day.nc")
    .uwnd.sortby("lat")
    .sel(lon=slice(100, 110), lat=slice(20, 30))
)
v_climate_data = (
    ecl.open_tutorial_dataset("vwnd_day_ltm_1991_2020_0to6day.nc")
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
            4.19551679e-05,
            5.21491046e-05,
            5.94482104e-05,
            5.67454483e-05,
            4.68051882e-05,
            4.56692292e-05,
            4.80743608e-05,
            5.52211832e-05,
            5.83449405e-05,
            5.41793069e-05,
            8.19941930e-05,
            7.55179976e-05,
            7.74269298e-05,
            8.24296946e-05,
            8.27685951e-05,
            1.18859462e-04,
            1.15524534e-04,
            1.20139966e-04,
            1.28560330e-04,
            1.31196055e-04,
            1.23843597e-04,
            1.29341672e-04,
            1.40635518e-04,
            1.52320858e-04,
            1.57412116e-04,
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
            -0.05699683,
            -0.04312915,
            -0.02416359,
            -0.01748616,
            -0.03120409,
            -0.02979126,
            -0.01563774,
            0.0014992,
            0.00404129,
            -0.01123914,
            0.00908052,
            0.03225594,
            0.04678952,
            0.03276326,
            -0.0086737,
            0.04798427,
            0.06669269,
            0.06484717,
            0.04368201,
            0.00017335,
            0.03583446,
            0.04942282,
            0.04089638,
            0.01673753,
            -0.01564289,
        ]
    )
    assert np.isclose(result_data1, refer_data1).all()


def test_calc_apparent_moisture_sink():
    result_data = ecl.calc_apparent_moisture_sink(
        u_data=uwnd_daily,
        v_data=vwnd_daily,
        omega_data=omega_daily,
        specific_humidity_data=q_daily,
        vertical_dim="level",
        vertical_dim_units="hPa",
        time_units="day",
        specific_humidity_units="kg/kg",
    )
    result_data1 = result_data.sel(level=500).isel(time=3).data.flatten()
    refer_data1 = np.array(
        [
            -0.07269812,
            -0.05038757,
            -0.01393634,
            -0.00042704,
            0.00794611,
            0.00750676,
            -0.04544461,
            -0.05117434,
            -0.02060636,
            -0.00413119,
            0.01671645,
            -0.01663402,
            -0.0420649,
            -0.02960256,
            0.00021111,
            -0.08287747,
            -0.02533122,
            -0.02680349,
            -0.03902617,
            -0.05393757,
            -0.20656966,
            -0.05570683,
            -0.02365906,
            -0.04624587,
            -0.09570991,
        ]
    )
    assert np.isclose(result_data1, refer_data1).all()


def test_calc_TN_wave_activity_horizontal_flux():
    result_data = ecl.calc_TN_wave_activity_horizontal_flux(
        z_prime_data=z_prime_data,
        u_climatology_data=u_climate_data_need,
        v_climatology_data=v_climate_data_need,
        vertical_dim="level",
        vertical_dim_units="hPa",
    )
    result_data1 = result_data.sel(level=500).isel(time=2).psi_p.data.flatten()
    result_data2 = result_data.sel(level=500).isel(time=2).fx.data.flatten()
    result_data3 = result_data.sel(level=500).isel(time=2).fy.data.flatten()

    refer_data1 = np.array(
        [
            -3634707.5334086,
            -3094413.24061066,
            -2652354.16253647,
            -1817353.7667043,
            -1080588.7385502,
            -2809504.61985641,
            -1887635.85238713,
            -1097462.72077178,
            -526782.13331078,
            43898.50968526,
            -2027268.91556854,
            -1192511.08311649,
            -318002.96168694,
            477004.45800343,
            993759.23593041,
            -1673561.50696147,
            -873162.53522288,
            -72763.54224161,
            691253.62651218,
            1491652.6548979,
            -1175946.24245749,
            -335984.64070214,
            436780.03814413,
            1377537.08965492,
            2284695.64047608,
        ]
    )
    refer_data2 = np.array(
        [
            0.00124622,
            0.00552013,
            0.00694739,
            0.00494649,
            0.00269781,
            0.00555633,
            0.00411022,
            0.00295194,
            0.002603,
            0.00275766,
            0.00563715,
            0.00532329,
            0.00485144,
            0.00392754,
            0.00328304,
            0.00468875,
            0.0045456,
            0.00441541,
            0.00432411,
            0.00416847,
            0.00423191,
            0.00464009,
            0.00520513,
            0.00620886,
            0.00751868,
        ]
    )
    refer_data3 = np.array(
        [
            0.01910575,
            0.01506972,
            0.00945703,
            0.00547807,
            0.00651758,
            0.00698939,
            0.0070297,
            0.00603156,
            0.00483476,
            0.0046966,
            0.00162786,
            0.00240627,
            0.0028132,
            0.00254405,
            0.00147791,
            0.00245604,
            0.00208707,
            0.00189897,
            0.00164943,
            0.00106359,
            0.00378404,
            0.00331469,
            0.00307501,
            0.00386894,
            0.00663514,
        ]
    )
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()


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
