"""
pytest for stat.py
"""

import pytest

import easyclimate as ecl
import numpy as np
from .util import round_sf_np_new
from datatree import DataTree

sst_data = ecl.tutorial.open_tutorial_dataset("mini_HadISST_sst").sst
sic_data_Barents_Sea = ecl.tutorial.open_tutorial_dataset("mini_HadISST_ice").sic
sic_data_Barents_Sea_12 = ecl.get_specific_months_data(sic_data_Barents_Sea, 12)


def test_calc_linregress_spatial():
    result_data = ecl.calc_linregress_spatial(
        sic_data_Barents_Sea_12.sel(lon=slice(34.5, 36.5), lat=slice(78.5, 80.5)),
        dim="time",
    ).compute()
    result_data1 = result_data["slope"].data
    result_data2 = result_data["intercept"].data
    result_data3 = result_data["rvalue"].data
    result_data4 = result_data["pvalue"].data
    result_data5 = result_data["stderr"].data
    result_data6 = result_data["intercept_stderr"].data

    refer_data1 = np.array(
        [
            [-0.01689814, -0.01618345, -0.01640629],
            [-0.01011993, -0.00922373, -0.0091192],
            [-0.00641115, -0.0054169, -0.00600519],
        ]
    )
    refer_data2 = np.array(
        [
            [1.05593576, 1.04461794, 1.04132889],
            [1.01817285, 1.01218153, 1.01741975],
            [0.89857147, 0.91509416, 0.94763013],
        ]
    )
    refer_data3 = np.array(
        [
            [-0.58339207, -0.57217978, -0.57376992],
            [-0.47457306, -0.4485609, -0.45343254],
            [-0.32601495, -0.28470031, -0.33127693],
        ]
    )
    refer_data4 = np.array(
        [
            [5.01586941e-05, 7.52887062e-05, 7.11385575e-05],
            [1.49647207e-03, 2.88846392e-03, 2.56361219e-03],
            [3.51172733e-02, 6.76380713e-02, 3.21088821e-02],
        ]
    )
    refer_data5 = np.array(
        [
            [0.00371969, 0.00366767, 0.00370284],
            [0.00296779, 0.00290584, 0.00283422],
            [0.00293946, 0.00288389, 0.00270435],
        ]
    )
    refer_data6 = np.array(
        [
            [0.08858534, 0.08734656, 0.08818416],
            [0.07067876, 0.06920341, 0.06749764],
            [0.07000405, 0.06868049, 0.06440477],
        ]
    )

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()
    assert np.isclose(result_data6, refer_data6).all()


def test_calc_linregress_spatial_datatree():
    data = sic_data_Barents_Sea_12.sel(
        lon=slice(34.5, 36.5), lat=slice(78.5, 80.5)
    ).to_dataset(name="sic")
    result_data = ecl.calc_linregress_spatial(data, dim="time").compute()
    assert isinstance(result_data, DataTree)


def test_calc_detrend_data():
    result_data = (
        ecl.calc_detrend_data(
            sic_data_Barents_Sea_12.sel(lon=slice(34.5, 36.5), lat=slice(78.5, 80.5)),
            time_dim="time",
        )
        .mean(dim=("lat", "lon"))
        .data
    )
    result_data = round_sf_np_new(result_data)
    refer_data = np.array(
        [
            -5.9e-02,
            -4.9e-02,
            -4.5e-02,
            -4.9e-01,
            -2.2e-02,
            -4.7e-03,
            2.0e-02,
            3.3e-02,
            1.3e-02,
            2.9e-02,
            -3.5e-04,
            2.5e-02,
            5.3e-02,
            3.0e-02,
            5.6e-02,
            1.5e-01,
            5.8e-02,
            1.5e-01,
            1.5e-01,
            1.5e-01,
            2.2e-02,
            1.2e-01,
            1.8e-01,
            1.5e-01,
            1.2e-01,
            -3.3e-01,
            -9.8e-02,
            1.6e-01,
            -1.2e-01,
            2.9e-01,
            1.3e-01,
            -4.2e-01,
            1.8e-01,
            2.7e-01,
            -2.0e-01,
            -6.2e-01,
            -2.3e-01,
            -6.0e-01,
            3.9e-01,
            5.8e-02,
            3.9e-01,
            -7.8e-02,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_ttestSpatialPattern_spatial():
    sic_data_Barents_Sea_12_spatial_mean = sic_data_Barents_Sea_12.sel(
        lon=slice(34.5, 36.5), lat=slice(78.5, 80.5)
    ).mean(dim=("lat", "lon"))
    test1 = sic_data_Barents_Sea_12_spatial_mean.isel(time=slice(0, 20))
    test2 = sic_data_Barents_Sea_12_spatial_mean.isel(time=slice(21, None))
    result_data = ecl.calc_ttestSpatialPattern_spatial(test1, test2)
    result_data1 = result_data["statistic"].data
    result_data2 = result_data["pvalue"].data
    refer_data1 = 3.43382768
    refer_data2 = 0.00142438
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()


def test_calc_levenetestSpatialPattern_spatial():
    sic_data_Barents_Sea_12_spatial_mean = sic_data_Barents_Sea_12.sel(
        lon=slice(34.5, 36.5), lat=slice(78.5, 80.5)
    ).mean(dim=("lat", "lon"))
    test1 = sic_data_Barents_Sea_12_spatial_mean.isel(time=slice(0, 20))
    test2 = sic_data_Barents_Sea_12_spatial_mean.isel(time=slice(21, None))
    result_data = ecl.calc_levenetestSpatialPattern_spatial(test1, test2)
    result_data1 = result_data["statistic"].data
    result_data2 = result_data["pvalue"].data
    refer_data1 = 13.79602081
    refer_data2 = 0.00063671
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()


def test_calc_levenetestSpatialPattern_spatial():
    sic_data_Barents_Sea_12_spatial_mean = sic_data_Barents_Sea_12.sel(
        lon=slice(34.5, 36.5), lat=slice(78.5, 80.5)
    ).mean(dim=("lat", "lon"))
    test1 = sic_data_Barents_Sea_12_spatial_mean.isel(time=slice(0, 20))
    test2 = sic_data_Barents_Sea_12_spatial_mean.isel(time=slice(21, None))
    result_data = ecl.calc_levenetestSpatialPattern_spatial(test1, test2)
    result_data1 = result_data["statistic"].data
    result_data2 = result_data["pvalue"].data
    refer_data1 = 13.79602081
    refer_data2 = 0.00063671
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()


def test_calc_skewness_spatial():
    sic_data_Barents_Sea_12_detrend = ecl.calc_detrend_data(
        sic_data_Barents_Sea_12.sel(lon=slice(34.5, 36.5), lat=slice(78.5, 80.5)),
        time_dim="time",
    )
    result_data = ecl.calc_skewness_spatial(sic_data_Barents_Sea_12_detrend, dim="time")
    result_data1 = result_data["skewness"].data
    result_data2 = result_data["pvalue"].data
    refer_data1 = np.array(
        [
            [-0.24333526, -0.32189173, -0.27430525],
            [-1.3397094, -1.5588326, -1.6165946],
            [-1.8677251, -2.209491, -2.330299],
        ]
    )
    refer_data2 = np.array(
        [
            [7.70400089e-01, 6.26686378e-01, 7.18524740e-01],
            [4.16622922e-04, 3.66448314e-05, 1.56112432e-05],
            [1.88511919e-06, 6.86471937e-08, 1.30767304e-08],
        ]
    )
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()


def test_calc_kurtosis_spatial():
    sic_data_Barents_Sea_12_detrend = ecl.calc_detrend_data(
        sic_data_Barents_Sea_12.sel(lon=slice(34.5, 36.5), lat=slice(78.5, 80.5)),
        time_dim="time",
    )
    result_data = ecl.calc_kurtosis_spatial(sic_data_Barents_Sea_12_detrend, dim="time")
    refer_data = np.array(
        [
            [2.7300231, 2.8368442, 2.7555804],
            [4.667509, 5.464381, 5.85138],
            [6.32586, 7.463421, 8.428185],
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_theilslopes_spatial():
    result_data = ecl.calc_theilslopes_spatial(
        sic_data_Barents_Sea_12.sel(lon=slice(34.5, 36.5), lat=slice(78.5, 80.5)),
        dim="time",
    ).compute()
    result_data1 = result_data["slope"].data.flatten()
    result_data2 = result_data["intercept"].data.flatten()
    result_data3 = result_data["low_slope"].data.flatten()
    result_data4 = result_data["high_slope"].data.flatten()

    refer_data1 = np.array(
        [
            -0.00999999,
            -0.00999999,
            -0.01,
            -0.0035,
            -0.00296296,
            -0.00258064,
            -0.0025,
            -0.00179487,
            -0.001875,
        ]
    )
    refer_data2 = np.array(
        [
            1.10499978,
            1.09499979,
            1.08499995,
            1.00675,
            0.99574073,
            0.99290321,
            0.91125007,
            0.92179486,
            0.93343744,
        ]
    )
    refer_data3 = np.array(
        [
            -0.01833333,
            -0.018,
            -0.019,
            -0.00853659,
            -0.0075,
            -0.00736842,
            -0.005625,
            -0.00416667,
            -0.00333334,
        ]
    )
    refer_data4 = np.array(
        [
            -0.003,
            -0.00285714,
            -0.00333333,
            -0.00090909,
            -0.00058823,
            -0.00058823,
            0.0,
            0.0,
            -0.000625,
        ]
    )

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
