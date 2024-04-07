"""
pytest for eof.py
"""

import pytest

import easyclimate as ecl
import xarray as xr
import numpy as np
import pandas as pd

data_time_series = xr.DataArray(
    np.array(
        [
            [57413.867, 57411.96],
            [57404.35, 57404.35],
            [57478.6, 57482.406],
            [57531.906, 57531.906],
            [57564.273, 57569.984],
            [57474.79, 57472.887],
            [57309.15, 57307.246],
            [57343.42, 57343.42],
            [57408.152, 57402.44],
            [57442.426, 57442.426],
            [57478.6, 57480.504],
            [57585.215, 57583.312],
            [57568.082, 57569.984],
            [57625.195, 57627.1],
            [57653.758, 57653.758],
            [57697.547, 57697.547],
            [57613.773, 57619.49],
            [57501.445, 57501.445],
            [57446.23, 57446.23],
            [57440.523, 57438.617],
            [57438.617, 57440.523],
            [57387.21, 57381.5],
            [57507.156, 57509.062],
            [57469.08, 57469.08],
            [57349.133, 57349.133],
            [57371.984, 57373.883],
            [57499.54, 57503.35],
            [57512.867, 57516.676],
            [57421.48, 57427.195],
            [57259.65, 57261.555],
        ]
    ),
    dims=("time", "lon"),
    coords={
        "time": pd.date_range("1982-01-01", periods=30, freq="ME"),
        "lon": np.array([100.125, 101.25]),
    },
)


def test_original_test():
    result_data = ecl.mk_test.original_test(data_time_series, dim="time")
    result_data1 = result_data["trend"].data
    result_data2 = result_data["h"].data
    result_data3 = result_data["p"].data
    result_data4 = result_data["z"].data
    result_data5 = result_data["Tau"].data
    result_data6 = result_data["s"].data
    result_data7 = result_data["var_s"].data
    result_data8 = result_data["slope"].data
    result_data9 = result_data["intercept"].data

    refer_data1 = np.array([0.0, 0.0])
    refer_data2 = np.array([0.0, 0.0])
    refer_data3 = np.array([0.62995903, 0.62995903])
    refer_data4 = np.array([-0.48178452, -0.48178452])
    refer_data5 = np.array([-0.06436782, -0.06436782])
    refer_data6 = np.array([-28.0, -28.0])
    refer_data7 = np.array([3140.66666667, 3140.66666667])
    refer_data8 = np.array([-1.34858333, -1.26911111])
    refer_data9 = np.array([57491.48945833, 57489.38561111])

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()
    assert np.isclose(result_data6, refer_data6).all()
    assert np.isclose(result_data7, refer_data7).all()
    assert np.isclose(result_data8, refer_data8).all()
    assert np.isclose(result_data9, refer_data9).all()


def test_hamed_rao_modification_test():
    result_data = ecl.mk_test.hamed_rao_modification_test(data_time_series, dim="time")
    result_data1 = result_data["trend"].data
    result_data2 = result_data["h"].data
    result_data3 = result_data["p"].data
    result_data4 = result_data["z"].data
    result_data5 = result_data["Tau"].data
    result_data6 = result_data["s"].data
    result_data7 = result_data["var_s"].data
    result_data8 = result_data["slope"].data
    result_data9 = result_data["intercept"].data

    refer_data1 = np.array([0.0, 0.0])
    refer_data2 = np.array([0.0, 0.0])
    refer_data3 = np.array([0.73706741, 0.73542695])
    refer_data4 = np.array([-0.33573938, -0.33791539])
    refer_data5 = np.array([-0.06436782, -0.06436782])
    refer_data6 = np.array([-28.0, -28.0])
    refer_data7 = np.array([6467.29937695, 6384.27507788])
    refer_data8 = np.array([-1.34858333, -1.26911111])
    refer_data9 = np.array([57491.48945833, 57489.38561111])

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()
    assert np.isclose(result_data6, refer_data6).all()
    assert np.isclose(result_data7, refer_data7).all()
    assert np.isclose(result_data8, refer_data8).all()
    assert np.isclose(result_data9, refer_data9).all()


def test_yue_wang_modification_test():
    result_data = ecl.mk_test.yue_wang_modification_test(data_time_series, dim="time")
    result_data1 = result_data["trend"].data
    result_data2 = result_data["h"].data
    result_data3 = result_data["p"].data
    result_data4 = result_data["z"].data
    result_data5 = result_data["Tau"].data
    result_data6 = result_data["s"].data
    result_data7 = result_data["var_s"].data
    result_data8 = result_data["slope"].data
    result_data9 = result_data["intercept"].data

    refer_data1 = np.array([0.0, 0.0])
    refer_data2 = np.array([0.0, 0.0])
    refer_data3 = np.array([0.52524504, 0.51986232])
    refer_data4 = np.array([-0.63528119, -0.64355765])
    refer_data5 = np.array([-0.06436782, -0.06436782])
    refer_data6 = np.array([-28.0, -28.0])
    refer_data7 = np.array([1806.32350299, 1760.16188858])
    refer_data8 = np.array([-1.34858333, -1.26911111])
    refer_data9 = np.array([57491.48945833, 57489.38561111])

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()
    assert np.isclose(result_data6, refer_data6).all()
    assert np.isclose(result_data7, refer_data7).all()
    assert np.isclose(result_data8, refer_data8).all()
    assert np.isclose(result_data9, refer_data9).all()


def test_pre_whitening_modification_test():
    result_data = ecl.mk_test.pre_whitening_modification_test(
        data_time_series, dim="time"
    )
    result_data1 = result_data["trend"].data
    result_data2 = result_data["h"].data
    result_data3 = result_data["p"].data
    result_data4 = result_data["z"].data
    result_data5 = result_data["Tau"].data
    result_data6 = result_data["s"].data
    result_data7 = result_data["var_s"].data
    result_data8 = result_data["slope"].data
    result_data9 = result_data["intercept"].data

    refer_data1 = np.array([0.0, 0.0])
    refer_data2 = np.array([0.0, 0.0])
    refer_data3 = np.array([0.22274021, 0.23730178])
    refer_data4 = np.array([-1.21927402, -1.1817579])
    refer_data5 = np.array([-0.16256158, -0.15763547])
    refer_data6 = np.array([-66.0, -64.0])
    refer_data7 = np.array([2842.0, 2842.0])
    refer_data8 = np.array([-1.34858333, -1.26911111])
    refer_data9 = np.array([57491.48945833, 57489.38561111])

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()
    assert np.isclose(result_data6, refer_data6).all()
    assert np.isclose(result_data7, refer_data7).all()
    assert np.isclose(result_data8, refer_data8).all()
    assert np.isclose(result_data9, refer_data9).all()


def test_trend_free_pre_whitening_modification_test():
    result_data = ecl.mk_test.trend_free_pre_whitening_modification_test(
        data_time_series, dim="time"
    )
    result_data1 = result_data["trend"].data
    result_data2 = result_data["h"].data
    result_data3 = result_data["p"].data
    result_data4 = result_data["z"].data
    result_data5 = result_data["Tau"].data
    result_data6 = result_data["s"].data
    result_data7 = result_data["var_s"].data
    result_data8 = result_data["slope"].data
    result_data9 = result_data["intercept"].data

    refer_data1 = np.array([0.0, 0.0])
    refer_data2 = np.array([0.0, 0.0])
    refer_data3 = np.array([0.11083858, 0.11948969])
    refer_data4 = np.array([-1.59443526, -1.55691913])
    refer_data5 = np.array([-0.21182266, -0.20689655])
    refer_data6 = np.array([-86.0, -84.0])
    refer_data7 = np.array([2842.0, 2842.0])
    refer_data8 = np.array([-1.34858333, -1.26911111])
    refer_data9 = np.array([57491.48945833, 57489.38561111])

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()
    assert np.isclose(result_data6, refer_data6).all()
    assert np.isclose(result_data7, refer_data7).all()
    assert np.isclose(result_data8, refer_data8).all()
    assert np.isclose(result_data9, refer_data9).all()


def test_seasonal_test():
    result_data = ecl.mk_test.seasonal_test(data_time_series, dim="time")
    result_data1 = result_data["trend"].data
    result_data2 = result_data["h"].data
    result_data3 = result_data["p"].data
    result_data4 = result_data["z"].data
    result_data5 = result_data["Tau"].data
    result_data6 = result_data["s"].data
    result_data7 = result_data["var_s"].data
    result_data8 = result_data["slope"].data
    result_data9 = result_data["intercept"].data

    refer_data1 = np.array([0.0, 0.0])
    refer_data2 = np.array([0.0, 0.0])
    refer_data3 = np.array([0.85010674, 0.85010674])
    refer_data4 = np.array([-0.18898224, -0.18898224])
    refer_data5 = np.array([-0.08333333, -0.08333333])
    refer_data6 = np.array([-2.0, -2.0])
    refer_data7 = np.array([28.0, 28.0])
    refer_data8 = np.array([-12.85125, -11.42425])
    refer_data9 = np.array([57487.46359375, 57484.78780208])

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()
    assert np.isclose(result_data6, refer_data6).all()
    assert np.isclose(result_data7, refer_data7).all()
    assert np.isclose(result_data8, refer_data8).all()
    assert np.isclose(result_data9, refer_data9).all()


def test_regional_test():
    result_data = ecl.mk_test.regional_test(data_time_series, dim="time")
    result_data1 = result_data["trend"].data
    result_data2 = result_data["h"].data
    result_data3 = result_data["p"].data
    result_data4 = result_data["z"].data
    result_data5 = result_data["Tau"].data
    result_data6 = result_data["s"].data
    result_data7 = result_data["var_s"].data
    result_data8 = result_data["slope"].data
    result_data9 = result_data["intercept"].data

    refer_data1 = np.array([0.0, 0.0])
    refer_data2 = np.array([0.0, 0.0])
    refer_data3 = np.array([0.62995903, 0.62995903])
    refer_data4 = np.array([-0.48178452, -0.48178452])
    refer_data5 = np.array([-0.06436782, -0.06436782])
    refer_data6 = np.array([-28.0, -28.0])
    refer_data7 = np.array([3140.66666667, 3140.66666667])
    refer_data8 = np.array([-1.34858333, -1.26911111])
    refer_data9 = np.array([57491.48945833, 57489.38561111])

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()
    assert np.isclose(result_data6, refer_data6).all()
    assert np.isclose(result_data7, refer_data7).all()
    assert np.isclose(result_data8, refer_data8).all()
    assert np.isclose(result_data9, refer_data9).all()


def test_correlated_seasonal_test():
    result_data = ecl.mk_test.correlated_seasonal_test(data_time_series, dim="time")
    result_data1 = result_data["trend"].data
    result_data2 = result_data["h"].data
    result_data3 = result_data["p"].data
    result_data4 = result_data["z"].data
    result_data5 = result_data["Tau"].data
    result_data6 = result_data["s"].data
    result_data7 = result_data["var_s"].data
    result_data8 = result_data["slope"].data
    result_data9 = result_data["intercept"].data

    refer_data1 = np.array([0.0, 0.0])
    refer_data2 = np.array([0.0, 0.0])
    refer_data3 = np.array([0.31731051, 0.31731051])
    refer_data4 = np.array([1.0, 1.0])
    refer_data5 = np.array([0.66666667, 0.66666667])
    refer_data6 = np.array([8.0, 8.0])
    refer_data7 = np.array([64.0, 64.0])
    refer_data8 = np.array([-12.85125, -11.42425])
    refer_data9 = np.array([57487.46359375, 57484.78780208])

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()
    assert np.isclose(result_data6, refer_data6).all()
    assert np.isclose(result_data7, refer_data7).all()
    assert np.isclose(result_data8, refer_data8).all()
    assert np.isclose(result_data9, refer_data9).all()


def test_sens_slope():
    result_data = ecl.mk_test.sens_slope(data_time_series, dim="time")
    result_data1 = result_data["slope"].data
    result_data2 = result_data["intercept"].data

    refer_data1 = np.array([-1.34858333, -1.26911111])
    refer_data2 = np.array([57491.48945833, 57489.38561111])

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()


def test_seasonal_sens_slope():
    result_data = ecl.mk_test.seasonal_sens_slope(data_time_series, dim="time")
    result_data1 = result_data["slope"].data
    result_data2 = result_data["intercept"].data

    refer_data1 = np.array([-12.85125, -11.42425])
    refer_data2 = np.array([57487.46359375, 57484.78780208])

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
