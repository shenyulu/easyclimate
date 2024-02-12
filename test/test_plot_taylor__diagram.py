"""
pytest for plot.taylor_diagrams.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

da_a = xr.DataArray(
    np.array([[1, 2, 3], [0.1, 0.2, 0.3], [3.2, 0.6, 1.8]]),
    dims=("lat", "time"),
    coords={
        "lat": np.array([-30, 0, 30]),
        "time": pd.date_range("2000-01-01", freq="D", periods=3),
    },
)
da_b = xr.DataArray(
    np.array([[0.2, 0.4, 0.6], [15, 10, 5], [3.2, 0.6, 1.8]]),
    dims=("lat", "time"),
    coords={
        "lat": np.array([-30, 0, 30]),
        "time": pd.date_range("2000-01-01", freq="D", periods=3),
    },
)
da_obs = (da_a + da_b) / 1.85

taylordiagrams_metadata = ecl.plot.calc_TaylorDiagrams_metadata(
    f=[da_a, da_b],
    r=[da_obs, da_obs],
    models_name=["f1", "f2"],
    weighted=True,
    normalized=True,
)

taylordiagrams_metadata_False = ecl.plot.calc_TaylorDiagrams_metadata(
    f=[da_a, da_b],
    r=[da_obs, da_obs],
    models_name=["f1", "f2"],
    weighted=True,
    normalized=False,
)


def test_calc_TaylorDiagrams_metadata1():
    result_data = ecl.plot.calc_TaylorDiagrams_metadata(
        f=[da_a, da_b],
        r=[da_obs, da_obs],
        models_name=["f1", "f2"],
        weighted=True,
        normalized=True,
    )
    result_data1 = list(result_data["item"].values)
    result_data2 = result_data["std"].values
    result_data3 = result_data["cc"].values
    result_data4 = result_data["centeredRMS"].values
    result_data5 = result_data["TSS"].values

    refer_data1 = ["Obs", "f1", "f2"]
    refer_data2 = np.array([1.0, 0.40462089, 2.05647001])
    refer_data3 = np.array([1.0, -0.42939816, 0.98408606])
    refer_data4 = np.array([0.0, 1.22931078, 1.08700596])
    refer_data5 = np.array([1.0020025, 0.00321028, 0.60040855])

    assert result_data1 == refer_data1
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()


def test_calc_TaylorDiagrams_metadata2():
    result_data = ecl.plot.calc_TaylorDiagrams_metadata(
        f=[da_a, da_b],
        r=[da_obs, da_obs],
        models_name=["f1", "f2"],
        weighted=False,
        normalized=False,
    )
    result_data1 = list(result_data["item"].values)
    result_data2 = result_data["std"].values
    result_data3 = result_data["cc"].values
    result_data4 = result_data["centeredRMS"].values
    result_data5 = result_data["TSS"].values

    refer_data1 = ["Obs", "f1", "f2"]
    refer_data2 = np.array([2.341526, 1.12754513, 4.87543604])
    refer_data3 = np.array([1.0, -0.38222542, 0.97689714])
    refer_data4 = np.array([0.0, 2.96182156, 2.63594059])
    refer_data5 = np.array([1.0020025, 0.00557519, 0.58269386])

    assert result_data1 == refer_data1
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()


@pytest.mark.mpl_image_compare
def test_draw_TaylorDiagrams_base1():
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    ecl.plot.draw_TaylorDiagrams_base(ax=ax, std_max=2.5)
    return fig


# @pytest.mark.mpl_image_compare
# def test_draw_TaylorDiagrams_base2():
#     fig, ax = plt.subplots(subplot_kw = {'projection': 'polar'})
#     ecl.plot.draw_TaylorDiagrams_base(std_max = 2.5, half_circle = False, normalized = True)
#     return fig


@pytest.mark.mpl_image_compare
def test_draw_TaylorDiagrams_base3():
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    ecl.plot.draw_TaylorDiagrams_base(
        ax=ax, std_max=2.5, half_circle=True, normalized=False
    )
    return fig


@pytest.mark.mpl_image_compare
def test_draw_TaylorDiagrams_base4():
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    ecl.plot.draw_TaylorDiagrams_base(
        ax=ax, std_max=2.5, half_circle=False, normalized=False
    )
    return fig


@pytest.mark.mpl_image_compare
def test_draw_TaylorDiagrams_base5():
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    ecl.plot.draw_TaylorDiagrams_base(
        ax=ax, std_max=2.5, half_circle=True, normalized=True
    )
    return fig


def test_draw_TaylorDiagrams_base6():
    fig, ax = plt.subplots()
    with pytest.raises(TypeError):
        ecl.plot.draw_TaylorDiagrams_base(ax=ax, std_max=2.5)
        assert 1 == 1


def test_draw_TaylorDiagrams_base7():
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    with pytest.raises(ValueError):
        ecl.plot.draw_TaylorDiagrams_base(ax=ax, std_min=-1, half_circle=True)
        assert 1 == 1


def test_draw_TaylorDiagrams_base8():
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    with pytest.raises(ValueError):
        ecl.plot.draw_TaylorDiagrams_base(ax=ax, std_min=-1, half_circle=False)
        assert 1 == 1


def test_draw_TaylorDiagrams_base9():
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    with pytest.raises(ValueError):
        ecl.plot.draw_TaylorDiagrams_base(
            ax=ax, std_min=-1, half_circle="This is a sample!"
        )
        assert 1 == 1


def test_draw_TaylorDiagrams_base9():
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    with pytest.raises(ValueError):
        ecl.plot.draw_TaylorDiagrams_base(ax=ax, std_max=-1, half_circle=False)
        assert 1 == 1


@pytest.mark.mpl_image_compare
def test_draw_TaylorDiagrams_metadata1():
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    ecl.plot.draw_TaylorDiagrams_base(ax=ax, std_max=2.5)
    ecl.plot.draw_TaylorDiagrams_metadata(
        taylordiagrams_metadata,
        ax=ax,
        marker_list=["o", "+", "*"],
        color_list=["black", "red", "green"],
        label_list=["1", "", "3"],
        legend_list=taylordiagrams_metadata["item"].to_list(),
        cc="cc",
        std="std",
    )
    return fig


@pytest.mark.mpl_image_compare
def test_draw_TaylorDiagrams_metadata2():
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    ecl.plot.draw_TaylorDiagrams_base(ax=ax, std_max=2.5)
    ecl.plot.draw_TaylorDiagrams_metadata(
        taylordiagrams_metadata_False,
        ax=ax,
        normalized=False,
        marker_list=["o", "+", "*"],
        color_list=["black", "red", "green"],
        label_list=["1", "", "3"],
        legend_list=taylordiagrams_metadata["item"].to_list(),
        cc="cc",
        std="std",
    )
    return fig


def test_draw_TaylorDiagrams_metadata3():
    fig, ax = plt.subplots()
    with pytest.raises(TypeError):
        ecl.plot.draw_TaylorDiagrams_metadata(
            taylordiagrams_metadata,
            ax=ax,
            marker_list=["o", "+", "*"],
            color_list=["black", "red", "green"],
            label_list=["1", "", "3"],
            legend_list=taylordiagrams_metadata["item"].to_list(),
            cc="cc",
            std="std",
            point_label_yoffset=[0.05, 0, 0.05],
            point_label_xoffset=[0.1, 0, 0],
        )
        assert 1 == 1
