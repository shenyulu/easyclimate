"""
pytest for wavelet_top.py
"""

import pytest

import easyclimate as ecl
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from .const_define import TEST_DATA_PATH

inputdata_nino3 = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_nino3_wavelet.nc"))
)["nino3"]
result_data_timeseries_wavelet_transform = ecl.filter.calc_timeseries_wavelet_transform(
    inputdata_nino3, dt=0.25
)


def test_timeseries_wavelet_transform1():
    refer_data = xr.open_dataset(
        str(Path(TEST_DATA_PATH, "test_output_waveletpytest1.nc"))
    )
    result_data = ecl.filter.calc_timeseries_wavelet_transform(
        inputdata_nino3,
        dt=0.25,
        mother="morlet",
        sigtest_wavelet="regular chi-square test",
        sigtest_global="time-average test",
    )
    result_data1 = result_data.power.data.flatten()
    result_data2 = result_data.sig.data.flatten()
    result_data3 = result_data.global_ws.data.flatten()
    result_data4 = result_data.global_signif.data.flatten()
    result_data5 = result_data.coi.data.flatten()
    result_data6 = result_data.coi_bottom.data.flatten()

    refer_data1 = refer_data.power.data.flatten()
    refer_data2 = refer_data.sig.data.flatten()
    refer_data3 = refer_data.global_ws.data.flatten()
    refer_data4 = refer_data.global_signif.data.flatten()
    refer_data5 = refer_data.coi.data.flatten()
    refer_data6 = refer_data.coi_bottom.data.flatten()

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()
    assert np.isclose(result_data6, refer_data6).all()


def test_timeseries_wavelet_transform2():
    refer_data = xr.open_dataset(
        str(Path(TEST_DATA_PATH, "test_output_waveletpytest2.nc"))
    )
    result_data = ecl.filter.calc_timeseries_wavelet_transform(
        inputdata_nino3,
        dt=0.25,
        mother="paul",
        sigtest_wavelet="regular chi-square test",
        sigtest_global="regular chi-square test",
    )
    result_data1 = result_data.power.data.flatten()
    result_data2 = result_data.sig.data.flatten()
    result_data3 = result_data.global_ws.data.flatten()
    result_data4 = result_data.global_signif.data.flatten()
    result_data5 = result_data.coi.data.flatten()
    result_data6 = result_data.coi_bottom.data.flatten()

    refer_data1 = refer_data.power.data.flatten()
    refer_data2 = refer_data.sig.data.flatten()
    refer_data3 = refer_data.global_ws.data.flatten()
    refer_data4 = refer_data.global_signif.data.flatten()
    refer_data5 = refer_data.coi.data.flatten()
    refer_data6 = refer_data.coi_bottom.data.flatten()

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()
    assert np.isclose(result_data6, refer_data6).all()


def test_timeseries_wavelet_transform3():
    refer_data = xr.open_dataset(
        str(Path(TEST_DATA_PATH, "test_output_waveletpytest3.nc"))
    )
    result_data = ecl.filter.calc_timeseries_wavelet_transform(
        inputdata_nino3,
        dt=0.25,
        mother="dog",
        sigtest_wavelet="regular chi-square test",
        sigtest_global="regular chi-square test",
    )
    result_data1 = result_data.power.data.flatten()
    result_data2 = result_data.sig.data.flatten()
    result_data3 = result_data.global_ws.data.flatten()
    result_data4 = result_data.global_signif.data.flatten()
    result_data5 = result_data.coi.data.flatten()
    result_data6 = result_data.coi_bottom.data.flatten()

    refer_data1 = refer_data.power.data.flatten()
    refer_data2 = refer_data.sig.data.flatten()
    refer_data3 = refer_data.global_ws.data.flatten()
    refer_data4 = refer_data.global_signif.data.flatten()
    refer_data5 = refer_data.coi.data.flatten()
    refer_data6 = refer_data.coi_bottom.data.flatten()

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()
    assert np.isclose(result_data6, refer_data6).all()


@pytest.mark.mpl_image_compare
def test_draw_global_wavelet_spectrum():
    fig, ax = plt.subplots()
    ecl.filter.draw_global_wavelet_spectrum(
        result_data_timeseries_wavelet_transform, ax=ax
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_draw_wavelet_transform():
    fig, ax = plt.subplots()
    ecl.filter.draw_wavelet_transform(result_data_timeseries_wavelet_transform, ax=ax)
    return fig
