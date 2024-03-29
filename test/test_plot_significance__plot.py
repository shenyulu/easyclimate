"""
pytest for plot.significance_plot.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from pathlib import Path
from .const_define import TEST_DATA_PATH

pvalue_data = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_plot_pvalue.nc"))
).pvalue.sortby("lat")


@pytest.mark.mpl_image_compare
def test_draw_significant_area_contourf1():
    fig, ax = plt.subplots()
    ecl.plot.draw_significant_area_contourf(pvalue_data, ax=ax)
    return fig


# @pytest.mark.mpl_image_compare
# def test_draw_significant_area_contourf2():
#     fig, ax = plt.subplots()
#     ecl.plot.draw_significant_area_contourf(pvalue_data, reverse_level_plot=True)
#     return fig


def test_draw_significant_area_contourf3():
    fig, ax = plt.subplots()

    with pytest.raises(ValueError):
        ecl.plot.draw_significant_area_contourf(
            pvalue_data, reverse_level_plot="This is a sample test!"
        )
        assert 1 == 1


def test_draw_significant_area_contourf4():
    fig, ax = plt.subplots()

    with pytest.raises(ValueError):
        ecl.plot.draw_significant_area_contourf(pvalue_data, thresh=-2)
        assert 1 == 1


def test_get_significance_point():
    sig_points = ecl.plot.get_significance_point(pvalue_data)
    result_data1 = sig_points.lat.values
    result_data2 = sig_points.lon.values
    refer_data1 = np.array(
        [
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -4.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -2.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -1.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            -0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            1.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            3.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
            4.5,
        ]
    )
    refer_data2 = np.array(
        [
            150.5,
            151.5,
            152.5,
            153.5,
            154.5,
            155.5,
            156.5,
            157.5,
            166.5,
            167.5,
            168.5,
            169.5,
            170.5,
            171.5,
            172.5,
            173.5,
            174.5,
            175.5,
            176.5,
            177.5,
            178.5,
            179.5,
            180.5,
            181.5,
            182.5,
            183.5,
            184.5,
            185.5,
            186.5,
            187.5,
            188.5,
            189.5,
            190.5,
            191.5,
            192.5,
            193.5,
            194.5,
            195.5,
            196.5,
            197.5,
            198.5,
            199.5,
            150.5,
            151.5,
            152.5,
            153.5,
            154.5,
            155.5,
            156.5,
            157.5,
            164.5,
            165.5,
            166.5,
            167.5,
            168.5,
            169.5,
            170.5,
            171.5,
            172.5,
            173.5,
            174.5,
            175.5,
            176.5,
            177.5,
            178.5,
            179.5,
            180.5,
            181.5,
            182.5,
            183.5,
            184.5,
            185.5,
            186.5,
            187.5,
            188.5,
            189.5,
            190.5,
            191.5,
            192.5,
            193.5,
            194.5,
            195.5,
            196.5,
            197.5,
            198.5,
            199.5,
            150.5,
            151.5,
            152.5,
            153.5,
            154.5,
            155.5,
            156.5,
            161.5,
            162.5,
            163.5,
            164.5,
            165.5,
            166.5,
            167.5,
            168.5,
            169.5,
            170.5,
            171.5,
            172.5,
            173.5,
            174.5,
            175.5,
            176.5,
            177.5,
            178.5,
            179.5,
            180.5,
            181.5,
            182.5,
            183.5,
            184.5,
            185.5,
            186.5,
            187.5,
            188.5,
            189.5,
            190.5,
            191.5,
            192.5,
            193.5,
            194.5,
            195.5,
            196.5,
            197.5,
            198.5,
            199.5,
            150.5,
            151.5,
            152.5,
            153.5,
            154.5,
            159.5,
            160.5,
            161.5,
            162.5,
            163.5,
            164.5,
            165.5,
            166.5,
            167.5,
            168.5,
            169.5,
            170.5,
            171.5,
            172.5,
            173.5,
            174.5,
            175.5,
            176.5,
            177.5,
            178.5,
            179.5,
            180.5,
            181.5,
            182.5,
            183.5,
            184.5,
            185.5,
            186.5,
            187.5,
            188.5,
            189.5,
            190.5,
            191.5,
            192.5,
            193.5,
            194.5,
            195.5,
            196.5,
            197.5,
            198.5,
            199.5,
            157.5,
            158.5,
            159.5,
            160.5,
            161.5,
            162.5,
            163.5,
            164.5,
            165.5,
            166.5,
            167.5,
            168.5,
            169.5,
            170.5,
            171.5,
            172.5,
            173.5,
            174.5,
            175.5,
            176.5,
            177.5,
            178.5,
            179.5,
            180.5,
            181.5,
            182.5,
            183.5,
            184.5,
            185.5,
            186.5,
            187.5,
            188.5,
            189.5,
            190.5,
            191.5,
            192.5,
            193.5,
            194.5,
            195.5,
            196.5,
            197.5,
            198.5,
            199.5,
            157.5,
            158.5,
            159.5,
            160.5,
            161.5,
            162.5,
            163.5,
            164.5,
            165.5,
            166.5,
            167.5,
            168.5,
            169.5,
            170.5,
            171.5,
            172.5,
            173.5,
            174.5,
            175.5,
            176.5,
            177.5,
            178.5,
            179.5,
            180.5,
            181.5,
            182.5,
            183.5,
            184.5,
            185.5,
            186.5,
            187.5,
            188.5,
            189.5,
            190.5,
            191.5,
            192.5,
            193.5,
            194.5,
            195.5,
            196.5,
            197.5,
            198.5,
            199.5,
            157.5,
            158.5,
            159.5,
            160.5,
            161.5,
            162.5,
            163.5,
            164.5,
            165.5,
            166.5,
            167.5,
            168.5,
            169.5,
            170.5,
            171.5,
            172.5,
            173.5,
            174.5,
            175.5,
            176.5,
            177.5,
            178.5,
            179.5,
            180.5,
            181.5,
            182.5,
            183.5,
            184.5,
            185.5,
            186.5,
            187.5,
            188.5,
            189.5,
            190.5,
            191.5,
            192.5,
            193.5,
            194.5,
            195.5,
            196.5,
            197.5,
            198.5,
            199.5,
            150.5,
            151.5,
            152.5,
            153.5,
            161.5,
            162.5,
            163.5,
            164.5,
            165.5,
            166.5,
            167.5,
            168.5,
            169.5,
            170.5,
            171.5,
            172.5,
            173.5,
            174.5,
            175.5,
            176.5,
            177.5,
            178.5,
            179.5,
            180.5,
            181.5,
            182.5,
            183.5,
            184.5,
            185.5,
            186.5,
            187.5,
            188.5,
            189.5,
            190.5,
            191.5,
            192.5,
            193.5,
            194.5,
            195.5,
            196.5,
            197.5,
            198.5,
            199.5,
            150.5,
            151.5,
            152.5,
            153.5,
            154.5,
            155.5,
            165.5,
            166.5,
            167.5,
            168.5,
            169.5,
            170.5,
            171.5,
            172.5,
            173.5,
            174.5,
            175.5,
            176.5,
            177.5,
            178.5,
            179.5,
            180.5,
            181.5,
            182.5,
            183.5,
            184.5,
            185.5,
            186.5,
            187.5,
            188.5,
            189.5,
            190.5,
            191.5,
            192.5,
            193.5,
            194.5,
            195.5,
            196.5,
            197.5,
            198.5,
            199.5,
            150.5,
            151.5,
            152.5,
            153.5,
            154.5,
            155.5,
            156.5,
            166.5,
            167.5,
            168.5,
            169.5,
            170.5,
            171.5,
            172.5,
            173.5,
            174.5,
            175.5,
            176.5,
            177.5,
            178.5,
            179.5,
            180.5,
            181.5,
            182.5,
            183.5,
            184.5,
            185.5,
            186.5,
            187.5,
            188.5,
            189.5,
            190.5,
            191.5,
            192.5,
            193.5,
            194.5,
            195.5,
            196.5,
            197.5,
            198.5,
            199.5,
        ]
    )
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()


@pytest.mark.mpl_image_compare
def test_draw_significant_area_scatter1():
    sig_points = ecl.plot.get_significance_point(pvalue_data)
    fig, ax = plt.subplots()
    ecl.plot.draw_significant_area_scatter(sig_points, s=8, ax=ax)
    return fig


# @pytest.mark.mpl_image_compare
# def test_draw_significant_area_scatter2():
#     sig_points = ecl.plot.get_significance_point(pvalue_data)
#     fig, ax = plt.subplots()
#     ecl.plot.draw_significant_area_scatter(sig_points, s = 8)
#     return fig
