"""
pytest for filter.spatial_pcf.py
"""

import pytest
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import easyclimate as ecl
from pathlib import Path
from .const_define import DOCS_DATA_PATH


@pytest.fixture(scope="module")
def pcf_result_fast():
    data = ecl.open_tutorial_dataset("uwnd_vwnd_hgt_equtorial_2021_2024").sortby("lat")

    # float64 -> float32
    uwnd = data.uwnd.astype("float32")
    vwnd = data.vwnd.astype("float32")
    hgt = data.hgt.astype("float32")

    result = ecl.filter.filter_2D_spatial_parabolic_cylinder_function(
        uwnd,
        vwnd,
        hgt,
        complex_dtype=np.complex64,
        real_dtype=np.float32,
    )
    return result


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_filter_2D_spatial_parabolic_cylinder_function1(pcf_result_fast):
    fig, ax = plt.subplots()
    pcf_result_fast.isel(time=500).sel(wave_type="kelvin").z.plot.contourf(
        ax=ax, levels=21, add_colorbar=False
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_filter_2D_spatial_parabolic_cylinder_function2(pcf_result_fast):
    fig, ax = plt.subplots()
    pcf_result_fast.isel(time=500).sel(wave_type="wmrg").z.plot.contourf(
        ax=ax, levels=21, add_colorbar=False
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_filter_2D_spatial_parabolic_cylinder_function3(pcf_result_fast):
    fig, ax = plt.subplots()
    pcf_result_fast.isel(time=500).sel(wave_type="r1").z.plot.contourf(
        ax=ax, levels=21, add_colorbar=False
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_filter_2D_spatial_parabolic_cylinder_function4(pcf_result_fast):
    fig, ax = plt.subplots()
    pcf_result_fast.isel(time=500).sel(wave_type="r2").z.plot.contourf(
        ax=ax, levels=21, add_colorbar=False
    )
    return fig
