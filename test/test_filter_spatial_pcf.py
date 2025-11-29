"""
pytest for filter.spatial_pcf.py
"""

import pytest
import xarray as xr
import matplotlib.pyplot as plt
import easyclimate as ecl
from pathlib import Path
from .const_define import DOCS_DATA_PATH

data = ecl.open_tutorial_dataset("uwnd_vwnd_hgt_equtorial_2021_2024").sortby("lat")

result = ecl.filter.filter_2D_spatial_parabolic_cylinder_function(
    data.uwnd, data.vwnd, data.hgt
)


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_filter_2D_spatial_parabolic_cylinder_function1():
    fig, ax = plt.subplots()
    result.isel(time=500).sel(wave_type="kelvin").z.plot.contourf(
        ax=ax, levels=21, add_colorbar=False
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_filter_2D_spatial_parabolic_cylinder_function2():
    fig, ax = plt.subplots()
    result.isel(time=500).sel(wave_type="wmrg").z.plot.contourf(
        ax=ax, levels=21, add_colorbar=False
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_filter_2D_spatial_parabolic_cylinder_function3():
    fig, ax = plt.subplots()
    result.isel(time=500).sel(wave_type="r1").z.plot.contourf(
        ax=ax, levels=21, add_colorbar=False
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_filter_2D_spatial_parabolic_cylinder_function4():
    fig, ax = plt.subplots()
    result.isel(time=500).sel(wave_type="r2").z.plot.contourf(
        ax=ax, levels=21, add_colorbar=False
    )
    return fig
