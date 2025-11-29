"""
pytest for filter/kf_filter.py
"""

import pytest
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import easyclimate as ecl
from pathlib import Path
from .const_define import DOCS_DATA_PATH

olr_data_interpolated = ecl.open_tutorial_dataset("olr_smooth_data").sortby("lat").olr

#
lf_result = ecl.filter.kf_filter_lf_wave(olr_data_interpolated, steps_per_day=1)
mjo_result = ecl.filter.kf_filter_mjo_wave(olr_data_interpolated, steps_per_day=1)
er_result = ecl.filter.kf_filter_er_wave(olr_data_interpolated, steps_per_day=1)
kelvin_result = ecl.filter.kf_filter_kelvin_wave(olr_data_interpolated, steps_per_day=1)
mt_result = ecl.filter.kf_filter_mt_wave(olr_data_interpolated, steps_per_day=1)
mrg_result = ecl.filter.kf_filter_mrg_wave(olr_data_interpolated, steps_per_day=1)
td_result = ecl.filter.kf_filter_td_wave(olr_data_interpolated, steps_per_day=1)

#
time1, time2 = "2017-12-01", "2018-02-28"
lon1, lon2, lats, latn = 39, 181, 5, 15

lf_result_ave = lf_result.sel(
    time=slice(time1, time2), lat=slice(lats, latn), lon=slice(lon1, lon2)
).mean(dim="lat")
mjo_result_ave = mjo_result.sel(
    time=slice(time1, time2), lat=slice(lats, latn), lon=slice(lon1, lon2)
).mean(dim="lat")
er_result_ave = er_result.sel(
    time=slice(time1, time2), lat=slice(lats, latn), lon=slice(lon1, lon2)
).mean(dim="lat")
kelvin_result_ave = kelvin_result.sel(
    time=slice(time1, time2), lat=slice(lats, latn), lon=slice(lon1, lon2)
).mean(dim="lat")
mt_result_ave = mt_result.sel(
    time=slice(time1, time2), lat=slice(lats, latn), lon=slice(lon1, lon2)
).mean(dim="lat")
mrg_result_ave = mrg_result.sel(
    time=slice(time1, time2), lat=slice(lats, latn), lon=slice(lon1, lon2)
).mean(dim="lat")
td_result_ave = td_result.sel(
    time=slice(time1, time2), lat=slice(lats, latn), lon=slice(lon1, lon2)
).mean(dim="lat")


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_kf_filter_lf_wave():
    fig, ax = plt.subplots(figsize=(8, 5))

    lf_result_ave.plot.contourf(
        ax=ax,
        yincrease=False,
        levels=np.linspace(-15, 15, 21),
        cbar_kwargs={"location": "right", "aspect": 30},
    )
    ecl.plot.set_lon_format_axis()
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_kf_filter_mjo_wave():
    fig, ax = plt.subplots(figsize=(8, 5))

    mjo_result_ave.plot.contourf(
        ax=ax,
        yincrease=False,
        levels=np.linspace(-15, 15, 21),
        cbar_kwargs={"location": "right", "aspect": 30},
    )
    ecl.plot.set_lon_format_axis()
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_kf_filter_er_wave():
    fig, ax = plt.subplots(figsize=(8, 5))

    er_result_ave.plot.contourf(
        ax=ax,
        yincrease=False,
        levels=np.linspace(-15, 15, 21),
        cbar_kwargs={"location": "right", "aspect": 30},
    )
    ecl.plot.set_lon_format_axis()
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_kf_filter_kelvin_wave():
    fig, ax = plt.subplots(figsize=(8, 5))

    kelvin_result_ave.plot.contourf(
        ax=ax,
        yincrease=False,
        levels=np.linspace(-15, 15, 21),
        cbar_kwargs={"location": "right", "aspect": 30},
    )
    ecl.plot.set_lon_format_axis()
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_kf_filter_mt_wave():
    fig, ax = plt.subplots(figsize=(8, 5))

    mt_result_ave.plot.contourf(
        ax=ax,
        yincrease=False,
        levels=np.linspace(-15, 15, 21),
        cbar_kwargs={"location": "right", "aspect": 30},
    )
    ecl.plot.set_lon_format_axis()
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_kf_filter_mrg_wave():
    fig, ax = plt.subplots(figsize=(8, 5))

    mrg_result_ave.plot.contourf(
        ax=ax,
        yincrease=False,
        levels=np.linspace(-15, 15, 21),
        cbar_kwargs={"location": "right", "aspect": 30},
    )
    ecl.plot.set_lon_format_axis()
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_kf_filter_td_wave():
    fig, ax = plt.subplots(figsize=(8, 5))

    td_result_ave.plot.contourf(
        ax=ax,
        yincrease=False,
        levels=np.linspace(-15, 15, 21),
        cbar_kwargs={"location": "right", "aspect": 30},
    )
    ecl.plot.set_lon_format_axis()
    return fig
