"""
pytest for plot/polar.py
"""

import pytest

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import easyclimate as ecl
import cartopy.crs as ccrs


@pytest.mark.mpl_image_compare(tolerance=20)
def test_draw_polar_basemap1():
    fig, ax = plt.subplots(
        figsize=(10 * 0.5, 12 * 0.5), subplot_kw={"projection": ccrs.NorthPolarStereo()}
    )
    ax.coastlines(edgecolor="black", linewidths=0.5)

    gl, meta = ecl.plot.draw_polar_basemap(
        lon_step=30,
        lat_step=10,
        lat_range=[50, 90],
        draw_labels=True,
        lat_label_lon=-30,
    )

    ecl.plot.set_polar_title("Mytest 1", meta, size=15)
    return fig


@pytest.mark.mpl_image_compare(tolerance=20)
def test_draw_polar_basemap2():
    fig, ax = plt.subplots(
        figsize=(10 * 0.5, 12 * 0.5), subplot_kw={"projection": ccrs.NorthPolarStereo()}
    )
    ax.coastlines(edgecolor="black", linewidths=0.5)

    gl, meta = ecl.plot.draw_polar_basemap(
        lon_step=30,
        lat_step=10,
        add_gridlines=False,
        lat_range=[50, 90],
        draw_labels=True,
    )

    ecl.plot.set_polar_title("Mytest 2", meta, ax=ax, size=18, loc="left")
    return fig


@pytest.mark.mpl_image_compare(tolerance=20)
def test_draw_polar_basemap3_southern_hemisphere_reversed_lat_range():
    fig, ax = plt.subplots(
        figsize=(10 * 0.5, 12 * 0.5), subplot_kw={"projection": ccrs.SouthPolarStereo()}
    )
    ax.coastlines(edgecolor="black", linewidths=0.5)

    gl, meta = ecl.plot.draw_polar_basemap(
        lon_step=30,
        lat_step=10,
        lat_range=[-50, -90],
        draw_labels=True,
        lat_label_lon=-30,
    )

    ecl.plot.set_polar_title("Mytest 3", meta, ax=ax, size=15)
    return fig
