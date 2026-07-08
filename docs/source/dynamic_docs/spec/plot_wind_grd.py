# -*- coding: utf-8 -*-
"""
Geographical Gradient (Spherical Harmonics)
==============================================


The spherical-harmonic gradient returns zonal and meridional derivatives of
a scalar field. This example calculates the gradient of zonal wind with
:py:func:`easyclimate.spec.calc_gradient <easyclimate.core.spec.wind.calc_gradient>`
and compares it with
:py:func:`easyclimate.spec.calc_gradient_rs <easyclimate.core.spec.wind.calc_gradient_rs>`.

"""
import cartopy.crs as ccrs
import xarray as xr
import matplotlib.pyplot as plt
import easyclimate as ecl


# %%
# Open the tutorial zonal and meridional wind components, combine them into one dataset, and select one 500 hPa time slice for the calculation.
u_data = ecl.open_tutorial_dataset("uwnd_2022_day5").uwnd
v_data = ecl.open_tutorial_dataset("vwnd_2022_day5").vwnd

uvdata = xr.Dataset()
uvdata["uwnd"] = u_data
uvdata["vwnd"] = v_data

uvdata_500_202201 = uvdata.sel(level=500).isel(time = 3)
uvdata_500_202201


# %%
# Prepare a scientific-notation formatter for the colorbars. Many spectral wind diagnostics have small physical units, so this keeps the labels readable.
import matplotlib.ticker as ticker
formatter = ticker.ScalarFormatter(useMathText=True, useOffset=True)
formatter.set_scientific(True)
formatter.set_powerlimits((0, 0))

# %%
# The returned dataset contains ``zonal_gradient`` and ``meridional_gradient`` on the input grid.
uvgrd_fp = ecl.spec.calc_gradient(
    data_input=uvdata_500_202201["uwnd"],
)

uvgrd_rs = ecl.spec.calc_gradient_rs(
    data_input=uvdata_500_202201["uwnd"],
)


# %%
# The first figure maps the zonal and meridional gradients from the Fortran-backed calculation.
fig, ax = plt.subplots(
    2, 1,
    figsize = (10, 10),
    subplot_kw={"projection": ccrs.Mercator(central_longitude=180)}
)

for axi in ax.flat:
    axi.coastlines()
    axi.gridlines(crs=ccrs.PlateCarree(), draw_labels=["bottom", "left"], alpha = 0)

axi = ax[0]
uvgrd_fp["zonal_gradient"].sortby("lat").sel(lat=slice(20, 80)).plot.contourf(
    levels=21,
    ax = axi,
    cbar_kwargs = {'location': 'bottom', 'format': formatter, 'pad': 0.1},
    transform = ccrs.PlateCarree(),
)

axi = ax[1]
uvgrd_fp["meridional_gradient"].sortby("lat").sel(lat=slice(20, 80)).plot.contourf(
    levels=21,
    ax = axi,
    cbar_kwargs = {'location': 'bottom', 'format': formatter, 'pad': 0.1},
    transform = ccrs.PlateCarree(),
)


# %%
# The final figure compares each gradient component from the two backends and plots the component-wise differences.
fig, ax = plt.subplots(2, 3, figsize = (15, 10))

uvgrd_fp["zonal_gradient"].sortby("lat").sel(lat=slice(20, 80)).plot.contourf(
    levels=21,
    ax = ax[0, 0],
    cbar_kwargs = {'location': 'bottom'},
)
ax[0, 0].set_title("Fortran")

uvgrd_rs["zonal_gradient"].sortby("lat").sel(lat=slice(20, 80)).plot.contourf(
    levels=21,
    ax = ax[0, 1],
    cbar_kwargs = {'location': 'bottom'},
)
ax[0, 1].set_title("Rust")

(uvgrd_fp["zonal_gradient"] - uvgrd_rs["zonal_gradient"]).sortby("lat").sel(lat=slice(20, 80)).plot(
    ax = ax[0, 2],
    cbar_kwargs = {'location': 'bottom'},
)
ax[0, 2].set_title("Diff: Fortran - Rust")

uvgrd_fp["meridional_gradient"].sortby("lat").sel(lat=slice(20, 80)).plot.contourf(
    levels=21,
    ax = ax[1, 0],
    cbar_kwargs = {'location': 'bottom'},
)
ax[1, 0].set_title("Fortran")

uvgrd_rs["meridional_gradient"].sortby("lat").sel(lat=slice(20, 80)).plot.contourf(
    levels=21,
    ax = ax[1, 1],
    cbar_kwargs = {'location': 'bottom'},
)
ax[1, 1].set_title("Rust")

(uvgrd_fp["meridional_gradient"] - uvgrd_rs["meridional_gradient"]).sortby("lat").sel(lat=slice(20, 80)).plot(
    ax = ax[1, 2],
    cbar_kwargs = {'location': 'bottom'},
)
ax[1, 2].set_title("Diff: Fortran - Rust")
