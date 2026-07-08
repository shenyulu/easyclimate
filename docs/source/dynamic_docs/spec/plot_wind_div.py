# -*- coding: utf-8 -*-
"""
Divergence (Spherical Harmonics)
==============================================


Horizontal divergence measures the horizontal outflow or inflow of the wind field.
This example uses
:py:func:`easyclimate.spec.calc_divergence <easyclimate.core.spec.wind.calc_divergence>`
and its Rust-backed counterpart
:py:func:`easyclimate.spec.calc_divergence_rs <easyclimate.core.spec.wind.calc_divergence_rs>`
on a single pressure-level wind field.

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
# The result is an :py:class:`xarray.DataArray <xarray.DataArray>` on the original latitude-longitude grid.
div_fp = ecl.spec.calc_divergence(
    u_data=uvdata_500_202201["uwnd"],
    v_data=uvdata_500_202201["vwnd"],
)

div_rs = ecl.spec.calc_divergence_rs(
    u_data=uvdata_500_202201["uwnd"],
    v_data=uvdata_500_202201["vwnd"],
)


# %%
# The first plot shows the Fortran-backed divergence field on a map projection over the Northern Hemisphere mid-to-high latitudes.
fig, ax = plt.subplots(
    figsize = (10, 5),
    subplot_kw={"projection": ccrs.Mercator(central_longitude=180)}
)
ax.coastlines()
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=["bottom", "left"], alpha = 0)

div_fp.sortby("lat").sel(lat=slice(20, 80)).plot.contourf(
    levels=21,
    cbar_kwargs = {'location': 'bottom', 'format': formatter, 'pad': 0.1},
    transform = ccrs.PlateCarree(),
)


# %%
# The final panel compares Fortran, Rust, and their difference. This is useful when validating backend consistency for the same input field.
fig, ax = plt.subplots(1, 3, figsize = (15, 5))

div_fp.sortby("lat").sel(lat=slice(20, 80)).plot.contourf(
    levels=21,
    ax = ax[0],
    cbar_kwargs = {'location': 'bottom', 'format': formatter},
)
ax[0].set_title("Fortran")

div_rs.sortby("lat").sel(lat=slice(20, 80)).plot.contourf(
    levels=21,
    ax = ax[1],
    cbar_kwargs = {'location': 'bottom', 'format': formatter},
)
ax[1].set_title("Rust")

(div_fp - div_rs).sortby("lat").sel(lat=slice(20, 80)).plot(
    ax = ax[2],
    cbar_kwargs = {'location': 'bottom'},
)
ax[2].set_title("Diff: Fortran - Rust")
