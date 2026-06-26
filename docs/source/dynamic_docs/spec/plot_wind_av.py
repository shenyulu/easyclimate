# -*- coding: utf-8 -*-
"""
Absolute Vorticity (Spherical Harmonics)
==============================================


Absolute vorticity is the sum of relative vorticity and planetary vorticity.
This example calculates absolute vorticity with
:py:func:`easyclimate.spec.calc_absolute_vorticity <easyclimate.core.spec.wind.calc_absolute_vorticity>`
and compares it with
:py:func:`easyclimate.spec.calc_absolute_vorticity_rs <easyclimate.core.spec.wind.calc_absolute_vorticity_rs>`.

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
# The output keeps the same horizontal grid as the input wind and includes both relative and planetary vorticity contributions.
av_fp = ecl.spec.calc_absolute_vorticity(
    u_data=uvdata_500_202201["uwnd"],
    v_data=uvdata_500_202201["vwnd"],
)

av_rs = ecl.spec.calc_absolute_vorticity_rs(
    u_data=uvdata_500_202201["uwnd"],
    v_data=uvdata_500_202201["vwnd"],
)


# %%
# The first plot shows the Fortran-backed absolute vorticity field. A scientific tick formatter is used because vorticity values are small.
fig, ax = plt.subplots(
    figsize = (10, 5),
    subplot_kw={"projection": ccrs.Mercator(central_longitude=180)}
)
ax.coastlines()
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=["bottom", "left"], alpha = 0)

av_fp.sortby("lat").sel(lat=slice(20, 80)).plot.contourf(
    levels=21,
    cmap = "jet",
    cbar_kwargs = {'location': 'bottom', 'format': formatter, 'pad': 0.1},
    transform = ccrs.PlateCarree(),
)


# %%
# The final panel compares Fortran, Rust, and their difference for absolute vorticity.
fig, ax = plt.subplots(1, 3, figsize = (15, 5))

av_fp.sortby("lat").sel(lat=slice(20, 80)).plot.contourf(
    levels=21,
    cmap = "jet",
    ax = ax[0],
    cbar_kwargs = {'location': 'bottom', 'format': formatter},
)
ax[0].set_title("Fortran")

av_rs.sortby("lat").sel(lat=slice(20, 80)).plot.contourf(
    levels=21,
    cmap = "jet",
    ax = ax[1],
    cbar_kwargs = {'location': 'bottom', 'format': formatter},
)
ax[1].set_title("Rust")

(av_fp - av_rs).sortby("lat").sel(lat=slice(20, 80)).plot(
    ax = ax[2],
    cbar_kwargs = {'location': 'bottom'},
)
ax[2].set_title("Diff: Fortran - Rust")
