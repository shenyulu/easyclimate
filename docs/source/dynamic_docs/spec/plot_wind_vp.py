# -*- coding: utf-8 -*-
"""
Velocity Potential (Spherical Harmonics)
==============================================


Velocity potential describes the irrotational, divergent part of the
horizontal wind. This example calculates it with
:py:func:`easyclimate.spec.calc_velocity_potential <easyclimate.core.spec.wind.calc_velocity_potential>`
and compares it with
:py:func:`easyclimate.spec.calc_velocity_potential_rs <easyclimate.core.spec.wind.calc_velocity_potential_rs>`.

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
# The returned velocity potential keeps the input coordinates and can be plotted directly with xarray.
vp_fp = ecl.spec.calc_velocity_potential(
    u_data=uvdata_500_202201["uwnd"],
    v_data=uvdata_500_202201["vwnd"],
)

vp_rs = ecl.spec.calc_velocity_potential_rs(
    u_data=uvdata_500_202201["uwnd"],
    v_data=uvdata_500_202201["vwnd"],
)


# %%
# The first plot shows the velocity potential from the Fortran-backed calculation on a map projection.
fig, ax = plt.subplots(
    figsize = (10, 5),
    subplot_kw={"projection": ccrs.Mercator(central_longitude=180)}
)
ax.coastlines()
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=["bottom", "left"], alpha = 0)

vp_fp.sortby("lat").sel(lat=slice(20, 80)).plot.contourf(
    levels=21,
    cbar_kwargs = {'location': 'bottom', 'format': formatter, 'pad': 0.1},
    transform = ccrs.PlateCarree(),
)


# %%
# The final panel compares Fortran, Rust, and their difference for the velocity potential.
fig, ax = plt.subplots(1, 3, figsize = (15, 5))

vp_fp.sortby("lat").sel(lat=slice(20, 80)).plot.contourf(
    levels=21,
    ax = ax[0],
    cbar_kwargs = {'location': 'bottom'},
)
ax[0].set_title("Fortran")

vp_rs.sortby("lat").sel(lat=slice(20, 80)).plot.contourf(
    levels=21,
    ax = ax[1],
    cbar_kwargs = {'location': 'bottom'},
)
ax[1].set_title("Rust")

(vp_fp - vp_rs).sortby("lat").sel(lat=slice(20, 80)).plot(
    ax = ax[2],
    cbar_kwargs = {'location': 'bottom'},
)
ax[2].set_title("Diff: Fortran - Rust")
