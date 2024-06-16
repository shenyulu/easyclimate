# -*- coding: utf-8 -*-
"""
Interpolation and Regriding
===================================

Before proceeding with all the steps, first import some necessary libraries and packages
"""
import easyclimate as ecl
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np

# %%
# Interpolation from points to grid
# ------------------------------------
# Open sample surface pressure data for the European region
data = ecl.open_tutorial_dataset("PressQFF_202007271200_872.csv")
print(data)

# %%
# :py:func:`easyclimate.interp.interp_point2mesh <easyclimate.interp.interp_point2mesh>` enables interpolation from site data to grid point data.
#
# .. seealso::
#
#     - https://github.com/MeteoSwiss/fast-barnes-py
#     - Zürcher, B. K.: Fast approximate Barnes interpolation: illustrated by Python-Numba implementation fast-barnes-py v1.0, Geosci. Model Dev., 16, 1697–1711, https://doi.org/10.5194/gmd-16-1697-2023, 2023.
meshdata = ecl.interp.interp_point2mesh(
    data,
    var_name="qff",
    grid_x=37.5,
    grid_y=75.0,
    point=[-26.0, 34.5],
    resolution=32,
    sigma=1,
)
meshdata

# %%
# Plotting interpolated grid point data and corresponding station locations
fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree(central_longitude=0)})

ax.gridlines(draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--")
ax.coastlines(edgecolor="black", linewidths=0.5)

# Draw interpolation results
meshdata.plot.contourf(
    ax=ax,
    transform=ccrs.PlateCarree(),
    cbar_kwargs={"location": "bottom"},
    cmap="RdBu_r",
    levels=21,
)

# Draw observation stations
ax.scatter(data["lon"], data["lat"], s=1, c="r", transform=ccrs.PlateCarree())

# %%
# Regriding
# ------------------------------------
# Reading example raw grid data
u_data = ecl.tutorial.open_tutorial_dataset("uwnd_202201_mon_mean").sortby("lat").uwnd
u_data

# %%
# Define the target grid (only for **latitude/longitude and regular grids**)
target_grid = xr.DataArray(
    dims=("lat", "lon"),
    coords={
        "lat": np.arange(-89, 89, 6) + 1 / 1.0,
        "lon": np.arange(-180, 180, 6) + 1 / 1.0,
    },
)

# %%
# :py:func:`easyclimate.interp.interp_point2mesh <easyclimate.interp.interp_point2mesh>` performs a regridding operation.
#
# .. seealso::
#
#   https://github.com/EXCITED-CO2/xarray-regrid
regriding_data = ecl.interp.interp_mesh2mesh(u_data, target_grid)
regriding_data

# %%
# Plotting differences before and after interpolation
fig, ax = plt.subplots(1, 2, figsize=(12, 5))

u_data.sel(level=500).isel(time=0).plot(ax=ax[0])
ax[0].set_title("Before", size=20)

regriding_data.sel(level=500).isel(time=0).plot(ax=ax[1])
ax[1].set_title("After", size=20)
