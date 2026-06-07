# -*- coding: utf-8 -*-
"""
ICON Cell-centered triangular plots
==============================================

This example draws ICON cell-centered scalar fields on the native triangular
mesh with
:py:func:`easyclimate.plot.icon.plot_cell_triangular <easyclimate.plot.icon.plot_cell_triangular>`.

Unlike contour plots, the values are rendered cell by cell, preserving the
native ICON triangular footprint in the selected domain.
"""
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import easyclimate as ecl

from easyclimate.plot.icon import plot_cell_triangular

# %%
# Open a sample ICON pressure-level file and inspect the dataset. The
# temperature field is stored on ICON cell centers.
data = xr.open_dataset("icon_native_pl_20080909T000000Z.nc")
data

# %%
# Draw cell-centered temperature directly on plain Matplotlib axes.
fig, ax = plt.subplots()

pc = plot_cell_triangular(
    data,
    data.temp.isel(time=0, plev=0),

    ax=ax,

    lon_min=110,
    lon_max=220,
    lat_min=25,
    lat_max=70,

    cmap="RdBu_r",
    cbar_kwargs={"location": "bottom", "aspect": 60},
)

ax.set_title("JW Wave 850hPa Temperature\n(plot_cell_triangular)")

# %%
# Pass a Cartopy GeoAxes and a geographic transform when the native ICON cell
# fills should be combined with coastlines and gridlines.
fig, ax = plt.subplots(
    subplot_kw={"projection": ccrs.PlateCarree(-120)}
)

pc = plot_cell_triangular(
    data,
    data.temp.isel(time=0, plev=0),

    ax=ax,
    transform=ccrs.PlateCarree(),

    lon_min=110,
    lon_max=220,
    lat_min=25,
    lat_max=70,

    cmap="RdBu_r",
    cbar_kwargs={"location": "bottom", "aspect": 60},
)

ax.set_title("JW Wave 850hPa Temperature\n(plot_cell_triangular)")
ax.coastlines(resolution="110m", linewidth=0.6, color="b")
ax.gridlines(draw_labels=True, alpha=0)

# %%
# Zoom into a smaller regional window without slicing the ICON dataset first.
fig, ax = plt.subplots(
    subplot_kw={"projection": ccrs.PlateCarree(-120)}
)

pc = plot_cell_triangular(
    data,
    data.temp.isel(time=0, plev=0),

    ax=ax,
    transform=ccrs.PlateCarree(),

    lon_min=110,
    lon_max=150,
    lat_min=45,
    lat_max=65,

    cmap="RdBu_r",
    cbar_kwargs={"location": "bottom", "aspect": 60},
)

ax.set_title("JW Wave 850hPa Temperature\n(plot_cell_triangular)")
ax.coastlines(resolution="50m", linewidth=0.6, color="b")
ax.gridlines(draw_labels=True, alpha=0)
