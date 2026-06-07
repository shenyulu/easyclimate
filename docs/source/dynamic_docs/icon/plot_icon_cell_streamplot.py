# -*- coding: utf-8 -*-
"""
ICON Cell-centered streamplots
==============================================

This example draws streamlines from ICON cell-centered vector winds with
:py:func:`easyclimate.plot.icon.plot_cell_streamplot <easyclimate.plot.icon.plot_cell_streamplot>`.

The ICON vectors are interpolated to a regular grid inside the requested
longitude-latitude window before calling Matplotlib's streamplot machinery.
"""
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import easyclimate as ecl

from easyclimate.plot.icon import plot_cell_streamplot

# %%
# Open a sample ICON model-level file and select the horizontal wind
# components stored on cell centers.
data = ecl.open_tutorial_dataset("icon_native_ml_20080909T000000Z")
u = data.u
v = data.v
u

# %%
# Draw streamlines over a regional window on plain Matplotlib axes.
fig, ax = plt.subplots()
q = plot_cell_streamplot(
    data, u, v,
    ax=ax,

    lon_min=140,
    lon_max=220,
    lat_min=25,
    lat_max=70,
)

ax.set_title("JW Wave Wind Field\n(height = 19, plot_cell_streamplot)")

# %%
# Draw streamlines on a shifted PlateCarree map. The helper keeps the visible
# window aligned with the requested longitude-latitude bounds.
# sphinx_gallery_thumbnail_number = -1
fig, ax = plt.subplots(
    subplot_kw={"projection": ccrs.PlateCarree(180)}
)
q = plot_cell_streamplot(
    data, u, v,
    ax=ax,

    lon_min=140,
    lon_max=220,
    lat_min=25,
    lat_max=70,

    transform=ccrs.PlateCarree(),
)

ax.set_title("JW Wave Wind Field\n(height = 19, plot_cell_streamplot)")
ax.coastlines(resolution="110m", linewidth=0.6, color="b")
ax.gridlines(draw_labels=True, alpha=0)
