# -*- coding: utf-8 -*-
"""
ICON Cell-centered quiver plots
==============================================

This example draws ICON cell-centered wind vectors with
:py:func:`easyclimate.plot.icon.plot_cell_quiver <easyclimate.plot.icon.plot_cell_quiver>`.

The helper thins irregular ICON cell vectors into longitude-latitude bins
before passing them to Matplotlib's quiver machinery.
"""
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import easyclimate as ecl

from easyclimate.plot.icon import plot_cell_quiver

# %%
# Open a sample ICON model-level file and select the horizontal wind
# components stored on cell centers.
data = ecl.open_tutorial_dataset("icon_native_ml_20080909T000000Z")
u = data.u
v = data.v
u

# %%
# Draw thinned vectors over a regional window on plain Matplotlib axes.
fig, ax = plt.subplots()
q = plot_cell_quiver(
    data, u, v,
    ax=ax,

    lon_min=140,
    lon_max=220,
    lat_min=25,
    lat_max=70,

    nx_bins=30,
    ny_bins=20,
    thin_method="nearest_center",
    scale=900,
    width = 0.004,
    headaxislength = 5,
)

ax.set_title("JW Wave Wind Field\n(height = 19, plot_cell_quiver)")

# %%
# Draw the same vectors on a shifted PlateCarree map. The ``transform`` keyword
# identifies the input wind locations as geographic coordinates.
# sphinx_gallery_thumbnail_number = -1
fig, ax = plt.subplots(
    subplot_kw={"projection": ccrs.PlateCarree(180)}
)
q = plot_cell_quiver(
    data, u, v,
    ax=ax,

    lon_min=140,
    lon_max=220,
    lat_min=25,
    lat_max=70,

    nx_bins=30,
    ny_bins=20,
    thin_method="nearest_center",
    scale=900,
    width = 0.004,
    headaxislength = 5,
    transform=ccrs.PlateCarree(),
)

ax.set_title("JW Wave Wind Field\n(height = 19, plot_cell_quiver)")
ax.coastlines(resolution="110m", linewidth=0.6, color="b")
ax.gridlines(draw_labels=True, alpha=0)
