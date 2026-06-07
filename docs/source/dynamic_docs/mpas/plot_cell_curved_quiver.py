# -*- coding: utf-8 -*-
"""
MPAS Cell-centered curved quiver plots
==============================================

This example draws MPAS cell-centered vector winds as curved arrows with
:py:func:`easyclimate.plot.mpas.plot_cell_curved_quiver <easyclimate.plot.mpas.plot_cell_curved_quiver>`.

The helper interpolates MPAS cell vectors to a regular plotting grid within the
requested longitude-latitude window.
"""
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import easyclimate as ecl

from easyclimate.plot.mpas import plot_cell_curved_quiver

# %%
# Open MPAS output and select the reconstructed zonal and meridional wind
# components.
data = ecl.open_tutorial_dataset("mpas_JWwave_T10_nVertLevels10")
U = data.uReconstructZonal
V = data.uReconstructMeridional

# %%
# Draw curved arrows for a midlatitude North Pacific domain. The ``nx`` and
# ``ny`` arguments control the interpolation grid resolution.
fig, ax = plt.subplots(subplot_kw= {"projection": ccrs.PlateCarree(-120)})

cq = plot_cell_curved_quiver(
    data,
    U,
    V,

    lon_min=190,
    lon_max=-60,
    lat_min=10,
    lat_max=80,

    nx=40,
    ny=20,
    density = 2,
    color = 'k',
)

ax.set_title("JW Wave Wind Field (plot_cell_curved_quiver)")
ax.coastlines(resolution="110m", linewidth=0.6, color = "brown")
ax.gridlines(draw_labels=True, alpha = 0)

# %%
# The same field can be drawn on a polar projection. ``regrid_shape`` controls
# the Cartopy-side regridding used for the projected vectors.

# sphinx_gallery_thumbnail_number = -2
fig, ax = plt.subplots(
    figsize = (6, 6),
    subplot_kw={"projection": ccrs.NorthPolarStereo(central_longitude=-90)}
)

ax.coastlines(color="brown", linewidths=0.5)

cq = plot_cell_curved_quiver(
    data,
    U,
    V,

    lon_min=0,
    lon_max=359,
    lat_min=30,
    lat_max=90,

    nx=20,
    ny=18,

    density = 2,
    regrid_shape = 26,
    title="",
)

gl, meta = ecl.plot.draw_polar_basemap(
    ax = ax,
    lon_step=30,
    lat_step=20,
    lat_range=[30, 90],
    draw_labels=True,
    gridlines_kwargs={"color": "grey", "alpha": 0.2, "linestyle": "--"},
    lat_label_lon=-60
)

ecl.plot.set_polar_title("JW Wave Wind Field", meta, size = 15)
