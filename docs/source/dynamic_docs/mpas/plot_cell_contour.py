# -*- coding: utf-8 -*-
"""
MPAS Cell-centered contour plots
==============================================

This example demonstrates filled and line contours for MPAS cell-centered
scalar fields.

:py:func:`easyclimate.plot.mpas.plot_cell_contourf <easyclimate.plot.mpas.plot_cell_contourf>`
draws filled contours from values located on ``nCells``. The companion
:py:func:`easyclimate.plot.mpas.plot_cell_contour <easyclimate.plot.mpas.plot_cell_contour>`
draws contour lines using the same MPAS cell coordinates.
"""
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import easyclimate as ecl

from easyclimate.plot.mpas import plot_cell_contourf, plot_cell_contour

# %%
# Open a sample MPAS output file. The ``divergence`` variable is defined at
# MPAS cell centers, so one time and vertical level are selected before
# plotting.
data = ecl.open_tutorial_dataset("mpas_JWwave_T10_nVertLevels10")
data

# %%
# Draw a filled contour plot over a longitude-latitude subset. The helper uses
# the MPAS cell coordinates and supports windows that cross the dateline.
plot_cell_contourf(
    data,
    data.divergence,
    cmap = "RdBu_r",
    vmax = 3e-6,
    vmin = -3e-6,
    levels = 21,

    lon_min=130,
    lon_max=320,
    lat_min=10,
    lat_max=80,

    cbar_kwargs = {'location': 'bottom', 'aspect': 60}
)

# %%
# Pass an existing Cartopy axes when the contour should be combined with coastlines,
# gridlines, and other map decorations.
fig, ax = plt.subplots(subplot_kw= {"projection": ccrs.PlateCarree(-120)})

plot_cell_contourf(
    data,
    data.divergence,
    transform = ccrs.PlateCarree(),
    cmap = "RdBu_r",
    vmax = 3e-6,
    vmin = -3e-6,
    levels = 21,
    ax = ax,

    lon_min=130,
    lon_max=320,
    lat_min=10,
    lat_max=80,

    cbar_kwargs = {'location': 'bottom', 'aspect': 60}
)

ax.set_title("JW Wave Divergence (plot_cell_contourf)")
ax.coastlines(resolution="110m", linewidth=0.6, color = "b")
ax.gridlines(draw_labels=True, alpha = 0)

# %%
# A smaller longitude-latitude window can be used to inspect regional structures
# without pre-slicing the MPAS mesh manually.
fig, ax = plt.subplots(subplot_kw= {"projection": ccrs.PlateCarree(160)})

plot_cell_contourf(
    data,
    data.divergence,
    transform = ccrs.PlateCarree(),
    cmap = "RdBu_r",
    levels = np.linspace(-1e-6, 1e-6, 21),
    ax = ax,

    lon_min=120,
    lon_max=200,
    lat_min=20,
    lat_max=70,

    cbar_kwargs = {'location': 'bottom', 'aspect': 60}
)

ax.set_title("Local JW Wave Divergence (plot_cell_contourf)")
ax.coastlines(resolution="50m", linewidth=0.6, color = "k")
ax.gridlines(draw_labels=True, alpha = 0)

# %%
# Use ``plot_cell_contour`` when contour lines and labels are preferred over a
# filled field.

# sphinx_gallery_thumbnail_number = -2
fig, ax = plt.subplots(subplot_kw= {"projection": ccrs.PlateCarree(120)})

cs = plot_cell_contour(
    data,
    data.divergence,
    levels=np.linspace(-1e-6, 1e-6, 9),
    colors="black",
    ax=ax,
    lon_min=150,
    lon_max=230,
    lat_min=10,
    lat_max=65,
    transform=ccrs.PlateCarree(),
)

# Matplotlib may emit empty/degenerate contour paths for some levels on
# triangulated MPAS meshes, and clabel() can fail on those paths.
has_labelable_segments = any(
    len(segment) > 1 for level_segments in cs.allsegs for segment in level_segments
)
if has_labelable_segments:
    ax.clabel(cs, inline=True, fontsize=8, fmt="%g")

ax.set_title("Local JW Wave Divergence (plot_cell_contour)")
ax.coastlines(resolution="50m", linewidth=0.6, color = "r")
ax.gridlines(draw_labels=True, alpha = 0)
