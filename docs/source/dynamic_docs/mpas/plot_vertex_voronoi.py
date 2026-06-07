# -*- coding: utf-8 -*-
"""
MPAS Vertex-centered Voronoi plots
==============================================

This example draws MPAS vertex-centered scalar data on the dual polygons around
each vertex with
:py:func:`easyclimate.plot.mpas.plot_vertex_voronoi <easyclimate.plot.mpas.plot_vertex_voronoi>`.

If cell-centered data are supplied, the helper can average neighboring cell
values to vertices before plotting.
"""
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import easyclimate as ecl

from easyclimate.plot.mpas import plot_vertex_voronoi

# %%
# Open the sample MPAS output and use a selected relative-vorticity field for the
# vertex-based examples.
data = ecl.open_tutorial_dataset("mpas_JWwave_T10_nVertLevels10")
data

# %%
# Draw vertex dual polygons over a broad region.
plot_vertex_voronoi(
    data,
    data.vorticity,
    lon_min=130,
    lon_max=320,
    lat_min=10,
    lat_max=80,
)

# %%
# The helper can also draw on an existing Cartopy axes.
fig, ax = plt.subplots(subplot_kw= {"projection": ccrs.PlateCarree(-120)})

plot_vertex_voronoi(
    data,
    data.vorticity,

    ax = ax,
    transform = ccrs.PlateCarree(),

    lon_min=130,
    lon_max=320,
    lat_min=10,
    lat_max=80,

    cbar_kwargs = {'location': 'bottom', 'aspect': 60}
)

ax.set_title("JW Wave Relative Vorticity (plot_vertex_voronoi)")
ax.coastlines(resolution="110m", linewidth=0.6, color = "b")
ax.gridlines(draw_labels=True, alpha = 0)

# %%
# Show polygon edges in a smaller region to inspect the dual mesh structure.

# sphinx_gallery_thumbnail_number = -2
fig, ax = plt.subplots(subplot_kw= {"projection": ccrs.PlateCarree(160)})

plot_vertex_voronoi(
    data,
    data.vorticity,

    transform = ccrs.PlateCarree(),
    ax = ax,

    lon_min=200,
    lon_max=270,
    lat_min=30,
    lat_max=70,

    linewidth = 0.1,
    edgecolor="grey",

    cbar_kwargs = {'location': 'bottom', 'aspect': 60}
)

ax.set_title("Local JW Wave Relative Vorticity (plot_cell_voronoi)")
ax.coastlines(resolution="50m", linewidth=0.6, color = "r")
ax.gridlines(draw_labels=True, alpha = 0)
