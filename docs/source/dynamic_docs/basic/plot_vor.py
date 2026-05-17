# -*- coding: utf-8 -*-
"""
Relative Vorticity
===================================

.. math::

    \\zeta = \\frac{\\partial v}{\\partial x} - \\frac{\\partial u}{\\partial y} + \\frac{u}{R} \\tan \\varphi

Before proceeding with all the steps, first import some necessary libraries and packages
"""
import easyclimate as ecl
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker

# %%
# Read sample data
udata = ecl.open_tutorial_dataset("uwnd_2022_day5")["uwnd"].sortby("lat")
vdata = ecl.open_tutorial_dataset("vwnd_2022_day5")["vwnd"].sortby("lat")
uvdata = xr.Dataset(data_vars = {"u": udata,"v": vdata})
uvdata

# %%
# Calculating vorticity using three types of functions
#
# - :py:func:`easyclimate.calc_vorticity <easyclimate.core.rvdv.calc_vorticity>`: Numpy Method.
# - :py:func:`easyclimate.calc_vorticity_ncl <easyclimate.core.rvdv.calc_vorticity_ncl>`: NCL Method.
# - :py:func:`easyclimate.calc_vorticity_rs <easyclimate.core.rvdv.calc_vorticity_rs>`: Rust Method (The calculation process is similar to that of NCL, but it is faster.).
vor_raw = ecl.calc_vorticity(
    uvdata.u,
    uvdata.v,
)

vor_ncl = ecl.calc_vorticity_ncl(
    uvdata.u,
    uvdata.v,
    cyclic_boundary_setting="nan"
)

vor_rust = ecl.calc_vorticity_rs(
    uvdata.u,
    uvdata.v,
    cyclic_boundary_setting="nan"
)
vor_rust

# %%
# Plot the results of the vorticity calculation
#
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker

lonlat_range = [70, 180, 0, 70]

proj = ccrs.PlateCarree()

fig, ax = plt.subplots(
    3, 1,
    figsize=(4.5, 9),
    subplot_kw={"projection": proj},
    constrained_layout=True
)

# Unify the color scale range
vmax = 1e-4
levels = np.linspace(-vmax, vmax, 21)
cmap = "RdBu_r"

# The three scenes that need to be drawn
fields = [
    (vor_raw.isel(time=0).sel(level=850),  "vor raw"),
    (vor_ncl.isel(time=0).sel(level=850),  "vor ncl"),
    (vor_rust.isel(time=0).sel(level=850), "vor rust"),
]

pcm = None

for i, (axi, (fld, title)) in enumerate(zip(ax.flat, fields)):
    axi.set_extent(lonlat_range, crs=ccrs.PlateCarree())

    # Coastlines and Base Maps
    axi.coastlines(resolution="50m", linewidth=0.8)

    # Grid lines
    gl = axi.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=0.5,
        color="gray",
        alpha=0.4,
        linestyle="--"
    )
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = True
    gl.bottom_labels = True

    # It looks cleaner if you keep the bottom labels on only the bottom row.
    if i < 2:
        gl.bottom_labels = False

    gl.xlocator = mticker.FixedLocator([80, 100, 120, 140, 160, 180])
    gl.ylocator = mticker.FixedLocator([0, 10, 20, 30, 40, 50, 60, 70])

    # Drawing
    pcm = fld.plot.contourf(
        ax=axi,
        transform=ccrs.PlateCarree(),
        levels=levels,
        cmap=cmap,
        add_colorbar=False,
        extend="both"
    )

    axi.set_title(title, fontsize=13)

# Shared colorbar
cbar = fig.colorbar(
    pcm,
    ax=ax,
    orientation="horizontal",
    shrink=0.9,
    pad=0.04,
    aspect=40
)
cbar.set_label("Vorticity [s$^{-1}$]")

# %%
# When comparing the results of different calculation methods, it is generally recommended to use the Rust-based method.
fig, ax = plt.subplots(2, 1, figsize = (5, 8), sharex=True)

(vor_ncl - vor_raw).isel(time = 0).sel(level = 850).plot(ax = ax[0])
ax[0].set_title("NCL - Numpy")
ax[0].set_xlabel("")
ax[0].set_ylabel("")

(vor_ncl - vor_rust).isel(time = 0).sel(level = 850).plot(ax = ax[1])
ax[1].set_title("NCL - Rust")
ax[1].set_xlabel("")
ax[1].set_ylabel("")
