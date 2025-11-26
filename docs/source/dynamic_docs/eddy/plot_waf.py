# -*- coding: utf-8 -*-
"""
.. _max_waf_example:

Wave Activity Horizontal Flux
============================================

In this tutorial, we'll explore how to visualize and analyze Rossby waves - those giant meanders in high-altitude winds that shape our weather patterns.
Think of them as the "weather rivers" flowing through our atmosphere!
"""
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import easyclimate as ecl
import matplotlib.ticker as ticker

formatter = ticker.ScalarFormatter(useMathText=True, useOffset=True)
formatter.set_powerlimits((0, 0))

# %%
z500_prime_data1 = xr.open_dataset("era5_daily_z500_prime_201411_N15.nc").z /9.8
z500_prime_data2 = xr.open_dataset("era5_daily_z500_prime_202411_N15.nc").z /9.8
u500_climatology_data = xr.open_dataset("era5_ymean_monthly_u500_199101_202012_N15.nc").u
v500_climatology_data = xr.open_dataset("era5_ymean_monthly_v500_199101_202012_N15.nc").v

# %%
# TN01 Wave Activity Horizontal Flux
# -----------------------------------------
tn01_result1 = ecl.calc_TN_wave_activity_horizontal_flux(
    z_prime_data = z500_prime_data1,
    u_climatology_data = u500_climatology_data,
    v_climatology_data = v500_climatology_data,
    vertical_dim = "level",
    vertical_dim_units = "hPa",
)
tn01_result1

# %%
tn01_mean_result1 = tn01_result1.mean(dim = "time").sortby("lat")
tn01_mean_result1

# %%
draw_shaded1 = tn01_mean_result1["psi_p"]
draw_shaded1 = ecl.plot.add_lon_cyclic(draw_shaded1, 6)
draw_quiver1 = tn01_mean_result1[["fx", "fy"]]

# %%
fig, ax = plt.subplots(
    figsize = (6, 6),
    subplot_kw={"projection": ccrs.NorthPolarStereo(central_longitude=140)}
)

ax.coastlines(edgecolor="black", linewidths=0.5)
ecl.plot.draw_Circlemap_PolarStereo(
    ax=ax,
    lon_step=30,
    lat_step=30,
    lat_range=[10, 90],
    draw_labels=True,
    gridlines_kwargs={"color": "grey", "alpha": 0.5, "linestyle": "--"},
)

fg1 = draw_shaded1.plot.contourf(
    ax = ax,
    levels = np.linspace(-1.5e7, 1.5e7, 21),
    cbar_kwargs = {'location': 'bottom', 'aspect': 40, 'format': formatter},
    transform=ccrs.PlateCarree(),
    zorder = 1
)

q = draw_quiver1.sel(lat = slice(10, None)).plot.quiver(
    ax = ax,
    x = "lon", y = "lat", u = "fx", v = "fy",
    transform = ccrs.PlateCarree(),
    scale = 600,
    regrid_shape = 18,
    headwidth = 5,
    width = 0.0045,
    headlength = 8,
    headaxislength = 8,
    minlength = 3,
    add_guide = False,
    zorder = 2,
)
qk = ax.quiverkey(
    q, 0, -0.1,
    50, "50", labelpos = "N",
    color = "k",
    coordinates = "axes",
    fontproperties = {"size": 12},
    labelsep = 0.05,
    transform = ccrs.PlateCarree(),
)
qk.set_zorder(3)

ax.set_title("TN01 WAF Analysis (500hPa, Nov. 2014)")

# %%
tn01_result2 = ecl.calc_TN_wave_activity_horizontal_flux(
    z_prime_data = z500_prime_data2,
    u_climatology_data = u500_climatology_data,
    v_climatology_data = v500_climatology_data,
    vertical_dim = "level",
    vertical_dim_units = "hPa",
)
tn01_mean_result2 = tn01_result2.mean(dim = "time").sortby("lat")
tn01_mean_result2

# %%
draw_shaded2 = tn01_mean_result2["psi_p"]
draw_shaded2 = ecl.plot.add_lon_cyclic(draw_shaded2, 6)
draw_quiver2 = tn01_mean_result2[["fx", "fy"]]

# %%
fig, ax = plt.subplots(
    figsize = (12, 6),
    subplot_kw={"projection": ccrs.Mercator(central_longitude=140)}
)
ax.set_extent([0, 360, 10, 85], crs = ccrs.PlateCarree())
ax.coastlines(edgecolor="black", linewidths=0.5)
ax.gridlines(draw_labels=["bottom", "left"], alpha=0)

fg1 = draw_shaded2.plot.contourf(
    ax = ax,
    levels = np.linspace(-1e7, 1e7, 21),
    cbar_kwargs = {'location': 'bottom', 'aspect': 80, 'shrink': 0.7, 'pad': 0.1, 'format': formatter},
    transform=ccrs.PlateCarree(),
    zorder = 1
)

q = draw_quiver2.sel(lat = slice(5, None)).plot.quiver(
    ax = ax,
    x = "lon", y = "lat", u = "fx", v = "fy",
    transform = ccrs.PlateCarree(),
    scale = 800,
    regrid_shape = 15,
    add_guide = False,
)
qk = ax.quiverkey(
    q, 0.95, 1.03,
    50, "50", labelpos = "N",
    color = "k",
    coordinates = "axes",
    fontproperties = {"size": 12},
    labelsep = 0.05,
    transform = ccrs.PlateCarree(),
)

ax.set_title("TN01 WAF Analysis (500hPa, Nov. 2024)")

# %%
# Plumb Wave Activity Horizontal Flux
# -----------------------------------------
pwaf_result = ecl.calc_Plumb_wave_activity_horizontal_flux(
    z_prime_data = z500_prime_data2,
    vertical_dim = "level",
    vertical_dim_units = "hPa",
)
pwaf_mean_result = pwaf_result.mean(dim = "time").sortby("lat")
pwaf_mean_result

# %%
draw_shaded3 = pwaf_mean_result["psi_p"]
draw_shaded3 = ecl.plot.add_lon_cyclic(draw_shaded3, 6)
draw_quiver3 = pwaf_mean_result[["fx", "fy"]]

# %%
fig, ax = plt.subplots(
    figsize = (12, 6),
    subplot_kw={"projection": ccrs.Mercator(central_longitude=140)}
)
ax.set_extent([0, 360, 10, 85], crs = ccrs.PlateCarree())
ax.coastlines(edgecolor="black", linewidths=0.5)
ax.gridlines(draw_labels=["bottom", "left"], alpha=0)

fg1 = draw_shaded3.plot.contourf(
    ax = ax,
    levels = np.linspace(-1e7, 1e7, 21),
    cbar_kwargs = {'location': 'bottom', 'aspect': 80, 'shrink': 0.7, 'pad': 0.1, 'format': formatter},
    transform=ccrs.PlateCarree(),
    zorder = 1
)

q = draw_quiver3.sel(lat = slice(5, None)).plot.quiver(
    ax = ax,
    x = "lon", y = "lat", u = "fx", v = "fy",
    transform = ccrs.PlateCarree(),
    scale = 800,
    regrid_shape = 15,
    add_guide = False,
)
qk = ax.quiverkey(
    q, 0.95, 1.03,
    50, "50", labelpos = "N",
    color = "k",
    coordinates = "axes",
    fontproperties = {"size": 12},
    labelsep = 0.05,
    transform = ccrs.PlateCarree(),
)

ax.set_title("Plumb WAF Analysis (500hPa, Nov. 2024)")
