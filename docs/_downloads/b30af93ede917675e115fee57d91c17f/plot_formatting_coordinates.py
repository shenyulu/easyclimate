# -*- coding: utf-8 -*-
"""
Formatting Coordinates
===================================

Before proceeding with all the steps, first import some necessary libraries and packages
"""
import easyclimate as ecl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import cartopy.crs as ccrs

#%%
# 
u_data = ecl.tutorial.open_tutorial_dataset('uwnd_202201_mon_mean').sortby('lat').uwnd
u_data

#%%
# 
draw_data1 = u_data.isel(time = 0).sel(level = 500)
draw_data1

#%%
# 
fig, ax  = plt.subplots(1, 1)

draw_data1.plot(
    ax = ax,
)

#%%
# 
fig, ax  = plt.subplots(1, 1)

draw_data1.plot(
    ax = ax,
)

ecl.plot.set_lon_format_axis(ax, axis = 'x')
ecl.plot.set_lat_format_axis(ax, axis = 'y')

#%%
# 
draw_data1_1 = u_data.isel(time = 0).sel(level = 500).sel(lon = slice(100, 110), lat = slice(20, 23))
draw_data1_1

#%%
# 
fig, ax = plt.subplots(1, 2, figsize = (15, 5))

for axi in ax.flat:
    draw_data1_1.plot(
        ax = axi,
        cmap = 'Reds'
    )
    axi.xaxis.set_major_locator(ticker.LinearLocator(5))

ecl.plot.set_lon_format_axis(ax[0], axis = 'x')
ecl.plot.set_lat_format_axis(ax[0], axis = 'y')

ecl.plot.set_lon_format_axis(ax[1], axis = 'x', dms = True)
ecl.plot.set_lat_format_axis(ax[1], axis = 'y', dms = True)

#%%
#
draw_data2 = u_data.isel(time = 0).mean(dim = 'lon')

#%%
#
fig, ax  = plt.subplots(1, 1)

draw_data2.plot.contourf(
    ax = ax,
    levels = 21,
    yincrease = False
)

#%%
#
fig, ax  = plt.subplots(1, 1)

draw_data2.plot.contourf(
    ax = ax,
    levels = 21,
    yincrease = False
)

ecl.plot.set_lat_format_axis(ax, axis = 'x')
ecl.plot.set_p_format_axis(ax, axis = 'y')

#%%
#
sic_data = ecl.tutorial.open_tutorial_dataset('mini_HadISST_ice').sic.isel(time = 10)
sic_data.plot.contourf(
    cmap = 'Blues',
    levels = 11
)

#%%
#
fig, ax = plt.subplots(subplot_kw = {'projection': ccrs.NorthPolarStereo()})

ax.coastlines(edgecolor = 'black', linewidths = 0.5)
ax.stock_img()

ecl.plot.draw_Circlemap_PolarStereo(
    ax = ax,
    lon_step = 30,
    lat_step = 10,
    lat_range = [50, 90],
    draw_labels = True,
    gridlines_kwargs = {'color': 'grey', 'alpha': 0.5, 'linestyle' : '--'}
)

sic_data.plot.contourf(
    cmap = 'Blues',
    levels = 11,
    transform = ccrs.PlateCarree()
)

#%%
#
fig, ax = plt.subplots(subplot_kw = {'projection': ccrs.NorthPolarStereo()})

ax.coastlines(edgecolor = 'black', linewidths = 0.5)
ax.stock_img()

ecl.plot.draw_Circlemap_PolarStereo(
    ax = ax, 
    lon_step = 30,
    lat_step = 10,
    lat_range = [50, 90],
    draw_labels = True,
    set_map_boundary_kwargs = {'north_pad': 0.3, 'south_pad': 0.4},
    gridlines_kwargs = {'color': 'grey', 'alpha': 0.5, 'linestyle' : '--'}
)

sic_data.plot(
    cmap = 'Blues',
    levels = 11,
    transform = ccrs.PlateCarree()
)