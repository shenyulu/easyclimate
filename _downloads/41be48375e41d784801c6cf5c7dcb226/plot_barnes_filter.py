# -*- coding: utf-8 -*-
"""
Barnes Filter
===================================

Barnes filter is a commonly used spatial filtering method that mainly uses two constants g and c to calculate Gaussian weights,
and performs spatial interpolation for each grid point, thus becoming a low-pass filter that filters out high-frequency fluctuations.
When using two different schemes of constant g and c schemes, both retain low-frequency fluctuations of different scales.
The difference between the filtering results of the two methods can result in mesoscale fluctuations.

.. raw:: html

    <iframe width="100%" height="422" src="https://player.bilibili.com/player.html?isOutside=true&aid=569810777&bvid=BV1tv4y1H79z&cid=1089298762&p=1" scrolling="no" border="0" frameborder="no" framespacing="0" allowfullscreen="true"></iframe>

.. seealso::

    - Maddox, R. A. (1980). An Objective Technique for Separating Macroscale and Mesoscale Features in Meteorological Data. Monthly Weather Review, 108(8), 1108-1121. https://journals.ametsoc.org/view/journals/mwre/108/8/1520-0493_1980_108_1108_aotfsm_2_0_co_2.xml
    - https://github.com/LinOuyang/pybarnes

Before proceeding with all the steps, first import some necessary libraries and packages
"""
import numpy as np
import cartopy.crs as ccrs
import easyclimate as ecl

# %%
# Open tuturial dataset
hgt_data = ecl.open_tutorial_dataset("hgt_2022_day5").hgt.sel(level = 1000)
hgt_data

# %%
# Filter dataset using :py:func:`easyclimate.filter.calc_barnes_lowpass <easyclimate.filter.calc_barnes_lowpass>`
hgt_data1 = ecl.filter.calc_barnes_lowpass(hgt_data)
hgt_data1

# %%
# Draw results and differences.
fig, ax = ecl.plot.quick_draw_spatial_basemap(3, 1, figsize=(5, 12), central_longitude=180)

axi = ax[0]
hgt_data.isel(time = 0).plot.contourf(
    ax = axi, levels = 21,
    transform = ccrs.PlateCarree(),
    cbar_kwargs={"location": "bottom", "aspect": 30, "label": ""},
)
axi.set_title("Raw data")

axi = ax[1]
hgt_data1.isel(time = 0).plot.contourf(
    ax = axi, levels = 21,
    transform = ccrs.PlateCarree(),
    cbar_kwargs={"location": "bottom", "aspect": 30, "label": ""},
)
axi.set_title("Filtered data")

axi = ax[2]
draw_dta = hgt_data.isel(time = 0) - hgt_data1.isel(time = 0)
draw_dta.plot.contourf(
    ax = axi, levels = np.linspace(-30, 30, 21),
    transform = ccrs.PlateCarree(),
    cbar_kwargs={"location": "bottom", "aspect": 30, "label": ""},
)
axi.set_title("Difference: Mesoscale fluctuations")
