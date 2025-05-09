# -*- coding: utf-8 -*-
"""
Onset and Retreat of Monsoon
=========================================================================================================

Before proceeding with all the steps, first import some necessary libraries and packages
"""
import xarray as xr
import cartopy.crs as ccrs
import easyclimate as ecl

# %%
# Open precipitable water example data (must be **daily frequency** data)
pw_data = ecl.open_tutorial_dataset("pr_wtr_eatm_2022").pr_wtr
pw_data

# %%
# Calculate the NPWI index using :py:func:`easyclimate.field.monsoon.calc_index_NPWI <easyclimate.field.monsoon.calc_index_NPWI>`
# .. seealso::
#
#     - Zeng, X., and E. Lu, 2004: Globally Unified Monsoon Onset and Retreat Indexes. J. Climate, 17, 2241–2248, https://journals.ametsoc.org/view/journals/clim/17/11/1520-0442_2004_017_2241_gumoar_2.0.co_2.xml.
#     - Tang Xu, Chen Baode, Liang Ping, Qian Weihong. Definition and features of the north edge of Asian summer monsoon. Acta Meteorologica Sinica (Chinese), 2009, (1): 83-89. doi: http://dx.doi.org/10.11676/qxxb2009.009
NPWI_index = ecl.field.monsoon.calc_index_NPWI(pw_data)
NPWI_index

# %%
# Separation of monsoon affected areas with :py:func:`easyclimate.field.monsoon.find_PW_monsoon_region <easyclimate.field.monsoon.find_PW_monsoon_region>`
PW_monsoon_region = ecl.field.monsoon.find_PW_monsoon_region(pw_data)
PW_monsoon_region

# %%
# Schematization of monsoon impact areas
fig, ax = ecl.plot.quick_draw_spatial_basemap(central_longitude = 200)
PW_monsoon_region.plot(
    ax = ax,
    cbar_kwargs = {'location': 'bottom'},
    transform = ccrs.PlateCarree(),
    levels = [0, 1],
    colors = ["grey"]
)

# %%
# Calculation of monsoon onset with :py:func:`easyclimate.field.monsoon.calc_NPWI_monsoon_onset <easyclimate.field.monsoon.calc_NPWI_monsoon_onset>`
monsoon_onset_date = ecl.field.monsoon.calc_NPWI_monsoon_onset(NPWI_index)
monsoon_onset_date

# %%
# Analyzing and mapping the monsoon onset time in the monsoon impact area
fig, ax = ecl.plot.quick_draw_spatial_basemap(central_longitude = 180)
monsoon_onset_date.where(PW_monsoon_region).plot(
    ax = ax,
    cbar_kwargs = {'location': 'bottom'},
    transform = ccrs.PlateCarree(),
    cmap = 'Reds'
)

# %%
# Calculation of monsoon retreat time with :py:func:`easyclimate.field.monsoon.calc_NPWI_monsoon_retreat <easyclimate.field.monsoon.calc_NPWI_monsoon_retreat>`
monsoon_retreat_date = ecl.field.monsoon.calc_NPWI_monsoon_retreat(NPWI_index, monsoon_onset_date)
monsoon_retreat_date

# %%
# Analyzing and mapping the monsoon retreat time in the monsoon impact area

# sphinx_gallery_thumbnail_number = -1
fig, ax = ecl.plot.quick_draw_spatial_basemap(central_longitude = 180)
monsoon_retreat_date.where(PW_monsoon_region).plot(
    ax = ax,
    cbar_kwargs = {'location': 'bottom'},
    transform = ccrs.PlateCarree(),
    cmap = 'Reds'
)
