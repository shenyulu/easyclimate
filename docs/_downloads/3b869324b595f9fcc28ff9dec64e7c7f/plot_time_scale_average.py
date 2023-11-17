# -*- coding: utf-8 -*-
"""
Time Scale Average
===================================

Before proceeding with all the steps, first import some necessary libraries and packages
"""
import easyclimate as ecl
import xarray as xr

#%%
#
sic_data_Barents_Sea = ecl.tutorial.open_tutorial_dataset('mini_HadISST_ice').sic
sic_data_Barents_Sea

#%%
# Mean States
# ------------------------------------
ecl.calc_climatological_mean(sic_data_Barents_Sea, dim = 'time')

#%%
#
ecl.calc_climatological_seasonal_mean(sic_data_Barents_Sea, dim = 'time')

#%%
#
ecl.calc_seasonal_cycle_mean(sic_data_Barents_Sea, dim = 'time')

#%%
# Remove Seasonal Cycle
# ------------------------------------
sic_data_Barents_Sea_remove_seasonal_cycle = ecl.remove_seasonal_cycle_mean(sic_data_Barents_Sea, dim = 'time')
sic_data_Barents_Sea_remove_seasonal_cycle

#%%
#
sic_data_Barents_Sea_remove_seasonal_cycle.mean(dim = ('lat', 'lon')).sel(time = slice('2010-01-01', '2015-12-31')).plot(
    figsize = (10, 3),
    marker = '.',
)

#%%
# Convert to the Month-mean State corresponding to Each Month
# ------------------------------------
everymonth_climatology = ecl.transfer_monmean2everymonthmean(sic_data_Barents_Sea)
everymonth_climatology;

#%%
#
everymonth_climatology.mean(dim = ('lat', 'lon')).sel(time = slice('2010-01-01', '2015-12-31')).plot(
    figsize = (10, 3),
    marker = '.',
)