# -*- coding: utf-8 -*-
"""
Time Scale Average
===================================

Before proceeding with all the steps, first import some necessary libraries and packages
"""
import easyclimate as ecl
import xarray as xr

# %%
# Example data are sea ice concentration (SIC) data for the Barents-Kara Sea (30°-90°E, 65°-85°N).
sic_data_Barents_Sea = ecl.tutorial.open_tutorial_dataset("mini_HadISST_ice").sic
sic_data_Barents_Sea

# %%
# Mean States
# ------------------------------------
# Solving for the overall climatological mean state was solved using :py:func:`easyclimate.calc_all_climatological_mean <easyclimate.calc_all_climatological_mean>`.
#
ecl.calc_all_climatological_mean(sic_data_Barents_Sea, dim="time")

# %%
# If the climate state is for each season, the results are solved using :py:func:`easyclimate.calc_seasonal_climatological_mean <easyclimate.calc_seasonal_climatological_mean>`.
ecl.calc_seasonal_climatological_mean(sic_data_Barents_Sea, dim="time")

# %%
# However, if the climate state is for each month, the results are solved using :py:func:`easyclimate.calc_seasonal_cycle_mean <easyclimate.calc_seasonal_cycle_mean>`.
ecl.calc_seasonal_cycle_mean(sic_data_Barents_Sea, dim="time")

# %%
# Remove Seasonal Cycle
# ------------------------------------
# :py:func:`easyclimate.remove_seasonal_cycle_mean <easyclimate.remove_seasonal_cycle_mean>` helps us to remove seasonal cycles (annual cycles) from the data in order to obtain monthly average anomalies.
sic_data_Barents_Sea_remove_seasonal_cycle = ecl.remove_seasonal_cycle_mean(
    sic_data_Barents_Sea, dim="time"
)
sic_data_Barents_Sea_remove_seasonal_cycle

# %%
# We can visualize the results by plotting.
sic_data_Barents_Sea_remove_seasonal_cycle.mean(dim=("lat", "lon")).sel(
    time=slice("2010-01-01", "2015-12-31")
).plot(
    figsize=(10, 3),
    marker=".",
)
