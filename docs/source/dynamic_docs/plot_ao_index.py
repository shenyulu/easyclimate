# -*- coding: utf-8 -*-
"""
Arctic Oscillation (AO) Index
===================================

The Arctic Oscillation (AO) Index (or Monthly Northern Hemisphere Annular Mode (NAM) Index) is a key metric used to describe
large-scale atmospheric variability in the Northern Hemisphere,
particularly influencing mid-to-high latitude weather patterns. It is defined by the leading mode of Empirical Orthogonal Function (EOF)
analysis of sea-level pressure (SLP) anomalies north of 20°N. The AO Index quantifies fluctuations
in atmospheric pressure between the Arctic and mid-latitudes, with positive and negative phases reflecting
distinct circulation patterns. In the positive phase, lower Arctic pressure and higher mid-latitude
pressure strengthen westerly winds, confining cold air to polar regions, often leading to milder
winters in North America and Europe. The negative phase, with higher Arctic pressure and weaker winds,
allows cold air to move southward, causing colder, stormier weather in these regions.

The AO's role in climate variability is significant, as it modulates temperature and precipitation, especially in winter.
The AO Index, typically derived from monthly or seasonal SLP data, reflects the strength of the polar vortex,
with positive values indicating a stronger vortex and negative values a weaker one.
It is closely linked to the North Atlantic Oscillation (NAO) due to shared variability patterns.

The AO's fluctuations are driven by internal atmospheric dynamics, stratospheric processes,
and external forcings like sea surface temperatures. Its teleconnections make it a critical factor in seasonal weather predictions and long-term climate modeling.
In a warming climate, Arctic amplification may alter AO dynamics, making its study essential for understanding future climate trends.

.. seealso::

    - Thompson, D. W. J., & Wallace, J. M. (1998). The Arctic oscillation signature in the wintertime geopotential height and temperature fields. Geophysical Research Letters, 25(9), 1297–1300. https://doi.org/10.1029/98gl00950
    - Fang, Z., Sun, X., Yang, X.-Q., & Zhu, Z. (2024). Interdecadal variations in the spatial pattern of the Arctic Oscillation Arctic center in wintertime. Geophysical Research Letters, 51, e2024GL111380. https://doi.org/10.1029/2024GL111380
    - Li, J., and J. X. L. Wang (2003), A modified zonal index and its physical sense, Geophys. Res. Lett., 30, 1632, doi: https://doi.org/10.1029/2003GL017441, 12.
    - Thompson, D. W. J. , & Wallace, J. M. . (1944). The Arctic oscillation signature in the wintertime geopotential height and temperature fields. Geophys. Res. Lett., doi: https://10.1029/98GL00950

Before proceeding with all the steps, first import some necessary libraries and packages
"""
import easyclimate as ecl
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# %%
# Load monthly mean sea level pressure (SLP) data for the Northern Hemisphere
# The data is read from a NetCDF file containing SLP monthly means
slp_data = xr.open_dataset("slp_monmean_NH.nc").slp

# %%
# Calculate seasonal mean for December-January-February (DJF)
# This aggregates the monthly data into winter seasonal means
slp_data_DJF_mean = ecl.calc_seasonal_mean(slp_data, extract_season = 'DJF')
slp_data_DJF_mean

# %%
# Remove the seasonal cycle mean to obtain anomalies
# This creates anomalies by subtracting the long-term seasonal mean
slp_data_DJF_anormaly = ecl.remove_seasonal_cycle_mean(slp_data_DJF_mean)
slp_data_DJF_anormaly

# %%
# Calculate Arctic Oscillation (AO) index using EOF method (Thompson & Wallace 1998)
# This performs EOF analysis on SLP anomalies north of 20°N to derive the AO index
index_ao = ecl.field.teleconnection.calc_index_AO_EOF_Thompson_Wallace_1998(slp_data_DJF_mean)
index_ao

# %%
# Calculate AO index using zonal mean SLP difference method (Li & Wang 2003)
# This computes the index as the normalized difference between 35°N and 65°N SLP
index_ao_point = ecl.field.teleconnection.calc_index_NAH_zonal_lat_Li_Wang_2003(slp_data_DJF_mean)
index_ao_point

# %%
# Calculate correlation between the two AO index calculation methods
# Shows how well the two different methods agree in capturing AO variability
np.corrcoef(index_ao_point, index_ao)

# %%
# Apply Gaussian filter to smooth the AO index time series
# Uses a 9-month window to highlight lower-frequency variations
index_ao_filtered = ecl.filter.calc_gaussian_filter(index_ao, window_length=9)
index_ao_filtered

# %%
# Perform linear regression between SLP anomalies and AO index
# Calculates the spatial pattern of SLP associated with AO variability
slp_reg_ao = ecl.calc_linregress_spatial(slp_data_DJF_anormaly, x = index_ao)
slp_reg_ao_rvalue = slp_reg_ao.rvalue
slp_reg_ao_pvalue = slp_reg_ao.pvalue

# %%
# Add cyclic point for plotting (avoids gap at 0/360° longitude)
slp_reg_ao_rvalue = ecl.plot.add_lon_cyclic(slp_reg_ao_rvalue, inter = 2.5)
slp_reg_ao_pvalue = ecl.plot.add_lon_cyclic(slp_reg_ao_pvalue, inter = 2.5)

# %%
# Prepare AO index data for plotting
# Converts time to year-only format for cleaner visualization
index_ao_bar = index_ao.copy(deep=True)
index_ao_bar['time'] = index_ao_bar['time'].dt.year.data

index_ao_filtered_bar = index_ao_filtered.copy(deep=True)
index_ao_filtered_bar['time'] = index_ao_filtered_bar['time'].dt.year.data

# %%
# Create time series plot of AO index
fig, ax = plt.subplots(figsize = (10, 4))

# Plot bar chart of annual AO index values
ecl.plot.bar_plot_with_threshold(index_ao_bar, ax = ax)
# Overlay smoothed time series
index_ao_filtered_bar.plot(color = 'k', lw = 2, ax = ax)

ax.set_title("Arctic Oscillation (AO) Index")

# %%
# Create spatial map of SLP regression pattern associated with AO
fig, ax = plt.subplots(subplot_kw={"projection": ccrs.NorthPolarStereo()})

# Configure map appearance
ax.coastlines(edgecolor="black", linewidths=0.5)

# Add polar map elements (grid lines, labels)
ecl.plot.draw_Circlemap_PolarStereo(
    ax=ax,
    lon_step=30,
    lat_step=10,
    lat_range=[20, 90],
    draw_labels=True,
    set_map_boundary_kwargs={"north_pad": 0.3, "south_pad": 0.4},
    gridlines_kwargs={"color": "grey", "alpha": 0.5, "linestyle": "--"},
)

# Plot regression coefficients (SLP pattern)
slp_reg_ao_rvalue.plot.contourf(
    cmap="RdBu_r",
    levels=11,
    transform=ccrs.PlateCarree(),
    cbar_kwargs = {'location': 'bottom'}
)
# Highlight statistically significant areas
ecl.plot.draw_significant_area_contourf(
    slp_reg_ao_pvalue,
    transform=ccrs.PlateCarree(),
    hatches = ".."
)
