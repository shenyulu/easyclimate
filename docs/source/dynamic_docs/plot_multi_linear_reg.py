# -*- coding: utf-8 -*-
"""
Multiple Variable Linear Regression
=========================================================================================================

This documentation demonstrates a comprehensive analysis of sea surface temperature (SST) variability
explained by two climate indices (AO and Niño 3.4) using multiple linear regression. The analysis
covers data preparation, index calculation, regression modeling, and visualization of results.

.. math::

    y = a_1 x_1 + a_2 x_2

Before proceeding with all the steps, first import some necessary libraries and packages
"""
import xarray as xr
import cartopy.crs as ccrs
import easyclimate as ecl

# %%
# Two time ranges are defined to account for seasonal analysis with different base periods
#
# - ``time_range``: Primary analysis period (1982-2020)
# - ``time_range_plus1``: Offset by one year for seasonal calculations
time_range = slice("1982-01-01", "2020-12-31")
time_range_plus1 = slice("1983-01-01", "2021-12-31")

# %%
# The Arctic Oscillation index is calculated using:
#
# 1. Sea level pressure (SLP) data from the Northern Hemisphere
# 2. Seasonal mean calculation for December-January-February (DJF)
# 3. EOF analysis following Thompson & Wallace (1998) methodology
#
# The resulting index is then subset to our analysis period.
slp_data = xr.open_dataset("slp_monmean_NH.nc").slp
slp_data_DJF_mean = ecl.calc_seasonal_mean(slp_data, extract_season = 'DJF')
index_ao = ecl.field.teleconnection.calc_index_AO_EOF_Thompson_Wallace_1998(slp_data_DJF_mean)
index_ao = index_ao.sel(time = time_range)
index_ao

# %%
# The Nino3.4 index is derived from:
#
# 1. Hadley Centre SST dataset
# 2. Seasonal mean for DJF period
# 3. Area averaging over the Niño 3.4 region (5°N-5°S, 170°W-120°W)
sst_data = ecl.open_tutorial_dataset("mini_HadISST_sst").sst
sst_data_DJF_mean = ecl.calc_seasonal_mean(sst_data, extract_season = 'DJF')
index_nino34 = ecl.field.air_sea_interaction.calc_index_nino34(sst_data_DJF_mean).sel(time = time_range)
index_nino34

# %%
# The dependent variable for our regression is prepared as:
#
# - Seasonal mean SST for September-October-November (SON)
# - Using the offset time range to examine potential lagged relationships
sst_data_SON_mean = ecl.calc_seasonal_mean(sst_data, extract_season = 'SON').sel(time = time_range_plus1)
sst_data_SON_mean

# %%
# The core analysis applies multiple linear regression to quantify how:
#
# - AO index (first predictor)
# - Niño 3.4 index (second predictor)
#
# jointly explain spatial patterns of SON SST variability.
#
# The function returns a dataset containing:
#
# - Regression coefficients (slopes) for each predictor
# - Intercept values
# - R-squared values (goodness of fit)
# - Statistical significance (p-values) for each parameter
result = ecl.calc_multiple_linear_regression_spatial(sst_data_SON_mean, [index_ao, index_nino34])
result

# %%
# The final visualization shows:
#
# - Top panel: Spatial pattern of AO influence on SON SST: Colors show regression coefficients, and Contours indicate statistically significant areas (p < 0.05)
# - Bottom panel: Spatial pattern of Niño 3.4 influence: Similar interpretation as top panel, and the central longitude is set to 200° for Pacific-centric viewing.
#
# Key interpretation points:
#
# - Positive coefficients indicate SST increases with positive phase of the index
# - Negative coefficients indicate inverse relationships
# - Non-significant areas suggest no robust statistical relationship
fig, ax = ecl.plot.quick_draw_spatial_basemap(nrows=2 ,figsize = (10, 5), central_longitude=200)

result.slopes.sel(coef = 0).plot(
    ax=ax[0],
    transform=ccrs.PlateCarree(),
    cbar_kwargs={"location": "bottom", "pad": 0.2, "aspect": 100, "shrink": 0.8},
)
ecl.plot.draw_significant_area_contourf(
    result.slopes_p.sel(coef = 0),
    ax = ax[0],
    transform=ccrs.PlateCarree()
)

result.slopes.sel(coef = 1).plot(
    ax=ax[1],
    transform=ccrs.PlateCarree(),
    cbar_kwargs={"location": "bottom", "pad": 0.2, "aspect": 100, "shrink": 0.8},
)
ecl.plot.draw_significant_area_contourf(
    result.slopes_p.sel(coef = 1),
    ax = ax[1],
    transform=ccrs.PlateCarree()
)
