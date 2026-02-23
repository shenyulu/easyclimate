# -*- coding: utf-8 -*-
"""
Spatial Detrend
===================================

Before proceeding with all the steps, first import some necessary libraries and packages
"""
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
import easyclimate as ecl

# %%
#
# The significance of removing spatial trends in meteorological data is profound,
# as it allows scientists to isolate and analyze the specific atmospheric phenomena
# they are interested in without the obfuscating influence of large-scale,
# systematic patterns. In essence, it is a fundamental pre-processing
# step that transforms raw data into a form primed for deeper
# scientific inquiry.
#
# To illustrate its significance, we will first generate sample data.
# coords
time = pd.date_range("2000-01-01", periods=50, freq="YS")
lat = np.linspace(-90, 90, 20)
lon = np.linspace(-180, 180, 30)

lat_da = xr.DataArray(lat, dims="lat", coords={"lat": lat})
lon_da = xr.DataArray(lon, dims="lon", coords={"lon": lon})

# deterministic components
base = 15.0
i = xr.DataArray(np.arange(len(time)), dims="time", coords={"time": time})
slope = 0.1
time_trend = slope * i

lat_pattern = -0.5 * abs(lat_da) / 90
lon_pattern = 0.01 * (lon_da / 180)

# raw noise
np.random.seed(42)
noise_raw = xr.DataArray(
    np.random.normal(0, 0.5, (len(time), len(lat), len(lon))),
    dims=["time", "lat", "lon"],
    coords={"time": time, "lat": lat, "lon": lon},
)

# make noise have NO linear component along time (per grid point)
noise_detr = noise_raw.reduce(signal.detrend, dim="time")  # removes (a*i+b)

# build data with a perfectly linear trend + detrended noise
data_trend = base + lat_pattern + lon_pattern + time_trend + noise_detr

# %%
# :py:func:`easyclimate.calc_detrend_spatial <easyclimate.calc_detrend_spatial>` can help us remove
# temporal trends from spatial data.
detr = ecl.calc_detrend_spatial(data_trend, "time")
detr

# %%
# Next, we consider the correct validation: detr(data_trend) should equal noise_detr.
err = np.max(np.abs(detr - noise_detr))
print("max abs error (should be ~1e-12~1e-10):", float(err))

# %%
# Here we also compared the difference before and after removing the trend at a specific point,
# which can be said to align with our expectations.
data_trend.isel(lat = 10, lon = 20).plot(label = "Trend data")
detr.isel(lat = 10, lon = 20).plot(label = "Detrend data")
plt.legend()

# %%
# To better handle large spatial datasets, easyclimate
# also provides :py:func:`easyclimate.calc_detrend_spatial_fast <easyclimate.calc_detrend_spatial_fast>`
# implementations for significantly faster solutions.
# These methods include: ``"scipy_reduce","scipy","numpy","rust","rust_chunked","rust_flexible"``.
# Their results demonstrate considerable robustness.
detr = ecl.calc_detrend_spatial_fast(data_trend, "time", method = "scipy_reduce")
# ✅ Correct validation: detr(data_trend) should equal noise_detr
err = np.max(np.abs(detr - noise_detr))
print("[scipy_reduce] max abs error (should be ~1e-12~1e-10):", float(err))

detr = ecl.calc_detrend_spatial_fast(data_trend, "time", method = "scipy")
# ✅ Correct validation: detr(data_trend) should equal noise_detr
err = np.max(np.abs(detr - noise_detr))
print("[scipy] max abs error (should be ~1e-12~1e-10):", float(err))

detr = ecl.calc_detrend_spatial_fast(data_trend, "time", method = "rust")
# ✅ Correct validation: detr(data_trend) should equal noise_detr
err = np.max(np.abs(detr - noise_detr))
print("[rust] max abs error (should be ~1e-12~1e-10):", float(err))

detr = ecl.calc_detrend_spatial_fast(data_trend, "time", method = "rust_chunked")
# ✅ Correct validation: detr(data_trend) should equal noise_detr
err = np.max(np.abs(detr - noise_detr))
print("[rust_chunked] max abs error (should be ~1e-12~1e-10):", float(err))

detr = ecl.calc_detrend_spatial_fast(data_trend, "time", method = "rust_flexible")
# ✅ Correct validation: detr(data_trend) should equal noise_detr
err = np.max(np.abs(detr - noise_detr))
print("[rust_flexible] max abs error (should be ~1e-12~1e-10):", float(err))
