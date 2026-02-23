"""
pytest for stat.py
"""

import pytest

import numpy as np
import pandas as pd
import xarray as xr
from scipy import signal
import easyclimate as ecl


@pytest.fixture
def sample_data():
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
    return [data_trend, noise_detr]


def calc_detrend_spatial_fast1(sample_data):
    data_trend = sample_data[0]
    noise_detr = sample_data[1]

    detr = ecl.calc_detrend_spatial_fast(data_trend, "time")
    err = np.max(np.abs(detr - noise_detr))

    assert err < 1e-6


def calc_detrend_spatial_fast2(sample_data):
    data_trend = sample_data[0]
    noise_detr = sample_data[1]

    detr = ecl.calc_detrend_spatial_fast(data_trend, "time", method="scipy_reduce")
    err = np.max(np.abs(detr - noise_detr))

    assert err < 1e-6


def calc_detrend_spatial_fast3(sample_data):
    data_trend = sample_data[0]
    noise_detr = sample_data[1]

    detr = ecl.calc_detrend_spatial_fast(data_trend, "time", method="scipy")
    err = np.max(np.abs(detr - noise_detr))

    assert err < 1e-6


def calc_detrend_spatial_fast4(sample_data):
    data_trend = sample_data[0]
    noise_detr = sample_data[1]

    detr = ecl.calc_detrend_spatial_fast(data_trend, "time", method="rust")
    err = np.max(np.abs(detr - noise_detr))

    assert err < 1e-6


def calc_detrend_spatial_fast5(sample_data):
    data_trend = sample_data[0]
    noise_detr = sample_data[1]

    detr = ecl.calc_detrend_spatial_fast(data_trend, "time", method="rust_chunked")
    err = np.max(np.abs(detr - noise_detr))

    assert err < 1e-6


def calc_detrend_spatial_fast6(sample_data):
    data_trend = sample_data[0]
    noise_detr = sample_data[1]

    detr = ecl.calc_detrend_spatial_fast(data_trend, "time", method="rust_flexible")
    err = np.max(np.abs(detr - noise_detr))

    assert err < 1e-6
