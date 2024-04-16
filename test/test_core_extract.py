"""
pytest for extract.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import pandas as pd
import xarray as xr

# data1
# For the `freq` alteration from `Y` to `YE`, please refer as follows:
# Reference: https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#timeseries-offset-aliases
# https://github.com/pandas-dev/pandas/issues/9586
timearray = pd.date_range("2020-01-01", freq="YE", periods=5)
data1 = xr.DataArray(np.arange(0, 5, 1), dims=("time"), coords={"time": timearray})

# data2
timearray = pd.date_range("2020-01-01", freq="ME", periods=5)
data2 = xr.DataArray(np.arange(0, 5, 1), dims=("time"), coords={"time": timearray})

# data3
timearray = pd.date_range("2020-01-01", freq="D", periods=5)
data3 = xr.DataArray(np.arange(0, 5, 1), dims=("time"), coords={"time": timearray})

# data4
timearray = pd.date_range("2020-01-01", freq="h", periods=5)
data4 = xr.DataArray(np.arange(0, 5, 1), dims=("time"), coords={"time": timearray})

# data5
timearray = pd.date_range("2020-01-01", freq="min", periods=5)
data5 = xr.DataArray(np.arange(0, 5, 1), dims=("time"), coords={"time": timearray})

# data6
timearray = pd.date_range("2020-01-01", freq="s", periods=5)
data6 = xr.DataArray(np.arange(0, 5, 1), dims=("time"), coords={"time": timearray})

# data7
timearray = pd.date_range("2020-01-01", freq="ms", periods=5)
data7 = xr.DataArray(np.arange(0, 5, 1), dims=("time"), coords={"time": timearray})

# data8
timearray = pd.date_range("2020-01-01", freq="ns", periods=5)
data8 = xr.DataArray(np.arange(0, 5, 1), dims=("time"), coords={"time": timearray})

# data9
timearray = pd.date_range("2020-01-01", freq="D", periods=14)
data9 = xr.DataArray(np.arange(0, 14, 1), dims=("time"), coords={"time": timearray})

# data10
timearray = pd.date_range("2020-01-01", freq="ME", periods=40)
data10 = xr.DataArray(np.arange(0, 40, 1), dims=("time"), coords={"time": timearray})

# data_sst_nino34
data_sst_nino34 = xr.DataArray(
    np.array(
        [
            1.0162432,
            0.5076949,
            -0.67902166,
            -0.70760345,
            0.11068994,
            1.3218967,
            -1.1079785,
            -0.8168568,
            0.14645328,
            0.63610655,
            0.64269495,
            0.3659633,
            0.3993566,
            -0.20599838,
            -0.4854034,
            1.2408715,
            -0.15124501,
            -1.1885,
            -0.89531666,
            -0.25720003,
            0.7296249,
            0.2966999,
            0.37730327,
            0.11304826,
            0.16522326,
            -0.5311017,
            -0.73997337,
            0.3690299,
            -0.45829165,
            -0.83581,
            -0.02034007,
            -0.22075672,
            0.2399766,
            1.5814567,
            0.37708327,
            -0.09190174,
            0.08657158,
            0.5442299,
            -0.35217008,
            -0.5423617,
            -0.9053801,
        ]
    ),
    dims="time",
    coords={"time": pd.date_range("1982-01-01", "2022-12-31", freq="YE")},
)


def test_get_specific_years_data():
    result_data = ecl.get_specific_years_data(data1, [2021, 2024]).data
    refer_data = np.array([1, 4])
    assert np.isclose(result_data, refer_data).all()


def test_get_specific_months_data():
    result_data = ecl.get_specific_months_data(data2, [2, 4]).data
    refer_data = np.array([1, 3])
    assert np.isclose(result_data, refer_data).all()


def test_get_specific_days_data():
    result_data = ecl.get_specific_days_data(data3, [1, 4]).data
    refer_data = np.array([0, 3])
    assert np.isclose(result_data, refer_data).all()


def test_get_specific_hours_data():
    result_data = ecl.get_specific_hours_data(data4, [2, 4]).data
    refer_data = np.array([2, 4])
    assert np.isclose(result_data, refer_data).all()


def test_get_specific_minutes_data():
    result_data = ecl.get_specific_minutes_data(data5, [0, 3]).data
    refer_data = np.array([0, 3])
    assert np.isclose(result_data, refer_data).all()


def test_get_specific_seconds_data():
    result_data = ecl.get_specific_seconds_data(data6, [1, 4]).data
    refer_data = np.array([1, 4])
    assert np.isclose(result_data, refer_data).all()


def test_get_specific_microseconds_data():
    result_data = ecl.get_specific_microseconds_data(data7, [1000, 2000]).data
    refer_data = np.array([1, 2])
    assert np.isclose(result_data, refer_data).all()


def test_get_specific_nanoseconds_data():
    result_data = ecl.get_specific_nanoseconds_data(data8, [0, 1, 4]).data
    refer_data = np.array([0, 1, 4])
    assert np.isclose(result_data, refer_data).all()


def test_get_specific_dayofweek_data():
    result_data = ecl.get_specific_dayofweek_data(data9, [3, 5]).data
    refer_data = np.array([1, 3, 8, 10])
    assert np.isclose(result_data, refer_data).all()


def test_get_yearmean_for_specific_months_data():
    result_data = ecl.get_yearmean_for_specific_months_data(data10, [10, 11]).data
    refer_data = np.array([9.5, 21.5, 33.5])
    assert np.isclose(result_data, refer_data).all()


def test_get_year_exceed_index_upper_bound():
    result_data = ecl.get_year_exceed_index_upper_bound(data_sst_nino34, 0.9)
    refer_data = np.array([1982, 1987, 1997, 2015])
    assert np.isclose(result_data, refer_data).all()


def test_get_year_exceed_index_lower_bound():
    result_data = ecl.get_year_exceed_index_lower_bound(data_sst_nino34, -0.9)
    refer_data = np.array([1988, 1999, 2022])
    assert np.isclose(result_data, refer_data).all()
