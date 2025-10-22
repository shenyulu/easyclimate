"""
pytest for core/normalized.py
"""

import pytest
import xarray as xr
import numpy as np
from easyclimate.core.normalized import (
    normalize_zscore,
    normalize_minmax,
    normalize_robust,
    normalize_mean,
)


# Shared fixtures
@pytest.fixture
def sample_da():
    time = np.array(["2020-01-01", "2020-01-02", "2020-01-03"], dtype="datetime64[D]")
    data = np.array(
        [
            [[1.0, 2.0], [3.0, 4.0]],
            [[5.0, 6.0], [7.0, 8.0]],
            [[9.0, 10.0], [11.0, 12.0]],
        ]
    )
    return xr.DataArray(
        data,
        dims=["time", "lat", "lon"],
        coords={"time": time, "lat": [0, 1], "lon": [0, 1]},
    )


@pytest.fixture
def sample_da_with_nan():
    time = np.array(["2020-01-01", "2020-01-02", "2020-01-03"], dtype="datetime64[D]")
    data = np.array(
        [
            [[1.0, 2.0], [3.0, np.nan]],
            [[5.0, 6.0], [7.0, 8.0]],
            [[9.0, 10.0], [11.0, 12.0]],
        ]
    )
    return xr.DataArray(
        data,
        dims=["time", "lat", "lon"],
        coords={"time": time, "lat": [0, 1], "lon": [0, 1]},
    )


@pytest.fixture
def sample_da_constant():
    time = np.array(["2020-01-01", "2020-01-02", "2020-01-03"], dtype="datetime64[D]")
    data = np.ones((3, 2, 2)) * 5.0
    return xr.DataArray(
        data,
        dims=["time", "lat", "lon"],
        coords={"time": time, "lat": [0, 1], "lon": [0, 1]},
    )


# Test class for normalize_zscore
class TestNormalizeZscore:
    def test_correctness(self, sample_da):
        result = normalize_zscore(sample_da, dim="time", ddof=1)
        # Check output type
        assert isinstance(result, xr.DataArray)
        # Check mean and std along time dimension
        assert np.allclose(result.mean(dim="time"), 0.0, atol=1e-6)
        assert np.allclose(result.std(dim="time", ddof=1), 1.0, atol=1e-6)
        # Check coordinates preservation
        assert result.dims == sample_da.dims
        assert result.coords["time"].equals(sample_da.coords["time"])

    # def test_invalid_dim(self, sample_da):
    #     with pytest.raises(ValueError, match=".*invalid_dim.*"):
    #         normalize_zscore(sample_da, dim="invalid_dim")

    def test_constant_data(self, sample_da_constant):
        result = normalize_zscore(sample_da_constant, dim="time", ddof=1)
        # Constant data should result in NaN due to zero standard deviation
        assert np.all(np.isnan(result))


# Test class for normalize_minmax
class TestNormalizeMinmax:
    def test_correctness(self, sample_da):
        result = normalize_minmax(sample_da, dim="time", feature_range=(0, 1))
        # Check output type
        assert isinstance(result, xr.DataArray)
        # Check range
        assert np.allclose(result.min(dim="time"), 0.0, atol=1e-6)
        assert np.allclose(result.max(dim="time"), 1.0, atol=1e-6)
        # Check coordinates preservation
        assert result.dims == sample_da.dims
        assert result.coords["time"].equals(sample_da.coords["time"])

    def test_custom_range(self, sample_da):
        result = normalize_minmax(sample_da, dim="time", feature_range=(-1, 1))
        assert np.allclose(result.min(dim="time"), -1.0, atol=1e-6)
        assert np.allclose(result.max(dim="time"), 1.0, atol=1e-6)

    # def test_invalid_dim(self, sample_da):
    #     with pytest.raises(ValueError, match=".*invalid_dim.*"):
    #         normalize_minmax(sample_da, dim="invalid_dim")

    def test_constant_data(self, sample_da_constant):
        result = normalize_minmax(sample_da_constant, dim="time")
        # Constant data should result in NaN due to zero range
        assert np.all(np.isnan(result))


# Test class for normalize_robust
class TestNormalizeRobust:
    def test_correctness(self, sample_da):
        result = normalize_robust(sample_da, dim="time", q_low=0.25, q_high=0.75)
        # Check output type
        assert isinstance(result, xr.DataArray)
        # Check median is approximately zero
        assert np.allclose(result.median(dim="time"), 0.0, atol=1e-6)
        # Check coordinates preservation
        assert result.dims == sample_da.dims
        assert result.coords["time"].equals(sample_da.coords["time"])

    # def test_invalid_dim(self, sample_da):
    #     with pytest.raises(ValueError, match=".*invalid_dim.*"):
    #         normalize_robust(sample_da, dim="invalid_dim")

    def test_constant_data(self, sample_da_constant):
        result = normalize_robust(sample_da_constant, dim="time")
        # Constant data should result in NaN due to zero IQR
        assert np.all(np.isnan(result))


# Test class for normalize_mean
class TestNormalizeMean:
    def test_correctness(self, sample_da):
        result = normalize_mean(sample_da, dim="time")
        # Check output type
        assert isinstance(result, xr.DataArray)
        # Check mean is approximately zero
        assert np.allclose(result.mean(dim="time"), 0.0, atol=1e-6)
        # Check coordinates preservation
        assert result.dims == sample_da.dims
        assert result.coords["time"].equals(sample_da.coords["time"])

    # def test_invalid_dim(self, sample_da):
    #     with pytest.raises(ValueError, match=".*invalid_dim.*"):
    #         normalize_mean(sample_da, dim="invalid_dim")

    def test_constant_data(self, sample_da_constant):
        result = normalize_mean(sample_da_constant, dim="time")
        # Constant data should result in NaN due to zero range
        assert np.all(np.isnan(result))
