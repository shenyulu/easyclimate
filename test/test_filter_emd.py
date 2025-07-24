"""
pytest for filter/emd.py
"""

import pytest
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
from easyclimate.filter.emd import filter_emd, filter_eemd


# Fixtures for common test data
@pytest.fixture
def sample_time_series():
    """Create a sample time series xarray DataArray for testing"""
    time = np.array([datetime(2020, 1, 1) + timedelta(days=i) for i in range(100)])
    data = np.sin(np.linspace(0, 10 * np.pi, 100)) + np.random.normal(0, 0.1, 100)
    return xr.DataArray(data, dims=["time"], coords={"time": time})


@pytest.fixture
def sample_time_array(sample_time_series):
    """Create a numeric time array for testing"""
    return np.arange(len(sample_time_series))


# Test filter_emd function
def test_filter_emd_basic(sample_time_series):
    """Test basic EMD functionality with default parameters"""
    result = filter_emd(sample_time_series, time_step="D")

    # Check return type
    assert isinstance(result, xr.Dataset)

    # Check input is preserved
    assert "input" in result
    assert np.allclose(result["input"], sample_time_series)

    # Check IMFs exist
    assert any("imf" in var for var in result.data_vars)

    # Check IMFs have same length as input
    for var in result.data_vars:
        if "imf" in var:
            assert len(result[var]) == len(sample_time_series)


def test_filter_emd_with_time_array(sample_time_series, sample_time_array):
    """Test EMD with custom time array"""
    result = filter_emd(sample_time_series, time_step="D", time_array=sample_time_array)

    # Check IMFs were created
    assert any("imf" in var for var in result.data_vars)


def test_filter_emd_parameters(sample_time_series):
    """Test EMD with different parameter combinations"""
    # Test with different spline kinds
    for spline in ["cubic", "linear", "akima"]:
        result = filter_emd(sample_time_series, time_step="D", spline_kind=spline)
        assert any("imf" in var for var in result.data_vars)

    # Test with max_imf limit
    result = filter_emd(sample_time_series, time_step="D", max_imf=3)
    assert len([var for var in result.data_vars if "imf" in var]) >= 3


def test_filter_emd_edge_cases(sample_time_series):
    """Test EMD edge cases"""
    # Test with very short time series
    short_data = sample_time_series.isel(time=slice(0, 10))
    result = filter_emd(short_data, time_step="D")
    assert isinstance(result, xr.Dataset)

    # Test with constant input (should still work but may produce fewer IMFs)
    const_data = sample_time_series.copy()
    const_data[:] = 1.0
    result = filter_emd(const_data, time_step="D")
    assert isinstance(result, xr.Dataset)


# Test filter_eemd function
def test_filter_eemd_basic(sample_time_series):
    """Test basic EEMD functionality with default parameters"""
    result = filter_eemd(sample_time_series, time_step="D")

    # Check return type
    assert isinstance(result, xr.Dataset)

    # Check input is preserved
    assert "input" in result
    assert np.allclose(result["input"], sample_time_series)

    # Check eIMFs exist
    assert any("eimf" in var for var in result.data_vars)

    # Check eIMFs have same length as input
    for var in result.data_vars:
        if "eimf" in var:
            assert len(result[var]) == len(sample_time_series)


def test_filter_eemd_with_noise_seed(sample_time_series):
    """Test EEMD reproducibility with noise seed"""
    result1 = filter_eemd(sample_time_series, time_step="D", noise_seed=42)
    result2 = filter_eemd(sample_time_series, time_step="D", noise_seed=42)

    # Results with same seed should be identical
    for var in result1.data_vars:
        if "eimf" in var:
            assert np.allclose(result1[var], result2[var])


def test_filter_eemd_parameters(sample_time_series):
    """Test EEMD with different parameter combinations"""
    # Test with different noise widths
    for noise_width in [0.01, 0.05, 0.1]:
        result = filter_eemd(
            sample_time_series,
            time_step="D",
            noise_width=noise_width,
            trials=10,  # Reduced for faster testing
        )
        assert any("eimf" in var for var in result.data_vars)

    # Test with parallel processing (if system supports it)
    result = filter_eemd(
        sample_time_series,
        time_step="D",
        parallel=True,
        processes=2,
        trials=10,  # Reduced for faster testing
    )
    assert any("eimf" in var for var in result.data_vars)


def test_filter_eemd_edge_cases(sample_time_series):
    """Test EEMD edge cases"""
    # Test with very short time series
    short_data = sample_time_series.isel(time=slice(0, 10))
    result = filter_eemd(short_data, time_step="D", trials=5)
    assert isinstance(result, xr.Dataset)

    # Test with constant input (should still work but may produce fewer eIMFs)
    const_data = sample_time_series.copy()
    const_data[:] = 1.0
    result = filter_eemd(const_data, time_step="D", trials=5)
    assert isinstance(result, xr.Dataset)


def test_filter_eemd_invalid_noise_seed(sample_time_series):
    """Test EEMD with invalid noise seed raises error"""
    with pytest.raises(ValueError):
        filter_eemd(
            sample_time_series,
            time_step="D",
            noise_seed="invalid",  # Not an int or None
        )


# Test both functions together
def test_emd_vs_eemd(sample_time_series):
    """Compare EMD and EEMD results (they should be different)"""
    emd_result = filter_emd(sample_time_series, time_step="D")
    eemd_result = filter_eemd(sample_time_series, time_step="D", trials=10)

    # They should both produce results but not identical
    assert len(emd_result.data_vars) > 1
    assert len(eemd_result.data_vars) > 1

    # The IMFs should be different between methods
    imf0 = emd_result["imf0"].values
    eimf0 = eemd_result["eimf0"].values
    assert not np.allclose(imf0, eimf0)
