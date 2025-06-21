"""
pytest for filter/gaussian_filter.py
"""

import pytest
import numpy as np
import xarray as xr
from scipy.ndimage import gaussian_filter1d
from easyclimate.filter import calc_gaussian_filter


@pytest.fixture
def sample_dataarray():
    """Create a sample xarray DataArray for testing."""
    time = np.arange(100)
    data = np.sin(time * 0.1) + np.random.normal(0, 0.1, size=100)
    da = xr.DataArray(
        data, dims=["time"], coords={"time": time}, attrs={"description": "test data"}
    )
    return da


def test_gaussian_filter_basic(sample_dataarray):
    """Test basic functionality of Gaussian filter."""
    window_length = 5
    result = calc_gaussian_filter(sample_dataarray, window_length)

    assert isinstance(result, xr.DataArray)
    assert result.shape == sample_dataarray.shape
    assert result.dims == sample_dataarray.dims
    assert np.all(result.values != sample_dataarray.values)  # Should be smoothed


def test_gaussian_filter_sigma_calculation(sample_dataarray):
    """Test automatic sigma calculation."""
    window_length = 5
    expected_sigma = window_length / np.sqrt(8 * np.log(2))

    result = calc_gaussian_filter(sample_dataarray, window_length)
    manual_result = calc_gaussian_filter(
        sample_dataarray, window_length, sigma=expected_sigma
    )

    assert np.allclose(result.values, manual_result.values)


def test_gaussian_filter_custom_sigma(sample_dataarray):
    """Test with custom sigma value."""
    window_length = 5
    custom_sigma = 2.0
    result = calc_gaussian_filter(sample_dataarray, window_length, sigma=custom_sigma)

    manual_smoothed = gaussian_filter1d(
        sample_dataarray.values,
        sigma=custom_sigma,
        axis=sample_dataarray.get_axis_num("time"),
    )

    assert np.allclose(result.values, manual_smoothed)


def test_gaussian_filter_different_dim():
    """Test filtering along a different dimension."""
    data = np.random.normal(0, 1, size=(10, 20))
    da = xr.DataArray(
        data, dims=["x", "y"], coords={"x": np.arange(10), "y": np.arange(20)}
    )

    result = calc_gaussian_filter(da, window_length=3, dim="y")
    manual_smoothed = gaussian_filter1d(
        da.values, sigma=3 / np.sqrt(8 * np.log(2)), axis=1
    )

    assert np.allclose(result.values, manual_smoothed)


def test_gaussian_filter_keep_attrs(sample_dataarray):
    """Test attribute preservation."""
    result = calc_gaussian_filter(sample_dataarray, window_length=5, keep_attrs=True)
    assert result.attrs == sample_dataarray.attrs

    result_no_attrs = calc_gaussian_filter(sample_dataarray, window_length=5)
    assert result_no_attrs.attrs == {}


def test_gaussian_filter_edge_cases():
    """Test edge cases like empty array or single value."""
    # Test empty array
    empty_da = xr.DataArray([], dims=["time"])
    result = calc_gaussian_filter(empty_da, window_length=1)
    assert len(result) == 0

    # Test single value
    single_da = xr.DataArray([1.0], dims=["time"])
    result = calc_gaussian_filter(single_da, window_length=1)
    assert np.isclose(result.values[0], 1.0)


def test_gaussian_filter_invalid_input(sample_dataarray):
    """Test invalid input handling."""

    with pytest.raises(ValueError):
        calc_gaussian_filter(sample_dataarray, window_length=-1)

    with pytest.raises(ValueError):
        calc_gaussian_filter(sample_dataarray, window_length=5, dim="invalid_dim")
