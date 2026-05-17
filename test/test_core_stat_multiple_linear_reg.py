"""
pytest for stat.py for ``calc_multiple_linear_regression_spatial``
"""

import pytest
import numpy as np
import xarray as xr
from typing import List
import warnings
from easyclimate.core.stat import calc_multiple_linear_regression_spatial


def test_basic_multiple_regression():
    """Test basic multiple linear regression with synthetic data."""
    # Create test data
    time = np.arange(10)
    lat = np.arange(3)
    lon = np.arange(4)

    # Create dependent variable (y = 2*x1 + 3*x2 + 5 + noise)
    x1 = xr.DataArray(
        np.random.rand(10, 3, 4),
        dims=("time", "lat", "lon"),
        coords={"time": time, "lat": lat, "lon": lon},
    )
    x2 = xr.DataArray(
        np.random.rand(10, 3, 4),
        dims=("time", "lat", "lon"),
        coords={"time": time, "lat": lat, "lon": lon},
    )
    noise = xr.DataArray(
        np.random.normal(0, 0.1, (10, 3, 4)),
        dims=("time", "lat", "lon"),
        coords={"time": time, "lat": lat, "lon": lon},
    )
    y = 2 * x1 + 3 * x2 + 5 + noise

    # Run regression
    result = calc_multiple_linear_regression_spatial(y, [x1, x2])

    # Check results
    assert "slopes" in result
    assert "intercept" in result
    assert "r_squared" in result
    assert "slopes_p" in result
    assert "intercept_p" in result

    # Check shapes
    assert result.slopes.shape == (3, 4, 2)  # 2 coefficients, 3 lats, 4 lons
    assert result.intercept.shape == (3, 4)
    assert result.r_squared.shape == (3, 4)
    assert result.slopes_p.shape == (3, 4, 2)
    assert result.intercept_p.shape == (3, 4)

    # Check approximate values (should be close to our coefficients)
    np.testing.assert_allclose(result.slopes.mean(dim=("lat", "lon")), [2, 3], rtol=0.2)
    np.testing.assert_allclose(result.intercept.mean(), 5, rtol=0.2)
    assert (result.r_squared > 0.9).all()  # Good fit expected


def test_single_predictor():
    """Test with just one predictor (simple linear regression)."""
    time = np.arange(20)
    lat = np.arange(2)
    lon = np.arange(2)

    x = xr.DataArray(
        np.random.rand(20, 2, 2),
        dims=("time", "lat", "lon"),
        coords={"time": time, "lat": lat, "lon": lon},
    )
    y = 1.5 * x + 2.0 + np.random.normal(0, 0.1, (20, 2, 2))

    result = calc_multiple_linear_regression_spatial(y, [x])

    assert result.slopes.shape == (2, 2, 1)
    np.testing.assert_allclose(result.slopes.mean(), 1.5, rtol=0.2)
    np.testing.assert_allclose(result.intercept.mean(), 2.0, rtol=0.2)


def test_time_dimension_name():
    """Test with custom time dimension name."""
    t = np.arange(10)
    lat = np.arange(2)

    x1 = xr.DataArray(
        np.random.rand(10, 2),
        dims=("custom_time", "lat"),
        coords={"custom_time": t, "lat": lat},
    )
    y = 2 * x1 + 3 + np.random.normal(0, 0.1, (10, 2))

    result = calc_multiple_linear_regression_spatial(y, [x1], dim="custom_time")

    assert result.slopes.shape == (2, 1)
    np.testing.assert_allclose(result.slopes.mean(), 2, rtol=0.2)


def test_time_coordinate_mismatch():
    """Test that warning is raised when time coordinates don't match."""
    t1 = np.arange(10)
    t2 = np.arange(5, 15)
    lat = np.arange(3)

    x = xr.DataArray(
        np.random.rand(10, 3), dims=("time", "lat"), coords={"time": t1, "lat": lat}
    )
    y = xr.DataArray(
        np.random.rand(10, 3), dims=("time", "lat"), coords={"time": t2, "lat": lat}
    )

    with pytest.warns(UserWarning):
        calc_multiple_linear_regression_spatial(y, [x])


def test_empty_input():
    """Test with empty input list of predictors."""
    time = np.arange(5)
    y = xr.DataArray(np.random.rand(5), dims=("time",), coords={"time": time})

    with pytest.raises(ValueError):
        calc_multiple_linear_regression_spatial(y, [])


def test_r_squared_range():
    """Test that R-squared values are between 0 and 1."""
    time = np.arange(10)
    lat = np.arange(2)

    x = xr.DataArray(
        np.random.rand(10, 2), dims=("time", "lat"), coords={"time": time, "lat": lat}
    )
    y = xr.DataArray(
        np.random.rand(10, 2), dims=("time", "lat"), coords={"time": time, "lat": lat}
    )

    result = calc_multiple_linear_regression_spatial(y, [x])

    assert (result.r_squared >= 0).all() and (result.r_squared <= 1).all()


def test_p_values_range():
    """Test that p-values are between 0 and 1."""
    time = np.arange(10)
    lat = np.arange(2)

    x = xr.DataArray(
        np.random.rand(10, 2), dims=("time", "lat"), coords={"time": time, "lat": lat}
    )
    y = 2 * x + 1 + np.random.normal(0, 0.1, (10, 2))

    result = calc_multiple_linear_regression_spatial(y, [x])

    assert (result.slopes_p >= 0).all() and (result.slopes_p <= 1).all()
    assert (result.intercept_p >= 0).all() and (result.intercept_p <= 1).all()


def test_missing_values_are_skipped_per_grid_point():
    """Test that NaNs in y or predictors do not make lstsq fail."""
    time = np.arange(8)
    lat = np.arange(2)

    x1_values = np.tile(time[:, None], (1, 2)).astype(float)
    x2_values = np.tile((time[:, None] ** 2), (1, 2)).astype(float)
    y_values = 2 * x1_values - 0.5 * x2_values + 4

    y_values[0, 0] = np.nan
    x1_values[1, 0] = np.nan
    x2_values[:, 1] = np.nan

    x1 = xr.DataArray(
        x1_values, dims=("time", "lat"), coords={"time": time, "lat": lat}
    )
    x2 = xr.DataArray(
        x2_values, dims=("time", "lat"), coords={"time": time, "lat": lat}
    )
    y = xr.DataArray(y_values, dims=("time", "lat"), coords={"time": time, "lat": lat})

    result = calc_multiple_linear_regression_spatial(y, [x1, x2])

    np.testing.assert_allclose(result.slopes.sel(lat=0), [2, -0.5], atol=1e-12)
    np.testing.assert_allclose(result.intercept.sel(lat=0), 4, atol=1e-12)
    assert np.isnan(result.slopes.sel(lat=1)).all()
    assert np.isnan(result.intercept.sel(lat=1))
