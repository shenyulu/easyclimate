"""
pytest for plot/bar.py
"""

import pytest

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from easyclimate.plot.bar import bar_plot_with_threshold


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_bar_plot_with_threshold1():
    np.random.seed(10)

    # Create sample data
    data = xr.DataArray(
        np.random.randn(10) * 2, dims=["x"], coords={"x": np.arange(10)}
    )

    fig, ax = plt.subplots()
    bar_plot_with_threshold(data, threshold=0)
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_bar_plot_with_threshold2():
    np.random.seed(10)

    # Create sample data
    data = xr.DataArray(
        np.random.randn(10) * 2, dims=["x"], coords={"x": np.arange(10)}
    )

    fig, ax = plt.subplots()
    bar_plot_with_threshold(data, threshold=0, ax=ax)
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_bar_plot_with_threshold3():
    np.random.seed(10)

    # Create sample data
    data = xr.DataArray(
        np.array([-0.0133197, -0.09011592, 0.03289608, -0.17159361, -0.50632054]),
        dims="time",
        coords={"time": pd.date_range("1981-12-01", freq="ME", periods=5)},
    )

    fig, ax = plt.subplots()
    bar_plot_with_threshold(data, threshold=0)
    return fig


def test_line_plot_with_threshold_input_validation():
    """Test input validation for line_plot_with_threshold function."""
    np.random.seed(10)

    # Create test cases

    # 1. Test with valid 1D input (should not raise an exception)
    valid_1d = xr.DataArray(
        np.random.rand(10), dims=["time"], coords={"time": np.arange(10)}
    )
    try:
        bar_plot_with_threshold(valid_1d)
    except ValueError:
        pytest.fail("Unexpected ValueError for valid 1D input")

    # 2. Test with 2D input (should raise ValueError)
    invalid_2d = xr.DataArray(
        np.random.rand(10, 5),
        dims=["time", "space"],
        coords={"time": np.arange(10), "space": np.arange(5)},
    )
    with pytest.raises(ValueError) as excinfo:
        bar_plot_with_threshold(invalid_2d)
    assert "Input DataArray must be 1-dimensional" in str(excinfo.value)

    # 3. Test with 0D input (should raise ValueError)
    invalid_0d = xr.DataArray(42)
    with pytest.raises(ValueError) as excinfo:
        bar_plot_with_threshold(invalid_0d)
    assert "Input DataArray must be 1-dimensional" in str(excinfo.value)

    # 4. Test with 3D input (should raise ValueError)
    invalid_3d = xr.DataArray(
        np.random.rand(10, 5, 3),
        dims=["time", "space", "level"],
        coords={"time": np.arange(10), "space": np.arange(5), "level": np.arange(3)},
    )
    with pytest.raises(ValueError) as excinfo:
        bar_plot_with_threshold(invalid_3d)
    assert "Input DataArray must be 1-dimensional" in str(excinfo.value)
