"""
pytest for plot/line.py
"""

import pytest

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from easyclimate.plot.line import line_plot_with_threshold


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_line_plot_with_threshold1():
    # Set random seed for reproducible test data
    np.random.seed(11)

    data = xr.DataArray(
        np.sin(np.linspace(0, 4 * np.pi, 100)) + 0.3 * np.random.randn(100),
        dims=["x"],
        coords={"x": np.linspace(0, 10, 100)},
    )

    fig, ax = plt.subplots(figsize=(10, 5))
    line_plot_with_threshold(
        data,
        threshold=0,
        pos_color="red",
        neg_color="blue",
        line_kwargs={"color": "black", "linewidth": 1.5},
        fill_kwargs={"alpha": 0.4},
    )
    ax.set_title("Line Plot with Threshold Shading")
    plt.grid(True, linestyle=":", alpha=0.5)
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_line_plot_with_threshold2():
    # Set random seed for reproducible test data
    np.random.seed(11)

    data = xr.DataArray(
        np.sin(np.linspace(0, 4 * np.pi, 100)) + 0.3 * np.random.randn(100),
        dims=["x"],
        coords={"x": np.linspace(0, 10, 100)},
    )

    fig, ax = plt.subplots(figsize=(10, 5))
    line_plot_with_threshold(
        data,
        threshold=0,
        pos_color="red",
        neg_color="blue",
        line_kwargs={"color": "black", "linewidth": 1.5},
        fill_kwargs={"alpha": 0.4},
        fill_neg_plot=False,
    )
    ax.set_title("Line Plot with Threshold Shading")
    plt.grid(True, linestyle=":", alpha=0.5)
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_line_plot_with_threshold3():
    # Set random seed for reproducible test data
    np.random.seed(11)

    data = xr.DataArray(
        np.sin(np.linspace(0, 4 * np.pi, 100)) + 0.3 * np.random.randn(100),
        dims=["x"],
        coords={"x": np.linspace(0, 10, 100)},
    )

    fig, ax = plt.subplots(figsize=(10, 5))
    line_plot_with_threshold(
        data,
        threshold=0,
        pos_color="red",
        neg_color="blue",
        line_kwargs={"color": "black", "linewidth": 1.5},
        fill_kwargs={"alpha": 0.4},
        fill_pos_plot=False,
    )
    ax.set_title("Line Plot with Threshold Shading")
    plt.grid(True, linestyle=":", alpha=0.5)
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_line_plot_with_threshold4():
    # Set random seed for reproducible test data
    np.random.seed(11)

    data = xr.DataArray(
        np.sin(np.linspace(0, 4 * np.pi, 100)) + 0.3 * np.random.randn(100),
        dims=["x"],
        coords={"x": np.linspace(0, 10, 100)},
    )

    fig, ax = plt.subplots(figsize=(10, 5))
    line_plot_with_threshold(
        data,
        threshold=0,
        pos_color="red",
        neg_color="blue",
        line_kwargs={"color": "black", "linewidth": 1.5},
        fill_kwargs={"alpha": 0.4},
        line_plot=False,
    )
    ax.set_title("Line Plot with Threshold Shading")
    plt.grid(True, linestyle=":", alpha=0.5)
    return fig


def test_line_plot_with_threshold_input_validation():
    """Test input validation for line_plot_with_threshold function."""

    # Create test cases

    # 1. Test with valid 1D input (should not raise an exception)
    valid_1d = xr.DataArray(
        np.random.rand(10), dims=["time"], coords={"time": np.arange(10)}
    )
    try:
        line_plot_with_threshold(valid_1d)
    except ValueError:
        pytest.fail("Unexpected ValueError for valid 1D input")

    # 2. Test with 2D input (should raise ValueError)
    invalid_2d = xr.DataArray(
        np.random.rand(10, 5),
        dims=["time", "space"],
        coords={"time": np.arange(10), "space": np.arange(5)},
    )
    with pytest.raises(ValueError) as excinfo:
        line_plot_with_threshold(invalid_2d)
    assert "Input must be 1-dimensional" in str(excinfo.value)

    # 3. Test with 0D input (should raise ValueError)
    invalid_0d = xr.DataArray(42)
    with pytest.raises(ValueError) as excinfo:
        line_plot_with_threshold(invalid_0d)
    assert "Input must be 1-dimensional" in str(excinfo.value)

    # 4. Test with 3D input (should raise ValueError)
    invalid_3d = xr.DataArray(
        np.random.rand(10, 5, 3),
        dims=["time", "space", "level"],
        coords={"time": np.arange(10), "space": np.arange(5), "level": np.arange(3)},
    )
    with pytest.raises(ValueError) as excinfo:
        line_plot_with_threshold(invalid_3d)
    assert "Input must be 1-dimensional" in str(excinfo.value)
