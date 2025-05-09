"""
Bar plot for xarray dataset
"""

import matplotlib.container
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import matplotlib

__all__ = ["bar_plot_with_threshold"]


def bar_plot_with_threshold(
    da: xr.DataArray,
    threshold: float = 0,
    pos_color: str = "red",
    neg_color: str = "blue",
    ax=None,
    **kwargs,
) -> matplotlib.container.BarContainer:
    """
    Plot a bar chart for a 1D xarray.DataArray with bars colored based on a threshold value.

    Parameters:
    -----------
    da : xarray.DataArray
        1-dimensional data array to plot
    threshold : float, optional
        Threshold value for color separation (default: 0)
    pos_color : str, optional
        Color for bars ≥ threshold (default: 'red')
    neg_color : str, optional
        Color for bars < threshold (default: 'blue')
    ax : matplotlib axes, optional
        Axes object to plot on (uses current axes if None)
    **kwargs :
        Additional arguments passed to plt.bar

    Returns:
    --------
    matplotlib.container.BarContainer
        The bar plot object

    Raises:
    -------
    ValueError
        If input DataArray is not 1-dimensional
    """
    # Verify 1D data
    if len(da.dims) != 1:
        raise ValueError("Input DataArray must be 1-dimensional")

    # Use current axes if none provided
    if ax is None:
        ax = plt.gca()

    # Extract coordinates and values
    x = da.coords[da.dims[0]].values
    y = da.values

    # Assign colors based on threshold
    colors = np.where(y >= threshold, pos_color, neg_color)

    # Create bar plot
    bars = ax.bar(x, y, color=colors, **kwargs)

    # Add threshold reference line
    ax.axhline(threshold, color="gray", linestyle="--", linewidth=0.8)

    return bars
