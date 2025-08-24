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
    width=0.8,
    threshold: float = 0,
    pos_color: str = "red",
    neg_color: str = "blue",
    ax=None,
    **kwargs,
) -> matplotlib.container.BarContainer:
    """
    Plot a bar chart with time for a 1D :py:class:`xarray.DataArray <xarray.DataArray>` with bars colored based on a threshold value.

    Parameters:
    -----------
    da : :py:class:`xarray.DataArray <xarray.DataArray>`
        1-dimensional data array to plot
    width: :py:class:`float <float>` or array-like, default: 0.8
        The width(s) of the bars.

        .. note::

            If x has units (e.g., datetime), then the width is converted to a multiple of the width relative to the difference units of the x values (e.g., time difference).

    threshold : :py:class:`float <float>`, optional
        Threshold value for color separation (default: 0)
    pos_color : :py:class:`str <str>`, optional
        Color for bars ≥ threshold (default: 'red')
    neg_color : :py:class:`str <str>`, optional
        Color for bars < threshold (default: 'blue')
    ax : matplotlib axes, optional
        Axes object to plot on (uses current axes if None)
    **kwargs :
        Additional arguments passed to plt.bar

    Returns:
    --------
    matplotlib.container.BarContainer
        The bar plot object

    .. seealso::

        :py:func:`matplotlib.pyplot.bar <matplotlib.pyplot.bar>`

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_ao_index.py
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

    # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.bar.html
    # https://stackoverflow.com/questions/59089739/matplotlib-datetime-x-axis-and-bar-widths
    if x.dtype == "datetime64[ns]":
        width_value = (x[1] - x[0]) * width
        kwargs.update({"width": width_value})
    else:
        kwargs.update({"width": width})

    # Assign colors based on threshold
    colors = np.where(y >= threshold, pos_color, neg_color)

    # Create bar plot
    bars = ax.bar(x, y, color=colors, **kwargs)

    # Add threshold reference line
    ax.axhline(threshold, color="gray", linestyle="--", linewidth=0.8)

    return bars
