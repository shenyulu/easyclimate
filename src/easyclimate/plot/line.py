"""
Bar plot for xarray dataset
"""

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

__all__ = ["line_plot_with_threshold"]


def line_plot_with_threshold(
    da: xr.DataArray,
    threshold: float = 0,
    pos_color: str = "red",
    neg_color: str = "blue",
    ax=None,
    line_plot: bool = True,
    fill_pos_plot: bool = True,
    fill_neg_plot: bool = True,
    line_kwargs=None,
    fill_kwargs=None,
) -> tuple:
    """
    Plot a line chart with proper shading at threshold crossings.

    Parameters:
    -----------
    da : :py:class:`xarray.DataArray <xarray.DataArray>`
        1-dimensional data array
    threshold : :py:class:`float <float>`, optional
        Color separation threshold (default: 0)
    pos_color : :py:class:`str <str>`, optional
        Color for values ≥ threshold (default: 'red')
    neg_color : :py:class:`str <str>`, optional
        Color for values < threshold (default: 'blue')
    ax : matplotlib axes, optional
        Axes to plot on (default: current axes)
    line_kwargs : :py:class:`dict <dict>`, optional
        Arguments for plt.plot
    fill_kwargs : :py:class:`dict <dict>`, optional
        Arguments for plt.fill_between

    Returns:
    --------
    tuple
        (line plot, fill objects)

    .. seealso::

        :py:func:`matplotlib.lines.Line2D <matplotlib.lines.Line2D>`

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_basic_statistical_analysis.py
        ./dynamic_docs/plot_corr_reg.py
    """
    # Input validation
    if len(da.dims) != 1:
        raise ValueError("Input must be 1-dimensional")

    ax = ax or plt.gca()
    line_kwargs = line_kwargs or {"color": "black", "linewidth": 1.5}
    fill_kwargs = fill_kwargs or {"alpha": 0.3}

    x = da.coords[da.dims[0]].values
    y = da.values

    # Find exact crossing points
    x_fine, y_fine = _interpolate_crossings(x, y, threshold)

    # Plot line and shaded regions
    if line_plot == True:
        line = ax.plot(x, y, **line_kwargs)
    else:
        line = None

    # Fill above threshold
    if fill_pos_plot == True:
        fill_pos = ax.fill_between(
            x_fine,
            y_fine,
            threshold,
            where=(y_fine >= threshold),
            facecolor=pos_color,
            **fill_kwargs,
        )
    else:
        fill_pos = None

    # Fill below threshold
    if fill_neg_plot == True:
        fill_neg = ax.fill_between(
            x_fine,
            y_fine,
            threshold,
            where=(y_fine < threshold),
            facecolor=neg_color,
            **fill_kwargs,
        )
    else:
        fill_neg = None

    ax.axhline(threshold, color="gray", linestyle="--", linewidth=0.8)
    return (line, (fill_pos, fill_neg))


def _interpolate_crossings(x, y, threshold):
    """Add interpolated points at threshold crossings"""
    x_fine = []
    y_fine = []

    for i in range(len(x) - 1):
        x_fine.append(x[i])
        y_fine.append(y[i])

        # Check if this segment crosses threshold
        if (y[i] < threshold and y[i + 1] > threshold) or (
            y[i] > threshold and y[i + 1] < threshold
        ):

            # Linear interpolation to find exact crossing
            t = (threshold - y[i]) / (y[i + 1] - y[i])
            x_cross = x[i] + t * (x[i + 1] - x[i])

            x_fine.append(x_cross)
            y_fine.append(threshold)

    # Add last point
    x_fine.append(x[-1])
    y_fine.append(y[-1])

    return np.array(x_fine), np.array(y_fine)
