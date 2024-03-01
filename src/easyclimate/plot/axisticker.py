"""
Quick processing of special axes
"""

import matplotlib
import cartopy.mpl.ticker as geoticker
import matplotlib.ticker as ticker


def set_lon_format_axis(ax: matplotlib.axes.Axes, axis: str = "x", **kwargs):
    """
    Setting the axes in longitude format.

    Parameters
    ----------
    ax : :py:class:`matplotlib.axes.Axes`
        The axes to which the boundary will be applied.
    axis: {'x', 'y'}, default: 'x'
        The axis to which the parameters are applied.
    **kwargs
        Additional keyword arguments to wrapped :py:func:`matplotlib.axis.Axis.set_major_formatter <matplotlib:matplotlib.axis.Axis.set_major_formatter>`.
    """
    if axis == "x":
        axis = ax.xaxis
    elif axis == "y":
        axis = ax.yaxis
    else:
        raise ValueError("`axis`  should be `x` or `y`.")

    axis.set_major_formatter(geoticker.LongitudeFormatter(**kwargs))


def set_lat_format_axis(ax: matplotlib.axes.Axes, axis: str = "y", **kwargs):
    """
    Setting the axes in latitude format.

    Parameters
    ----------
    ax : :py:class:`matplotlib.axes.Axes`
        The axes to which the boundary will be applied.
    axis: {'x', 'y'}, default: 'y'
        The axis to which the parameters are applied.
    **kwargs
        Additional keyword arguments to wrapped :py:func:`matplotlib.axis.Axis.set_major_formatter <matplotlib:matplotlib.axis.Axis.set_major_formatter>`.
    """
    if axis == "x":
        axis = ax.xaxis
    elif axis == "y":
        axis = ax.yaxis
    else:
        raise ValueError("`axis`  should be `x` or `y`.")

    axis.set_major_formatter(geoticker.LatitudeFormatter(**kwargs))


def set_p_format_axis(
    ax: matplotlib.axes.Axes,
    axis: str = "y",
    axis_limits: tuple = (1000, 100),
    ticker_step: float = 100,
):
    """
    Setting the axes in logarithmic vertical barometric pressure format.

    Parameters
    ----------
    ax : :py:class:`matplotlib.axes.Axes`
        The axes to which the boundary will be applied.
    axis: {'x', 'y'}, default: 'y'
        The axis to which the parameters are applied.
    axis_limits: :py:class:`tuple`, default `(1000, 100)`.
        Assuming that the distribution of coordinates exhibits an isotropic series distribution,
        this item sets the maximum value (near surface air pressure) and the minimum value (near overhead air pressure).
    ticker_step: :py:class:`float`, default `100`.
        Assuming an isotropic series of coordinate distributions, the term sets the tolerance.
    """
    if axis == "x":
        axis = ax.xaxis
        ax.set_xscale("log")
        ax.set_xlim(axis_limits)
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.xaxis.set_major_locator(ticker.MultipleLocator(ticker_step))
    elif axis == "y":
        axis = ax.yaxis
        ax.set_yscale("log")
        ax.set_ylim(axis_limits)
        ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.yaxis.set_major_locator(ticker.MultipleLocator(ticker_step))
    else:
        raise ValueError("`axis`  should be `x` or `y`.")
