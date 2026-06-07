"""
Quick processing of geographic axes labels
"""

import matplotlib
import matplotlib.pyplot as plt

__all__ = ["add_geolatitude_label", "add_geolongitude_label"]


def add_geolatitude_label(
    ax: matplotlib.axes.Axes = None,
    **kwargs,
):
    """
    Add a latitude label to a Cartopy PlateCarree axes.

    Parameters
    ----------
    ax : :py:class:`matplotlib.axes.Axes`
        The axes to which the label will be applied.
    **kwargs
        Additional keyword arguments to wrapped :py:meth:`matplotlib.axes.Axes.text <matplotlib:matplotlib.axes.Axes.text>`.
    """
    if ax is None:
        ax = plt.gca()

    kwargs.setdefault("x", -0.07)
    kwargs.setdefault("y", 0.55)
    kwargs.setdefault("s", "")
    kwargs.setdefault("va", "bottom")
    kwargs.setdefault("ha", "center")
    kwargs.setdefault("rotation", "vertical")
    kwargs.setdefault("rotation_mode", "anchor")

    return ax.text(**kwargs)


def add_geolongitude_label(
    ax: matplotlib.axes.Axes = None,
    **kwargs,
):
    """
    Add a longitude label to a Cartopy PlateCarree axes.

    Parameters
    ----------
    ax : :py:class:`matplotlib.axes.Axes`
        The axes to which the label will be applied.
    **kwargs
        Additional keyword arguments to wrapped :py:meth:`matplotlib.axes.Axes.text <matplotlib:matplotlib.axes.Axes.text>`.
    """
    if ax is None:
        ax = plt.gca()

    kwargs.setdefault("x", 0.5)
    kwargs.setdefault("y", -0.2)
    kwargs.setdefault("s", "")
    kwargs.setdefault("va", "bottom")
    kwargs.setdefault("ha", "center")
    kwargs.setdefault("rotation", "vertical")
    kwargs.setdefault("rotation_mode", "anchor")

    return ax.text(**kwargs)
