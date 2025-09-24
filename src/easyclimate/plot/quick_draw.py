"""
The quick drawing function
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib
import matplotlib.patches as patches
import numpy as np
import xarray as xr
from ..core.utility import transfer_xarray_lon_from180TO360

__all__ = ["quick_draw_spatial_basemap", "quick_draw_rectangular_box"]


def quick_draw_spatial_basemap(
    nrows: int = 1,
    ncols: int = 1,
    figsize=None,
    central_longitude: float = 0.0,
    draw_labels: str | bool | list | dict = ["bottom", "left"],
    gridlines_color: str = "grey",
    gridlines_alpha: float = 0.5,
    gridlines_linestyle: str = "--",
    coastlines_edgecolor: str = "black",
    coastlines_kwargs: dict = {"lw": 0.5},
):
    """
    Create geographical and spatial base map.

    Parameters
    ----------
    nrows, ncols :py:class:`int <int>`, default: 1
        Number of rows/columns of the subplot grid.
    figsize: (:py:class:`float <float>`, :py:class:`float <float>`)
        Width, height in inches.
    central_longitude: :py:class:`float <float>`, default: 0.
        The central longitude for `cartopy.crs.PlateCarree` projection.
    draw_labels: :py:class:`str <str>` | :py:class:`bool <bool>` | :py:class:`list <list>` | :py:class:`dict <dict>`, default: ["bottom", "left"].
        Toggle whether to draw labels. For finer control, attributes of Gridliner may be modified individually.

        - string: `"x"` or `"y"` to only draw labels of the respective coordinate in the CRS.
        - list: Can contain the side identifiers and/or coordinate types to select which ones to draw. For all labels one would use `["x", "y", "top", "bottom", "left", "right", "geo"]`.
        - dict: The keys are the side identifiers `("top", "bottom", "left", "right")` and the values are the coordinates `("x", "y")`; this way you can precisely decide what kind of label to draw and where. For x labels on the bottom and y labels on the right you could pass in `{"bottom": "x", "left": "y"}`.

        Note that, by default, x and y labels are not drawn on left/right and top/bottom edges respectively unless explicitly requested.
    gridlines_color: :py:class:`str <str>`, default: `grey`.
        The parameter `color` for `ax.gridlines`.
    gridlines_alpha: :py:class:`float <float>`, default: `0.5`.
        The parameter `alpha` for `ax.gridlines`.
    gridlines_linestyle: :py:class:`str <str>`, default: `"--"`.
        The parameter `linestyle` for `ax.gridlines`.
    coastlines_edgecolor: :py:class:`str <str>`, default: `"black"`.
        The parameter `color` for `ax.coastlines`.
    coastlines_kwargs: :py:class:`float <float>`, default: ``{"lw": 0.5}``.
        The kwargs for `ax.coastlines`.

    Returns
    -------
    - fig: :py:class:`Figure <matplotlib:matplotlib.figure.Figure>`
    - ax: :py:class:`Axes <matplotlib:matplotlib.axes.Axes>` or array of Axes: `ax` can be either a single :py:class:`Axes <matplotlib:matplotlib.axes.Axes>` object, or an array of Axes objects if more than one subplot was created. The dimensions of the resulting array can be controlled with the squeeze keyword.

    .. seealso::
        :py:func:`matplotlib.pyplot.subplots <matplotlib:matplotlib.pyplot.subplots>`
        :py:func:`cartopy.mpl.geoaxes.GeoAxes.gridlines <cartopy:cartopy.mpl.geoaxes.GeoAxes.gridlines>`
        :py:func:`cartopy.mpl.geoaxes.GeoAxes.coastlines <cartopy:cartopy.mpl.geoaxes.GeoAxes.coastlines>`
    """

    fig, ax = plt.subplots(
        figsize=figsize,
        nrows=nrows,
        ncols=ncols,
        subplot_kw={
            "projection": ccrs.PlateCarree(central_longitude=central_longitude)
        },
    )

    if nrows == 1 and ncols == 1:
        ax.gridlines(
            draw_labels=draw_labels,
            color=gridlines_color,
            alpha=gridlines_alpha,
            linestyle=gridlines_linestyle,
        )
        ax.coastlines(color=coastlines_edgecolor, **coastlines_kwargs)
    else:
        for axi in ax.flat:
            axi.gridlines(
                draw_labels=draw_labels,
                color=gridlines_color,
                alpha=gridlines_alpha,
                linestyle=gridlines_linestyle,
            )
            axi.coastlines(color=coastlines_edgecolor, **coastlines_kwargs)
    return fig, ax


def quick_draw_rectangular_box(
    lon1: float,
    lon2: float,
    lat1: float,
    lat2: float,
    ax: matplotlib.axes.Axes = None,
    **patches_kwargs,
):
    """
    Create geographical rectangular box.

    Parameters
    ----------
    lon1, lon2: :py:class:`float <float>`.
        Rectangular box longitude point. The applicable value should be between -180 :math:`^\\circ` and 360 :math:`^\\circ`.
        `lon1` and `lon2` must have a certain difference, should not be equal,
        do not strictly require the size relationship between `lon1` and `lon2`.
    lat1, lat2: :py:class:`float <float>`.
        Rectangular box latitude point. The applicable value should be between -90 :math:`^\\circ` and 90 :math:`^\\circ`.
        `lat1` and `lat2` must have a certain difference, should not be equal,
        do not strictly require the size relationship between `lat1` and `lat2`.
    ax : :py:class:`matplotlib.axes.Axes`, optional.
        Axes on which to plot. By default, use the current axes. Mutually exclusive with `size` and `figsize`.
    **patches_kwargs:
        Patch properties. see more in :py:class:`matplotlib.patches.Patch <matplotlib.patches.Patch>`

    .. seealso::
        :py:class:`matplotlib.patches.Rectangle <matplotlib.patches.Rectangle>`
    """
    # Get Axes
    if ax == None:
        ax = plt.gca()
    else:
        pass

    if lon1 > 360 or lon2 > 360 or lon1 < -180 or lon2 < -180:
        raise ValueError(
            "`lon1` or `lon2` should remain between -180 degree to 360 degree."
        )
    if lat1 > 90 or lat2 > 90 or lat1 < -90 or lat2 < -90:
        raise ValueError(
            "`lat1` or `lat2` should remain between -90 degree to 90 degree."
        )

    if lon1 < 0 or lon2 < 0:
        # Transfer from -180-180 degree to 0-360 degree
        data_raw = np.array([lon1, lon2])
        data_transfered = xr.DataArray(data_raw, dims="lon", coords={"lon": data_raw})
        data_transfered = transfer_xarray_lon_from180TO360(data_transfered)
        lon1 = data_transfered["lon"].data[0]
        lon2 = data_transfered["lon"].data[1]

    width = lon2 - lon1
    height = lat2 - lat1

    if width < 0:
        tmp = lon2
        lon2 = lon1
        lon1 = tmp
        width = np.abs(width)
    elif width == 0:
        raise ValueError("`lon1` and `lon2` should not be same!")

    if height < 0:
        tmp = lat2
        lat2 = lat1
        lat1 = tmp
        height = np.abs(height)
    elif height == 0:
        raise ValueError("`lat1` and `lat2` should not be same!")

    rect = patches.Rectangle((lon1, lat1), width, height, **patches_kwargs)
    ax.add_patch(rect)
