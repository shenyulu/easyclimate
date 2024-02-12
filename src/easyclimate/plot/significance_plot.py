"""
Mapping areas of significance
"""

from __future__ import annotations

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd


def draw_significant_area_contourf(
    p_value: xr.DataArray,
    thresh: float = 0.05,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    ax: matplotlib.axes.Axes = None,
    hatches: str = "...",
    hatch_colors: str = "k",
    reverse_level_plot: bool = False,
    **kwargs,
) -> matplotlib.contour.QuadContourSet:
    """
    Draw significant area by :py:func:`matplotlib.axes.Axes.contourf<matplotlib.axes.Axes.contourf>`.

    Parameters
    ----------
    p_value: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The p value data.
    thresh: :py:class:`float <float>`.
        The threshold value.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    ax : :py:class:`matplotlib.axes.Axes`, optional.
        Axes on which to plot. By default, use the current axes. Mutually exclusive with `size` and `figsize`.
    hatches: :py:class:`list[str]`, default: `...`
        A list of cross hatch patterns to use on the filled areas. If None, no hatching will be added to the contour. Hatching is supported in the PostScript, PDF, SVG and Agg backends only.
    hatch_colors: :py:class:`str <str>`, default: `k`.
        The colors of the hatches.
    reverse_level_plot: :py:class:`bool<bool>`, default: `False`.
        Whether to reverse the drawing area.
    **kwargs, optional:
        Additional keyword arguments to :py:func:`xarray.plot.contourf<xarray.plot.contourf>`.

    Returns
    -------
    :py:class:`matplotlib.contour.QuadContourSet<matplotlib.contour.QuadContourSet>`.
    """
    if ax == None:
        ax = plt.gca()

    if 0 < thresh < 1:
        pass
    else:
        raise ValueError("The parameter `thresh` should be between 0 and 1.")

    if reverse_level_plot == False:
        hatches_value = [hatches, None]
    elif reverse_level_plot == True:
        hatches_value = [None, hatches]
    else:
        raise ValueError("The parameter `hatches_value` should be bool type.")

    cs = p_value.plot.contourf(
        x=lon_dim,
        y=lat_dim,
        ax=ax,
        levels=[0, thresh],
        hatches=hatches_value,
        colors="none",
        add_colorbar=False,
        **kwargs,
    )

    # Change hatch pattern color
    # see: https://github.com/matplotlib/matplotlib/issues/2789/
    #
    # https://matplotlib.org/stable/api/prev_api_changes/api_changes_3.8.0.html#contourset-collections
    #
    # For each level, we set the color of its hatch
    for i, collection in enumerate(cs.collections):
        collection.set_edgecolor(hatch_colors[i % len(hatch_colors)])

    # Doing this also colors in the box around each level
    # We can remove the colored line around the levels by setting the linewidth to 0
    for collection in cs.collections:
        collection.set_linewidth(0.0)


def get_significance_point(
    p_value: xr.DataArray,
    thresh: float = 0.05,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
) -> pd.DataFrame:
    """
    Obtain longitude and latitude array values that meet the conditions within the threshold from a two-dimensional array of p-values

    Parameters
    ----------
    p_value: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The p value data.
    thresh: :py:class:`float <float>`.
        The threshold value.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

    Returns
    -------
    :py:class:`pandas.DataFrame <pandas.DataFrame>`.
    """
    index = np.where(p_value < thresh)
    point_lat = p_value[lat_dim][index[0]].data
    point_lon = p_value[lon_dim][index[1]].data
    return pd.DataFrame({lat_dim: point_lat, lon_dim: point_lon})


def draw_significant_area_scatter(
    significant_points_dataframe: pd.DataFrame,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    ax: matplotlib.axes.Axes = None,
    **kwargs,
):
    """
    Draw significant area by :py:func:`matplotlib.axes.Axes.scatter<matplotlib.axes.Axes.scatter>`.

    Parameters
    ----------
    significant_points_dataframe: :py:class:`pandas.DataFrame <pandas.DataFrame>`.
        The data contains the significant points, which is obtained by the :py:func:`easyclimate.plot.get_significance_point <easyclimate.plot.get_significance_point>`.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    ax : :py:class:`matplotlib.axes.Axes`, optional
        Axes on which to plot. By default, use the current axes. Mutually exclusive with `size` and `figsize`.
    **kwargs, optional:
        Additional keyword arguments to :py:func:`matplotlib.axes.Axes.scatter<matplotlib.axes.Axes.scatter>`.

        .. attention::
            You must specify `kwargs = {'transform': ccrs.PlateCarree()}` (`import cartopy.crs as ccrs`) in the cartopy `GeoAxes` or `GeoAxesSubplot`, otherwise projection errors may occur.
    """
    if ax == None:
        ax = plt.gca()

    ax.scatter(
        significant_points_dataframe[lon_dim].values,
        significant_points_dataframe[lat_dim].values,
        **kwargs,
    )
