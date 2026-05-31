"""
Grid-line plots for MPAS Voronoi meshes.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import cartopy.crs as ccrs

from .common import (
    _rad2deg,
    _is_cartopy_axis,
    wrap_lon_180,
    unwrap_edge_lon,
    edge_near_lon_window,
)

__all__ = ["plot_voronoi_grid"]


def plot_voronoi_grid(
    data,
    verticesOnEdge="verticesOnEdge",
    lonVertex="lonVertex",
    latVertex="latVertex",
    lon_min=None,
    lon_max=None,
    lat_min=None,
    lat_max=None,
    ax=None,
    **linekwags,
):
    """
    Plot MPAS Voronoi mesh edges.

    Parameters
    ----------
    data : :py:class:`xarray.Dataset <xarray.Dataset>`
        MPAS dataset containing ``verticesOnEdge``, ``lonVertex``, and
        ``latVertex``. Used as the fallback source when these variables are not
        provided explicitly.
    verticesOnEdge : :py:class:`str <str>`, optional
        Name of the edge-to-vertex connectivity variable using MPAS 1-based
        vertex indices.
    lonVertex, latVertex : :py:class:`str <str>`, optional
        Names of the vertex longitude and latitude variables in radians.
    lon_min, lon_max, lat_min, lat_max : :py:class:`float <float>`, optional
        Optional subset window in degrees. Cross-dateline longitude windows are
        supported.
    ax : :py:class:`matplotlib.axes.Axes <matplotlib.axes.Axes>`, optional
        Axes on which to draw. By default, use the current axes.
    **linekwags
        Additional keyword arguments passed to
        :py:class:`matplotlib.collections.LineCollection`.

    Returns
    -------
    :py:class:`matplotlib.collections.LineCollection <matplotlib.collections.LineCollection>`
        Line collection containing selected Voronoi mesh edges.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        voronoi_grid
    """
    linekwags = {
        "colors": "grey",
        "linewidths": 0.4,
        **linekwags,
    }

    verticesOnEdge = data[verticesOnEdge].values.astype(int) - 1

    lonVertex = wrap_lon_180(_rad2deg(data[lonVertex].values))
    latVertex = _rad2deg(data[latVertex].values)

    segments = []

    use_subset = (
        lon_min is not None
        and lon_max is not None
        and lat_min is not None
        and lat_max is not None
    )

    for e in range(verticesOnEdge.shape[0]):
        v0, v1 = verticesOnEdge[e, :]

        if v0 < 0 or v1 < 0:
            continue

        x0 = lonVertex[v0]
        y0 = latVertex[v0]
        x1 = lonVertex[v1]
        y1 = latVertex[v1]

        # Locally unwrap edges crossing the dateline so they draw as short
        # segments instead of long lines across the globe.
        x0u, x1u = unwrap_edge_lon(x0, x1)

        if use_subset:
            if max(y0, y1) < lat_min or min(y0, y1) > lat_max:
                continue

            if not edge_near_lon_window(x0, x1, lon_min, lon_max):
                continue

        segments.append([(x0u, y0), (x1u, y1)])

    if ax is None:
        ax = plt.gca()

    cartopy_axis = _is_cartopy_axis(ax)

    if cartopy_axis:
        linekwags.setdefault("transform", ccrs.PlateCarree())
    else:
        linekwags.pop("transform", None)

    lc = LineCollection(segments, **linekwags)

    ax.add_collection(lc)

    if segments:
        xy = np.asarray(segments, dtype=float)
        x = xy[:, :, 0]
        y = xy[:, :, 1]

        x_min = float(np.nanmin(x))
        x_max = float(np.nanmax(x))
        y_min = float(np.nanmin(y))
        y_max = float(np.nanmax(y))

        dx = x_max - x_min
        dy = y_max - y_min

        if dx == 0:
            dx = 1.0
        if dy == 0:
            dy = 1.0

        margin = 0.05
        xlim = np.asarray(
            [x_min - margin * dx, x_max + margin * dx],
            dtype=float,
        )
        ylim = np.asarray(
            [y_min - margin * dy, y_max + margin * dy],
            dtype=float,
        )

        if not (np.all(np.isfinite(xlim)) and np.all(np.isfinite(ylim))):
            return lc

        if cartopy_axis:
            extent = [
                max(float(xlim[0]), -180.0),
                min(float(xlim[1]), 180.0),
                max(float(ylim[0]), -90.0),
                min(float(ylim[1]), 90.0),
            ]

            if extent[0] < extent[1] and extent[2] < extent[3]:
                try:
                    ax.set_extent(extent, crs=ccrs.PlateCarree())
                except ValueError:
                    pass
        else:
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)

    return lc
