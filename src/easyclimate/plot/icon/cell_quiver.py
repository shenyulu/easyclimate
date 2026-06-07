"""
Quiver plots for ICON cell winds.
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from .common import (
    _is_cartopy_axis,
    _prepare_vector_points,
    _set_cartopy_extent,
    thin_vectors_by_bins,
)


__all__ = ["plot_cell_quiver"]


def plot_cell_quiver(
    ds,
    u_da,
    v_da,
    *,
    ax=None,
    # extent
    lon_min=None,
    lon_max=None,
    lat_min=None,
    lat_max=None,
    center_lon=None,
    # coordinate
    input_radians=True,
    lon_name_cell="clon",
    lat_name_cell="clat",
    # thinning
    nx_bins=35,
    ny_bins=20,
    thin_method="nearest_center",
    use_thinning_cache=True,
    # quiver key
    add_quiverkey=True,
    quiverkey_value=None,
    quiverkey_label=None,
    quiverkey_x=0.88,
    quiverkey_y=1.04,
    # cartopy
    transform=None,
    # figure
    title=None,
    **quiver_kwargs,
):
    """
    Plot cell-centered vector wind on an ICON mesh.

    ICON cell-centered vectors are thinned into regular longitude-latitude bins
    before they are passed to ``Axes.quiver``. This keeps vector density
    readable on native ICON triangular meshes.

    Parameters
    ----------
    ds : :py:class:`xarray.Dataset <xarray.Dataset>`
        ICON dataset containing cell coordinate variables.
    u_da, v_da : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`numpy.ndarray <numpy.ndarray>`
        Cell-centered zonal/eastward and meridional/northward vector
        components. Each input must be one-dimensional on ``ncells`` after
        selecting time or vertical dimensions.
    ax : :py:class:`matplotlib.axes.Axes <matplotlib.axes.Axes>`, optional
        Axes on which to draw. By default, use the current axes.
    lon_min, lon_max, lat_min, lat_max : :py:class:`float <float>`
        Plot extent in degrees. These bounds are required so vector thinning is
        restricted to the requested ICON domain.
    center_lon : :py:class:`float <float>`, optional
        Longitude center used for wrapping. If None, infer it from the extent or
        valid cell centers.
    input_radians : :py:class:`bool <bool>`, default: True
        If True, ICON longitude and latitude variables are interpreted as
        radians and converted to degrees.
    lon_name_cell, lat_name_cell : :py:class:`str <str>`
        ICON cell coordinate variable names.
    nx_bins, ny_bins : :py:class:`int <int>`, default: 35, 20
        Number of longitude and latitude bins used to thin vectors.
    thin_method : {"nearest_center", "max_speed", "mean"}, default: "nearest_center"
        Rule used to choose or aggregate one vector per bin.
    use_thinning_cache : :py:class:`bool <bool>`, default: True
        If True, cache coordinate-only bin metadata for repeated calls on the
        same ICON mesh and extent.
    add_quiverkey : :py:class:`bool <bool>`, default: True
        If True, add a reference quiver key.
    quiverkey_value : :py:class:`float <float>`, optional
        Reference vector magnitude for the quiver key. If None, use the 75th
        percentile of thinned vector speeds.
    quiverkey_label : :py:class:`str <str>`, optional
        Label for the quiver key. If None, build a metres-per-second label from
        ``quiverkey_value``.
    quiverkey_x, quiverkey_y : :py:class:`float <float>`, default: 0.88, 1.04
        Quiver-key position in axes coordinates.
    transform : :py:class:`cartopy.crs.CRS <cartopy.crs.CRS>`, optional
        Coordinate reference system of the input coordinates for Cartopy axes.
    title : :py:class:`str <str>`, optional
        Axes title. If None, use a default wind title.
    **quiver_kwargs
        Additional keyword arguments passed to ``Axes.quiver``.

    Returns
    -------
    :py:class:`matplotlib.quiver.Quiver <matplotlib.quiver.Quiver>`
        Quiver object returned by ``Axes.quiver``.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        icon_cell_quiver
    """
    default_quiver_kwargs = dict(
        scale=None,
        scale_units=None,
        width=0.0025,
        color="black",
        alpha=0.9,
        pivot="middle",
    )
    default_quiver_kwargs.update(quiver_kwargs)
    quiver_kwargs = default_quiver_kwargs

    geo = _prepare_vector_points(
        ds,
        u_da,
        v_da,
        lon_min=lon_min,
        lon_max=lon_max,
        lat_min=lat_min,
        lat_max=lat_max,
        center_lon=center_lon,
        input_radians=input_radians,
        lon_name_cell=lon_name_cell,
        lat_name_cell=lat_name_cell,
    )

    thinning_cache_key = None
    if use_thinning_cache:
        thinning_cache_key = (
            "icon_cell_quiver_bins",
            id(ds),
            lon_name_cell,
            lat_name_cell,
            input_radians,
            geo["center_lon"],
            geo["lon_min_plot"],
            geo["lon_max_plot"],
            geo["lat_min_plot"],
            geo["lat_max_plot"],
            int(nx_bins),
            int(ny_bins),
            int(geo["lon"].size),
        )

    xq, yq, uq, vq, sq = thin_vectors_by_bins(
        geo["x"],
        geo["y"],
        geo["u"],
        geo["v"],
        lon_min=geo["lon_min_plot"],
        lon_max=geo["lon_max_plot"],
        lat_min=geo["lat_min_plot"],
        lat_max=geo["lat_max_plot"],
        nx_bins=nx_bins,
        ny_bins=ny_bins,
        method=thin_method,
        cache_key=thinning_cache_key,
        use_cache=use_thinning_cache,
    )

    if xq.size == 0:
        raise RuntimeError("No vectors left after thinning.")

    if ax is None:
        ax = plt.gca()

    cartopy_axis = _is_cartopy_axis(ax)
    q_kwargs = dict(quiver_kwargs)

    if q_kwargs.get("scale") is None:
        q_kwargs.pop("scale", None)
    if q_kwargs.get("scale_units") is None:
        q_kwargs.pop("scale_units", None)

    if cartopy_axis:
        if transform is None:
            transform = ccrs.PlateCarree()

        _set_cartopy_extent(
            ax,
            geo["lon_min_plot"],
            geo["lon_max_plot"],
            geo["lat_min_plot"],
            geo["lat_max_plot"],
        )
        q_kwargs["transform"] = transform
    else:
        ax.set_xlim(geo["lon_min_plot"], geo["lon_max_plot"])
        ax.set_ylim(geo["lat_min_plot"], geo["lat_max_plot"])
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

    q = ax.quiver(xq, yq, uq, vq, **q_kwargs)

    if add_quiverkey:
        if quiverkey_value is None:
            quiverkey_value = float(np.nanpercentile(sq, 75))

        if quiverkey_label is None:
            quiverkey_label = f"{quiverkey_value:.2g} m s$^{{-1}}$"

        ax.quiverkey(
            q,
            X=quiverkey_x,
            Y=quiverkey_y,
            U=quiverkey_value,
            label=quiverkey_label,
            labelpos="E",
        )

    if title is None:
        title = "Cell-centered wind"

    ax.set_title(title)

    return q
