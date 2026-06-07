"""
Cell-centered barb plotting utilities for ICON meshes.
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from .common import (
    _is_cartopy_axis,
    _prepare_vector_points,
    _set_cartopy_extent,
    thin_vectors_by_bins,
)


__all__ = ["plot_cell_barbs"]


def plot_cell_barbs(
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
    # cartopy
    transform=None,
    # figure
    title=None,
    **barbs_kwargs,
):
    """
    Plot cell-centered vector wind as wind barbs on an ICON mesh.

    ICON cell-centered vectors are thinned into regular longitude-latitude bins
    before they are passed to ``Axes.barbs``. This keeps barb density readable
    on native ICON triangular meshes.

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
    transform : :py:class:`cartopy.crs.CRS <cartopy.crs.CRS>`, optional
        Coordinate reference system of the input coordinates for Cartopy axes.
    title : :py:class:`str <str>`, optional
        Axes title. If None, use a default wind-barb title.
    **barbs_kwargs
        Additional keyword arguments passed to ``Axes.barbs``.

    Returns
    -------
    :py:class:`matplotlib.quiver.Barbs <matplotlib.quiver.Barbs>`
        Barb container returned by ``Axes.barbs``.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        icon_cell_barbs
    """
    default_barbs_kwargs = dict(
        length=5,
        linewidth=0.6,
        color="black",
        alpha=0.9,
        pivot="middle",
    )
    default_barbs_kwargs.update(barbs_kwargs)
    barbs_kwargs = default_barbs_kwargs

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
            "icon_cell_barbs_bins",
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

    xb, yb, ub, vb, _ = thin_vectors_by_bins(
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

    if xb.size == 0:
        raise RuntimeError("No vectors left after thinning.")

    if ax is None:
        ax = plt.gca()

    cartopy_axis = _is_cartopy_axis(ax)
    b_kwargs = dict(barbs_kwargs)

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
        b_kwargs["transform"] = transform
    else:
        ax.set_xlim(geo["lon_min_plot"], geo["lon_max_plot"])
        ax.set_ylim(geo["lat_min_plot"], geo["lat_max_plot"])
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

    b = ax.barbs(xb, yb, ub, vb, **b_kwargs)

    if title is None:
        title = "Cell-centered wind barbs"

    ax.set_title(title)

    return b
