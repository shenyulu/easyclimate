"""
Cell-centered barb plotting utilities for MPAS meshes.
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from .common import (
    _adjust_global_lon_window_for_regrid,
    _infer_center_lon,
    _is_cartopy_axis,
    _lon_wrap_to_center,
    _normalize_lon_window,
    _vector_to_projection_grid,
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
    lon_name_cell="lonCell",
    lat_name_cell="latCell",
    # thinning
    nx_bins=35,
    ny_bins=20,
    thin_method="nearest_center",
    use_thinning_cache=True,
    # cartopy
    transform=None,
    regrid_shape=None,
    regrid_global_epsilon=1.0e-6,
    # figure
    title=None,
    **barbs_kwargs,
):
    """
    Plot cell-centered vector wind as wind barbs on an MPAS mesh.

    Parameters
    ----------
    ds : :py:class:`xarray.Dataset <xarray.Dataset>`
        MPAS dataset containing cell coordinate variables.
    u_da, v_da : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`numpy.ndarray <numpy.ndarray>`
        Cell-centered zonal/eastward and meridional/northward vectors, 1D on
        nCells after selection.
    ax : :py:class:`matplotlib.axes.Axes <matplotlib.axes.Axes>`, optional
        Axes on which to draw. By default, use the current axes.
    lon_min, lon_max, lat_min, lat_max : :py:class:`float <float>`
        Plot extent in degrees. These bounds are required so vector thinning is
        robust on variable-resolution meshes.
    center_lon : :py:class:`float <float>`, optional
        Longitude center used for wrapping. If None, infer it from the extent or
        valid cell centers.
    input_radians : :py:class:`bool <bool>`, default: True
        If True, MPAS longitude and latitude variables are interpreted as
        radians and converted to degrees.
    lon_name_cell, lat_name_cell : :py:class:`str <str>`
        MPAS cell coordinate variable names.
    nx_bins, ny_bins : :py:class:`int <int>`, default: 35, 20
        Number of longitude and latitude bins used to thin vectors.
    thin_method : {"nearest_center", "max_speed", "mean"}, default: "nearest_center"
        Rule used to choose or aggregate one vector per bin.
    use_thinning_cache : :py:class:`bool <bool>`, default: True
        If True, cache coordinate-only bin metadata for repeated calls.
    transform : :py:class:`cartopy.crs.CRS <cartopy.crs.CRS>`, optional
        Coordinate reference system of the input coordinates for Cartopy axes.
    regrid_shape : :py:class:`int <int>` or (:py:class:`int <int>`, :py:class:`int <int>`) or :data:`None`, optional
        Regrid vectors onto a regular grid in the target map projection before
        plotting. This is only supported for Cartopy GeoAxes.
    regrid_global_epsilon : :py:class:`float <float>`, default: 1.0e-6
        Small longitude shrink applied to full-width global regridding windows.
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

        cell_barbs
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

    if lon_min is None or lon_max is None or lat_min is None or lat_max is None:
        raise ValueError(
            "Please provide lon_min, lon_max, lat_min, lat_max. "
            "This makes thinning robust for variable-resolution meshes."
        )

    u = np.asarray(u_da.values if hasattr(u_da, "values") else u_da)
    v = np.asarray(v_da.values if hasattr(v_da, "values") else v_da)

    if u.ndim != 1 or v.ndim != 1:
        raise ValueError(
            f"u_da and v_da must be 1D on nCells. Got {u.shape}, {v.shape}."
        )

    if u.shape != v.shape:
        raise ValueError(f"u_da and v_da shape mismatch: {u.shape}, {v.shape}")

    lon_raw = np.asarray(ds[lon_name_cell].values)
    lat = np.asarray(ds[lat_name_cell].values)

    if input_radians:
        lon_raw = np.degrees(lon_raw)
        lat = np.degrees(lat)

    if lon_raw.size != u.size:
        raise ValueError(f"nCells mismatch: lonCell={lon_raw.size}, vector={u.size}")

    valid = np.isfinite(u) & np.isfinite(v)

    if center_lon is None:
        center_lon = _infer_center_lon(
            lon_min=lon_min,
            lon_max=lon_max,
            lon_cell_raw=lon_raw,
            valid=valid,
        )

    lon = _lon_wrap_to_center(lon_raw, center_lon)
    lon_min_plot, lon_max_plot = _normalize_lon_window(lon_min, lon_max)

    coord_valid = np.isfinite(lon) & np.isfinite(lat)
    mask = (
        coord_valid
        & (lon >= lon_min_plot)
        & (lon <= lon_max_plot)
        & (lat >= lat_min)
        & (lat <= lat_max)
    )

    x = lon[mask]
    y = lat[mask]
    uu = u[mask]
    vv = v[mask]

    if x.size == 0:
        raise RuntimeError("No valid cell vectors selected. Check extent.")

    thinning_cache_key = None
    if use_thinning_cache:
        thinning_cache_key = (
            "cell_barbs_bins",
            id(ds),
            lon_name_cell,
            lat_name_cell,
            input_radians,
            center_lon,
            lon_min_plot,
            lon_max_plot,
            float(lat_min),
            float(lat_max),
            int(nx_bins),
            int(ny_bins),
            int(lon.size),
        )

    xb, yb, ub, vb, _ = thin_vectors_by_bins(
        x,
        y,
        uu,
        vv,
        lon_min=lon_min_plot,
        lon_max=lon_max_plot,
        lat_min=lat_min,
        lat_max=lat_max,
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

        extent_lon_min = lon_min_plot
        extent_lon_max = lon_max_plot

        if regrid_shape is not None:
            extent_lon_min, extent_lon_max, _ = _adjust_global_lon_window_for_regrid(
                lon_min_plot,
                lon_max_plot,
                epsilon=regrid_global_epsilon,
            )

        ax.set_extent(
            [extent_lon_min, extent_lon_max, lat_min, lat_max],
            crs=ccrs.PlateCarree(),
        )

        if regrid_shape is not None:
            xb, yb, ub, vb = _vector_to_projection_grid(
                transform,
                ax.projection,
                regrid_shape,
                xb,
                yb,
                ub,
                vb,
                target_extent=ax.get_extent(crs=ax.projection),
                lon_min=lon_min_plot,
                lon_max=lon_max_plot,
                lat_min=lat_min,
                lat_max=lat_max,
                center_lon=center_lon,
                cyclic=True,
            )
        else:
            b_kwargs["transform"] = transform
    else:
        if regrid_shape is not None:
            raise ValueError(
                "`regrid_shape` is only supported when `ax` is a cartopy GeoAxes."
            )

        ax.set_xlim(lon_min_plot, lon_max_plot)
        ax.set_ylim(lat_min, lat_max)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

    b = ax.barbs(xb, yb, ub, vb, **b_kwargs)

    if title is None:
        title = "Cell-centered wind barbs"

    ax.set_title(title)

    return b
