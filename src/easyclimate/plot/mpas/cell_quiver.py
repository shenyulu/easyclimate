"""
Quiver plots for MPAS cell winds.
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
    lon_name_cell="lonCell",
    lat_name_cell="latCell",
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
    regrid_shape=None,
    regrid_global_epsilon=1.0e-6,
    # figure
    title=None,
    **quiver_kwargs,
):
    """
    Plot cell-centered vector wind on MPAS mesh.

    Parameters
    ----------
    ds : :py:class:`xarray.Dataset <xarray.Dataset>`
        MPAS dataset containing cell coordinate variables.
    u_da, v_da : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`numpy.ndarray <numpy.ndarray>`
        Cell-centered zonal/eastward and meridional/northward vector
        components. Each input must be one-dimensional on ``nCells`` after
        selection.
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
    regrid_shape : :py:class:`int <int>` or (:py:class:`int <int>`, :py:class:`int <int>`) or :data:`None`, optional
        Regrid vectors onto a regular grid in the target map projection before
        plotting. This is only supported for Cartopy GeoAxes.
    regrid_global_epsilon : :py:class:`float <float>`, default: 1.0e-6
        Small longitude shrink applied to full-width global regridding windows.
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

        cell_quiver
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

    if lon_min is None or lon_max is None or lat_min is None or lat_max is None:
        raise ValueError(
            "Please provide lon_min, lon_max, lat_min, lat_max. "
            "This makes thinning robust for variable-resolution meshes."
        )

    # -----------------------------
    # Read vectors
    # -----------------------------
    u = np.asarray(u_da.values if hasattr(u_da, "values") else u_da)
    v = np.asarray(v_da.values if hasattr(v_da, "values") else v_da)

    if u.ndim != 1 or v.ndim != 1:
        raise ValueError(
            f"u_da and v_da must be 1D on nCells. Got {u.shape}, {v.shape}."
        )

    if u.shape != v.shape:
        raise ValueError(f"u_da and v_da shape mismatch: {u.shape}, {v.shape}")

    # -----------------------------
    # Coordinates
    # -----------------------------
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

    # -----------------------------
    # Thin vectors
    # -----------------------------
    thinning_cache_key = None

    if use_thinning_cache:
        thinning_cache_key = (
            "cell_quiver_bins",
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

    xq, yq, uq, vq, sq = thin_vectors_by_bins(
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

    if xq.size == 0:
        raise RuntimeError("No vectors left after thinning.")

    # -----------------------------
    # Axis
    # -----------------------------
    if ax is None:
        ax = plt.gca()

    fig = ax.figure

    cartopy_axis = _is_cartopy_axis(ax)

    q_kwargs = dict(quiver_kwargs)

    if q_kwargs.get("scale") is None:
        q_kwargs.pop("scale", None)

    if q_kwargs.get("scale_units") is None:
        q_kwargs.pop("scale_units", None)

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
            xq, yq, uq, vq = _vector_to_projection_grid(
                transform,
                ax.projection,
                regrid_shape,
                xq,
                yq,
                uq,
                vq,
                target_extent=ax.get_extent(crs=ax.projection),
                lon_min=lon_min_plot,
                lon_max=lon_max_plot,
                lat_min=lat_min,
                lat_max=lat_max,
                center_lon=center_lon,
                cyclic=True,
            )
        else:
            q_kwargs["transform"] = transform

    else:
        if regrid_shape is not None:
            raise ValueError(
                "`regrid_shape` is only supported when `ax` is a cartopy GeoAxes."
            )

        ax.set_xlim(lon_min_plot, lon_max_plot)
        ax.set_ylim(lat_min, lat_max)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

    # -----------------------------
    # Quiver
    # -----------------------------
    q = ax.quiver(
        xq,
        yq,
        uq,
        vq,
        **q_kwargs,
    )

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
