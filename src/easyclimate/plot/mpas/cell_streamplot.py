"""
Cell-centered streamline plotting utilities for MPAS meshes.
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from .common import (
    _adjust_global_lon_window_for_regrid,
    _apply_tri_interp_weights,
    _get_cached_tri_interp_weights,
    _get_cached_triangulation,
    _infer_center_lon,
    _is_cartopy_axis,
    _lon_wrap_to_center,
    _normalize_lon_window,
    _vector_to_projection_grid,
)


__all__ = ["plot_cell_streamplot"]


def plot_cell_streamplot(
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
    # interpolation grid
    nx=160,
    ny=100,
    interpolation_padding=None,
    interpolation_padding_cells=2.0,
    min_valid_fraction=0.02,
    use_triangulation_cache=True,
    use_interp_weight_cache=True,
    # cartopy
    transform=None,
    project_to_map=True,
    regrid_shape=None,
    regrid_global_epsilon=1.0e-6,
    # figure
    title=None,
    **streamplot_kwargs,
):
    """
    Plot cell-centered vector wind as streamlines on an MPAS mesh.

    The MPAS cell-centered irregular vectors are first interpolated onto a
    regular lon-lat grid using Matplotlib triangular linear interpolation, then
    passed to ``Axes.streamplot`` / Cartopy GeoAxes ``streamplot``.

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
        Plot extent in degrees. These bounds are required to define the regular
        interpolation grid.
    center_lon : :py:class:`float <float>`, optional
        Longitude center used for wrapping. If None, infer it from the extent or
        valid cell centers.
    input_radians : :py:class:`bool <bool>`, default: True
        If True, MPAS longitude and latitude variables are interpreted as
        radians and converted to degrees.
    lon_name_cell, lat_name_cell : :py:class:`str <str>`
        MPAS cell coordinate variable names.
    nx, ny : :py:class:`int <int>`, default: 160, 100
        Number of longitude and latitude points in the regular interpolation
        grid.
    interpolation_padding : :py:class:`float <float>` or (:py:class:`float <float>`, :py:class:`float <float>`) or :data:`None`, optional
        Extra longitude/latitude degrees used only when selecting source cells
        for triangular interpolation.  Padding keeps the requested output grid
        inside the interpolation hull, reducing blank margins near plot edges.
        If None, padding is inferred from ``interpolation_padding_cells`` and
        the regular grid spacing.

    interpolation_padding_cells : :py:class:`float <float>`, optional
        Number of regular-grid cells used as automatic interpolation padding.

    project_to_map : :py:class:`bool <bool>`, optional
        For Cartopy GeoAxes, transform the regular lon-lat grid and vector
        components to the target map projection before calling streamplot.
        This avoids Cartopy's internal streamplot regridding from clipping a
        shifted PlateCarree window near the map edge.

    regrid_shape : :py:class:`int <int>` or (:py:class:`int <int>`, :py:class:`int <int>`) or :data:`None`, optional
        Regrid vectors onto a regular grid in the target map projection before
        plotting. This is only supported for Cartopy GeoAxes. If provided, this
        projection-grid path takes precedence over ``project_to_map``.

    use_triangulation_cache : :py:class:`bool <bool>`, optional
        Cache the geometry-only triangulation for repeated calls with the same
        dataset, extent, padding, and interpolation grid. This avoids rebuilding
        Qhull triangulation and TriFinder state when plotting multiple levels or
        time slices on the same MPAS mesh.

    use_interp_weight_cache : :py:class:`bool <bool>`, optional
        Cache target-grid barycentric interpolation weights. This avoids repeated
        ``LinearTriInterpolator`` setup and plane-coefficient calculations for
        repeated calls with the same geometry and interpolation grid.
    transform : :py:class:`cartopy.crs.CRS <cartopy.crs.CRS>`, optional
        Coordinate reference system of the input coordinates for Cartopy axes.
    min_valid_fraction : :py:class:`float <float>`, default: 0.02
        Minimum fraction of valid interpolated grid points required before
        plotting.
    title : :py:class:`str <str>`, optional
        Axes title. If None, use a default streamline title.
    **streamplot_kwargs
        Additional keyword arguments passed to ``Axes.streamplot``.

    Returns
    -------
    :py:class:`matplotlib.streamplot.StreamplotSet <matplotlib.streamplot.StreamplotSet>`
        Streamplot set returned by ``Axes.streamplot``.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        cell_streamplot
    """

    default_streamplot_kwargs = dict(
        density=1.5,
        color="black",
        linewidth=1.0,
        arrowsize=1.0,
    )
    default_streamplot_kwargs.update(streamplot_kwargs)
    streamplot_kwargs = default_streamplot_kwargs

    if lon_min is None or lon_max is None or lat_min is None or lat_max is None:
        raise ValueError(
            "Please provide lon_min, lon_max, lat_min, lat_max. "
            "Streamplot needs an explicit regular interpolation grid."
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

    lon_grid = np.linspace(lon_min_plot, lon_max_plot, int(nx))
    lat_grid = np.linspace(float(lat_min), float(lat_max), int(ny))
    xx, yy = np.meshgrid(lon_grid, lat_grid)

    if interpolation_padding is None:
        dlon_grid = (lon_max_plot - lon_min_plot) / max(int(nx) - 1, 1)
        dlat_grid = (float(lat_max) - float(lat_min)) / max(int(ny) - 1, 1)
        pad_lon = float(interpolation_padding_cells) * dlon_grid
        pad_lat = float(interpolation_padding_cells) * dlat_grid
    else:
        try:
            pad_lon, pad_lat = interpolation_padding
        except TypeError:
            pad_lon = pad_lat = interpolation_padding

        pad_lon = float(pad_lon)
        pad_lat = float(pad_lat)

    coord_valid = np.isfinite(lon) & np.isfinite(lat)

    mask = (
        valid
        & coord_valid
        & (lon >= lon_min_plot - pad_lon)
        & (lon <= lon_max_plot + pad_lon)
        & (lat >= float(lat_min) - pad_lat)
        & (lat <= float(lat_max) + pad_lat)
    )

    x = lon[mask]
    y = lat[mask]
    uu = u[mask]
    vv = v[mask]

    if x.size < 3:
        raise RuntimeError(
            "At least 3 valid cell vectors are required for streamplot interpolation."
        )

    if ax is None:
        ax = plt.gca()

    cartopy_axis = _is_cartopy_axis(ax)
    sp_kwargs = dict(streamplot_kwargs)

    if regrid_shape is not None and not cartopy_axis:
        raise ValueError(
            "`regrid_shape` is only supported when `ax` is a cartopy GeoAxes."
        )

    needs_lonlat_grid_interp = not (cartopy_axis and regrid_shape is not None)

    if needs_lonlat_grid_interp:
        triangulation_cache_key = None

        if use_triangulation_cache:
            triangulation_cache_key = (
                "cell_streamplot_triangulation",
                id(ds),
                lon_name_cell,
                lat_name_cell,
                input_radians,
                float(center_lon),
                float(lon_min_plot),
                float(lon_max_plot),
                float(lat_min),
                float(lat_max),
                float(pad_lon),
                float(pad_lat),
                int(nx),
                int(ny),
                int(lon.size),
                int(x.size),
            )

        triang = _get_cached_triangulation(
            x,
            y,
            cache_key=triangulation_cache_key,
            use_cache=use_triangulation_cache,
            precompute_trifinder=True,
        )

        interp_weight_cache_key = None

        if use_interp_weight_cache and triangulation_cache_key is not None:
            interp_weight_cache_key = (
                "cell_streamplot_interp_weights",
                triangulation_cache_key,
            )

        interp_weights = _get_cached_tri_interp_weights(
            triang,
            x,
            y,
            xx,
            yy,
            cache_key=interp_weight_cache_key,
            use_cache=use_interp_weight_cache,
        )

        u_interp = _apply_tri_interp_weights(uu, interp_weights)
        v_interp = _apply_tri_interp_weights(vv, interp_weights)

        u_grid = np.ma.masked_invalid(u_interp)
        v_grid = np.ma.masked_invalid(v_interp)
        valid_fraction = np.ma.count(u_grid) / u_grid.size

        if valid_fraction < min_valid_fraction:
            raise RuntimeError(
                "Too few regular-grid vectors were interpolated. Check extent or input data."
            )

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
            grid_x, grid_y, u_plot, v_plot = _vector_to_projection_grid(
                transform,
                ax.projection,
                regrid_shape,
                x,
                y,
                uu,
                vv,
                target_extent=ax.get_extent(crs=ax.projection),
                lon_min=lon_min_plot,
                lon_max=lon_max_plot,
                lat_min=lat_min,
                lat_max=lat_max,
                center_lon=center_lon,
                cyclic=True,
            )

            lon_grid = grid_x[0, :]
            lat_grid = grid_y[:, 0]
            u_grid = np.ma.masked_invalid(u_plot)
            v_grid = np.ma.masked_invalid(v_plot)
        elif project_to_map:
            target_xyz = ax.projection.transform_points(transform, xx, yy)
            target_x = target_xyz[..., 0]
            target_y = target_xyz[..., 1]

            rectilinear = (
                np.all(np.isfinite(target_x))
                and np.all(np.isfinite(target_y))
                and np.allclose(target_x, target_x[0, :][None, :])
                and np.allclose(target_y, target_y[:, 0][:, None])
            )

            if rectilinear:
                u_filled = np.ma.filled(u_grid, np.nan)
                v_filled = np.ma.filled(v_grid, np.nan)
                u_plot, v_plot = ax.projection.transform_vectors(
                    transform,
                    xx,
                    yy,
                    u_filled,
                    v_filled,
                )

                x_plot = target_x[0, :]
                y_plot = target_y[:, 0]

                if x_plot[1] < x_plot[0]:
                    x_plot = x_plot[::-1]
                    u_plot = u_plot[:, ::-1]
                    v_plot = v_plot[:, ::-1]

                if y_plot[1] < y_plot[0]:
                    y_plot = y_plot[::-1]
                    u_plot = u_plot[::-1, :]
                    v_plot = v_plot[::-1, :]

                u_grid = np.ma.masked_invalid(u_plot)
                v_grid = np.ma.masked_invalid(v_plot)
                lon_grid = x_plot
                lat_grid = y_plot
            else:
                sp_kwargs["transform"] = transform
        else:
            sp_kwargs["transform"] = transform
    else:
        ax.set_xlim(lon_min_plot, lon_max_plot)
        ax.set_ylim(lat_min, lat_max)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

    sp = ax.streamplot(lon_grid, lat_grid, u_grid, v_grid, **sp_kwargs)

    if title is None:
        title = "Cell-centered wind streamlines"

    ax.set_title(title)

    return sp
