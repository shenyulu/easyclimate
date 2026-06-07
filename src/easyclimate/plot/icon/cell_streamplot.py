"""
Cell-centered streamline plotting utilities for ICON meshes.
"""

import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs

from .common import (
    _interpolate_vectors_to_grid,
    _is_cartopy_axis,
    _normalize_lon_window,
    _prepare_vector_points,
    _set_cartopy_extent,
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
    lon_name_cell="clon",
    lat_name_cell="clat",
    # interpolation grid
    nx=160,
    ny=100,
    interpolation_padding=None,
    interpolation_padding_cells=2.0,
    min_valid_fraction=0.02,
    max_triangle_edge=None,
    use_triangulation_cache=True,
    # cartopy
    transform=None,
    # figure
    title=None,
    **streamplot_kwargs,
):
    """
    Plot cell-centered vector wind as streamlines on an ICON mesh.

    ICON cell-centered irregular vectors are first interpolated onto a regular
    lon-lat grid with SciPy linear interpolation, then passed to
    ``Axes.streamplot`` / Cartopy GeoAxes ``streamplot``.

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
        Plot extent in degrees. These bounds are required to define the regular
        interpolation grid.
    center_lon : :py:class:`float <float>`, optional
        Longitude center used for wrapping. If None, infer it from the extent or
        valid cell centers.
    input_radians : :py:class:`bool <bool>`, default: True
        If True, ICON longitude and latitude variables are interpreted as
        radians and converted to degrees.
    lon_name_cell, lat_name_cell : :py:class:`str <str>`
        ICON cell coordinate variable names.
    nx, ny : :py:class:`int <int>`, default: 160, 100
        Number of longitude and latitude points in the regular interpolation
        grid.
    interpolation_padding : :py:class:`float <float>` or (:py:class:`float <float>`, :py:class:`float <float>`) or :data:`None`, optional
        Extra longitude/latitude degrees used only when selecting source ICON
        cells for interpolation. Padding keeps the requested output grid inside
        the interpolation hull, reducing blank margins near plot edges. If None,
        padding is inferred from ``interpolation_padding_cells`` and the regular
        grid spacing.
    interpolation_padding_cells : :py:class:`float <float>`, default: 2.0
        Number of regular-grid cells used as automatic interpolation padding.
    min_valid_fraction : :py:class:`float <float>`, default: 0.02
        Minimum fraction of valid interpolated grid points required before
        plotting.
    max_triangle_edge : :py:class:`float <float>`, optional
        Maximum accepted source-cell separation for interpolation support.
        Longer source triangles are masked when a triangulation-based cache path
        is used.
    use_triangulation_cache : :py:class:`bool <bool>`, default: True
        If True, cache geometry-only interpolation metadata for repeated calls
        on the same ICON mesh and grid.
    transform : :py:class:`cartopy.crs.CRS <cartopy.crs.CRS>`, optional
        Coordinate reference system of the input coordinates for Cartopy axes.
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

        icon_cell_streamplot
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

    lon_min_norm, lon_max_norm = _normalize_lon_window(lon_min, lon_max)
    lat_min_norm = float(lat_min)
    lat_max_norm = float(lat_max)

    if interpolation_padding is None:
        dlon_grid = (lon_max_norm - lon_min_norm) / max(int(nx) - 1, 1)
        dlat_grid = (lat_max_norm - lat_min_norm) / max(int(ny) - 1, 1)
        pad_lon = float(interpolation_padding_cells) * dlon_grid
        pad_lat = float(interpolation_padding_cells) * dlat_grid
    else:
        try:
            pad_lon, pad_lat = interpolation_padding
        except TypeError:
            pad_lon = pad_lat = interpolation_padding

        pad_lon = float(pad_lon)
        pad_lat = float(pad_lat)

    pad_lon = max(pad_lon, 0.0)
    pad_lat = max(pad_lat, 0.0)

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
        padding_lon=pad_lon,
        padding_lat=pad_lat,
    )

    lon_min_plot = geo["lon_min_plot"]
    lon_max_plot = geo["lon_max_plot"]
    lat_min_plot = geo["lat_min_plot"]
    lat_max_plot = geo["lat_max_plot"]

    lon_min_interp = lon_min_plot - pad_lon
    lon_max_interp = lon_max_plot + pad_lon
    lat_min_interp = lat_min_plot - pad_lat
    lat_max_interp = lat_max_plot + pad_lat

    cache_key = (
        "icon_cell_streamplot_triangulation",
        id(ds),
        lon_name_cell,
        lat_name_cell,
        input_radians,
        geo["center_lon"],
        lon_min_plot,
        lon_max_plot,
        lat_min_plot,
        lat_max_plot,
        lon_min_interp,
        lon_max_interp,
        lat_min_interp,
        lat_max_interp,
        pad_lon,
        pad_lat,
        int(nx),
        int(ny),
        int(geo["lon"].size),
        int(geo["x"].size),
        max_triangle_edge,
    )

    lon_grid, lat_grid, xx, yy, u_grid, v_grid = _interpolate_vectors_to_grid(
        geo["x"],
        geo["y"],
        geo["u"],
        geo["v"],
        lon_min=lon_min_interp,
        lon_max=lon_max_interp,
        lat_min=lat_min_interp,
        lat_max=lat_max_interp,
        nx=nx,
        ny=ny,
        max_triangle_edge=max_triangle_edge,
        cache_key=cache_key,
        use_triangulation_cache=use_triangulation_cache,
        min_valid_fraction=min_valid_fraction,
    )

    if ax is None:
        ax = plt.gca()

    cartopy_axis = _is_cartopy_axis(ax)

    if cartopy_axis:
        if transform is None:
            transform = ccrs.PlateCarree()

        _set_cartopy_extent(
            ax,
            lon_min_plot,
            lon_max_plot,
            lat_min_plot,
            lat_max_plot,
        )

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

            lon_grid = x_plot
            lat_grid = y_plot
            u_grid = np.ma.masked_invalid(u_plot)
            v_grid = np.ma.masked_invalid(v_plot)
        else:
            streamplot_kwargs["transform"] = transform
    else:
        ax.set_xlim(lon_min_plot, lon_max_plot)
        ax.set_ylim(lat_min_plot, lat_max_plot)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

    sp = ax.streamplot(lon_grid, lat_grid, u_grid, v_grid, **streamplot_kwargs)

    if title is None:
        title = "Cell-centered wind streamlines"

    ax.set_title(title)

    return sp
