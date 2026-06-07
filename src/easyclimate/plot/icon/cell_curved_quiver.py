"""
Curved quiver plots for ICON cell winds.
"""

from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs

from ..modplot import CurvedQuiverplotSet, velovect
from .common import (
    _interpolate_vectors_to_grid,
    _is_cartopy_axis,
    _normalize_lon_window,
    _prepare_vector_points,
    _set_cartopy_extent,
)


__all__ = ["plot_cell_curved_quiver"]


def plot_cell_curved_quiver(
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
    project_to_map=True,
    # curved quiver
    density=1,
    linewidth=None,
    color=None,
    cmap=None,
    norm=None,
    arrowsize=1,
    arrowstyle="-|>",
    zorder=None,
    start_points=None,
    integration_direction="both",
    grains=15,
    broken_streamlines=True,
    ref_magnitude: float | None = None,
    ref_length: float | None = None,
    min_frac_length: float = 0.0,
    length_norm: Literal["reference", "max", "percentile"] = "reference",
    mask_density: int | tuple[int, int] = 10,
    line_start_stride: int = 1,
    arrow_stride: int = 1,
    min_distance: float = 0.0,
    arrow_head_ratio: float = 1.0,
    arrow_position: float = 1.0,
    # figure
    title=None,
) -> CurvedQuiverplotSet:
    """
    Plot ICON cell-centered vector wind as curved quiver trajectories.

    ICON cell-centered irregular vectors are first interpolated onto a regular
    grid, then drawn with the same curved-quiver engine used by
    :py:func:`easyclimate.plot.curved_quiver_plot.curved_quiver`.

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
        cells for interpolation. If None, padding is inferred from
        ``interpolation_padding_cells`` and the regular grid spacing.
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
    project_to_map : :py:class:`bool <bool>`, default: True
        For Cartopy GeoAxes, transform the regular lon-lat grid and vector
        components to the target projection before drawing when the transformed
        grid remains rectilinear.
    density, linewidth, color, cmap, norm, arrowsize, arrowstyle
        Curved-quiver style options passed to ``velovect``.
    zorder, start_points, integration_direction, grains, broken_streamlines
        Curved-quiver trajectory options passed to ``velovect``.
    ref_magnitude, ref_length, min_frac_length, length_norm
        Options controlling how vector magnitude maps to curved-vector length.
    mask_density, line_start_stride, arrow_stride, min_distance
        Options controlling seed-point density and collision masking.
    arrow_head_ratio, arrow_position
        Arrow-head size and placement controls.
    title : :py:class:`str <str>`, optional
        Axes title. If None, use a default curved-quiver title.

    Returns
    -------
    :py:class:`easyclimate.plot.modplot.CurvedQuiverplotSet <easyclimate.plot.modplot.CurvedQuiverplotSet>`
        Container with line and arrow artists returned by ``velovect``.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        icon_cell_curved_quiver
    """
    if lon_min is None or lon_max is None or lat_min is None or lat_max is None:
        raise ValueError(
            "Please provide lon_min, lon_max, lat_min, lat_max. "
            "Curved quiver needs an explicit regular interpolation grid."
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

    cache_key = (
        "icon_cell_curved_quiver_triangulation",
        id(ds),
        lon_name_cell,
        lat_name_cell,
        input_radians,
        geo["center_lon"],
        lon_min_plot,
        lon_max_plot,
        lat_min_plot,
        lat_max_plot,
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
        lon_min=lon_min_plot,
        lon_max=lon_max_plot,
        lat_min=lat_min_plot,
        lat_max=lat_max_plot,
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
    plot_transform = None

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

        if project_to_map:
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
                plot_transform = ax.transData
            else:
                plot_transform = transform._as_mpl_transform(ax)
        else:
            plot_transform = transform._as_mpl_transform(ax)
    else:
        ax.set_xlim(lon_min_plot, lon_max_plot)
        ax.set_ylim(lat_min_plot, lat_max_plot)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

    obj = velovect(
        ax,
        lon_grid,
        lat_grid,
        u_grid,
        v_grid,
        density=density,
        linewidth=linewidth,
        color=color,
        cmap=cmap,
        norm=norm,
        arrowsize=arrowsize,
        arrowstyle=arrowstyle,
        transform=plot_transform,
        zorder=zorder,
        start_points=start_points,
        integration_direction=integration_direction,
        grains=grains,
        broken_streamlines=broken_streamlines,
        ref_magnitude=ref_magnitude,
        ref_length=ref_length,
        min_frac_length=min_frac_length,
        length_norm=length_norm,
        mask_density=mask_density,
        line_start_stride=line_start_stride,
        arrow_stride=arrow_stride,
        min_distance=min_distance,
        arrow_head_ratio=arrow_head_ratio,
        glyph_mode=False,
        arrow_position=arrow_position,
    )

    if title is None:
        title = "Cell-centered curved wind"

    ax.set_title(title)

    return obj
