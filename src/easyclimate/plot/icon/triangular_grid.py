"""
Grid-line plots for ICON native triangular meshes.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import cartopy.crs as ccrs

from .common import (
    _rad2deg,
    _lon_wrap_to_center,
    _normalize_lon_window,
    _infer_center_lon,
    _infer_extent_from_points,
    _is_cartopy_axis,
    _set_cartopy_extent,
)


__all__ = ["plot_triangular_grid"]


def _split_segments_at_lon_seam(segments, center_lon, max_lon_span):
    """
    Split line segments crossing the longitude seam into short edge stubs.
    """
    if segments.size == 0:
        return segments

    segments = np.asarray(segments, dtype=float)
    lon0 = segments[:, 0, 0]
    lon1 = segments[:, 1, 0]
    dlon = lon1 - lon0
    crossing = np.abs(dlon) > max_lon_span

    if not np.any(crossing):
        return segments

    kept = segments[~crossing]
    split_segments = []
    seam_low = float(center_lon) - 180.0
    seam_high = float(center_lon) + 180.0

    for segment in segments[crossing]:
        x0, y0 = segment[0]
        x1, y1 = segment[1]

        if x0 - x1 > max_lon_span:
            x1_cont = x1 + 360.0
            t = (seam_high - x0) / (x1_cont - x0)
            y_seam = y0 + t * (y1 - y0)
            split_segments.append([[x0, y0], [seam_high, y_seam]])
            split_segments.append([[seam_low, y_seam], [x1, y1]])
        elif x1 - x0 > max_lon_span:
            x1_cont = x1 - 360.0
            t = (seam_low - x0) / (x1_cont - x0)
            y_seam = y0 + t * (y1 - y0)
            split_segments.append([[x0, y0], [seam_low, y_seam]])
            split_segments.append([[seam_high, y_seam], [x1, y1]])

    if not split_segments:
        return kept

    if kept.size == 0:
        return np.asarray(split_segments, dtype=float)

    return np.concatenate([kept, np.asarray(split_segments, dtype=float)], axis=0)


def plot_triangular_grid(
    ds,
    *,
    lon_min=None,
    lon_max=None,
    lat_min=None,
    lat_max=None,
    ax=None,
    transform=None,
    auto_extent=True,
    extent_margin=0.05,
    cell_margin=0.05,
    global_plot=False,
    center_lon=None,
    input_radians=True,
    skip_large_lon_jump=True,
    max_polygon_lon_span=180.0,
    lon_name_cell="clon",
    lat_name_cell="clat",
    lon_bnds_name="clon_bnds",
    lat_bnds_name="clat_bnds",
    deduplicate_edges=True,
    project_to_map=None,
    **linekwargs,
):
    """
    Plot ICON native triangular mesh edges.

    Parameters
    ----------
    ds : :py:class:`xarray.Dataset <xarray.Dataset>`
        ICON dataset containing triangular cell bounds.
    lon_min, lon_max, lat_min, lat_max : :py:class:`float <float>`, optional
        Optional plot extent in degrees. Cross-dateline longitude windows are
        supported.
    ax : :py:class:`matplotlib.axes.Axes <matplotlib.axes.Axes>`, optional
        Axes on which to draw. By default, use the current axes.
    transform : :py:class:`cartopy.crs.CRS <cartopy.crs.CRS>`, optional
        Coordinate reference system of the input coordinates for Cartopy axes.
    auto_extent : :py:class:`bool <bool>`, default: True
        If True and no explicit extent is supplied, infer the extent from ICON
        cell centers.
    extent_margin : :py:class:`float <float>`, default: 0.05
        Fractional margin added to the automatically inferred extent.
    cell_margin : :py:class:`float <float>`, default: 0.05
        Fractional padding used when selecting triangles near the visible
        extent.
    global_plot : :py:class:`bool <bool>`, default: False
        If True, draw a global extent centered on ``center_lon``.
    center_lon : :py:class:`float <float>`, optional
        Longitude center used for wrapping. If None, infer it from the extent or
        valid cell centers.
    input_radians : :py:class:`bool <bool>`, default: True
        If True, ICON longitude and latitude variables are interpreted as
        radians and converted to degrees.
    skip_large_lon_jump : :py:class:`bool <bool>`, default: True
        If True, skip triangles with very large wrapped longitude spans.
    max_polygon_lon_span : :py:class:`float <float>`, default: 180.0
        Maximum accepted wrapped triangle longitude span in degrees.
    lon_name_cell, lat_name_cell, lon_bnds_name, lat_bnds_name : :py:class:`str <str>`
        ICON coordinate and triangular-bound variable names.
    deduplicate_edges : :py:class:`bool <bool>`, default: True
        If True, remove duplicate shared cell edges before drawing.
    project_to_map : :py:class:`bool <bool>` or {"force"}, optional
        If True or None on Cartopy axes, use the fast endpoint-projection path
        only for global-scale plots and use Cartopy's exact path transform for
        regional plots. If False, always let Cartopy transform each path itself.
        If "force", always project all edge endpoints to the target map
        projection in one vectorized operation before creating the collection.
        The forced path avoids Cartopy's per-path Shapely projection overhead,
        but may differ from Cartopy's exact transform on regional projections.
    **linekwargs
        Additional keyword arguments passed to
        :py:class:`matplotlib.collections.LineCollection`.

    Returns
    -------
    :py:class:`matplotlib.collections.LineCollection <matplotlib.collections.LineCollection>`
        Line collection containing selected ICON triangular mesh edges.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        icon_triangular_grid
    """
    linekwargs = {
        "colors": "grey",
        "linewidths": 0.4,
        **linekwargs,
    }

    if ax is None:
        ax = plt.gca()

    cartopy_axis = _is_cartopy_axis(ax)
    projection = getattr(ax, "projection", None)

    if transform is None and cartopy_axis:
        transform = ccrs.PlateCarree()

    lon_bnds_raw = np.asarray(ds[lon_bnds_name].values)
    lat_bnds = np.asarray(ds[lat_bnds_name].values)

    if lon_bnds_raw.ndim != 2 or lat_bnds.ndim != 2:
        raise ValueError("ICON bound variables must be two-dimensional.")

    if lon_bnds_raw.shape != lat_bnds.shape:
        raise ValueError(
            "`lon_bnds_name` and `lat_bnds_name` must have matching shapes."
        )

    if input_radians:
        lon_bnds_raw = _rad2deg(lon_bnds_raw)
        lat_bnds = _rad2deg(lat_bnds)

    n_cells = lon_bnds_raw.shape[0]

    lon_cell_raw = np.asarray(ds[lon_name_cell].values)
    lat_cell = np.asarray(ds[lat_name_cell].values)

    if input_radians:
        lon_cell_raw = _rad2deg(lon_cell_raw)
        lat_cell = _rad2deg(lat_cell)

    if lon_cell_raw.size != n_cells or lat_cell.size != n_cells:
        raise ValueError("`lon_name_cell` and `lat_name_cell` must have length ncells.")

    center_valid = np.isfinite(lon_cell_raw) & np.isfinite(lat_cell)

    if center_lon is None:
        if cartopy_axis and projection is not None:
            center_lon = float(
                getattr(projection, "proj4_params", {}).get("lon_0", 0.0)
            )
        else:
            center_lon = _infer_center_lon(
                lon_min=lon_min,
                lon_max=lon_max,
                lon_cell=lon_cell_raw,
                valid=center_valid,
            )

    lon_cell = _lon_wrap_to_center(lon_cell_raw, center_lon)
    lon_bnds = lon_cell[:, np.newaxis] + _lon_wrap_to_center(
        lon_bnds_raw - lon_cell_raw[:, np.newaxis],
        0.0,
    )

    user_has_extent = (
        lon_min is not None
        and lon_max is not None
        and lat_min is not None
        and lat_max is not None
    )

    if global_plot:
        lon_min_plot = center_lon - 180.0
        lon_max_plot = center_lon + 180.0
        lat_min_plot = -90.0
        lat_max_plot = 90.0
    elif user_has_extent:
        lon_min_plot, lon_max_plot = _normalize_lon_window(lon_min, lon_max)
        lat_min_plot = float(lat_min)
        lat_max_plot = float(lat_max)
    elif auto_extent:
        lon_min_plot, lon_max_plot, lat_min_plot, lat_max_plot = (
            _infer_extent_from_points(
                lon_cell[center_valid],
                lat_cell[center_valid],
                margin=extent_margin,
            )
        )
    else:
        lon_min_plot = lon_max_plot = None
        lat_min_plot = lat_max_plot = None

    candidate = (
        center_valid
        & np.all(np.isfinite(lon_bnds), axis=1)
        & np.all(np.isfinite(lat_bnds), axis=1)
    )

    if skip_large_lon_jump:
        lon_span = np.nanmax(lon_bnds, axis=1) - np.nanmin(lon_bnds, axis=1)
        candidate &= lon_span <= max_polygon_lon_span

    if lon_min_plot is not None and not global_plot:
        dlon_plot = lon_max_plot - lon_min_plot
        dlat_plot = lat_max_plot - lat_min_plot
        lon_pad = cell_margin * dlon_plot
        lat_pad = cell_margin * dlat_plot

        candidate &= (
            (lon_cell >= lon_min_plot - lon_pad)
            & (lon_cell <= lon_max_plot + lon_pad)
            & (lat_cell >= lat_min_plot - lat_pad)
            & (lat_cell <= lat_max_plot + lat_pad)
        )

    selected_lons = lon_bnds[candidate]
    selected_lats = lat_bnds[candidate]

    edge_lons = np.stack(
        [selected_lons, np.roll(selected_lons, -1, axis=1)],
        axis=2,
    )
    edge_lats = np.stack(
        [selected_lats, np.roll(selected_lats, -1, axis=1)],
        axis=2,
    )
    segments = np.stack([edge_lons, edge_lats], axis=-1).reshape(-1, 2, 2)

    if deduplicate_edges and segments.size:
        edge_keys = np.round(segments, decimals=12)
        start = edge_keys[:, 0, :]
        end = edge_keys[:, 1, :]
        swap = (start[:, 0] > end[:, 0]) | (
            (start[:, 0] == end[:, 0]) & (start[:, 1] > end[:, 1])
        )

        if np.any(swap):
            edge_keys = edge_keys.copy()
            edge_keys[swap, 0, :] = end[swap]
            edge_keys[swap, 1, :] = start[swap]

        _, unique_index = np.unique(
            edge_keys.reshape(edge_keys.shape[0], -1),
            axis=0,
            return_index=True,
        )
        segments = segments[np.sort(unique_index)]

    if segments.size == 0:
        raise RuntimeError(
            "No ICON triangular grid edges selected. Check extent or coordinate units."
        )

    if project_to_map not in {None, True, False, "force"}:
        raise ValueError("`project_to_map` must be None, True, False, or 'force'.")

    if project_to_map in {None, True}:
        if lon_min_plot is not None:
            plot_lon_span = lon_max_plot - lon_min_plot
            plot_lat_span = lat_max_plot - lat_min_plot
        else:
            plot_lon_span = np.nanmax(segments[:, :, 0]) - np.nanmin(segments[:, :, 0])
            plot_lat_span = np.nanmax(segments[:, :, 1]) - np.nanmin(segments[:, :, 1])

        project_to_map = bool(
            cartopy_axis
            and transform is not None
            and (global_plot or plot_lon_span >= 300.0 or plot_lat_span >= 150.0)
        )
    elif project_to_map == "force":
        project_to_map = True

    if cartopy_axis and transform is not None and project_to_map:
        transform_center_lon = float(
            getattr(transform, "proj4_params", {}).get("lon_0", center_lon)
        )
        segments = segments.copy()
        segments[:, :, 0] = _lon_wrap_to_center(
            segments[:, :, 0],
            transform_center_lon,
        )
        segments = _split_segments_at_lon_seam(
            segments,
            transform_center_lon,
            max_polygon_lon_span,
        )

        if segments.size == 0:
            raise RuntimeError(
                "No ICON triangular grid edges remain after seam clipping."
            )

        projected = projection.transform_points(
            transform,
            segments[:, :, 0].ravel(),
            segments[:, :, 1].ravel(),
        )[:, :2]
        segments = projected.reshape(segments.shape)
        finite_segments = np.all(np.isfinite(segments), axis=(1, 2))
        segments = segments[finite_segments]
        linekwargs.pop("transform", None)
    elif cartopy_axis and transform is not None:
        linekwargs.setdefault("transform", transform)
    else:
        linekwargs.pop("transform", None)

    if segments.size == 0:
        raise RuntimeError("No ICON triangular grid edges remain after map projection.")

    lc = LineCollection(segments, **linekwargs)
    ax.add_collection(lc)

    if cartopy_axis:
        if global_plot:
            ax.set_global()
        elif lon_min_plot is not None:
            _set_cartopy_extent(
                ax,
                lon_min_plot,
                lon_max_plot,
                lat_min_plot,
                lat_max_plot,
            )
    else:
        if lon_min_plot is not None:
            ax.set_xlim(lon_min_plot, lon_max_plot)
            ax.set_ylim(lat_min_plot, lat_max_plot)
        else:
            ax.autoscale_view()

        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

    return lc
