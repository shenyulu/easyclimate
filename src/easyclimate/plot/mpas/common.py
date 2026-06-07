"""
Shared MPAS plotting helpers.
"""

import numpy as np
import cartopy.crs as ccrs


_VERTEX_DUAL_GEOMETRY_CACHE = {}
_CELL_POLYGON_GEOMETRY_CACHE = {}
_TRIANGULATION_CACHE = {}
_VECTOR_BIN_CACHE = {}
_TRI_INTERP_WEIGHT_CACHE = {}


__all__ = [
    "_rad2deg",
    "_lon_wrap_to_center",
    "_infer_center_lon",
    "_infer_extent_from_points",
    "_is_cartopy_axis",
    "_normalize_lon_window",
    "_adjust_global_lon_window_for_regrid",
    "_infer_center_lon_from_window",
    "_mask_bad_triangles",
    "_infer_default_cmap",
    "_read_vertex_contour_data",
    "_prepare_vertex_contour_geometry",
    "_resolve_contour_levels",
    "_decorate_vertex_contour_axis",
    "_get_vertex_dual_geometry",
    "_get_cell_polygon_geometry",
    "_get_cached_triangulation",
    "_get_cached_tri_interp_weights",
    "_apply_tri_interp_weights",
    "_get_cached_vector_bin_groups",
    "_vector_to_projection_grid",
    "thin_vectors_by_bins",
    "wrap_lon_180",
    "unwrap_edge_lon",
    "lon_in_range_cyclic",
    "edge_near_lon_window",
]


def _rad2deg(x):
    """
    Convert radians to degrees.
    """
    return np.degrees(x)


def wrap_lon_180(lon):
    """
    Wrap longitudes to the [-180, 180) interval.
    """
    return (lon + 180.0) % 360.0 - 180.0


def unwrap_edge_lon(lon0, lon1):
    """
    Locally unwrap one edge so ``lon1`` follows the short path from ``lon0``.

    For example, an edge from 179 degrees to -179 degrees is represented as
    179 degrees to 181 degrees, avoiding a long segment across the whole map.
    """
    dlon = (lon1 - lon0 + 180.0) % 360.0 - 180.0
    lon1_unwrapped = lon0 + dlon
    return lon0, lon1_unwrapped


def lon_in_range_cyclic(lon, lon_min, lon_max):
    """
    Test whether a longitude falls inside a cyclic longitude window.

    Supports both regular windows, such as 90 to 160 degrees, and dateline
    crossing windows, such as 90 to -60 degrees.
    """
    lon = lon % 360.0
    lon_min = lon_min % 360.0
    lon_max = lon_max % 360.0

    if lon_min <= lon_max:
        return lon_min <= lon <= lon_max

    return lon >= lon_min or lon <= lon_max


def edge_near_lon_window(lon0, lon1, lon_min, lon_max):
    """
    Test whether an edge is near a cyclic longitude window.

    For short MPAS mesh edges, checking both endpoints is usually sufficient.
    """
    return lon_in_range_cyclic(lon0, lon_min, lon_max) or lon_in_range_cyclic(
        lon1, lon_min, lon_max
    )


def _lon_wrap_to_center(lon, center_lon):
    """
    Wrap longitude to [center_lon - 180, center_lon + 180).

    Example
    -------
    center_lon = -125:
        170E -> -190
        240E -> -120
        -60  -> -60
    """
    lon = np.asarray(lon)
    return ((lon - center_lon + 180.0) % 360.0) + center_lon - 180.0


def _infer_center_lon(
    lon_min=None,
    lon_max=None,
    lon_cell=None,
    lon_cell_raw=None,
    valid=None,
):
    """
    Infer longitude wrapping center.

    Priority
    --------
    1. If lon_min/lon_max are provided, use the midpoint of the normalized
       continuous longitude window.
    2. Otherwise use circular mean of valid cell longitudes.
    3. Fallback to 0.
    """
    if lon_min is not None and lon_max is not None:
        lon_min_n, lon_max_n = _normalize_lon_window(lon_min, lon_max)
        return 0.5 * (lon_min_n + lon_max_n)

    if lon_cell is None and lon_cell_raw is not None:
        lon_cell = lon_cell_raw

    if lon_cell is None:
        return 0.0

    lon = np.asarray(lon_cell)

    if valid is not None:
        lon = lon[valid]

    lon = lon[np.isfinite(lon)]

    if lon.size == 0:
        return 0.0

    lon_rad = np.deg2rad(lon)
    mean_angle = np.arctan2(
        np.nanmean(np.sin(lon_rad)),
        np.nanmean(np.cos(lon_rad)),
    )

    return float(np.rad2deg(mean_angle))


def _infer_extent_from_points(lon, lat, margin=0.05):
    """
    Infer extent from valid lon/lat points.
    """
    lon = np.asarray(lon)
    lat = np.asarray(lat)

    m = np.isfinite(lon) & np.isfinite(lat)
    lon = lon[m]
    lat = lat[m]

    if lon.size == 0:
        raise RuntimeError("Cannot infer extent: no valid lon/lat points.")

    lon_min = float(np.nanmin(lon))
    lon_max = float(np.nanmax(lon))
    lat_min = float(np.nanmin(lat))
    lat_max = float(np.nanmax(lat))

    dlon = lon_max - lon_min
    dlat = lat_max - lat_min

    if dlon == 0:
        dlon = 1.0
    if dlat == 0:
        dlat = 1.0

    lon_min -= margin * dlon
    lon_max += margin * dlon
    lat_min -= margin * dlat
    lat_max += margin * dlat

    lat_min = max(lat_min, -90.0)
    lat_max = min(lat_max, 90.0)

    return lon_min, lon_max, lat_min, lat_max


def _is_cartopy_axis(ax):
    """
    Weak check for Cartopy GeoAxes.
    """
    return hasattr(ax, "projection")


def _normalize_lon_window(lon_min, lon_max):
    """
    Normalize a longitude window into a continuous increasing interval.

    Examples
    --------
    90, -60   -> 90, 300
    -190, -60 -> -190, -60
    170, -120 -> 170, 240
    """
    lon_min = float(lon_min)
    lon_max = float(lon_max)

    while lon_max <= lon_min:
        lon_max += 360.0

    return lon_min, lon_max


def _adjust_global_lon_window_for_regrid(lon_min, lon_max, *, epsilon=1.0e-6):
    """
    Slightly shrink a global longitude window for Cartopy vector regridding.

    Cartopy's GeoAxes quiver/barbs ``regrid_shape`` path builds a regular grid
    in the target projection.  A fully closed 360-degree source extent contains
    both sides of the same cyclic seam, which can create artificial interpolation
    across the projection cut.  Keeping the requested center but making the
    interval open by a tiny epsilon avoids duplicate seam samples while leaving
    the visible extent unchanged for practical plotting purposes.
    """
    lon_min = float(lon_min)
    lon_max = float(lon_max)
    width = lon_max - lon_min

    if width < 360.0 - epsilon:
        return lon_min, lon_max, False

    half_epsilon = 0.5 * float(epsilon)
    return lon_min + half_epsilon, lon_max - half_epsilon, True


def _infer_center_lon_from_window(lon_min, lon_max):
    """
    Infer the longitude wrapping center from a longitude window.
    """
    lon_min_n, lon_max_n = _normalize_lon_window(lon_min, lon_max)
    return 0.5 * (lon_min_n + lon_max_n)


def _mask_bad_triangles(triang, x, y, max_edge=None):
    """
    Mask triangles with too-long edges.
    This helps remove artificial long triangles near dateline or domain edges.
    """
    if max_edge is None:
        return None

    triangles = triang.triangles

    xtri = x[triangles]
    ytri = y[triangles]

    d01 = np.hypot(xtri[:, 0] - xtri[:, 1], ytri[:, 0] - ytri[:, 1])
    d12 = np.hypot(xtri[:, 1] - xtri[:, 2], ytri[:, 1] - ytri[:, 2])
    d20 = np.hypot(xtri[:, 2] - xtri[:, 0], ytri[:, 2] - ytri[:, 0])

    bad = (d01 > max_edge) | (d12 > max_edge) | (d20 > max_edge)
    return bad


def _infer_default_cmap(values, cmap=None, sequential="viridis", diverging="RdBu_r"):
    """
    Infer a default colormap similar to xarray.plot.

    If the user supplied `cmap`, return it unchanged.  Otherwise use a
    diverging colormap when finite data contain both negative and positive
    values, and a sequential colormap for one-sided data.
    """
    if cmap is not None:
        return cmap

    values = np.asarray(values)
    values = values[np.isfinite(values)]

    if values.size == 0:
        return sequential

    if np.nanmin(values) < 0.0 and np.nanmax(values) > 0.0:
        return diverging

    return sequential


def _read_vertex_contour_data(da):
    """
    Read one-dimensional vertex data and metadata.
    """
    if hasattr(da, "values"):
        z = np.asarray(da.values)
        varname = getattr(da, "name", None)
        units = da.attrs.get("units", "")
        long_name = da.attrs.get("long_name", varname or "value")
    else:
        z = np.asarray(da)
        units = ""
        long_name = "value"

    if z.ndim != 1:
        raise ValueError(
            f"`da` must be 1D on nVertices after selection, got shape {z.shape}."
        )

    return z, units, long_name


def _prepare_vertex_contour_geometry(
    ds,
    da,
    *,
    lon_min=None,
    lon_max=None,
    lat_min=None,
    lat_max=None,
    auto_extent=True,
    extent_margin=0.05,
    contour_margin=0.05,
    global_plot=False,
    center_lon=None,
    input_radians=True,
    lon_name_vertex="lonVertex",
    lat_name_vertex="latVertex",
):
    """
    Prepare wrapped vertex coordinates and extent metadata for contour plots.
    """
    z, units, long_name = _read_vertex_contour_data(da)

    lon_vertex_raw = np.asarray(ds[lon_name_vertex].values)
    lat_vertex = np.asarray(ds[lat_name_vertex].values)

    if input_radians:
        lon_vertex_raw = _rad2deg(lon_vertex_raw)
        lat_vertex = _rad2deg(lat_vertex)

    n_vertices = lon_vertex_raw.size

    if z.size != n_vertices:
        raise ValueError(
            f"`da` length does not match nVertices: "
            f"len(da)={z.size}, nVertices={n_vertices}"
        )

    valid = np.isfinite(z)

    if center_lon is None:
        center_lon = _infer_center_lon(
            lon_min=lon_min,
            lon_max=lon_max,
            lon_cell=lon_vertex_raw,
            valid=valid,
        )

    lon_vertex = _lon_wrap_to_center(lon_vertex_raw, center_lon)

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
                lon_vertex[valid],
                lat_vertex[valid],
                margin=extent_margin,
            )
        )
    else:
        lon_min_plot = lon_max_plot = None
        lat_min_plot = lat_max_plot = None

    mask = np.isfinite(lon_vertex) & np.isfinite(lat_vertex) & valid

    if lon_min_plot is not None and not global_plot:
        dlon_plot = lon_max_plot - lon_min_plot
        dlat_plot = lat_max_plot - lat_min_plot
        lon_pad = contour_margin * dlon_plot
        lat_pad = contour_margin * dlat_plot

        mask &= (
            (lon_vertex >= lon_min_plot - lon_pad)
            & (lon_vertex <= lon_max_plot + lon_pad)
            & (lat_vertex >= lat_min_plot - lat_pad)
            & (lat_vertex <= lat_max_plot + lat_pad)
        )

    idx = np.flatnonzero(mask)
    x = lon_vertex[idx]
    y = lat_vertex[idx]
    zz = z[idx]

    if x.size < 3:
        raise RuntimeError("Not enough valid vertices to draw contours.")

    return {
        "z": z,
        "x": x,
        "y": y,
        "zz": zz,
        "idx": idx,
        "units": units,
        "long_name": long_name,
        "center_lon": center_lon,
        "lon_min_plot": lon_min_plot,
        "lon_max_plot": lon_max_plot,
        "lat_min_plot": lat_min_plot,
        "lat_max_plot": lat_max_plot,
    }


def _resolve_contour_levels(
    zz,
    levels,
    *,
    vmin=None,
    vmax=None,
    symmetric=False,
    percentile=98,
    extend=None,
):
    """
    Resolve contour levels, color limits, and extension mode.
    """
    levels_is_scalar = np.isscalar(levels)

    if (vmin is None or vmax is None) and levels_is_scalar:
        if symmetric:
            vmax_auto = np.nanpercentile(np.abs(zz), percentile)
            vmin_auto = -vmax_auto
        else:
            vmin_auto = np.nanpercentile(zz, 100.0 - percentile)
            vmax_auto = np.nanpercentile(zz, percentile)

        if vmin is None:
            vmin = vmin_auto
        if vmax is None:
            vmax = vmax_auto

    contour_levels = levels
    if levels_is_scalar and vmin is not None and vmax is not None:
        n_levels = int(levels)
        if n_levels < 2:
            raise ValueError("`levels` must be at least 2 when it is an integer.")
        contour_levels = np.linspace(vmin, vmax, n_levels)
    elif not levels_is_scalar:
        contour_levels = np.asarray(levels, dtype=float)
        if contour_levels.ndim != 1 or contour_levels.size < 2:
            raise ValueError("`levels` must be a 1D array with at least 2 values.")
        if vmin is None:
            vmin = np.nanmin(contour_levels)
        if vmax is None:
            vmax = np.nanmax(contour_levels)

    if extend is None:
        level_min = np.nanmin(contour_levels)
        level_max = np.nanmax(contour_levels)
        has_under = np.nanmin(zz) < level_min
        has_over = np.nanmax(zz) > level_max

        if has_under and has_over:
            extend = "both"
        elif has_under:
            extend = "min"
        elif has_over:
            extend = "max"
        else:
            extend = "neither"

    return contour_levels, vmin, vmax, extend


def _decorate_vertex_contour_axis(
    ax,
    *,
    cartopy_axis,
    global_plot=False,
    lon_min_plot=None,
    lon_max_plot=None,
    lat_min_plot=None,
    lat_max_plot=None,
    xlabel="Longitude",
    ylabel="Latitude",
    aspect="auto",
):
    """
    Apply extent and axis labels for vertex contour plots.
    """
    if cartopy_axis:
        if global_plot:
            ax.set_global()
        elif lon_min_plot is not None:
            ax.set_extent(
                [lon_min_plot, lon_max_plot, lat_min_plot, lat_max_plot],
                crs=ccrs.PlateCarree(),
            )
    else:
        if lon_min_plot is not None:
            ax.set_xlim(lon_min_plot, lon_max_plot)
            ax.set_ylim(lat_min_plot, lat_max_plot)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        if aspect is not None:
            ax.set_aspect(aspect, adjustable="box")


def thin_vectors_by_bins(
    x,
    y,
    u,
    v,
    *,
    lon_min,
    lon_max,
    lat_min,
    lat_max,
    nx_bins=35,
    ny_bins=20,
    method="nearest_center",
    cache_key=None,
    use_cache=True,
):
    """
    Thin vectors by lon-lat bins.

    method:
        nearest_center: choose vector closest to each bin center
        max_speed: choose strongest vector in each bin
        mean: average all vectors in each bin
    """

    if method not in {"nearest_center", "max_speed", "mean"}:
        raise ValueError("method must be 'nearest_center', 'max_speed', or 'mean'.")

    x = np.asarray(x)
    y = np.asarray(y)
    u = np.asarray(u)
    v = np.asarray(v)

    groups = _get_cached_vector_bin_groups(
        x,
        y,
        lon_min=lon_min,
        lon_max=lon_max,
        lat_min=lat_min,
        lat_max=lat_max,
        nx_bins=nx_bins,
        ny_bins=ny_bins,
        cache_key=cache_key,
        use_cache=use_cache,
    )

    valid = groups["valid"] & np.isfinite(u) & np.isfinite(v)

    x_valid = x[valid]
    y_valid = y[valid]
    u_valid = u[valid]
    v_valid = v[valid]

    if x_valid.size == 0:
        empty = np.asarray([], dtype=float)
        return empty, empty, empty, empty, empty

    bin_id = groups["bin_id"][valid]

    if method == "mean":
        n_bins = nx_bins * ny_bins
        count = np.bincount(bin_id, minlength=n_bins)
        used = count > 0

        x_out = (
            np.bincount(bin_id, weights=x_valid, minlength=n_bins)[used] / count[used]
        )
        y_out = (
            np.bincount(bin_id, weights=y_valid, minlength=n_bins)[used] / count[used]
        )
        u_out = (
            np.bincount(bin_id, weights=u_valid, minlength=n_bins)[used] / count[used]
        )
        v_out = (
            np.bincount(bin_id, weights=v_valid, minlength=n_bins)[used] / count[used]
        )

    else:
        if method == "nearest_center":
            dist2 = groups["dist2"][valid]
            order = np.lexsort((dist2, bin_id))
        else:
            speed = np.hypot(u_valid, v_valid)
            order = np.lexsort((-speed, bin_id))

        sorted_bins = bin_id[order]
        first = np.r_[True, sorted_bins[1:] != sorted_bins[:-1]]
        selected = order[first]

        x_out = x_valid[selected]
        y_out = y_valid[selected]
        u_out = u_valid[selected]
        v_out = v_valid[selected]

    speed_out = np.hypot(u_out, v_out)

    return x_out, y_out, u_out, v_out, speed_out


def _get_cached_vector_bin_groups(
    x,
    y,
    *,
    lon_min,
    lon_max,
    lat_min,
    lat_max,
    nx_bins=35,
    ny_bins=20,
    cache_key=None,
    use_cache=True,
):
    """
    Build or retrieve coordinate-only bin metadata for vector thinning.

    The cache intentionally stores only coordinate-derived arrays, so repeated
    plots over the same mesh, extent, wrapping center, and bin layout can reuse
    bin assignment and distance-to-bin-center work while still accepting new
    vector values for every call.
    """
    if cache_key is not None and use_cache:
        cached = _VECTOR_BIN_CACHE.get(cache_key)
        if cached is not None:
            return cached

    x = np.asarray(x)
    y = np.asarray(y)

    lon_edges = np.linspace(lon_min, lon_max, nx_bins + 1)
    lat_edges = np.linspace(lat_min, lat_max, ny_bins + 1)

    ix = np.searchsorted(lon_edges, x, side="right") - 1
    iy = np.searchsorted(lat_edges, y, side="right") - 1

    valid = (
        (ix >= 0)
        & (ix < nx_bins)
        & (iy >= 0)
        & (iy < ny_bins)
        & np.isfinite(x)
        & np.isfinite(y)
    )

    bin_id = np.full(x.shape, -1, dtype=np.intp)
    bin_id[valid] = iy[valid] * nx_bins + ix[valid]

    lon_centers = 0.5 * (lon_edges[:-1] + lon_edges[1:])
    lat_centers = 0.5 * (lat_edges[:-1] + lat_edges[1:])

    dist2 = np.full(x.shape, np.inf, dtype=float)
    dist2[valid] = (x[valid] - lon_centers[ix[valid]]) ** 2 + (
        y[valid] - lat_centers[iy[valid]]
    ) ** 2

    groups = {
        "valid": valid,
        "bin_id": bin_id,
        "dist2": dist2,
    }

    if cache_key is not None and use_cache:
        _VECTOR_BIN_CACHE[cache_key] = groups

    return groups


def _vector_to_projection_grid(
    src_crs,
    target_proj,
    regrid_shape,
    x,
    y,
    u,
    v,
    *,
    target_extent,
    lon_min=None,
    lon_max=None,
    lat_min=None,
    lat_max=None,
    center_lon=None,
    cyclic=True,
):
    """
    Regrid lon-lat vectors onto a regular grid in target projection coordinates.

    This is a local replacement for Cartopy's vector regridding path used by
    GeoAxes quiver/barbs.  It builds the regular grid in target-projection data
    coordinates, transforms those grid points back to source lon-lat, wraps the
    longitudes consistently with the source data, interpolates in the source
    lon-lat space, and finally rotates vector components into the target
    projection frame.

    The important difference from Cartopy's default implementation is that the
    source-space interpolation is done after cyclic longitude wrapping and, for
    global windows, optional ±360-degree point duplication.  This avoids the
    artificial seam interpolation that can appear near 180 degrees in polar
    full-longitude plots.
    """
    try:
        from scipy.interpolate import griddata
    except ImportError as exc:
        raise ImportError("Projection vector regridding requires scipy.") from exc

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)

    try:
        nx, ny = regrid_shape
    except TypeError:
        nx = ny = regrid_shape

    nx = int(nx)
    ny = int(ny)

    if nx <= 1 or ny <= 1:
        raise ValueError("`regrid_shape` must define at least a 2 by 2 grid.")

    x0, x1, y0, y1 = target_extent
    grid_x, grid_y = np.meshgrid(
        np.linspace(x0, x1, nx),
        np.linspace(y0, y1, ny),
    )

    source_xyz = src_crs.transform_points(target_proj, grid_x, grid_y)
    source_x = source_xyz[..., 0]
    source_y = source_xyz[..., 1]

    if center_lon is None:
        if lon_min is not None and lon_max is not None:
            lon_min_n, lon_max_n = _normalize_lon_window(lon_min, lon_max)
            center_lon = 0.5 * (lon_min_n + lon_max_n)
        else:
            center_lon = 0.0

    x_wrapped = _lon_wrap_to_center(x, center_lon)
    source_x_wrapped = _lon_wrap_to_center(source_x, center_lon)

    valid = np.isfinite(x_wrapped) & np.isfinite(y) & np.isfinite(u) & np.isfinite(v)

    x_src = x_wrapped[valid]
    y_src = y[valid]
    u_src = u[valid]
    v_src = v[valid]

    if x_src.size < 3:
        raise RuntimeError("At least 3 valid vectors are required for regridding.")

    if cyclic:
        x_src = np.concatenate([x_src - 360.0, x_src, x_src + 360.0])
        y_src = np.tile(y_src, 3)
        u_src = np.tile(u_src, 3)
        v_src = np.tile(v_src, 3)

    grid_valid = np.isfinite(source_x_wrapped) & np.isfinite(source_y)

    if lon_min is not None and lon_max is not None:
        lon_min_n, lon_max_n = _normalize_lon_window(lon_min, lon_max)
        grid_valid &= (source_x_wrapped >= lon_min_n) & (source_x_wrapped <= lon_max_n)

    if lat_min is not None:
        grid_valid &= source_y >= float(lat_min)

    if lat_max is not None:
        grid_valid &= source_y <= float(lat_max)

    points = np.column_stack([x_src, y_src])
    interp_points = np.column_stack(
        [
            source_x_wrapped.ravel(),
            source_y.ravel(),
        ]
    )

    u_grid = griddata(points, u_src, interp_points, method="linear").reshape(
        grid_x.shape
    )
    v_grid = griddata(points, v_src, interp_points, method="linear").reshape(
        grid_x.shape
    )

    u_grid = np.where(grid_valid, u_grid, np.nan)
    v_grid = np.where(grid_valid, v_grid, np.nan)

    u_target, v_target = target_proj.transform_vectors(
        src_crs,
        source_x_wrapped,
        source_y,
        u_grid,
        v_grid,
    )

    u_target = np.where(grid_valid, u_target, np.nan)
    v_target = np.where(grid_valid, v_target, np.nan)

    return grid_x, grid_y, u_target, v_target


def _get_vertex_dual_geometry(
    cells_on_vertex,
    lon_cell,
    lat_cell,
    *,
    cache_key=None,
    use_cache=True,
):
    """
    Build or retrieve cached vertex-dual polygon geometry.

    The cache stores per-vertex cell indices, polygon bounding boxes, and
    geometry-valid masks.  It intentionally does not store data values or
    extent-dependent selections.
    """
    cells_on_vertex = np.asarray(cells_on_vertex, dtype=int)
    lon_cell = np.asarray(lon_cell)
    lat_cell = np.asarray(lat_cell)

    if cache_key is not None and use_cache:
        cached = _VERTEX_DUAL_GEOMETRY_CACHE.get(cache_key)
        if cached is not None:
            return cached

    valid_cids = cells_on_vertex >= 0
    safe_cids = np.where(valid_cids, cells_on_vertex, 0)

    lon_poly = lon_cell[safe_cids]
    lat_poly = lat_cell[safe_cids]

    finite_coord = np.isfinite(lon_poly) & np.isfinite(lat_poly)
    finite_poly = valid_cids & finite_coord
    n_valid = valid_cids.sum(axis=1)
    geometry_valid = (n_valid >= 3) & np.where(valid_cids, finite_coord, True).all(
        axis=1
    )

    lon_masked = np.where(valid_cids, lon_poly, np.nan)
    lat_masked = np.where(valid_cids, lat_poly, np.nan)

    lon_lo = np.nanmin(lon_masked, axis=1)
    lon_hi = np.nanmax(lon_masked, axis=1)
    lat_lo = np.nanmin(lat_masked, axis=1)
    lat_hi = np.nanmax(lat_masked, axis=1)

    geometry = {
        "safe_cids": safe_cids,
        "valid_cids": valid_cids,
        "n_valid": n_valid,
        "valid": geometry_valid,
        "lon_lo": lon_lo,
        "lon_hi": lon_hi,
        "lat_lo": lat_lo,
        "lat_hi": lat_hi,
    }

    if cache_key is not None and use_cache:
        _VERTEX_DUAL_GEOMETRY_CACHE[cache_key] = geometry

    return geometry


def _get_cell_polygon_geometry(
    vertices_on_cell,
    n_edges_on_cell,
    lon_vertex,
    lat_vertex,
    *,
    cache_key=None,
    use_cache=True,
):
    """
    Build or retrieve cached native cell polygon geometry.

    The geometry is independent of plotted values and extents.  It stores safe
    vertex indices, per-cell polygon bounding boxes, and valid-geometry masks.
    """
    vertices_on_cell = np.asarray(vertices_on_cell, dtype=int)
    n_edges_on_cell = np.asarray(n_edges_on_cell, dtype=int)
    lon_vertex = np.asarray(lon_vertex)
    lat_vertex = np.asarray(lat_vertex)

    if cache_key is not None and use_cache:
        cached = _CELL_POLYGON_GEOMETRY_CACHE.get(cache_key)
        if cached is not None:
            return cached

    max_edges = vertices_on_cell.shape[1]
    edge_index = np.arange(max_edges)[None, :]
    valid_vids = (vertices_on_cell >= 0) & (edge_index < n_edges_on_cell[:, None])
    safe_vids = np.where(valid_vids, vertices_on_cell, 0)

    lon_poly = lon_vertex[safe_vids]
    lat_poly = lat_vertex[safe_vids]

    finite_coord = np.isfinite(lon_poly) & np.isfinite(lat_poly)
    geometry_valid = (n_edges_on_cell >= 3) & np.where(
        valid_vids, finite_coord, True
    ).all(axis=1)

    lon_masked = np.where(valid_vids, lon_poly, np.nan)
    lat_masked = np.where(valid_vids, lat_poly, np.nan)

    lon_lo = np.nanmin(lon_masked, axis=1)
    lon_hi = np.nanmax(lon_masked, axis=1)
    lat_lo = np.nanmin(lat_masked, axis=1)
    lat_hi = np.nanmax(lat_masked, axis=1)

    geometry = {
        "safe_vids": safe_vids,
        "n_edges": n_edges_on_cell,
        "valid": geometry_valid,
        "lon_lo": lon_lo,
        "lon_hi": lon_hi,
        "lat_lo": lat_lo,
        "lat_hi": lat_hi,
    }

    if cache_key is not None and use_cache:
        _CELL_POLYGON_GEOMETRY_CACHE[cache_key] = geometry

    return geometry


def _get_cached_triangulation(
    x,
    y,
    *,
    mask=None,
    max_edge=None,
    cache_key=None,
    use_cache=True,
    precompute_trifinder=False,
):
    """
    Build or retrieve a cached Matplotlib Triangulation.

    This is useful when repeated contour plots use the same wrapped grid,
    extent subset, and triangle quality settings.
    """
    if cache_key is not None and use_cache:
        cached = _TRIANGULATION_CACHE.get(cache_key)
        if cached is not None:
            return cached

    import matplotlib.tri as mtri

    triang = mtri.Triangulation(x, y)

    tri_mask = _mask_bad_triangles(
        triang,
        x,
        y,
        max_edge=max_edge,
    )

    if tri_mask is not None:
        triang.set_mask(tri_mask)

    if precompute_trifinder:
        triang.get_trifinder()

    if cache_key is not None and use_cache:
        _TRIANGULATION_CACHE[cache_key] = triang

    return triang


def _get_cached_tri_interp_weights(
    triang,
    x,
    y,
    xi,
    yi,
    *,
    cache_key=None,
    use_cache=True,
):
    """
    Build or retrieve barycentric interpolation weights for fixed target points.

    The returned mapping depends only on triangulation geometry and target grid,
    not on the interpolated data values.  Reusing it avoids repeated
    ``LinearTriInterpolator`` setup, plane-coefficient calculation, and
    point-location work for repeated streamplot calls on the same grid.
    """
    if cache_key is not None and use_cache:
        cached = _TRI_INTERP_WEIGHT_CACHE.get(cache_key)
        if cached is not None:
            return cached

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xi = np.asarray(xi, dtype=float)
    yi = np.asarray(yi, dtype=float)

    tri_index = triang.get_trifinder()(xi, yi)
    valid = tri_index >= 0

    triangles = triang.triangles
    vertex_indices = np.zeros(tri_index.shape + (3,), dtype=np.intp)
    weights = np.full(tri_index.shape + (3,), np.nan, dtype=float)

    if np.any(valid):
        tri_valid = tri_index[valid]
        verts = triangles[tri_valid]

        x0 = x[verts[:, 0]]
        y0 = y[verts[:, 0]]
        x1 = x[verts[:, 1]]
        y1 = y[verts[:, 1]]
        x2 = x[verts[:, 2]]
        y2 = y[verts[:, 2]]
        xp = xi[valid]
        yp = yi[valid]

        denom = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
        finite = np.isfinite(denom) & (denom != 0.0)

        w0 = np.full(tri_valid.shape, np.nan, dtype=float)
        w1 = np.full(tri_valid.shape, np.nan, dtype=float)
        w2 = np.full(tri_valid.shape, np.nan, dtype=float)

        w0[finite] = (
            (y1[finite] - y2[finite]) * (xp[finite] - x2[finite])
            + (x2[finite] - x1[finite]) * (yp[finite] - y2[finite])
        ) / denom[finite]
        w1[finite] = (
            (y2[finite] - y0[finite]) * (xp[finite] - x2[finite])
            + (x0[finite] - x2[finite]) * (yp[finite] - y2[finite])
        ) / denom[finite]
        w2[finite] = 1.0 - w0[finite] - w1[finite]

        vertex_indices[valid] = verts
        weights[valid, 0] = w0
        weights[valid, 1] = w1
        weights[valid, 2] = w2
        valid = valid.copy()
        valid[valid] = finite

    mapping = {
        "vertex_indices": vertex_indices,
        "weights": weights,
        "valid": valid,
    }

    if cache_key is not None and use_cache:
        _TRI_INTERP_WEIGHT_CACHE[cache_key] = mapping

    return mapping


def _apply_tri_interp_weights(values, mapping):
    """
    Apply cached barycentric interpolation weights to one data field.
    """
    values = np.asarray(values, dtype=float)
    vertex_indices = mapping["vertex_indices"]
    weights = mapping["weights"]
    valid = mapping["valid"]

    out = np.full(valid.shape, np.nan, dtype=float)

    if np.any(valid):
        vals = values[vertex_indices[valid]]
        good = np.isfinite(vals).all(axis=1) & np.isfinite(weights[valid]).all(axis=1)
        interp = np.full(vals.shape[0], np.nan, dtype=float)
        interp[good] = np.sum(vals[good] * weights[valid][good], axis=1)
        out[valid] = interp

    return out
