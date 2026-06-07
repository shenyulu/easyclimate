"""
Shared ICON plotting helpers.
"""

import numpy as np
import matplotlib.tri as mtri
import cartopy.crs as ccrs


_TRIANGULATION_CACHE = {}
_VECTOR_BIN_CACHE = {}


__all__ = [
    "_rad2deg",
    "_lon_wrap_to_center",
    "_normalize_lon_window",
    "_infer_center_lon",
    "_infer_extent_from_points",
    "_is_cartopy_axis",
    "_set_cartopy_extent",
    "_infer_default_cmap",
    "_read_cell_data",
    "_prepare_cell_geometry",
    "_resolve_contour_levels",
    "_mask_bad_triangles",
    "_get_cached_triangulation",
    "_read_cell_vectors",
    "_prepare_vector_points",
    "thin_vectors_by_bins",
    "_interpolate_vectors_to_grid",
]


def _rad2deg(x):
    """
    Convert radians to degrees.
    """
    return np.degrees(x)


def _lon_wrap_to_center(lon, center_lon):
    """
    Wrap longitude to [center_lon - 180, center_lon + 180).
    """
    lon = np.asarray(lon)
    return ((lon - center_lon + 180.0) % 360.0) + center_lon - 180.0


def _normalize_lon_window(lon_min, lon_max):
    """
    Normalize a longitude window into a continuous increasing interval.
    """
    lon_min = float(lon_min)
    lon_max = float(lon_max)

    while lon_max <= lon_min:
        lon_max += 360.0

    return lon_min, lon_max


def _infer_center_lon(lon_min=None, lon_max=None, lon_cell=None, valid=None):
    """
    Infer longitude wrapping center.
    """
    if lon_min is not None and lon_max is not None:
        lon_min_n, lon_max_n = _normalize_lon_window(lon_min, lon_max)
        return 0.5 * (lon_min_n + lon_max_n)

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

    if dlon == 0.0:
        dlon = 1.0
    if dlat == 0.0:
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


def _wrap_lon_180(lon):
    """
    Wrap longitude to the [-180, 180) interval.
    """
    return (lon + 180.0) % 360.0 - 180.0


def _set_cartopy_extent(ax, lon_min, lon_max, lat_min, lat_max):
    """
    Set Cartopy extent for continuous longitude windows.

    ``GeoAxes.set_extent`` can expand windows such as 110E to 220E to a global
    extent on shifted PlateCarree axes.  For PlateCarree projections, convert
    the geographic longitudes into the target projection's native x coordinate
    and set the extent in that coordinate system.
    """
    projection = getattr(ax, "projection", None)

    if projection is not None and projection.__class__.__name__ == "PlateCarree":
        lon_0 = float(getattr(projection, "proj4_params", {}).get("lon_0", 0.0))
        x_min = _wrap_lon_180(float(lon_min) - lon_0)
        x_max = _wrap_lon_180(float(lon_max) - lon_0)

        while x_max <= x_min:
            x_max += 360.0

        ax.set_extent([x_min, x_max, lat_min, lat_max], crs=projection)
    else:
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())


def _infer_default_cmap(values, cmap=None, sequential="viridis", diverging="RdBu_r"):
    """
    Infer a default colormap similar to xarray.plot.
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


def _read_cell_data(da):
    """
    Read one-dimensional ICON cell data and metadata.
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
            f"`da` must be 1D on ncells after selection, got shape {z.shape}."
        )

    return z, units, long_name


def _prepare_cell_geometry(
    ds,
    da,
    *,
    lon_min=None,
    lon_max=None,
    lat_min=None,
    lat_max=None,
    auto_extent=True,
    extent_margin=0.05,
    global_plot=False,
    center_lon=None,
    input_radians=True,
    lon_name_cell="clon",
    lat_name_cell="clat",
):
    """
    Prepare wrapped ICON cell centers and extent metadata.
    """
    z, units, long_name = _read_cell_data(da)

    lon_cell_raw = np.asarray(ds[lon_name_cell].values)
    lat_cell = np.asarray(ds[lat_name_cell].values)

    if input_radians:
        lon_cell_raw = _rad2deg(lon_cell_raw)
        lat_cell = _rad2deg(lat_cell)

    n_cells = lon_cell_raw.size

    if z.size != n_cells:
        raise ValueError(
            f"`da` length does not match ncells: len(da)={z.size}, ncells={n_cells}"
        )

    valid = np.isfinite(z)

    if center_lon is None:
        center_lon = _infer_center_lon(
            lon_min=lon_min,
            lon_max=lon_max,
            lon_cell=lon_cell_raw,
            valid=valid,
        )

    lon_cell = _lon_wrap_to_center(lon_cell_raw, center_lon)

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
                lon_cell[valid],
                lat_cell[valid],
                margin=extent_margin,
            )
        )
    else:
        lon_min_plot = lon_max_plot = None
        lat_min_plot = lat_max_plot = None

    return {
        "z": z,
        "units": units,
        "long_name": long_name,
        "lon_cell_raw": lon_cell_raw,
        "lat_cell": lat_cell,
        "lon_cell": lon_cell,
        "valid": valid,
        "center_lon": center_lon,
        "lon_min_plot": lon_min_plot,
        "lon_max_plot": lon_max_plot,
        "lat_min_plot": lat_min_plot,
        "lat_max_plot": lat_max_plot,
    }


def _resolve_contour_levels(
    values,
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
    values = np.asarray(values)
    levels_is_scalar = np.isscalar(levels)

    if (vmin is None or vmax is None) and levels_is_scalar:
        if symmetric:
            vmax_auto = np.nanpercentile(np.abs(values), percentile)
            vmin_auto = -vmax_auto
        else:
            vmin_auto = np.nanpercentile(values, 100.0 - percentile)
            vmax_auto = np.nanpercentile(values, percentile)

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
        has_under = np.nanmin(values) < level_min
        has_over = np.nanmax(values) > level_max

        if has_under and has_over:
            extend = "both"
        elif has_under:
            extend = "min"
        elif has_over:
            extend = "max"
        else:
            extend = "neither"

    return contour_levels, vmin, vmax, extend


def _mask_bad_triangles(triang, x, y, max_edge=None):
    """
    Mask triangles with too-long edges.
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


def _get_cached_triangulation(
    x,
    y,
    *,
    max_edge=None,
    cache_key=None,
    use_cache=True,
    precompute_trifinder=False,
):
    """
    Build or reuse a triangulation for ICON cell-center contouring.
    """
    if use_cache and cache_key is not None and cache_key in _TRIANGULATION_CACHE:
        return _TRIANGULATION_CACHE[cache_key]

    triang = mtri.Triangulation(x, y)
    mask = _mask_bad_triangles(triang, x, y, max_edge=max_edge)

    if mask is not None:
        triang.set_mask(mask)

    if precompute_trifinder:
        triang.get_trifinder()

    if use_cache and cache_key is not None:
        _TRIANGULATION_CACHE[cache_key] = triang

    return triang


def _read_cell_vectors(u_da, v_da):
    """
    Read one-dimensional ICON cell vector components.
    """
    u = np.asarray(u_da.values if hasattr(u_da, "values") else u_da)
    v = np.asarray(v_da.values if hasattr(v_da, "values") else v_da)

    if u.ndim != 1 or v.ndim != 1:
        raise ValueError(
            f"u_da and v_da must be 1D on ncells. Got {u.shape}, {v.shape}."
        )

    if u.shape != v.shape:
        raise ValueError(f"u_da and v_da shape mismatch: {u.shape}, {v.shape}")

    return u, v


def _prepare_vector_points(
    ds,
    u_da,
    v_da,
    *,
    lon_min,
    lon_max,
    lat_min,
    lat_max,
    center_lon=None,
    input_radians=True,
    lon_name_cell="clon",
    lat_name_cell="clat",
    padding_lon=0.0,
    padding_lat=0.0,
):
    """
    Read, wrap, and subset ICON cell vectors for an explicit plot window.
    """
    if lon_min is None or lon_max is None or lat_min is None or lat_max is None:
        raise ValueError("Please provide lon_min, lon_max, lat_min, lat_max.")

    u, v = _read_cell_vectors(u_da, v_da)

    lon_raw = np.asarray(ds[lon_name_cell].values)
    lat = np.asarray(ds[lat_name_cell].values)

    if input_radians:
        lon_raw = _rad2deg(lon_raw)
        lat = _rad2deg(lat)

    if lon_raw.size != u.size:
        raise ValueError(f"ncells mismatch: clon={lon_raw.size}, vector={u.size}")

    valid = np.isfinite(u) & np.isfinite(v)

    if center_lon is None:
        center_lon = _infer_center_lon(
            lon_min=lon_min,
            lon_max=lon_max,
            lon_cell=lon_raw,
            valid=valid,
        )

    lon = _lon_wrap_to_center(lon_raw, center_lon)
    lon_min_plot, lon_max_plot = _normalize_lon_window(lon_min, lon_max)

    coord_valid = np.isfinite(lon) & np.isfinite(lat)
    mask = (
        valid
        & coord_valid
        & (lon >= lon_min_plot - float(padding_lon))
        & (lon <= lon_max_plot + float(padding_lon))
        & (lat >= float(lat_min) - float(padding_lat))
        & (lat <= float(lat_max) + float(padding_lat))
    )

    x = lon[mask]
    y = lat[mask]
    uu = u[mask]
    vv = v[mask]

    if x.size == 0:
        raise RuntimeError("No valid ICON cell vectors selected. Check extent.")

    return {
        "x": x,
        "y": y,
        "u": uu,
        "v": vv,
        "lon": lon,
        "lat": lat,
        "center_lon": center_lon,
        "lon_min_plot": lon_min_plot,
        "lon_max_plot": lon_max_plot,
        "lat_min_plot": float(lat_min),
        "lat_max_plot": float(lat_max),
    }


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

    groups = {"valid": valid, "bin_id": bin_id, "dist2": dist2}

    if cache_key is not None and use_cache:
        _VECTOR_BIN_CACHE[cache_key] = groups

    return groups


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


def _interpolate_vectors_to_grid(
    x,
    y,
    u,
    v,
    *,
    lon_min,
    lon_max,
    lat_min,
    lat_max,
    nx=160,
    ny=100,
    max_triangle_edge=None,
    cache_key=None,
    use_triangulation_cache=True,
    min_valid_fraction=0.02,
):
    """
    Interpolate irregular ICON cell vectors onto a regular lon-lat grid.
    """
    try:
        from scipy.interpolate import griddata
    except ImportError as exc:
        raise ImportError("ICON vector interpolation requires scipy.") from exc

    if len(x) < 3:
        raise RuntimeError("At least 3 valid ICON cell vectors are required.")

    lon_grid = np.linspace(float(lon_min), float(lon_max), int(nx))
    lat_grid = np.linspace(float(lat_min), float(lat_max), int(ny))
    xx, yy = np.meshgrid(lon_grid, lat_grid)

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)

    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(u) & np.isfinite(v)
    points = np.column_stack([x[valid], y[valid]])
    values_u = u[valid]
    values_v = v[valid]

    if points.shape[0] < 3:
        raise RuntimeError("At least 3 finite ICON vectors are required.")

    _, unique_idx = np.unique(points, axis=0, return_index=True)
    points = points[unique_idx]
    values_u = values_u[unique_idx]
    values_v = values_v[unique_idx]

    if points.shape[0] < 3:
        raise RuntimeError("At least 3 unique ICON vector locations are required.")

    interp_points = np.column_stack([xx.ravel(), yy.ravel()])
    u_interp = griddata(points, values_u, interp_points, method="linear").reshape(
        xx.shape
    )
    v_interp = griddata(points, values_v, interp_points, method="linear").reshape(
        xx.shape
    )

    u_grid = np.ma.masked_invalid(u_interp)
    v_grid = np.ma.masked_invalid(v_interp)

    valid_fraction = min(np.ma.count(u_grid), np.ma.count(v_grid)) / u_grid.size

    if valid_fraction < min_valid_fraction:
        raise RuntimeError(
            "Too few regular-grid vectors were interpolated. Check extent or input data."
        )

    return lon_grid, lat_grid, xx, yy, u_grid, v_grid
