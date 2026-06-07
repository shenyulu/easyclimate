"""
Voronoi-style plots for MPAS vertices.
"""

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.collections import PolyCollection
import cartopy.crs as ccrs

from .common import (
    _rad2deg,
    _lon_wrap_to_center,
    _infer_center_lon,
    _infer_extent_from_points,
    _is_cartopy_axis,
    _normalize_lon_window,
    _infer_default_cmap,
    _get_vertex_dual_geometry,
)


__all__ = ["plot_vertex_voronoi"]


def plot_vertex_voronoi(
    ds,
    da,
    *,
    ax=None,
    transform=None,
    # extent control
    lon_min=None,
    lon_max=None,
    lat_min=None,
    lat_max=None,
    auto_extent=True,
    extent_margin=0.05,
    vertex_margin=0.05,
    global_plot=False,
    # coordinate control
    center_lon=None,
    input_radians=True,
    skip_large_lon_jump=True,
    max_polygon_lon_span=180.0,
    # color control
    cmap=None,
    vmin=None,
    vmax=None,
    symmetric=None,
    percentile=98,
    # mesh style
    edgecolor="none",
    linewidth=0.0,
    # map decorations
    add_colorbar=True,
    cbar_kwargs=None,
    title=None,
    # plain matplotlib labels
    xlabel="Longitude",
    ylabel="Latitude",
    aspect="auto",
    # MPAS mesh variable names
    cells_on_vertex_name="cellsOnVertex",
    lon_name_cell="lonCell",
    lat_name_cell="latCell",
    lon_name_vertex="lonVertex",
    lat_name_vertex="latVertex",
    use_geometry_cache=True,
):
    """
    Plot MPAS vertex-centered scalar data on dual polygons.

    Each vertex value is drawn on the dual polygon formed by the surrounding
    cell centers. If ``da`` is cell-centered with length ``nCells``, values are
    first averaged to vertices using ``cellsOnVertex``.

    Parameters
    ----------
    ds : :py:class:`xarray.Dataset <xarray.Dataset>`
        MPAS dataset containing vertex, cell, and connectivity variables.
    da : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`numpy.ndarray <numpy.ndarray>`
        Vertex-centered values on ``nVertices`` or cell-centered values on
        ``nCells``. Cell-centered values are averaged to vertices before
        plotting.
    ax : :py:class:`matplotlib.axes.Axes <matplotlib.axes.Axes>`, optional
        Axes on which to draw. By default, use the current axes.
    transform : :py:class:`cartopy.crs.CRS <cartopy.crs.CRS>`, optional
        Coordinate reference system of the input coordinates for Cartopy axes.
    lon_min, lon_max, lat_min, lat_max : :py:class:`float <float>`, optional
        Plot extent in degrees. Cross-dateline longitude windows are supported.
    auto_extent : :py:class:`bool <bool>`, default: True
        If True and no explicit extent is supplied, infer the extent from valid
        vertices.
    extent_margin : :py:class:`float <float>`, default: 0.05
        Fractional margin added to the automatically inferred extent.
    vertex_margin : :py:class:`float <float>`, default: 0.05
        Fractional padding used when selecting dual polygons near the visible
        extent.
    global_plot : :py:class:`bool <bool>`, default: False
        If True, draw a global extent centered on ``center_lon``.
    center_lon : :py:class:`float <float>`, optional
        Longitude center used for wrapping. If None, infer it from the extent or
        valid vertices.
    input_radians : :py:class:`bool <bool>`, default: True
        If True, MPAS longitude and latitude variables are interpreted as
        radians and converted to degrees.
    skip_large_lon_jump : :py:class:`bool <bool>`, default: True
        If True, skip polygons with very large wrapped longitude spans.
    max_polygon_lon_span : :py:class:`float <float>`, default: 180.0
        Maximum accepted wrapped polygon longitude span in degrees.
    cmap : :py:class:`str <str>` or :py:class:`matplotlib.colors.Colormap <matplotlib.colors.Colormap>`, optional
        Colormap used for polygon values.
    vmin, vmax : :py:class:`float <float>`, optional
        Color limits. If omitted, infer them from selected values.
    symmetric : :py:class:`bool <bool>`, optional
        If True, infer symmetric color limits around zero. If None, default to
        True for vertex-dual plots.
    percentile : :py:class:`float <float>`, default: 98
        Percentile used for automatic color-limit inference.
    edgecolor, linewidth
        Polygon edge style passed to :py:class:`matplotlib.collections.PolyCollection`.
    add_colorbar : :py:class:`bool <bool>`, default: True
        If True, add a colorbar for the polygon collection.
    cbar_kwargs : :py:class:`dict <dict>`, optional
        Keyword arguments passed to ``Figure.colorbar``.
    title : :py:class:`str <str>`, optional
        Axes title. If None, use the data long name.
    xlabel, ylabel : :py:class:`str <str>`
        Axis labels used for plain Matplotlib axes.
    aspect : :py:class:`str <str>` or :py:class:`float <float>`, default: "auto"
        Aspect setting used for plain Matplotlib axes.
    cells_on_vertex_name : :py:class:`str <str>`
        MPAS variable giving neighboring cells for each vertex.
    lon_name_cell, lat_name_cell, lon_name_vertex, lat_name_vertex : :py:class:`str <str>`
        MPAS coordinate variable names.
    use_geometry_cache : :py:class:`bool <bool>`, default: True
        If True, cache vertex-dual polygon geometry for repeated calls.

    Returns
    -------
    :py:class:`matplotlib.collections.PolyCollection <matplotlib.collections.PolyCollection>`
        Polygon collection added to the axes.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        vertex_voronoi
    """

    if hasattr(da, "values"):
        z = np.asarray(da.values)
        varname = da.name
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

    cells_on_vertex = ds[cells_on_vertex_name].values.astype(int) - 1

    lon_cell_raw = np.asarray(ds[lon_name_cell].values)
    lat_cell = np.asarray(ds[lat_name_cell].values)

    lon_vertex_raw = np.asarray(ds[lon_name_vertex].values)
    lat_vertex = np.asarray(ds[lat_name_vertex].values)

    if input_radians:
        lon_cell_raw = _rad2deg(lon_cell_raw)
        lat_cell = _rad2deg(lat_cell)
        lon_vertex_raw = _rad2deg(lon_vertex_raw)
        lat_vertex = _rad2deg(lat_vertex)

    n_cells = lon_cell_raw.size
    n_vertices = lon_vertex_raw.size

    if z.size == n_cells:
        z_cell = z
        valid_cids = cells_on_vertex >= 0
        safe_cids = np.where(valid_cids, cells_on_vertex, 0)
        vals = z_cell[safe_cids]
        finite_vals = valid_cids & np.isfinite(vals)

        sums = np.where(finite_vals, vals, 0.0).sum(axis=1)
        counts = finite_vals.sum(axis=1)

        z = np.full(n_vertices, np.nan, dtype=float)
        np.divide(sums, counts, out=z, where=counts > 0)
    elif z.size != n_vertices:
        raise ValueError(
            f"`da` length does not match nVertices or nCells: "
            f"len(da)={z.size}, nVertices={n_vertices}, nCells={n_cells}"
        )

    valid = np.isfinite(z)

    if center_lon is None:
        center_lon = _infer_center_lon(
            lon_min=lon_min,
            lon_max=lon_max,
            lon_cell=lon_vertex_raw,
            valid=valid,
        )

    lon_cell = _lon_wrap_to_center(lon_cell_raw, center_lon)
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

    # ---- axis ----
    if ax is None:
        ax = plt.gca()

    fig = ax.figure

    cartopy_axis = _is_cartopy_axis(ax)

    if transform is None and cartopy_axis:
        transform = ccrs.PlateCarree()

    geometry_cache_key = (
        "vertex_dual",
        id(ds),
        cells_on_vertex_name,
        lon_name_cell,
        lat_name_cell,
        center_lon,
        input_radians,
    )
    geometry = _get_vertex_dual_geometry(
        cells_on_vertex,
        lon_cell,
        lat_cell,
        cache_key=geometry_cache_key,
        use_cache=use_geometry_cache,
    )

    safe_cids = geometry["safe_cids"]
    n_valid = geometry["n_valid"]
    candidate = geometry["valid"] & np.isfinite(z)

    if skip_large_lon_jump:
        candidate &= (geometry["lon_hi"] - geometry["lon_lo"]) <= max_polygon_lon_span

    if lon_min_plot is not None and not global_plot:
        dlon_plot = lon_max_plot - lon_min_plot
        dlat_plot = lat_max_plot - lat_min_plot
        lon_pad = vertex_margin * dlon_plot
        lat_pad = vertex_margin * dlat_plot

        # Select by polygon overlap with a padded window instead of only
        # requiring the vertex point to be inside the visible extent.  This
        # mirrors `plot_cell_voronoi` and keeps partial dual cells along
        # plot boundaries, avoiding blank strips near the edges.
        candidate &= ~(
            (geometry["lon_hi"] < lon_min_plot - lon_pad)
            | (geometry["lon_lo"] > lon_max_plot + lon_pad)
            | (geometry["lat_hi"] < lat_min_plot - lat_pad)
            | (geometry["lat_lo"] > lat_max_plot + lat_pad)
        )

    polys = []
    values = []

    for v in np.flatnonzero(candidate):
        val = z[v]
        cids = safe_cids[v, : n_valid[v]]
        lons = lon_cell[cids]
        lats = lat_cell[cids]

        polys.append(np.column_stack([lons, lats]))
        values.append(val)

    values = np.asarray(values)

    if len(polys) == 0:
        raise RuntimeError(
            "No vertex dual triangles selected. Check extent or cellsOnVertex."
        )

    cmap = _infer_default_cmap(values, cmap=cmap)

    if symmetric is None:
        symmetric = True

    if vmin is None or vmax is None:
        if symmetric:
            vmax_auto = np.nanpercentile(np.abs(values), percentile)
            vmin_auto = -vmax_auto
        else:
            vmin_auto = np.nanpercentile(values, 100 - percentile)
            vmax_auto = np.nanpercentile(values, percentile)

        if vmin is None:
            vmin = vmin_auto
        if vmax is None:
            vmax = vmax_auto

    pc_kwargs = dict(
        cmap=cmap,
        edgecolor=edgecolor,
        linewidth=linewidth,
    )

    if cartopy_axis and transform is not None:
        pc_kwargs["transform"] = transform

    pc = PolyCollection(polys, closed=True, **pc_kwargs)

    pc.set_array(values)
    pc.set_clim(vmin, vmax)

    ax.add_collection(pc)

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
        else:
            ax.autoscale_view()

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        if aspect is not None:
            ax.set_aspect(aspect, adjustable="box")

    if title is None:
        title = str(long_name)

    ax.set_title(title)

    if add_colorbar:
        if cbar_kwargs is None:
            cbar_kwargs = {}

        cb = fig.colorbar(pc, ax=ax, **cbar_kwargs)
        if units:
            cb.set_label(f"{long_name} ({units})")
        else:
            cb.set_label(str(long_name))

    return pc
