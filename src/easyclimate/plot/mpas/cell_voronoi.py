"""
Voronoi-style plots for MPAS cells.
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
    _get_cell_polygon_geometry,
)


__all__ = ["plot_cell_voronoi"]


def plot_cell_voronoi(
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
    cell_margin=0.05,
    global_plot=False,
    # coordinate control
    input_radians=True,
    center_lon=None,
    skip_large_lon_jump=True,
    max_polygon_lon_span=180.0,
    # color control
    cmap=None,
    vmin=None,
    vmax=None,
    symmetric=False,
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
    lon_name_cell="lonCell",
    lat_name_cell="latCell",
    lon_name_vertex="lonVertex",
    lat_name_vertex="latVertex",
    vertices_on_cell_name="verticesOnCell",
    n_edges_on_cell_name="nEdgesOnCell",
    use_geometry_cache=True,
):
    """
    Plot MPAS cell-centered scalar field on native polygon cells.

    The scalar field is rendered on each native MPAS cell polygon. The function
    works with both plain Matplotlib axes and Cartopy GeoAxes.

    Parameters
    ----------
    ds : :py:class:`xarray.Dataset <xarray.Dataset>`
        Dataset containing MPAS mesh geometry.
    da : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`numpy.ndarray <numpy.ndarray>`
        Cell-centered scalar field. Must be 1D on nCells after selection.
    ax : :py:class:`matplotlib.axes.Axes <matplotlib.axes.Axes>`, optional
        Axes on which to draw. By default, use the current axes.
    projection : :py:class:`cartopy.crs.CRS <cartopy.crs.CRS>`, optional
        Reserved for compatibility with older call sites.
    transform : :py:class:`cartopy.crs.CRS <cartopy.crs.CRS>`, optional
        Coordinate reference system of the input coordinates for Cartopy axes.
    lon_min, lon_max, lat_min, lat_max : :py:class:`float <float>`, optional
        Plot extent. Cross-dateline windows are supported, for example:
        ``lon_min=-190, lon_max=-60`` describes the continuous longitude window
        from 170E to 60W.
    auto_extent : :py:class:`bool <bool>`, default: True
        If True and no explicit extent is supplied, infer the extent from valid
        cell centers.
    extent_margin : :py:class:`float <float>`, default: 0.05
        Fractional margin added to the automatically inferred extent.
    cell_margin : :py:class:`float <float>`, default: 0.05
        Extra fraction of the plot extent used to select cells outside the
        visible window. This helps avoid missing partial cells along boundaries.
    global_plot : :py:class:`bool <bool>`, default: False
        If True, draw a global extent centered on ``center_lon``.
    input_radians : :py:class:`bool <bool>`, default: True
        If True, MPAS longitude and latitude variables are interpreted as
        radians and converted to degrees.
    center_lon : :py:class:`float <float>`, optional
        Longitude wrapping center. If None, inferred from lon_min/lon_max
        or valid cell centers.
    skip_large_lon_jump : :py:class:`bool <bool>`, default: True
        If True, skip polygons with very large wrapped longitude spans.
    max_polygon_lon_span : :py:class:`float <float>`, default: 180.0
        Maximum accepted wrapped polygon longitude span in degrees.
    cmap : :py:class:`str <str>` or :py:class:`matplotlib.colors.Colormap <matplotlib.colors.Colormap>`, optional
        Colormap used for cell values.
    vmin, vmax : :py:class:`float <float>`, optional
        Color limits. If omitted, infer them from selected values.
    symmetric : :py:class:`bool <bool>`, default: False
        If True, infer symmetric color limits around zero.
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
    figsize : (:py:class:`float <float>`, :py:class:`float <float>`), default: (8, 6)
        Reserved for compatibility with older call sites.
    xlabel, ylabel : :py:class:`str <str>`
        Axis labels used for plain Matplotlib axes.
    aspect : :py:class:`str <str>` or :py:class:`float <float>`, default: "auto"
        Aspect setting used for plain Matplotlib axes.
    lon_name_cell, lat_name_cell, lon_name_vertex, lat_name_vertex : :py:class:`str <str>`
        MPAS coordinate variable names.
    vertices_on_cell_name, n_edges_on_cell_name : :py:class:`str <str>`
        MPAS connectivity variable names.
    use_geometry_cache : :py:class:`bool <bool>`, default: True
        If True, cache native cell polygon geometry for repeated calls.

    Returns
    -------
    ax : :py:class:`matplotlib.axes.Axes <matplotlib.axes.Axes>`
        Axes containing the plot.
    pc : :py:class:`matplotlib.collections.PolyCollection <matplotlib.collections.PolyCollection>`
        Polygon collection added to the axes.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        cell_voronoi
    """

    # ------------------------------------------------------------
    # Read data first
    # ------------------------------------------------------------
    if hasattr(da, "values"):
        values_all = np.asarray(da.values)

        varname = getattr(da, "name", None)
        units = da.attrs.get("units", "")
        long_name = da.attrs.get("long_name", varname or "value")
    else:
        values_all = np.asarray(da)
        varname = "value"
        units = ""
        long_name = "value"

    if values_all.ndim != 1:
        raise ValueError(
            f"`da` must be 1D on nCells after selection, got shape {values_all.shape}. "
            "Example: ds['divergence'].isel(Time=0, nVertLevels=0)"
        )

    valid_data = np.isfinite(values_all)

    # ------------------------------------------------------------
    # Read mesh geometry
    # ------------------------------------------------------------
    vertices_on_cell = ds[vertices_on_cell_name].values.astype(int) - 1
    n_edges_on_cell = ds[n_edges_on_cell_name].values.astype(int)

    lon_vertex_raw = np.asarray(ds[lon_name_vertex].values)
    lat_vertex = np.asarray(ds[lat_name_vertex].values)

    lon_cell_raw = np.asarray(ds[lon_name_cell].values)
    lat_cell = np.asarray(ds[lat_name_cell].values)

    if input_radians:
        lon_vertex_raw = _rad2deg(lon_vertex_raw)
        lat_vertex = _rad2deg(lat_vertex)
        lon_cell_raw = _rad2deg(lon_cell_raw)
        lat_cell = _rad2deg(lat_cell)

    n_cells = lon_cell_raw.size

    if values_all.size != n_cells:
        raise ValueError(
            f"`da` length does not match nCells: "
            f"len(da)={values_all.size}, nCells={n_cells}"
        )

    # ------------------------------------------------------------
    # Infer longitude wrapping center
    # ------------------------------------------------------------
    if center_lon is None:
        center_lon = _infer_center_lon(
            lon_min=lon_min,
            lon_max=lon_max,
            lon_cell=lon_cell_raw,
            valid=valid_data,
        )

    lon_vertex = _lon_wrap_to_center(lon_vertex_raw, center_lon)
    lon_cell = _lon_wrap_to_center(lon_cell_raw, center_lon)

    # ------------------------------------------------------------
    # 4. Decide plotting extent
    # ------------------------------------------------------------
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
                lon_cell[valid_data],
                lat_cell[valid_data],
                margin=extent_margin,
            )
        )
    else:
        lon_min_plot = lon_max_plot = None
        lat_min_plot = lat_max_plot = None

    # ------------------------------------------------------------
    # Create axis
    # ------------------------------------------------------------
    if ax is None:
        ax = plt.gca()

    fig = ax.figure

    cartopy_axis = _is_cartopy_axis(ax)

    if transform is None and cartopy_axis:
        transform = ccrs.PlateCarree()

    # ------------------------------------------------------------
    # Build cell polygons
    # ------------------------------------------------------------
    geometry_cache_key = (
        "cell_polygon",
        id(ds),
        vertices_on_cell_name,
        n_edges_on_cell_name,
        lon_name_vertex,
        lat_name_vertex,
        center_lon,
        input_radians,
    )
    geometry = _get_cell_polygon_geometry(
        vertices_on_cell,
        n_edges_on_cell,
        lon_vertex,
        lat_vertex,
        cache_key=geometry_cache_key,
        use_cache=use_geometry_cache,
    )

    candidate = geometry["valid"] & valid_data

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

    if skip_large_lon_jump:
        candidate &= (geometry["lon_hi"] - geometry["lon_lo"]) <= max_polygon_lon_span

    safe_vids = geometry["safe_vids"]
    n_edges = geometry["n_edges"]

    polys = []
    values = []

    for c in np.flatnonzero(candidate):
        val = values_all[c]
        vids = safe_vids[c, : n_edges[c]]
        lons = lon_vertex[vids]
        lats = lat_vertex[vids]
        polys.append(np.column_stack([lons, lats]))
        values.append(val)

    values = np.asarray(values)

    if len(polys) == 0:
        raise RuntimeError(
            "No cells selected. Check extent, longitude wrapping, or data validity."
        )

    cmap = _infer_default_cmap(values, cmap=cmap)

    # ------------------------------------------------------------
    # Color limits
    # ------------------------------------------------------------
    if vmin is None or vmax is None:
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

    # ------------------------------------------------------------
    # Draw PolyCollection
    # ------------------------------------------------------------
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

    # ------------------------------------------------------------
    # Extent and decorations
    # ------------------------------------------------------------
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

    # ------------------------------------------------------------
    # Colorbar and title
    # ------------------------------------------------------------
    if add_colorbar:
        if cbar_kwargs is None:
            cbar_kwargs = {}

        cb = fig.colorbar(
            pc,
            ax=ax,
            **cbar_kwargs,
        )

        if units:
            cb.set_label(f"{long_name} ({units})")
        else:
            cb.set_label(str(long_name))

    if title is None:
        title = str(long_name)

    ax.set_title(title)

    return ax, pc
