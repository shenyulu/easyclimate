"""
Contour plots for MPAS cell data.
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from .voronoi_extract import extract_cell_latlon
from .common import _get_cached_triangulation, _infer_default_cmap, _is_cartopy_axis

__all__ = [
    "plot_cell_contourf",
    "plot_cell_contour",
]


def plot_cell_contourf(
    ds,
    da,
    *,
    ax=None,
    # extent
    lon_min=None,
    lon_max=None,
    lat_min=None,
    lat_max=None,
    auto_extent=True,
    extent_margin=0.05,
    contour_margin=0.05,
    global_plot=False,
    center_lon=None,
    # contour
    levels=10,
    cmap=None,
    vmin=None,
    vmax=None,
    symmetric=False,
    percentile=98,
    extend=None,
    # triangle quality control
    max_triangle_edge=None,
    # coordinate
    input_radians=True,
    # cartopy
    transform=None,
    # decoration
    add_colorbar=True,
    cbar_kwargs=None,
    title=None,
    # MPAS names
    lon_name_cell="lonCell",
    lat_name_cell="latCell",
    lon_name_vertex="lonVertex",
    lat_name_vertex="latVertex",
    use_triangulation_cache=True,
):
    """
    Draw filled contour for MPAS cell-centered scalar field.

    The scalar values are drawn from MPAS cell centers using
    :py:class:`matplotlib.tri.Triangulation`.

    Parameters
    ----------
    ds : :py:class:`xarray.Dataset <xarray.Dataset>`
        MPAS dataset containing cell and vertex coordinate variables.
    da : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`numpy.ndarray <numpy.ndarray>`
        Cell-centered scalar values. The data must be one-dimensional on
        ``nCells`` after any time or vertical-level selection.
    ax : :py:class:`matplotlib.axes.Axes <matplotlib.axes.Axes>`, optional
        Axes on which to draw. By default, use the current axes.
    lon_min, lon_max, lat_min, lat_max : :py:class:`float <float>`, optional
        Plot extent in degrees. Cross-dateline longitude windows are supported.
    auto_extent : :py:class:`bool <bool>`, default: True
        If True and no explicit extent is supplied, infer the extent from valid
        cell centers.
    extent_margin : :py:class:`float <float>`, default: 0.05
        Fractional margin added to the automatically inferred extent.
    contour_margin : :py:class:`float <float>`, default: 0.05
        Fractional padding used when selecting source cells for contouring.
    global_plot : :py:class:`bool <bool>`, default: False
        If True, draw a global extent centered on ``center_lon``.
    center_lon : :py:class:`float <float>`, optional
        Longitude center used for wrapping. If None, infer it from the extent or
        valid cell centers.
    levels : :py:class:`int <int>` or array-like, default: 10
        Number of contour levels, or explicit level values.
    cmap : :py:class:`str <str>` or :py:class:`matplotlib.colors.Colormap <matplotlib.colors.Colormap>`, optional
        Colormap used for filled contours.
    vmin, vmax : :py:class:`float <float>`, optional
        Color/data limits used by contour levels and color mapping. If one side
        is omitted, infer it from selected finite values.
    symmetric : :py:class:`bool <bool>`, default: False
        If True, infer symmetric limits around zero when ``vmin`` or ``vmax`` is
        omitted.
    percentile : :py:class:`float <float>`, default: 98
        Percentile used for automatic color-limit inference.
    extend : {"neither", "both", "min", "max"}, optional
        Colorbar extension mode. If None, infer it from data and levels.
    max_triangle_edge : :py:class:`float <float>`, optional
        Maximum accepted edge length for triangles. Longer triangles are masked.
    input_radians : :py:class:`bool <bool>`, default: True
        If True, MPAS longitude and latitude variables are interpreted as
        radians and converted to degrees.
    transform : :py:class:`cartopy.crs.CRS <cartopy.crs.CRS>`, optional
        Coordinate reference system of the input coordinates for Cartopy axes.
    add_colorbar : :py:class:`bool <bool>`, default: True
        If True, add a colorbar for the filled contour set.
    cbar_kwargs : :py:class:`dict <dict>`, optional
        Keyword arguments passed to ``Figure.colorbar``.
    title : :py:class:`str <str>`, optional
        Axes title. If None, use the data long name.
    lon_name_cell, lat_name_cell, lon_name_vertex, lat_name_vertex : :py:class:`str <str>`
        MPAS coordinate variable names.
    use_triangulation_cache : :py:class:`bool <bool>`, default: True
        If True, cache the triangulation for repeated calls on the same mesh and
        extent.

    Returns
    -------
    :py:class:`matplotlib.contour.QuadContourSet <matplotlib.contour.QuadContourSet>`
        Filled contour set returned by ``Axes.tricontourf``.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        cell_contour
    """

    # ---- data ----
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
            f"`da` must be 1D on nCells after selection, got shape {z.shape}."
        )

    # ---- lon/lat extraction, including dateline handling ----
    geo = extract_cell_latlon(
        ds,
        da=da,
        lon_min=lon_min,
        lon_max=lon_max,
        lat_min=lat_min,
        lat_max=lat_max,
        auto_extent=auto_extent,
        extent_margin=extent_margin,
        global_plot=global_plot,
        center_lon=center_lon,
        input_radians=input_radians,
        lon_name_cell=lon_name_cell,
        lat_name_cell=lat_name_cell,
        lon_name_vertex=lon_name_vertex,
        lat_name_vertex=lat_name_vertex,
    )

    lon_cell = geo["lon_cell"]
    lat_cell = geo["lat_cell"]

    lon_min_plot = geo["lon_min_plot"]
    lon_max_plot = geo["lon_max_plot"]
    lat_min_plot = geo["lat_min_plot"]
    lat_max_plot = geo["lat_max_plot"]
    center_lon = geo["center_lon"]

    # ---- subset ----
    mask = np.isfinite(lon_cell) & np.isfinite(lat_cell) & np.isfinite(z)

    if lon_min_plot is not None and not global_plot:
        dlon_plot = lon_max_plot - lon_min_plot
        dlat_plot = lat_max_plot - lat_min_plot
        lon_pad = contour_margin * dlon_plot
        lat_pad = contour_margin * dlat_plot

        mask &= (
            (lon_cell >= lon_min_plot - lon_pad)
            & (lon_cell <= lon_max_plot + lon_pad)
            & (lat_cell >= lat_min_plot - lat_pad)
            & (lat_cell <= lat_max_plot + lat_pad)
        )

    idx = np.flatnonzero(mask)
    x = lon_cell[idx]
    y = lat_cell[idx]
    zz = z[idx]

    if x.size < 3:
        raise RuntimeError("Not enough valid cell centers to draw contours.")

    cmap = _infer_default_cmap(zz, cmap=cmap)

    # ---- color/data limits and contour levels ----
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

    # ---- axis ----
    if ax is None:
        ax = plt.gca()

    fig = ax.figure

    cartopy_axis = _is_cartopy_axis(ax)

    if transform is None and cartopy_axis:
        transform = ccrs.PlateCarree()

    # ---- triangulation ----
    triang_cache_key = (
        "cell_contour_triangulation",
        id(ds),
        lon_name_cell,
        lat_name_cell,
        center_lon,
        input_radians,
        lon_min_plot,
        lon_max_plot,
        lat_min_plot,
        lat_max_plot,
        contour_margin,
        global_plot,
        max_triangle_edge,
        tuple(idx.tolist()),
    )
    triang = _get_cached_triangulation(
        x,
        y,
        max_edge=max_triangle_edge,
        cache_key=triang_cache_key,
        use_cache=use_triangulation_cache,
    )

    # ---- draw filled contour ----
    kwargs = {}

    if cartopy_axis and transform is not None:
        kwargs["transform"] = transform

    cs = ax.tricontourf(
        triang,
        zz,
        levels=contour_levels,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        extend=extend,
        **kwargs,
    )

    # ---- extent and decoration ----
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

        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_aspect("auto")

    if add_colorbar:
        if cbar_kwargs is None:
            cbar_kwargs = {}

        cb = fig.colorbar(cs, ax=ax, **cbar_kwargs)

        if units:
            cb.set_label(f"{long_name} ({units})")
        else:
            cb.set_label(str(long_name))

    if title is None:
        title = str(long_name)

    ax.set_title(title)

    return cs


def plot_cell_contour(
    ds,
    da,
    *,
    ax=None,
    # extent
    lon_min=None,
    lon_max=None,
    lat_min=None,
    lat_max=None,
    auto_extent=True,
    extent_margin=0.05,
    contour_margin=0.05,
    global_plot=False,
    center_lon=None,
    # contour
    levels=10,
    colors=None,
    linewidths=0.8,
    linestyles=None,
    cmap=None,
    vmin=None,
    vmax=None,
    symmetric=False,
    percentile=98,
    extend=None,
    # triangle quality control
    max_triangle_edge=None,
    # coordinate
    input_radians=True,
    # cartopy
    transform=None,
    # decoration
    title=None,
    # MPAS names
    lon_name_cell="lonCell",
    lat_name_cell="latCell",
    lon_name_vertex="lonVertex",
    lat_name_vertex="latVertex",
    use_triangulation_cache=True,
):
    """
    Draw contour lines for MPAS cell-centered scalar field.

    Labels are not added automatically. Add them outside this function with
    ``ax.clabel(cs, inline=True, fontsize=8, fmt="%g")``.

    Parameters
    ----------
    ds : :py:class:`xarray.Dataset <xarray.Dataset>`
        MPAS dataset containing cell and vertex coordinate variables.
    da : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`numpy.ndarray <numpy.ndarray>`
        Cell-centered scalar values. The data must be one-dimensional on
        ``nCells`` after any time or vertical-level selection.
    ax : :py:class:`matplotlib.axes.Axes <matplotlib.axes.Axes>`, optional
        Axes on which to draw. By default, use the current axes.
    lon_min, lon_max, lat_min, lat_max : :py:class:`float <float>`, optional
        Plot extent in degrees. Cross-dateline longitude windows are supported.
    auto_extent : :py:class:`bool <bool>`, default: True
        If True and no explicit extent is supplied, infer the extent from valid
        cell centers.
    extent_margin : :py:class:`float <float>`, default: 0.05
        Fractional margin added to the automatically inferred extent.
    contour_margin : :py:class:`float <float>`, default: 0.05
        Fractional padding used when selecting source cells for contouring.
    global_plot : :py:class:`bool <bool>`, default: False
        If True, draw a global extent centered on ``center_lon``.
    center_lon : :py:class:`float <float>`, optional
        Longitude center used for wrapping. If None, infer it from the extent or
        valid cell centers.
    levels : :py:class:`int <int>` or array-like, default: 10
        Number of contour levels, or explicit level values.
    colors : color or :py:class:`collections.abc.Sequence <collections.abc.Sequence>` of colors, optional
        Fixed line colors passed to ``Axes.tricontour``. If omitted, use
        ``cmap``.
    linewidths, linestyles
        Contour line style options passed to ``Axes.tricontour``.
    cmap : :py:class:`str <str>` or :py:class:`matplotlib.colors.Colormap <matplotlib.colors.Colormap>`, optional
        Colormap used when ``colors`` is not supplied.
    vmin, vmax : :py:class:`float <float>`, optional
        Data limits used by contour levels and color mapping.
    symmetric : :py:class:`bool <bool>`, default: False
        If True, infer symmetric limits around zero when ``vmin`` or ``vmax`` is
        omitted.
    percentile : :py:class:`float <float>`, default: 98
        Percentile used for automatic data-limit inference.
    extend : {"neither", "both", "min", "max"}, optional
        Contour extension mode. If None, infer it from data and levels.
    max_triangle_edge : :py:class:`float <float>`, optional
        Maximum accepted edge length for triangles. Longer triangles are masked.
    input_radians : :py:class:`bool <bool>`, default: True
        If True, MPAS longitude and latitude variables are interpreted as
        radians and converted to degrees.
    transform : :py:class:`cartopy.crs.CRS <cartopy.crs.CRS>`, optional
        Coordinate reference system of the input coordinates for Cartopy axes.
    title : :py:class:`str <str>`, optional
        Axes title. If None, use the data long name.
    lon_name_cell, lat_name_cell, lon_name_vertex, lat_name_vertex : :py:class:`str <str>`
        MPAS coordinate variable names.
    use_triangulation_cache : :py:class:`bool <bool>`, default: True
        If True, cache the triangulation for repeated calls on the same mesh and
        extent.

    Returns
    -------
    :py:class:`matplotlib.contour.QuadContourSet <matplotlib.contour.QuadContourSet>`
        Contour set returned by ``Axes.tricontour``.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        cell_contour
    """

    # ---- data ----
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
            f"`da` must be 1D on nCells after selection, got shape {z.shape}."
        )

    # ---- lon/lat extraction, including dateline handling ----
    geo = extract_cell_latlon(
        ds,
        da=da,
        lon_min=lon_min,
        lon_max=lon_max,
        lat_min=lat_min,
        lat_max=lat_max,
        auto_extent=auto_extent,
        extent_margin=extent_margin,
        global_plot=global_plot,
        center_lon=center_lon,
        input_radians=input_radians,
        lon_name_cell=lon_name_cell,
        lat_name_cell=lat_name_cell,
        lon_name_vertex=lon_name_vertex,
        lat_name_vertex=lat_name_vertex,
    )

    lon_cell = geo["lon_cell"]
    lat_cell = geo["lat_cell"]

    lon_min_plot = geo["lon_min_plot"]
    lon_max_plot = geo["lon_max_plot"]
    lat_min_plot = geo["lat_min_plot"]
    lat_max_plot = geo["lat_max_plot"]

    # ---- subset ----
    mask = np.isfinite(lon_cell) & np.isfinite(lat_cell) & np.isfinite(z)

    if lon_min_plot is not None and not global_plot:
        dlon_plot = lon_max_plot - lon_min_plot
        dlat_plot = lat_max_plot - lat_min_plot
        lon_pad = contour_margin * dlon_plot
        lat_pad = contour_margin * dlat_plot

        mask &= (
            (lon_cell >= lon_min_plot - lon_pad)
            & (lon_cell <= lon_max_plot + lon_pad)
            & (lat_cell >= lat_min_plot - lat_pad)
            & (lat_cell <= lat_max_plot + lat_pad)
        )

    idx = np.flatnonzero(mask)
    x = lon_cell[idx]
    y = lat_cell[idx]
    zz = z[idx]

    if x.size < 3:
        raise RuntimeError("Not enough valid cell centers to draw contours.")

    cmap = _infer_default_cmap(zz, cmap=cmap)

    # ---- data limits and contour levels ----
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

    # ---- axis ----
    if ax is None:
        ax = plt.gca()

    cartopy_axis = _is_cartopy_axis(ax)

    if transform is None and cartopy_axis:
        transform = ccrs.PlateCarree()

    # ---- triangulation ----
    triang_cache_key = (
        "cell_contour_triangulation",
        id(ds),
        lon_name_cell,
        lat_name_cell,
        geo["center_lon"],
        input_radians,
        lon_min_plot,
        lon_max_plot,
        lat_min_plot,
        lat_max_plot,
        contour_margin,
        global_plot,
        max_triangle_edge,
        tuple(idx.tolist()),
    )
    triang = _get_cached_triangulation(
        x,
        y,
        max_edge=max_triangle_edge,
        cache_key=triang_cache_key,
        use_cache=use_triangulation_cache,
    )

    # ---- draw contour lines ----
    kwargs = {}

    if cartopy_axis and transform is not None:
        kwargs["transform"] = transform

    contour_kwargs = dict(
        linewidths=linewidths,
        linestyles=linestyles,
        vmin=vmin,
        vmax=vmax,
        extend=extend,
    )

    if colors is None:
        contour_kwargs["cmap"] = cmap
    else:
        contour_kwargs["colors"] = colors

    cs = ax.tricontour(
        triang,
        zz,
        levels=contour_levels,
        **contour_kwargs,
        **kwargs,
    )

    # ---- extent and decoration ----
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

        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_aspect("auto")

    if title is None:
        title = str(long_name)

    ax.set_title(title)

    return cs
