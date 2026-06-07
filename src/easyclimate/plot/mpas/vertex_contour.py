"""
Contour plots for MPAS vertex data.
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from .common import (
    _decorate_vertex_contour_axis,
    _get_cached_triangulation,
    _infer_default_cmap,
    _is_cartopy_axis,
    _prepare_vertex_contour_geometry,
    _resolve_contour_levels,
)


__all__ = [
    "plot_vertex_contourf",
    "plot_vertex_contour",
]


def plot_vertex_contourf(
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
    # triangle quality/cache control
    max_triangle_edge=None,
    use_triangulation_cache=True,
    # coordinate
    input_radians=True,
    # cartopy
    transform=None,
    # decoration
    add_colorbar=True,
    cbar_kwargs=None,
    title=None,
    xlabel="Longitude",
    ylabel="Latitude",
    aspect="auto",
    # MPAS names
    lon_name_vertex="lonVertex",
    lat_name_vertex="latVertex",
):
    """
    Draw filled contour for MPAS vertex-centered scalar field.

    Parameters
    ----------
    ds : :py:class:`xarray.Dataset <xarray.Dataset>`
        MPAS dataset containing vertex coordinate variables.
    da : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`numpy.ndarray <numpy.ndarray>`
        Vertex-centered scalar values. The data must be one-dimensional on
        ``nVertices`` after any time or vertical-level selection.
    ax : :py:class:`matplotlib.axes.Axes <matplotlib.axes.Axes>`, optional
        Axes on which to draw. By default, use the current axes.
    lon_min, lon_max, lat_min, lat_max : :py:class:`float <float>`, optional
        Plot extent in degrees. Cross-dateline longitude windows are supported.
    auto_extent : :py:class:`bool <bool>`, default: True
        If True and no explicit extent is supplied, infer the extent from valid
        vertices.
    extent_margin : :py:class:`float <float>`, default: 0.05
        Fractional margin added to the automatically inferred extent.
    contour_margin : :py:class:`float <float>`, default: 0.05
        Fractional padding used when selecting source vertices for contouring.
    global_plot : :py:class:`bool <bool>`, default: False
        If True, draw a global extent centered on ``center_lon``.
    center_lon : :py:class:`float <float>`, optional
        Longitude center used for wrapping. If None, infer it from the extent or
        valid vertices.
    levels : :py:class:`int <int>` or array-like, default: 10
        Number of contour levels, or explicit level values.
    cmap : :py:class:`str <str>` or :py:class:`matplotlib.colors.Colormap <matplotlib.colors.Colormap>`, optional
        Colormap used for filled contours.
    vmin, vmax : :py:class:`float <float>`, optional
        Data limits used by contour levels and color mapping.
    symmetric : :py:class:`bool <bool>`, default: False
        If True, infer symmetric limits around zero when ``vmin`` or ``vmax`` is
        omitted.
    percentile : :py:class:`float <float>`, default: 98
        Percentile used for automatic data-limit inference.
    extend : {"neither", "both", "min", "max"}, optional
        Colorbar extension mode. If None, infer it from data and levels.
    max_triangle_edge : :py:class:`float <float>`, optional
        Maximum accepted edge length for triangles. Longer triangles are masked.
    use_triangulation_cache : :py:class:`bool <bool>`, default: True
        If True, cache the triangulation for repeated calls on the same mesh and
        extent.
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
    xlabel, ylabel : :py:class:`str <str>`
        Axis labels used for plain Matplotlib axes.
    aspect : :py:class:`str <str>` or :py:class:`float <float>`, default: "auto"
        Aspect setting used for plain Matplotlib axes.
    lon_name_vertex, lat_name_vertex : :py:class:`str <str>`
        MPAS vertex coordinate variable names.

    Returns
    -------
    :py:class:`matplotlib.contour.QuadContourSet <matplotlib.contour.QuadContourSet>`
        Filled contour set returned by ``Axes.tricontourf``.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        vertex_contour
    """
    geo = _prepare_vertex_contour_geometry(
        ds,
        da,
        lon_min=lon_min,
        lon_max=lon_max,
        lat_min=lat_min,
        lat_max=lat_max,
        auto_extent=auto_extent,
        extent_margin=extent_margin,
        contour_margin=contour_margin,
        global_plot=global_plot,
        center_lon=center_lon,
        input_radians=input_radians,
        lon_name_vertex=lon_name_vertex,
        lat_name_vertex=lat_name_vertex,
    )

    x = geo["x"]
    y = geo["y"]
    zz = geo["zz"]
    idx = geo["idx"]

    cmap = _infer_default_cmap(zz, cmap=cmap)
    contour_levels, vmin, vmax, extend = _resolve_contour_levels(
        zz,
        levels,
        vmin=vmin,
        vmax=vmax,
        symmetric=symmetric,
        percentile=percentile,
        extend=extend,
    )

    if ax is None:
        ax = plt.gca()

    fig = ax.figure
    cartopy_axis = _is_cartopy_axis(ax)

    if transform is None and cartopy_axis:
        transform = ccrs.PlateCarree()

    triang_cache_key = (
        "vertex_contour_triangulation",
        id(ds),
        lon_name_vertex,
        lat_name_vertex,
        geo["center_lon"],
        input_radians,
        geo["lon_min_plot"],
        geo["lon_max_plot"],
        geo["lat_min_plot"],
        geo["lat_max_plot"],
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

    _decorate_vertex_contour_axis(
        ax,
        cartopy_axis=cartopy_axis,
        global_plot=global_plot,
        lon_min_plot=geo["lon_min_plot"],
        lon_max_plot=geo["lon_max_plot"],
        lat_min_plot=geo["lat_min_plot"],
        lat_max_plot=geo["lat_max_plot"],
        xlabel=xlabel,
        ylabel=ylabel,
        aspect=aspect,
    )

    if add_colorbar:
        if cbar_kwargs is None:
            cbar_kwargs = {}

        cb = fig.colorbar(cs, ax=ax, **cbar_kwargs)
        if geo["units"]:
            cb.set_label(f"{geo['long_name']} ({geo['units']})")
        else:
            cb.set_label(str(geo["long_name"]))

    if title is None:
        title = str(geo["long_name"])

    ax.set_title(title)

    return cs


def plot_vertex_contour(
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
    # triangle quality/cache control
    max_triangle_edge=None,
    use_triangulation_cache=True,
    # coordinate
    input_radians=True,
    # cartopy
    transform=None,
    # decoration
    title=None,
    xlabel="Longitude",
    ylabel="Latitude",
    aspect="auto",
    # MPAS names
    lon_name_vertex="lonVertex",
    lat_name_vertex="latVertex",
):
    """
    Draw contour lines for MPAS vertex-centered scalar field.

    Parameters
    ----------
    ds : :py:class:`xarray.Dataset <xarray.Dataset>`
        MPAS dataset containing vertex coordinate variables.
    da : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`numpy.ndarray <numpy.ndarray>`
        Vertex-centered scalar values. The data must be one-dimensional on
        ``nVertices`` after any time or vertical-level selection.
    ax : :py:class:`matplotlib.axes.Axes <matplotlib.axes.Axes>`, optional
        Axes on which to draw. By default, use the current axes.
    lon_min, lon_max, lat_min, lat_max : :py:class:`float <float>`, optional
        Plot extent in degrees. Cross-dateline longitude windows are supported.
    auto_extent : :py:class:`bool <bool>`, default: True
        If True and no explicit extent is supplied, infer the extent from valid
        vertices.
    extent_margin : :py:class:`float <float>`, default: 0.05
        Fractional margin added to the automatically inferred extent.
    contour_margin : :py:class:`float <float>`, default: 0.05
        Fractional padding used when selecting source vertices for contouring.
    global_plot : :py:class:`bool <bool>`, default: False
        If True, draw a global extent centered on ``center_lon``.
    center_lon : :py:class:`float <float>`, optional
        Longitude center used for wrapping. If None, infer it from the extent or
        valid vertices.
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
    use_triangulation_cache : :py:class:`bool <bool>`, default: True
        If True, cache the triangulation for repeated calls on the same mesh and
        extent.
    input_radians : :py:class:`bool <bool>`, default: True
        If True, MPAS longitude and latitude variables are interpreted as
        radians and converted to degrees.
    transform : :py:class:`cartopy.crs.CRS <cartopy.crs.CRS>`, optional
        Coordinate reference system of the input coordinates for Cartopy axes.
    title : :py:class:`str <str>`, optional
        Axes title. If None, use the data long name.
    xlabel, ylabel : :py:class:`str <str>`
        Axis labels used for plain Matplotlib axes.
    aspect : :py:class:`str <str>` or :py:class:`float <float>`, default: "auto"
        Aspect setting used for plain Matplotlib axes.
    lon_name_vertex, lat_name_vertex : :py:class:`str <str>`
        MPAS vertex coordinate variable names.

    Returns
    -------
    :py:class:`matplotlib.contour.QuadContourSet <matplotlib.contour.QuadContourSet>`
        Contour set returned by ``Axes.tricontour``.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        vertex_contour
    """
    geo = _prepare_vertex_contour_geometry(
        ds,
        da,
        lon_min=lon_min,
        lon_max=lon_max,
        lat_min=lat_min,
        lat_max=lat_max,
        auto_extent=auto_extent,
        extent_margin=extent_margin,
        contour_margin=contour_margin,
        global_plot=global_plot,
        center_lon=center_lon,
        input_radians=input_radians,
        lon_name_vertex=lon_name_vertex,
        lat_name_vertex=lat_name_vertex,
    )

    x = geo["x"]
    y = geo["y"]
    zz = geo["zz"]
    idx = geo["idx"]

    cmap = _infer_default_cmap(zz, cmap=cmap)
    contour_levels, vmin, vmax, extend = _resolve_contour_levels(
        zz,
        levels,
        vmin=vmin,
        vmax=vmax,
        symmetric=symmetric,
        percentile=percentile,
        extend=extend,
    )

    if ax is None:
        ax = plt.gca()

    cartopy_axis = _is_cartopy_axis(ax)

    if transform is None and cartopy_axis:
        transform = ccrs.PlateCarree()

    triang_cache_key = (
        "vertex_contour_triangulation",
        id(ds),
        lon_name_vertex,
        lat_name_vertex,
        geo["center_lon"],
        input_radians,
        geo["lon_min_plot"],
        geo["lon_max_plot"],
        geo["lat_min_plot"],
        geo["lat_max_plot"],
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

    _decorate_vertex_contour_axis(
        ax,
        cartopy_axis=cartopy_axis,
        global_plot=global_plot,
        lon_min_plot=geo["lon_min_plot"],
        lon_max_plot=geo["lon_max_plot"],
        lat_min_plot=geo["lat_min_plot"],
        lat_max_plot=geo["lat_max_plot"],
        xlabel=xlabel,
        ylabel=ylabel,
        aspect=aspect,
    )

    if title is None:
        title = str(geo["long_name"])

    ax.set_title(title)

    return cs
