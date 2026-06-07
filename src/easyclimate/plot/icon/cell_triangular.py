"""
Triangular cell plots for ICON native grids.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import cartopy.crs as ccrs

from .common import (
    _rad2deg,
    _lon_wrap_to_center,
    _is_cartopy_axis,
    _infer_default_cmap,
    _prepare_cell_geometry,
    _set_cartopy_extent,
)


__all__ = ["plot_cell_triangular", "plot_icon_triangles"]


def plot_cell_triangular(
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
    center_lon=None,
    input_radians=True,
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
    # ICON mesh variable names
    lon_name_cell="clon",
    lat_name_cell="clat",
    lon_bnds_name="clon_bnds",
    lat_bnds_name="clat_bnds",
):
    """
    Plot ICON cell-centered scalar data on native triangular cells.

    Parameters
    ----------
    ds : :py:class:`xarray.Dataset <xarray.Dataset>`
        ICON dataset containing cell centers and triangular cell bounds.
    da : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`numpy.ndarray <numpy.ndarray>`
        Cell-centered scalar values on ``ncells``. Select time and vertical
        dimensions before calling, for example
        ``ds["temp"].isel(time=0, plev=0)``.
    ax : :py:class:`matplotlib.axes.Axes <matplotlib.axes.Axes>`, optional
        Axes on which to draw. By default, use the current axes.
    transform : :py:class:`cartopy.crs.CRS <cartopy.crs.CRS>`, optional
        Coordinate reference system of the input coordinates for Cartopy axes.
    lon_min, lon_max, lat_min, lat_max : :py:class:`float <float>`, optional
        Plot extent in degrees. Cross-dateline longitude windows are supported.
    auto_extent : :py:class:`bool <bool>`, default: True
        If True and no explicit extent is supplied, infer the extent from valid
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
    cmap : :py:class:`str <str>` or :py:class:`matplotlib.colors.Colormap <matplotlib.colors.Colormap>`, optional
        Colormap used for triangle values.
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
    xlabel, ylabel : :py:class:`str <str>`
        Axis labels used for plain Matplotlib axes.
    aspect : :py:class:`str <str>` or :py:class:`float <float>`, default: "auto"
        Aspect setting used for plain Matplotlib axes.
    lon_name_cell, lat_name_cell, lon_bnds_name, lat_bnds_name : :py:class:`str <str>`
        ICON coordinate and triangular-bound variable names.

    Returns
    -------
    :py:class:`matplotlib.collections.PolyCollection <matplotlib.collections.PolyCollection>`
        Polygon collection added to the axes.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        icon_cell_triangular
    """

    lon_bnds_raw = np.asarray(ds[lon_bnds_name].values)
    lat_bnds = np.asarray(ds[lat_bnds_name].values)

    if input_radians:
        lon_bnds_raw = _rad2deg(lon_bnds_raw)
        lat_bnds = _rad2deg(lat_bnds)

    geo = _prepare_cell_geometry(
        ds,
        da,
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
    )

    z = geo["z"]
    units = geo["units"]
    long_name = geo["long_name"]
    valid = geo["valid"]
    center_lon = geo["center_lon"]
    lon_cell = geo["lon_cell"]
    lat_cell = geo["lat_cell"]
    lon_min_plot = geo["lon_min_plot"]
    lon_max_plot = geo["lon_max_plot"]
    lat_min_plot = geo["lat_min_plot"]
    lat_max_plot = geo["lat_max_plot"]

    n_cells = lon_cell.size

    if lon_bnds_raw.shape[0] != n_cells or lat_bnds.shape[0] != n_cells:
        raise ValueError(
            "`lon_bnds_name` and `lat_bnds_name` must have ncells as their "
            "first dimension."
        )

    lon_bnds = _lon_wrap_to_center(lon_bnds_raw, center_lon)

    # ---- axis ----
    if ax is None:
        ax = plt.gca()

    fig = ax.figure

    cartopy_axis = _is_cartopy_axis(ax)

    if transform is None and cartopy_axis:
        transform = ccrs.PlateCarree()

    candidate = (
        valid
        & np.all(np.isfinite(lon_bnds), axis=1)
        & np.all(np.isfinite(lat_bnds), axis=1)
        & np.isfinite(lon_cell)
        & np.isfinite(lat_cell)
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

    polys = []
    values = []

    for c in np.flatnonzero(candidate):
        val = z[c]
        lons = lon_bnds[c]
        lats = lat_bnds[c]
        polys.append(np.column_stack([lons, lats]))
        values.append(val)

    values = np.asarray(values)

    if len(polys) == 0:
        raise RuntimeError(
            "No ICON triangular cells selected. Check extent or coordinate units."
        )

    cmap = _infer_default_cmap(values, cmap=cmap)

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
