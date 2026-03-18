"""
Polar region mapping
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.path as mpath
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from typing import Literal

__all__ = ["draw_polar_basemap", "set_polar_title"]


def draw_polar_basemap(
    *,
    lat_range: list,
    add_gridlines: bool = True,
    lon_step: float = None,
    lat_step: float = None,
    ax=None,
    draw_labels: bool = True,
    lat_label_lon: float | None = None,
    lon_label_lat: float | None = None,
    lon_label_pad_deg: float = -4.0,
    set_map_boundary_kwargs: dict | None = None,
    gridlines_kwargs: dict | None = None,
    lat_label_kwargs: dict | None = None,
    lon_label_kwargs: dict | None = None,
):
    """
    Draw a polar stereographic circular basemap with optional gridlines and
    manually positioned coordinate labels.

    This function is designed for north and south polar stereographic maps.
    It first clips the map boundary into a circular polar view, then optionally
    adds longitude and latitude gridlines, and finally places readable manual
    latitude and longitude labels around the polar boundary.

    Parameters
    ----------
    lat_range: list[float] | tuple[float, float]
        Latitude range used to define the circular polar boundary, such as
        ``[50, 90]`` for the Northern Hemisphere or ``[-90, -50]`` for the
        Southern Hemisphere.
    add_gridlines: :py:class:`bool <bool>`, default: ``True``
        Whether to draw longitude and latitude gridlines.
    lon_step: float, optional
        Interval of longitude gridlines and longitude labels in degrees. This
        parameter must be specified when `add_gridlines=True` or
        `draw_labels=True`.
    lat_step: :py:class:`float <float>`, optional
        Interval of latitude gridlines and latitude labels in degrees. This
        parameter must be specified when `add_gridlines=True` or when manual
        latitude labels are enabled.
    ax: cartopy.mpl.geoaxes.GeoAxesSubplot, optional
        Target axes with a polar stereographic projection. (uses current axes if None)
    draw_labels: :py:class:`bool <bool>`, default: ``True``
        Whether to manually draw latitude and longitude labels.
    lat_label_lon: :py:class:`float <float>` | None, optional
        Longitude at which latitude labels are placed.
    lon_label_lat: :py:class:`float <float>` | None, optional
        Latitude at which longitude labels are placed. If not provided, it is
        inferred from the lower edge of `lat_range` plus `lon_label_pad_deg`.
    lon_label_pad_deg: :py:class:`float <float>`, default: ``-4.0``
        Padding, in degrees latitude, used to shift longitude labels relative
        to the outer latitude boundary.
    set_map_boundary_kwargs: :py:class:`dict <dict>` | None, optional
        Additional keyword arguments passed to when constructing the circular boundary.

        - north_pad : :py:class:`int <int>`
            A constant to be added to the second entry in lat_range. Use this
            if the northern edge of the plot is cut off. Defaults to 0.
        - south_pad : :py:class:`int <int>`
            A constant to be subtracted from the first entry in lat_range. Use
            this if the southern edge of the plot is cut off. Defaults to 0.
        - east_pad : :py:class:`int <int>`
            A constant to be added to the second entry in lon_range. Use this
            if the eastern edge of the plot is cut off. Defaults to 0.
        - west_pad : :py:class:`int <int>`
            A constant to be subtracted from the first entry in lon_range. Use
            this if the western edge of the plot is cut off. Defaults to 0.

    gridlines_kwargs: :py:class:`dict <dict>` | None, optional
        Additional keyword arguments passed to
        :py:func:`cartopy.mpl.geoaxes.GeoAxes.gridlines <cartopy.mpl.geoaxes.GeoAxes.gridlines>`
        when drawing gridlines.
    lat_label_kwargs: :py:class:`dict <dict>` | None, optional
        Additional text style options for latitude labels.
    lon_label_kwargs: :py:class:`dict <dict>` | None, optional
        Additional text style options for longitude labels.

    Returns
    -------
    gl: :py:class:`cartopy.mpl.gridliner.Gridliner <cartopy.mpl.gridliner.Gridliner>` or None
        The gridliner object returned by Cartopy. If `add_gridlines=False`,
        ``None`` is returned.
    meta: :py:class:`dict <dict>`
        Metadata dictionary describing the polar label layout. It can be passed
        directly to :py:func:`set_polar_title <set_polar_title>`
        to help position titles above the circular map.

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_formatting_coordinates.py
        ./dynamic_docs/plot_da_bbo.py
        ./dynamic_docs/plot_curve_quiver.py
        ./dynamic_docs/plot_waf.py.py
    """

    if set_map_boundary_kwargs is None:
        set_map_boundary_kwargs = {"north_pad": 0.3, "south_pad": 0.4}
    if gridlines_kwargs is None:
        gridlines_kwargs = {"color": "grey", "alpha": 0.5, "linestyle": "--"}
    if lat_label_kwargs is None:
        lat_label_kwargs = {}
    if lon_label_kwargs is None:
        lon_label_kwargs = {}

    if ax is None:
        ax = plt.gca()

    if not isinstance(ax.projection, (ccrs.NorthPolarStereo, ccrs.SouthPolarStereo)):
        raise TypeError("Projection must be NorthPolarStereo or SouthPolarStereo.")

    lat_min = float(min(lat_range))
    lat_max = float(max(lat_range))
    boundary_lat_range = [lat_min, lat_max]
    is_northern_hemisphere = isinstance(ax.projection, ccrs.NorthPolarStereo)
    outer_boundary_lat = lat_min if is_northern_hemisphere else lat_max

    # 1) circular boundary
    theta = np.linspace(0, 2 * np.pi, 361)
    circle_center = np.array([0.5, 0.5])
    circle_radius = 0.5
    circular_path = mpath.Path(
        np.column_stack([np.sin(theta), np.cos(theta)]) * circle_radius + circle_center
    )
    ax.set_boundary(circular_path, transform=ax.transAxes)

    extent_kwargs = dict(set_map_boundary_kwargs)
    north_pad = float(extent_kwargs.pop("north_pad", 0.0))
    south_pad = float(extent_kwargs.pop("south_pad", 0.0))
    east_pad = float(extent_kwargs.pop("east_pad", 0.0))
    west_pad = float(extent_kwargs.pop("west_pad", 0.0))
    if extent_kwargs:
        raise TypeError(
            "Unsupported set_map_boundary_kwargs for draw_polar_basemap: "
            f"{', '.join(extent_kwargs.keys())}"
        )

    ax.set_extent(
        [
            -180 - west_pad,
            180 + east_pad,
            lat_min - south_pad,
            lat_max + north_pad,
        ],
        crs=ccrs.PlateCarree(),
    )

    gl = None
    if add_gridlines:
        if lon_step is None or lat_step is None:
            raise ValueError(
                "If add_gridlines=True, lon_step and lat_step must be specified."
            )

        xlocs = mticker.FixedLocator(np.arange(-180, 180, lon_step))
        ylocs = mticker.FixedLocator(np.arange(lat_min, lat_max + 0.1, lat_step))

        gl = ax.gridlines(
            crs=ccrs.PlateCarree(),
            draw_labels=False,
            xlocs=xlocs,
            ylocs=ylocs,
            **gridlines_kwargs,
        )

    # 2) manual labels
    pc = ccrs.PlateCarree()
    if draw_labels:
        lat_fmt = LatitudeFormatter()
        lon_fmt = LongitudeFormatter(zero_direction_label=False)

        # 2.1 latitude labels along fixed longitude
        if lat_label_lon is not None:
            lats = np.arange(lat_min, lat_max + 0.1, lat_step)
            default_lat_bbox = dict(
                facecolor="white", alpha=0.8, edgecolor="none", pad=1.0
            )

            for lat in lats:
                if lat <= lat_min or lat >= lat_max:
                    continue
                ax.text(
                    lat_label_lon,
                    lat,
                    lat_fmt(lat),
                    transform=pc,
                    ha=lat_label_kwargs.get("ha", "center"),
                    va=lat_label_kwargs.get("va", "center"),
                    fontsize=lat_label_kwargs.get("fontsize", 10),
                    bbox=lat_label_kwargs.get("bbox", default_lat_bbox),
                    zorder=lat_label_kwargs.get("zorder", 20),
                )

        # 2.2 longitude labels around boundary
        base_lat = outer_boundary_lat

        if lon_label_lat is None:
            lat_pad = float(lon_label_pad_deg)
            lon_label_lat = (
                base_lat + lat_pad if is_northern_hemisphere else base_lat - lat_pad
            )
        else:
            lon_label_lat = float(lon_label_lat)

        lons = np.arange(-180, 180, lon_step)
        default_lon_bbox = dict(facecolor="white", alpha=0.8, edgecolor="none", pad=1.0)
        fontsize = lon_label_kwargs.get("fontsize", 10)

        for lon in lons:
            # rotation computed in projection plane
            x, y = ax.projection.transform_point(lon, lon_label_lat, pc)
            ang = np.degrees(np.arctan2(y, x))
            rot = ang - 90.0

            # keep readable
            if rot > 90:
                rot -= 180
            if rot < -90:
                rot += 180

            ax.text(
                lon,
                lon_label_lat,
                lon_fmt(lon),
                transform=pc,
                rotation=rot,
                rotation_mode="anchor",
                ha=lon_label_kwargs.get("ha", "center"),
                va=lon_label_kwargs.get("va", "center"),
                fontsize=fontsize,
                bbox=lon_label_kwargs.get("bbox", default_lon_bbox),
                zorder=lon_label_kwargs.get("zorder", 20),
            )

    meta = dict(
        lon_step=lon_step,
        lon_label_lat=lon_label_lat,
        lat_range=boundary_lat_range,
        lon_label_pad_deg=lon_label_pad_deg,
    )
    return gl, meta


import numpy as np
import cartopy.crs as ccrs


def set_polar_title(
    title: str,
    meta: dict,
    ax=None,
    *,
    loc: Literal["left", "center", "right"] = "center",
    dy: float | None = None,
    xpad: float = 0.02,
    auto_dy_factor: float = 1.2,
    **text_kwargs,
):
    """
    Add a title above a polar circular map using the layout metadata returned
    by :py:func:`draw_polar_basemap <draw_polar_basemap>`.

    This function estimates the uppermost position of manually placed longitude
    labels on a polar map, then places the title slightly above them so that
    the title does not overlap with the circular boundary labels.

    Parameters
    ----------
    title: :py:class:`str <str>`
        Title text.
    meta: :py:class:`dict <dict>`
        Metadata dictionary returned by :py:func:`draw_polar_basemap <draw_polar_basemap>`.
        It must at least contain `lon_step` and `lon_label_lat`.
    ax: cartopy.mpl.geoaxes.GeoAxesSubplot
        Polar stereographic axes on which the title will be drawn.
    loc: {"left", "center", "right"}, default: ``"center"``
        Horizontal alignment of the title.
    dy: :py:class:`float <float>` | None, optional
        Vertical offset in axes coordinates relative to the topmost longitude
        label position. If ``None``, the offset is estimated automatically from
        the font size.
    xpad: :py:class:`float <float>`, default: ``0.02``
        Horizontal padding used when `loc` is ``"left"`` or ``"right"``.
    auto_dy_factor: :py:class:`float <float>`, default: ``1.2``
        Multiplicative factor used when automatically converting font height
        into a safe vertical title offset.
    **text_kwargs
        Additional keyword arguments passed to
        :py:func:`matplotlib.axes.Axes.text <matplotlib.axes.Axes.text>`.

    Returns
    -------
    :py:class:`matplotlib.text.Text <matplotlib.text.Text>`
        The created title text artist.

    .. note::
        If `dy` is ``None``, the function automatically computes a reasonable
        vertical offset based on the rendered font size, which is usually the most
        robust option for polar circular maps.

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_formatting_coordinates.py
        ./dynamic_docs/plot_waf.py.py
    """
    if ax is None:
        ax = plt.gca()

    lon_step = meta["lon_step"]
    lon_label_lat = meta["lon_label_lat"]

    pc = ccrs.PlateCarree()
    lons = np.arange(-180, 180, lon_step)

    # Find the highest y-coordinate of the Axes with the longitude label.
    ys = []
    for lon in lons:
        x_proj, y_proj = ax.projection.transform_point(lon, lon_label_lat, pc)
        x_disp, y_disp = ax.transData.transform((x_proj, y_proj))
        x_ax, y_ax = ax.transAxes.inverted().transform((x_disp, y_disp))
        ys.append(y_ax)

    y_top = np.nanmax(ys)

    # Automatic dy calculation
    if dy is None:
        # Obtain font size (priority: text_kwargs, then rcParams)
        fontsize = text_kwargs.get(
            "fontsize",
            (
                ax.figure._cachedRenderer.points_to_pixels(10)  # fallback
                if hasattr(ax.figure, "_cachedRenderer")
                else 12
            ),
        )

        # More stable method: Use matplotlib's "points -> pixels"
        fig = ax.figure
        renderer = fig.canvas.get_renderer()

        # Font height (points -> pixels)
        fontsize_pt = text_kwargs.get("fontsize", 12)
        fontsize_px = fontsize_pt * fig.dpi / 72.0

        # Additional safety factor
        offset_px = fontsize_px * auto_dy_factor

        # Move y_top to the display coordinates
        x_disp, y_disp = ax.transAxes.transform((0, y_top))

        # Add pixel offset
        y_disp_new = y_disp + offset_px

        # Return to Axes coordinates
        _, y_title = ax.transAxes.inverted().transform((x_disp, y_disp_new))

    else:
        y_title = y_top + dy

    # Position X
    loc = loc.lower()
    if loc == "left":
        x_title, ha = xpad, "left"
    elif loc == "right":
        x_title, ha = 1 - xpad, "right"
    else:
        x_title, ha = 0.5, "center"

    defaults = dict(
        transform=ax.transAxes,
        ha=ha,
        va="bottom",
        clip_on=False,
        zorder=50,
    )

    defaults.update(text_kwargs)

    return ax.text(x_title, y_title, title, **defaults)
