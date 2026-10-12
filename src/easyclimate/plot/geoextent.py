"""Version-independent geographic map extents."""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from ._geo import _central_longitude, _set_cartopy_extent

__all__ = ["set_geoextent"]


def set_geoextent(extents, *, ax=None, crs=None):
    """Set a map's geographic window across Cartopy versions.

    Parameters
    ----------
    extents : sequence of four floats
        ``[lon_min, lon_max, lat_min, lat_max]`` in degrees. Longitudes
        may use either the -180--180 or 0--360 convention. A decreasing
        longitude pair, such as ``[170, -170]``, crosses the dateline.
        A span of 360 degrees, such as ``[0, 360]`` or ``[-180, 180]``,
        requests all longitudes. Equal longitude endpoints are rejected.
    ax : cartopy.mpl.geoaxes.GeoAxes, optional
        Target map axes. Defaults to the current axes.
    crs : cartopy.crs.CRS, optional
        Coordinate system of ``extents``. Defaults to ``PlateCarree()``.
        Shifted PlateCarree coordinates are supported. For other coordinate
        systems, delegates to ``ax.set_extent`` without longitude handling.

    Returns
    -------
    cartopy.mpl.geoaxes.GeoAxes
        The modified axes.

    Notes
    -----
    PlateCarree and Mercator windows are set in native projection coordinates
    to avoid coincident full-longitude edges and changes to Cartopy's
    PlateCarree central-meridian parameters. A window crossing the target
    projection's seam must show both map edges, so its bounding box spans
    the map width. Choose a central longitude inside the window to avoid this.
    Full-longitude bounds on other projections are computed from sampled
    geographic points; polar stereographic bounds are computed analytically.
    This function sets the view limits, without adding a circular boundary.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import cartopy.crs as ccrs
    >>> import easyclimate as ecl
    >>> fig, ax = plt.subplots(subplot_kw={"projection": ccrs.Mercator(140)})
    >>> ax = ecl.plot.set_geoextent([0, 360, 10, 85], ax=ax)
    """
    if ax is None:
        ax = plt.gca()
    if not hasattr(ax, "projection") or not hasattr(ax, "set_extent"):
        raise TypeError("ax must be a Cartopy GeoAxes.")
    values = np.asarray(extents, dtype=float)
    if values.shape != (4,) or not np.all(np.isfinite(values)):
        raise ValueError("extents must contain four finite values.")
    if crs is None:
        crs = ccrs.PlateCarree()
    if not isinstance(crs, ccrs.PlateCarree):
        ax.set_extent(values, crs=crs)
        return ax

    lon_min, lon_max, lat_min, lat_max = values
    if not -90 <= lat_min < lat_max <= 90:
        raise ValueError("Latitude bounds must satisfy -90 <= lat_min < lat_max <= 90.")
    span = lon_max - lon_min
    if span == 0 or abs(span) > 360:
        raise ValueError(
            "Longitude bounds must define a nonzero window of at most 360 degrees."
        )
    center = _central_longitude(crs)
    _set_cartopy_extent(ax, lon_min + center, lon_max + center, lat_min, lat_max)
    return ax
