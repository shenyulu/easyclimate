"""
The Eurasian pattern (EU) pattern
"""

import xarray as xr
from ...core.utility import sort_ascending_latlon_coordinates
from ...core.variability import remove_seasonal_cycle_mean

__all__ = ["calc_index_EU_Wallace_Gutzler_1981"]


def calc_index_EU_Wallace_Gutzler_1981(
    z_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    time_dim: str = "time",
) -> xr.DataArray:
    """
    The calculation of monthly mean Eurasian pattern (EU) index using Pointwise method following Wallace and Gutzler (1981):

    .. math::
        \\mathrm{EU = - \\frac{1}{4} Z^*(55^{\\circ}N, 20^{\\circ}E) + \\frac{1}{2} Z^*(55^{\\circ}N, 75^{\\circ}E) - \\frac{1}{4} Z^*(40^{\\circ}N, 145^{\\circ}E)}

    where :math:`Z^*` denotes monthly mean 500 hPa height anomaly.

    Parameters
    ----------
    z_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The monthly geopotential height. The 500hPa layer is recommended.
    time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
        The time range of seasonal cycle means to be calculated. The default value is the entire time range.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.

    Returns
    -------
    The monthly mean EU index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - `Wallace, J. M., & Gutzler, D. S. (1981). Teleconnections in the Geopotential Height Field during the Northern Hemisphere Winter. Monthly Weather Review, 109(4), 784-812. <https://journals.ametsoc.org/view/journals/mwre/109/4/1520-0493_1981_109_0784_titghf_2_0_co_2.xml>`__
    """
    z_monthly_data = sort_ascending_latlon_coordinates(
        z_monthly_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    # anomaly
    z_anomaly_data = remove_seasonal_cycle_mean(
        z_monthly_data, dim=time_dim, time_range=time_range
    )

    # Scandinavia and Poland: Z*(55°N,20°E)
    part1 = z_anomaly_data.sel(lat=55, lon=20, method="nearest")
    # Siberia: Z*(55°N,75°E)
    part2 = z_anomaly_data.sel(lat=55, lon=75, method="nearest")
    # Japan: Z*(40°N,145°E)
    part3 = z_anomaly_data.sel(lat=40, lon=145, method="nearest")
    index_EU = -0.25 * part1 + 0.5 * part2 - 0.25 * part3

    # Normalized
    index_normalized_std = index_EU.sel({time_dim: time_range}).std(dim=time_dim).data
    return (index_EU / index_normalized_std).drop_vars("month")
