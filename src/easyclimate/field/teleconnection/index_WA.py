"""
The western Atlantic (WA) pattern
"""

import xarray as xr
from ...core.utility import sort_ascending_latlon_coordinates
from ...core.variability import remove_seasonal_cycle_mean

__all__ = ["calc_index_WA_Wallace_Gutzler_1981"]


def calc_index_WA_Wallace_Gutzler_1981(
    z_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    time_dim: str = "time",
) -> xr.DataArray:
    """
    The calculation of monthly mean western Atlantic (WA) index using Pointwise method following Wallace and Gutzler (1981):

    .. math::
        \\mathrm{WA = \\frac{1}{2} [Z^*(55^{\\circ}N, 55^{\\circ}W) - Z^*(30^{\\circ}N, 55^{\\circ}W)] }

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
    The monthly mean WA index (:py:class:`xarray.DataArray<xarray.DataArray>`).

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

    # Z*(55°N,55°W)
    part1 = z_anomaly_data.sel(lat=55, lon=305, method="nearest")
    # Z*(30°N,55°W)
    part2 = z_anomaly_data.sel(lat=30, lon=305, method="nearest")
    index_WA = 0.5 * (part1 - part2)

    # Normalized
    index_normalized_std = index_WA.sel({time_dim: time_range}).std(dim=time_dim).data
    return (index_WA / index_normalized_std).drop_vars("month")
