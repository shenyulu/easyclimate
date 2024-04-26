"""
The eastern Atlantic (EA) pattern
"""

import xarray as xr
from ...core.utility import sort_ascending_latlon_coordinates
from ...core.variability import remove_seasonal_cycle_mean

__all__ = ["calc_index_EA_Wallace_Gutzler_1981"]


def calc_index_EA_Wallace_Gutzler_1981(
    z_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    time_dim: str = "time",
) -> xr.DataArray:
    """
    The calculation of monthly mean eastern Atlantic (EA) index using Pointwise method following Wallace and Gutzler (1981):

    .. math::
        \\mathrm{EA = \\frac{1}{2} Z^*(55^{\\circ}N, 20^{\\circ}W) - \\frac{1}{4} Z^*(25^{\\circ}N, 25^{\\circ}W) - \\frac{1}{4} Z^*(50^{\\circ}N, 40^{\\circ}E)}

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
    The monthly mean EA index (:py:class:`xarray.DataArray<xarray.DataArray>`).

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

    # Z*(55°N,20°W)
    part1 = z_anomaly_data.sel(lat=55, lon=340, method="nearest")
    # Z*(25°N,25°W)
    part2 = z_anomaly_data.sel(lat=25, lon=335, method="nearest")
    # Z*(50°N,40°E)
    part3 = z_anomaly_data.sel(lat=50, lon=320, method="nearest")
    index_PNA = 0.5 * part1 - 0.25 * part2 - 0.25 * part3

    # Normalized
    index_normalized_std = index_PNA.sel({time_dim: time_range}).std(dim=time_dim).data
    return (index_PNA / index_normalized_std).drop_vars("month")
