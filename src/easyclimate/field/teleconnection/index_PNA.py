"""
Pacific/North American (PNA) Index

The Pacific/North American teleconnection pattern (PNA) is one of the most prominent modes of low-frequency
variability in the Northern Hemisphere extratropics. The positive phase of the PNA pattern features above-average
heights in the vicinity of Hawaii and over the intermountain region of North America, and below-average heights
located south of the Aleutian Islands and over the southeastern United States. The PNA pattern is associated with
strong fluctuations in the strength and location of the East Asian jet stream. The positive phase is associated
with an enhanced East Asian jet stream and with an eastward shift in the jet exit region toward the western
United States. The negative phase is associated with a westward retraction of that jet stream toward eastern Asia,
blocking activity over the high latitudes of the North pacific, and a strong split-flow configuration over the central North Pacific.

The positive phase of the PNA pattern is associated with above-average temperatures over western Canada
and the extreme western United States, and below-average temperatures across the south-central and
southeastern U.S. The PNA tends to have little impact on surface temperature variability over
North America during summer. The associated precipitation anomalies include above-average
totals in the Gulf of Alaska extending into the Pacific Northwestern United States,
and below-average totals over the upper Midwestern United States.

Although the PNA pattern is a natural internal mode of climate variability,
it is also strongly influenced by the El Niño/ Southern Oscillation (ENSO) phenomenon.
The positive phase of the PNA pattern tends to be associated with Pacific warm episodes (El Niño),
and the negative phase tends to be associated with Pacific cold episodes (La Niña).

https://www.cpc.ncep.noaa.gov/data/teledoc/pna.shtml
"""

import xarray as xr
from ...core.utility import sort_ascending_latlon_coordinates
from ...core.variability import remove_seasonal_cycle_mean
from ...core.eof import get_REOF_model, calc_REOF_analysis

__all__ = [
    "calc_index_PNA_modified_pointwise",
    "calc_index_PNA_Wallace_Gutzler_1981",
    "calc_index_PNA_NH_REOF",
]


def calc_index_PNA_modified_pointwise(
    z_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    time_dim: str = "time",
) -> xr.DataArray:
    """
    The calculation of monthly mean PNA index is constructed by following modified pointwise method:

    .. math::
        \\mathrm{PNA = Z^*(15^{\\circ}N-25^{\\circ}N,180-140^{\\circ}W)-Z^*(40^{\\circ}N-50^{\\circ}N,180-140^{\\circ}W)+Z^*(45^{\\circ}N-60^{\\circ}N,125^{\\circ}W-105^{\\circ}W)-Z^*(25^{\\circ}N-35^{\\circ}N,90^{\\circ}W-70^{\\circ}W)}

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
    The monthly mean PNA index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/month_pna_index2.shtml
    """
    z_monthly_data = sort_ascending_latlon_coordinates(
        z_monthly_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    # anomaly
    z_anomaly_data = remove_seasonal_cycle_mean(
        z_monthly_data, dim=time_dim, time_range=time_range
    )

    # Z*(15°N-25°N,180-140°W)
    part1 = z_anomaly_data.sel(lat=slice(15, 25), lon=slice(180, 220)).mean(
        dim=(lat_dim, lon_dim)
    )
    # Z*(40°N-50°N,180-140°W)
    part2 = z_anomaly_data.sel(lat=slice(40, 50), lon=slice(180, 220)).mean(
        dim=(lat_dim, lon_dim)
    )
    # Z*(45°N-60°N,125°W-105°W)
    part3 = z_anomaly_data.sel(lat=slice(45, 60), lon=slice(235, 255)).mean(
        dim=(lat_dim, lon_dim)
    )
    # Z*(25°N-35°N,90°W-70°W)
    part4 = z_anomaly_data.sel(lat=slice(25, 35), lon=slice(270, 290)).mean(
        dim=(lat_dim, lon_dim)
    )
    # Summation
    index_PNA = (part1 - part2 + part3 - part4) / 4

    # Normalized
    index_normalized_std = index_PNA.sel({time_dim: time_range}).std(dim=time_dim).data
    result = (index_PNA / index_normalized_std).drop_vars("month")
    return result


def calc_index_PNA_Wallace_Gutzler_1981(
    z_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    time_dim: str = "time",
) -> xr.DataArray:
    """
    The calculation of monthly mean PNA index using Pointwise method following Wallace and Gutzler (1981):

    .. math::
        \\mathrm{PNA = Z^*(20^{\\circ}N,160^{\\circ}W)-Z^*(45^{\\circ}N,165^{\\circ}W)+Z^*(55^{\\circ}N,115^{\\circ}W)-Z^*(30^{\\circ}N,85^{\\circ}W)}

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
    The monthly mean PNA index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - `Wallace, J. M., & Gutzler, D. S. (1981). Teleconnections in the Geopotential Height Field during the Northern Hemisphere Winter. Monthly Weather Review, 109(4), 784-812. <https://journals.ametsoc.org/view/journals/mwre/109/4/1520-0493_1981_109_0784_titghf_2_0_co_2.xml>`__
    - https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/month_pna_index2.shtml
    """
    z_monthly_data = sort_ascending_latlon_coordinates(
        z_monthly_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    # anomaly
    z_anomaly_data = remove_seasonal_cycle_mean(
        z_monthly_data, dim=time_dim, time_range=time_range
    )

    # Z*(20°N,160°W)
    part1 = z_anomaly_data.sel(lat=20, lon=200, method="nearest")
    # Z*(45°N,165°W)
    part2 = z_anomaly_data.sel(lat=45, lon=195, method="nearest")
    # Z*(55°N,115°W)
    part3 = z_anomaly_data.sel(lat=55, lon=245, method="nearest")
    # Z*(30°N,85°W)
    part4 = z_anomaly_data.sel(lat=30, lon=275, method="nearest")
    index_PNA = (part1 - part2 + part3 - part4) / 4

    # Normalized
    index_normalized_std = index_PNA.sel({time_dim: time_range}).std(dim=time_dim).data
    result = (index_PNA / index_normalized_std).drop_vars("month")
    return result


def calc_index_PNA_NH_REOF(
    z_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    lat_range: slice = slice(20, 85),
    time_dim: str = "time",
    random_state=None,
    solver="auto",
    solver_kwargs={},
) -> xr.DataArray:
    """
    The calculation of monthly mean PNA index using rotated empirical orthogonal functions (REOFs) method over the entire Northern Hemisphere:

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
    lat_range: :py:class:`slice <slice>`, default: `slice(20, 85)`.
        The latitude range of computation using REOFs over the Northern Hemisphere. The default value is from :math:`\\mathrm{20^{\\circ}N}` to :math:`\\mathrm{85^{\\circ}N}`.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.
    random_state: :py:class:`int<int>`, default `None`.
        Seed for the random number generator.
    solver: {"auto", "full", "randomized"}, default: "auto".
        Solver to use for the REOFs computation.
    solver_kwargs: :py:class:`dict<dict>`, default `{}`.
        Additional keyword arguments to be passed to the REOFs solver.

    Returns
    -------
    The monthly mean PNA index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - `Rodionov, S., & Assel, R. (2001). A new look at the Pacific/North American index. Geophysical Research Letters, 28(8), 1519-1522. <https://doi.org/10.1029/2000GL012185>`__
    - `Soulard, N., Lin, H. The spring relationship between the Pacific-North American pattern and the North Atlantic Oscillation. Clim Dyn 48, 619–629 (2017). <https://doi.org/10.1007/s00382-016-3098-3>`__
    - https://www.ncei.noaa.gov/access/monitoring/pna/
    - https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/history/method.shtml

    .. seealso::
        :py:func:`get_REOF_model <easyclimate.core.eof.get_REOF_model>`
    """
    z_monthly_data = sort_ascending_latlon_coordinates(
        z_monthly_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    # anomaly
    z_monthly_data_NH = z_monthly_data.sel(lat=lat_range)
    z_anomaly_data_NH = remove_seasonal_cycle_mean(
        z_monthly_data_NH, dim=time_dim, time_range=time_range
    )

    # REOF
    z_REOF_model = get_REOF_model(
        z_anomaly_data_NH,
        lat_dim=lat_dim,
        lon_dim=lon_dim,
        time_dim=time_dim,
        n_modes=2,
        random_state=random_state,
        solver=solver,
        solver_kwargs=solver_kwargs,
    )
    z_REOF_result = calc_REOF_analysis(z_REOF_model)
    index_PNA = z_REOF_result["PC"].sel(mode=2)

    # Normalized
    index_normalized_std = index_PNA.sel({time_dim: time_range}).std(dim=time_dim).data
    result = (index_PNA / index_normalized_std).drop_vars("month")
    return result
