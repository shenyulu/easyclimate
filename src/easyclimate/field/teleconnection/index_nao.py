"""
North Atlantic Oscillation (NAO) Index

The North Atlantic Oscillation (NAO) index is based on the surface sea-level pressure difference between the Subtropical (Azores) High and the Subpolar Low. The positive phase of the NAO reflects below-normal heights and pressure across the high latitudes of the North Atlantic and above-normal heights and pressure over the central North Atlantic, the eastern United States and western Europe. The negative phase reflects an opposite pattern of height and pressure anomalies over these regions. Both phases of the NAO are associated with basin-wide changes in the intensity and location of the North Atlantic jet stream and storm track, and in large-scale modulations of the normal patterns of zonal and meridional heat and moisture transport, which in turn results in changes in temperature and precipitation patterns often extending from eastern North America to western and central Europe.

Strong positive phases of the NAO tend to be associated with above-normal temperatures in the eastern United States and across northern Europe and below-normal temperatures in Greenland and oftentimes across southern Europe and the Middle East. They are also associated with above-normal precipitation over northern Europe and Scandinavia and below-normal precipitation over southern and central Europe. Opposite patterns of temperature and precipitation anomalies are typically observed during strong negative phases of the NAO. During particularly prolonged periods dominated by one particular phase of the NAO, abnormal height and temperature patterns are also often seen extending well into central Russia and north-central Siberia. The NAO exhibits considerable interseasonal and interannual variability, and prolonged periods (several months) of both positive and negative phases of the pattern are common.

https://www.ncei.noaa.gov/access/monitoring/nao/
"""

import xarray as xr
from ...core.utility import sort_ascending_latlon_coordinates
from ...core.variability import remove_seasonal_cycle_mean
from ...core.eof import get_REOF_model, calc_REOF_analysis

__all__ = [
    "calc_index_NAO_NH_REOF",
]


def calc_index_NAO_NH_REOF(
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
    The calculation of monthly mean NAO index using rotated empirical orthogonal functions (REOFs) method:

    Parameters
    ----------
    z_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The monthly geopotential height.
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
        Solver to use for the SVD computation.
    solver_kwargs: :py:class:`dict<dict>`, default `{}`.
        Additional keyword arguments to be passed to the SVD solver.

    Returns
    -------
    The monthly mean NAO index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - `Soulard, N., Lin, H. The spring relationship between the Pacific-North American pattern and the North Atlantic Oscillation. Clim Dyn 48, 619–629 (2017). <https://doi.org/10.1007/s00382-016-3098-3>`__
    - https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/history/method.shtml
    - https://www.ncei.noaa.gov/access/monitoring/nao

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
    index_NAO = z_REOF_result["PC"].sel(mode=1)

    # Normalized
    index_normalized_std = index_NAO.sel({time_dim: time_range}).std(dim=time_dim).data
    return (index_NAO / index_normalized_std).drop_vars(["month", "mode"])
