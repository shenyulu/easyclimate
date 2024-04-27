"""
Circumglobal teleconnection pattern (CGT) Index
"""

import xarray as xr
from ...core.utility import sort_ascending_latlon_coordinates
from ...core.variability import remove_seasonal_cycle_mean
from ...core.eof import get_EOF_model, calc_EOF_analysis

__all__ = [
    "calc_index_CGT_1point_Ding_Wang_2005",
    "calc_index_CGT_NH_EOF2_Ding_Wang_2005",
]


def calc_index_CGT_1point_Ding_Wang_2005(
    z200_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    time_dim: str = "time",
) -> xr.DataArray:
    """
    The calculation of monthly mean circumglobal teleconnection pattern (CGT) index is constructed by following method:

        Z200 anomalies averaged over the area (:math:`\\mathrm{ 35 ^{\\circ}N - 40 ^{\\circ}N; 60^{\\circ}-70^{\\circ}E }`).

    Parameters
    ----------
    z200_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The monthly 200-hPa geopotential height.
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
    The monthly mean CGT index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - `Ding, Q., & Wang, B. (2005). Circumglobal Teleconnection in the Northern Hemisphere Summer. Journal of Climate, 18(17), 3483-3505. <https://doi.org/10.1175/JCLI3473.1>`__
    - `Liu Y. Relationship between the Silk Road and Circumglobal Teleconnection Patterns on the Interannual and Interdecadal Timescales. Atmosphere. 2023; 14(11):1626. <https://doi.org/10.3390/atmos14111626>`__

    """
    z200_monthly_data = sort_ascending_latlon_coordinates(
        z200_monthly_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    # anomaly
    z_anomaly_data = remove_seasonal_cycle_mean(
        z200_monthly_data, dim=time_dim, time_range=time_range
    )

    # Z200*(35°N-40°N, 60-70°E)
    index_CGT = z_anomaly_data.sel(lat=slice(35, 40), lon=slice(60, 70)).mean(
        dim=(lat_dim, lon_dim)
    )

    # Normalized
    index_normalized_std = index_CGT.sel({time_dim: time_range}).std(dim=time_dim).data
    result = (index_CGT / index_normalized_std).drop_vars("month")
    return result


def calc_index_CGT_NH_EOF2_Ding_Wang_2005(
    z200_monthly_data: xr.DataArray,
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
    The calculation of monthly mean circumglobal teleconnection pattern (CGT) index using empirical orthogonal functions (EOFs) method over the entire Northern Hemisphere:

    Parameters
    ----------
    z200_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The monthly 200-hPa geopotential height.
    time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
        The time range of seasonal cycle means to be calculated. The default value is the entire time range.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    lat_range: :py:class:`slice <slice>`, default: `slice(20, 85)`.
        The latitude range of computation using EOFs over the Northern Hemisphere. The default value is from :math:`\\mathrm{20^{\\circ}N}` to :math:`\\mathrm{85^{\\circ}N}`.
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
    The monthly mean CGT index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - `Ding, Q., & Wang, B. (2005). Circumglobal Teleconnection in the Northern Hemisphere Summer. Journal of Climate, 18(17), 3483-3505. <https://doi.org/10.1175/JCLI3473.1>`__
    - `Liu Y. Relationship between the Silk Road and Circumglobal Teleconnection Patterns on the Interannual and Interdecadal Timescales. Atmosphere. 2023; 14(11):1626. <https://doi.org/10.3390/atmos14111626>`__

    .. seealso::
        :py:func:`get_EOF_model <easyclimate.core.eof.get_EOF_model>`
    """
    z200_monthly_data = sort_ascending_latlon_coordinates(
        z200_monthly_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    # anomaly
    z_monthly_data_NH = z200_monthly_data.sel(lat=lat_range)
    z_anomaly_data_NH = remove_seasonal_cycle_mean(
        z_monthly_data_NH, dim=time_dim, time_range=time_range
    )

    # REOF
    z_REOF_model = get_EOF_model(
        z_anomaly_data_NH,
        lat_dim=lat_dim,
        lon_dim=lon_dim,
        time_dim=time_dim,
        n_modes=2,
        random_state=random_state,
        solver=solver,
        solver_kwargs=solver_kwargs,
    )
    z_REOF_result = calc_EOF_analysis(z_REOF_model)
    index_PNA = z_REOF_result["PC"].sel(mode=2)

    # Normalized
    index_normalized_std = index_PNA.sel({time_dim: time_range}).std(dim=time_dim).data
    result = (index_PNA / index_normalized_std).drop_vars("month")
    return result
