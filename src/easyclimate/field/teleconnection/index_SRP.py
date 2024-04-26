"""
Silk Road pattern (SRP) Index
"""

import xarray as xr
from ...core.utility import sort_ascending_latlon_coordinates
from ...core.variability import remove_seasonal_cycle_mean
from ...core.eof import get_EOF_model, calc_EOF_analysis

__all__ = [
    "calc_index_SRP_EOF1_Yasui_Watanabe_2010",
    "calc_index_SRP_EOF1_Kosaka_2009",
    "calc_index_SRP_EOF1_Chen_Huang_2012",
    "calc_index_SRP_EOF1_Sato_Takahashi_2006",
    "calc_index_SRP_1point_Lu_2002",
]


def calc_index_SRP_EOF1_Yasui_Watanabe_2010(
    v200_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    lat_range: slice = slice(20, 60),
    lon_range: slice = slice(0, 150),
    time_dim: str = "time",
    random_state=None,
    solver="auto",
    solver_kwargs={},
) -> xr.DataArray:
    """
    The calculation of monthly mean SRP index using empirical orthogonal functions (EOFs) method based on Yasui and Watanabe (2010):

        EOF1 of V200 over (:math:`\\mathrm{20 ^{\\circ}N - 60 ^{\\circ}N; 0 ^{\\circ} - 150 ^{\\circ}E}`).

    Parameters
    ----------
    v200_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The monthly 200-hPa meridional wind.
    time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
        The time range of seasonal cycle means to be calculated. The default value is the entire time range.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    lat_range: :py:class:`slice <slice>`, default: `slice(20, 60)`.
        The latitude range of computation using EOFs over the region. The default value is from :math:`\\mathrm{20^{\\circ}N}` to :math:`\\mathrm{60^{\\circ}N}`.
    lon_range: :py:class:`slice <slice>`, default: `slice(0, 150)`.
        The longitude range of computation using EOFs over the region. The default value is from :math:`\\mathrm{0^{\\circ}}` to :math:`\\mathrm{150^{\\circ}E}`.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.
    random_state: :py:class:`int<int>`, default `None`.
        Seed for the random number generator.
    solver: {"auto", "full", "randomized"}, default: "auto".
        Solver to use for the EOFs computation.
    solver_kwargs: :py:class:`dict<dict>`, default `{}`.
        Additional keyword arguments to be passed to the EOFs solver.

    Returns
    -------
    The monthly mean SRP index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - `Yasui, S., & Watanabe, M. (2010). Forcing Processes of the Summertime Circumglobal Teleconnection Pattern in a Dry AGCM. Journal of Climate, 23(8), 2093-2114. <https://doi.org/10.1175/2009JCLI3323.1>`__
    - `Liu Y. Relationship between the Silk Road and Circumglobal Teleconnection Patterns on the Interannual and Interdecadal Timescales. Atmosphere. 2023; 14(11):1626. <https://doi.org/10.3390/atmos14111626>`__

    .. seealso::
        :py:func:`get_EOF_model <easyclimate.core.eof.get_EOF_model>`
    """
    v200_monthly_data = sort_ascending_latlon_coordinates(
        v200_monthly_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    # anomaly
    v200_monthly_data_area = v200_monthly_data.sel(lat=lat_range, lon=lon_range)
    v200_monthly_data_area = remove_seasonal_cycle_mean(
        v200_monthly_data_area, dim=time_dim, time_range=time_range
    )

    # EOF
    v_EOF_model = get_EOF_model(
        v200_monthly_data_area,
        lat_dim=lat_dim,
        lon_dim=lon_dim,
        time_dim=time_dim,
        n_modes=1,
        random_state=random_state,
        solver=solver,
        solver_kwargs=solver_kwargs,
    )
    v_EOF_result = calc_EOF_analysis(v_EOF_model)
    index_SRP = v_EOF_result["PC"].sel(mode=1)

    # Normalized
    index_normalized_std = index_SRP.sel({time_dim: time_range}).std(dim=time_dim).data
    result = (index_SRP / index_normalized_std).drop_vars(["month", "mode"])
    result.name = "SRPI"
    return result


def calc_index_SRP_EOF1_Kosaka_2009(
    v200_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    lat_range: slice = slice(20, 60),
    lon_range: slice = slice(30, 130),
    time_dim: str = "time",
    random_state=None,
    solver="auto",
    solver_kwargs={},
) -> xr.DataArray:
    """
    The calculation of monthly mean SRP index using empirical orthogonal functions (EOFs) method based on Kosaka et al. (2009):

        EOF1 of V200 over (:math:`\\mathrm{20 ^{\\circ}N - 60 ^{\\circ}N; 30 ^{\\circ} - 130 ^{\\circ}E}`).

    Parameters
    ----------
    v200_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The monthly 200-hPa meridional wind.
    time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
        The time range of seasonal cycle means to be calculated. The default value is the entire time range.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    lat_range: :py:class:`slice <slice>`, default: `slice(20, 60)`.
        The latitude range of computation using EOFs over the region. The default value is from :math:`\\mathrm{20^{\\circ}N}` to :math:`\\mathrm{60^{\\circ}N}`.
    lon_range: :py:class:`slice <slice>`, default: `slice(30, 130)`.
        The longitude range of computation using EOFs over the region. The default value is from :math:`\\mathrm{30^{\\circ}}` to :math:`\\mathrm{130^{\\circ}E}`.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.
    random_state: :py:class:`int<int>`, default `None`.
        Seed for the random number generator.
    solver: {"auto", "full", "randomized"}, default: "auto".
        Solver to use for the EOFs computation.
    solver_kwargs: :py:class:`dict<dict>`, default `{}`.
        Additional keyword arguments to be passed to the EOFs solver.

    Returns
    -------
    The monthly mean SRP index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - `Kosaka, Y., Nakamura, H., Watanabe, M., & Kimoto, M. (2009). Analysis on the dynamics of a wave-like teleconnection pattern along the summertime Asian jet based on a reanalysis dataset and climate model simulations. Journal of the Meteorological Society of Japan. Ser. II, 87(3), 561-580. <https://doi.org/10.2151/jmsj.87.561>`__
    - `Liu Y. Relationship between the Silk Road and Circumglobal Teleconnection Patterns on the Interannual and Interdecadal Timescales. Atmosphere. 2023; 14(11):1626. <https://doi.org/10.3390/atmos14111626>`__

    .. seealso::
        :py:func:`get_EOF_model <easyclimate.core.eof.get_EOF_model>`
    """
    v200_monthly_data = sort_ascending_latlon_coordinates(
        v200_monthly_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    # anomaly
    v200_monthly_data_area = v200_monthly_data.sel(lat=lat_range, lon=lon_range)
    v200_monthly_data_area = remove_seasonal_cycle_mean(
        v200_monthly_data_area, dim=time_dim, time_range=time_range
    )

    # EOF
    v_EOF_model = get_EOF_model(
        v200_monthly_data_area,
        lat_dim=lat_dim,
        lon_dim=lon_dim,
        time_dim=time_dim,
        n_modes=1,
        random_state=random_state,
        solver=solver,
        solver_kwargs=solver_kwargs,
    )
    v_EOF_result = calc_EOF_analysis(v_EOF_model)
    index_SRP = v_EOF_result["PC"].sel(mode=1)

    # Normalized
    index_normalized_std = index_SRP.sel({time_dim: time_range}).std(dim=time_dim).data
    result = (index_SRP / index_normalized_std).drop_vars(["month", "mode"])
    result.name = "SRPI"
    return result


def calc_index_SRP_EOF1_Chen_Huang_2012(
    v200_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    lat_range: slice = slice(30, 60),
    lon_range: slice = slice(30, 130),
    time_dim: str = "time",
    random_state=None,
    solver="auto",
    solver_kwargs={},
) -> xr.DataArray:
    """
    The calculation of monthly mean SRP index using empirical orthogonal functions (EOFs) method based on Chen and Huang (2009):

        EOF1 of V200 over (:math:`\\mathrm{30 ^{\\circ}N - 60 ^{\\circ}N; 30 ^{\\circ} - 130 ^{\\circ}E}`).

    Parameters
    ----------
    v200_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The monthly 200-hPa meridional wind.
    time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
        The time range of seasonal cycle means to be calculated. The default value is the entire time range.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    lat_range: :py:class:`slice <slice>`, default: `slice(20, 60)`.
        The latitude range of computation using EOFs over the region. The default value is from :math:`\\mathrm{30^{\\circ}N}` to :math:`\\mathrm{60^{\\circ}N}`.
    lon_range: :py:class:`slice <slice>`, default: `slice(30, 130)`.
        The longitude range of computation using EOFs over the region. The default value is from :math:`\\mathrm{30^{\\circ}}` to :math:`\\mathrm{130^{\\circ}E}`.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.
    random_state: :py:class:`int<int>`, default `None`.
        Seed for the random number generator.
    solver: {"auto", "full", "randomized"}, default: "auto".
        Solver to use for the EOFs computation.
    solver_kwargs: :py:class:`dict<dict>`, default `{}`.
        Additional keyword arguments to be passed to the EOFs solver.

    Returns
    -------
    The monthly mean SRP index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - `Chen, G., & Huang, R. (2012). Excitation Mechanisms of the Teleconnection Patterns Affecting the July Precipitation in Northwest China. Journal of Climate, 25(22), 7834-7851. <https://doi.org/10.1175/JCLI-D-11-00684.1>`__
    - `Liu Y. Relationship between the Silk Road and Circumglobal Teleconnection Patterns on the Interannual and Interdecadal Timescales. Atmosphere. 2023; 14(11):1626. <https://doi.org/10.3390/atmos14111626>`__

    .. seealso::
        :py:func:`get_EOF_model <easyclimate.core.eof.get_EOF_model>`
    """
    v200_monthly_data = sort_ascending_latlon_coordinates(
        v200_monthly_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    # anomaly
    v200_monthly_data_area = v200_monthly_data.sel(lat=lat_range, lon=lon_range)
    v200_monthly_data_area = remove_seasonal_cycle_mean(
        v200_monthly_data_area, dim=time_dim, time_range=time_range
    )

    # EOF
    v_EOF_model = get_EOF_model(
        v200_monthly_data_area,
        lat_dim=lat_dim,
        lon_dim=lon_dim,
        time_dim=time_dim,
        n_modes=1,
        random_state=random_state,
        solver=solver,
        solver_kwargs=solver_kwargs,
    )
    v_EOF_result = calc_EOF_analysis(v_EOF_model)
    index_SRP = v_EOF_result["PC"].sel(mode=1)

    # Normalized
    index_normalized_std = index_SRP.sel({time_dim: time_range}).std(dim=time_dim).data
    result = (index_SRP / index_normalized_std).drop_vars(["month", "mode"])
    result.name = "SRPI"
    return result


def calc_index_SRP_EOF1_Sato_Takahashi_2006(
    v200_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    lat_range: slice = slice(30, 60),
    lon_range: slice = slice(80, 200),
    time_dim: str = "time",
    random_state=None,
    solver="auto",
    solver_kwargs={},
) -> xr.DataArray:
    """
    The calculation of monthly mean SRP index using empirical orthogonal functions (EOFs) method based on Sato and Takahashi (2006):

        EOF1 of V200 over (:math:`\\mathrm{30 ^{\\circ}N - 60 ^{\\circ}N; 80 ^{\\circ}E - 160 ^{\\circ}W}`).

    Parameters
    ----------
    v200_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The monthly 200-hPa meridional wind.
    time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
        The time range of seasonal cycle means to be calculated. The default value is the entire time range.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    lat_range: :py:class:`slice <slice>`, default: `slice(30, 60)`.
        The latitude range of computation using EOFs over the region. The default value is from :math:`\\mathrm{30^{\\circ}N}` to :math:`\\mathrm{60^{\\circ}N}`.
    lon_range: :py:class:`slice <slice>`, default: `slice(80, 200)`.
        The longitude range of computation using EOFs over the region. The default value is from :math:`\\mathrm{80^{\\circ}E}` to :math:`\\mathrm{160^{\\circ}W}`.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.
    random_state: :py:class:`int<int>`, default `None`.
        Seed for the random number generator.
    solver: {"auto", "full", "randomized"}, default: "auto".
        Solver to use for the EOFs computation.
    solver_kwargs: :py:class:`dict<dict>`, default `{}`.
        Additional keyword arguments to be passed to the EOFs solver.

    Returns
    -------
    The monthly mean SRP index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - `Chen, G., & Huang, R. (2012). Excitation Mechanisms of the Teleconnection Patterns Affecting the July Precipitation in Northwest China. Journal of Climate, 25(22), 7834-7851. <https://doi.org/10.1175/JCLI-D-11-00684.1>`__
    - `Liu Y. Relationship between the Silk Road and Circumglobal Teleconnection Patterns on the Interannual and Interdecadal Timescales. Atmosphere. 2023; 14(11):1626. <https://doi.org/10.3390/atmos14111626>`__

    .. seealso::
        :py:func:`get_EOF_model <easyclimate.core.eof.get_EOF_model>`
    """
    v200_monthly_data = sort_ascending_latlon_coordinates(
        v200_monthly_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    # anomaly
    v200_monthly_data_area = v200_monthly_data.sel(lat=lat_range, lon=lon_range)
    v200_monthly_data_area = remove_seasonal_cycle_mean(
        v200_monthly_data_area, dim=time_dim, time_range=time_range
    )

    # EOF
    v_EOF_model = get_EOF_model(
        v200_monthly_data_area,
        lat_dim=lat_dim,
        lon_dim=lon_dim,
        time_dim=time_dim,
        n_modes=1,
        random_state=random_state,
        solver=solver,
        solver_kwargs=solver_kwargs,
    )
    v_EOF_result = calc_EOF_analysis(v_EOF_model)
    index_SRP = v_EOF_result["PC"].sel(mode=1)

    # Normalized
    index_normalized_std = index_SRP.sel({time_dim: time_range}).std(dim=time_dim).data
    result = (index_SRP / index_normalized_std).drop_vars(["month", "mode"])
    result.name = "SRPI"
    return result


def calc_index_SRP_1point_Lu_2002(
    v200_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    time_dim: str = "time",
) -> xr.DataArray:
    """
    The calculation of monthly mean SRP index based on Lu et al. (2002):

        V200 at :math:`\\mathrm{42.5 ^{\\circ}N, 105 ^{\\circ}E}`.

    Parameters
    ----------
    v200_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The monthly 200-hPa meridional wind.
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
    The monthly mean SRP index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - `Lu, R.-Y., Oh, J.-H., & Kim, B.-J. (2002). A teleconnection pattern in upper-level meridional wind over the North African and Eurasian continent in summer. Tellus A: Dynamic Meteorology and Oceanography, 54(1), 44-55. <https://doi.org/10.3402/tellusa.v54i1.12122>`__
    - `Liu Y. Relationship between the Silk Road and Circumglobal Teleconnection Patterns on the Interannual and Interdecadal Timescales. Atmosphere. 2023; 14(11):1626. <https://doi.org/10.3390/atmos14111626>`__
    """
    v200_monthly_data = sort_ascending_latlon_coordinates(
        v200_monthly_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    # anomaly
    v200_anomaly_data = remove_seasonal_cycle_mean(
        v200_monthly_data, dim=time_dim, time_range=time_range
    )

    # V200 (42.5°N, 105°E)
    index_SRP = v200_anomaly_data.sel(lat=42.5, lon=105, method="nearest")

    # Normalized
    index_normalized_std = index_SRP.sel({time_dim: time_range}).std(dim=time_dim).data
    result = (index_SRP / index_normalized_std).drop_vars("month")
    result.name = "SRPI"
    return result
