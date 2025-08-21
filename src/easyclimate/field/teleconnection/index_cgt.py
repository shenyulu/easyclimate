"""
Circumglobal teleconnection pattern (CGT) Index
"""

import xarray as xr
from ...core.utility import sort_ascending_latlon_coordinates
from ...core.variability import remove_seasonal_cycle_mean
from ...core.yearstat import calc_yearly_climatological_mean
from ...core.extract import get_specific_months_data
from ...core.stat import calc_detrend_spatial
from ...core.normalized import normalize_zscore
from ...filter.spectrum import filter_fourier_harmonic_analysis
from typing import Literal

__all__ = [
    "calc_index_CGT_1point_Ding_Wang_2005",
    "calc_index_CGT_NH_Ding_Wang_2005",
]


def calc_index_CGT_1point_Ding_Wang_2005(
    z200_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    time_dim: str = "time",
    normalized: bool = True,
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
    normalized: :py:class:`bool <bool>`, default `True`, optional.
        Whether to standardize the index based on standard deviation over `time_range`.

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
    index_CGT = z_anomaly_data.sel(
        {lat_dim: slice(35, 40), lon_dim: slice(60, 70)}
    ).mean(dim=(lat_dim, lon_dim))

    # Normalized
    if normalized == True:
        index_normalized_std = (
            index_CGT.sel({time_dim: time_range}).std(dim=time_dim).data
        )
        result = (index_CGT / index_normalized_std).drop_vars("month")
        return result
    elif normalized == False:
        return index_CGT


def calc_index_CGT_NH_Ding_Wang_2005(
    z200_monthly_data: xr.DataArray,
    output_freq: Literal["monthly", "seasonally"],
    time_range: slice = slice(None, None),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    lat_range: slice = slice(20, 85),
    time_dim: str = "time",
) -> xr.DataArray:
    """
    The calculation of monthly mean circumglobal teleconnection pattern (CGT) index using Ding, Q., & Wang, B. (2005) method

    Parameters
    ----------
    z200_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The monthly 200-hPa geopotential height.
    output_freq: :py:class:`str <str>`.
        The output frequency. Optional values are `monthly`, `seasonally`.
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

    Returns
    -------
    The monthly/seasonally mean CGT index (:py:class:`xarray.DataArray<xarray.DataArray>`).

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
    z_monthly_data_NH = z200_monthly_data.sel({lat_dim: lat_range})
    z_anomaly_data_NH = remove_seasonal_cycle_mean(
        z_monthly_data_NH, dim=time_dim, time_range=time_range
    )
    # JJAS
    z_anomaly_data_NH_JJAS = get_specific_months_data(
        z_anomaly_data_NH, month_array=[6, 7, 8, 9], dim=time_dim
    )

    if output_freq == "monthly":
        # Detrend
        z_anomaly_data_NH_JJAS_detrend = calc_detrend_spatial(
            z_anomaly_data_NH_JJAS, time_dim=time_dim
        )
        # Fourier Harmonic Analysis
        z_anomaly_data_NH_JJAS_detrend_filtered = filter_fourier_harmonic_analysis(
            da=z_anomaly_data_NH_JJAS_detrend,
            time_dim=time_dim,
            period_bounds=(None, 8),
            filter_type="highpass",
            sampling_interval=1 / 4,
            apply_window=True,
        )
    elif output_freq == "seasonally":
        # Detrend
        z_anomaly_data_NH_JJAS_yearly = calc_yearly_climatological_mean(
            z_anomaly_data_NH_JJAS, dim=time_dim
        )
        z_anomaly_data_NH_JJAS_yearly_detrend = calc_detrend_spatial(
            z_anomaly_data_NH_JJAS_yearly, time_dim=time_dim
        )
        # Fourier Harmonic Analysis
        z_anomaly_data_NH_JJAS_detrend_filtered = filter_fourier_harmonic_analysis(
            da=z_anomaly_data_NH_JJAS_yearly_detrend,
            time_dim=time_dim,
            period_bounds=(None, 8),
            filter_type="highpass",
            sampling_interval=1,
            apply_window=True,
        )
    else:
        raise ValueError(
            f"`output_freq` should be `monthly` or `seasonally`, not {output_freq}"
        )

    index_cgt = normalize_zscore(
        z_anomaly_data_NH_JJAS_detrend_filtered.sel(
            lon=slice(60, 70), lat=slice(35, 40)
        ).mean(dim=("lat", "lon"))
    )

    index_cgt.attrs = dict()
    index_cgt.attrs["output_freq"] = output_freq
    return index_cgt
