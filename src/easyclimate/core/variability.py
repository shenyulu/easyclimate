"""
This module calculate climate variability
"""

import numpy as np
import xarray as xr
from typing import Literal
from scipy.fft import rfft, irfft
from .extract import get_specific_months_data

__all__ = [
    "calc_all_climatological_mean",
    "calc_seasonal_climatological_mean",
    "calc_seasonal_cycle_mean",
    "calc_seasonal_cycle_std",
    "calc_seasonal_cycle_var",
    "calc_seasonal_mean",
    "remove_seasonal_cycle_mean",
    "calc_monthly_climatological_std_without_seasonal_cycle_mean",
    "calc_monthly_climatological_var_without_seasonal_cycle_mean",
    "smooth_daily_annual_cycle",
    "calc_daily_annual_cycle_mean",
    "calc_daily_annual_cycle_std",
    "calc_daily_annual_cycle_var",
    "remove_smooth_daily_annual_cycle_mean",
    "calc_horizontal_wind_components_std",
    "calc_windspeed_dataset",
    "calc_windspeed_dataarray",
    "populate_monmean2everymon",
    "populate_daymean2everyday",
    "calc_daily_climatological_anomaly",
    "remove_low_frequency_signal",
]


def calc_all_climatological_mean(
    data_input: xr.DataArray | xr.Dataset, dim: str = "time", **kwargs
) -> xr.DataArray:
    """
    Calculation of the climatological mean over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
    dim: :py:class:`str <str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return data_input.mean(dim=dim, **kwargs)


def calc_seasonal_climatological_mean(
    data_input: xr.DataArray | xr.Dataset, dim: str = "time", **kwargs
) -> xr.DataArray:
    """
    Calculation of the seasonal climatological mean over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
    dim: :py:class:`str <str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return data_input.groupby(data_input[dim].dt.season).mean(dim=dim, **kwargs)


def calc_seasonal_cycle_mean(
    data_input: xr.DataArray | xr.Dataset, dim: str = "time", **kwargs
) -> xr.DataArray:
    """
    Calculation of the seasonal cycle means over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

    .. caution:: `data_input` must be **monthly** data.

    dim: :py:class:`str <str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return data_input.groupby(data_input[dim].dt.month).mean(dim=dim, **kwargs)


def calc_seasonal_cycle_std(
    data_input: xr.DataArray | xr.Dataset, dim: str = "time", **kwargs
) -> xr.DataArray:
    """
    Calculation of the seasonal cycle standard deviation over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

    .. caution:: `data_input` must be **monthly** data.

    dim: :py:class:`str <str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating standard deviation on this object's data.
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return data_input.groupby(data_input[dim].dt.month).std(dim=dim, **kwargs)


def calc_seasonal_cycle_var(
    data_input: xr.DataArray | xr.Dataset, dim: str = "time", **kwargs
) -> xr.DataArray:
    """
    Calculation of the seasonal cycle standard deviation over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

    .. caution:: `data_input` must be **monthly** data.

    dim: :py:class:`str <str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating variance on this object's data.
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return data_input.groupby(data_input[dim].dt.month).var(dim=dim, **kwargs)


def calc_seasonal_mean(
    data_input: xr.DataArray | xr.Dataset,
    dim: str = "time",
    extract_season: Literal["DJF", "MAM", "JJA", "SON", None] = None,
    **kwargs,
) -> xr.DataArray:
    """
    Calculation of the seasonal means per year over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

    .. caution:: `data_input` must be **monthly** data.

    dim: :py:class:`str <str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    extract_season: :py:class:`list <list>`, e.g., one or multiple items from `['DJF', 'MAM', 'JJA', 'SON']`. default: None.
        Extraction seasons. A variety of seasons can be placed in it.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray <xarray.DataArray>`.

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_ao_index.py
        ./dynamic_docs/plot_oceanic_front.py
        ./dynamic_docs/plot_multi_linear_reg.py
    """
    result_seasonal_mean = data_input.resample({dim: "QS-DEC"}).mean(dim=dim)

    extract_season_list = []
    if extract_season is not None:
        if "DJF" in extract_season:
            extract_season_list.append(12)
        if "MAM" in extract_season:
            extract_season_list.append(3)
        if "JJA" in extract_season:
            extract_season_list.append(6)
        if "SON" in extract_season:
            extract_season_list.append(9)
        return get_specific_months_data(result_seasonal_mean, extract_season_list)
    else:
        return result_seasonal_mean


def remove_seasonal_cycle_mean(
    data_input: xr.DataArray | xr.Dataset,
    dim: str = "time",
    time_range: slice = slice(None, None),
) -> xr.DataArray:
    """
    Remove of the seasonal cycle means over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
        The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

    .. caution:: `data_input` must be **monthly** data.

    dim: :py:class:`str <str>`.
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

    time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
        The time range of seasonal cycle means to be calculated. The default value is the entire time range.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_ao_index.py
        ./dynamic_docs/plot_basic_statistical_analysis.py
        ./dynamic_docs/plot_da_bbo.py
        ./dynamic_docs/plot_multieof.py
        ./dynamic_docs/plot_ocean_mix_layer.py
        ./dynamic_docs/plot_time_scale_average.py
        ./dynamic_docs/plot_corr_reg.py
    """
    gb = data_input.groupby(data_input[dim].dt.month)
    data_input_mean = data_input.sel({dim: time_range})
    gb_mean = data_input_mean.groupby(data_input_mean[dim].dt.month)
    return gb - gb_mean.mean(dim=dim)


def calc_monthly_climatological_std_without_seasonal_cycle_mean(
    data_input: xr.DataArray | xr.Dataset, dim: str = "time", **kwargs
) -> xr.DataArray:
    """
    Calculate the standard deviation of monthly data anomalies over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

    .. caution:: `data_input` must be **monthly** data.

    dim: :py:class:`str <str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating standard deviation on this object's data.
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return remove_seasonal_cycle_mean(data_input, dim=dim).std(dim=dim, **kwargs)


def calc_monthly_climatological_var_without_seasonal_cycle_mean(
    data_input: xr.DataArray | xr.Dataset, dim: str = "time", **kwargs
) -> xr.DataArray:
    """
    Calculate the variance of monthly data anomalies over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The data of :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>` to be calculated.

    .. caution:: `data_input` must be **monthly** data.

    dim: :py:class:`str <str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating variance on this object's data.
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return remove_seasonal_cycle_mean(data_input, dim=dim).var(dim=dim, **kwargs)


def smooth_daily_annual_cycle(
    daily_annual_cycle_data: xr.DataArray,
    harmonics_num: int = 3,
    time_dim: str = "dayofyear",
) -> xr.DataArray:
    """
    Calculates a smooth mean daily annual cycle for an array nominally dimensioned.

    Parameters
    ----------
    daily_annual_cycle_data : :py:class:`xarray.DataArray<xarray.DataArray>`
        The input data array with time as the first dimension. The time dimension should be named as specified by `time_dim`.
    harmonics_num : int, optional
        The number of harmonics to retain in the FFT. Default is 3.
    time_dim : str, optional
        The name of the time dimension in the DataArray. Default is "dayofyear".

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`

        The smoothed daily cycle data.

    .. seealso::

        - https://www.ncl.ucar.edu/Document/Functions/Contributed/smthClmDayTLL.shtml

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_smooth_daily_cycle.py
    """
    # move time dimension to 1st
    dims_order = daily_annual_cycle_data.dims
    daily_annual_cycle_data = daily_annual_cycle_data.transpose(time_dim, ...)

    # get time dimension size
    nt = daily_annual_cycle_data[time_dim].shape[0]

    cf = rfft(
        daily_annual_cycle_data.values, axis=0
    )  # xarray.DataArray -> numpy.ndarray
    cf[harmonics_num, :, :] = 0.5 * cf[harmonics_num, :, :]  # mini-taper.
    cf[harmonics_num + 1 :, :, :] = 0.0  # set all higher coef to 0.0
    icf = irfft(cf, n=nt, axis=0)  # reconstructed series

    # create `xarray.DataArray` and transpose original dimensions order
    clmDaySmth = daily_annual_cycle_data.copy(data=icf, deep=True)
    clmDaySmth = clmDaySmth.transpose(*dims_order)
    return clmDaySmth


def calc_daily_annual_cycle_mean(
    data_input: xr.DataArray | xr.Dataset,
    dim: str = "time",
    **kwargs,
) -> xr.DataArray:
    """
    Calculation of the daily means per year over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

    .. caution::

        `data_input` must be **daily** or **hourly** data.
        At least one year of time range must be included in the `data_input`.

    dim: :py:class:`str <str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed to the mean function.

    Returns
    -------
    :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`.

    .. caution::

        - For complete coverage, the data should span at least one full year.
        - If the data is sub-daily (e.g., hourly), the mean is taken over all sub-daily time points for each day of the year.
    """
    result_daily_cycle_mean = data_input.groupby(data_input[dim].dt.dayofyear).mean(
        dim=dim, **kwargs
    )

    return result_daily_cycle_mean


def calc_daily_annual_cycle_std(
    data_input: xr.DataArray | xr.Dataset,
    dim: str = "time",
    **kwargs,
) -> xr.DataArray:
    """
    Calculation of the daily standard deviation per year over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

    .. caution::

        `data_input` must be **daily** or **hourly** data.
        At least one year of time range must be included in the `data_input`.

    dim: :py:class:`str <str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed to the standard deviation function.

    Returns
    -------
    :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`.

    .. caution::

        - For complete coverage, the data should span at least one full year.
        - If the data is sub-daily (e.g., hourly), the std is taken over all sub-daily time points for each day of the year.
    """
    result_daily_cycle_std = data_input.groupby(data_input[dim].dt.dayofyear).std(
        dim=dim, **kwargs
    )

    return result_daily_cycle_std


def calc_daily_annual_cycle_var(
    data_input: xr.DataArray | xr.Dataset,
    dim: str = "time",
    **kwargs,
) -> xr.DataArray:
    """
    Calculation of the daily variance per year over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

    .. caution::

        `data_input` must be **daily** or **hourly** data.
        At least one year of time range must be included in the `data_input`.

    dim: :py:class:`str <str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed to the variance function.

    Returns
    -------
    :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`.

    .. caution::

        - For complete coverage, the data should span at least one full year.
        - If the data is sub-daily (e.g., hourly), the var is taken over all sub-daily time points for each day of the year.
    """
    result_daily_cycle_var = data_input.groupby(data_input[dim].dt.dayofyear).var(
        dim=dim, **kwargs
    )

    return result_daily_cycle_var


def remove_smooth_daily_annual_cycle_mean(
    data_input: xr.DataArray,
    daily_cycle_mean_time_range: slice = slice(None, None),
    extract_time_range: slice = slice(None, None),
    harmonics_num: int = 3,
    dim: str = "time",
):
    """
    Removes the smooth daily annual cycle mean from the input data.

    This function first calculates the daily cycle mean over a specified time range,
    smooths that mean using a specified number of harmonics, and then subtracts this
    smoothed cycle from the input data over another specified time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        The input data from which to remove the smooth daily annual cycle mean.
    daily_cycle_mean_time_range : slice, optional
        The time range used to compute the daily annual cycle mean. Default is all time.
    extract_time_range : slice, optional
        The time range from which to extract the data and remove the daily annual cycle. Default is all time.
    harmonics_num : int, optional
        The number of harmonics to use in smoothing the daily annual cycle mean. Default is 3.
    dim : str, optional
        The name of the time dimension. Default is "time".

    Returns
    -------
    :py:class:`xarray.DataArray <xarray.DataArray>`
        The input data with the smooth daily cycle mean removed.
    """
    day_clm = calc_daily_annual_cycle_mean(
        data_input.sel({dim: daily_cycle_mean_time_range}), dim=dim
    )
    day_clm = day_clm.compute()
    day_clm_sm = smooth_daily_annual_cycle(day_clm, harmonics_num, "dayofyear")
    extract_input = data_input.sel({dim: extract_time_range})
    smoothed_data = extract_input.groupby(extract_input[dim].dt.dayofyear) - day_clm_sm
    return smoothed_data


def calc_horizontal_wind_components_std(
    uv_dataset: xr.Dataset, u_dim="u", v_dim="v", time_dim="time", ddof=0
) -> xr.Dataset:
    """
    Calculate the standard deviation of vector wind speed and direction.

    The standard deviation of vector wind speed

    .. math::
        \\sigma_s = [U^2 \\sigma_u^2 + V^2 \\sigma_v^2 + 2 U V \\sigma_{uv}]^{1/2} S^{-1},

    The standard deviation of vector wind direction

    .. math::
        \\sigma_d = [V^2 \\sigma_u^2 + U^2 \\sigma_v^2 + 2 U V \\sigma_{uv}]^{1/2} S^{-2},

    Where time mean of :math:`u` is :math:`U = n^{-1} \\sum u_i`, time mean of :math:`v` is :math:`V = n^{-1} \\sum v_i`,
    time variance of :math:`u` is :math:`\\sigma_u^2 = n^{-1} \\sum u_{i}^{2} - U^2`,
    time variance of :math:`v` is :math:`\\sigma_v^2 = n^{-1} \\sum v_{i}^{2} - V^2`,
    time covariance of :math:`u`, :math:`v` is :math:`\\sigma_{uv} = n^{-1} \\sum u_i v_i - UV`,
    vector mean wind speed is :math:`S = (U^2 + V^2)^{1/2}`.

    Parameters
    ----------
    uv_dataset : :py:class:`xarray.Dataset<xarray.Dataset>`
        :py:class:`xarray.Dataset<xarray.Dataset>` data containing zonal and meridional wind components.
    u_dim: :py:class:`str <str>`, default: `u`
        Variable name for the u velocity (in x direction).
    v_dim: :py:class:`str <str>`, default: `v`
        Variable name for the v velocity (in y direction).
    time_dim : :py:class:`str <str>`, default: `time`
        Dimension(s) over which to apply. By default is applied over the `time` dimension.
    ddof : :py:class:`int <int>`, default: 1
        If `ddof=1`, covariance is normalized by `N-1`, giving an unbiased estimate, else normalization is by `N`.

    Returns
    -------
    :py:class:`xarray.Dataset<xarray.Dataset>`
        - sigma_s: the standard deviation of vector wind speed.
        - sigma_d: the standard deviation of vector wind direction.

    Reference
    --------------
    - G. R. Ackermann. (1983). Means and Standard Deviations of Horizontal Wind Components. https://doi.org/10.1175/1520-0450(1983)022%3C0959:MASDOH%3E2.0.CO;2
    """
    U = uv_dataset[u_dim].mean(dim=time_dim)
    V = uv_dataset[v_dim].mean(dim=time_dim)
    S = np.hypot(U, V)
    # D = np.arctan(U / V)

    sigma2_u = uv_dataset[u_dim].var(dim=time_dim)
    sigma2_v = uv_dataset[v_dim].var(dim=time_dim)
    sigma_uv = xr.cov(uv_dataset[u_dim], uv_dataset[v_dim], dim=time_dim, ddof=ddof)

    sigma_s = (U**2 * sigma2_u + V**2 * sigma2_v + 2 * U * V * sigma_uv) ** (
        1 / 2
    ) * S ** (-1)
    sigma_d = (V**2 * sigma2_u + U**2 * sigma2_v - 2 * U * V * sigma_uv) ** (
        1 / 2
    ) * S ** (-2)

    return uv_dataset.assign({"sigma_s": sigma_s, "sigma_d": sigma_d})


def calc_windspeed_dataset(
    uv_dataset: xr.Dataset, u_dim: str = "u", v_dim: str = "v"
) -> xr.Dataset:
    """
    Calculate the horizontal wind speed from zonal and meridional wind components in a :py:class:`xarray.Dataset<xarray.Dataset>`.

    The wind speed is computed as the magnitude of the horizontal wind vector:

    .. math::

        S = \\sqrt{u^2 + v^2},

    where :math:`u` is the zonal wind component and :math:`v` is the meridional wind component.

    Parameters
    ----------
    uv_dataset : :py:class:`xarray.Dataset<xarray.Dataset>`
        :py:class:`xarray.Dataset<xarray.Dataset>` containing zonal and meridional wind components.
    u_dim : :py:class:`str <str>`, default: `u`
        Variable name for the zonal wind component (in x direction).
    v_dim : :py:class:`str <str>`, default: `v`
        Variable name for the meridional wind component (in y direction).

    Returns
    -------
    :py:class:`xarray.Dataset<xarray.Dataset>`
        A copy of the input dataset with an additional variable `speed` containing the wind speed.

    Examples
    --------
    >>> ds = xr.Dataset({"u": (("time",), [1, 2, 3]), "v": (("time",), [4, 5, 6])})
    >>> ds_with_speed = calc_windspeed_dataset(ds, u_dim="u", v_dim="v")
    >>> print(ds_with_speed["speed"])
    <xarray.DataArray 'speed' (time: 3)> Size: 24B
    array([4.12310563, 5.38516481, 6.70820393])
    Dimensions without coordinates: time
    """
    ds_ = uv_dataset.copy(deep=True)
    ds_["speed"] = np.sqrt(uv_dataset[u_dim] ** 2 + uv_dataset[v_dim] ** 2)
    return ds_


def calc_windspeed_dataarray(
    u_data: xr.DataArray, v_data: xr.DataArray
) -> xr.DataArray:
    """
    Calculate the horizontal wind speed from zonal and meridional wind components in :py:class:`xarray.DataArray<xarray.DataArray>`.

    The wind speed is computed as the magnitude of the horizontal wind vector:

    .. math::

        S = \\sqrt{u^2 + v^2},

    where :math:`u` is the zonal wind component and :math:`v` is the meridional wind component.

    Parameters
    ----------
    u_data : :py:class:`xarray.DataArray<xarray.DataArray>`
        :py:class:`xarray.DataArray<xarray.DataArray>` containing the zonal wind component (in x direction).
    v_data : :py:class:`xarray.DataArray<xarray.DataArray>`
        :py:class:`xarray.DataArray<xarray.DataArray>` containing the meridional wind component (in y direction).

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`
        A :py:class:`xarray.DataArray<xarray.DataArray>` containing the wind speed.

    Examples
    --------
    >>> u = xr.DataArray([1, 2, 3], dims="time")
    >>> v = xr.DataArray([4, 5, 6], dims="time")
    >>> speed = calc_windspeed_dataarray(u, v)
    >>> print(speed)
    <xarray.DataArray (time: 3)> Size: 24B
    array([4.12310563, 5.38516481, 6.70820393])
    Dimensions without coordinates: time
    """
    ds_ = np.sqrt(u_data**2 + v_data**2)
    ds_.attrs = dict()
    return ds_


def populate_monmean2everymon(
    data_monthly: xr.DataArray,
    data_climatology_monthly_data: xr.DataArray = None,
    time_dim: str = "time",
) -> xr.DataArray:
    """
    Populate the data of each month using the monthly mean state of the `data_monthly` or given dataset.

    Parameters
    ----------
    data_monthly: :py:class:`xarray.DataArray<xarray.DataArray>`.
        :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
    data_climatology_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`, default `None`.
        The monthly climatology dataset. If it is `None`, the climatology is derived from `data_monthly`.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    if data_climatology_monthly_data is None:
        time_step_all = data_monthly[time_dim].shape[0]
        month_int = data_monthly.time.dt.month
        month_climate = data_monthly.groupby(data_monthly[time_dim].dt.month).mean(
            dim=time_dim
        )
        data_monthly_empty = xr.full_like(data_monthly, fill_value=np.nan)

        for time_step in np.arange(0, time_step_all):
            time_step_month = month_int.isel({time_dim: time_step}).data
            data_monthly_empty[{time_dim: time_step}] = month_climate.sel(
                month=time_step_month
            )

        return data_monthly_empty

    else:
        result_data = xr.full_like(data_monthly, fill_value=np.nan)
        time_length = result_data[time_dim].shape[0]
        time_month = result_data[time_dim].dt.month.data
        climate_data_month_index = data_climatology_monthly_data[time_dim].dt.month.data

        for time_item in np.arange(time_length):
            # Target month index
            time_month_item = time_month[time_item]
            # Correspond month index in the climate data
            month_index = np.transpose(
                np.nonzero(climate_data_month_index == time_month_item)
            )
            month_index = month_index.item()

            result_data[{time_dim: time_item}] = data_climatology_monthly_data.isel(
                {time_dim: month_index}
            )

        return result_data


def populate_daymean2everyday(
    data_daily: xr.DataArray,
    data_climatology_daily_data: xr.DataArray = None,
    time_dim: str = "time",
) -> xr.DataArray:
    """
    Populate the data of each day using the daily mean state of the `data_daily` or given dataset.

    Parameters
    ----------
    - data_daily: :py:class:`xarray.DataArray<xarray.DataArray>`.
        :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
    - data_climatology_daily_data: :py:class:`xarray.DataArray<xarray.DataArray>`, default `None`.
        The daily climatology dataset. If it is `None`, the climatology is derived from `data_monthly`.
    - time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    if data_climatology_daily_data is None:
        time_step_all = data_daily[time_dim].shape[0]
        day_int = data_daily.time.dt.day
        day_climate = data_daily.groupby(data_daily[time_dim].dt.day).mean(dim=time_dim)
        data_daily_empty = xr.full_like(data_daily, fill_value=np.nan)

        for time_step in np.arange(0, time_step_all):
            time_step_day = day_int.isel({time_dim: time_step}).data
            data_daily_empty[{time_dim: time_step}] = day_climate.sel(day=time_step_day)

        return data_daily_empty

    else:
        result_data = xr.full_like(data_daily, fill_value=np.nan)
        time_length = result_data[time_dim].shape[0]
        time_dayofyear = result_data[time_dim].dt.dayofyear.data
        climate_data_dayofyear_index = data_climatology_daily_data[
            time_dim
        ].dt.dayofyear.data

        for time_item in np.arange(time_length):
            # Target dayofyear index
            time_dayofyear_item = time_dayofyear[time_item]
            # Correspond dayofyear index in the climate data
            dayofyear_index = np.transpose(
                np.nonzero(climate_data_dayofyear_index == time_dayofyear_item)
            )
            dayofyear_index = dayofyear_index.item()

            result_data[{time_dim: time_item}] = data_climatology_daily_data.isel(
                {time_dim: dayofyear_index}
            )

        return result_data


def calc_daily_climatological_anomaly(
    data_daily: xr.DataArray | xr.Dataset,
    data_climatology_daily_data: xr.DataArray | xr.Dataset,
    timd_dim="time",
) -> xr.DataArray | xr.Dataset:
    """
    Calulate daily anomaly using the given dataset of climatological mean state .

    - data_daily: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
        :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
    - data_climatology_daily_data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
        The daily climatology dataset.
    - time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
    """
    data_daily_anomaly = (
        data_daily.groupby(data_daily[timd_dim].dt.dayofyear)
        - data_climatology_daily_data.groupby(
            data_climatology_daily_data[timd_dim].dt.dayofyear
        ).mean()
    )
    data_daily_anomaly = data_daily_anomaly.drop_vars("dayofyear")
    return data_daily_anomaly


def remove_low_frequency_signal(
    da: xr.DataArray, window: int = 120, center: bool = False, time_dim: str = "time"
) -> xr.DataArray:
    """
    Remove low-frequency signal by subtracting the running mean from a time series.

    This function removes the effect of interannual variability by subtracting the
    running mean of the specified window (default 120 days), as described in Wheeler
    and Hendon (2004). The method is commonly used in the context of the Madden-Julian
    Oscillation (MJO) index calculation for monitoring and prediction.

    Parameters
    ----------
    da : :py:class:`xarray.DataArray<xarray.DataArray>`
        Input time series data array with a time dimension.
    window : :py:class:`int <int>`, optional
        Size of the moving average window in days (default is 120).
    center : :py:class:`bool <bool>`, optional
        If ``True``, the moving average is centered (mean of window around each point).
        If ``False``, the moving average is trailing (mean of last window days).
        Default is ``False``.
    time_dim : :py:class:`str <str>`, optional
        Name of the time dimension in the input DataArray (default is "time").

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`
        The input data array with the low-frequency signal (running mean) removed.

    References
    ----------
    - Wheeler, M. C., & Hendon, H. H. (2004). An All-Season Real-Time Multivariate MJO Index: Development of an Index for Monitoring and Prediction. Monthly Weather Review, 132(8), 1917-1932. https://journals.ametsoc.org/view/journals/mwre/132/8/1520-0493_2004_132_1917_aarmmi_2.0.co_2.xml

    Examples
    --------
    >>> import xarray as xr
    >>> da = xr.DataArray([...], dims=['time'], coords={'time': [...]})
    >>> result = remove_low_frequency_signal(da, window=120, center=False, time_dim='time')
    """
    # Calculate moving average over the specified time dimension
    # Using rolling mean with center=False for trailing window
    rolling_mean = da.rolling({time_dim: window}, center=center).mean()
    result = da - rolling_mean

    return result
