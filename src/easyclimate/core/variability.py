"""
This module calculate climate variability
"""
import numpy as np
import xarray as xr

__all__ =  ["calc_all_climatological_mean", "calc_seasonal_climatological_mean", 
            "calc_seasonal_cycle_mean", "calc_seasonal_cycle_std", "calc_seasonal_cycle_var",
            "remove_seasonal_cycle_mean", "calc_monthly_climatological_std_without_seasonal_cycle_mean",
            "calc_monthly_climatological_var_without_seasonal_cycle_mean",
            "calc_horizontal_wind_components_std", "transfer_monmean2everymonthmean",
            "mapping_daily_climatological_mean2every_day", "mapping_monthly_climatological_mean2every_month",
            "calc_daily_climatological_anomaly"]

def calc_all_climatological_mean(
    data_input: xr.DataArray | xr.Dataset, 
    dim = 'time', 
    **kwargs
) -> xr.DataArray:
    """
    Calculation of the climatological mean over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
    dim : :py:class:`str<python.str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data. 
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return data_input.mean(dim = dim, **kwargs)

def calc_seasonal_climatological_mean(
    data_input: xr.DataArray | xr.Dataset,
    dim = 'time',
    **kwargs
) -> xr.DataArray:
    """
    Calculation of the seasonal climatological mean over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
    dim : :py:class:`str<python.str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data. 
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return data_input.groupby(data_input[dim].dt.season).mean(dim = dim, **kwargs)

def calc_seasonal_cycle_mean(
    data_input: xr.DataArray | xr.Dataset,
    dim = 'time',
    **kwargs
) -> xr.DataArray:
    """
    Calculation of the seasonal cycle means over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

    .. caution:: `data_input` must be **monthly** data.

    dim : :py:class:`str<python.str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data. 
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return data_input.groupby(data_input[dim].dt.month).mean(dim = dim, **kwargs)

def calc_seasonal_cycle_std(
    data_input: xr.DataArray | xr.Dataset,
    dim = 'time',
    **kwargs
) -> xr.DataArray:
    """
    Calculation of the seasonal cycle standard deviation over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

    .. caution:: `data_input` must be **monthly** data.

    dim : :py:class:`str<python.str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating standard deviation on this object's data. 
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return data_input.groupby(data_input[dim].dt.month).std(dim = dim, **kwargs)

def calc_seasonal_cycle_var(
    data_input: xr.DataArray | xr.Dataset,
    dim = 'time',
    **kwargs
) -> xr.DataArray:
    """
    Calculation of the seasonal cycle standard deviation over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

    .. caution:: `data_input` must be **monthly** data.

    dim : :py:class:`str<python.str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating variance on this object's data. 
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return data_input.groupby(data_input[dim].dt.month).var(dim = dim, **kwargs)

def remove_seasonal_cycle_mean(
    data_input: xr.DataArray | xr.Dataset, 
    dim = 'time', 
    **kwargs
) -> xr.DataArray:
    """
    Remove of the seasonal cycle means over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

    .. caution:: `data_input` must be **monthly** data.

    dim : :py:class:`str<python.str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data. 
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    gb = data_input.groupby(data_input[dim].dt.month)
    return gb - gb.mean(dim = dim)

def calc_monthly_climatological_std_without_seasonal_cycle_mean(
    data_input: xr.DataArray | xr.Dataset,
    dim = 'time',
    **kwargs
) -> xr.DataArray:
    """
    Calculate the standard deviation of monthly data anomalies over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

    .. caution:: `data_input` must be **monthly** data.

    dim : :py:class:`str<python.str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating standard deviation on this object's data. 
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return remove_seasonal_cycle_mean(data_input, dim = dim).std(dim = dim, **kwargs)

def calc_monthly_climatological_var_without_seasonal_cycle_mean(
    data_input: xr.DataArray | xr.Dataset,
    dim = 'time',
    **kwargs
) -> xr.DataArray:
    """
    Calculate the variance of monthly data anomalies over the entire time range.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

    .. caution:: `data_input` must be **monthly** data.

    dim : :py:class:`str<python.str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating variance on this object's data. 
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return remove_seasonal_cycle_mean(data_input, dim = dim).var(dim = dim, **kwargs)

def calc_horizontal_wind_components_std(
    uv_dataset: xr.Dataset,
    u = 'u',
    v = 'v',
    time_dim = 'time',
    ddof = 0
) -> xr.Dataset:
    '''
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
    u : str, default: `u`
        Variable name for the u velocity (in x direction).
    v : str, default: `v`
        Variable name for the v velocity (in y direction).
    time_dim : str, default: `time`
        Dimension(s) over which to apply. By default is applied over the `time` dimension.
    ddof : int, default: 1
        If `ddof=1`, covariance is normalized by `N-1`, giving an unbiased estimate, else normalization is by `N`.

    Returns
    -------
    :py:class:`xarray.Dataset<xarray.Dataset>`
        - sigma_s: the standard deviation of vector wind speed.
        - sigma_d: the standard deviation of vector wind direction.
    
    Reference
    --------------
    G. R. Ackermann. (1983). Means and Standard Deviations of Horizontal Wind Components. 
    Website: https://doi.org/10.1175/1520-0450(1983)022%3C0959:MASDOH%3E2.0.CO;2    
    '''
    U = uv_dataset[u].mean(dim = time_dim)
    V = uv_dataset[v].mean(dim = time_dim)
    S = np.hypot(U, V)
    # D = np.arctan(U / V)

    sigma2_u = uv_dataset[u].var(dim = time_dim)
    sigma2_v = uv_dataset[v].var(dim = time_dim)
    sigma_uv = xr.cov(uv_dataset[u], uv_dataset[v], dim = time_dim, ddof = ddof)

    sigma_s = (U**2 * sigma2_u + V**2 * sigma2_v + 2 * U * V * sigma_uv)**(1/2) * S**(-1)
    sigma_d = (V**2 * sigma2_u + U**2 * sigma2_v - 2 * U * V * sigma_uv)**(1/2) * S**(-2)

    return uv_dataset.assign({'sigma_s': sigma_s, 'sigma_d': sigma_d})

def transfer_monmean2everymonthmean(
    data_input: xr.DataArray, 
    time_dim: str = 'time'
) -> xr.DataArray:
    """
    Convert to the month-mean state corresponding to each month.

    Parameters
    ----------
    - data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.    
    """
    time_step_all = data_input[time_dim].shape[0]
    month_int = data_input.time.dt.month
    month_climate = data_input.groupby(time_dim + '.month').mean(dim = time_dim)
    data_input_empty = xr.full_like(data_input, fill_value = np.nan)

    for time_step in np.arange(0, time_step_all):
        time_step_month = month_int.isel(time = time_step).data
        data_input_empty[{time_dim: time_step}] = month_climate.sel(month = time_step_month)

    return data_input_empty

def mapping_daily_climatological_mean2every_day(
    data_daily: xr.DataArray,
    data_climatology_daily_data: xr.DataArray,
    timd_dim: str = 'time',
) -> xr.DataArray:
    """
    
    """
    result_data = xr.full_like(data_daily, fill_value=np.nan)
    time_length = result_data[timd_dim].shape[0]
    time_dayofyear = result_data[timd_dim].dt.dayofyear.data
    climate_data_dayofyear_index = data_climatology_daily_data[timd_dim].dt.dayofyear.data

    for time_item in np.arange(time_length):
        # Target dayofyear index
        time_dayofyear_item = time_dayofyear[time_item]
        # Correspond dayofyear index in the climate data
        dayofyear_index = np.transpose(np.nonzero(climate_data_dayofyear_index == time_dayofyear_item))
        dayofyear_index = int(dayofyear_index)

        result_data[{timd_dim: time_item}] = data_climatology_daily_data.isel({timd_dim: dayofyear_index})   

    return result_data

def mapping_monthly_climatological_mean2every_month(
    data_monthly: xr.DataArray,
    data_climatology_monthly_data: xr.DataArray,
    timd_dim: str = 'time',
) -> xr.DataArray:
    """
    
    """
    result_data = xr.full_like(data_monthly, fill_value=np.nan)
    time_length = result_data[timd_dim].shape[0]
    time_month = result_data[timd_dim].dt.month.data
    climate_data_month_index = data_climatology_monthly_data[timd_dim].dt.month.data

    for time_item in np.arange(time_length):
        # Target month index
        time_month_item = time_month[time_item]
        # Correspond month index in the climate data
        month_index = np.transpose(np.nonzero(climate_data_month_index == time_month_item))
        month_index = int(month_index)

        result_data[{timd_dim: time_item}] = data_climatology_monthly_data.isel({timd_dim: month_index})   

    return result_data

def calc_daily_climatological_anomaly(
    data_daily: xr.DataArray | xr.Dataset,
    data_climatology_daily_data: xr.DataArray | xr.Dataset,
    timd_dim = 'time',
) -> xr.DataArray:
    """
    
    """
    data_daily_anomaly = data_daily.groupby(data_daily[timd_dim].dt.dayofyear) - data_climatology_daily_data.groupby(data_climatology_daily_data[timd_dim].dt.dayofyear).mean()
    data_daily_anomaly = data_daily_anomaly.drop('dayofyear')
    return data_daily_anomaly