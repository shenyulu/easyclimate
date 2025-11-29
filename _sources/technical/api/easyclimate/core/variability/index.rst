easyclimate.core.variability
============================

.. py:module:: easyclimate.core.variability

.. autoapi-nested-parse::

   This module calculate climate variability



Functions
---------

.. autoapisummary::

   easyclimate.core.variability.calc_all_climatological_mean
   easyclimate.core.variability.calc_seasonal_climatological_mean
   easyclimate.core.variability.calc_seasonal_cycle_mean
   easyclimate.core.variability.calc_seasonal_cycle_std
   easyclimate.core.variability.calc_seasonal_cycle_var
   easyclimate.core.variability.calc_seasonal_mean
   easyclimate.core.variability.remove_seasonal_cycle_mean
   easyclimate.core.variability.calc_monthly_climatological_std_without_seasonal_cycle_mean
   easyclimate.core.variability.calc_monthly_climatological_var_without_seasonal_cycle_mean
   easyclimate.core.variability.smooth_daily_annual_cycle
   easyclimate.core.variability.calc_daily_annual_cycle_mean
   easyclimate.core.variability.calc_daily_annual_cycle_std
   easyclimate.core.variability.calc_daily_annual_cycle_var
   easyclimate.core.variability.remove_smooth_daily_annual_cycle_mean
   easyclimate.core.variability.calc_horizontal_wind_components_std
   easyclimate.core.variability.calc_windspeed_dataset
   easyclimate.core.variability.calc_windspeed_dataarray
   easyclimate.core.variability.populate_monmean2everymon
   easyclimate.core.variability.populate_daymean2everyday
   easyclimate.core.variability.calc_daily_climatological_anomaly
   easyclimate.core.variability.remove_low_frequency_signal


Module Contents
---------------

.. py:function:: calc_all_climatological_mean(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', **kwargs) -> xarray.DataArray

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


.. py:function:: calc_seasonal_climatological_mean(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', **kwargs) -> xarray.DataArray

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


.. py:function:: calc_seasonal_cycle_mean(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', **kwargs) -> xarray.DataArray

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


.. py:function:: calc_seasonal_cycle_std(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', **kwargs) -> xarray.DataArray

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


.. py:function:: calc_seasonal_cycle_var(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', **kwargs) -> xarray.DataArray

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


.. py:function:: calc_seasonal_mean(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', extract_season: Literal['seasonly', 'DJF', 'MAM', 'JJA', 'SON', 'JJAS', 'OND'] = 'seasonly', **kwargs) -> xarray.DataArray

   Calculation of the seasonal means per year over the entire time range.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

   .. caution:: `data_input` must be **monthly** data.

   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   extract_season: :py:class:`list <list>`. default: "seasonly".
       Extraction seasons.

       - "seasonly": Calculate by meteorological seasons (DJF, MAM, JJA, SON)
       - custom seasons: Calculate climatology by custom seasons, e.g., `'DJF'`, `'MAM'`, `'JJA'`, `'SON'`, ``'JJAS'``, ``'OND'``.

   Returns
   -------
   :py:class:`xarray.DataArray <xarray.DataArray>`.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_ao_index.py
       ./dynamic_docs/plot_oceanic_front.py
       ./dynamic_docs/plot_multi_linear_reg.py


.. py:function:: remove_seasonal_cycle_mean(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', time_range: slice = slice(None, None)) -> xarray.DataArray

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


.. py:function:: calc_monthly_climatological_std_without_seasonal_cycle_mean(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', **kwargs) -> xarray.DataArray

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


.. py:function:: calc_monthly_climatological_var_without_seasonal_cycle_mean(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', **kwargs) -> xarray.DataArray

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


.. py:function:: smooth_daily_annual_cycle(daily_annual_cycle_data: xarray.DataArray, harmonics_num: int = 3, time_dim: str = 'dayofyear') -> xarray.DataArray

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


.. py:function:: calc_daily_annual_cycle_mean(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', **kwargs) -> xarray.DataArray

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


.. py:function:: calc_daily_annual_cycle_std(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', **kwargs) -> xarray.DataArray

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


.. py:function:: calc_daily_annual_cycle_var(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', **kwargs) -> xarray.DataArray

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


.. py:function:: remove_smooth_daily_annual_cycle_mean(data_input: xarray.DataArray, daily_cycle_mean_time_range: slice = slice(None, None), extract_time_range: slice = slice(None, None), harmonics_num: int = 3, dim: str = 'time')

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


.. py:function:: calc_horizontal_wind_components_std(uv_dataset: xarray.Dataset, u_dim='u', v_dim='v', time_dim='time', ddof=0) -> xarray.Dataset

   Calculate the standard deviation of vector wind speed and direction.

   The standard deviation of vector wind speed

   .. math::
       \sigma_s = [U^2 \sigma_u^2 + V^2 \sigma_v^2 + 2 U V \sigma_{uv}]^{1/2} S^{-1},

   The standard deviation of vector wind direction

   .. math::
       \sigma_d = [V^2 \sigma_u^2 + U^2 \sigma_v^2 + 2 U V \sigma_{uv}]^{1/2} S^{-2},

   Where time mean of :math:`u` is :math:`U = n^{-1} \sum u_i`, time mean of :math:`v` is :math:`V = n^{-1} \sum v_i`,
   time variance of :math:`u` is :math:`\sigma_u^2 = n^{-1} \sum u_{i}^{2} - U^2`,
   time variance of :math:`v` is :math:`\sigma_v^2 = n^{-1} \sum v_{i}^{2} - V^2`,
   time covariance of :math:`u`, :math:`v` is :math:`\sigma_{uv} = n^{-1} \sum u_i v_i - UV`,
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


.. py:function:: calc_windspeed_dataset(uv_dataset: xarray.Dataset, u_dim: str = 'u', v_dim: str = 'v') -> xarray.Dataset

   Calculate the horizontal wind speed from zonal and meridional wind components in a :py:class:`xarray.Dataset<xarray.Dataset>`.

   The wind speed is computed as the magnitude of the horizontal wind vector:

   .. math::

       S = \sqrt{u^2 + v^2},

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


.. py:function:: calc_windspeed_dataarray(u_data: xarray.DataArray, v_data: xarray.DataArray) -> xarray.DataArray

   Calculate the horizontal wind speed from zonal and meridional wind components in :py:class:`xarray.DataArray<xarray.DataArray>`.

   The wind speed is computed as the magnitude of the horizontal wind vector:

   .. math::

       S = \sqrt{u^2 + v^2},

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


.. py:function:: populate_monmean2everymon(data_monthly: xarray.DataArray, data_climatology_monthly_data: xarray.DataArray = None, time_dim: str = 'time') -> xarray.DataArray

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


.. py:function:: populate_daymean2everyday(data_daily: xarray.DataArray, data_climatology_daily_data: xarray.DataArray = None, time_dim: str = 'time') -> xarray.DataArray

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


.. py:function:: calc_daily_climatological_anomaly(data_daily: xarray.DataArray | xarray.Dataset, data_climatology_daily_data: xarray.DataArray | xarray.Dataset, timd_dim='time') -> xarray.DataArray | xarray.Dataset

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


.. py:function:: remove_low_frequency_signal(da: xarray.DataArray, window: int = 120, center: bool = False, time_dim: str = 'time') -> xarray.DataArray

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


