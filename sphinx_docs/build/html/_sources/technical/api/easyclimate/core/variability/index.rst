:py:mod:`easyclimate.core.variability`
======================================

.. py:module:: easyclimate.core.variability

.. autoapi-nested-parse::

   This module calculate climate variability



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.core.variability.calc_climatological_mean
   easyclimate.core.variability.calc_climatological_seasonal_mean
   easyclimate.core.variability.calc_seasonal_cycle_mean
   easyclimate.core.variability.calc_seasonal_cycle_std
   easyclimate.core.variability.calc_seasonal_cycle_var
   easyclimate.core.variability.remove_seasonal_cycle_mean
   easyclimate.core.variability.calc_climate_monthly_std
   easyclimate.core.variability.calc_climate_monthly_var
   easyclimate.core.variability.calc_horizontal_wind_components_std



.. py:function:: calc_climatological_mean(data_input: xarray.DataArray, dim='time', **kwargs) -> xarray.DataArray

   Calculation of the climatological mean over the entire time range.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   dim : :py:class:`str<python.str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data. 
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: calc_climatological_seasonal_mean(data_input: xarray.DataArray, dim='time', **kwargs) -> xarray.DataArray

   Calculation of the seasonal climatological mean over the entire time range.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   dim : :py:class:`str<python.str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data. 
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: calc_seasonal_cycle_mean(data_input: xarray.DataArray, dim='time', **kwargs) -> xarray.DataArray

   Calculation of the seasonal cycle means over the entire time range.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
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


.. py:function:: calc_seasonal_cycle_std(data_input: xarray.DataArray, dim='time', **kwargs) -> xarray.DataArray

   Calculation of the seasonal cycle standard deviation over the entire time range.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
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


.. py:function:: calc_seasonal_cycle_var(data_input: xarray.DataArray, dim='time', **kwargs) -> xarray.DataArray

   Calculation of the seasonal cycle standard deviation over the entire time range.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
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


.. py:function:: remove_seasonal_cycle_mean(data_input: xarray.DataArray, dim='time', **kwargs) -> xarray.DataArray

   Remove of the seasonal cycle means over the entire time range.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
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


.. py:function:: calc_climate_monthly_std(data_input: xarray.DataArray, dim='time', **kwargs) -> xarray.DataArray

   Calculate the standard deviation of monthly data anomalies over the entire time range.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
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


.. py:function:: calc_climate_monthly_var(data_input: xarray.DataArray, dim='time', **kwargs) -> xarray.DataArray

   Calculate the variance of monthly data anomalies over the entire time range.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
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


.. py:function:: calc_horizontal_wind_components_std(uv_dataset: xarray.Dataset, u='u', v='v', time_dim='time', ddof=0) -> xarray.Dataset

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


