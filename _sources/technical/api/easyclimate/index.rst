easyclimate
===========

.. py:module:: easyclimate


Submodules
----------

.. toctree::
   :maxdepth: 1

   /technical/api/easyclimate/accessor/index
   /technical/api/easyclimate/backend/index
   /technical/api/easyclimate/core/index
   /technical/api/easyclimate/field/index
   /technical/api/easyclimate/filter/index
   /technical/api/easyclimate/interp/index
   /technical/api/easyclimate/physics/index
   /technical/api/easyclimate/plot/index
   /technical/api/easyclimate/version/index
   /technical/api/easyclimate/wrf/index


Attributes
----------

.. autoapisummary::

   easyclimate.__version__


Classes
-------

.. autoapisummary::

   easyclimate.DataNode


Functions
---------

.. autoapisummary::

   easyclimate.show_versions
   easyclimate.calc_gradient
   easyclimate.calc_lon_gradient
   easyclimate.calc_lat_gradient
   easyclimate.calc_lon_laplacian
   easyclimate.calc_lat_laplacian
   easyclimate.calc_lon_lat_mixed_derivatives
   easyclimate.calc_p_gradient
   easyclimate.calc_time_gradient
   easyclimate.calc_delta_pressure
   easyclimate.calc_p_integral
   easyclimate.calc_top2surface_integral
   easyclimate.calc_laplacian
   easyclimate.calc_divergence
   easyclimate.calc_vorticity
   easyclimate.calc_geostrophic_wind
   easyclimate.calc_geostrophic_wind_vorticity
   easyclimate.calc_horizontal_water_flux
   easyclimate.calc_vertical_water_flux
   easyclimate.calc_water_flux_top2surface_integral
   easyclimate.calc_divergence_watervaporflux
   easyclimate.calc_divergence_watervaporflux_top2surface_integral
   easyclimate.calc_u_advection
   easyclimate.calc_v_advection
   easyclimate.calc_p_advection
   easyclimate.calc_eady_growth_rate
   easyclimate.calc_apparent_heat_source
   easyclimate.calc_total_diabatic_heating
   easyclimate.calc_apparent_moisture_sink
   easyclimate.calc_Plumb_wave_activity_horizontal_flux
   easyclimate.calc_TN_wave_activity_horizontal_flux
   easyclimate.calc_EP_horizontal_flux
   easyclimate.get_specific_years_data
   easyclimate.get_specific_months_data
   easyclimate.get_specific_days_data
   easyclimate.get_specific_hours_data
   easyclimate.get_specific_minutes_data
   easyclimate.get_specific_seconds_data
   easyclimate.get_specific_microseconds_data
   easyclimate.get_specific_nanoseconds_data
   easyclimate.get_specific_dayofweek_data
   easyclimate.get_yearmean_for_specific_months_data
   easyclimate.get_year_exceed_index_upper_bound
   easyclimate.get_year_exceed_index_lower_bound
   easyclimate.get_time_exceed_index_upper_bound
   easyclimate.get_time_exceed_index_lower_bound
   easyclimate.open_muliti_dataset
   easyclimate.calc_linregress_spatial
   easyclimate.calc_detrend_spatial
   easyclimate.calc_corr_spatial
   easyclimate.calc_leadlag_corr_spatial
   easyclimate.calc_multiple_linear_regression_spatial
   easyclimate.calc_ttestSpatialPattern_spatial
   easyclimate.calc_windmask_ttestSpatialPattern_spatial
   easyclimate.calc_levenetestSpatialPattern_spatial
   easyclimate.calc_skewness_spatial
   easyclimate.calc_kurtosis_spatial
   easyclimate.calc_theilslopes_spatial
   easyclimate.calc_lead_lag_correlation_coefficients
   easyclimate.calc_timeseries_correlations
   easyclimate.calc_non_centered_corr
   easyclimate.calc_pattern_corr
   easyclimate.calc_all_climatological_mean
   easyclimate.calc_seasonal_climatological_mean
   easyclimate.calc_seasonal_cycle_mean
   easyclimate.calc_seasonal_cycle_std
   easyclimate.calc_seasonal_cycle_var
   easyclimate.calc_seasonal_mean
   easyclimate.remove_seasonal_cycle_mean
   easyclimate.calc_monthly_climatological_std_without_seasonal_cycle_mean
   easyclimate.calc_monthly_climatological_var_without_seasonal_cycle_mean
   easyclimate.smooth_daily_annual_cycle
   easyclimate.calc_daily_annual_cycle_mean
   easyclimate.calc_daily_annual_cycle_std
   easyclimate.calc_daily_annual_cycle_var
   easyclimate.remove_smooth_daily_annual_cycle_mean
   easyclimate.calc_horizontal_wind_components_std
   easyclimate.calc_windspeed_dataset
   easyclimate.calc_windspeed_dataarray
   easyclimate.populate_monmean2everymon
   easyclimate.populate_daymean2everyday
   easyclimate.calc_daily_climatological_anomaly
   easyclimate.remove_low_frequency_signal
   easyclimate.calc_yearly_climatological_mean
   easyclimate.calc_yearly_climatological_sum
   easyclimate.calc_yearly_climatological_std
   easyclimate.calc_yearly_climatological_var
   easyclimate.calc_yearly_climatological_max
   easyclimate.calc_yearly_climatological_min
   easyclimate.open_tutorial_dataset
   easyclimate.open_datanode


Package Contents
----------------

.. py:data:: __version__
   :value: '2025.11.0'


.. py:function:: show_versions() -> str

.. py:function:: calc_gradient(data_input: xarray.DataArray | xarray.Dataset, dim: str, varargs: int = 1, edge_order: int = 2) -> xarray.DataArray | xarray.Dataset

   Compute the gradient along the coordinate `dim` direction.

   The gradient is computed using **second order accurate central differences** in the interior points
   and either first or second order accurate one-sides (forward or backwards) differences at the boundaries.
   The returned gradient hence has the same shape as the input array.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
        The spatio-temporal data to be calculated.
   dim : :py:class:`str <str>`.
       Dimension(s) over which to apply gradient. By default gradient is applied over the `time` dimension.
   varargs: :py:class:`list <list>` of scalar or array, optional
       Spacing between f values. Default unitary spacing for all dimensions. Spacing can be specified using:

       1. Single scalar to specify a sample distance for all dimensions.
       2. N scalars to specify a constant sample distance for each dimension. i.e. :math:`\mathrm{d}x, \mathrm{d}y, \mathrm{d}z, ...`
       3. N arrays to specify the coordinates of the values along each dimension of F.
          The length of the array must match the size of the corresponding dimension.
       4. Any combination of N scalars/arrays with the meaning of 2. and 3.

   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.

   Returns
   -------
   The gradient along the coordinate `dim` direction (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`numpy.gradient <numpy:numpy.gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_lon_gradient(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon', lat_dim: str = 'lat', min_dx: float = 1.0, edge_order: int = 2, R: float = 6370000) -> xarray.DataArray | xarray.Dataset

   Calculate the gradient along the longitude.

   .. math::
       \frac{\partial F}{\partial x} = \frac{1}{R \cos\varphi} \cdot \frac{\partial F}{\partial \lambda}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   min_dx: :py:class:`float <float>`, default: `1.0`.
       The minimum acceptable value of `dx`, below which parts will set `nan` to avoid large computational errors.
       The unit is m. You can set it to a negative value in order to remove this benefit.
   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 2.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The gradient along the longitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_lat_gradient(data_input: xarray.DataArray | xarray.Dataset, lat_dim: str = 'lat', min_dy: float = 1.0, edge_order: int = 2, R: float = 6370000) -> xarray.DataArray | xarray.Dataset

   Calculate the gradient along the latitude.

   .. math::
       \frac{\partial F}{\partial y} = \frac{1}{R} \cdot \frac{\partial F}{\partial \varphi}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   min_dy: :py:class:`float <float>`, default: `1.0`.
       The minimum acceptable value of `dy`, below which parts will set `nan` to avoid large computational errors.
       The unit is m. You can set it to a negative value in order to remove this benefit.
   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The gradient along the latitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_lon_laplacian(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon', lat_dim: str = 'lat', min_dx2: float = 1000000000.0, edge_order: int = 2, R: float = 6370000) -> xarray.DataArray | xarray.Dataset

   Calculation of the second-order partial derivative term (Laplace term) along longitude.

   .. math::
       \frac{\partial^2 F}{\partial x^2} = \frac{1}{(R \cos\varphi)^2} \cdot \frac{\partial^2 F}{\partial \lambda^2}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   min_dx2: :py:class:`float <float>`, default: `1e9`.
       The minimum acceptable value of :math:`(\mathrm{d}x)^2`, below which parts will set `nan` to avoid large computational errors.
       The unit is m. You can set it to a negative value in order to remove this benefit.
   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The second-order partial derivative term (Laplace term) along longitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_lat_laplacian(data_input: xarray.DataArray | xarray.Dataset, lat_dim: str = 'lat', min_dy2: float = 1.0, edge_order: int = 2, R: float = 6370000) -> xarray.DataArray | xarray.Dataset

   Calculation of the second-order partial derivative term (Laplace term) along latitude.

   .. math::
       \frac{\partial^2 F}{\partial y^2} = \frac{1}{R^2} \cdot \frac{\partial^2 F}{\partial \varphi^2}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   min_dy2: :py:class:`float <float>`, default: `1.0`.
       The minimum acceptable value of :math:`(\mathrm{d}y)^2`, below which parts will set `nan` to avoid large computational errors.
       The unit is m. You can set it to a negative value in order to remove this benefit.
   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The second-order partial derivative term (Laplace term) along latitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_lon_lat_mixed_derivatives(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon', lat_dim: str = 'lat', min_dxdy: float = 10000000000.0, edge_order: int = 2, R: float = 6370000) -> xarray.DataArray | xarray.Dataset

   Calculation of second-order mixed partial derivative terms along longitude and latitude.

   .. math::
       \frac{\partial^2 F}{\partial x \partial y} = \frac{1}{R^2 \cos\varphi} \cdot \frac{\partial^2 F}{\partial \lambda \partial \varphi}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   min_dxdy: :py:class:`float <float>`, default: `1e10`.
       The minimum acceptable value of :math:`\mathrm{d}x\mathrm{d}y`, below which parts will set `nan` to avoid large computational errors.
       The unit is m. You can set it to a negative value in order to remove this benefit.
   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The second-order mixed partial derivative terms along longitude and latitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_p_gradient(data_input: xarray.DataArray, vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar']) -> xarray.DataArray

   Calculate the gradient along the barometric pressure direction in the p-coordinate system.

   .. math::
       \frac{\partial F}{\partial p}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.

   Returns
   -------
   The gradient along the barometric pressure direction in the p-coordinate system (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>`


.. py:function:: calc_time_gradient(data_input: xarray.DataArray, time_units: str, time_dim: str = 'time') -> xarray.DataArray

   Calculate the gradient along the time direction.

   .. math::
       \frac{\partial F}{\partial t}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   time_units: :py:class:`str <str>`.
       The unit corresponding to the time dimension value. Optional values are `seconds`, `months`, `years` and so on.
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name.

   Returns
   -------
   The gradient along the time direction (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. caution:: The units for partial derivative of `time` are :math:`\mathrm{s^{-1}}`.

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>`


.. py:function:: calc_delta_pressure(data_input: xarray.DataArray, surface_pressure_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: str, surface_pressure_data_units: str) -> xarray.DataArray

   Calculates the pressure layer thickness (delta pressure) of a constant
   pressure level coordinate system.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The spatio-temporal data to be calculated.
   surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Mean surface sea level pressure.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   surface_pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.

   Returns
   -------
   The pressure layer thickness (delta pressure) of a constant pressure level coordinate system (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       - :py:func:`geocat.comp.meteorology.delta_pressure <geocat-comp:geocat.comp.meteorology.delta_pressure>`
       - `dpres_plevel - NCL <https://www.ncl.ucar.edu/Document/Functions/Built-in/dpres_plevel.shtml>`__


.. py:function:: calc_p_integral(data_input: xarray.DataArray, vertical_dim: str, normalize: bool = True) -> xarray.DataArray

   Calculate the vertical integral along the barometric pressure direction in the p-coordinate system.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The spatio-temporal data to be calculated.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   normalize: :py:class:`bool<bool>`, default: `True`.
       Whether or not the integral results are averaged over the entire layer.

   Returns
   -------
   The vertical integral along the barometric pressure direction in the p-coordinate system (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. attention::
       This method ignores the effect of topography, so it applies to altitudes **above 900hPa** and is **NOT applicable to the Tibetan Plateau region**.
       For a fully accurate vertical integration, please use the :py:func:`calc_top2surface_integral <calc_top2surface_integral>` function to calculate,
       but the speed of the calculation is slightly slowed down.


.. py:function:: calc_top2surface_integral(data_input: xarray.DataArray, surface_pressure_data: xarray.DataArray, vertical_dim: str, surface_pressure_data_units: str, vertical_dim_units: str, method: Literal['Boer-vibeta', 'Trenberth-vibeta'] = 'Trenberth-vibeta', normalize: bool = True) -> xarray.DataArray

   Calculate the vertical integral in the p-coordinate system from the ground to the zenith along the barometric pressure direction.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The spatio-temporal data to be calculated.
   surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Mean surface sea level pressure.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   surface_pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   method: :py:class:`str <str>`, default: `'Trenberth-vibeta'`.
       vertical integration method. Optional values are `Boer-vibeta`, `'Trenberth-vibeta'`.

       .. note::
           The trapezoidal rule of integration is exactly equivalent to

           .. math::
               I = \sum_{j=1,2J-1,2} (\beta M)_j \Delta p_j,

           where Kevin E. Trenberth (1991) define

           .. math::
               \beta_j = \left\lbrace
               \begin{array}{ll}
               1, & \mathrm{if} \ p_{j-1} < p_s,\\
               0, & \mathrm{if} \ p_{j+1} > p_s ,\\
               \frac{p_s - p_{j+1}}{p_{j-1} - p_{j+1}}, & \mathrm{if}  \ p_{j-1} > p_s > p_{j+1}.
               \end{array}
               \right.

           While G. J. Boer (1982) define :math:`\beta = 0, 1` only.

   normalize: :py:class:`bool<bool>`, default: `True`.
       Whether or not the integral results are averaged over the entire layer.

   Returns
   -------
   The vertical integral in the p-coordinate system from the ground to the zenith along the barometric pressure direction. (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - `Boer, G. J., 1982: Diagnostic Equations in Isobaric Coordinates. Mon. Wea. Rev., 110, 1801–1820, <https://doi.org/10.1175/1520-0493(1982)110%3C1801:DEIIC%3E2.0.CO;2>`__
   - `Trenberth, K. E., 1991: Climate Diagnostics from Global Analyses: Conservation of Mass in ECMWF Analyses. J. Climate, 4, 707–722, <https://doi.org/10.1175/1520-0442(1991)004%3C0707:CDFGAC%3E2.0.CO;2>`__

   .. seealso::
       - `vibeta - NCL <https://www.ncl.ucar.edu/Document/Functions/Built-in/vibeta.shtml>`__
       - `dpres_plevel - NCL <https://www.ncl.ucar.edu/Document/Functions/Built-in/dpres_plevel.shtml>`__


.. py:function:: calc_laplacian(data_input: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', R: float = 6370000, spherical_coord: bool = True) -> xarray.DataArray

   Calculate the horizontal Laplace term.

   rectangular coordinates

   .. math::
       \nabla^2 F = \frac{\partial^2 F}{\partial x^2} + \frac{\partial^2 F}{\partial y^2}

   Spherical coordinates

   .. math::
       \nabla^2 F = \frac{\partial^2 F}{\partial x^2} + \frac{\partial^2 F}{\partial y^2} - \frac{1}{R} \frac{\partial F}{\partial y} \tan \varphi

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The spatio-temporal data to be calculated.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.
   spherical_coord: :py:class:`bool <bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates.

   Returns
   -------
   The horizontal Laplace term. (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_divergence(u_data: xarray.DataArray, v_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', R: float = 6370000, spherical_coord=True) -> xarray.DataArray

   Calculate the horizontal divergence term.

   rectangular coordinates

   .. math::
       \mathrm{D} = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y}

   Spherical coordinates

   .. math::
       \mathrm{D} = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} - \frac{v}{R} \tan \varphi

   Parameters
   ----------
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.
   spherical_coord: :py:class:`bool<bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates.

   Returns
   -------
   The horizontal divergence term. (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_vorticity(u_data: xarray.DataArray, v_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', R: float = 6370000, spherical_coord: bool = True) -> xarray.DataArray

   Calculate the horizontal relative vorticity term.

   rectangular coordinates

   .. math::
       \zeta = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}

   Spherical coordinates

   .. math::
       \zeta = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y} + \frac{u}{R} \tan \varphi

   Parameters
   ----------
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.
   spherical_coord: :py:class:`bool<bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates.

   Returns
   -------
   The horizontal relative vorticity term. (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_geostrophic_wind(z_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', omega: float = 7.292e-05, g: float = 9.8, R: float = 6370000) -> xarray.DataArray

   Calculate the geostrophic wind.

   .. math::
       u_g = - \frac{g}{f} \frac{\partial H}{\partial y}

   .. math::
       v_g = \frac{g}{f} \frac{\partial H}{\partial x}

   Parameters
   ----------
   z_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric geopotential height.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   omega: :py:class:`float <float>`, default: `7.292e-5`.
       The angular speed of the earth.
   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The geostrophic wind term. (:py:class:`xarray.DataArray<xarray.DataArray>`).
       - ug
       - vg

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_geostrophic_wind_vorticity(z_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', spherical_coord: bool = True, omega: float = 7.292e-05, g: float = 9.8, R: float = 6370000) -> xarray.DataArray

   Calculate the geostrophic vorticity.

   Rectangular coordinates

   .. math::
       \zeta_g = \frac{\partial v_g}{\partial x} - \frac{\partial u_g}{\partial y}

   Spherical coordinates

   .. math::
       \zeta_g = \frac{\partial v_g}{\partial x} - \frac{\partial u_g}{\partial y} + \frac{u_g}{R} \tan \varphi

   Parameters
   ----------
   z_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric geopotential height.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   spherical_coord: :py:class:`bool<bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates.
   omega: :py:class:`float <float>`, default: `7.292e-5`.
       The angular speed of the earth.
   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The geostrophic vorticity term. (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_horizontal_water_flux(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, g: float = 9.8) -> xarray.Dataset

   Calculate horizontal water vapor flux at each vertical level.

   .. math::
       \frac{1}{g} q \mathbf{V} = \frac{1}{g} (u q\ \mathbf{i} + vq\ \mathbf{j})

   Parameters
   ----------
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.

   Returns
   -------
   The water vapor flux. (:py:class:`xarray.Dataset<xarray.Dataset>`).

   - :math:`qu`: zonal water vapor flux.
   - :math:`qv`: meridional water vapor flux.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_vertical_water_flux(specific_humidity_data: xarray.DataArray, omega_data: xarray.DataArray, g: float = 9.8) -> xarray.DataArray

   Calculate vertical water vapor flux.

   .. math::
       -\omega \frac{q}{g}

   Parameters
   ----------
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   omega_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The vertical velocity data (:math:`\frac{\mathrm{d} p}{\mathrm{d} t}`).
   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.

   Returns
   -------
   The vertical water flux. (:py:class:`xarray.DataArray <xarray.DataArray>`).

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_water_flux_top2surface_integral(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, surface_pressure_data: xarray.DataArray, surface_pressure_data_units: str, vertical_dim: str, vertical_dim_units: str, method: Literal['Boer-vibeta', 'Trenberth-vibeta'] = 'Trenberth-vibeta', g: float = 9.8) -> xarray.DataArray

   Calculate the water vapor flux across the vertical level.

   Parameters
   ----------
   specific_humidity: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Mean surface sea level pressure.
   surface_pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   method: :py:class:`str <str>`, default: `'Trenberth-vibeta'`.
       Vertical integration method. Optional values are `Boer-vibeta`, `'Trenberth-vibeta'`.
   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.

   Returns
   -------
   The water vapor flux. (:py:class:`xarray.Dataset<xarray.Dataset>`).

   - :math:`qu`: zonal water vapor flux.
   - :math:`qv`: meridional water vapor flux.

   .. seealso::
       :py:func:`calc_top2surface_integral <calc_top2surface_integral>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_divergence_watervaporflux(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, specific_humidity_data_units: str, spherical_coord: bool = True, lon_dim: str = 'lon', lat_dim: str = 'lat', g: float = 9.8, R: float = 6370000) -> xarray.DataArray

   Calculate water vapor flux divergence at each vertical level.

   .. math::
       \nabla \left( \frac{1}{g} q \mathbf{V} \right) = \frac{1}{g} \nabla \cdot \left( q \mathbf{V} \right)


   Parameters
   ----------
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   specific_humidity_data_units: :py:class:`str <str>`.
       The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
   spherical_coord: :py:class:`bool<bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The water vapor flux divergence. (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_divergence_watervaporflux_top2surface_integral(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, surface_pressure_data: xarray.DataArray, vertical_dim: str, specific_humidity_data_units: str, surface_pressure_data_units: str, vertical_dim_units: str, spherical_coord: bool = True, lon_dim: str = 'lon', lat_dim: str = 'lat', method: Literal['Boer-vibeta', 'Trenberth-vibeta'] = 'Trenberth-vibeta', g: float = 9.8, R: float = 6370000) -> xarray.DataArray

   Calculate water vapor flux divergence across the vertical level.

   Parameters
   ----------
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Mean surface sea level pressure.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   specific_humidity_data_units: :py:class:`str <str>`.
       The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
   surface_pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   spherical_coord: :py:class:`bool<bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The water vapor flux divergence. (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_u_advection(u_data: xarray.DataArray, temper_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', min_dx: float = 1.0, edge_order: int = 2, R: float = 6370000) -> xarray.DataArray

   Calculate zonal temperature advection at each vertical level.

   .. math::
       -u \frac{\partial T}{\partial x}

   Parameters
   ----------
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   min_dx: :py:class:`float <float>`, default: `1.0`.
       The minimum acceptable value of `dx`, below which parts will set `nan` to avoid large computational errors.
       The unit is m. You can set it to a negative value in order to remove this benefit.
   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The zonal temperature advection. (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_v_advection(v_data: xarray.DataArray, temper_data: xarray.DataArray, lat_dim: str = 'lat', min_dy: float = 1.0, edge_order: int = 2, R: float = 6370000) -> xarray.DataArray

   Calculate meridional temperature advection at each vertical level.

   .. math::
       -v \frac{\partial T}{\partial y}

   Parameters
   ----------
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The meridional temperature advection. (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_p_advection(omega_data: xarray.DataArray, temper_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: str) -> xarray.DataArray

   Calculate vertical temperature transport at each vertical level.

   .. math::
       -\omega \frac{\partial T}{\partial p}

   Parameters
   ----------
   omega: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The vertical velocity data (:math:`\frac{\mathrm{d} p}{\mathrm{d} t}`).
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.

   Returns
   -------
   The vertical temperature transport. (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_eady_growth_rate(u_daily_data: xarray.DataArray, z_daily_data: xarray.DataArray, temper_daily_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], lat_dim='lat', g=9.8) -> xarray.Dataset

   Calculate the maximum Eady growth rate.

   .. math::
       \sigma = 0.3098 \frac{f}{N} \frac{\mathrm{d} U}{\mathrm{d} z}

   .. caution::
       Eady growth rate (EGR) is a non-linear quantity. Hence, `calc_eady_growth_rate` should **NOT** be **directly applied to monthly means** variables.
       If a monthly climatology of EGR is desired, the EGR values at the high frequency temporal time steps should be calculated;
       then, use calculate monthly mean.

   Parameters
   ----------
   u_daily_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind daily data.
   z_daily_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Daily atmospheric geopotential height.

   .. attention:: The unit of `z_daily_data` should be **meters**, NOT :math:`\mathrm{m^2 \cdot s^2}` which is the unit used in the representation of potential energy.

   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Daily air temperature.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.

   Returns
   -------
   The maximum Eady growth rate (:py:class:`xarray.Dataset<xarray.Dataset>`).

   - `eady_growth_rate`: The maximum Eady growth rate.
   - `dudz`: :math:`\frac{\mathrm{d} U}{\mathrm{d} z}`
   - `brunt_vaisala_frequency`: Brunt-väisälä frequency.

   .. seealso::
       - `eady_growth_rate -NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/eady_growth_rate.shtml>`__
       - `瞬变涡旋诊断量 <https://renqlsysu.github.io/2020/02/16/wave_activity_flux/>`__


.. py:function:: calc_apparent_heat_source(u_data: xarray.DataArray, v_data: xarray.DataArray, omega_data: xarray.DataArray, temper_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], time_units: str, lon_dim='lon', lat_dim='lat', time_dim='time', c_p=1005.7) -> xarray.DataArray

   Calculate the apparent heat source.

   .. math::
       Q_1 = C_p \frac{T}{\theta} \left( \frac{\partial \theta}{\partial t} + u \frac{\partial \theta}{\partial x} + v \frac{\partial \theta}{\partial y} + \omega \frac{\partial \theta}{\partial p} \right)

   Parameters
   ----------
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   omega_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The vertical velocity data (:math:`\frac{\mathrm{d} p}{\mathrm{d} t}`).
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   time_units: :py:class:`str <str>`.
       The unit corresponding to the time dimension value. Optional values are `seconds`, `months`, `years` and so on.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   time_dim: :py:class:`str <str>`.
       The time coordinate dimension name.
   c_p: :py:class:`float <float>`, default: `1005.7`.
       The specific heat at constant pressure of dry air.

       .. note::
           `specific heat capacity - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Specific_heat_capacity>`__

   Returns
   -------
   The apparent heat source (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::
       - `Yanai, M., & Tomita, T. (1998). Seasonal and Interannual Variability of Atmospheric Heat Sources and Moisture Sinks as Determined from NCEP–NCAR Reanalysis, Journal of Climate, 11(3), 463-482. <https://journals.ametsoc.org/view/journals/clim/11/3/1520-0442_1998_011_0463_saivoa_2.0.co_2.xml>`__
       - `Ling, J., & Zhang, C. (2013). Diabatic Heating Profiles in Recent Global Reanalyses, Journal of Climate, 26(10), 3307-3325. <https://doi.org/10.1175/JCLI-D-12-00384.1>`__


.. py:function:: calc_total_diabatic_heating(u_data: xarray.DataArray, v_data: xarray.DataArray, omega_data: xarray.DataArray, temper_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], time_units: str, lat_dim='lat', lon_dim='lon', time_dim='time', c_p=1005.7) -> xarray.DataArray

   Calculate the total diabatic heating.

   Calculated in exactly the same way as for the apparent heat source.

   Parameters
   ----------
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   omega_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The vertical velocity data (:math:`\frac{\mathrm{d} p}{\mathrm{d} t}`).
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   time_units: :py:class:`str <str>`.
       The unit corresponding to the time dimension value. Optional values are `seconds`, `months`, `years` and so on.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name.
   c_p: :py:class:`float <float>`, default: `1005.7` (:math:`\mathrm{J \cdot kg^{-1} \cdot K^{-1}}`).
       The specific heat at constant pressure of dry air.

       .. note::
           `specific heat capacity - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Specific_heat_capacity>`__

   Returns
   -------
   The total diabatic heating (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::
       :py:func:`calc_apparent_heat_source <calc_apparent_heat_source>`


.. py:function:: calc_apparent_moisture_sink(u_data: xarray.DataArray, v_data: xarray.DataArray, omega_data: xarray.DataArray, specific_humidity_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], time_units: str, specific_humidity_data_units: str, lon_dim='lon', lat_dim='lat', time_dim='time', latent_heat_of_condensation=2501000.0) -> xarray.DataArray

   Calculate the apparent moisture sink.

   .. math::
       Q_2 = -L \left( \frac{\partial q}{\partial t} + u \frac{\partial q}{\partial x} + v \frac{\partial q}{\partial y} + \omega \frac{\partial q}{\partial p}  \right)

   Parameters
   ----------
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   omega_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The vertical velocity data (:math:`\frac{\mathrm{d} p}{\mathrm{d} t}`).
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   time_units: :py:class:`str <str>`.
       The unit corresponding to the time dimension value. Optional values are `seconds`, `months`, `years` and so on.
   specific_humidity_data_units: :py:class:`str <str>`.
       The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name.
   latent_heat_of_condensation: :py:class:`float <float>`, default: `2.5008e6` (:math:`\mathrm{J \cdot kg^{-1}}`).
       Latent heat of condensation of water at 0°C.

       .. note::
           - `latent heat - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Latent_heat>`__
           - `Latent heat - Wikipedia <https://en.wikipedia.org/wiki/Latent_heat>`__

   Returns
   -------
   The apparent moisture sink (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::
       - `Yanai, M., & Tomita, T. (1998). Seasonal and Interannual Variability of Atmospheric Heat Sources and Moisture Sinks as Determined from NCEP–NCAR Reanalysis, Journal of Climate, 11(3), 463-482. <https://journals.ametsoc.org/view/journals/clim/11/3/1520-0442_1998_011_0463_saivoa_2.0.co_2.xml>`__
       - `HAO Lisheng, MA Ning, HE Liye. Circulation anomalies characteritics of the abnormal drought and high temperature event in the middle and lower reaches of the Yangtze River in summer of 2022[J]. Arid Meteorology, 2022, 40(5): 721-732 <https://doi.org/10.11755/j.issn.1006-7639(2022)-05-0721>`__


.. py:function:: calc_Plumb_wave_activity_horizontal_flux(z_prime_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], lon_dim='lon', lat_dim='lat', omega=7.292e-05, g=9.8, R=6370000) -> xarray.Dataset

   Calculate Plumb wave activity horizontal flux.

   Parameters
   ----------
   z_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The anormaly of atmospheric geopotential height.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   omega: :py:class:`float <float>`, default: `7.292e-5`.
       The angular speed of the earth.
   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The Plumb wave activity horizontal flux (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::
       - `Plumb, R. A., 1985: On the Three-Dimensional Propagation of Stationary Waves. J. Atmos. Sci., 42, 217–229 <https://journals.ametsoc.org/view/journals/atsc/42/3/1520-0469_1985_042_0217_ottdpo_2_0_co_2.xml>`__


.. py:function:: calc_TN_wave_activity_horizontal_flux(z_prime_data: xarray.DataArray, u_climatology_data: xarray.DataArray, v_climatology_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], lon_dim: str = 'lon', lat_dim: str = 'lat', omega: float = 7.292e-05, g: float = 9.8, R: float = 6370000) -> xarray.DataArray

   Calculate TN wave activity horizontal flux.

   .. math::
       \mathbf{W_h} = \frac{p\cos\varphi}{2\lvert \mathbf{U_c} \rvert}\begin{pmatrix}
                             \frac{U_c}{R^2 \cos^2 \varphi} \left[ \left( \frac{\partial \psi'}{\partial \lambda} \right)^2 - \psi'\frac{\partial^2 \psi'}{\partial \lambda^2} \right] + \frac{V_c}{R^2 \cos \varphi} \left[ \frac{\partial \psi'}{\partial \lambda} \frac{\partial \psi'}{\partial \varphi} - \psi' \frac{\partial^2 \psi'}{\partial \lambda \partial \varphi} \right] \\
                             \frac{U_c}{R^2 \cos \varphi} \left[ \frac{\partial \psi'}{\partial \lambda} \frac{\partial \psi'}{\partial \varphi} - \psi' \frac{\partial^2 \psi'}{\partial \lambda \partial \varphi} \right] + \frac{V_c}{R^2} \left[ \left( \frac{\partial \psi'}{\partial \varphi} \right)^2 - \psi'\frac{\partial^2 \psi'}{\partial \varphi^2} \right] \\
                              \end{pmatrix}

   Parameters
   ----------
   z_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The anormaly of atmospheric geopotential height.
   u_climatology_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The climatology of zonal wind data.
   v_climatology_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The climatology of meridional wind data.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   omega: :py:class:`float <float>`, default: `7.292e-5`.
       The angular speed of the earth.
   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The TN wave activity horizontal flux (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::
       - http://www.atmos.rcast.u-tokyo.ac.jp/nishii/programs/index.html
       - https://github.com/laishenggx/T-N_Wave-Activity-Flux


.. py:function:: calc_EP_horizontal_flux(u_prime_data: xarray.DataArray, v_prime_data: xarray.DataArray, time_dim: str = 'time', lat_dim: str = 'lat') -> xarray.Dataset

   Calculate horizontal Eliassen–Palm Flux.

   Parameters
   ----------
   u_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The anormaly of zonal wind data.
   v_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The anormaly of meridional wind data.
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The Eliassen–Palm Flux (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::
       - https://www.ncl.ucar.edu/Applications/EPflux.shtml
       - https://renqlsysu.github.io/2020/02/16/wave_activity_flux/


.. py:function:: get_specific_years_data(data_input: xarray.DataArray | xarray.Dataset, year_array: np.array(int) | List[int], dim: str = 'time') -> xarray.DataArray | xarray.Dataset

   Slicing and extracting the part of the data containing the specified year based on an array of given integer years.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   year_array: :py:class:`list[int]`
       Year(s) to be extracted.
   dim : :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_months_data(data_input: xarray.DataArray | xarray.Dataset, month_array: numpy.array, dim: str = 'time') -> xarray.DataArray | xarray.Dataset

   Slicing and extracting the part of the data containing the specified year based on an array of given integer months.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   month_array: :py:class:`list[int]`
       Month(s) to be extracted.
   dim : :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_basic_statistical_analysis.py


.. py:function:: get_specific_days_data(data_input: xarray.DataArray | xarray.Dataset, day_array: np.array(int) | List[int], dim: str = 'time') -> xarray.DataArray | xarray.Dataset

   Slicing and extracting the part of the data containing the specified year based on an array of given integer days.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   day_array: :py:class:`list[int]`
       Days(s) to be extracted.
   dim : :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_hours_data(data_input: xarray.DataArray | xarray.Dataset, hour_array: np.array(int) | List[int], dim: str = 'time') -> xarray.DataArray | xarray.Dataset

   Slicing and extracting the part of the data containing the specified year based on an array of given integer hours.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   hour_array: :py:class:`list[int]`
       Hour(s) to be extracted.
   dim : :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_minutes_data(data_input: xarray.DataArray | xarray.Dataset, minute_array: np.array(int) | List[int], dim: str = 'time') -> xarray.DataArray | xarray.Dataset

   Slicing and extracting the part of the data containing the specified year based on an array of given integer minutes.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   minute_array: :py:class:`list[int]`
       Minute(s) to be extracted.
   dim : :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_seconds_data(data_input: xarray.DataArray | xarray.Dataset, second_array: np.array(int) | List[int], dim: str = 'time') -> xarray.DataArray | xarray.Dataset

   Slicing and extracting the part of the data containing the specified year based on an array of given integer seconds.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   second_array: :py:class:`list[int]`
       Second(s) to be extracted.
   dim : :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_microseconds_data(data_input: xarray.DataArray | xarray.Dataset, microsecond_array: np.array(int) | List[int], dim: str = 'time') -> xarray.DataArray | xarray.Dataset

   Slicing and extracting the part of the data containing the specified year based on an array of given integer microseconds.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   microsecond_array: :py:class:`list[int]`
       Microsecond(s) to be extracted.
   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_nanoseconds_data(data_input: xarray.DataArray | xarray.Dataset, nanosecond_array: np.array(int) | List[int], dim: str = 'time') -> xarray.DataArray | xarray.Dataset

   Slicing and extracting the part of the data containing the specified year based on an array of given integer nanoseconds.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   nanosecond_array: :py:class:`list[int]`
       Nanosecond(s) to be extracted.
   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_dayofweek_data(data_input: xarray.DataArray | xarray.Dataset, dayofweek_array: np.array(int) | List[int], dim: str = 'time') -> xarray.DataArray | xarray.Dataset

   Slicing and extracting the part of the data containing the specified year based on an array of given integer dayofweek.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   dayofweek_array: :py:class:`list[int]`
       The days of the week to be extracted.

       The integer numbers correspond to the days of the week as follows.

   +-------------------+-------------------+
   | Day of the week   | Integer numbers   |
   +===================+===================+
   |      Monday       |         0         |
   +-------------------+-------------------+
   |      Tuesday      |         1         |
   +-------------------+-------------------+
   |      Wednesday    |         2         |
   +-------------------+-------------------+
   |      Thursday     |         3         |
   +-------------------+-------------------+
   |      Friday       |         4         |
   +-------------------+-------------------+
   |      Saturday     |         5         |
   +-------------------+-------------------+
   |      Sunday       |         6         |
   +-------------------+-------------------+

   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_yearmean_for_specific_months_data(data_input: xarray.DataArray | xarray.Dataset, month_array: np.array(int) | List[int], dim: str = 'time', **kwargs) -> xarray.DataArray | xarray.Dataset

   Get the annual average of certain months.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   month_array: :py:class:`list[int]`
       Month(s) to be extracted.
   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_year_exceed_index_upper_bound(data_input: xarray.DataArray, thresh: float, time_dim: str = 'time') -> numpy.array

   Extract the years under the specified threshold (upper bound) in the annual average index (one-dimensional data with only a `time` dimension).

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The one-dimensional data with only a `time` dimension.
   thresh: :py:class:`float <float>`.
       The threshold value.
   time_dim: :py:class:`str <str>`.
       The time coordinate dimension name.

   Returns
   -------
   :py:class:`numpy.array <numpy:numpy.array>`.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_basic_statistical_analysis.py


.. py:function:: get_year_exceed_index_lower_bound(data_input: xarray.DataArray, thresh: float, time_dim: str = 'time') -> numpy.array

   Extract the years under the specified threshold (lower bound) in the annual average index (one-dimensional data with only a `time` dimension).

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The one-dimensional data with only a `time` dimension.
   thresh: :py:class:`float <float>`.
       The threshold value.
   time_dim: :py:class:`str <str>`.
       The time coordinate dimension name.

   Returns
   -------
   :py:class:`numpy.array <numpy:numpy.array>`.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_basic_statistical_analysis.py


.. py:function:: get_time_exceed_index_upper_bound(data_input: xarray.DataArray, thresh: float, time_dim: str = 'time') -> numpy.array

   Extract the time under the specified threshold (upper bound) in the annual average index (one-dimensional data with only a `time` dimension).

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The one-dimensional data with only a `time` dimension.
   thresh: :py:class:`float <float>`.
       The threshold value.
   time_dim: :py:class:`str <str>`.
       The time coordinate dimension name.

   Returns
   -------
   Time array.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_da_bbo.py


.. py:function:: get_time_exceed_index_lower_bound(data_input: xarray.DataArray, thresh: float, time_dim: str = 'time') -> numpy.array

   Extract the time under the specified threshold (lower bound) in the annual average index (one-dimensional data with only a `time` dimension).

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The one-dimensional data with only a `time` dimension.
   thresh: :py:class:`float <float>`.
       The threshold value.
   time_dim: :py:class:`str <str>`.
       The time coordinate dimension name.

   Returns
   -------
   Time array.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_da_bbo.py


.. py:function:: open_muliti_dataset(files: str, dim: str, **kwargs) -> xarray.Dataset

   Open multiple netCDF files without the need for xarray's necessary dimension checks

   Parameters
   ----------
   - ver1: Version number 1
   - ver2: Version number 2

   Returns
   -------
   :py:class:`int <int>`.

   .. note::
       If `ver1<ver2`, return `-1`; If `ver1=ver2`, return `0`; If `ver1>ver2`, return `1`.

   Examples
   --------

   .. code:: python

       >>> import easyclimate as ecl
       >>> result = assert_compared_version("10.12.2.6.5", "10.12.2.6")
       >>> print(result)
       1

   .. note::
       - https://medium.com/pangeo/accessing-netcdf-and-grib-file-collections-as-cloud-native-virtual-datasets-using-kerchunk-625a2d0a9191
       - https://github.com/fsspec/kerchunk/issues/240


.. py:function:: calc_linregress_spatial(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', x: numpy.array = None, alternative: str = 'two-sided', returns_type: {'dataset_returns', 'dataset_vars'} = 'dataset_returns') -> xarray.Dataset | xarray.DataTree

   Calculate a linear least-squares regression (**trend**) for spatial data of time.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>` to be regression.
   dim: :py:class:`str <str>`, default `time`.
       Dimension(s) over which to apply linregress. By default linregress is applied over the `time` dimension.
   x: :py:class:`numpy.array <numpy.array>`
       Independent variable. If None, use `np.arange(len(data_input['time'].shape[0]))` instead.
   returns_type: :py:class:`str <str>`, default `'dataset_returns'`.
       Return data type.

   Returns
   -------
   result : ``LinregressResult`` Dataset
       The return Dataset have following data_var:

       **slope**: :py:class:`float <float>`
           Slope of the regression line.
       **intercept**: :py:class:`float <float>`
           Intercept of the regression line.
       **rvalue**: :py:class:`float <float>`
           The Pearson correlation coefficient. The square of ``rvalue``
           is equal to the coefficient of determination.
       **pvalue**: :py:class:`float <float>`
           The p-value for a hypothesis test whose null hypothesis is
           that the slope is zero, using Wald Test with t-distribution of
           the test statistic. See `alternative` above for alternative
           hypotheses.
       **stderr**: :py:class:`float <float>`
           Standard error of the estimated slope (gradient), under the
           assumption of residual normality.
       **intercept_stderr**: :py:class:`float <float>`
           Standard error of the estimated intercept, under the assumption
           of residual normality.

   .. seealso::
       :py:func:`scipy.stats.linregress <scipy:scipy.stats.linregress>`.


.. py:function:: calc_detrend_spatial(data_input: xarray.DataArray | xarray.Dataset, time_dim: str = 'time') -> xarray.DataArray | xarray.DataTree

   Remove linear trend along axis from data.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`
        The spatio-temporal data of :py:class:`xarray.DataArray<xarray.DataArray>` to be detrended.
   time_dim: :py:class:`str <str>`
       Dimension(s) over which to detrend. By default dimension is applied over the `time` dimension.

   Returns
   -------
   - :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::
       :py:func:`scipy.signal.detrend <scipy:scipy.signal.detrend>`.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_basic_statistical_analysis.py


.. py:function:: calc_corr_spatial(data_input: xarray.DataArray, x: xarray.DataArray | numpy.ndarray, time_dim: str = 'time', method: Literal['scipy', 'xarray'] = 'xarray') -> xarray.Dataset

   Calculate Pearson correlation coefficients and corresponding p-values between spatial data
   and a time series using ``scipy.stats.pearsonr``.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       Input spatial data with dimensions ``(time, ...)``.

       .. note::

           NaN values are automatically skipped in calculations.

   x : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`numpy.ndarray<numpy.ndarray>`
       Time series data with dimension ``(time,)``. Must have the same length as data_input's time dimension.

       .. note::

           NaN values are automatically skipped in calculations.

   time_dim: :py:class:`str <str>`
       Dimension(s) over which to detrend. By default dimension is applied over the `time` dimension.
   method : {'scipy', 'xarray'}, optional
       Method used to compute correlations:

       - 'scipy': Uses :py:func:`scipy.stats.pearsonr<scipy:scipy.stats.pearsonr>` for direct calculation
       - 'xarray': Uses xarray's built-in correlation with t-test conversion (faster)

       Default is 'xarray'.

   Returns
   -------
   reg_coeff, corr & pvalue (:py:class:`xarray.Dataset<xarray.Dataset>`)

   reg_coeff: :py:class:`xarray.DataArray<xarray.DataArray>`
       Regression coefficient, in units of ``data_input`` per standard deviation of the index.

   corr : :py:class:`xarray.DataArray<xarray.DataArray>`
       Pearson correlation coefficients with dimensions.
       Values range from -1 to 1 where:

       - 1: perfect positive correlation
       - -1: perfect negative correlation
       - 0: no correlation

   pvalue : :py:class:`xarray.DataArray<xarray.DataArray>`
       Two-tailed p-values with dimensions.
       Small p-values (<0.05) indicate statistically significant correlations.

   Examples
   --------
   >>> data_input = xr.DataArray(np.random.rand(10, 3, 4),
   ...                           dims=['time', 'lat', 'lon'],
   ...                           coords={'time': pd.date_range('2000-01-01', periods=10)})
   >>> x = xr.DataArray(np.random.rand(10), dims=['time'])
   >>> corr_dataset = ecl.calc_corr_spatial(data_input, x)

   .. seealso::
       :py:func:`scipy.stats.pearsonr<scipy:scipy.stats.pearsonr>`:
       The underlying correlation function used for calculations.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_corr_reg.py


.. py:function:: calc_leadlag_corr_spatial(data_input: xarray.DataArray, x: xarray.DataArray | numpy.ndarray, leadlag_array: numpy.array | List[int], time_dim: str = 'time', method: Literal['scipy', 'xarray'] = 'xarray')

   Calculate Pearson correlation coefficients and corresponding p-values between spatial data
   and a time series with specified lead or lag shifts, using ``scipy.stats.pearsonr`` or xarray methods.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       Input spatial data with dimensions ``(time, ...)`` representing spatial fields over time.

       .. note::
           NaN values are automatically skipped in calculations.

   x : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`numpy.ndarray<numpy.ndarray>`
       Time series data with dimension ``(time,)``. Must have the same length as ``data_input``'s time dimension.

       .. note::
           NaN values are automatically skipped in calculations.

   leadlag_array : :py:class:`numpy.ndarray<numpy.ndarray>` or :py:class:`List[int]<list>`
       Array or list of integers specifying the lead or lag shifts (in time steps) to apply to the time series `x`
       relative to `data_input`.

       - **Positive values** indicate a **lag**: the time series `x` is shifted forward in time (e.g., a value of +2 means `x` is delayed by 2 time steps relative to `data_input`).
       - **Negative values** indicate a **lead**: the time series `x` is shifted backward in time (e.g., a value of -2 means `x` leads `data_input` by 2 time steps).
       - A value of **0** means no shift (synchronous correlation).

       Example: If ``leadlag_array = [-2, 0, 2]``, correlations are computed for :math:`x` leading by 2 time steps, no shift, and lagging by 2 time steps, respectively.

   time_dim : :py:class:`str<str>`
       Name of the time dimension in `data_input` and `x`. Default is `"time"`.

   method : {'scipy', 'xarray'}, optional
       Method used to compute correlations:

       - `'scipy'`: Uses :py:func:`scipy.stats.pearsonr<scipy:scipy.stats.pearsonr>` for direct calculation, which may be more precise but slower.
       - `'xarray'`: Uses xarray's built-in correlation function with t-test conversion, which is typically faster.

       Default is `'xarray'`.

   Returns
   -------
   result : :py:class:`xarray.Dataset<xarray.Dataset>`
       Dataset containing two variables:

       - **corr** : :py:class:`xarray.DataArray<xarray.DataArray>`
           Pearson correlation coefficients with dimensions ``(leadlag, ...)``.
           Values range from -1 to 1, where:
           - 1 indicates a perfect positive correlation.
           - -1 indicates a perfect negative correlation.
           - 0 indicates no correlation.

       - **pvalue** : :py:class:`xarray.DataArray<xarray.DataArray>`
           Two-tailed p-values with dimensions ``(leadlag, ...)``.
           Small p-values (<0.05) indicate statistically significant correlations.

   Notes
   -----
   - The function iterates over each lead/lag value in `leadlag_array`, computes the correlation between the shifted `x` and `data_input`, and concatenates results along a new `leadlag` dimension.
   - Shifting `x` may introduce NaN values at the edges of the time series, which are handled automatically during correlation calculations.
   - Ensure `data_input` and `x` have compatible time dimensions to avoid errors.

   Examples
   --------
   >>> import xarray as xr
   >>> import numpy as np
   >>> data = xr.DataArray(np.random.rand(100, 10, 10), dims=["time", "lat", "lon"])
   >>> ts = xr.DataArray(np.random.rand(100), dims=["time"])
   >>> leadlag = [-2, 0, 2]
   >>> result = calc_leadlag_corr_spatial(data, ts, leadlag, time_dim="time", method="xarray")
   >>> print(result)
   Processing leadlag: 2 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
   <xarray.Dataset> Size: 7kB
   Dimensions:    (leadlag: 3, lat: 10, lon: 10)
   Coordinates:
   * leadlag    (leadlag) int64 24B -2 0 2
   Dimensions without coordinates: lat, lon
   Data variables:
       reg_coeff  (leadlag, lat, lon) float64 2kB 0.006322 0.002647 ... -0.02781
       corr       (leadlag, lat, lon) float64 2kB 0.02141 0.00894 ... -0.09169
       pvalue     (leadlag, lat, lon) float64 2kB 0.8326 0.9297 ... 0.3053 0.3643


.. py:function:: calc_multiple_linear_regression_spatial(y_data: xarray.DataArray, x_datas: List[xarray.DataArray], dim='time') -> xarray.Dataset

   Apply multiple linear regression to dataset across spatial dimensions.

   .. math::

       y = a_1 x_1 + a_2 x_2 + \cdots

   Parameters
   -----------
   y_data : :py:class:`xarray.DataArray<xarray.DataArray>`
       Dependent variable with dimensions, each with dimensions ``(time,)``.
   x_datas : :py:class:`list <list>` of :py:class:`xarray.DataArray<xarray.DataArray>`
       List of independent variables, each with dimensions ``(time,)``.
   dim : :py:class:`str <str>`, optional
       Time dimension name (default: ``'time'``)

   Returns
   --------
   :py:class:`xarray.Dataset <xarray.Dataset>`
       :py:class:`xarray.Dataset <xarray.Dataset>` containing regression results with:

       - slopes: slope coefficients for each predictor ``(coef, lat, lon)``
       - intercept: intercept values ``(lat, lon)``
       - r_squared: coefficient of determination ``(lat, lon)``
       - slopes_p: p-values for slope coefficients ``(coef, lat, lon)``
       - intercept_p: p-values for intercept ``(lat, lon)``

   Raises
   -------
   ValueError
       If the time coordinates of input variables don't match.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_multi_linear_reg.py


.. py:function:: calc_ttestSpatialPattern_spatial(data_input1: xarray.DataArray, data_input2: xarray.DataArray, dim: str = 'time', equal_var: bool = True, alternative: Literal['two-sided', 'less', 'greater'] = 'two-sided', method: Literal['scipy', 'xarray'] = 'xarray') -> xarray.Dataset

   Calculate the T-test for the means of two independent sptial samples along with other axis (i.e. 'time') of scores.

   Parameters
   ----------
   data_input1: :py:class:`xarray.DataArray<xarray.DataArray>`
        The first spatio-temporal data of xarray DataArray to be calculated.
   data_input2: :py:class:`xarray.DataArray<xarray.DataArray>`
        The second spatio-temporal data of xarray DataArray to be calculated.

   .. note::
       - The order of `data_input1` and `data_input2` has no effect on the calculation result.
       - The non-time dimensions of the two data sets must be exactly the same, and the dimensionality values must be arranged in the same order (ascending or descending).

   dim: :py:class:`str <str>`
       Dimension(s) over which to apply the test. By default the test is applied over the `time` dimension.

   equal_var: :py:class:`bool <bool>`
       If True (default), perform a standard independent 2 sample test that assumes equal population variances (see https://en.wikipedia.org/wiki/T-test#Independent_two-sample_t-test).
       If False, perform Welch’s t-test, which does not assume equal population variance (see https://en.wikipedia.org/wiki/Welch%27s_t-test).

   alternative : {'two-sided', 'less', 'greater'}, optional
       Defines the alternative hypothesis.
       The following options are available (default is 'two-sided'):

       - 'two-sided': the means of the distributions underlying the samples are unequal.
       - 'less': the mean of the distribution underlying the first sample is less than the mean of the distribution underlying the second sample.
       - 'greater': the mean of the distribution underlying the first sample is greater than the mean of the distribution underlying the second sample.

   method : {'scipy', 'xarray'}, optional
       Method used to compute correlations:

       - 'scipy': Uses :py:func:`scipy.stats.ttest_ind<scipy:scipy.stats.ttest_ind>` for direct calculation
       - 'xarray': Uses xarray's built-in method to calculate (faster)

       Default is 'xarray'.

   Returns
   -------
   - **statistic**, **pvalue**: :py:class:`xarray.Dataset<xarray.Dataset>`.

   .. seealso::
       :py:func:`scipy.stats.ttest_ind <scipy:scipy.stats.ttest_ind>`.


.. py:function:: calc_windmask_ttestSpatialPattern_spatial(data_input1: xarray.Dataset, data_input2: xarray.Dataset, dim: str = 'time', u_dim: str = 'u', v_dim: str = 'v', mask_method: Literal['or', 'and'] = 'or', thresh: float = 0.05, equal_var: bool = True, alternative: Literal['two-sided', 'less', 'greater'] = 'two-sided', method: Literal['scipy', 'xarray'] = 'xarray')

   Generate a significance mask for T-tests on the means of two independent spatial zonal (u) and meridional (v) wind samples,
   aggregated over the specified dimension (default 'time').

   Parameters
   ----------
   data_input1 : :py:class:`xarray.Dataset`
        The first spatio-temporal data of xarray Dataset to be calculated. It is necessary to include the zonal wind component (u_dim) and the meridional wind component (v_dim).
   data_input2 : :py:class:`xarray.Dataset`
        The second spatio-temporal data of xarray Dataset to be calculated. It is necessary to include the zonal wind component (u_dim) and the meridional wind component (v_dim).

   .. note::
       - The order of `data_input1` and `data_input2` has no effect on the calculation result.
       - The non-time dimensions of the two data sets must be exactly the same, and the dimensionality values must be arranged in the same order (ascending or descending).

   dim : :py:class:`str`, default: `time`
       Dimension(s) over which to apply the test. By default the test is applied over the `time` dimension.
   u_dim : :py:class:`str`, default: `u`
       Variable name for the u velocity (zonal wind, in x direction).
   v_dim : :py:class:`str`, default: `v`
       Variable name for the v velocity (meridional wind, in y direction).
   mask_method : Literal["or", "and"], default: "or"
       Method to combine the significance masks for u and v components:

       - "or": A grid point is considered significant if either the u or v component is significant (p <= thresh).
       - "and": A grid point is considered significant if both the u and v components are significant (p <= thresh).

   thresh : :py:class:`float`, default: 0.05
       The significance level (alpha) for the p-value threshold used to create the mask.
   equal_var : :py:class:`bool`, default: True
       If True (default), perform a standard independent 2 sample test that assumes equal population variances (see https://en.wikipedia.org/wiki/T-test#Independent_two-sample_t-test).
       If False, perform Welch’s t-test, which does not assume equal population variance (see https://en.wikipedia.org/wiki/Welch%27s_t-test).

   alternative : {'two-sided', 'less', 'greater'}, optional=
       Defines the alternative hypothesis.
       The following options are available (default is 'two-sided'):

       - 'two-sided': the means of the distributions underlying the samples are unequal.
       - 'less': the mean of the distribution underlying the first sample is less than the mean of the distribution underlying the second sample.
       - 'greater': the mean of the distribution underlying the first sample is greater than the mean of the distribution underlying the second sample.

   method : {'scipy', 'xarray'}, optional
       Method used to compute t-tests:

       - 'scipy': Uses :py:func:`scipy.stats.ttest_ind<scipy:scipy.stats.ttest_ind>` for direct calculation
       - 'xarray': Uses xarray's built-in method to calculate (faster)

       Default is 'xarray'.

   Returns
   -------
   masked_pvalue : :py:class:`xarray.DataArray`
       A boolean mask indicating significant regions (True where p <= thresh, combined via mask_method for u and v).

   .. seealso::
       :py:func:`scipy.stats.ttest_ind <scipy:scipy.stats.ttest_ind>`.


.. py:function:: calc_levenetestSpatialPattern_spatial(data_input1: xarray.DataArray, data_input2: xarray.DataArray, dim: str = 'time', center: {'mean', 'median', 'trimmed'} = 'median', proportiontocut: float = 0.05) -> xarray.Dataset

   Perform Levene test for equal variances of two independent sptial samples along with other axis (i.e. 'time') of scores.

   The Levene test tests the null hypothesis that all input samples are from populations with equal variances.
   Levene's test is an alternative to Bartlett's test in the case where there are significant deviations from normality.

   Parameters
   ----------
   data_input1: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The first spatio-temporal data of xarray DataArray to be calculated.
   data_input2: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The second spatio-temporal data of xarray DataArray to be calculated.

   .. note::
       - The order of `data_input1` and `data_input2` has no effect on the calculation result.
       - The non-time dimensions of the two data sets must be exactly the same, and the dimensionality values must be arranged in the same order (ascending or descending).

   dim: :py:class:`str <str>`.
       Dimension(s) over which to apply the test. By default the test is applied over the `time` dimension.
   center: {'mean', 'median', 'trimmed'}, default `'median'`.
       Which function of the data to use in the test.

       .. note::

           Three variations of Levene’s test are possible. The possibilities and their recommended usages are:

           - median: Recommended for skewed (non-normal) distributions.
           - mean: Recommended for symmetric, moderate-tailed distributions.
           - trimmed: Recommended for heavy-tailed distributions.

           The test version using the mean was proposed in the original article of Levene (Levene, H., 1960) while the median and trimmed mean have been studied by Brown and Forsythe (Brown, M. B. and Forsythe, A. B., 1974), sometimes also referred to as Brown-Forsythe test.


   proportiontocut: :py:class:`float <float>`, default `0.05`.
       When center is `'trimmed'`, this gives the proportion of data points to cut from each end (See :py:func:`scipy.stats.trim_mean <scipy:scipy.stats.trim_mean>`).

   Returns
   -------
   - **statistic**, **pvalue**: :py:class:`xarray.Dataset<xarray.Dataset>`.

   Reference
   --------------
   - Levene, H. (1960). In Contributions to Probability and Statistics: Essays in Honor of Harold Hotelling, I. Olkin et al. eds., Stanford University Press, pp. 278-292.
   - Morton B. Brown & Alan B. Forsythe (1974) Robust Tests for the Equality of Variances, Journal of the American Statistical Association, 69:346, 364-367, DOI: https://doi.org/10.1080/01621459.1974.10482955

   .. seealso::
       :py:func:`scipy.stats.levene <scipy:scipy.stats.levene>`.


.. py:function:: calc_skewness_spatial(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time') -> xarray.Dataset | xarray.DataTree

   Calculate the skewness of the spatial field on the time axis and its significance test.

   The :math:`k` th statistical moment about the mean is given by

   .. math::
       m_k = \sum_{i=1}^{N} \frac{(x_i-\bar{x})^k}{N}

   where :math:`x_i` is the :math:`i` th observation, :math:`\bar{x}` the mean and :math:`N` the number of observations.

   One definition of the coefficient of skewness is

   .. math::
       a_3 = \frac{m_3}{(m_2)^{3/2}}

   Skewness is a measure of the asymmetry of a distribution and is zero for a normal distribution. If the longer wing of a distribution
   occurs for values of :math:`x` higher than the mean, that distribution is said to have positive skewness. If thelonger wing occurs for
   values of :math:`x` lower than the mean, the distribution is said to have negative skewness.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`
        The spatio-temporal data of xarray DataArray to be calculated.
   dim: :py:class:`str <str>`
       Dimension(s) over which to apply skewness. By default skewness is applied over the `time` dimension.

   Returns
   -------
   - **skewness**, **pvalue**: :py:class:`xarray.Dataset<xarray.Dataset>`.

   Reference
   --------------
   White, G. H. (1980). Skewness, Kurtosis and Extreme Values of
   Northern Hemisphere Geopotential Heights, Monthly Weather Review, 108(9), 1446-1455.
   Website: https://journals.ametsoc.org/view/journals/mwre/108/9/1520-0493_1980_108_1446_skaevo_2_0_co_2.xml

   .. seealso::
       :py:func:`scipy.stats.skew <scipy:scipy.stats.skew>`, :py:func:`scipy.stats.normaltest <scipy:scipy.stats.normaltest>`.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_basic_statistical_analysis.py


.. py:function:: calc_kurtosis_spatial(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time') -> xarray.DataArray | xarray.DataTree

   Calculate the kurtosis of the spatial field on the time axis and its significance test.

   The :math:`k` th statistical moment about the mean is given by

   .. math::
       m_k = \sum_{i=1}^{N} \frac{(x_i-\bar{x})^k}{N}

   where :math:`x_i` is the :math:`i` th observation, :math:`\bar{x}` the mean and :math:`N` the number of observations.

   The coefficient of kurtosis is defined by

   .. math::
       a_4 = \frac{m_4}{(m_2)^{2}}

   The kurtosis of a normal distribution is 3. If a distribution has a large central region which is flatter than a normal distribution
   with the same mean and variance, it has a kurtosis of less than 3. If the distribution has a central maximum more peaked and with
   longer wings than the equivalent normal distribution, its kurtosis is higher than 3 (Brooks and Carruthers, 1954).
   Extreme departures from the mean will cause very high values of kurtosis. Consequently, high kurtosis has been used as
   an indicator of bad data (Craddock and Flood, 1969). For the same reason, high values of kurtosis can be a result of one or two
   extreme events in a period of several years.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`
        The spatio-temporal data of xarray DataArray to be calculated.
   dim: :py:class:`str <str>`
       Dimension(s) over which to apply kurtosis. By default kurtosis is applied over the `time` dimension.

   Returns
   -------
   - kurtosis: :py:class:`xarray.DataArray<xarray.DataArray>`.

   Reference
   --------------
   White, G. H. (1980). Skewness, Kurtosis and Extreme Values of
   Northern Hemisphere Geopotential Heights, Monthly Weather Review, 108(9), 1446-1455.
   Website: https://journals.ametsoc.org/view/journals/mwre/108/9/1520-0493_1980_108_1446_skaevo_2_0_co_2.xml

   Køie, M., Brooks, C.E., & Carruthers, N. (1954). Handbook of Statistical Methods in Meteorology. Oikos, 4, 202.

   Craddock, J.M. and Flood, C.R. (1969), Eigenvectors for representing the 500 mb geopotential
   surface over the Northern Hemisphere. Q.J.R. Meteorol. Soc., 95: 576-593.
   doi: https://doi.org/10.1002/qj.49709540510

   .. seealso::
       :py:func:`scipy.stats.kurtosis <scipy:scipy.stats.kurtosis>`.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_basic_statistical_analysis.py


.. py:function:: calc_theilslopes_spatial(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', x=None, alpha: float = 0.95, method: {'joint', 'separate'} = 'separate', returns_type: {'dataset_returns', 'dataset_vars'} = 'dataset_returns') -> xarray.Dataset | xarray.DataTree

   Computes the Theil-Sen estimator.

   Theilslopes implements a method for robust linear regression. It computes the slope as the median of all slopes between paired values.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>` to be regression.
   dim: :py:class:`str <str>`, default `time`.
       Dimension(s) over which to apply linregress. By default linregress is applied over the `time` dimension.
   x: :py:class:`numpy.array <numpy.array>`
       Independent variable. If None, use `np.arange(len(data_input['time'].shape[0]))` instead.
   alpha: :py:class:`float <float>`, default 0.95.
       Confidence degree between 0 and 1. Default is 95% confidence. Note that alpha is symmetric around 0.5, i.e. both 0.1 and 0.9 are interpreted as "find the 90% confidence interval".
   method: {'joint', 'separate'}, default `'separate'`.
       Method to be used for computing estimate for intercept. Following methods are supported,

       - *joint*: Uses `np.median(y - slope * x)` as intercept.
       - *separate*: Uses `np.median(y) - slope * np.median(x)` as intercept.

   returns_type: :py:class:`str <str>`, default `'dataset_returns'`.
       Return data type.

   Returns
   -------
   result : ``TheilslopesResult`` Dataset
       The return Dataset have following data_var:

       **slope**: :py:class:`float <float>`
           Theil slope.
       **intercept**: :py:class:`float <float>`
           Intercept of the Theil line.
       **low_slope**: :py:class:`float <float>`
           Lower bound of the confidence interval on `slope`.
       **high_slope**: :py:class:`float <float>`
           Upper bound of the confidence interval on `slope`.

   .. seealso::
       :py:func:`scipy.stats.theilslopes <scipy:scipy.stats.theilslopes>`.


.. py:function:: calc_lead_lag_correlation_coefficients(pcs: dict, pairs: List[tuple], max_lag: int) -> xarray.Dataset

   Compute lead-lag correlation coefficients for specified pairs of indexes.

   This function calculates the cross-correlation between pairs of time series (e.g., MJO/BSISO principal components PC1 vs. PC2)
   to determine their lead-lag relationships. The correlation coefficients are computed for a range
   of lags, and the maximum correlation and corresponding lag are stored as attributes in the output
   dataset.

   Parameters
   ----------
   pcs : py:class:`dict <dict>`
       Dictionary mapping PC names (e.g., 'PC1', 'PC2') to :py:class:`xarray.DataArray<xarray.DataArray>` objects, where each
       DataArray represents a principal component time series with a ``'time'`` dimension.
   pairs : py:class:`list <list>` of tuples
       List of tuples, where each tuple contains (pair_name, pc_name1, pc_name2). For example,
       [('PC1_vs_PC2', 'PC1', 'PC2'), ('PC3_vs_PC4', 'PC3', 'PC4')].
   max_lag : py:class:`int <int>`
       Maximum lag (in time steps) to consider for the cross-correlation. The function computes
       correlations for lags from ``-max_lag`` to ``+max_lag``.

   Returns
   -------
   corr_da : :py:class:`xarray.Dataset<xarray.Dataset>`
       Dataset containing correlation coefficients for each pair, with a 'lag' dimension.
       Each variable (e.g., 'PC1_vs_PC2') has attributes ``'max_correlation' (float)`` and
       ``'lag_at_max_correlation' (int)`` indicating the maximum correlation and the lag at which it occurs.

   Notes
   -----

   - Positive lags indicate that the first PC leads the second; negative lags indicate the opposite.
   - The input PCs should have no missing values. Use ``fillna`` or ``interpolate_na`` if needed.
   - The correlation coefficients are normalized to range between :math:`-1` and :math:`1`.


.. py:function:: calc_timeseries_correlations(da: dict[str, xarray.DataArray] | list[xarray.DataArray], dim: str = 'time') -> xarray.DataArray

   Calculate the correlation matrix between multiple DataArray time series.

   This function calculates pairwise correlations between time series in the input DataArrays
   using the specified correlation method along the given dimension. The output is a symmetric
   correlation matrix stored as an xarray DataArray with dimensions (var1, var2).

   Parameters
   ----------
   da : :py:class:`dict[str, xarray.DataArray<xarray.DataArray>]` or :py:class:`list[xarray.DataArray<xarray.DataArray>]`.
       A dictionary with names as keys and DataArrays as values, or a list of DataArrays.
       Each DataArray must contain the specified dimension.
   dim : :py:class:`str <str>`, default: `time`.
       The dimension along which to compute correlations. All DataArrays must have this dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.
       A DataArray containing the correlation matrix with dimensions (var1, var2).
       Coordinates are set to the names of the input time series.

   Examples
   --------
   >>> time = pd.date_range('2020-01-01', '2020-12-31', freq='D')
   >>> da1 = xr.DataArray(np.random.randn(len(time)), dims='time', coords={'time': time}, name='series1')
   >>> da2 = xr.DataArray(da1 * 0.5 + np.random.randn(len(time)) * 0.5, dims='time', coords={'time': time}, name='series2')
   >>> data = {'series1': da1, 'series2': da2}
   >>> corr_matrix = calc_timeseries_correlations(data, method='pearson')
   >>> print(corr_matrix)


.. py:function:: calc_non_centered_corr(data_input1, data_input2, dim=None)

   Compute the non-centered (uncentered) correlation coefficient between two xarray DataArrays.
   This is equivalent to the cosine similarity, calculated as the sum of the product of the two arrays
   divided by the product of their L2 norms (Euclidean norms), without subtracting the means.

   The formula is:

   .. math::

       r = \frac{\sum (x \cdot y)}{\sqrt{\sum x^2} \cdot \sqrt{\sum y^2}}

   Parameters
   ----------
   data_input1 : :py:class:`xarray.DataArray`
       The first input data array to be correlated.
   data_input2 : :py:class:`xarray.DataArray`
       The second input data array to be correlated.

   .. note::
       - Both inputs must be xarray DataArray objects.
       - The arrays must have compatible shapes: if `dim` is specified, it must be a shared dimension;
         if `dim` is None, all dimensions are flattened into a single vector for computation.
       - The result is set to 0 where the denominator (product of norms) is zero to avoid division by zero.

   dim : :py:class:`str` or None, optional
       Dimension(s) over which to compute the correlation. If None (default), the arrays are flattened
       across all dimensions into a single vector before computation. If a string, the correlation is
       computed along the specified dimension, preserving other dimensions.

   Returns
   -------
   corr : :py:class:`xarray.DataArray`
       The non-centered correlation coefficient, with the same dimensions as the input arrays
       (or broadcasted appropriately).

   .. seealso::
       :py:func:`scipy.spatial.distance.cosine <scipy:scipy.spatial.distance.cosine>`
       (for the related cosine distance metric).

   Examples
   --------
   Compute correlation along a dimension:

   >>> import xarray as xr
   >>> import numpy as np
   >>> import easyclimate as ecl
   >>> da1 = xr.DataArray(np.array([[1, 2], [3, 4]]), dims=['x', 'y'])
   >>> da2 = xr.DataArray(np.array([[2, 3], [4, 5]]), dims=['x', 'y'])
   >>> corr = ecl.calc_non_centered_corr(da1, da2, dim='y')
   >>> print(corr)
   <xarray.DataArray (x: 2)> Size: 16B
   array([0.99227788, 0.99951208])

   Flatten and compute scalar correlation:

   >>> corr_flat = calc_non_centered_corr(da1, da2)
   >>> print(corr_flat)
   Dimensions without coordinates: x
   <xarray.DataArray ()> Size: 8B
   array(0.99380799)


.. py:function:: calc_pattern_corr(data_input1: xarray.DataArray, data_input2: xarray.DataArray, time_dim: str = 'time')

   Compute the pattern correlation (non-centered) between two xarray DataArrays over their common
   spatial dimensions. This is useful for comparing spatial patterns, such as in climate data.

   It uses the non-centered correlation formula:

   .. math::

       r = \frac{\sum (x \cdot y)}{\sqrt{\sum x^2} \cdot \sqrt{\sum y^2}}

   where the summation is over the stacked spatial (pattern) dimensions.

   The spatial pattern dimensions are automatically detected as the intersection of the input dimensions,
   excluding 'time' (if present). Both inputs are stacked along these pattern dimensions into a temporary
   'pattern' dimension, and the non-centered correlation is computed along it.

   - If both inputs lack 'time', the result is a scalar.
   - If one input has 'time' and the other does not, the result preserves the 'time' dimension from the timed input.
   - Broadcasting occurs automatically for compatible shapes.

   Parameters
   ----------
   data_input1 : :py:class:`xarray.DataArray`
       The first input data array (e.g., spatial pattern or time series of patterns).
   data_input2 : :py:class:`xarray.DataArray`
       The second input data array (must have compatible spatial dimensions).
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name.

   Returns
   -------
   corr : :py:class:`xarray.DataArray` or scalar
       The pattern correlation coefficient(s). Dimensions match the non-spatial dimensions of the inputs
       (e.g., 'time' if present in one input).

   .. note::
       - Assumes inputs have compatible shapes and the only differing dimension is 'time'.
       - Equivalent to cosine similarity over the spatial pattern.
       - For zero-norm cases, correlation is set to 0.

   .. seealso::
       - `pattern_cor -NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/pattern_cor.shtml>`__
       - :py:func:`calc_non_centered_corr`

   Examples
   --------
   Scalar correlation between two spatial patterns:

   >>> import xarray as xr
   >>> import numpy as np
   >>> import easyclimate as ecl
   >>> # Create a random number generator with a fixed seed.
   >>> rng = np.random.default_rng(42)
   >>> pat1 = xr.DataArray(rng.random((2, 3)), dims=['lat', 'lon'])
   >>> pat2 = xr.DataArray(rng.random((2, 3)), dims=['lat', 'lon'])
   >>> pcc = ecl.calc_pattern_corr(pat1, pat2)
   >>> print(pcc)
   <xarray.DataArray 'pcc' ()> Size: 8B
   array(0.85730639)
   Attributes:
       long_name:  Pattern Correlation Coefficient (non-centered)
       units:      dimensionless

   Time series correlation (one with time):

   >>> # Create a random number generator with a fixed seed.
   >>> rng = np.random.default_rng(42)
   >>> time = xr.DataArray(np.arange(4), dims=['time'])
   >>> timed_pat = xr.DataArray(rng.random((4, 2, 3)), dims=['time', 'lat', 'lon'])
   >>> pcc_time = ecl.calc_pattern_corr(timed_pat, pat2)
   >>> print(pcc_time)
   <xarray.DataArray 'pcc' (time: 4)> Size: 32B
   array([0.85730639, 1.        , 0.78188174, 0.88162673])
   Dimensions without coordinates: time
   Attributes:
       long_name:  Pattern Correlation Coefficient (non-centered)
       units:      dimensionless


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


.. py:function:: calc_seasonal_mean(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', extract_season: Literal['DJF', 'MAM', 'JJA', 'SON', None] = None, **kwargs) -> xarray.DataArray

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


.. py:function:: calc_yearly_climatological_mean(data_input: xarray.DataArray, dim: str = 'time', **kwargs)

   Calculate yearly mean.

   For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same year it is:

   .. math::
       o(t, x) = \mathrm{mean} \left \lbrace i(t', x), t_1 < t' \leqslant t_n \right\rbrace

   .. tip::
       This function uses :py:meth:`xarray.DataArray.groupby<xarray:xarray.DataArray.groupby>` to implement the calculation.
       To substantially improve the performance of GroupBy operations, particularly with dask `install the flox package <https://flox.readthedocs.io/>`_.
       flox `extends Xarray's in-built GroupBy capabilities <https://flox.readthedocs.io/en/latest/xarray.html>`_ by allowing
       grouping by multiple variables, and lazy grouping by dask arrays. If installed, Xarray will automatically use flox by default.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

       .. note::

           The recommended frequence of the `data_input` is monthly.

   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` with `time` dimension type of :py:class:`numpy.datetime64<numpy.datetime64>`.

   .. seealso::
       :py:func:`numpy.mean <numpy:numpy.mean>`, :py:func:`dask.array.mean <dask:dask.array.mean>`,
       :py:meth:`xarray.DataArray.mean <xarray:xarray.DataArray.mean>`,
       :py:meth:`xarray.core.groupby.DataArrayGroupBy.mean <xarray:xarray.core.groupby.DataArrayGroupBy.mean>`.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_basic_statistical_analysis.py


.. py:function:: calc_yearly_climatological_sum(data_input: xarray.DataArray, dim: str = 'time', **kwargs)

   Calculate yearly sum.

   For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same year it is:

   .. math::
       o(t, x) = \mathrm{sum} \left \lbrace i(t', x), t_1 < t' \leqslant t_n \right\rbrace

   .. tip::
       This function uses :py:meth:`xarray.DataArray.groupby<xarray:xarray.DataArray.groupby>` to implement the calculation.
       To substantially improve the performance of GroupBy operations, particularly with dask `install the flox package <https://flox.readthedocs.io/>`_.
       flox `extends Xarray's in-built GroupBy capabilities <https://flox.readthedocs.io/en/latest/xarray.html>`_ by allowing
       grouping by multiple variables, and lazy grouping by dask arrays. If installed, Xarray will automatically use flox by default.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

       .. note::

           The recommended frequence of the `data_input` is monthly.

   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating sum on this object's data.
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` with `time` dimension type of :py:class:`numpy.datetime64<numpy.datetime64>`.

   .. seealso::
       :py:func:`numpy.sum <numpy:numpy.sum>`, :py:func:`dask.array.sum <dask:dask.array.sum>`,
       :py:meth:`xarray.DataArray.sum <xarray:xarray.DataArray.sum>`,
       :py:meth:`xarray.core.groupby.DataArrayGroupBy.sum <xarray:xarray.core.groupby.DataArrayGroupBy.sum>`.


.. py:function:: calc_yearly_climatological_std(data_input: xarray.DataArray, dim: str = 'time', **kwargs)

   Calculate yearly standard deviation.

   For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same year it is:

   .. math::
       o(t, x) = \mathrm{std} \left \lbrace i(t', x), t_1 < t' \leqslant t_n \right\rbrace

   .. tip::
       This function uses :py:meth:`xarray.DataArray.groupby<xarray:xarray.DataArray.groupby>` to implement the calculation.
       To substantially improve the performance of GroupBy operations, particularly with dask `install the flox package <https://flox.readthedocs.io/>`_.
       flox `extends Xarray's in-built GroupBy capabilities <https://flox.readthedocs.io/en/latest/xarray.html>`_ by allowing
       grouping by multiple variables, and lazy grouping by dask arrays. If installed, Xarray will automatically use flox by default.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

       .. note::

           The recommended frequence of the `data_input` is monthly.

   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating std on this object's data.
       These could include dask-specific kwargs like split_every.

   .. note::
       The parameter `ddof` is `Delta Degrees of Freedom`: the divisor used in the calculation is `N - ddof`,
       where `N` represents the number of elements. If the data needs to be Normalize by `(n-1)`, then `ddof=1`.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` with `time` dimension type of :py:class:`numpy.datetime64<numpy.datetime64>`.

   .. seealso::
       :py:func:`numpy.std <numpy:numpy.std>`, :py:func:`dask.array.std <dask:dask.array.std>`,
       :py:meth:`xarray.DataArray.std <xarray:xarray.DataArray.std>`,
       :py:meth:`xarray.core.groupby.DataArrayGroupBy.std <xarray:xarray.core.groupby.DataArrayGroupBy.std>`.


.. py:function:: calc_yearly_climatological_var(data_input: xarray.DataArray, dim: str = 'time', **kwargs)

   Calculate yearly standard deviation.

   For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same year it is:

   .. math::
       o(t, x) = \mathrm{var} \left \lbrace i(t', x), t_1 < t' \leqslant t_n \right\rbrace

   .. tip::
       This function uses :py:meth:`xarray.DataArray.groupby<xarray:xarray.DataArray.groupby>` to implement the calculation.
       To substantially improve the performance of GroupBy operations, particularly with dask `install the flox package <https://flox.readthedocs.io/>`_.
       flox `extends Xarray's in-built GroupBy capabilities <https://flox.readthedocs.io/en/latest/xarray.html>`_ by allowing
       grouping by multiple variables, and lazy grouping by dask arrays. If installed, Xarray will automatically use flox by default.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

       .. note::

           The recommended frequence of the `data_input` is monthly.

   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating var on this object's data.
       These could include dask-specific kwargs like split_every.

   .. note::
       The parameter `ddof` is `Delta Degrees of Freedom`: the divisor used in the calculation is `N - ddof`,
       where `N` represents the number of elements. If the data needs to be Normalize by `(n-1)`, then `ddof=1`.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` with `time` dimension type of :py:class:`numpy.datetime64<numpy.datetime64>`.

   .. seealso::
       :py:func:`numpy.var <numpy:numpy.var>`, :py:func:`dask.array.var <dask:dask.array.var>`,
       :py:meth:`xarray.DataArray.var <xarray:xarray.DataArray.var>`,
       :py:meth:`xarray.core.groupby.DataArrayGroupBy.var <xarray:xarray.core.groupby.DataArrayGroupBy.var>`.


.. py:function:: calc_yearly_climatological_max(data_input: xarray.DataArray, dim: str = 'time', **kwargs)

   Calculate yearly standard deviation.

   For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same year it is:

   .. math::
       o(t, x) = \mathrm{max} \left \lbrace i(t', x), t_1 < t' \leqslant t_n \right\rbrace

   .. tip::
       This function uses :py:meth:`xarray.DataArray.groupby<xarray:xarray.DataArray.groupby>` to implement the calculation.
       To substantially improve the performance of GroupBy operations, particularly with dask `install the flox package <https://flox.readthedocs.io/>`_.
       flox `extends Xarray's in-built GroupBy capabilities <https://flox.readthedocs.io/en/latest/xarray.html>`_ by allowing
       grouping by multiple variables, and lazy grouping by dask arrays. If installed, Xarray will automatically use flox by default.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

       .. note::

           The recommended frequence of the `data_input` is monthly.

   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating max on this object's data.
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` with `time` dimension type of :py:class:`numpy.datetime64<numpy.datetime64>`.

   .. seealso::
       :py:func:`numpy.maximum <numpy:numpy.maximum>`, :py:func:`dask.array.max <dask:dask.array.max>`,
       :py:meth:`xarray.DataArray.max <xarray:xarray.DataArray.max>`,
       :py:meth:`xarray.core.groupby.DataArrayGroupBy.max <xarray:xarray.core.groupby.DataArrayGroupBy.max>`.


.. py:function:: calc_yearly_climatological_min(data_input: xarray.DataArray, dim: str = 'time', **kwargs)

   Calculate yearly standard deviation.

   For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same year it is:

   .. math::
       o(t, x) = \mathrm{min} \left \lbrace i(t', x), t_1 < t' \leqslant t_n \right\rbrace

   .. tip::
       This function uses :py:meth:`xarray.DataArray.groupby<xarray:xarray.DataArray.groupby>` to implement the calculation.
       To substantially improve the performance of GroupBy operations, particularly with dask `install the flox package <https://flox.readthedocs.io/>`_.
       flox `extends Xarray's in-built GroupBy capabilities <https://flox.readthedocs.io/en/latest/xarray.html>`_ by allowing
       grouping by multiple variables, and lazy grouping by dask arrays. If installed, Xarray will automatically use flox by default.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

       .. note::

           The recommended frequence of the `data_input` is monthly.

   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating min on this object's data.
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` with `time` dimension type of :py:class:`numpy.datetime64<numpy.datetime64>`.

   .. seealso::
       :py:func:`numpy.minimum <numpy:numpy.minimum>`, :py:func:`dask.array.min <dask:dask.array.min>`,
       :py:meth:`xarray.DataArray.min <xarray:xarray.DataArray.min>`,
       :py:meth:`xarray.core.groupby.DataArrayGroupBy.min <xarray:xarray.core.groupby.DataArrayGroupBy.min>`.


.. py:function:: open_tutorial_dataset(name: str, cache: bool = True, cache_dir: None | str | os.PathLike = None, progressbar: bool = False, *, engine: xarray.backends.api.T_Engine = None, **kws) -> xarray.Dataset

   Open a dataset from the online repository (requires internet).

   If a local copy is found then always use that to avoid network traffic.

   Available datasets:

   * ``"air_202201_mon_mean"``: 2m air temperature of the NCEP reanalysis subset
   * ``"hgt_202201_mon_mean"``: Geopotential height of the NCEP reanalysis subset
   * ``"precip_202201_mon_mean"``: Precipitation of the NCEP reanalysis subset
   * ``"pressfc_202201_mon_mean"``: Mean sea surface pressure of the NCEP reanalysis subset
   * ``"shum_202201_mon_mean"``: Absolute humidity of the NCEP reanalysis subset
   * ``"uwnd_202201_mon_mean"``: Zonal wind of the NCEP reanalysis subset
   * ``"vwnd_202201_mon_mean"``: Meridional wind of the NCEP reanalysis subset
   * ``"omega_202201_mon_mean"``: Vertical velocity of the NCEP reanalysis subset
   * ``"mini_HadISST_ice"``: Hadley Centre Sea Ice and Sea Surface Temperature data set (HadISST) subset
   * ``"PressQFF_202007271200_872"``: Observational data from European stations (from https://github.com/EXCITED-CO2/xarray-regrid)
   * ``"pr_wtr_eatm_2022"``: Precipitable water of the NCEP reanalysis subset in the 2022
   * ``"sst_mnmean_oisst"``: NOAA Optimum Interpolation (OI) SST V2 (from https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html)

   Parameters
   ----------
   name : :py:class:`str <str>`
       Name of the file containing the dataset.
       e.g. 'air_202201_mon_mean'
   cache_dir : path-like, optional
       The directory in which to search for and write cached data.
   cache : dim: :py:class:`bool <bool>`, optional
       If True, then cache data locally for use on subsequent calls
   progressbar: :py:class:`bool <bool>`, default `False`.
       If True, will print a progress bar of the download to standard error (stderr).
   **kws : :py:class:`dict <dict>`, optional
       Passed to xarray.open_dataset

   Returns
   -------
   :py:class:`xarray.Dataset<xarray.Dataset>`

   Reference
   --------------
   - Kalnay et al.,The NCEP/NCAR 40-year reanalysis project, Bull. Amer. Meteor. Soc., 77, 437-470, 1996
   - Rayner, N. A.; Parker, D. E.; Horton, E. B.; Folland, C. K.; Alexander, L. V.; Rowell, D. P.; Kent, E. C.; Kaplan, A. (2003) Global analyses of sea surface temperature, sea ice, and night marine air temperature since the late nineteenth century J. Geophys. Res.Vol. 108, No. D14, 4407 10.1029/2002JD002670  (pdf ~9Mb)

   .. seealso::
       - :py:func:`xarray.tutorial.load_dataset<xarray.tutorial.load_dataset>`
       - :py:func:`xarray.open_dataset<xarray.open_dataset>`
       - :py:func:`xarray.load_dataset<xarray.load_dataset>`


.. py:class:: DataNode(name='root')

   A hierarchical data structure that dynamically manages attributes and nested nodes.

   The :py:class:`DataNode <DataNode>` class provides a flexible way to organize and access data in a tree-like
   structure. It supports automatic creation of nested nodes, path-style access,
   and rich HTML representation in Jupyter environments.

   Parameters
   ----------
   name : :py:class:`str <str>`, optional
       The name of the root node (default: ``"root"``).


   .. py:attribute:: _attributes


   .. py:attribute:: name
      :value: 'root'



   .. py:method:: __getattr__(key)

      Dynamically access or create node attributes.

      Automatically creates nested :py:class:`DataNode <DataNode>` for non-existent attributes,
      while filtering out special IPython/Jupyter attribute requests.

      Parameters
      ----------
      key : :py:class:`str <str>`
          The attribute name to access.

      Returns
      -------
      Any
          The requested attribute or a new :py:class:`DataNode <DataNode>` if attribute doesn't exist.

      Raises
      ------
      AttributeError
          If the attribute is a special IPython/Jupyter attribute.



   .. py:method:: __setattr__(key, value)

      Set node attributes while protecting internal attributes.

      Parameters
      ----------
      key : :py:class:`str <str>`
          The attribute name to set.
      value : Any
          The value to assign to the attribute.



   .. py:method:: __getitem__(key)

      Access attributes using path-style notation (e.g., "path/to/attribute").

      Parameters
      ----------
      key : :py:class:`str <str>`
          The attribute path to access, with '/' separators for nested nodes.

      Returns
      -------
      Any
          The value at the specified path.

      Raises
      ------
      KeyError
          If any part of the path doesn't exist.



   .. py:method:: __contains__(key)


   .. py:method:: __setitem__(key, value)

      Set attributes using path-style notation, creating intermediate nodes as needed.

      Parameters
      ----------
      key : :py:class:`str <str>`
          The attribute path to set, with '/' separators for nested nodes.
      value : Any
          The value to assign at the specified path.



   .. py:method:: _repr_html_()

      Generate HTML representation for Jupyter notebooks.

      Returns
      -------
      :py:class:`str <str>`
          HTML string representing the node and its contents.



   .. py:method:: _format_html()

      Generate complete HTML representation including styles and scripts.

      Returns
      -------
      :py:class:`str <str>`
          Complete HTML document as a string.



   .. py:method:: _format_node_html(node, level=0, parent_id=None)

      Recursively generate HTML for a node and its children.

      Parameters
      ----------
      node : :py:class:`DataNode <DataNode>`
          The node to format.
      level : :py:class:`int <int>`, optional
          Current nesting level (default: 0).
      parent_id : :py:class:`str <str>`, optional
          ID of parent node for DOM construction (default: None).

      Returns
      -------
      :py:class:`str <str>`
          HTML string representing the node.



   .. py:method:: _format_value(value)

      Format values for display, truncating long sequences.

      Parameters
      ----------
      value : Any
          The value to format.

      Returns
      -------
      :py:class:`str <str>`
          Formatted string representation of the value.



   .. py:method:: format_tree(level=0, html=False, is_last_list=None)

      Generate a tree-structured representation of the node.

      Parameters
      ----------
      level : :py:class:`int <int>`, optional
          Current indentation level (default: 0).
      html : :py:class:`bool <bool>`, optional
          Whether to generate HTML output (default: False).
      is_last_list : :py:class:`list <list>` of :py:class:`bool <bool>`, optional
          Track position in hierarchy for proper indentation (default: None).

      Returns
      -------
      :py:class:`str <str>`
          Formatted tree representation.



   .. py:method:: __repr__()


   .. py:method:: _repr_pretty_(p, cycle)

      Support the IPython of pretty printing



   .. py:method:: _repr_mimebundle_(include=None, exclude=None)


   .. py:method:: __dir__()


   .. py:method:: to_zarr(filepath: Union[str, pathlib.Path], zarr_format: Literal[2, 3] = 2)

      Save the DataNode and its contents to a Zarr storage format.

      Parameters
      ----------
      filepath : Union[str, Path]
          Directory path to save the data.
      zarr_format : Literal[2, 3], optional
          Zarr storage format version (default: 2).



   .. py:method:: load(filepath: Union[str, pathlib.Path])
      :classmethod:


      Load a DataNode from a Zarr storage directory.

      Parameters
      ----------
      filepath : Union[str, Path]
          Directory path containing the saved DataNode.

      Returns
      -------
      :py:class:`DataNode <DataNode>`
          The reconstructed DataNode with all its attributes.



.. py:function:: open_datanode(filepath: str) -> DataNode

   Load a DataNode object from native location.

   This function provides a convenient way to load a DataNode that was previously saved
   using the ``DataNode.to_zarr()`` method.

   Parameters
   ----------
   filepath : :py:class:`str <str>`
       The path to the directory containing the saved DataNode data.
       This should be the same path used with DataNode.to_zarr().

   Returns
   -------
   :py:class:`DataNode <DataNode>`
       The loaded DataNode object with all its attributes and nested structure.

   Examples
   --------
   >>> node = open_datanode("path/to/saved_node")
   >>> node.some_attribute  # Access attributes as usual

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_multieof.py


