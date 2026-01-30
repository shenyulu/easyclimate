easyclimate.core
================

.. py:module:: easyclimate.core


Submodules
----------

.. toctree::
   :maxdepth: 1

   /technical/api/easyclimate/core/datanode/index
   /technical/api/easyclimate/core/diff/index
   /technical/api/easyclimate/core/eddy/index
   /technical/api/easyclimate/core/eof/index
   /technical/api/easyclimate/core/extract/index
   /technical/api/easyclimate/core/mk_test/index
   /technical/api/easyclimate/core/normalized/index
   /technical/api/easyclimate/core/read/index
   /technical/api/easyclimate/core/spharm/index
   /technical/api/easyclimate/core/stat/index
   /technical/api/easyclimate/core/stats/index
   /technical/api/easyclimate/core/tutorial/index
   /technical/api/easyclimate/core/units/index
   /technical/api/easyclimate/core/utility/index
   /technical/api/easyclimate/core/variability/index
   /technical/api/easyclimate/core/windspharm/index


Classes
-------

.. autoapisummary::

   easyclimate.core.DataNode


Functions
---------

.. autoapisummary::

   easyclimate.core.calc_gradient
   easyclimate.core.calc_dx_gradient
   easyclimate.core.calc_dlon_radian_gradient
   easyclimate.core.calc_dlon_degree_gradient
   easyclimate.core.calc_dy_gradient
   easyclimate.core.calc_dlat_radian_gradient
   easyclimate.core.calc_dlat_degree_gradient
   easyclimate.core.calc_dx_laplacian
   easyclimate.core.calc_dy_laplacian
   easyclimate.core.calc_dxdy_mixed_derivatives
   easyclimate.core.calc_p_gradient
   easyclimate.core.calc_time_gradient
   easyclimate.core.calc_delta_pressure
   easyclimate.core.calc_p_integral
   easyclimate.core.calc_top2surface_integral
   easyclimate.core.calc_dxdy_laplacian
   easyclimate.core.calc_divergence
   easyclimate.core.calc_vorticity
   easyclimate.core.calc_geostrophic_wind
   easyclimate.core.calc_geostrophic_wind_vorticity
   easyclimate.core.calc_horizontal_water_flux
   easyclimate.core.calc_vertical_water_flux
   easyclimate.core.calc_water_flux_top2surface_integral
   easyclimate.core.calc_divergence_watervaporflux
   easyclimate.core.calc_divergence_watervaporflux_top2surface_integral
   easyclimate.core.calc_u_advection
   easyclimate.core.calc_v_advection
   easyclimate.core.calc_p_advection
   easyclimate.core.calc_shear_stretch_deform
   easyclimate.core.calc_eady_growth_rate
   easyclimate.core.calc_apparent_heat_source
   easyclimate.core.calc_total_diabatic_heating
   easyclimate.core.calc_apparent_moisture_sink
   easyclimate.core.calc_Plumb_wave_activity_horizontal_flux
   easyclimate.core.calc_TN_wave_activity_horizontal_flux
   easyclimate.core.calc_EP_horizontal_flux
   easyclimate.core.calc_monthly_rossby_wave_source
   easyclimate.core.get_specific_years_data
   easyclimate.core.get_specific_months_data
   easyclimate.core.get_specific_days_data
   easyclimate.core.get_specific_hours_data
   easyclimate.core.get_specific_minutes_data
   easyclimate.core.get_specific_seconds_data
   easyclimate.core.get_specific_microseconds_data
   easyclimate.core.get_specific_nanoseconds_data
   easyclimate.core.get_specific_dayofweek_data
   easyclimate.core.get_yearmean_for_specific_months_data
   easyclimate.core.get_year_exceed_index_upper_bound
   easyclimate.core.get_year_exceed_index_lower_bound
   easyclimate.core.get_time_exceed_index_upper_bound
   easyclimate.core.get_time_exceed_index_lower_bound
   easyclimate.core.open_muliti_dataset
   easyclimate.core.calc_linregress_spatial
   easyclimate.core.calc_detrend_spatial
   easyclimate.core.calc_corr_spatial
   easyclimate.core.calc_leadlag_corr_spatial
   easyclimate.core.calc_multiple_linear_regression_spatial
   easyclimate.core.calc_ttestSpatialPattern_spatial
   easyclimate.core.calc_windmask_ttestSpatialPattern_spatial
   easyclimate.core.calc_levenetestSpatialPattern_spatial
   easyclimate.core.calc_skewness_spatial
   easyclimate.core.calc_kurtosis_spatial
   easyclimate.core.calc_theilslopes_spatial
   easyclimate.core.calc_lead_lag_correlation_coefficients
   easyclimate.core.calc_timeseries_correlations
   easyclimate.core.calc_non_centered_corr
   easyclimate.core.calc_pattern_corr
   easyclimate.core.remove_sst_trend
   easyclimate.core.calc_detrend_spatial_fast
   easyclimate.core.calc_monthly_mean
   easyclimate.core.calc_monthly_sum
   easyclimate.core.calc_monthly_std
   easyclimate.core.calc_monthly_var
   easyclimate.core.calc_monthly_max
   easyclimate.core.calc_monthly_min
   easyclimate.core.calc_yearly_mean
   easyclimate.core.calc_yearly_sum
   easyclimate.core.calc_yearly_std
   easyclimate.core.calc_yearly_var
   easyclimate.core.calc_yearly_max
   easyclimate.core.calc_yearly_min
   easyclimate.core.calc_all_climatological_mean
   easyclimate.core.calc_seasonal_climatological_mean
   easyclimate.core.calc_seasonal_cycle_mean
   easyclimate.core.calc_seasonal_cycle_std
   easyclimate.core.calc_seasonal_cycle_var
   easyclimate.core.calc_seasonal_mean
   easyclimate.core.remove_seasonal_cycle_mean
   easyclimate.core.calc_monthly_climatological_std_without_seasonal_cycle_mean
   easyclimate.core.calc_monthly_climatological_var_without_seasonal_cycle_mean
   easyclimate.core.smooth_daily_annual_cycle
   easyclimate.core.calc_daily_annual_cycle_mean
   easyclimate.core.calc_daily_annual_cycle_std
   easyclimate.core.calc_daily_annual_cycle_var
   easyclimate.core.remove_smooth_daily_annual_cycle_mean
   easyclimate.core.calc_horizontal_wind_components_std
   easyclimate.core.calc_windspeed_dataset
   easyclimate.core.calc_windspeed_dataarray
   easyclimate.core.populate_monmean2everymon
   easyclimate.core.populate_daymean2everyday
   easyclimate.core.calc_daily_climatological_anomaly
   easyclimate.core.remove_low_frequency_signal
   easyclimate.core.open_tutorial_dataset
   easyclimate.core.open_datanode
   easyclimate.core.transfer_units_coeff
   easyclimate.core.transfer_data_multiple_units
   easyclimate.core.transfer_data_difference_units
   easyclimate.core.transfer_data_temperature_units
   easyclimate.core.transfer_data_units


Package Contents
----------------

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
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 2.

   Returns
   -------
   The gradient along the coordinate `dim` direction (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`numpy.gradient <numpy:numpy.gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_dx_gradient(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon', lat_dim: str = 'lat', min_dx: float = 1.0, edge_order: int = 2, R: float = 6371200.0) -> xarray.DataArray | xarray.Dataset

   Calculate the zonal gradient of the input data in physical units (meters).

   This function computes the partial derivative :math:`\partial F / \partial x`,
   where :math:`x` is the eastward distance along a parallel (in meters). It is
   the full physical zonal gradient on a sphere, given by:

   .. math::

       \frac{\partial F}{\partial x} = \frac{1}{R \cos \varphi} \cdot \frac{\partial F}{\partial \lambda}

   where :math:`R` is the Earth's radius, :math:`\varphi` is latitude, and
   :math:`\lambda` is longitude in radians. This is essential for dynamical
   calculations like advection or wave propagation in atmospheric/oceanic models.

   The computation uses finite differences along the longitude dimension:
   :math:`\partial F / \partial x = (\partial F / \partial i) / (\partial x / \partial i)`,
   where :math:`i` is the grid index. Longitude is assumed in degrees and
   converted to radians; latitude is broadcasted for the cosine factor.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated. If a Dataset, the gradient is applied to all
       data variables.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default, the gradient is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. Used to compute the cosine factor; broadcasted if necessary.
   min_dx: :py:class:`float <float>`, default: `1.0` (:math:`\mathrm{m}`).
       The minimum acceptable value of `dx` (zonal spacing in meters), below which the output
       is set to NaN to avoid numerical instabilities from very small grid spacings.
       Set to a negative value to disable this check.
   edge_order: {1, 2}, optional
       Order of the finite difference used at the boundaries. 1 uses first-order accurate
       one-sided differences; 2 uses second-order accurate one-sided differences. Default: 2.
   R: :py:class:`float <float>`, default: `6370000` (:math:`\mathrm{m}`).
       Radius of the Earth in meters (approximate mean radius). Can be adjusted for specific models.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The zonal gradient :math:`\partial F / \partial x`, with the same shape and coordinates
       as the input. Units are those of the input data divided by meters (e.g., if F is in K,
       output is K/m). Invalid regions (dx < min_dx) are NaN.

   .. seealso::
       - :py:func:`calc_gradient <calc_gradient>`
       - :py:func:`calc_dlon_radian_gradient <calc_dlon_radian_gradient>`
       - :py:func:`calc_dy_gradient <calc_dy_gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_dlon_radian_gradient(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon', edge_order: int = 2) -> xarray.DataArray | xarray.Dataset

   Calculate the gradient of the input data with respect to longitude in radians.

   This function computes the partial derivative :math:`\partial F / \partial \lambda`,
   where :math:`\lambda` is the longitude in radians. It is useful for spherical coordinate
   calculations, such as in wave activity flux (WAF) formulations, where angular gradients
   must be in radians for consistency with trigonometric functions and the Earth's radius.

   The computation uses finite differences:

   .. math::

       \frac{\partial F}{\partial \lambda} = \frac{\partial F}{\partial i } / \frac{\partial \lambda}{\partial i},

   where :math:`i` is the grid index along the longitude dimension. Longitude values are
   assumed to be in degrees initially and converted to radians for the denominator.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated. If a Dataset, the gradient is applied to all data variables.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default, the gradient is applied over the `lon` dimension.
   edge_order: {1, 2}, optional
       Order of the finite difference used at the boundaries. 1 uses first-order accurate
       one-sided differences; 2 uses second-order accurate one-sided differences. Default: 2.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The gradient :math:`\partial F / \partial \lambda` (in radians), with the same shape
       and coordinates as the input. Units are inherited from the input data divided by radians
       (e.g., if F is in K, output is K/rad).


.. py:function:: calc_dlon_degree_gradient(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon', edge_order: int = 2) -> xarray.DataArray | xarray.Dataset

   Calculate the gradient of the input data with respect to longitude in degrees.

   This function computes the partial derivative :math:`\partial F / \partial \lambda`,
   where :math:`\lambda` is the longitude in degrees. It is suitable for general-purpose
   gradient calculations where degree units are preferred for interpretability.

   The computation uses finite differences:

   .. math::

       \frac{\partial F}{\partial \lambda} = \frac{\partial F}{\partial i} / \frac{\partial \lambda}{\partial i},

   where :math:`i` is the grid index along the longitude dimension. Longitude values remain
   in degrees for the denominator.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated. If a Dataset, the gradient is applied to all
       data variables.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default, the gradient is applied over the `lon` dimension.
   edge_order: {1, 2}, optional
       Order of the finite difference used at the boundaries. 1 uses first-order accurate
       one-sided differences; 2 uses second-order accurate one-sided differences. Default: 2.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The gradient :math:`\partial F / \partial \lambda` (in degrees), with the same shape
       and coordinates as the input. Units are inherited from the input data divided by degrees
       (e.g., if F is in K, output is K/deg).

   .. seealso::
       - :py:func:`calc_gradient <calc_gradient>`
       - :py:func:`calc_dlon_radian_gradient <calc_dlon_radian_gradient>`
       - :py:func:`calc_dlat_degree_gradient <calc_dlat_degree_gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_dy_gradient(data_input: xarray.DataArray | xarray.Dataset, lat_dim: str = 'lat', min_dy: float = 1.0, edge_order: int = 2, R: float = 6371200.0) -> xarray.DataArray | xarray.Dataset

   Calculate the meridional gradient of the input data in physical units (meters).

   This function computes the partial derivative :math:`\partial F / \partial y`,
   where :math:`y` is the northward distance along a meridian (in meters). It is
   the full physical meridional gradient on a sphere, given by:

   .. math::

       \frac{\partial F}{\partial y} = \frac{1}{R} \cdot \frac{\partial F}{\partial \varphi}

   where :math:`R` is the Earth's radius and :math:`\varphi` is latitude in radians.
   This is essential for dynamical calculations like advection or vorticity in models.

   The computation uses finite differences along the latitude dimension:
   :math:`\partial F / \partial y = (\partial F / \partial j) / (\partial y / \partial j)`,
   where :math:`j` is the grid index. Latitude is assumed in degrees and
   converted to radians.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated. If a Dataset, the gradient is applied to all
       data variables.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default, the gradient is applied over the `lat` dimension.
   min_dy: :py:class:`float <float>`, default: `1.0`.
       The minimum acceptable value of `dy` (meridional spacing in meters), below which the output
       is set to NaN to avoid numerical instabilities from very small grid spacings.
       Set to a negative value to disable this check. Unit: meters.
   edge_order: {1, 2}, optional
       Order of the finite difference used at the boundaries. 1 uses first-order accurate
       one-sided differences; 2 uses second-order accurate one-sided differences. Default: 2.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth in meters (approximate mean radius). Can be adjusted for specific models.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The meridional gradient :math:`\partial F / \partial y`, with the same shape and coordinates
       as the input. Units are those of the input data divided by meters (e.g., if F is in K,
       output is K/m). Invalid regions (dy < min_dy) are NaN.

   .. seealso::
       - :py:func:`calc_gradient <calc_gradient>`
       - :py:func:`calc_dlat_radian_gradient <calc_dlat_radian_gradient>`
       - :py:func:`calc_dx_gradient <calc_dx_gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_dlat_radian_gradient(data_input: xarray.DataArray | xarray.Dataset, lat_dim: str = 'lat', edge_order: int = 2) -> xarray.DataArray | xarray.Dataset

   Calculate the gradient of the input data with respect to latitude in radians.

   This function computes the partial derivative :math:`\partial F / \partial \phi`,
   where :math:`\phi` is the latitude in radians. It is useful for spherical coordinate
   calculations, such as in wave activity flux (WAF) formulations.

   The computation uses finite differences:

   .. math::

       \frac{\partial F}{\partial \phi} = \frac{\partial F}{\partial j} / \frac{\partial \phi}{\partial j},

   where :math:`j` is the grid index along the latitude dimension. Latitude values are
   assumed to be in degrees initially and converted to radians for the denominator.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated. If a Dataset, the gradient is applied to all
       data variables.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default, the gradient is applied over the `lat` dimension.
   edge_order: {1, 2}, optional
       Order of the finite difference used at the boundaries. 1 uses first-order accurate
       one-sided differences; 2 uses second-order accurate one-sided differences. Default: 2.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The gradient :math:`\partial F / \partial \phi` (in radians), with the same shape
       and coordinates as the input. Units are inherited from the input data divided by radians
       (e.g., if F is in K, output is K/rad).

   .. seealso::
       - :py:func:`calc_gradient <calc_gradient>`
       - :py:func:`calc_dlat_degree_gradient <calc_dlat_degree_gradient>`
       - :py:func:`calc_dlon_radian_gradient <calc_dlon_radian_gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_dlat_degree_gradient(data_input: xarray.DataArray | xarray.Dataset, lat_dim: str = 'lat', edge_order: int = 2) -> xarray.DataArray | xarray.Dataset

   Calculate the gradient of the input data with respect to latitude in degrees.

   This function computes the partial derivative :math:`\partial F / \partial \phi`,
   where :math:`\phi` is the latitude in degrees. It is suitable for general-purpose
   gradient calculations where degree units are preferred.

   The computation uses finite differences:

   .. math::

       \frac{\partial F}{\partial \phi} = \frac{\partial F}{\partial j} / \frac{\partial \phi}{\partial j},

   where :math:`j` is the grid index along the latitude dimension. Latitude values remain
   in degrees for the denominator.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated. If a Dataset, the gradient is applied to all
       data variables.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default, the gradient is applied over the `lat` dimension.
   edge_order: {1, 2}, optional
       Order of the finite difference used at the boundaries. 1 uses first-order accurate
       one-sided differences; 2 uses second-order accurate one-sided differences. Default: 2.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The gradient :math:`\partial F / \partial \phi` (in degrees), with the same shape
       and coordinates as the input. Units are inherited from the input data divided by degrees
       (e.g., if F is in K, output is K/deg).

   .. seealso::
       - :py:func:`calc_gradient <calc_gradient>`
       - :py:func:`calc_dlat_radian_gradient <calc_dlat_radian_gradient>`
       - :py:func:`calc_dlon_degree_gradient <calc_dlon_degree_gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_dx_laplacian(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon', lat_dim: str = 'lat', min_dx2: float = 1000000000.0, edge_order: int = 2, R: float = 6371200.0) -> xarray.DataArray | xarray.Dataset

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


.. py:function:: calc_dy_laplacian(data_input: xarray.DataArray | xarray.Dataset, lat_dim: str = 'lat', min_dy2: float = 1.0, edge_order: int = 2, R: float = 6371200.0) -> xarray.DataArray | xarray.Dataset

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


.. py:function:: calc_dxdy_mixed_derivatives(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon', lat_dim: str = 'lat', min_dxdy: float = 10000000000.0, edge_order: int = 2, R: float = 6371200.0) -> xarray.DataArray | xarray.Dataset

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


.. py:function:: calc_delta_pressure(data_input: xarray.DataArray, surface_pressure_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], surface_pressure_data_units: Literal['hPa', 'Pa', 'mbar']) -> xarray.DataArray

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

   Examples
   --------
   The results in :py:func:`geocat.comp.meteorology.delta_pressure <geocat-comp:geocat.comp.meteorology.delta_pressure>`:

   >>> from geocat.comp.meteorology import delta_pressure
   >>> dp = delta_pressure(
   ...     pressure_lev= np.array([1000.,925.,850.,700.,600.,500., 400.,300.,250.,200.,150.,100., 70.,50.,30.,20.,10.]),
   ...     surface_pressure = np.array([1013]),
   ... )
   >>> print(dp)
   [[ 50.5  75.  112.5 125.  100.  100.  100.   75.   50.   50.   50.   40.
      25.   20.   15.   10.    5. ]]

   For comparison, the results in :py:func:`easyclimate.calc_delta_pressure <calc_delta_pressure>`:

   >>> temp_sample = xr.DataArray(
   ...     np.array([[292.,285.,283.,277.,270.,260., 250.,235.,225.,215.,207.,207., 213.,220.,225.,228.,230.]]),
   ...     dims = ("lat", "plev"),
   ...     coords = {"plev": np.array([1000.,925.,850.,700.,600.,500., 400.,300.,250.,200.,150.,100., 70.,50.,30.,20.,10.]),
   ...             "lat": np.array([0])}
   ... )
   >>> dp = ecl.calc_delta_pressure(
   ...     data_input = temp_sample,
   ...     surface_pressure_data = xr.DataArray([1013], dims = "lat"),
   ...     vertical_dim = "plev",
   ...     surface_pressure_data_units = "Pa",
   ...     vertical_dim_units = "Pa",
   ... ).transpose("lat", "plev")
   >>> print(dp)
   <xarray.DataArray 'plev' (lat: 1, plev: 17)> Size: 136B
   array([[ 50.5,  75. , 112.5, 125. , 100. , 100. , 100. ,  75. ,  50. ,
           50. ,  50. ,  40. ,  25. ,  20. ,  15. ,  10. ,   5. ]])
   Coordinates:
   * lat      (lat) int64 8B 0
   * plev     (plev) float64 136B 1e+03 925.0 850.0 700.0 ... 50.0 30.0 20.0 10.0


.. py:function:: calc_p_integral(data_input: xarray.DataArray, vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], normalize: bool = False) -> xarray.DataArray

   Calculate the vertical integral along the barometric pressure direction in the p-coordinate system.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The spatio-temporal data to be calculated.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   normalize: :py:class:`bool<bool>`, default: `True`.
       Whether or not the integral results are averaged over the entire layer.

   Returns
   -------
   The vertical integral along the barometric pressure direction in the p-coordinate system (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. attention::
       This method ignores the effect of topography, so it applies to altitudes **above 900hPa** and is **NOT applicable to the Tibetan Plateau region**.
       For a fully accurate vertical integration, please use the :py:func:`calc_top2surface_integral <calc_top2surface_integral>` function to calculate,
       but the speed of the calculation is slightly slowed down.


.. py:function:: calc_top2surface_integral(data_input: xarray.DataArray, surface_pressure_data: xarray.DataArray, vertical_dim: str, surface_pressure_data_units: Literal['hPa', 'Pa', 'mbar'], vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], method: Literal['Boer1982', 'Trenberth1991', 'vibeta-ncl'] = 'vibeta-ncl', normalize: bool = False) -> xarray.DataArray

   Calculate the vertical integral in the p-coordinate system from the ground to the zenith along the barometric pressure direction.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The spatio-temporal data to be calculated.
   surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Surface level pressure.

   .. warning::

       This parameter only accepts local pressure (slp) and must **NOT** be substituted
       with mean sea level pressure (msl). There is a fundamental difference in their physical meaning and numerical
       characteristics—local pressure reflects atmospheric pressure at the actual elevation of local area,
       while mean sea level pressure is a theoretical value adjusted to sea level.

   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   surface_pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   method: {"Boer1982", "Trenberth1991", "vibeta-ncl"}, default: `vibeta-ncl`.
       vertical integration method. Optional values are ``Boer1982``, ``Trenberth1991`` and ``vibeta-ncl``.

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


.. py:function:: calc_dxdy_laplacian(data_input: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', R: float = 6371200.0, spherical_coord: bool = True) -> xarray.DataArray

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


.. py:function:: calc_divergence(u_data: xarray.DataArray, v_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', R: float = 6371200.0, spherical_coord=True, cyclic_boundary: bool = False, method: Literal['easyclimate', 'uv2dv_cfd-ncl'] = 'uv2dv_cfd-ncl') -> xarray.DataArray

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
   R: :py:class:`float <float>`, default: `6371200.0`.
       Radius of the Earth.
   spherical_coord: :py:class:`bool<bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates. The parameter is applicable only when ``method = easyclimate``.
   cyclic_boundary: :py:class:`bool <bool>`, default: `False`.
       If True, assume cyclic (periodic) boundaries in longitude. The parameter is applicable only when ``method = ddvfidf-ncl``.
   method: {"easyclimate", "ddvfidf-ncl"}, default: `ddvfidf-ncl`.
       The method to calculate horizontal divergence term. Optional values are ``easyclimate`` and ``ddvfidf-ncl``.

   Returns
   -------
   The horizontal divergence term. (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::

       - https://www.ncl.ucar.edu/Document/Functions/Built-in/uv2dv_cfd.shtml
       - Howard B. Bluestein. (1992). Synoptic-Dynamic Meteorology in Midlatitudes: Principles of Kinematics and Dynamics, Vol. 1. p113-114

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_vorticity(u_data: xarray.DataArray, v_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', R: float = 6371200.0, spherical_coord: bool = True, cyclic_boundary: bool = False, method: Literal['easyclimate', 'uv2vr_cfd-ncl'] = 'uv2vr_cfd-ncl') -> xarray.DataArray

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
   cyclic_boundary: :py:class:`bool <bool>`, default: `False`.
       If True, assume cyclic (periodic) boundaries in longitude. The parameter is applicable only when ``method = uv2vr_cfd-ncl``.
   method: {"easyclimate", "uv2vr_cfd-ncl"}, default: `uv2vr_cfd-ncl`.
       The method to calculate horizontal divergence term. Optional values are ``easyclimate`` and ``uv2vr_cfd-ncl``.

   Returns
   -------
   The horizontal relative vorticity term. (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::

       - https://www.ncl.ucar.edu/Document/Functions/Built-in/uv2vr_cfd.shtml
       - Howard B. Bluestein. (1992). Synoptic-Dynamic Meteorology in Midlatitudes: Principles of Kinematics and Dynamics, Vol. 1. p113-114

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_geostrophic_wind(z_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', omega: float = 7.292e-05, g: float = 9.8, R: float = 6371200.0) -> xarray.DataArray

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


.. py:function:: calc_geostrophic_wind_vorticity(z_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', spherical_coord: bool = True, omega: float = 7.292e-05, g: float = 9.8, R: float = 6371200.0, cyclic_boundary: bool = False, method: Literal['easyclimate', 'uv2vr_cfd-ncl'] = 'uv2vr_cfd-ncl') -> xarray.DataArray

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
   cyclic_boundary: :py:class:`bool <bool>`, default: `False`.
       If True, assume cyclic (periodic) boundaries in longitude. The parameter is applicable only when ``method = uv2vr_cfd-ncl``.
   method: {"easyclimate", "uv2vr_cfd-ncl"}, default: `uv2vr_cfd-ncl`.
       The method to calculate horizontal divergence term. Optional values are ``easyclimate`` and ``uv2vr_cfd-ncl``.

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


.. py:function:: calc_water_flux_top2surface_integral(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, surface_pressure_data: xarray.DataArray, surface_pressure_data_units: Literal['hPa', 'Pa', 'mbar'], specific_humidity_data_units: Literal['kg/kg', 'g/kg', 'g/g'], vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], method: Literal['Boer1982', 'Trenberth1991', 'vibeta-ncl'] = 'vibeta-ncl', g: float = 9.8) -> xarray.DataArray

   Calculate the water vapor flux across the vertical level.

   .. math::

       \frac{1}{g} \int_0^{p_s} (q\mathbf{v}),dp

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
   specific_humidity_data_units: :py:class:`str <str>`.
       The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   method: :py:class:`str <str>`, default: `'Trenberth1991'`.
       Vertical integration method. Optional values are `Boer1982`, `'Trenberth1991'`.
   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.

   Returns
   -------
   The water vapor flux. (:py:class:`xarray.Dataset<xarray.Dataset>`, :math:`\mathrm{kg \cdot m^-1 \cdot s^-1 }`).

   - :math:`qu`: zonal water vapor flux.
   - :math:`qv`: meridional water vapor flux.

   .. seealso::
       :py:func:`calc_top2surface_integral <calc_top2surface_integral>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_divergence_watervaporflux(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, specific_humidity_data_units: Literal['kg/kg', 'g/kg', 'g/g'], spherical_coord: bool = True, cyclic_boundary: bool = False, method: Literal['easyclimate', 'uv2dv_cfd-ncl'] = 'uv2dv_cfd-ncl', lon_dim: str = 'lon', lat_dim: str = 'lat', g: float = 9.8, R: float = 6371200.0) -> xarray.DataArray

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
       Whether or not to compute the horizontal Laplace term in spherical coordinates. The parameter is applicable only when ``method = easyclimate``.
   cyclic_boundary: :py:class:`bool <bool>`, default: `False`.
       If True, assume cyclic (periodic) boundaries in longitude. The parameter is applicable only when ``method = ddvfidf-ncl``.
   method: {"easyclimate", "ddvfidf-ncl"}, default: `ddvfidf-ncl`.
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


.. py:function:: calc_divergence_watervaporflux_top2surface_integral(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, surface_pressure_data: xarray.DataArray, vertical_dim: str, specific_humidity_data_units: Literal['kg/kg', 'g/kg', 'g/g'], surface_pressure_data_units: Literal['hPa', 'Pa', 'mbar'], vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], spherical_coord: bool = True, cyclic_boundary: bool = False, lon_dim: str = 'lon', lat_dim: str = 'lat', integral_method: Literal['Boer1982', 'Trenberth1991', 'vibeta-ncl'] = 'vibeta-ncl', div_method: Literal['easyclimate', 'uv2dv_cfd-ncl'] = 'uv2dv_cfd-ncl', g: float = 9.8, R: float = 6371200.0) -> xarray.DataArray

   Calculate water vapor flux divergence across the vertical level.

   .. math::

       \nabla \cdot \frac{1}{g} \int_0^{p_s} (q\mathbf{v}),dp

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
       Whether or not to compute the horizontal Laplace term in spherical coordinates. The parameter is applicable only when ``method = easyclimate``.
   cyclic_boundary: :py:class:`bool <bool>`, default: `False`.
       If True, assume cyclic (periodic) boundaries in longitude. The parameter is applicable only when ``method = ddvfidf-ncl``.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   integral_method: {"Boer1982", "Trenberth1991", "vibeta-ncl"}, default: `vibeta-ncl`.
       The vertical integration method. Optional values are ``Boer1982``, ``Trenberth1991`` and ``vibeta-ncl``.

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

   div_method: {"easyclimate", "ddvfidf-ncl"}, default: `ddvfidf-ncl`.
       The method to calculate horizontal divergence term. Optional values are ``easyclimate`` and ``ddvfidf-ncl``.

   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The water vapor flux divergence. (:py:class:`xarray.DataArray<xarray.DataArray>`, :math:`\mathrm{kg \cdot m^-2 \cdot s^-1 }`).

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_u_advection(u_data: xarray.DataArray, temper_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', min_dx: float = 1.0, edge_order: int = 2, R: float = 6371200.0) -> xarray.DataArray

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


.. py:function:: calc_v_advection(v_data: xarray.DataArray, temper_data: xarray.DataArray, lat_dim: str = 'lat', min_dy: float = 1.0, edge_order: int = 2, R: float = 6371200.0) -> xarray.DataArray

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


.. py:function:: calc_p_advection(omega_data: xarray.DataArray, temper_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar']) -> xarray.DataArray

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


.. py:function:: calc_shear_stretch_deform(u_data: xarray.DataArray, v_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', edge_order: {1, 2} = 2, R: float = 6371200.0)

   - https://www.ncl.ucar.edu/Document/Functions/Contributed/shear_stretch_deform_cfd.shtml
   - Spensberger, C., & Spengler, T. (2014). A New Look at Deformation as a Diagnostic for Large-Scale Flow. Journal of the Atmospheric Sciences, 71(11), 4221-4234. https://doi.org/10.1175/JAS-D-14-0108.1


.. py:function:: calc_eady_growth_rate(u_daily_data: xarray.DataArray, z_daily_data: xarray.DataArray, temper_daily_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], lat_dim: str = 'lat') -> xarray.Dataset

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

   Returns
   -------
   The maximum Eady growth rate (:py:class:`xarray.Dataset<xarray.Dataset>`).

   - `eady_growth_rate`: The maximum Eady growth rate.
   - `dudz`: :math:`\frac{\mathrm{d} U}{\mathrm{d} z}`
   - `brunt_vaisala_frequency`: Brunt-väisälä frequency.

   .. seealso::
       - Eady, E. T. (1949). Long Waves and Cyclone Waves. Tellus, 1(3), 33–52. https://doi.org/10.3402/tellusa.v1i3.8507, https://www.tandfonline.com/doi/abs/10.3402/tellusa.v1i3.8507
       - Lindzen, R. S. , & Farrell, B. (1980). A Simple Approximate Result for the Maximum Growth Rate of Baroclinic Instabilities. Journal of Atmospheric Sciences, 37(7), 1648-1654. https://journals.ametsoc.org/view/journals/atsc/37/7/1520-0469_1980_037_1648_asarft_2_0_co_2.xml
       - Simmonds, I., and E.-P. Lim (2009), Biases in the calculation of Southern Hemisphere mean baroclinic eddy growth rate, Geophys. Res. Lett., 36, L01707, https://doi.org/10.1029/2008GL036320.
       - Sloyan, B. M., and T. J. O'Kane (2015), Drivers of decadal variability in the Tasman Sea, J. Geophys. Res. Oceans, 120, 3193–3210, https://doi.org/10.1002/2014JC010550.
       - `eady_growth_rate -NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/eady_growth_rate.shtml>`__
       - `瞬变涡旋诊断量 <https://renqlsysu.github.io/2020/02/16/wave_activity_flux/>`__

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_egr.py


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


.. py:function:: calc_TN_wave_activity_horizontal_flux(z_prime_data: xarray.DataArray | None, u_climatology_data: xarray.DataArray, v_climatology_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], psi_prime_data: xarray.DataArray | None = None, lon_dim: str = 'lon', lat_dim: str = 'lat', time_dim: str = 'time', omega: float = 7.292e-05, g: float = 9.8, R: float = 6371200.0) -> xarray.DataArray

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

   .. attention:: The unit of `z_prime_data` should be **meters**, NOT :math:`\mathrm{m^2 \cdot s^2}` which is the unit used in the representation of potential energy.

   u_climatology_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The climatology of zonal wind data.
   v_climatology_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The climatology of meridional wind data.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   psi_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Perturbation stream function. Geostrophic stream function perturbation calculated from geopotential height anomaly.
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

   Reference
   --------------
   - Takaya, K., & Nakamura, H. (2001). A Formulation of a Phase-Independent Wave-Activity Flux for Stationary and Migratory Quasigeostrophic Eddies on a Zonally Varying Basic Flow. Journal of the Atmospheric Sciences, 58(6), 608-627. https://journals.ametsoc.org/view/journals/atsc/58/6/1520-0469_2001_058_0608_afoapi_2.0.co_2.xml

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


.. py:function:: calc_monthly_rossby_wave_source(u_data: xarray.DataArray, v_data: xarray.DataArray, u_climatology_data: xarray.DataArray, v_climatology_data: xarray.DataArray, lat_dim: str = 'lat', omega: float = 7.292e-05, R: float = 6371200.0) -> xarray.DataArray

   Calculate the Rossby wave source following Sardeshmukh and Hoskins (1988).

   The Rossby wave source (RWS) represents the forcing term for stationary Rossby waves
   and is widely used to diagnose atmospheric teleconnections. It is decomposed into
   five terms representing different physical mechanisms:

   .. math::

       S' = -\nabla \cdot (\mathbf{v}_\chi \zeta)' = \text{term1a} + \text{term1b} + \text{term2a} + \text{term2b} + \text{term3} + \text{term4} + \text{term5}

   where:

   - **term1a**: :math:`-\bar{\zeta} \nabla \cdot \mathbf{v}'_\chi`
     (stretching of climatological vorticity by divergent wind anomalies)
   - **Term1b**: :math:`-f \nabla \cdot \mathbf{v'_\chi}`
     (Stretching of planetary vorticity by anomalous divergence)
   - **term2a**: :math:`-\mathbf{v}'_\chi \cdot \nabla\bar{\zeta}`
     (advection of climatological vorticity by divergent wind anomalies)
   - **term2b**: :math:`-\beta v'_\chi`
     (planetary vorticity advection, where :math:`\beta = \partial f/\partial y`)
   - **term3**: :math:`-\zeta' \nabla \cdot \bar{\mathbf{v}}_\chi`
     (stretching of vorticity anomalies by climatological divergent wind)
   - **term4**: :math:`-\bar{\mathbf{v}}_\chi \cdot \nabla\zeta'`
     (advection of vorticity anomalies by climatological divergent wind)
   - **term5**: :math:`-\nabla \cdot (\mathbf{v}'_\chi \zeta')`
     (nonlinear term: divergence of transient vorticity flux)

   where :math:`\mathbf{v}_\chi = (u_\chi, v_\chi)` is the divergent (irrotational)
   wind component, :math:`\zeta` is relative vorticity, overbar denotes climatology,
   and prime denotes anomaly.

   .. note::
       - Anomalies are computed as deviations from the monthly climatology.
       - The planetary term (term6) represents the beta effect, which is important for Rossby wave propagation, especially in the tropics and subtropics.
       - Terms 1, 2, and 6 typically dominate the Rossby wave source.
       - Positive RWS indicates a cyclonic wave source; negative indicates anticyclonic.

   Parameters
   ----------
   u_data : :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{m/s}`)
       Zonal wind component. Must contain a 'time' dimension with monthly or
       higher frequency data. Expected dimensions: ``(time, lat, lon)`` or ``(time, level, lat, lon)``.
   v_data : :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{m/s}`)
       Meridional wind component. Must have the same dimensions as `u_data`.
   u_climatology_data : :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{m/s}`)
       Monthly climatology of zonal wind. Expected dimensions: ``(time, lat, lon)``
       or ``(time, level, lat, lon)``, where time/month ranges from 1 to 12.
   v_climatology_data : :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{m/s}`)
       Monthly climatology of meridional wind. Must have the same dimensions
       as `u_climatology_data`.

   Returns
   -------
   result : :py:class:`xarray.Dataset<xarray.Dataset>`
       Dataset containing the Rossby wave source and its five decomposed terms:

       - **RWS**: Total Rossby wave source (sum of all terms)
       - **term1a**: Vorticity stretching by anomalous divergence
       - **term1b**: Stretching of planetary vorticity by anomalous divergence
       - **term2a**: Vorticity advection by anomalous divergent wind
       - **term2b**: Planetary vorticity advection (beta effect)
       - **term3**: Anomalous vorticity stretching by climatological divergence
       - **term4**: Anomalous vorticity advection by climatological divergent wind
       - **term5**: Nonlinear transient eddy term


       All variables have units of :math:`\mathrm{s^{-2}}` and retain the spatial/temporal dimensions
       of the input data.

   .. seealso::

       - Sardeshmukh, P. D., & Hoskins, B. J. (1988). The generation of global rotational flow by steady idealized tropical divergence. Journal of the Atmospheric Sciences, 45(7), 1228-1251. https://journals.ametsoc.org/view/journals/atsc/45/7/1520-0469_1988_045_1228_tgogrf_2_0_co_2.xml, https://doi.org/10.1175/1520-0469(1988)045<1228:TGOGRF>2.0.CO;2

   Examples
   --------
   >>> # Calculate monthly climatology
   >>> u_clim = u.groupby('time.month').mean('time')
   >>> v_clim = v.groupby('time.month').mean('time')
   >>>
   >>> # Compute Rossby wave source
   >>> rws_result = calc_rossby_wave_source(u, v, u_clim, v_clim)
   >>>
   >>> # Visualize the dominant terms
   >>> rws_result['RWS'].sel(time='2015-12').plot()
   >>> (rws_result['term1'] + rws_result['term2']).sel(time='2015-12').plot()


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


.. py:function:: remove_sst_trend(ssta: xarray.DataArray, spatial_dims: list[str] = ['lat', 'lon']) -> xarray.DataArray

   Remove the global mean SST anomaly trend from the SST anomaly field
   for EOF/PC analysis and so on.

   This function computes the deviation of the SST anomaly at each grid point
   :math:`(x, y)` from the global-mean SST anomaly for the same time step :math:`t`, following
   the approach described in Zhang et al. (1997).
   The resulting field is:

   .. math::

       \mathrm{SSTA}^*_{x,y,t} = \mathrm{SSTA}_{x,y,t} - [\mathrm{SSTA}]_t

   where :math:`[\mathrm{SSTA}]_t` is the spatial mean over the specified dimensions
   for time :math:`t`. This removes the dominant **global warming mode** trend, avoiding
   unrealistic orthogonality constraints in subsequent PC analysis.

   Parameters
   ----------
   ssta : :py:class:`xarray.DataArray <xarray.DataArray>`
       The SST anomaly field with dimensions including 'time' and the spatial dimensions
       (e.g., 'lat', 'lon').
   spatial_dims : :py:class:`list <list>` of :py:class:`str <str>`, optional
       The spatial dimensions over which to compute the global mean. Default: `['lat', 'lon']`.

   Returns
   -------
   :py:class:`xarray.DataArray <xarray.DataArray>`
       The deviation SST anomaly field with the same shape and coordinates as the input.

   Reference
   --------------
   - Zhang, Y., Wallace, J. M., & Battisti, D. S. (1997). ENSO-like Interdecadal Variability: 1900-93. Journal of Climate, 10(5), 1004-1020. https://journals.ametsoc.org/view/journals/clim/10/5/1520-0442_1997_010_1004_eliv_2.0.co_2.xml

   Examples
   --------
   >>> import xarray as xr
   >>> # Load or create an xarray Dataset with SST anomalies
   >>> ds = xr.open_dataset('path_to_sst_data.nc')
   >>> # Assume 'ssta' is the variable name for SST anomalies
   >>> ssta_dev = remove_sst_tendency(ds['ssta'])
   >>> print(ssta_dev)
   <xarray.DataArray 'ssta' (time: 120, lat: 180, lon: 360)>
   Coordinates:
     * time     (time) datetime64[ns] 1940-01-01 ... 2024-12-01
     * lat      (lat) float32 -89.5 -88.5 -87.5 ... 87.5 88.5 89.5
     * lon      (lon) float32 0.5 1.5 2.5 ... 357.5 358.5 359.5


.. py:function:: calc_detrend_spatial_fast(data_input: xarray.DataArray, time_dim: str = 'time', min_valid_fraction: float = 0.5, method: Literal['scipy_reduce', 'scipy', 'numpy', 'rust', 'rust_chunked', 'rust_flexible', 'auto'] = 'auto', **kwargs) -> xarray.DataArray

   Remove linear trend along time dimension from spatio-temporal data.

   Supports multiple computation methods with optional automatic selection.

   Parameters
   ----------
   data_input : xr.DataArray
       The spatio-temporal data to be detrended.
   time_dim : str, default "time"
       Name of the time dimension.
   min_valid_fraction : float, default 0.5
       Minimum fraction of valid (finite) values required for detrending.
       Grid points with fewer valid values will be set to NaN.
   method : str, default 'auto'
       Computation method to use:
       - 'scipy_reduce': Simplified version using scipy.signal.detrend
       - 'scipy': Optimized version using scipy.signal.detrend
       - 'numpy': Manual numpy vectorized implementation
       - 'rust': High-performance Rust backend
       - 'rust_chunked': Rust backend with chunked processing (for large datasets)
       - 'rust_flexible': Rust backend with flexible dimension handling
       - 'auto': Automatically selects the best available method
   **kwargs : dict
       Additional arguments passed to specific methods:
       - chunk_size: int (for 'rust_chunked' method)
       - use_chunked: bool (for 'rust_chunked' method)

   Returns
   -------
   xr.DataArray
       Detrended data with same shape and coordinates as input.

   Raises
   ------
   TypeError
       If data_input is not an xarray.DataArray.
   ValueError
       If input parameters are invalid.
   ImportError
       If Rust method is selected but Rust backend is not available.

   Examples
   --------
   >>> import xarray as xr
   >>> import numpy as np
   >>>
   >>> # Create sample data
   >>> data = xr.DataArray(
   ...     np.random.randn(100, 50, 100),
   ...     dims=['time', 'lat', 'lon']
   ... )
   >>>
   >>> # Using scipy method
   >>> result1 = calc_detrend_spatial_fast(data, method='scipy')
   >>>
   >>> # Using numpy method
   >>> result2 = calc_detrend_spatial_fast(data, method='numpy')
   >>>
   >>> # Using Rust method (if available)
   >>> try:
   >>>     result3 = calc_detrend_spatial_fast(data, method='rust')
   >>> except ImportError:
   >>>     print("Rust backend not available")
   >>>
   >>> # Automatic method selection
   >>> result4 = calc_detrend_spatial_fast(data, method='auto')

   Notes
   -----
   - 'scipy_reduce': Simplest method but less robust with NaN values
   - 'scipy': Optimized scipy version with better special value handling
   - 'numpy': Manual vectorized implementation, typically 2-3x faster than scipy
   - 'rust': High-performance Rust implementation, typically 10-50x faster than numpy
   - 'rust_chunked': Rust chunked version for large memory datasets
   - 'rust_flexible': Rust version with flexible dimension ordering
   - 'auto': Uses 'rust' if available, otherwise 'numpy'


.. py:function:: calc_monthly_mean(data_input: xarray.DataArray, dim: str = 'time', **kwargs)

   Calculate monthly mean.

   For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same month it is:

   .. math::
       o(t, x) = \mathrm{mean} \left \lbrace i(t', x), t_1 < t' \leqslant t_n \right\rbrace

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

       .. note::

           The recommended frequence of the `data_input` is daily.

   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.

   Examples
   --------
   >>> import xarray as xr
   >>> import numpy as np
   >>> import pandas as pd
   >>> import easyclimate as ecl
   >>> # Create sample data with daily frequency
   >>> time_index = pd.date_range('2020-01-01', '2020-03-31', freq='D')
   >>> rng = np.random.default_rng(42)
   >>> data = rng.random((len(time_index), 3, 3))
   >>> da = xr.DataArray(data, dims=['time', 'x', 'y'], coords={'time': time_index})
   >>> # Calculate monthly mean
   >>> monthly_mean = ecl.calc_monthly_mean(da)
   >>> print(monthly_mean)
   <xarray.DataArray (time: 3, x: 3, y: 3)> Size: 216B
   array([[[0.50635175, 0.45767908, 0.50271707],
           [0.52091523, 0.44830133, 0.44293946],
           [0.47944591, 0.53314083, 0.48073062]],
       [[0.53431127, 0.48259521, 0.47464862],
           [0.41070456, 0.51619935, 0.4872374 ],
           [0.6009132 , 0.43963445, 0.6028882 ]],
       [[0.48363606, 0.60589154, 0.42622008],
           [0.47519641, 0.46989711, 0.45327877],
           [0.44193025, 0.45050389, 0.641573  ]]])
   Coordinates:
   * time     (time) datetime64[ns] 24B 2020-01-01 2020-02-01 2020-03-01
   Dimensions without coordinates: x, y

   .. seealso::
       :py:func:`numpy.mean <numpy:numpy.mean>`, :py:func:`dask.array.mean <dask:dask.array.mean>`,
       :py:meth:`xarray.DataArray.mean <xarray:xarray.DataArray.mean>`,
       :py:meth:`xarray.core.groupby.DataArrayGroupBy.mean <xarray:xarray.core.groupby.DataArrayGroupBy.mean>`.


.. py:function:: calc_monthly_sum(data_input: xarray.DataArray, dim: str = 'time', **kwargs)

   Calculate monthly sum.

   For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same month it is:

   .. math::
       o(t, x) = \mathrm{sum} \left \lbrace i(t', x), t_1 < t' \leqslant t_n \right\rbrace

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

       .. note::

           The recommended frequence of the `data_input` is daily.

   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.

   Examples
   --------
   >>> import xarray as xr
   >>> import numpy as np
   >>> import pandas as pd
   >>> import easyclimate as ecl
   >>> # Create sample data with daily frequency
   >>> time_index = pd.date_range('2020-01-01', '2020-03-31', freq='D')
   >>> rng = np.random.default_rng(42)
   >>> data = rng.random((len(time_index), 3, 3))
   >>> da = xr.DataArray(data, dims=['time', 'x', 'y'], coords={'time': time_index})
   >>> # Calculate monthly sum
   >>> monthly_sum = ecl.calc_monthly_sum(da)
   >>> print(monthly_sum)
   <xarray.DataArray (time: 3, x: 3, y: 3)> Size: 216B
   array([[[15.69690421, 14.18805154, 15.58422905],
           [16.14837203, 13.89734131, 13.7311233 ],
           [14.86282316, 16.52736588, 14.90264935]],
       [[15.49502686, 13.99526102, 13.76481005],
           [11.91043229, 14.96978113, 14.1298847 ],
           [17.42648281, 12.7493991 , 17.48375789]],
       [[14.99271788, 18.78263759, 13.21282259],
           [14.73108859, 14.56681052, 14.05164186],
           [13.69983774, 13.96562059, 19.88876291]]])
   Coordinates:
   * time     (time) datetime64[ns] 24B 2020-01-01 2020-02-01 2020-03-01
   Dimensions without coordinates: x, y

   .. seealso::
       :py:func:`numpy.sum <numpy:numpy.sum>`, :py:func:`dask.array.sum <dask:dask.array.sum>`,
       :py:meth:`xarray.DataArray.sum <xarray:xarray.DataArray.sum>`,
       :py:meth:`xarray.core.groupby.DataArrayGroupBy.sum <xarray:xarray.core.groupby.DataArrayGroupBy.sum>`.


.. py:function:: calc_monthly_std(data_input: xarray.DataArray, dim: str = 'time', **kwargs)

   Calculate monthly standard deviation.

   For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same month it is:

   .. math::
       o(t, x) = \mathrm{std} \left \lbrace i(t', x), t_1 < t' \leqslant t_n \right\rbrace

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

       .. note::

           The recommended frequence of the `data_input` is daily.

   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.

   Examples
   --------
   >>> import xarray as xr
   >>> import numpy as np
   >>> import pandas as pd
   >>> import easyclimate as ecl
   >>> # Create sample data with daily frequency
   >>> time_index = pd.date_range('2020-01-01', '2020-03-31', freq='D')
   >>> rng = np.random.default_rng(42)
   >>> data = rng.random((len(time_index), 3, 3))
   >>> da = xr.DataArray(data, dims=['time', 'x', 'y'], coords={'time': time_index})
   >>> # Calculate monthly std
   >>> monthly_std = ecl.calc_monthly_std(da)
   >>> print(monthly_std)
   <xarray.DataArray (time: 3, x: 3, y: 3)> Size: 216B
   array([[[0.30528844, 0.29905447, 0.26472868],
           [0.23456056, 0.30879525, 0.29333846],
           [0.26139562, 0.306974  , 0.27987361]],
       [[0.30196879, 0.24783961, 0.26078164],
           [0.2708643 , 0.3012602 , 0.29801453],
           [0.24816804, 0.33863555, 0.25623523]],
       [[0.29399346, 0.31595077, 0.30336434],
           [0.31117807, 0.3130123 , 0.28909393],
           [0.27104435, 0.26864038, 0.22912052]]])
   Coordinates:
   * time     (time) datetime64[ns] 24B 2020-01-01 2020-02-01 2020-03-01
   Dimensions without coordinates: x, y

   .. seealso::
       :py:func:`numpy.std <numpy:numpy.std>`, :py:func:`dask.array.std <dask:dask.array.std>`,
       :py:meth:`xarray.DataArray.std <xarray:xarray.DataArray.std>`,
       :py:meth:`xarray.core.groupby.DataArrayGroupBy.std <xarray:xarray.core.groupby.DataArrayGroupBy.std>`.


.. py:function:: calc_monthly_var(data_input: xarray.DataArray, dim: str = 'time', **kwargs)

   Calculate monthly variance.

   For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same month it is:

   .. math::
       o(t, x) = \mathrm{var} \left \lbrace i(t', x), t_1 < t' \leqslant t_n \right\rbrace

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

       .. note::

           The recommended frequence of the `data_input` is daily.

   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.

   Examples
   --------
   >>> import xarray as xr
   >>> import numpy as np
   >>> import pandas as pd
   >>> import easyclimate as ecl
   >>> # Create sample data with daily frequency
   >>> time_index = pd.date_range('2020-01-01', '2020-03-31', freq='D')
   >>> rng = np.random.default_rng(42)
   >>> data = rng.random((len(time_index), 3, 3))
   >>> da = xr.DataArray(data, dims=['time', 'x', 'y'], coords={'time': time_index})
   >>> # Calculate monthly var
   >>> monthly_var = ecl.calc_monthly_var(da)
   >>> print(monthly_var)
   <xarray.DataArray (time: 3, x: 3, y: 3)> Size: 216B
   array([[[0.09320103, 0.08943358, 0.07008127],
           [0.05501865, 0.09535451, 0.08604745],
           [0.06832767, 0.09423304, 0.07832924]],
       [[0.09118515, 0.06142447, 0.06800706],
           [0.07336747, 0.09075771, 0.08881266],
           [0.06158738, 0.11467404, 0.0656565 ]],
       [[0.08643216, 0.09982489, 0.09202992],
           [0.09683179, 0.0979767 , 0.0835753 ],
           [0.07346504, 0.07216766, 0.05249621]]])
   Coordinates:
   * time     (time) datetime64[ns] 24B 2020-01-01 2020-02-01 2020-03-01
   Dimensions without coordinates: x, y

   .. seealso::
       :py:func:`numpy.var <numpy:numpy.var>`, :py:func:`dask.array.var <dask:dask.array.var>`,
       :py:meth:`xarray.DataArray.var <xarray:xarray.DataArray.var>`,
       :py:meth:`xarray.core.groupby.DataArrayGroupBy.var <xarray:xarray.core.groupby.DataArrayGroupBy.var>`.


.. py:function:: calc_monthly_max(data_input: xarray.DataArray, dim: str = 'time', **kwargs)

   Calculate monthly maximum.

   For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same month it is:

   .. math::
       o(t, x) = \mathrm{max} \left \lbrace i(t', x), t_1 < t' \leqslant t_n \right\rbrace

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

       .. note::

           The recommended frequence of the `data_input` is daily.

   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.

   Examples
   --------
   >>> import xarray as xr
   >>> import numpy as np
   >>> import pandas as pd
   >>> import easyclimate as ecl
   >>> # Create sample data with daily frequency
   >>> time_index = pd.date_range('2020-01-01', '2020-03-31', freq='D')
   >>> rng = np.random.default_rng(42)
   >>> data = rng.random((len(time_index), 3, 3))
   >>> da = xr.DataArray(data, dims=['time', 'x', 'y'], coords={'time': time_index})
   >>> # Calculate monthly max
   >>> monthly_max = ecl.calc_monthly_max(da)
   >>> print(monthly_max)
   <xarray.DataArray (time: 3, x: 3, y: 3)> Size: 216B
   <xarray.DataArray (time: 3, x: 3, y: 3)> Size: 216B
   array([[[0.96189766, 0.95855921, 0.93604357],
           [0.97182643, 0.97069802, 0.97562235],
           [0.99237556, 0.96623191, 0.91601185]],
       [[0.95119466, 0.88414571, 0.85053368],
           [0.94602445, 0.99910473, 0.99546447],
           [0.98002718, 0.98663154, 0.99703466]],
       [[0.989133  , 0.99874337, 0.99308458],
           [0.99032166, 0.98180595, 0.92746046],
           [0.99758004, 0.91879368, 0.99470175]]])
   Coordinates:
   * time     (time) datetime64[ns] 24B 2020-01-01 2020-02-01 2020-03-01
   Dimensions without coordinates: x, y

   .. seealso::
       :py:func:`numpy.max <numpy:numpy.max>`, :py:func:`dask.array.max <dask:dask.array.max>`,
       :py:meth:`xarray.DataArray.max <xarray:xarray.DataArray.max>`,
       :py:meth:`xarray.core.groupby.DataArrayGroupBy.max <xarray:xarray.core.groupby.DataArrayGroupBy.max>`.


.. py:function:: calc_monthly_min(data_input: xarray.DataArray, dim: str = 'time', **kwargs)

   Calculate monthly minimum.

   For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same month it is:

   .. math::
       o(t, x) = \mathrm{min} \left \lbrace i(t', x), t_1 < t' \leqslant t_n \right\rbrace

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

       .. note::

           The recommended frequence of the `data_input` is daily.

   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.

   Examples
   --------
   >>> import xarray as xr
   >>> import numpy as np
   >>> import pandas as pd
   >>> import easyclimate as ecl
   >>> # Create sample data with daily frequency
   >>> time_index = pd.date_range('2020-01-01', '2020-03-31', freq='D')
   >>> rng = np.random.default_rng(42)
   >>> data = rng.random((len(time_index), 3, 3))
   >>> da = xr.DataArray(data, dims=['time', 'x', 'y'], coords={'time': time_index})
   >>> # Calculate monthly min
   >>> monthly_min = ecl.calc_monthly_min(da)
   >>> print(monthly_min)
   <xarray.DataArray (time: 3, x: 3, y: 3)> Size: 216B
   array([[[0.02280387, 0.02485949, 0.05338193],
           [0.02271207, 0.01783678, 0.02161208],
           [0.00736227, 0.04161417, 0.02114802]],
       [[0.0289995 , 0.04347506, 0.0401513 ],
           [0.01072764, 0.09172101, 0.01468284],
           [0.07205915, 0.01230269, 0.00542983]],
       [[0.0040076 , 0.0165798 , 0.00166071],
           [0.04896371, 0.01903415, 0.00123306],
           [0.04737402, 0.00450012, 0.22825288]]])
   Coordinates:
   * time     (time) datetime64[ns] 24B 2020-01-01 2020-02-01 2020-03-01
   Dimensions without coordinates: x, y

   .. seealso::
       :py:func:`numpy.min <numpy:numpy.min>`, :py:func:`dask.array.min <dask:dask.array.min>`,
       :py:meth:`xarray.DataArray.min <xarray:xarray.DataArray.min>`,
       :py:meth:`xarray.core.groupby.DataArrayGroupBy.min <xarray:xarray.core.groupby.DataArrayGroupBy.min>`.


.. py:function:: calc_yearly_mean(data_input: xarray.DataArray, dim: str = 'time', **kwargs)

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


.. py:function:: calc_yearly_sum(data_input: xarray.DataArray, dim: str = 'time', **kwargs)

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


.. py:function:: calc_yearly_std(data_input: xarray.DataArray, dim: str = 'time', **kwargs)

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


.. py:function:: calc_yearly_var(data_input: xarray.DataArray, dim: str = 'time', **kwargs)

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


.. py:function:: calc_yearly_max(data_input: xarray.DataArray, dim: str = 'time', **kwargs)

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


.. py:function:: calc_yearly_min(data_input: xarray.DataArray, dim: str = 'time', **kwargs)

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


.. py:function:: transfer_units_coeff(input_units: str, output_units: str) -> float

   Compute the unit conversion factor from input units to output units.

   Parameters
   ----------
   input_units : :py:class:`str <str>`
       The input unit string (e.g., 'm').
   output_units : :py:class:`str <str>`
       The output unit string (e.g., 'km').

   Returns
   -------
   :py:class:`float <float>`
       The conversion factor to multiply the input values by to obtain output values.

   Example
   -------
   >>> import easyclimate as ecl
   >>> result1 = ecl.transfer_units_coeff("m/s", "km/h")
   >>> result2 = ecl.transfer_units_coeff("hPa", "mbar")
   >>> result3 = ecl.transfer_units_coeff("mm/day", "m/month")
   >>> print(result1, result2, result3)
   3.6 1.0 0.0304375


.. py:function:: transfer_data_multiple_units(input_data: xarray.DataArray | xarray.Dataset, input_units: str, output_units: str) -> xarray.DataArray | xarray.Dataset

   Convert data units for multiplicative transitions (e.g., m to km).

   Parameters
   ----------
   input_data : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`
       The input data with units to convert.
   input_units : :py:class:`str <str>`
       The input unit string. Below is a table of supported temperature unit standards and aliases:
   output_units : :py:class:`str <str>`
       The output unit string.

   Returns
   -------
   :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`
       The converted data with updated 'units' attribute.

   Example
   -------
   >>> import xarray as xr
   >>> import numpy as np
   >>> import easyclimate as ecl
   >>> ds = xr.DataArray(
   ...     np.array([[3e-5, 5.4e-5], [5.2, -75.5]]),
   ...     dims=("lon", "lat"),
   ...     coords={"lon": np.array([160, 70]), "lat": np.array([87.5, -87.5])}
   ... )
   >>> result = ecl.transfer_data_multiple_units(ds, "mm/day", "m/day")
   >>> print(result)
   <xarray.DataArray (lon: 2, lat: 2)> Size: 32B
   array([[ 3.00e-08,  5.40e-08],
       [ 5.20e-03, -7.55e-02]])
   Coordinates:
   * lon      (lon) int64 16B 160 70
   * lat      (lat) float64 16B 87.5 -87.5
   Attributes:
       units:    m/day


.. py:function:: transfer_data_difference_units(input_data: xarray.DataArray | xarray.Dataset, input_units: str, output_units: str) -> xarray.DataArray | xarray.Dataset

   Convert data units for difference-based transitions (e.g., adjustments requiring offset).

   Parameters
   ----------
   input_data : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`
       The input data with units to convert.
   input_units : :py:class:`str <str>`
       The input unit string.
   output_units : :py:class:`str <str>`
       The output unit string.

   Returns
   -------
   :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`
       The converted data with updated 'units' attribute.

   Example
   -------
   >>> import xarray as xr
   >>> import easyclimate as ecl
   >>> result = ecl.transfer_data_difference_units(xr.DataArray(15), "celsius", "kelvin")
   >>> print(result)
   <xarray.DataArray ()> Size: 8B
   array(288.15)
   Attributes:
       units:    kelvin


.. py:function:: transfer_data_temperature_units(input_data: xarray.DataArray | xarray.Dataset, input_units: str, output_units: str) -> xarray.DataArray | xarray.Dataset

   Convert temperature data from one unit to another, supporting aliases.

   Parameters
   ----------
   input_data : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`
       The input temperature data to convert.
   input_units : :py:class:`str <str>`
       The input temperature unit (e.g., 'degC', 'K', '摄氏度'). The available values are as follows:

   .. tab-set::

       .. tab-item:: Kelvin scale

           K, kelvin, degK, 开氏度, ケルビン

       .. tab-item:: Celsius scale

           degC, celsius, degC, 摄氏度, セルシウス度

       .. tab-item:: Fahrenheit scale

           degF, fahrenheit, 华氏度, カ氏度, 華氏度, ファーレンハイト度

       .. tab-item:: Rankine scale

           degR, rankine, 兰氏度, ランキン度

       .. tab-item:: Réaumur scale

           degRe, reaumur, 列氏度, レオミュール度

   output_units : :py:class:`str <str>`
       The output temperature unit (e.g., 'degF', 'kelvin'). The available values are as follows:

   .. tab-set::

       .. tab-item:: Kelvin scale

           K, kelvin, degK, 开氏度, ケルビン

       .. tab-item:: Celsius scale

           degC, celsius, degC, 摄氏度, セルシウス度

       .. tab-item:: Fahrenheit scale

           degF, fahrenheit, 华氏度, カ氏度, 華氏度, ファーレンハイト度

       .. tab-item:: Rankine scale

           degR, rankine, 兰氏度, ランキン度

       .. tab-item:: Réaumur scale

           degRe, reaumur, 列氏度, レオミュール度

   Returns
   -------
   :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`
       The converted temperature data with updated 'units' attribute.

   Example
   -------
   >>> import xarray as xr
   >>> import easyclimate as ecl
   >>> result = ecl.transfer_data_temperature_units(
   ...     xr.DataArray([104, 100, 92, 92, 86, 80, 80, 60, 30]), "摄氏度", "カ氏度"
   ... )
   >>> print(result)
   <xarray.DataArray (dim_0: 9)> Size: 72B
   array([219.2, 212. , 197.6, 197.6, 186.8, 176. , 176. , 140. ,  86. ])
   Dimensions without coordinates: dim_0
   Attributes:
       units:    カ氏度

   .. seealso::
       - https://en.wikipedia.org/wiki/Kelvin
       - https://en.wikipedia.org/wiki/Celsius
       - https://en.wikipedia.org/wiki/Fahrenheit
       - https://en.wikipedia.org/wiki/Rankine_scale
       - https://en.wikipedia.org/wiki/R%C3%A9aumur_scale
       - :py:func:`transfer_data_units <transfer_data_units>`


.. py:function:: transfer_data_units(input_data: xarray.DataArray | xarray.Dataset, input_units: str, output_units: str) -> xarray.DataArray | xarray.Dataset

   Convert data units for any type of transition using Pint.

   .. warning::

       Does NOT support ``dask``.

   Parameters
   ----------
   input_data : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`
       The input data with units to convert.
   input_units : :py:class:`str <str>`
       The input unit string.
   output_units : :py:class:`str <str>`
       The output unit string.

   Returns
   -------
   :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`
       The converted data with updated 'units' attribute.

   Example
   -------
   >>> import xarray as xr
   >>> import easyclimate as ecl
   >>> result = ecl.transfer_data_units(
   ...     xr.DataArray([104, 100, 92, 92, 86, 80, 80, 60, 30]), "degC", "degF"
   ... )
   >>> print(result)
   <xarray.DataArray (dim_0: 9)> Size: 72B
   array([219.2, 212. , 197.6, 197.6, 186.8, 176. , 176. , 140. ,  86. ])
   Dimensions without coordinates: dim_0
   Attributes:
       units:    degF

   .. seealso::
       :py:func:`transfer_data_temperature_units <transfer_data_temperature_units>`


