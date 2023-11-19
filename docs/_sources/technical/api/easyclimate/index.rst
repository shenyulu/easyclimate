:py:mod:`easyclimate`
=====================

.. py:module:: easyclimate


Subpackages
-----------
.. toctree::
   :titlesonly:
   :maxdepth: 3

   core/index.rst
   filter/index.rst
   index/index.rst
   interp/index.rst
   ocean/index.rst
   plot/index.rst
   windspharm/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   easyclimate.BarnesFilter



Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.open_tutorial_dataset
   easyclimate.calc_gradient
   easyclimate.calc_p_gradient
   easyclimate.transfer_deg2rad
   easyclimate.transfer_units_coeff
   easyclimate.calc_brunt_vaisala_frequency_atm
   easyclimate.get_coriolis_parameter
   easyclimate.get_potential_temperature
   easyclimate.calc_static_stability
   easyclimate.find_dims_axis
   easyclimate.transfer_inf2nan
   easyclimate.transfer_data_units
   easyclimate.get_weighted_spatial_data
   easyclimate.generate_dataset_dispatcher
   easyclimate.calc_lon_gradient
   easyclimate.calc_lat_gradient
   easyclimate.calc_lon_laplacian
   easyclimate.calc_lat_laplacian
   easyclimate.calc_lon_lat_mixed_derivatives
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
   easyclimate.transfer_dFdp2dFdz
   easyclimate.calc_eady_growth_rate
   easyclimate.calc_apparent_heat_source
   easyclimate.calc_total_diabatic_heating
   easyclimate.calc_apparent_moisture_sink
   easyclimate.calc_Plumb_wave_activity_horizontal_flux
   easyclimate.calc_TN_wave_activity_horizontal_flux
   easyclimate.calc_TN_wave_activity_3D_flux
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
   easyclimate.calc_yearmean
   easyclimate.calc_yearsum
   easyclimate.calc_yearstd
   easyclimate.calc_yearvar
   easyclimate.calc_yearmax
   easyclimate.calc_yearmin
   easyclimate.assert_compared_version
   easyclimate.transfer_int2datetime
   easyclimate.transfer_datetime2int
   easyclimate.transfer_monmean2everymonthmean
   easyclimate.get_compress_xarraydata
   easyclimate.sort_ascending_latlon_coordinates
   easyclimate.generate_datatree_dispatcher
   easyclimate.transfer_xarray_lon_from180TO360
   easyclimate.transfer_xarray_lon_from360TO180
   easyclimate.module_available
   easyclimate.open_muliti_dataset
   easyclimate.calc_linregress_spatial
   easyclimate.calc_detrend_data
   easyclimate.calc_ttestSpatialPattern_spatial
   easyclimate.calc_skewness_spatial
   easyclimate.calc_kurtosis_spatial
   easyclimate.calc_climatological_mean
   easyclimate.calc_climatological_seasonal_mean
   easyclimate.calc_seasonal_cycle_mean
   easyclimate.calc_seasonal_cycle_std
   easyclimate.calc_seasonal_cycle_var
   easyclimate.remove_seasonal_cycle_mean
   easyclimate.calc_climate_monthly_std
   easyclimate.calc_climate_monthly_var
   easyclimate.calc_horizontal_wind_components_std
   easyclimate.get_EOF_model
   easyclimate.calc_EOF_analysis
   easyclimate.get_EOF_projection
   easyclimate.save_EOF_model
   easyclimate.load_EOF_model
   easyclimate.get_MCA_model
   easyclimate.calc_MCA_analysis
   easyclimate.get_MCA_projection
   easyclimate.save_MCA_model
   easyclimate.load_MCA_model
   easyclimate.calc_index_NPWI
   easyclimate.find_PW_monsoon_region
   easyclimate.cal_NPWI_monsoon_onset
   easyclimate.cal_NPWI_monsoon_detreat
   easyclimate.get_weighted_spatial_data
   easyclimate.sort_ascending_latlon_coordinates
   easyclimate.calc_intensity_STFZ
   easyclimate.calc_intensity_SAFZ
   easyclimate.calc_location_STFZ
   easyclimate.calc_location_SAFZ
   easyclimate.calc_location_line_STFZ
   easyclimate.calc_location_line_SAFZ
   easyclimate.transfer_xarray_lon_from180TO360
   easyclimate.generate_dataset_dispatcher
   easyclimate.interp_mesh2mesh
   easyclimate.interp_point2mesh
   easyclimate.interp_point2mesh_S2
   easyclimate.field_grids
   easyclimate.find_dims_axis
   easyclimate.calc_butter_bandpass
   easyclimate.calc_butter_lowpass
   easyclimate.calc_butter_highpass



.. py:function:: open_tutorial_dataset(name: str, cache: bool = True, cache_dir: None | str | os.PathLike = None, *, engine: xarray.backends.api.T_Engine = None, **kws) -> xarray.Dataset

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
   * ``"mini_HadISST_ice"``: Hadley Centre Sea Ice and Sea Surface Temperature data set (HadISST) subset
   * ``"PressQFF_202007271200_872"``: Observational data from European stations (from https://github.com/EXCITED-CO2/xarray-regrid)


   Parameters
   ----------
   name : str
       Name of the file containing the dataset.
       e.g. 'air_202201_mon_mean'
   cache_dir : path-like, optional
       The directory in which to search for and write cached data.
   cache : bool, optional
       If True, then cache data locally for use on subsequent calls
   **kws : dict, optional
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


.. py:function:: calc_gradient(data_input: xr.DataArray | xr.Dataset, dim: str, varargs=1, edge_order=2) -> xr.DataArray | xr.Dataset

   Compute the gradient along the coordinate `dim` direction.

   The gradient is computed using **second order accurate central differences** in the interior points 
   and either first or second order accurate one-sides (forward or backwards) differences at the boundaries. 
   The returned gradient hence has the same shape as the input array.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
   dim : str
       Dimension(s) over which to apply gradient. By default gradient is applied over the `time` dimension.
   varargs: list of scalar or array, optional
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


.. py:function:: calc_p_gradient(data_input: xarray.DataArray, vertical_dim: str, vertical_dim_units: str) -> xarray.DataArray

   Calculate the gradient along the barometric pressure direction in the p-coordinate system.

   .. math::
       \frac{\partial F}{\partial p}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.

   Returns
   -------
   The gradient along the barometric pressure direction in the p-coordinate system (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>`


.. py:function:: transfer_deg2rad(ds: xarray.DataArray) -> xarray.DataArray

   Convert Degrees to Radians.

   Parameters
   ----------
   - ds: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Degrees data.

   Returns
   -------
   - Radians data.: :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: transfer_units_coeff(input_units: str, output_units: str) -> float

   Unit conversion factor


.. py:function:: calc_brunt_vaisala_frequency_atm(potential_temperature_data: xarray.DataArray, z_data: xarray.DataArray, vertical_dim: str, g=9.8) -> xarray.DataArray

   Calculation of the Brunt-väisälä frequency for the vertical atmosphere.

   .. math::
       N = \left( \frac{g}{\theta} \frac{\mathrm{d}\theta}{\mathrm{d}z} \right)^\frac{1}{2}

   Parameters
   ----------
   potential_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Vertical atmospheric potential temperature.
   z_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Vertical atmospheric geopotential height.

   .. attention:: The unit of `z_data` should be **meters**, NOT :math:`\mathrm{m^2 \cdot s^2}` which is the unit used in the representation of potential energy.

   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   g: :py:class:`float<python.float>`, default: `9.8`.
       The acceleration of gravity.

   Returns
   -------
   Brunt-väisälä frequency (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - `Brunt-väisälä frequency - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Brunt-v%C3%A4is%C3%A4l%C3%A4_frequency>`__

   .. seealso::
       - `brunt_vaisala_frequency — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.brunt_vaisala_frequency.html>`__
       - `brunt_vaisala_atm - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/brunt_vaisala_atm.shtml>`__
       


.. py:function:: get_coriolis_parameter(lat_data, omega=7.292e-05) -> xarray.DataArray

   Calculate the Coriolis parameter at each point.

   .. math::
       f = 2 \Omega \sin(\phi)

   Parameters
   ----------
   lat_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Latitude at each point.
   omega: :py:class:`float<python.float>`, default: `7.292e-5`.
       The angular speed of the earth.

   Returns
   -------
   Corresponding Coriolis force at each point (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - `Coriolis parameter - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Coriolis_parameter>`__

   .. seealso::
       - `coriolis_parameter — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.coriolis_parameter.html>`__
       - `coriolis_param - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/coriolis_param.shtml>`__


.. py:function:: get_potential_temperature(temper_data: xarray.DataArray, vertical_dim: str, vertical_units: str, kappa=287 / 1005.7) -> xarray.DataArray

   Calculate the potential temperature.

   Uses the Poisson equation to calculation the potential temperature given pressure and temperature.

   .. math::
       \theta = T \left( \frac{p_0}{p} \right) ^\kappa

   Parameters
   ----------
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   kappa: :py:class:`float<python.float>`, default: `287/1005.7`.
       Poisson constant :math:`\kappa`.

       .. note::
           `Poisson constant - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Poisson_constant>`__

   Returns
   -------
   Potential temperature corresponding to the temperature and pressure (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - `Potential temperature - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Potential_temperature>`__
   - `Potential-temperature.pdf <http://weatherclimatelab.mit.edu/wp-content/uploads/2018/02/Potential-temperature.pdf>`__
   - `大气位温、相当位温、饱和相当位温、静力稳定度 <https://renqlsysu.github.io/2019/10/23/potential_temperature/>`__

   .. seealso::
       - `potential_temperature — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.potential_temperature.html>`__
       - `pot_temp - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/pot_temp.shtml>`__


.. py:function:: calc_static_stability(temper_data: xarray.DataArray, vertical_dim: str, vertical_units: str) -> xarray.DataArray

   Calculate the static stability within a vertical profile.

   .. math::
       \sigma = - T \frac{\partial \ln \theta}{\partial p}

   Parameters
   ----------
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.

   Returns
   -------
   Static stability (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Howard B. Bluestein. (1992). Synoptic-Dynamic Meteorology in Midlatitudes: Principles of Kinematics and Dynamics, Vol. 1

   .. seealso::
       - `static_stability - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/static_stability.shtml>`__
       - `static_stability — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.static_stability.html>`__
       - `Static stability parameters · Issue #2535 · Unidata/MetPy <https://github.com/Unidata/MetPy/issues/2535>`__


.. py:function:: find_dims_axis(data: xarray.DataArray, dim: str) -> int

   Find the index of `dim` in the xarray DataArray.

   Parameters
   ----------
   - data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   - dim : :py:class:`str<python.str>`
       Dimension(s) over which to find axis.

   Returns
   -------
   :py:class:`int<python.int>`.


.. py:function:: transfer_inf2nan(ds: xarray.DataArray) -> xarray.DataArray

   Convert `np.inf` in `ds` to `np.nan`, respectively.

   Parameters
   ----------
   - ds: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Data include `np.inf`.

   Returns
   -------
   - Data include `np.nan`.: :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: transfer_data_units(input_data: xr.DataArray | xr.Dataset, input_units: str, output_units: str) -> xr.DataArray | xr.Dataset

   Data unit conversion


.. py:function:: get_weighted_spatial_data(data_input: xarray.DataArray, lat_dim: str = 'lat', lon_dim: str = 'lon', method: str = 'cos_lat') -> xarray.DataArray

   Get the area-weighting data.

   Parameters
   ----------
   - data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   - lat_dim: :py:class:`str<python.str>`.
       Latitude dimension over which to apply. By default is applied over the `lat` dimension.
   - lon_dim: :py:class:`str<python.str>`.
       Longitude dimension over which to apply. By default is applied over the `lon` dimension.
   - method: {`'cos_lat'`, `'area'`}.
       area-weighting methods.

       1. `'cos_lat'`: weighting data by the cosine of latitude.
       2. `'area'`: weighting data by area, where you weight each data point by the area of each grid cell.

   .. Caution:: 
       - `data_input` must be **regular lonlat grid**.
       - If you are calculating global average temperature just on land, 
         then you need to mask out the ocean in your area dataset at first.

   .. seealso::
       - `The Correct Way to Average the Globe (Why area-weighting your data is important) <https://towardsdatascience.com/the-correct-way-to-average-the-globe-92ceecd172b7>`__.
       - Kevin Cowtan, Peter Jacobs, Peter Thorne, Richard Wilkinson, 
         Statistical analysis of coverage error in simple global temperature estimators, 
         Dynamics and Statistics of the Climate System, Volume 3, Issue 1, 2018, dzy003, https://doi.org/10.1093/climsys/dzy003.


.. py:function:: generate_dataset_dispatcher(func)

   Function Dispensers: Iterate over the variables in the `xarray.Dataset` data using a function that only supports `xarray.DataArray` data


.. py:function:: calc_lon_gradient(data_input: xr.DataArray | xr.Dataset, lon_dim='lon', lat_dim='lat', min_dx=1.0, edge_order=2, R=6370000) -> xr.DataArray | xr.Dataset

   Calculate the gradient along the longitude.

   .. math::
       \frac{\partial F}{\partial x} = \frac{1}{R \cos\varphi} \cdot \frac{\partial F}{\partial \lambda}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   min_dx: :py:class:`float<python.float>`, default: `1.0`.
       The minimum acceptable value of `dx`, below which parts will set `nan` to avoid large computational errors. 
       The unit is m. You can set it to a negative value in order to remove this benefit.
   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
   R: :py:class:`float<python.float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The gradient along the longitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>`


.. py:function:: calc_lat_gradient(data_input: xr.DataArray | xr.Dataset, lat_dim='lat', min_dy=1.0, edge_order=2, R=6370000) -> xr.DataArray | xr.Dataset

   Calculate the gradient along the latitude.

   .. math::
       \frac{\partial F}{\partial y} = \frac{1}{R} \cdot \frac{\partial F}{\partial \varphi}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   min_dy: :py:class:`float<python.float>`, default: `1.0`.
       The minimum acceptable value of `dy`, below which parts will set `nan` to avoid large computational errors. 
       The unit is m. You can set it to a negative value in order to remove this benefit.
   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
   R: :py:class:`float<python.float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The gradient along the latitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>`  


.. py:function:: calc_lon_laplacian(data_input: xr.DataArray | xr.Dataset, lon_dim='lon', lat_dim='lat', min_dx2=1000000000.0, edge_order=2, R=6370000) -> xr.DataArray | xr.Dataset

   Calculation of the second-order partial derivative term (Laplace term) along longitude.

   .. math::
       \frac{\partial^2 F}{\partial x^2} = \frac{1}{(R \cos\varphi)^2} \cdot \frac{\partial^2 F}{\partial \lambda^2}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   min_dx2: :py:class:`float<python.float>`, default: `1e9`.
       The minimum acceptable value of :math:`(\mathrm{d}x)^2`, below which parts will set `nan` to avoid large computational errors. 
       The unit is m. You can set it to a negative value in order to remove this benefit.
   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
   R: :py:class:`float<python.float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The second-order partial derivative term (Laplace term) along longitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>` 


.. py:function:: calc_lat_laplacian(data_input: xr.DataArray | xr.Dataset, lat_dim='lat', min_dy2=1.0, edge_order=2, R=6370000) -> xr.DataArray | xr.Dataset

   Calculation of the second-order partial derivative term (Laplace term) along latitude.

   .. math::
       \frac{\partial^2 F}{\partial y^2} = \frac{1}{R^2} \cdot \frac{\partial^2 F}{\partial \varphi^2}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.            
   min_dy2: :py:class:`float<python.float>`, default: `1.0`.
       The minimum acceptable value of :math:`(\mathrm{d}y)^2`, below which parts will set `nan` to avoid large computational errors. 
       The unit is m. You can set it to a negative value in order to remove this benefit.
   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
   R: :py:class:`float<python.float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The second-order partial derivative term (Laplace term) along latitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>`             


.. py:function:: calc_lon_lat_mixed_derivatives(data_input: xr.DataArray | xr.Dataset, lon_dim='lon', lat_dim='lat', min_dxdy=10000000000.0, edge_order=2, R=6370000) -> xr.DataArray | xr.Dataset

   Calculation of second-order mixed partial derivative terms along longitude and latitude.

   .. math::
       \frac{\partial^2 F}{\partial x \partial y} = \frac{1}{R^2 \cos\varphi} \cdot \frac{\partial^2 F}{\partial \lambda \partial \varphi}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   min_dxdy: :py:class:`float<python.float>`, default: `1e10`.
       The minimum acceptable value of :math:`\mathrm{d}x\mathrm{d}y`, below which parts will set `nan` to avoid large computational errors. 
       The unit is m. You can set it to a negative value in order to remove this benefit.
   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
   R: :py:class:`float<python.float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The second-order mixed partial derivative terms along longitude and latitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>` 


.. py:function:: calc_time_gradient(data_input: xarray.DataArray, time_units: str, time_dim='time') -> xarray.DataArray

   Calculate the gradient along the time direction.

   .. math::
       \frac{\partial F}{\partial t}    

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   time_units: :py:class:`str<python.str>`.
       The unit corresponding to the time dimension value. Optional values are `seconds`, `months`, `years` and so on.
   time_dim: :py:class:`str<python.str>`, default: `time`.
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
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   surface_pressure_data_units: :py:class:`str<python.str>`.
       The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.

   Returns
   -------
   The pressure layer thickness (delta pressure) of a constant pressure level coordinate system (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       - :py:func:`geocat.comp.meteorology.delta_pressure <geocat-comp:geocat.comp.meteorology.delta_pressure>`
       - `dpres_plevel - NCL <https://www.ncl.ucar.edu/Document/Functions/Built-in/dpres_plevel.shtml>`__


.. py:function:: calc_p_integral(data_input: xarray.DataArray, vertical_dim: str, normalize=True) -> xarray.DataArray

   Calculate the vertical integral along the barometric pressure direction in the p-coordinate system.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The spatio-temporal data to be calculated.
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   normalize: :py:class:`bool<python.bool>`, default: `True`.
       Whether or not the integral results are averaged over the entire layer.

   Returns
   -------
   The vertical integral along the barometric pressure direction in the p-coordinate system (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. attention::
       This method ignores the effect of topography, so it applies to altitudes **above 900hPa** and is **NOT applicable to the Tibetan Plateau region**. 
       For a fully accurate vertical integration, please use the :py:func:`calc_top2surface_integral <calc_top2surface_integral>` function to calculate, 
       but the speed of the calculation is slightly slowed down.


.. py:function:: calc_top2surface_integral(data_input: xarray.DataArray, surface_pressure_data: xarray.DataArray, vertical_dim: str, surface_pressure_data_units: str, vertical_dim_units: str, method='Trenberth-vibeta', normalize=True) -> xarray.DataArray

   Calculate the vertical integral in the p-coordinate system from the ground to the zenith along the barometric pressure direction.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The spatio-temporal data to be calculated.
   surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Mean surface sea level pressure.    
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   surface_pressure_data_units: :py:class:`str<python.str>`.
       The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   method: :py:class:`str<python.str>`, default: `'Trenberth-vibeta'`.
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

   normalize: :py:class:`bool<python.bool>`, default: `True`.
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


.. py:function:: calc_laplacian(data_input: xarray.DataArray, lon_dim='lon', lat_dim='lat', R=6370000, spherical_coord=True) -> xarray.DataArray

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
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   R: :py:class:`float<python.float>`, default: `6370000`.
       Radius of the Earth.
   spherical_coord: :py:class:`bool<python.bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates.

   Returns
   -------
   The horizontal Laplace term. (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_divergence(u_data: xarray.DataArray, v_data: xarray.DataArray, lon_dim='lon', lat_dim='lat', R=6370000, spherical_coord=True) -> xarray.DataArray

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
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   R: :py:class:`float<python.float>`, default: `6370000`.
       Radius of the Earth.
   spherical_coord: :py:class:`bool<python.bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates.

   Returns
   -------
   The horizontal divergence term. (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_vorticity(u_data: xarray.DataArray, v_data: xarray.DataArray, lon_dim='lon', lat_dim='lat', R=6370000, spherical_coord=True) -> xarray.DataArray

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
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   R: :py:class:`float<python.float>`, default: `6370000`.
       Radius of the Earth.
   spherical_coord: :py:class:`bool<python.bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates.

   Returns
   -------
   The horizontal relative vorticity term. (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_geostrophic_wind(z_data: xarray.DataArray, lon_dim='lon', lat_dim='lat', omega=7.292e-05, g=9.8, R=6370000) -> xarray.DataArray

   Calculate the geostrophic wind.

   .. math::
       u_g = - \frac{g}{f} \frac{\partial H}{\partial y}

   .. math::
       v_g = \frac{g}{f} \frac{\partial H}{\partial x}

   Parameters
   ----------
   z_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric geopotential height.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   omega: :py:class:`float<python.float>`, default: `7.292e-5`.
       The angular speed of the earth.
   g: :py:class:`float<python.float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float<python.float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The geostrophic wind term. (:py:class:`xarray.DataArray<xarray.DataArray>`).
       - ug
       - vg


.. py:function:: calc_geostrophic_wind_vorticity(z_data: xarray.DataArray, lon_dim='lon', lat_dim='lat', spherical_coord=True, omega=7.292e-05, g=9.8, R=6370000) -> xarray.DataArray

   Calculate the geostrophic vorticity.

   rectangular coordinates

   .. math::
       \zeta_g = \frac{\partial v_g}{\partial x} - \frac{\partial u_g}{\partial y}

   Spherical coordinates

   .. math::
       \zeta_g = \frac{\partial v_g}{\partial x} - \frac{\partial u_g}{\partial y} + \frac{u_g}{R} \tan \varphi

   Parameters
   ----------
   z_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric geopotential height.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   spherical_coord: :py:class:`bool<python.bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates.
   omega: :py:class:`float<python.float>`, default: `7.292e-5`.
       The angular speed of the earth.
   g: :py:class:`float<python.float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float<python.float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The geostrophic vorticity term. (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_horizontal_water_flux(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, g=9.8) -> xarray.Dataset

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
   g: :py:class:`float<python.float>`, default: `9.8`.
       The acceleration of gravity.    

   Returns
   -------
   The water vapor flux. (:py:class:`xarray.Dataset<xarray.Dataset>`).

   - :math:`qu`: zonal water vapor flux.
   - :math:`qv`: meridional water vapor flux.


.. py:function:: calc_vertical_water_flux(specific_humidity_data: xarray.DataArray, omega_data: xarray.DataArray, g=9.8) -> xarray.Dataset

   Calculate vertical water vapor flux.

   .. math::
       -\omega \frac{q}{g}

   Parameters
   ----------
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   omega_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The vertical velocity data (:math:`\frac{\mathrm{d} p}{\mathrm{d} t}`).
   g: :py:class:`float<python.float>`, default: `9.8`.
       The acceleration of gravity.  


.. py:function:: calc_water_flux_top2surface_integral(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, surface_pressure_data: xarray.DataArray, surface_pressure_data_units: str, vertical_dim: str, vertical_dim_units: str, method='Trenberth-vibeta', g=9.8) -> xarray.DataArray

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
   surface_pressure_data_units: :py:class:`str<python.str>`.
       The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   method: :py:class:`str<python.str>`, default: `'Trenberth-vibeta'`.
       Vertical integration method. Optional values are `Boer-vibeta`, `'Trenberth-vibeta'`.
   g: :py:class:`float<python.float>`, default: `9.8`.
       The acceleration of gravity.    

   Returns
   -------
   The water vapor flux. (:py:class:`xarray.Dataset<xarray.Dataset>`).

   - :math:`qu`: zonal water vapor flux.
   - :math:`qv`: meridional water vapor flux.

   .. seealso::
       :py:func:`calc_top2surface_integral <calc_top2surface_integral>` 


.. py:function:: calc_divergence_watervaporflux(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, specific_humidity_units: str, spherical_coord=True, lon_dim='lon', lat_dim='lat', g=9.8, R=6370000) -> xarray.DataArray

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
   specific_humidity_units: :py:class:`str<python.str>`.
       The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
   spherical_coord: :py:class:`bool<python.bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   g: :py:class:`float<python.float>`, default: `9.8`.
       The acceleration of gravity. 
   R: :py:class:`float<python.float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The water vapor flux divergence. (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_divergence_watervaporflux_top2surface_integral(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, surface_pressure_data: xarray.DataArray, vertical_dim: str, specific_humidity_units: str, surface_pressure_data_units: str, vertical_dim_units: str, spherical_coord=True, lon_dim='lon', lat_dim='lat', method='Trenberth-vibeta', g=9.8, R=6370000) -> xarray.DataArray

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
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   specific_humidity_units: :py:class:`str<python.str>`.
       The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
   surface_pressure_data_units: :py:class:`str<python.str>`.
       The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   spherical_coord: :py:class:`bool<python.bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   g: :py:class:`float<python.float>`, default: `9.8`.
       The acceleration of gravity. 
   R: :py:class:`float<python.float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The water vapor flux divergence. (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_u_advection(u_data: xarray.DataArray, temper_data: xarray.DataArray, lon_dim='lon', lat_dim='lat') -> xarray.DataArray

   Calculate zonal temperature advection at each vertical level.

   .. math::
       -u \frac{\partial T}{\partial x}

   Parameters
   ----------
   u: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The zonal temperature advection. (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_v_advection(v_data: xarray.DataArray, temper_data: xarray.DataArray, lon_dim='lon', lat_dim='lat') -> xarray.DataArray

   Calculate meridional temperature advection at each vertical level.

   .. math::
       -v \frac{\partial T}{\partial y}

   Parameters
   ----------
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The meridional temperature advection. (:py:class:`xarray.DataArray<xarray.DataArray>`).


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
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.

   Returns
   -------
   The vertical temperature transport. (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: transfer_dFdp2dFdz(dFdp_data: xr.DataArray | xr.Dataset, rho_d: float = 1292.8, g: float = 9.8)

   The transformation relationship between the z coordinate system and the p coordinate system.

   .. math::
       \frac{\partial F}{\partial z} = \frac{\partial F}{\partial p} \frac{\partial p}{\partial z} = - \rho g \frac{\partial F}{\partial p}


.. py:function:: calc_eady_growth_rate(u_daily_data: xarray.DataArray, z_daily_data: xarray.DataArray, temper_daily_data: xarray.DataArray, vertical_dim: str, vertical_units: str, lat_dim='lat', g=9.8) -> xarray.Dataset

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
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Daily air temperature.
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   g: :py:class:`float<python.float>`, default: `9.8`.
       The acceleration of gravity. 

   Returns
   -------
   The maximum Eady growth rate. (:py:class:`xarray.Dataset<xarray.Dataset>`).

   - `eady_growth_rate`: The maximum Eady growth rate.
   - `dudz`: :math:`\frac{\mathrm{d} U}{\mathrm{d} z}`
   - `brunt_vaisala_frequency`: Brunt-väisälä frequency.

   .. seealso::
       - `eady_growth_rate -NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/eady_growth_rate.shtml>`__
       - `瞬变涡旋诊断量 <https://renqlsysu.github.io/2020/02/16/wave_activity_flux/>`__


.. py:function:: calc_apparent_heat_source(u_data: xarray.DataArray, v_data: xarray.DataArray, omega_data: xarray.DataArray, temper_data: xarray.DataArray, vertical_dim: str, vertical_units: str, time_units: str, lon_dim='lon', lat_dim='lat', time_dim='time', c_p=1005.7) -> xarray.DataArray

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
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   time_units: :py:class:`str<python.str>`.
       The unit corresponding to the time dimension value. Optional values are `seconds`, `months`, `years` and so on.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   time_dim: :py:class:`str<python.str>`.
       The time coordinate dimension name.
   c_p: :py:class:`float<python.float>`, default: `1005.7`.
       The specific heat at constant pressure of dry air.

       .. note::
           `specific heat capacity - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Specific_heat_capacity>`__

   Returns
   -------
   The apparent heat source (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::
       - `Yanai, M., & Tomita, T. (1998). Seasonal and Interannual Variability of Atmospheric Heat Sources and Moisture Sinks as Determined from NCEP–NCAR Reanalysis, Journal of Climate, 11(3), 463-482. <https://journals.ametsoc.org/view/journals/clim/11/3/1520-0442_1998_011_0463_saivoa_2.0.co_2.xml>`__
       - `Ling, J., & Zhang, C. (2013). Diabatic Heating Profiles in Recent Global Reanalyses, Journal of Climate, 26(10), 3307-3325. <https://doi.org/10.1175/JCLI-D-12-00384.1>`__


.. py:function:: calc_total_diabatic_heating(u_data: xarray.DataArray, v_data: xarray.DataArray, omega_data: xarray.DataArray, temper_data: xarray.DataArray, vertical_dim: str, vertical_units: str, time_units: str, lat_dim='lat', lon_dim='lon', time_dim='time', c_p=1005.7) -> xarray.DataArray

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
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   time_units: :py:class:`str<python.str>`.
       The unit corresponding to the time dimension value. Optional values are `seconds`, `months`, `years` and so on.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   time_dim: :py:class:`str<python.str>`, default: `time`.
       The time coordinate dimension name.
   c_p: :py:class:`float<python.float>`, default: `1005.7` (:math:`\mathrm{J \cdot kg^{-1} \cdot K^{-1}}`).
       The specific heat at constant pressure of dry air.

       .. note::
           `specific heat capacity - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Specific_heat_capacity>`__

   Returns
   -------
   The total diabatic heating (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::
       :py:func:`calc_apparent_heat_source <calc_apparent_heat_source>`


.. py:function:: calc_apparent_moisture_sink(u_data: xarray.DataArray, v_data: xarray.DataArray, omega_data: xarray.DataArray, specific_humidity_data: xarray.DataArray, vertical_dim: str, vertical_units: str, time_units: str, specific_humidity_units: str, lon_dim='lon', lat_dim='lat', time_dim='time', latent_heat_of_condensation=2501000.0)

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
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   time_units: :py:class:`str<python.str>`.
       The unit corresponding to the time dimension value. Optional values are `seconds`, `months`, `years` and so on.
   specific_humidity_units: :py:class:`str<python.str>`.
       The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   time_dim: :py:class:`str<python.str>`, default: `time`.
       The time coordinate dimension name.
   latent_heat_of_condensation: :py:class:`float<python.float>`, default: `2.5008e6` (:math:`\mathrm{J \cdot kg^{-1}}`).
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


.. py:function:: calc_Plumb_wave_activity_horizontal_flux(z_prime_data: xarray.DataArray, vertical_dim: str, vertical_units: str, lon_dim='lon', lat_dim='lat', omega=7.292e-05, g=9.8, R=6370000)

   Calculate Plumb wave activity horizontal flux.

   Parameters
   ----------
   z_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The anormaly of atmospheric geopotential height.
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   omega: :py:class:`float<python.float>`, default: `7.292e-5`.
       The angular speed of the earth.
   g: :py:class:`float<python.float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float<python.float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The Plumb wave activity horizontal flux (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::
       - `Plumb, R. A., 1985: On the Three-Dimensional Propagation of Stationary Waves. J. Atmos. Sci., 42, 217–229 <https://journals.ametsoc.org/view/journals/atsc/42/3/1520-0469_1985_042_0217_ottdpo_2_0_co_2.xml>`__


.. py:function:: calc_TN_wave_activity_horizontal_flux(z_prime_data, u_climatology_data, v_climatology_data, vertical_dim, vertical_dim_units, lon_dim='lon', lat_dim='lat', omega=7.292e-05, g=9.8, R=6370000)

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
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   omega: :py:class:`float<python.float>`, default: `7.292e-5`.
       The angular speed of the earth.
   g: :py:class:`float<python.float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float<python.float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The TN wave activity horizontal flux (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::
       - http://www.atmos.rcast.u-tokyo.ac.jp/nishii/programs/index.html
       - http://500hpa.cn/pyinmet/tnflux/
       - http://tytd.gx.cn/exchange/tnflux/
       - https://github.com/laishenggx/T-N_Wave-Activity-Flux


.. py:function:: calc_TN_wave_activity_3D_flux(z_prime_data, u_climatology_data, v_climatology_data, temper_data, vertical_dim, vertical_units, z_data=None, lon_dim='lon', lat_dim='lat', omega=7.292e-05, g=9.8, R=6370000, scale_height=8000, kappa=287 / 1005.7, method='practical_height')

   Calculate TN wave activity 3D flux.

   Parameters
   ----------
   z_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The anormaly of atmospheric geopotential height.
   u_climatology_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The climatology of zonal wind data.
   v_climatology_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The climatology of meridional wind data.
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   z_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric geopotential height.
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   omega: :py:class:`float<python.float>`, default: `7.292e-5`.
       The angular speed of the earth.
   g: :py:class:`float<python.float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float<python.float>`, default: `6370000`.
       Radius of the Earth.
   scale_height: :py:class:`float<python.float>`, default: `8000`.
       Scale height.
   kappa: :py:class:`float<python.float>`, default: `287/1005.7`.
       Poisson constant :math:`\kappa`.

       .. note::
           `Poisson constant - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Poisson_constant>`__

   method: :py:class:`str<python.str>`, default: `'practical_height'`.
       The calculation method of :math:`\mathrm{d}z`. Optional values are `'practical_height'`, `'scale_height'`.

   Returns
   -------
   The TN wave activity 3D flux (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_EP_horizontal_flux(u_prime_data, v_prime_data, time_dim='time', lat_dim='lat')

   Calculate horizontal Eliassen–Palm Flux.

   Parameters
   ----------
   u_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The anormaly of zonal wind data.
   v_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The anormaly of meridional wind data.
   time_dim: :py:class:`str<python.str>`, default: `time`.
       The time coordinate dimension name.        
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The Eliassen–Palm Flux (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::
       - https://www.ncl.ucar.edu/Applications/EPflux.shtml
       - https://renqlsysu.github.io/2020/02/16/wave_activity_flux/


.. py:function:: get_specific_years_data(data_input: easyclimate.core.yearstat.xr.DataArray, year_array: easyclimate.core.yearstat.np.array, dim='time') -> easyclimate.core.yearstat.xr.DataArray

   Slicing and extracting the part of the data containing the specified year based on an array of given integer years.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   year_array: numpy.array containing int data-type objects
       Year(s) to be extracted.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_months_data(data_input: easyclimate.core.yearstat.xr.DataArray, month_array: easyclimate.core.yearstat.np.array, dim='time') -> easyclimate.core.yearstat.xr.DataArray

   Slicing and extracting the part of the data containing the specified year based on an array of given integer months.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   month_array: numpy.array containing int data-type objects
       Month(s) to be extracted.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_days_data(data_input: easyclimate.core.yearstat.xr.DataArray, day_array: easyclimate.core.yearstat.np.array, dim='time') -> easyclimate.core.yearstat.xr.DataArray

   Slicing and extracting the part of the data containing the specified year based on an array of given integer days.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   day_array: numpy.array containing int data-type objects
       Days(s) to be extracted.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_hours_data(data_input: easyclimate.core.yearstat.xr.DataArray, hour_array: easyclimate.core.yearstat.np.array, dim='time') -> easyclimate.core.yearstat.xr.DataArray

   Slicing and extracting the part of the data containing the specified year based on an array of given integer hours.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   hour_array: numpy.array containing int data-type objects
       Hour(s) to be extracted.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_minutes_data(data_input: easyclimate.core.yearstat.xr.DataArray, minute_array: easyclimate.core.yearstat.np.array, dim='time') -> easyclimate.core.yearstat.xr.DataArray

   Slicing and extracting the part of the data containing the specified year based on an array of given integer minutes.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   minute_array: numpy.array containing int data-type objects
       Minute(s) to be extracted.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_seconds_data(data_input: easyclimate.core.yearstat.xr.DataArray, second_array: easyclimate.core.yearstat.np.array, dim='time') -> easyclimate.core.yearstat.xr.DataArray

   Slicing and extracting the part of the data containing the specified year based on an array of given integer seconds.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   second_array: numpy.array containing int data-type objects
       Second(s) to be extracted.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_microseconds_data(data_input: easyclimate.core.yearstat.xr.DataArray, microsecond_array: easyclimate.core.yearstat.np.array, dim='time') -> easyclimate.core.yearstat.xr.DataArray

   Slicing and extracting the part of the data containing the specified year based on an array of given integer microseconds.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   microsecond_array: numpy.array containing int data-type objects
       Microsecond(s) to be extracted.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_nanoseconds_data(data_input: easyclimate.core.yearstat.xr.DataArray, nanosecond_array: easyclimate.core.yearstat.np.array, dim='time') -> easyclimate.core.yearstat.xr.DataArray

   Slicing and extracting the part of the data containing the specified year based on an array of given integer nanoseconds.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   nanosecond_array: numpy.array containing int data-type objects
       Nanosecond(s) to be extracted.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_dayofweek_data(data_input: easyclimate.core.yearstat.xr.DataArray, dayofweek_array: easyclimate.core.yearstat.np.array, dim='time') -> easyclimate.core.yearstat.xr.DataArray

   Slicing and extracting the part of the data containing the specified year based on an array of given integer dayofweek.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   dayofweek_array: numpy.array containing int data-type objects
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

   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_yearmean_for_specific_months_data(data_input: easyclimate.core.yearstat.xr.DataArray, month_array: easyclimate.core.yearstat.np.array, dim='time', kwargs=None) -> easyclimate.core.yearstat.xr.DataArray

   Get the annual average of certain months.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   month_array: numpy.array containing int data-type objects
       Month(s) to be extracted.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data. 
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_year_exceed_index_upper_bound(data_input: easyclimate.core.yearstat.xr.DataArray, thresh: float, time_dim: str = 'time') -> easyclimate.core.yearstat.np.array

   Extract the years under the specified threshold (upper bound) in the annual average index (one-dimensional data with only a `time` dimension).

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The one-dimensional data with only a `time` dimension.
   thresh: :py:class:`float<python.float>`.
       The threshold value.
   time_dim: :py:class:`str<python.str>`.
       The time coordinate dimension name.

   Returns
   -------
   :py:class:`numpy.array <numpy:numpy.array>`.


.. py:function:: get_year_exceed_index_lower_bound(data_input: easyclimate.core.yearstat.xr.DataArray, thresh: float, time_dim: str = 'time') -> easyclimate.core.yearstat.np.array

   Extract the years under the specified threshold (lower bound) in the annual average index (one-dimensional data with only a `time` dimension).

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The one-dimensional data with only a `time` dimension.
   thresh: :py:class:`float<python.float>`.
       The threshold value.
   time_dim: :py:class:`str<python.str>`.
       The time coordinate dimension name.

   Returns
   -------
   :py:class:`numpy.array <numpy:numpy.array>`.


.. py:function:: calc_yearmean(data_input, dim='time', **kwargs)

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
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   dim : str
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


.. py:function:: calc_yearsum(data_input, dim='time', **kwargs)

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
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   dim : str
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


.. py:function:: calc_yearstd(data_input, dim='time', **kwargs)

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
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   dim : str
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


.. py:function:: calc_yearvar(data_input, dim='time', **kwargs)

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
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   dim : str
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


.. py:function:: calc_yearmax(data_input, dim='time', **kwargs)

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
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   dim : str
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


.. py:function:: calc_yearmin(data_input, dim='time', **kwargs)

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
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   dim : str
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


.. py:function:: assert_compared_version(ver1: float, ver2: float) -> int

   Compare python library versions.

   .. attention::
       - Only for incoming version numbers without alphabetic characters.
       - Based on this method, the version number comparison should result in the following `"10.12.2.6.5">"10.12.2.6"`.

   Parameters
   ----------
   - ver1: Version number 1
   - ver2: Version number 2

   Returns
   -------
   :py:class:`int<python.int>`.

   .. note::
       If `ver1<ver2`, return `-1`; If `ver1=ver2`, return `0`; If `ver1>ver2`, return `1`.

   Examples
   --------

   .. code:: python

       >>> import easyclimate as ecl
       >>> result = ecl.assert_compared_version("10.12.2.6.5", "10.12.2.6")
       >>> print(result)
       1


.. py:function:: transfer_int2datetime(data: numpy.array) -> numpy.datetime64

   Convert a numpy array of years of type integer to `np.datetime64` type.

   Parameters
   ----------
   - data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

   Examples
   --------

   .. code:: python

       >>> import easyclimate as ecl
       >>> import numpy as np
       >>> intyear = np.array([2054, 2061, 2062, 2067, 2071, 2075, 2076, 2078, 2085, 2089, 2096])
       >>> ecl.transfer_int2datetime(intyear)
       array(['2054-01-01T00:00:00.000000000', '2061-01-01T00:00:00.000000000',
              '2062-01-01T00:00:00.000000000', '2067-01-01T00:00:00.000000000',
              '2071-01-01T00:00:00.000000000', '2075-01-01T00:00:00.000000000',
              '2076-01-01T00:00:00.000000000', '2078-01-01T00:00:00.000000000',
              '2085-01-01T00:00:00.000000000', '2089-01-01T00:00:00.000000000',
              '2096-01-01T00:00:00.000000000'], dtype='datetime64[ns]')

   .. seealso::
       `Python(pandas)整数类型数据转换为时间类型 <https://www.jianshu.com/p/d12d95fbc90c>`__.


.. py:function:: transfer_datetime2int(ds: xarray.DataArray) -> xarray.DataArray

   Convert `np.datetime64` type with years and days to `year` and `day` coordinates.

   Parameters
   ----------
   - data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

   .. seealso::
       `Function in xarray to regroup monthly data into months and # of years <https://github.com/pydata/xarray/discussions/5119>`__.


.. py:function:: transfer_monmean2everymonthmean(data_input: xarray.DataArray, time_dim: str = 'time') -> xarray.DataArray

   Convert to the month-mean state corresponding to each month.

   Parameters
   ----------
   - data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.    


.. py:function:: get_compress_xarraydata(data: xr.DataArray | xr.Dataset, complevel: int) -> xr.DataArray | xr.Dataset

   Export compressible netCDF files from xarray data (:py:class:`xarray.DataArray<xarray.DataArray>`, :py:class:`xarray.Dataset<xarray.Dataset>`)


.. py:function:: sort_ascending_latlon_coordinates(data: xr.DataArray | xr.Dataset, lat_dim: str = 'lat', lon_dim: str = 'lon') -> xr.DataArray | xr.Dataset

   Sort the dimensions `lat`, `lon` in ascending order.


.. py:function:: generate_datatree_dispatcher(func)

   Function Dispensers: Iterate over the variables in the `xarray.Dataset` data using a function that only supports `xarray.DataArray` data


.. py:function:: transfer_xarray_lon_from180TO360(data_input: xr.DataArray | xr.Dataset, lon_dim: str = 'lon') -> xr.DataArray | xr.Dataset

   Longitude conversion -180-180 to 0-360.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.

   .. seealso::
       :py:func:`transfer_xarray_lon_from360TO180 <transfer_xarray_lon_from360TO180>`


.. py:function:: transfer_xarray_lon_from360TO180(data_input: xr.DataArray | xr.Dataset, lon_dim: str = 'lon') -> xr.DataArray | xr.Dataset

   Longitude conversion 0-360 to -180-180.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.

   .. seealso::
       :py:func:`transfer_xarray_lon_from180TO360 <transfer_xarray_lon_from180TO360>`


.. py:function:: module_available(module: str) -> bool

   Checks whether a module is installed without importing it.

   Use this for a lightweight check and lazy imports.

   Parameters
   ----------
   module : str
       Name of the module.

   Returns
   -------
   available : bool
       Whether the module is installed.


.. py:function:: open_muliti_dataset(files: str, dim: str, **kwargs) -> xarray.Dataset

   Compare python library versions.

   .. attention::
       - Only for incoming version numbers without alphabetic characters.
       - Based on this method, the version number comparison should result in the following `"10.12.2.6.5">"10.12.2.6"`.

   Parameters
   ----------
   - ver1: Version number 1
   - ver2: Version number 2

   Returns
   -------
   :py:class:`int<python.int>`.

   .. note::
       If `ver1<ver2`, return `-1`; If `ver1=ver2`, return `0`; If `ver1>ver2`, return `1`.

   Examples
   --------

   .. code:: python
       >>> import easyclimate as ecl
       >>> result = assert_compared_version("10.12.2.6.5", "10.12.2.6")
       >>> print(result)
       1

   .. todo::
       - https://medium.com/pangeo/accessing-netcdf-and-grib-file-collections-as-cloud-native-virtual-datasets-using-kerchunk-625a2d0a9191
       - https://github.com/fsspec/kerchunk/issues/240


.. py:function:: calc_linregress_spatial(data_input, dim='time', x=None, alternative='two-sided', returns_type='dataset_returns', engine='scipy_linregress')

   Calculate a linear least-squares regression for spatial data of time.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>` to be regression.
   dim : str
       Dimension(s) over which to apply linregress. By default linregress is applied over the `time` dimension.
   x : numpy.array
   returns_type: str

   Returns
   -------
   result : ``LinregressResult`` Dataset
       The return Dataset have following data_var:

       slope : float
           Slope of the regression line.
       intercept : float
           Intercept of the regression line.
       rvalue : float
           The Pearson correlation coefficient. The square of ``rvalue``
           is equal to the coefficient of determination.
       pvalue : float
           The p-value for a hypothesis test whose null hypothesis is
           that the slope is zero, using Wald Test with t-distribution of
           the test statistic. See `alternative` above for alternative
           hypotheses.
       stderr : float
           Standard error of the estimated slope (gradient), under the
           assumption of residual normality.
       intercept_stderr : float
           Standard error of the estimated intercept, under the assumption
           of residual normality.

   .. seealso::
       :py:func:`scipy.stats.linregress <scipy:scipy.stats.linregress>`.


.. py:function:: calc_detrend_data(data_input, time_dim='time')

   Remove linear trend along axis from data.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        The spatio-temporal data of :py:class:`xarray.DataArray<xarray.DataArray>` to be detrended.
   dim : `str`
       Dimension(s) over which to detrend. By default dimension is applied over the `time` dimension.

   Returns
   -------
   - :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::
       :py:func:`scipy.signal.detrend <scipy:scipy.signal.detrend>`.


.. py:function:: calc_ttestSpatialPattern_spatial(data_input1, data_input2, dim='time')

   Calculate the T-test for the means of two independent sptial samples along with other axis (i.e. 'time') of scores.

   Parameters
   ----------
   data_input1 : :py:class:`xarray.DataArray<xarray.DataArray>`
        The first spatio-temporal data of xarray DataArray to be calculated.
   data_input2 : :py:class:`xarray.DataArray<xarray.DataArray>`
        The second spatio-temporal data of xarray DataArray to be calculated.

   .. note::
       The order of `data_input1` and `data_input2` has no effect on the calculation result.

   dim : `str`
       Dimension(s) over which to apply skewness. By default skewness is applied over the `time` dimension.

   Returns
   -------
   - statistic, pvalue: :py:class:`xarray.Dataset<xarray.Dataset>`.

   .. seealso::
       :py:func:`scipy.stats.ttest_ind <scipy:scipy.stats.ttest_ind>`.


.. py:function:: calc_skewness_spatial(data_input, dim='time')

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
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        The spatio-temporal data of xarray DataArray to be calculated.
   dim : str
       Dimension(s) over which to apply skewness. By default skewness is applied over the `time` dimension.

   Returns
   -------
   - skewness, pvalue: :py:class:`xarray.Dataset<xarray.Dataset>`.

   Reference
   --------------
   White, G. H. (1980). Skewness, Kurtosis and Extreme Values of 
   Northern Hemisphere Geopotential Heights, Monthly Weather Review, 108(9), 1446-1455. 
   Website: https://journals.ametsoc.org/view/journals/mwre/108/9/1520-0493_1980_108_1446_skaevo_2_0_co_2.xml

   .. seealso::
       :py:func:`scipy.stats.skew <scipy:scipy.stats.skew>`, :py:func:`scipy.stats.normaltest <scipy:scipy.stats.normaltest>`.


.. py:function:: calc_kurtosis_spatial(data_input, dim='time')

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
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        The spatio-temporal data of xarray DataArray to be calculated.
   dim : str
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


.. py:function:: get_EOF_model(data_input: xarray.DataArray, time_dim: str = 'time', n_modes=10, standardize=False, use_coslat=False, use_weights=False, weights=None, solver='auto', **solver_kwargs)

       
       


.. py:function:: calc_EOF_analysis(model: xeofs.models.eof.EOF, mini_summary=False)


.. py:function:: get_EOF_projection(model: xeofs.models.eof.EOF, data: xarray.DataArray)


.. py:function:: save_EOF_model(model: xeofs.models.eof.EOF, path: save_EOF_model.str)

       
       


.. py:function:: load_EOF_model(path: str)

       
       


.. py:function:: get_MCA_model(data_input1: xarray.DataArray, data_input2: xarray.DataArray, time_dim: str = 'time', n_modes=10, standardize=False, use_coslat=False, use_weights=False, weights1=None, weights2=None, n_pca_modes=None, solver='auto', **solver_kwargs)

       
       


.. py:function:: calc_MCA_analysis(model: xeofs.models.mca.MCA, correction=None, alpha=0.05, mini_summary=False)


.. py:function:: get_MCA_projection(model: xeofs.models.mca.MCA, **kwargs)

       
       


.. py:function:: save_MCA_model(model: xeofs.models.mca.MCA, path: save_MCA_model.str)

       
       


.. py:function:: load_MCA_model(path: str)

       
       


.. py:function:: calc_index_NPWI(precipitable_water_daily: xarray.DataArray, time_dim='time') -> xarray.DataArray

   Calculate the normalized precipitable water index (NPWI).

   .. math::
       \mathrm{NPWI} = \frac{\mathrm{PW} - \mathrm{PW_{min}}}{\mathrm{PW_{max}} - \mathrm{PW_{min}}}

   Parameters
   ----------
   precipitable_water_daily: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Daily precipitable water data. 
   time_dim: :py:class:`str<python.str>`.
       The time coordinate dimension name.

   Returns
   -------
   Normalized precipitable water index (NPWI).

   Reference
   --------------
   Zeng, X., and E. Lu, 2004: Globally Unified Monsoon Onset and Retreat Indexes. J. Climate, 17, 2241–2248, https://doi.org/10.1175/1520-0442(2004)017<2241:GUMOAR>2.0.CO;2.


.. py:function:: find_PW_monsoon_region(precipitable_water_daily: xarray.DataArray, time_dim='time') -> xarray.DataArray

   The refined monsoon regions.

   .. note::
       To refine the definition of monsoon regions on a grid- cell-by-cell basis, 
       we first compute the 10-yr-averaged monthly PW over each cell. 
       Then we obtain the maximum monthly PW during the three summer months [e.g., June–August in the Northern Hemisphere (NH), denoted as PWw], 
       and the maximum monthly PW during the three winter months (e.g., December–February for the NH, denoted as PWc). 
       The refined monsoon regions are simply defined as grid cells that are within the monsoon 
       regions given in the above studies and have a difference between PWw and PWc greater than 12 mm. 
       Initially we have also tried to use the annual maximum and minimum monthly PW values.

   Parameters
   ----------
   precipitable_water_daily: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Daily precipitable water data.
   time_dim: :py:class:`str<python.str>`.
       The time coordinate dimension name.

   Reference
   --------------
   Zeng, X., and E. Lu, 2004: Globally Unified Monsoon Onset and Retreat Indexes. J. Climate, 17, 2241–2248, https://doi.org/10.1175/1520-0442(2004)017<2241:GUMOAR>2.0.CO;2.


.. py:function:: cal_NPWI_monsoon_onset(NPWI, thresh=0.618, consecutive_days=3, n=7, lon_dim='lon', lat_dim='lat', time_dim='time') -> xarray.DataArray

   Calculate the summer monsoon onset date.

   The summer monsoon onset date for grid cell G is defined as the first day (:math:`d`) 
   when NWPI is greater than the Golden Ratio (0.618) for three consecutive days
   in seven of the nine cells centered at cell G in day :math:`d` or (:math:`d \pm 1`).

   .. note::
       If one or more of the nine grids are undefined, for example, at the edge of monsoon regions, 
       the required number of seven is correspondingly reduced. 
       For instance, if only seven grid cells are defined, the required number is five.

   Parameters
   ----------
   NPWI: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Normalized precipitable water index (NPWI). 

       .. attention::
           It must include three dimensions: `time`, `longitude`, and `latitude`.

   thresh: :py:class:`float<python.float>`, default: `0.618`.
       Golden Ratio value for the threshold value.
   consecutive_days: :py:class:`int<python.int>`, default: `3`.
       Consecutive days values.
   n: :py:class:`int<python.int>`, default: `7`.
       :math:`n` of the nine cells centered at cell G in day :math:`d` or (:math:`d \pm 1`).
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   time_dim: :py:class:`str<python.str>`.
       The time coordinate dimension name.

   Returns
   -------
   Summer monsoon onset date.

   Reference
   --------------
   Zeng, X., and E. Lu, 2004: Globally Unified Monsoon Onset and Retreat Indexes. J. Climate, 17, 2241–2248, https://doi.org/10.1175/1520-0442(2004)017<2241:GUMOAR>2.0.CO;2.


.. py:function:: cal_NPWI_monsoon_detreat(NPWI, monsoon_onset_date, thresh=0.618, consecutive_days=3, n=7, lon_dim='lon', lat_dim='lat', time_dim='time') -> xarray.DataArray

   Calculate the summer monsoon retreat date.

   The summer monsoon retreat date for grid cell G is defined as the first day (:math:`d`) 
   when NWPI is less than the Golden Ratio (0.618) for three consecutive days
   in seven of the nine cells centered at cell G in day :math:`d` or (:math:`d \pm 1`).

   .. note::
       If one or more of the nine grids are undefined, for example, at the edge of monsoon regions, 
       the required number of seven is correspondingly reduced.
       For instance, if only seven grid cells are defined, the required number is five.

   Parameters
   ----------
   NPWI: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Normalized precipitable water index (NPWI). 

       .. attention::
           It must include three dimensions: `time`, `longitude`, and `latitude`.

   monsoon_onset_date: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Summer monsoon onset date. The results is generated by :py:func:`easyclimate.index.cal_NPWI_monsoon_onset <easyclimate.index.cal_NPWI_monsoon_onset>`.
   thresh: :py:class:`float<python.float>`, default: `0.618`.
       Golden Ratio value for the threshold value.
   consecutive_days: :py:class:`int<python.int>`, default: `3`.
       Consecutive days values.
   n: :py:class:`int<python.int>`, default: `7`.
       :math:`n` of the nine cells centered at cell G in day :math:`d` or (:math:`d \pm 1`).
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   time_dim: :py:class:`str<python.str>`.
       The time coordinate dimension name.

   Returns
   -------
   Summer monsoon retreat date.

   Reference
   --------------
   Zeng, X., and E. Lu, 2004: Globally Unified Monsoon Onset and Retreat Indexes. J. Climate, 17, 2241–2248, https://doi.org/10.1175/1520-0442(2004)017<2241:GUMOAR>2.0.CO;2.


.. py:function:: get_weighted_spatial_data(data_input: xarray.DataArray, lat_dim: str = 'lat', lon_dim: str = 'lon', method: str = 'cos_lat') -> xarray.DataArray

   Get the area-weighting data.

   Parameters
   ----------
   - data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   - lat_dim: :py:class:`str<python.str>`.
       Latitude dimension over which to apply. By default is applied over the `lat` dimension.
   - lon_dim: :py:class:`str<python.str>`.
       Longitude dimension over which to apply. By default is applied over the `lon` dimension.
   - method: {`'cos_lat'`, `'area'`}.
       area-weighting methods.

       1. `'cos_lat'`: weighting data by the cosine of latitude.
       2. `'area'`: weighting data by area, where you weight each data point by the area of each grid cell.

   .. Caution:: 
       - `data_input` must be **regular lonlat grid**.
       - If you are calculating global average temperature just on land, 
         then you need to mask out the ocean in your area dataset at first.

   .. seealso::
       - `The Correct Way to Average the Globe (Why area-weighting your data is important) <https://towardsdatascience.com/the-correct-way-to-average-the-globe-92ceecd172b7>`__.
       - Kevin Cowtan, Peter Jacobs, Peter Thorne, Richard Wilkinson, 
         Statistical analysis of coverage error in simple global temperature estimators, 
         Dynamics and Statistics of the Climate System, Volume 3, Issue 1, 2018, dzy003, https://doi.org/10.1093/climsys/dzy003.


.. py:function:: sort_ascending_latlon_coordinates(data: xr.DataArray | xr.Dataset, lat_dim: str = 'lat', lon_dim: str = 'lon') -> xr.DataArray | xr.Dataset

   Sort the dimensions `lat`, `lon` in ascending order.


.. py:function:: calc_intensity_STFZ(data_sst_DtDy: xarray.DataArray, criteria=0.45 * 1e-05, lat_range=[24, 32], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray

   Calculate the intensity of the subtropical frontal zone (STFZ). 
   The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

   .. math::
       \mathrm{ITS} = \sum_{i=1}^{N} \frac{G_i}{N}

   where :math:`G_i` is the value of zonally-averaged SST meridional gradient that is no less than 
   an empirically-given critical value (here, :math:`0.45 \times 10^{-5} \mathrm{km^{-1}}` for STFZ) at the :math:`i`-th latitudinal grid point within the zone, 
   and :math:`N` is the number of total grid points that satisfy the criteria above.

   Parameters
   ----------
   data_sst_DtDy : :py:class:`xarray.DataArray<xarray.DataArray>` 
       The SST meridional gradient data.
   criteria: :py:class:`float<python.float>`, default: `0.45 *1e-5`.
       Empirically-given critical value.
   lat_range: :py:class:`list<python.list>`, default: `[24, 32]`.
       The latitude range of the oceanic frontal zone.
   lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
       The longitude range of the oceanic frontal zone.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   The intensity of the STFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017), 
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts 
   and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766


.. py:function:: calc_intensity_SAFZ(data_sst_DtDy: xarray.DataArray, criteria=0.8 * 1e-05, lat_range=[36, 44], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray

   Calculate the intensity of the subarctic frontal zone (SAFZ). 
   The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

   .. math::
       \mathrm{ITS} = \sum_{i=1}^{N} \frac{G_i}{N}

   where :math:`G_i` is the value of zonally-averaged SST meridional gradient that is no less than 
   an empirically-given critical value (here, :math:`0.80 \times 10^{-5} \mathrm{km^{-1}}` for SAFZ) at the :math:`i`-th latitudinal grid point within the zone, 
   and :math:`N` is the number of total grid points that satisfy the criteria above.

   Parameters
   ----------
   data_sst_DtDy : :py:class:`xarray.DataArray<xarray.DataArray>` 
       The SST meridional gradient data.
   criteria: :py:class:`float<python.float>`, default: `0.80 *1e-5`.
       Empirically-given critical value.
   lat_range: :py:class:`list<python.list>`, default: `[36, 44]`.
       The latitude range of the oceanic frontal zone.
   lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
       The longitude range of the oceanic frontal zone.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   The intensity of the SAFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017), 
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts 
   and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766


.. py:function:: calc_location_STFZ(data_sst_DtDy: xarray.DataArray, criteria=0.45 * 1e-05, lat_range=[24, 32], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray

   Calculate the location index of the subtropical frontal zone (STFZ). 
   The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

   .. math::
       \mathrm{LCT} = \sum_{i=1}^{N} (G_i \times \mathrm{LAT}_i) / \sum_{i=1}^{N} G_i

   where :math:`\mathrm{LAT}_i` is the latitude at the :math:`i`-th grid point within the front zone. 
   Obviously, this definition reflects a weighted-average of :math:`\mathrm{LAT}_i` with respect to :math:`G_i`, 
   indicating that the location of a front is mainly determined by larger SST meridional gradients within the frontal zone.

   Parameters
   ----------
   data_sst_DtDy : :py:class:`xarray.DataArray<xarray.DataArray>` 
       The SST meridional gradient data.
   criteria: :py:class:`float<python.float>`, default: `0.45 *1e-5`.
       Empirically-given critical value.
   lat_range: :py:class:`list<python.list>`, default: `[24, 32]`.
       The latitude range of the oceanic frontal zone.
   lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
       The longitude range of the oceanic frontal zone.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   The location index of the STFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017), 
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts 
   and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766


.. py:function:: calc_location_SAFZ(data_sst_DtDy: xarray.DataArray, criteria=0.8 * 1e-05, lat_range=[36, 44], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray

   Calculate the location index of the subarctic frontal zone (SAFZ).
   The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

   .. math::
       \mathrm{LCT} = \sum_{i=1}^{N} (G_i \times \mathrm{LAT}_i) / \sum_{i=1}^{N} G_i

   where :math:`\mathrm{LAT}_i` is the latitude at the :math:`i`-th grid point within the front zone. 
   Obviously, this definition reflects a weighted-average of :math:`\mathrm{LAT}_i` with respect to :math:`G_i`, 
   indicating that the location of a front is mainly determined by larger SST meridional gradients within the frontal zone.

   Parameters
   ----------
   data_sst_DtDy : :py:class:`xarray.DataArray<xarray.DataArray>` 
       The SST meridional gradient data.
   criteria: :py:class:`float<python.float>`, default: `0.80 *1e-5`.
       Empirically-given critical value.
   lat_range: :py:class:`list<python.list>`, default: `[36, 44]`.
       The latitude range of the oceanic frontal zone.
   lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
       The longitude range of the oceanic frontal zone.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   The location index of the SAFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017), 
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts 
   and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766


.. py:function:: calc_location_line_STFZ(data_sst_DtDy: xarray.DataArray, criteria=0.45 * 1e-05, lat_range=[24, 32], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray

   Calculate the location of the subtropical frontal zone (STFZ). 
   The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

   .. math::
       \mathrm{LCT} = (\sum_{i=1}^{N} (G_i \times \mathrm{LAT}_i)) / G_i

   where :math:`\mathrm{LAT}_i` is the latitude at the :math:`i`-th grid point within the front zone. 
   Obviously, this definition reflects a weighted-average of :math:`\mathrm{LAT}_i` with respect to :math:`G_i`, 
   indicating that the location of a front is mainly determined by larger SST meridional gradients within the frontal zone.

   Parameters
   ----------
   data_sst_DtDy : :py:class:`xarray.DataArray<xarray.DataArray>` 
       The SST meridional gradient data.
   criteria: :py:class:`float<python.float>`, default: `0.45 *1e-5`.
       Empirically-given critical value.
   lat_range: :py:class:`list<python.list>`, default: `[24, 32]`.
       The latitude range of the oceanic frontal zone.
   lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
       The longitude range of the oceanic frontal zone.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   The location of the STFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017), 
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts 
   and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766    


.. py:function:: calc_location_line_SAFZ(data_sst_DtDy: xarray.DataArray, criteria=0.8 * 1e-05, lat_range=[36, 44], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray

   Calculate the location of the subarctic frontal zone (SAFZ).
   The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

   .. math::
       \mathrm{LCT} = (\sum_{i=1}^{N} (G_i \times \mathrm{LAT}_i)) / G_i

   where :math:`\mathrm{LAT}_i` is the latitude at the :math:`i`-th grid point within the front zone. 
   Obviously, this definition reflects a weighted-average of :math:`\mathrm{LAT}_i` with respect to :math:`G_i`, 
   indicating that the location of a front is mainly determined by larger SST meridional gradients within the frontal zone.

   Parameters
   ----------
   data_sst_DtDy : :py:class:`xarray.DataArray<xarray.DataArray>` 
       The SST meridional gradient data.
   criteria: :py:class:`float<python.float>`, default: `0.80 *1e-5`.
       Empirically-given critical value.
   lat_range: :py:class:`list<python.list>`, default: `[36, 44]`.
       The latitude range of the oceanic frontal zone.
   lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
       The longitude range of the oceanic frontal zone.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   The location of the SAFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017), 
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts 
   and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766    


.. py:function:: transfer_xarray_lon_from180TO360(data_input: xr.DataArray | xr.Dataset, lon_dim: str = 'lon') -> xr.DataArray | xr.Dataset

   Longitude conversion -180-180 to 0-360.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.

   .. seealso::
       :py:func:`transfer_xarray_lon_from360TO180 <transfer_xarray_lon_from360TO180>`


.. py:function:: generate_dataset_dispatcher(func)

   Function Dispensers: Iterate over the variables in the `xarray.Dataset` data using a function that only supports `xarray.DataArray` data


.. py:function:: interp_mesh2mesh(data_input: xr.DataArray | xr.Dataset, target_grid: xr.DataArray | xr.Dataset, lon_dim: str = 'lon', lat_dim: str = 'lat', method: str = 'linear')

   Regridding regular or lat-lon grid data.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   target_grid: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       Target grid to be regridding.

       :py:class:`xarray.DataArray<xarray.DataArray>` version sample

       .. code:: python

           target_grid = xr.DataArray(
               dims=('lat', 'lon'),
               coords={'lat': np.arange(-89, 89, 3) + 1 / 1.0, 'lon': np.arange(-180, 180, 3) + 1 / 1.0}
           )

       :py:class:`xarray.Dataset<xarray.Dataset>` version sample

       .. code:: python

           target_grid = xr.Dataset()
           target_grid['lat'] = np.arange(-89, 89, 3) + 1 / 1.0
           target_grid['lon'] = np.arange(-180, 180, 3) + 1 / 1.0

   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   method: :py:class:`str<python.str>`, default: `linear`.
       The methods of regridding.

       - `linear`: linear, bilinear, or higher dimensional linear interpolation.
       - `nearest`: nearest-neighbor regridding.
       - `cubic`: cubic spline regridding.
       - `conservative`: conservative regridding.

   Reference
   --------------
   https://github.com/EXCITED-CO2/xarray-regrid


.. py:function:: interp_point2mesh(data: pandas.DataFrame, var_name: str, point: list[int], grid_x: float, grid_y: float, resolution: float, sigma: float, lon_dim_name='lon', lat_dim_name='lat', method='optimized_convolution', num_iter=4, min_weight=0.001) -> xarray.DataArray

   Computes the Barnes interpolation for observation values `var_name` taken at sample
   points `data` using Gaussian weights for the width parameter `sigma`.
   The underlying grid embedded in a Euclidean space is given with start point
   `point`, regular x-direction length `grid_x` (degree), regular y-direction length `grid_y` (degree),
   and resolution `resolution`.

   Barnes interpolation is a method that is widely used in geospatial sciences like meteorology 
   to remodel data values recorded at irregularly distributed points into a representative 
   analytical field. It is defined as

   .. math::
       f(\boldsymbol{x})=\frac{\sum_{k=1}^N f_k\cdot w_k(\boldsymbol{x})}{\sum_{k=1}^N w_k(\boldsymbol{x})}

   with Gaussian weights

   .. math::
       w_k(\boldsymbol{x})=\text{e}^{-\frac{1}{2\sigma^2}\left\|x-\boldsymbol{x}_k\right\|^2}

   Naive computation of Barnes interpolation leads to an algorithmic complexity of :math:`O(N \times W \times H)`, 
   where :math:`N` is the number of sample points and :math:`W \times H` the size of the underlying grid.

   For sufficiently large :math:`n` (in general in the range from 3 to 6) a good approximation of 
   Barnes interpolation with a reduced complexity :math:`O(N + W \times H)` can be obtained by the convolutional expression

   .. math::
       f(\boldsymbol{x})\approx \frac{ (\sum_{k=1}^{N}f_k\cdot\delta_{\boldsymbol{x}_k}) *  ( r_n^{*n[x]}(x)\cdot r_n^{*n[y]}(y) )   }{ ( \sum_{k=1}^{N} \delta_{\boldsymbol{x}_k}  ) *  (  r_{n}^{*n[x]}(x)\cdot r_{n}^{*n[y]}(y)  )   }

   where :math:`\delta` is the Dirac impulse function and :math:`r(.)` an elementary rectangular function of a specific length that depends on :math:`\sigma` and :math:`n`.

   - data : :py:class:`pandas.DataFrame<pandas.DataFrame>`
       Longitude and latitude grid point discrete data. There should be a similar structure as follows

       +------------+------------+-----------+
       |    lon     |    lat     |    qff    |
       +============+============+===========+
       |   -3.73    |   56.33    |   995.1   |
       +------------+------------+-----------+
       |    2.64    |   47.05    |  1012.5   |
       +------------+------------+-----------+
       |    ...     |   ...      |   ...     |
       +------------+------------+-----------+

       .. note:: 
           Data points should contain longitude (`lon`), latitude (`lat`) and data variables (the above data variable name is `qff`).

   - var_name: :py:class:`str<python.str>`
       The name of the data variable. This should match the one in the parameter `data`.
   - lat_dim_name: :py:class:`str<python.str>`.
       Latitude dimension name. This should match the one in the parameter `data`. By default is `lat`.
   - lon_dim_name: :py:class:`str<python.str>`.
       Longitude dimension name. This should match the one in the parameter `data`. By default is `lon`.
   - point : numpy ndarray
       A 1-dimensional array of size 2 containing the coordinates of the
       start point of the grid to be used.
   - grid_x : :py:class:`int<python.int>`.
       Length in degrees in the x-direction of the interpolated rectangular grid.
   - grid_y : :py:class:`int<python.int>`.
       Length in degrees in the y-direction of the interpolated rectangular grid.
   - resolution: float
       Grid resolution. The distance between regular grid points is the reciprocal of the value. Common values: 4.0, 8.0, 16.0, 32.0, 64.0.
   - sigma : float
       The Gaussian width parameter to be used. Common values: 0.25, 0.5, 1.0, 2.0, 4.0.
   - method : {'optimized_convolution', 'convolution', 'radius', 'naive'}, default: 'optimized_convolution'.
       Designates the Barnes interpolation method to be used. The possible
       implementations that can be chosen are 'naive' for the straightforward
       implementation (algorithm A from paper), 'radius' to consider only sample
       points within a specific radius of influence, both with an algorithmic
       complexity of :math:`O(N \times W \times H)`.
       The choice 'convolution' implements algorithm B specified in the paper
       and 'optimized_convolution' is its optimization by appending tail values
       to the rectangular kernel. The latter two algorithms reduce the complexity
       down to :math:`O(N + W \times H)`.
       The default is 'optimized_convolution'.
   - num_iter : int, optional
       The number of performed self-convolutions of the underlying rect-kernel.
       Applies only if method is 'optimized_convolution' or 'convolution'.
       The default is 4. Applies only to Convol interpolations: one of 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50.
   - min_weight : float, optional
       Choose radius of influence such that Gaussian weight of considered sample
       points is greater than `min_weight`.
       Applies only if method is 'radius'. Recommended values are 0.001 and less.
       The default is 0.001, which corresponds to a radius of 3.717 * sigma.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::   
       - https://github.com/MeteoSwiss/fast-barnes-py
       - Zürcher, B. K.: Fast approximate Barnes interpolation: illustrated by Python-Numba implementation fast-barnes-py v1.0, Geosci. Model Dev., 16, 1697–1711, https://doi.org/10.5194/gmd-16-1697-2023, 2023.


.. py:function:: interp_point2mesh_S2(data: pandas.DataFrame, var_name: str, point: list[int], grid_x: float, grid_y: float, resolution: float, sigma: float, lon_dim_name='lon', lat_dim_name='lat', method='optimized_convolution_S2', num_iter=4, resample=True) -> xarray.DataArray

   Computes the Barnes interpolation for observation values `var_name` taken at sample
   points `data` using Gaussian weights for the width parameter `sigma`.

   The underlying grid embedded on the unit sphere :math:`S^2` and thus inherits the
   spherical distance measure (taken in degrees). The grid is given by the start
   point `point`, regular x-direction length `grid_x` (degree), regular y-direction length `grid_y` (degree),
   and resolution `resolution`.

   Parameters
   ----------
   - data : :py:class:`pandas.DataFrame<pandas.DataFrame>`
       Longitude and latitude grid point discrete data. There should be a similar structure as follows

       +------------+------------+-----------+
       |    lon     |    lat     |    qff    |
       +============+============+===========+
       |   -3.73    |   56.33    |   995.1   |
       +------------+------------+-----------+
       |    2.64    |   47.05    |  1012.5   |
       +------------+------------+-----------+
       |    ...     |   ...      |   ...     |
       +------------+------------+-----------+

       .. note:: 
           Data points should contain longitude (`lon`), latitude (`lat`) and data variables (the above data variable name is `qff`).

   - var_name: :py:class:`str<python.str>`
       The name of the data variable. This should match the one in the parameter `data`.
   - lat_dim_name: :py:class:`str<python.str>`.
       Latitude dimension name. This should match the one in the parameter `data`. By default is `lat`.
   - lon_dim_name: :py:class:`str<python.str>`.
       Longitude dimension name. This should match the one in the parameter `data`. By default is `lon`.
   - point : numpy ndarray
       A 1-dimensional array of size 2 containing the coordinates of the
       start point of the grid to be used.
   - grid_x : :py:class:`int<python.int>`.
       Length in degrees in the x-direction of the interpolated rectangular grid.
   - grid_y : :py:class:`int<python.int>`.
       Length in degrees in the y-direction of the interpolated rectangular grid.
   - resolution: float
       Grid resolution. The distance between regular grid points is the reciprocal of the value. Common values: 4.0, 8.0, 16.0, 32.0, 64.0.
   - sigma : float
       The Gaussian width parameter to be used. Common values: 0.25, 0.5, 1.0, 2.0, 4.0.
   - method : {'optimized_convolution_S2', 'naive_S2'}, default: 'optimized_convolution_S2'.
       Designates the Barnes interpolation method to be used. The possible
       implementations that can be chosen are 'naive_S2' for the straightforward
       implementation (algorithm A from the paper) with an algorithmic complexity
       of :math:`O(N \times W \times H)`.
       The choice 'optimized_convolution_S2' implements the optimized algorithm B
       specified in the paper by appending tail values to the rectangular kernel.
       The latter algorithm has a reduced complexity of :math:`O(N + W \times H)`.
       The default is 'optimized_convolution_S2'.
   - num_iter : int, optional, default: 4.
       The number of performed self-convolutions of the underlying rect-kernel.
       Applies only if method is 'optimized_convolution_S2'.
       The default is 4. Applies only to Convol interpolations: one of 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50.
   - resample : bool, optional, default: `True`.
       Specifies whether to resample Lambert grid field to lonlat grid.
       Applies only if method is 'optimized_convolution_S2'.
       The default is True.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::   
       - https://github.com/MeteoSwiss/fast-barnes-py
       - Zürcher, B. K.: Fast approximate Barnes interpolation: illustrated by Python-Numba implementation fast-barnes-py v1.0, Geosci. Model Dev., 16, 1697–1711, https://doi.org/10.5194/gmd-16-1697-2023, 2023.


.. py:function:: field_grids(data, grids)

   For each point of the grids, its nearby region are chosen.

   Parameters
   ----------
   data : array
       An N-dimensional array.
   grids : int or 2-size tuple
       The total number of grid points of the x and y directions.

   Returns
   -------
   D : ndarray
       A ndarray with extra two dimensions. The last two dimensions are
       each points's nearby region.


.. py:class:: BarnesFilter(data_arr, lon=None, lat=None, radius_degree=10)

   The Barnes method performs grid point interpolation by selecting appropriate
   filtering parameters *c* and *g* to filter out shortwave noise in the original field,
   making the analysis results stable and smooth. In addition, it can form a bandpass filter
   to separate various sub weather scales that affect weather processes according to actual needs,
   achieving the purpose of scale separation.

   Reference:
   DOI : https://doi.org/10.1175/1520-0493(1980)108<1108:AOTFSM>2.0.CO;2

   .. py:method:: __convert_data(data)


   .. py:method:: __calculate_distance(lon, lat)


   .. py:method:: __lowpass(g=0.3, c=150000)


   .. py:method:: lowpass(g=0.3, c=150000)

      Selecting different parameters *g* and *c*
      will result in different filtering characteristics.

      Reference:
      DOI : https://doi.org/10.1175/1520-0493(1980)108<1108:AOTFSM>2.0.CO;2

      Parameters
      ----------
      g : float, generally between (0, 1]
          Constant parameter.
      c : int
          Constant parameter. When *c* takes a larger value, the filter function converges
          at a larger wavelength, and the response function slowly approaches the maximum value,
          which means that high-frequency fluctuations have been filtered out.

      Returns
      -------
      data_vars : array
          Data field after filtering out high-frequency fluctuations



   .. py:method:: bandpass(g1=0.3, c1=30000, g2=0.3, c2=150000)

      Select two different filtering schemes 1 and 2, and perform the filtering separately.
      And then perform the difference, that means *scheme1 - scheme2*.
      The mesoscale fluctuations are thus preserved.

      Parameters
      ----------
      g1 : float, generally between (0, 1]
          Constant parameter of scheme1.
      c1 : int
          Constant parameterof scheme1.
      g2 : float, generally between (0, 1]
          Constant parameter of scheme2.
      c2 : int
          Constant parameterof scheme2.

      Returns
      -------
      data_vars : array
          Mesoscale wave field filtered out from raw data



.. py:function:: find_dims_axis(data: xarray.DataArray, dim: str) -> int

   Find the index of `dim` in the xarray DataArray.

   Parameters
   ----------
   - data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   - dim : :py:class:`str<python.str>`
       Dimension(s) over which to find axis.

   Returns
   -------
   :py:class:`int<python.int>`.


.. py:function:: calc_butter_bandpass(data: xarray.DataArray, sampling_frequency: int, period: list, N=3, dim='time') -> xarray.DataArray

   Butterworth bandpass filter.

   Parameters
   ----------
   - data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The array of data to be filtered.
   - sampling_frequency: int.
       Data sampling frequency. If it is daily data with only one time level record per day, 
       then the parameter value is 1; If a day contains four time level records, the parameter value is 4.
   - period: list.
       The time period interval of the bandpass filter to be acquired. 
       If we are obtaining a 3-10 day bandpass filter, the value of this parameter is `[3, 10]`. 
       Note that the units of this parameter should be consistent with `sampling_frequency`.
   - N: int.
       The order of the filter. Default is 3.
   - dim: str.
       Dimension(s) over which to apply bandpass filter. By default gradient is applied over the `time` dimension.

   .. seealso::
       :py:func:`scipy.signal.butter <scipy:scipy.signal.butter>`, :py:func:`scipy.signal.sosfiltfilt <scipy:scipy.signal.sosfiltfilt>`


.. py:function:: calc_butter_lowpass(data: xarray.DataArray, sampling_frequency: int, period: int, N=3, dim='time') -> xarray.DataArray

   Butterworth lowpass filter.

   Parameters
   ----------
   - data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The array of data to be filtered.
   - sampling_frequency: int.
       Data sampling frequency. If it is daily data with only one time level record per day, 
       then the parameter value is 1; If a day contains four time level records, the parameter value is 4.
   - period: list.
       The low-pass filtering time period, above which the signal (low frequency signal) will pass. 
       If you are getting a 10-day low-pass filter, the value of this parameter is `10`. 
       Note that the units of this parameter should be consistent with `sampling_frequency`.
   - N: int.
       The order of the filter. Default is 3.
   - dim: str.
       Dimension(s) over which to apply lowpass filter. By default gradient is applied over the `time` dimension.

   .. seealso::
       :py:func:`scipy.signal.butter <scipy:scipy.signal.butter>`, :py:func:`scipy.signal.sosfiltfilt <scipy:scipy.signal.sosfiltfilt>`


.. py:function:: calc_butter_highpass(data: xarray.DataArray, sampling_frequency: int, period: int, N=3, dim='time') -> xarray.DataArray

   Butterworth highpass filter.

   Parameters
   ----------
   - data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The array of data to be filtered.
   - sampling_frequency: int.
       Data sampling frequency. If it is daily data with only one time level record per day, 
       then the parameter value is 1; If a day contains four time level records, the parameter value is 4.
   - period: list.
       The high-pass filtering time period below which the signal (high-frequency signal) will pass. 
       If you are obtaining a 10-day high-pass filter, the value of this parameter is `10`. 
       Note that the units of this parameter should be consistent with `sampling_frequency`.
   - N: int.
       The order of the filter. Default is 3.
   - dim: str.
       Dimension(s) over which to apply highpass filter. By default gradient is applied over the `time` dimension.

   .. seealso::
       :py:func:`scipy.signal.butter <scipy:scipy.signal.butter>`, :py:func:`scipy.signal.sosfiltfilt <scipy:scipy.signal.sosfiltfilt>`


