easyclimate.field.ocean
=======================

.. py:module:: easyclimate.field.ocean


Submodules
----------

.. toctree::
   :maxdepth: 1

   /technical/api/easyclimate/field/ocean/mixlayer/index
   /technical/api/easyclimate/field/ocean/oceanic_front/index
   /technical/api/easyclimate/field/ocean/stability/index
   /technical/api/easyclimate/field/ocean/thermal/index


Functions
---------

.. autoapisummary::

   easyclimate.field.ocean.calc_gradient
   easyclimate.field.ocean.calc_u_advection
   easyclimate.field.ocean.calc_v_advection
   easyclimate.field.ocean.calc_mixed_layer_depth
   easyclimate.field.ocean.calc_MLD_depth_weighted
   easyclimate.field.ocean.calc_MLD_temper_tendency
   easyclimate.field.ocean.get_data_within_MLD
   easyclimate.field.ocean.get_temper_within_MLD
   easyclimate.field.ocean.get_data_average_within_MLD
   easyclimate.field.ocean.get_temper_average_within_MLD
   easyclimate.field.ocean.calc_MLD_average_horizontal_advection
   easyclimate.field.ocean.calc_MLD_average_vertical_advection
   easyclimate.field.ocean.calc_ocean_surface_heat_flux
   easyclimate.field.ocean.get_weighted_spatial_data
   easyclimate.field.ocean.sort_ascending_latlon_coordinates
   easyclimate.field.ocean.calc_intensity_STFZ
   easyclimate.field.ocean.calc_intensity_SAFZ
   easyclimate.field.ocean.calc_location_STFZ
   easyclimate.field.ocean.calc_location_SAFZ
   easyclimate.field.ocean.calc_location_line_STFZ
   easyclimate.field.ocean.calc_location_line_SAFZ
   easyclimate.field.ocean.find_dims_axis
   easyclimate.field.ocean.calc_N2_from_temp_salt
   easyclimate.field.ocean.calc_potential_density_from_temp_salt
   easyclimate.field.ocean.calc_gradient
   easyclimate.field.ocean.calc_seawater_thermocline_depth
   easyclimate.field.ocean.calc_Dx_depth
   easyclimate.field.ocean.calc_D14_depth
   easyclimate.field.ocean.calc_D17_depth
   easyclimate.field.ocean.calc_D20_depth
   easyclimate.field.ocean.calc_D26_depth
   easyclimate.field.ocean.calc_D28_depth


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
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.

   Returns
   -------
   The gradient along the coordinate `dim` direction (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`numpy.gradient <numpy:numpy.gradient>`


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


.. py:function:: calc_mixed_layer_depth(seawater_temperature_data: xarray.DataArray, seawater_practical_salinity_data: xarray.DataArray, criterion: {'temperature', 'density', 'pdvar'} = 'pdvar', depth_dim: str = 'depth', lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.DataArray

   Calculate the mixed layer depth.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{^\circ C}`)
       Vertical seawater temperature data.
   seawater_practical_salinity_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{PSU}`)
       Vertical seawater salinity data (practical salinity).
   criterion: {'temperature', 'density', 'pdvar'}, default `'pdvar'`.
       Mixed layer depth criteria

       - **temperature** : Computed based on constant temperature difference criterion, i.e., :math:`CT(0) - T[mld] = 0.5 \mathrm{^\circ C}`.
       - **density** : Computed based on the constant potential density difference criterion, i.e., :math:`pd[0] - pd[mld] = 0.125` in sigma units.
       - **pdvar** : Computed based on variable potential density criterion :math:`pd[0] - pd[mld] = var(T[0], S[0])`, where var is a variable potential density difference which corresponds to constant temperature difference of :math:`0.5 \mathrm{^\circ C}`.

   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The mixed layer depth (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::

       - https://github.com/pyoceans/oceans
       - https://pyoceans.github.io/python-oceans/ocfis.html#oceans.ocfis.mld


.. py:function:: calc_MLD_depth_weighted(seawater_temperature_data: xarray.DataArray | xarray.Dataset, mixed_layer_depth: xarray.DataArray, depth_dim: str = 'depth') -> xarray.DataArray | xarray.Dataset

   Calculate the weights of mixed layer depth. The weights required by the ocean model under non-uniform distribution grids in the depth direction.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{^\circ C}`).
       Vertical seawater temperature data.
   mixed_layer_depth: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{m}`).
       The mixed layer depth.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   The weights of the mixed layer depth (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).


.. py:function:: calc_MLD_temper_tendency(seawater_temperature_anomaly_data: xarray.DataArray, mixed_layer_depth: xarray.DataArray, depth_weight: xarray.DataArray, depth_dim='depth', time_dim='time') -> xarray.DataArray

   Calculate the tendency of the mixing layer temperature.

   Parameters
   ----------
   seawater_temperature_anomaly_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{^\circ C}`).
       The anomaly of the vertical seawater temperature data.
   mixed_layer_depth: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{m}`).
       The mixed layer depth.
   depth_weight: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The weights of the mixed layer depth. The data is generated by the :py:class:`easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted <easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted>`
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name.

   Returns
   -------
   The weights of the mixed layer depth (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: get_data_within_MLD(data_input: xarray.DataArray, mixed_layer_depth: xarray.DataArray, depth_dim: str = 'depth') -> xarray.DataArray

   Obtain data within the mixed layer.

   .. caution::

       This function sets the data outside the mixing layer as missing values, i.e. `np.nan`,
       but it does not calculate the average value for the data inside the mixing layer.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The spatio-temporal data to be calculated.
   mixed_layer_depth: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{m}`).
       The mixed layer depth.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   The data within the mixed layer (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: get_temper_within_MLD(seawater_temperature_data: xarray.DataArray, mixed_layer_depth: xarray.DataArray, depth_dim: str = 'depth') -> xarray.DataArray

   Obtain seawater temperature data within the mixing layer.

   .. caution::

       This function sets the data outside the mixing layer as missing values, i.e. `np.nan`,
       but it does not calculate the average value for the data inside the mixing layer.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{^\circ C}`)
       Vertical seawater temperature data.
   mixed_layer_depth: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{m}`).
       The mixed layer depth.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   The seawater temperature data within the mixed layer (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: get_data_average_within_MLD(data_input: xarray.DataArray, mixed_layer_depth: xarray.DataArray, depth_weight: xarray.DataArray, depth_dim: str = 'depth') -> xarray.DataArray

   Obtain averaged data within the mixed layer.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`
        The spatio-temporal data to be calculated.
   mixed_layer_depth: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{m}`).
       The mixed layer depth.
   depth_weight: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The weights of the mixed layer depth. The data is generated by the :py:class:`easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted <easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted>`
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   The averaged data within the mixed layer (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: get_temper_average_within_MLD(seawater_temperature_data: xarray.DataArray, mixed_layer_depth: xarray.DataArray, depth_weight: xarray.DataArray, depth_dim: str = 'depth') -> xarray.DataArray

   Obtain averaged seawater temperature data within the mixing layer.

   .. caution::

       This function sets the data outside the mixing layer as missing values, i.e. `np.nan`,
       but it does not calculate the average value for the data inside the mixing layer.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{^\circ C}`)
       Vertical seawater temperature data.
   mixed_layer_depth: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{m}`).
       The mixed layer depth.
   depth_weight: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The weights of the mixed layer depth. The data is generated by the :py:class:`easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted <easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted>`
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   The averaged seawater temperature data within the mixed layer (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_MLD_average_horizontal_advection(u_monthly_data: xarray.DataArray, v_monthly_data: xarray.DataArray, seawater_temperature_data: xarray.DataArray, mixed_layer_depth: xarray.DataArray, depth_weight: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', depth_dim: str = 'depth', min_dx: float = 1.0, min_dy: float = 1.0, edge_order: int = 2, R: float = 6370000) -> xarray.DataArray

   Obtain the average horizontal advection within the mixed layer

   Parameters
   ----------
   u_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{m \cdot s^{-1}}`).
       The monthly ocean current data.
   v_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{m \cdot s^{-1}}`).
       The monthly meridional ocean current data.
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{^\circ C}`).
       Vertical seawater temperature data.
   mixed_layer_depth: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{m}`).
       The mixed layer depth.
   depth_weight: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The weights of the mixed layer depth. The data is generated by the :py:class:`easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted <easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted>`
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.
   min_dx: :py:class:`float <float>`, default: `1.0`.
       The minimum acceptable value of `dx`, below which parts will set `nan` to avoid large computational errors.
       The unit is m. You can set it to a negative value in order to remove this benefit.
   min_dy: :py:class:`float <float>`, default: `1.0`.
       The minimum acceptable value of `dy`, below which parts will set `nan` to avoid large computational errors.
       The unit is m. You can set it to a negative value in order to remove this benefit.
   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The average horizontal advection within the mixed layer (:math:`\mathrm{^\circ C} \cdot \mathrm{month}^{-1}`, :py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_MLD_average_vertical_advection(w_monthly_data: xarray.DataArray, seawater_temperature_data: xarray.DataArray, mixed_layer_depth: xarray.DataArray, depth_weight: xarray.DataArray, depth_dim: str = 'depth') -> xarray.DataArray

   Obtain the average vertical advection within the mixed layer.

   Parameters
   ----------
   w_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{m \cdot s^{-1}}`).
       The monthly vertical ocean current data.
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{^\circ C}`).
       Vertical seawater temperature data.
   mixed_layer_depth: :py:class:`xarray.DataArray<xarray.DataArray>`(:math:`\mathrm{m}`).
       The mixed layer depth.
   depth_weight: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The weights of the mixed layer depth. The data is generated by the :py:class:`easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted <easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted>`
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   The average vertical advection within the mixed layer (:math:`\mathrm{^\circ C} \cdot \mathrm{month}^{-1}`, :py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_ocean_surface_heat_flux(qnet_monthly_anomaly_data: xarray.DataArray, mixed_layer_depth: xarray.DataArray | float, rho_0: float = 1000, c_p: float = 4000) -> xarray.DataArray

   Obtain ocean surface heat flux.

   Parameters
   ----------
   qnet_monthly_anomaly_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{W \cdot m^{-2}}`).
       The monthly anomaly of the downward net heat flux at the ocean surface.
   mixed_layer_depth: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{m}`).
       The mixed layer depth.
   rho_0: :py:class:`float <float>`, default: `1000` (:math:`\mathrm{kg \cdot m^{-3}}`).
       The density of water.
   c_p: :py:class:`float <float>`, default: `4000` (:math:`\mathrm{J \cdot kg \cdot K^{-1}}`).
       The specific heat of water.

   Returns
   -------
   The ocean surface heat flux (:math:`\mathrm{^\circ C} \cdot \mathrm{month}^{-1}`, :py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Nnamchi, H., Li, J., Kucharski, F. et al. Thermodynamic controls of the Atlantic Niño. Nat Commun 6, 8895 (2015). https://doi.org/10.1038/ncomms9895


.. py:function:: get_weighted_spatial_data(data_input: xarray.DataArray, lat_dim: str = 'lat', lon_dim: str = 'lon', method: str = 'cos_lat') -> xarray.DataArray

   Get the area-weighting data.

   Parameters
   ----------
   - data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   - lat_dim: :py:class:`str <str>`.
       Latitude dimension over which to apply. By default is applied over the `lat` dimension.
   - lon_dim: :py:class:`str <str>`.
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


.. py:function:: sort_ascending_latlon_coordinates(data: xarray.DataArray | xarray.Dataset, lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray | xarray.Dataset

   Sort the dimensions `lat`, `lon` in ascending order.


.. py:function:: calc_intensity_STFZ(sst_DtDy_data: xarray.DataArray, criteria=0.45 * 1e-05, lat_range=[24, 32], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray

   Calculate the intensity of the subtropical frontal zone (STFZ).
   The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

   .. math::
       \mathrm{ITS} = \sum_{i=1}^{N} \frac{G_i}{N}

   where :math:`G_i` is the value of zonally-averaged SST meridional gradient that is no less than
   an empirically-given critical value (here, :math:`0.45 \times 10^{-5} \mathrm{km^{-1}}` for STFZ) at the :math:`i`-th latitudinal grid point within the zone,
   and :math:`N` is the number of total grid points that satisfy the criteria above.

   Parameters
   ----------
   sst_DtDy_data : :py:class:`xarray.DataArray<xarray.DataArray>`
       The SST meridional gradient data.
   criteria: :py:class:`float <float>`, default: `0.45 *1e-5`.
       Empirically-given critical value.
   lat_range: :py:class:`list<python.list>`, default: `[24, 32]`.
       The latitude range of the oceanic frontal zone.
   lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
       The longitude range of the oceanic frontal zone.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   The intensity of the STFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017),
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts
   and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766


.. py:function:: calc_intensity_SAFZ(sst_DtDy_data: xarray.DataArray, criteria=0.8 * 1e-05, lat_range=[36, 44], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray

   Calculate the intensity of the subarctic frontal zone (SAFZ).
   The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

   .. math::
       \mathrm{ITS} = \sum_{i=1}^{N} \frac{G_i}{N}

   where :math:`G_i` is the value of zonally-averaged SST meridional gradient that is no less than
   an empirically-given critical value (here, :math:`0.80 \times 10^{-5} \mathrm{km^{-1}}` for SAFZ) at the :math:`i`-th latitudinal grid point within the zone,
   and :math:`N` is the number of total grid points that satisfy the criteria above.

   Parameters
   ----------
   sst_DtDy_data : :py:class:`xarray.DataArray<xarray.DataArray>`
       The SST meridional gradient data.
   criteria: :py:class:`float <float>`, default: `0.80 *1e-5`.
       Empirically-given critical value.
   lat_range: :py:class:`list<python.list>`, default: `[36, 44]`.
       The latitude range of the oceanic frontal zone.
   lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
       The longitude range of the oceanic frontal zone.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   The intensity of the SAFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017),
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts
   and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766


.. py:function:: calc_location_STFZ(sst_DtDy_data: xarray.DataArray, criteria=0.45 * 1e-05, lat_range=[24, 32], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray

   Calculate the location index of the subtropical frontal zone (STFZ).
   The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

   .. math::
       \mathrm{LCT} = \sum_{i=1}^{N} (G_i \times \mathrm{LAT}_i) / \sum_{i=1}^{N} G_i

   where :math:`\mathrm{LAT}_i` is the latitude at the :math:`i`-th grid point within the front zone.
   Obviously, this definition reflects a weighted-average of :math:`\mathrm{LAT}_i` with respect to :math:`G_i`,
   indicating that the location of a front is mainly determined by larger SST meridional gradients within the frontal zone.

   Parameters
   ----------
   sst_DtDy_data : :py:class:`xarray.DataArray<xarray.DataArray>`
       The SST meridional gradient data.
   criteria: :py:class:`float <float>`, default: `0.45 *1e-5`.
       Empirically-given critical value.
   lat_range: :py:class:`list<python.list>`, default: `[24, 32]`.
       The latitude range of the oceanic frontal zone.
   lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
       The longitude range of the oceanic frontal zone.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   The location index of the STFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017),
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts
   and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766


.. py:function:: calc_location_SAFZ(sst_DtDy_data: xarray.DataArray, criteria=0.8 * 1e-05, lat_range=[36, 44], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray

   Calculate the location index of the subarctic frontal zone (SAFZ).
   The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

   .. math::
       \mathrm{LCT} = \sum_{i=1}^{N} (G_i \times \mathrm{LAT}_i) / \sum_{i=1}^{N} G_i

   where :math:`\mathrm{LAT}_i` is the latitude at the :math:`i`-th grid point within the front zone.
   Obviously, this definition reflects a weighted-average of :math:`\mathrm{LAT}_i` with respect to :math:`G_i`,
   indicating that the location of a front is mainly determined by larger SST meridional gradients within the frontal zone.

   Parameters
   ----------
   sst_DtDy_data : :py:class:`xarray.DataArray<xarray.DataArray>`
       The SST meridional gradient data.
   criteria: :py:class:`float <float>`, default: `0.80 *1e-5`.
       Empirically-given critical value.
   lat_range: :py:class:`list<python.list>`, default: `[36, 44]`.
       The latitude range of the oceanic frontal zone.
   lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
       The longitude range of the oceanic frontal zone.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   The location index of the SAFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017),
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts
   and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766


.. py:function:: calc_location_line_STFZ(sst_DtDy_data: xarray.DataArray, criteria=0.45 * 1e-05, lat_range=[24, 32], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray

   Calculate the location of the subtropical frontal zone (STFZ).
   The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

   .. math::
       \mathrm{LCT} = (\sum_{i=1}^{N} (G_i \times \mathrm{LAT}_i)) / G_i

   where :math:`\mathrm{LAT}_i` is the latitude at the :math:`i`-th grid point within the front zone.
   Obviously, this definition reflects a weighted-average of :math:`\mathrm{LAT}_i` with respect to :math:`G_i`,
   indicating that the location of a front is mainly determined by larger SST meridional gradients within the frontal zone.

   Parameters
   ----------
   sst_DtDy_data : :py:class:`xarray.DataArray<xarray.DataArray>`
       The SST meridional gradient data.
   criteria: :py:class:`float <float>`, default: `0.45 *1e-5`.
       Empirically-given critical value.
   lat_range: :py:class:`list<python.list>`, default: `[24, 32]`.
       The latitude range of the oceanic frontal zone.
   lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
       The longitude range of the oceanic frontal zone.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   The location of the STFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017),
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts
   and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766


.. py:function:: calc_location_line_SAFZ(sst_DtDy_data: xarray.DataArray, criteria=0.8 * 1e-05, lat_range=[36, 44], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray

   Calculate the location of the subarctic frontal zone (SAFZ).
   The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

   .. math::
       \mathrm{LCT} = (\sum_{i=1}^{N} (G_i \times \mathrm{LAT}_i)) / G_i

   where :math:`\mathrm{LAT}_i` is the latitude at the :math:`i`-th grid point within the front zone.
   Obviously, this definition reflects a weighted-average of :math:`\mathrm{LAT}_i` with respect to :math:`G_i`,
   indicating that the location of a front is mainly determined by larger SST meridional gradients within the frontal zone.

   Parameters
   ----------
   sst_DtDy_data : :py:class:`xarray.DataArray<xarray.DataArray>`
       The SST meridional gradient data.
   criteria: :py:class:`float <float>`, default: `0.80 *1e-5`.
       Empirically-given critical value.
   lat_range: :py:class:`list<python.list>`, default: `[36, 44]`.
       The latitude range of the oceanic frontal zone.
   lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
       The longitude range of the oceanic frontal zone.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   The location of the SAFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017),
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts
   and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766


.. py:function:: find_dims_axis(data: xarray.DataArray, dim: str) -> int

   Find the index of `dim` in the xarray DataArray.

   Parameters
   ----------
   - data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   - dim : :py:class:`str <str>`
       Dimension(s) over which to find axis.

   Returns
   -------
   :py:class:`int <int>`.


.. py:function:: calc_N2_from_temp_salt(seawater_temperature_data: xarray.DataArray, seawater_practical_salinity_data: xarray.DataArray, time_dim: str | None, depth_dim: str = 'depth', lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.Dataset

   Calculate the frequency of seawater buoyancy.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{^\circ C}`)
       Vertical seawater temperature data.
   seawater_practical_salinity_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{PSU}`)
       Vertical seawater salinity data (practical salinity).
   time_dim: :py:class:`str <str>` or `None`, default: `time`.
       The time coordinate dimension name.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The frequency of seawater buoyancy (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::

       - http://www.teos-10.org/pubs/gsw/html/gsw_Nsquared.html
       - http://www.teos-10.org/pubs/gsw/html/gsw_contents.html


.. py:function:: calc_potential_density_from_temp_salt(seawater_temperature_data: xarray.DataArray, seawater_practical_salinity_data: xarray.DataArray, time_dim: str | None, depth_dim: str = 'depth', lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.Dataset

   Calculate the potential density of seawater.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{^\circ C}`)
       Vertical seawater temperature data.
   seawater_practical_salinity_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{PSU}`)
       Vertical seawater salinity data (practical salinity).
   time_dim: :py:class:`str <str>` or `None`, default: `time`.
       The time coordinate dimension name.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The potential density of seawater (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::

       - http://www.teos-10.org/pubs/gsw/html/gsw_Nsquared.html
       - http://www.teos-10.org/pubs/gsw/html/gsw_contents.html


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


.. py:function:: calc_seawater_thermocline_depth(seawater_temperature_data: xarray.DataArray, depth_dim: str = 'depth') -> xarray.DataArray

   Caculate thermocline depth of ocean temperature.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\mathrm{^\circ C}`).
       ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: calc_Dx_depth(seawater_temperature_data: xarray.DataArray, value: float, depth_dim: str = 'depth') -> xarray.DataArray

   Caculate `value` depth of ocean temperature.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\mathrm{^\circ C}`).
       ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
   value: float.
       The depth of ocean temperature to be calculated.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: calc_D14_depth(seawater_temperature_data: xarray.DataArray, value: float = 14, depth_dim: str = 'depth') -> xarray.DataArray

   Caculate 14m depth of ocean temperature.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\mathrm{^\circ C}`).
       ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
   value: float (:math:`\mathrm{m}`).
       The depth of ocean temperature to be calculated.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: calc_D17_depth(seawater_temperature_data: xarray.DataArray, value: float = 17, depth_dim: str = 'depth') -> xarray.DataArray

   Caculate 17m depth of ocean temperature.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\mathrm{^\circ C}`).
       ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
   value: float (:math:`\mathrm{m}`).
       The depth of ocean temperature to be calculated.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: calc_D20_depth(seawater_temperature_data: xarray.DataArray, value: float = 20, depth_dim: str = 'depth') -> xarray.DataArray

   Caculate 20m depth of ocean temperature.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\mathrm{^\circ C}`).
       ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
   value: float (:math:`\mathrm{m}`).
       The depth of ocean temperature to be calculated.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: calc_D26_depth(seawater_temperature_data: xarray.DataArray, value: float = 26, depth_dim: str = 'depth') -> xarray.DataArray

   Caculate 26m depth of ocean temperature.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\mathrm{^\circ C}`).
       ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
   value: float (:math:`\mathrm{m}`).
       The depth of ocean temperature to be calculated.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: calc_D28_depth(seawater_temperature_data: xarray.DataArray, value: float = 28, depth_dim: str = 'depth') -> xarray.DataArray

   Caculate 28m depth of ocean temperature.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\mathrm{^\circ C}`).
       ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
   value: float (:math:`\mathrm{m}`).
       The depth of ocean temperature to be calculated.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


