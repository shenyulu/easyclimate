easyclimate
===========

.. py:module:: easyclimate


Subpackages
-----------

.. toctree::
   :maxdepth: 1

   /technical/api/easyclimate/core/index
   /technical/api/easyclimate/field/index
   /technical/api/easyclimate/filter/index
   /technical/api/easyclimate/interp/index
   /technical/api/easyclimate/plot/index


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

   easyclimate.calc_brunt_vaisala_frequency_atm
   easyclimate.get_coriolis_parameter
   easyclimate.calc_potential_temperature
   easyclimate.calc_virtual_temperature_Hobbs2006
   easyclimate.calc_virtual_temperature
   easyclimate.calc_static_stability
   easyclimate.calc_dewpoint
   easyclimate.calc_mixing_ratio
   easyclimate.calc_vapor_pressure
   easyclimate.calc_saturation_vapor_pressure
   easyclimate.calc_saturation_mixing_ratio
   easyclimate.transfer_mixing_ratio_2_specific_humidity
   easyclimate.transfer_specific_humidity_2_mixing_ratio
   easyclimate.transfer_dewpoint_2_specific_humidity
   easyclimate.transfer_specific_humidity_2_dewpoint
   easyclimate.transfer_dewpoint_2_relative_humidity
   easyclimate.transfer_mixing_ratio_2_relative_humidity
   easyclimate.transfer_specific_humidity_2_relative_humidity
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
   easyclimate.calc_ttestSpatialPattern_spatial
   easyclimate.calc_levenetestSpatialPattern_spatial
   easyclimate.calc_skewness_spatial
   easyclimate.calc_kurtosis_spatial
   easyclimate.calc_theilslopes_spatial
   easyclimate.calc_all_climatological_mean
   easyclimate.calc_seasonal_climatological_mean
   easyclimate.calc_seasonal_cycle_mean
   easyclimate.calc_seasonal_cycle_std
   easyclimate.calc_seasonal_cycle_var
   easyclimate.calc_seasonal_mean
   easyclimate.remove_seasonal_cycle_mean
   easyclimate.calc_monthly_climatological_std_without_seasonal_cycle_mean
   easyclimate.calc_monthly_climatological_var_without_seasonal_cycle_mean
   easyclimate.calc_horizontal_wind_components_std
   easyclimate.populate_monmean2everymon
   easyclimate.populate_daymean2everyday
   easyclimate.calc_daily_climatological_anomaly
   easyclimate.calc_yearly_climatological_mean
   easyclimate.calc_yearly_climatological_sum
   easyclimate.calc_yearly_climatological_std
   easyclimate.calc_yearly_climatological_var
   easyclimate.calc_yearly_climatological_max
   easyclimate.calc_yearly_climatological_min
   easyclimate.open_tutorial_dataset


Package Contents
----------------

.. py:class:: DataNode(name='root')

   .. py:attribute:: _attributes


   .. py:attribute:: name


   .. py:method:: __getattr__(key)


   .. py:method:: __setattr__(key, value)


   .. py:method:: __getitem__(key)


   .. py:method:: __setitem__(key, value)


   .. py:method:: format_tree(level=0, html=False)


   .. py:method:: __repr__()


.. py:function:: calc_brunt_vaisala_frequency_atm(potential_temperature_data: xarray.DataArray, z_data: xarray.DataArray, vertical_dim: str, g: float = 9.8) -> xarray.DataArray

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

   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   g: :py:class:`float <float>`, default: `9.8`.
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



.. py:function:: get_coriolis_parameter(lat_data: xarray.DataArray | numpy.array, omega: float = 7.292e-05) -> xarray.DataArray | numpy.array

   Calculate the Coriolis parameter at each point.

   .. math::
       f = 2 \Omega \sin(\phi)

   Parameters
   ----------
   lat_data: :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`numpy.array <numpy.array>`.
       Latitude at each point.
   omega: :py:class:`float <float>`, default: `7.292e-5`.
       The angular speed of the earth.

   Returns
   -------
   Corresponding Coriolis force at each point (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`numpy.array <numpy.array>`).

   Reference
   --------------
   - `Coriolis parameter - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Coriolis_parameter>`__

   .. seealso::
       - `coriolis_parameter — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.coriolis_parameter.html>`__
       - `coriolis_param - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/coriolis_param.shtml>`__


.. py:function:: calc_potential_temperature(temper_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: str, kappa: float = 287 / 1005.7) -> xarray.DataArray

   Calculate the potential temperature.

   Uses the Poisson equation to calculation the potential temperature given pressure and temperature.

   .. math::
       \theta = T \left( \frac{p_0}{p} \right) ^\kappa

   Parameters
   ----------
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   kappa: :py:class:`float <float>`, default: `287/1005.7`.
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


.. py:function:: calc_virtual_temperature_Hobbs2006(temper_data: xarray.DataArray, specific_humidity_data: xarray.DataArray, specific_humidity_data_units: Literal['kg/kg', 'g/g', 'g/kg'], epsilon: float = 0.6219569100577033) -> xarray.DataArray

   Calculate virtual temperature.

   The virtual temperature (:math:`T_v`) is the temperature at which dry air would have the same density as the moist air, at a given pressure.
   In other words, two air samples with the same virtual temperature have the same density, regardless of their actual temperature or relative humidity.
   The virtual temperature is always greater than  the absolute air temperature.

   This calculation must be given an air parcel's temperature and mixing ratio. The implementation uses the formula outlined in [Hobbs2006] pg.67 & 80.

   .. math::
       T_v = T \frac{\text{q} + \epsilon}{\epsilon\,(1 + \text{q})}

   where :math:`\epsilon \approx 0.622` when the mixing ratio (specific humidity) :math:`q` is expressed in :math:`\mathrm{g \ g^{-1}}`.

   Parameters
   ----------
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   specific_humidity_data_units: :py:class:`str <str>`.
       The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/g`, `g/kg` and so on.
   epsilon: :py:class:`float <float>`.
       The molecular weight ratio, which is molecular weight of the constituent gas to that assumed for air. Defaults to the ratio for water vapor to dry air. (:math:`\epsilon \approx 0.622`)

   Reference
   --------------
   - Hobbs, P. V., and J. M. Wallace, 2006: Atmospheric Science: An Introductory Survey. 2nd ed. Academic Press, 504 pp. https://www.sciencedirect.com/book/9780127329512/atmospheric-science
   - Doswell , C. A., and E. N. Rasmussen, 1994: The Effect of Neglecting the Virtual Temperature Correction on CAPE Calculations. Wea. Forecasting, 9, 625–629, https://doi.org/10.1175/1520-0434(1994)009<0625:TEONTV>2.0.CO;2.
   - https://en.wikipedia.org/wiki/Virtual_temperature
   - https://glossary.ametsoc.org/wiki/Virtual_temperature

   .. seealso::
       - `virtual_temperature — MetPy <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.virtual_temperature.html>`__
       - `temp_virtual - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/temp_virtual.shtml>`__


.. py:function:: calc_virtual_temperature(temper_data: xarray.DataArray, specific_humidity_data: xarray.DataArray, specific_humidity_data_units: Literal['kg/kg', 'g/g', 'g/kg'], epsilon: float = 0.608) -> xarray.DataArray

   Calculate virtual temperature.

   The virtual temperature (:math:`T_v`) is the temperature at which dry air would have the same density as the moist air, at a given pressure.
   In other words, two air samples with the same virtual temperature have the same density, regardless of their actual temperature or relative humidity.
   The virtual temperature is always greater than  the absolute air temperature.

   .. math::
       T_v = T(1+ \epsilon q)

   where :math:`\epsilon = 0.608` when the mixing ratio (specific humidity) :math:`q` is expressed in :math:`\mathrm{g \cdot g^{-1}}`.

   Parameters
   ----------
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   specific_humidity_data_units: :py:class:`str <str>`.
       The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/g`, `g/kg` and so on.
   epsilon: :py:class:`float <float>`.
       A constant.

   Reference
   --------------
   - Doswell , C. A., and E. N. Rasmussen, 1994: The Effect of Neglecting the Virtual Temperature Correction on CAPE Calculations. Wea. Forecasting, 9, 625–629, https://doi.org/10.1175/1520-0434(1994)009<0625:TEONTV>2.0.CO;2.
   - https://en.wikipedia.org/wiki/Virtual_temperature
   - https://glossary.ametsoc.org/wiki/Virtual_temperature

   .. seealso::
       - `virtual_temperature — MetPy <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.virtual_temperature.html>`__
       - `temp_virtual - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/temp_virtual.shtml>`__


.. py:function:: calc_static_stability(temper_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: str) -> xarray.DataArray

   Calculate the static stability within a vertical profile.

   .. math::
       \sigma = - T \frac{\partial \ln \theta}{\partial p}

   Parameters
   ----------
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
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


.. py:function:: calc_dewpoint(vapor_pressure_data: xarray.DataArray, vapor_pressure_data_units: Literal['hPa', 'Pa']) -> xarray.DataArray

   Calculate the ambient dewpoint given the vapor pressure.

   Parameters
   ----------
   vapor_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Water vapor partial pressure.
   total_pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `total_pressure_data` value. Optional values are `hPa`, `Pa`.

   Returns
   -------
   The dew point (:py:class:`xarray.DataArray<xarray.DataArray>`), degrees Celsius.

   .. seealso::
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.dewpoint.html


.. py:function:: calc_mixing_ratio(partial_pressure_data: xarray.DataArray, total_pressure_data: xarray.DataArray, molecular_weight_ratio: float = 0.6219569100577033) -> xarray.DataArray

   Calculate the mixing ratio of a gas.

   This calculates mixing ratio given its partial pressure and the total pressure of the air.
   There are no required units for the input arrays, other than that they have the same units.

   Parameters
   ----------
   partial_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Partial pressure of the constituent gas.
   total_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Total air pressure.
   molecular_weight_ratio : :py:class:`float <float>`, optional.
       The ratio of the molecular weight of the constituent gas to that assumed for air.
       Defaults to the ratio for water vapor to dry air (:math:`\epsilon\approx0.622`).

   .. note::
       The units of `partial_pressure_data` and `total_pressure_data` should be the same.

   Returns
   -------
   The mixing ratio (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless (e.g. Kg/Kg or g/g).

   .. seealso::
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.mixing_ratio.html


.. py:function:: calc_vapor_pressure(pressure_data: xarray.DataArray, mixing_ratio_data: xarray.DataArray, pressure_data_units: Literal['hPa', 'Pa'] = None, epsilon: float = 0.6219569100577033) -> xarray.DataArray

   Parameters
   ----------
   pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The pressure data set.
   mixing_ratio_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The mixing ratio of a gas.
   epsilon: :py:class:`float <float>`.
       The molecular weight ratio, which is molecular weight of the constituent gas to that assumed for air.
       Defaults to the ratio for water vapor to dry air. (:math:`\epsilon \approx 0.622`)
   pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.

   Returns
   -------
   The water vapor (partial) pressure (:py:class:`xarray.DataArray<xarray.DataArray>`), units according to `pressure_data_units`.

   .. seealso::
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.vapor_pressure.html


.. py:function:: calc_saturation_vapor_pressure(temperature_data: xarray.DataArray, temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit']) -> xarray.DataArray

   Calculate the saturation water vapor (partial) pressure.

   Parameters
   ----------
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric temperature.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.

   Returns
   -------
   The mixing ratio (:py:class:`xarray.DataArray<xarray.DataArray>`), hPa.

   .. seealso::
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.saturation_vapor_pressure.html


.. py:function:: calc_saturation_mixing_ratio(total_pressure_data: xarray.DataArray, temperature_data: xarray.DataArray, temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], total_pressure_data_units: Literal['hPa', 'Pa']) -> xarray.DataArray

   Calculate the saturation mixing ratio of water vapor.

   This calculation is given total atmospheric pressure and air temperature.

   Parameters
   ----------
   total_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Total atmospheric pressure.
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric temperature.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
   total_pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `total_pressure_data` value. Optional values are `hPa`, `Pa`.

   Returns
   -------
   The saturation mixing ratio (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless.

   .. seealso::
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.saturation_mixing_ratio.html


.. py:function:: transfer_mixing_ratio_2_specific_humidity(mixing_ratio_data: xarray.DataArray) -> xarray.DataArray

   Calculate the specific humidity from mixing ratio.

   Parameters
   ----------
   mixing_ratio_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The mixing ratio of a gas.

   Returns
   -------
   The specific humidity (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless (e.g. Kg/Kg or g/g).

   .. seealso::
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.specific_humidity_from_mixing_ratio.html


.. py:function:: transfer_specific_humidity_2_mixing_ratio(specific_humidity_data: xarray.DataArray, specific_humidity_data_units: Literal['kg/kg', 'g/g', 'g/kg']) -> xarray.DataArray

   Calculate the mixing ratio from specific humidity.

   Parameters
   ----------
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The Specific humidity of air.
   specific_humidity_data_units: :py:class:`str <str>`.
       The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/g`, `g/kg` and so on.

   Returns
   -------
   The mixing ratio (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless.

   .. seealso::
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.mixing_ratio_from_specific_humidity.html


.. py:function:: transfer_dewpoint_2_specific_humidity(dewpoint_data: xarray.DataArray, pressure_data: xarray.DataArray, dewpoint_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], pressure_data_units: Literal['hPa', 'Pa']) -> xarray.DataArray

   Calculate the specific humidity from the dewpoint temperature and pressure.

   Parameters
   ----------
   dewpoint_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The dewpoint temperature.
   pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The pressure data set.
   dewpoint_data_units: :py:class:`str <str>`.
       The unit corresponding to `dewpoint_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
   pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.

   Returns
   -------
   The specific humidity (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless.

   .. seealso::
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.specific_humidity_from_dewpoint.html


.. py:function:: transfer_specific_humidity_2_dewpoint(specific_humidity_data: xarray.DataArray, pressure_data: xarray.DataArray, specific_humidity_data_units: Literal['kg/kg', 'g/g', 'g/kg'], pressure_data_units: Literal['hPa', 'Pa'], epsilon: float = 0.6219569100577033) -> xarray.DataArray

   Calculate the dewpoint from specific humidity and pressure.

   Parameters
   ----------
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The pressure data set.
   specific_humidity_data_units: :py:class:`str <str>`.
       The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/g`, `g/kg` and so on.
   pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
   epsilon: :py:class:`float <float>`.
       The molecular weight ratio, which is molecular weight of the constituent gas to that assumed for air.
       Defaults to the ratio for water vapor to dry air. (:math:`\epsilon \approx 0.622`)

   Returns
   -------
   The dewpoint (:py:class:`xarray.DataArray<xarray.DataArray>`), degrees Celsius.

   .. seealso::
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.dewpoint_from_specific_humidity.html


.. py:function:: transfer_dewpoint_2_relative_humidity(temperature_data: xarray.DataArray, dewpoint_data: xarray.DataArray, temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], dewpoint_data_units: Literal['celsius', 'kelvin', 'fahrenheit']) -> xarray.DataArray

   Calculate the relative humidity from dewpoint.

   Uses temperature and dewpoint to calculate relative humidity as the ratio of vapor pressure to saturation vapor pressures.

   Parameters
   ----------
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric temperature.
   dewpoint_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The dewpoint temperature.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
   dewpoint_data_units: :py:class:`str <str>`.
       The unit corresponding to `dewpoint_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.

   Returns
   -------
   The relative humidity (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless.

   .. seealso::
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.relative_humidity_from_dewpoint.html


.. py:function:: transfer_mixing_ratio_2_relative_humidity(pressure_data: xarray.DataArray, temperature_data: xarray.DataArray, mixing_ratio_data: xarray.DataArray, pressure_data_units: Literal['hPa', 'Pa'], temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], epsilon: float = 0.6219569100577033) -> xarray.DataArray

   Calculate the relative humidity from mixing ratio, temperature, and pressure.

   Parameters
   ----------
   pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The pressure data set.
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric temperature.
   mixing_ratio_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The mixing ratio of a gas.
   pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
   epsilon: :py:class:`float <float>`.
       The molecular weight ratio, which is molecular weight of the constituent gas to that assumed for air.
       Defaults to the ratio for water vapor to dry air. (:math:`\epsilon \approx 0.622`)

   Returns
   -------
   The relative humidity (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless.

   .. seealso::
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.relative_humidity_from_mixing_ratio.html


.. py:function:: transfer_specific_humidity_2_relative_humidity(pressure_data: xarray.DataArray, temperature_data: xarray.DataArray, specific_humidity_data: xarray.DataArray, pressure_data_units: Literal['hPa', 'Pa'], temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], specific_humidity_data_units: Literal['kg/kg', 'g/g', 'g/kg']) -> xarray.DataArray

   Calculate the relative humidity from specific humidity, temperature, and pressure.

   Parameters
   ----------
   pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The pressure data set.
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric temperature.
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
   specific_humidity_data_units: :py:class:`str <str>`.
       The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/g`, `g/kg` and so on.

   Returns
   -------
   The relative humidity (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless.

   .. seealso::
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.relative_humidity_from_specific_humidity.html


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


.. py:function:: calc_p_gradient(data_input: xarray.DataArray, vertical_dim: str, vertical_dim_units: str) -> xarray.DataArray

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


.. py:function:: calc_top2surface_integral(data_input: xarray.DataArray, surface_pressure_data: xarray.DataArray, vertical_dim: str, surface_pressure_data_units: str, vertical_dim_units: str, method: str = 'Trenberth-vibeta', normalize: bool = True) -> xarray.DataArray

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


.. py:function:: calc_geostrophic_wind_vorticity(z_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', spherical_coord: bool = True, omega: float = 7.292e-05, g: float = 9.8, R: float = 6370000) -> xarray.DataArray

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


.. py:function:: calc_water_flux_top2surface_integral(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, surface_pressure_data: xarray.DataArray, surface_pressure_data_units: str, vertical_dim: str, vertical_dim_units: str, method: str = 'Trenberth-vibeta', g: float = 9.8) -> xarray.DataArray

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


.. py:function:: calc_divergence_watervaporflux_top2surface_integral(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, surface_pressure_data: xarray.DataArray, vertical_dim: str, specific_humidity_data_units: str, surface_pressure_data_units: str, vertical_dim_units: str, spherical_coord: bool = True, lon_dim: str = 'lon', lat_dim: str = 'lat', method: str = 'Trenberth-vibeta', g: float = 9.8, R: float = 6370000) -> xarray.DataArray

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


.. py:function:: calc_eady_growth_rate(u_daily_data: xarray.DataArray, z_daily_data: xarray.DataArray, temper_daily_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: str, lat_dim='lat', g=9.8) -> xarray.Dataset

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
   The maximum Eady growth rate. (:py:class:`xarray.Dataset<xarray.Dataset>`).

   - `eady_growth_rate`: The maximum Eady growth rate.
   - `dudz`: :math:`\frac{\mathrm{d} U}{\mathrm{d} z}`
   - `brunt_vaisala_frequency`: Brunt-väisälä frequency.

   .. seealso::
       - `eady_growth_rate -NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/eady_growth_rate.shtml>`__
       - `瞬变涡旋诊断量 <https://renqlsysu.github.io/2020/02/16/wave_activity_flux/>`__


.. py:function:: calc_apparent_heat_source(u_data: xarray.DataArray, v_data: xarray.DataArray, omega_data: xarray.DataArray, temper_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: str, time_units: str, lon_dim='lon', lat_dim='lat', time_dim='time', c_p=1005.7) -> xarray.DataArray

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


.. py:function:: calc_total_diabatic_heating(u_data: xarray.DataArray, v_data: xarray.DataArray, omega_data: xarray.DataArray, temper_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: str, time_units: str, lat_dim='lat', lon_dim='lon', time_dim='time', c_p=1005.7) -> xarray.DataArray

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


.. py:function:: calc_apparent_moisture_sink(u_data: xarray.DataArray, v_data: xarray.DataArray, omega_data: xarray.DataArray, specific_humidity_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: str, time_units: str, specific_humidity_data_units: str, lon_dim='lon', lat_dim='lat', time_dim='time', latent_heat_of_condensation=2501000.0) -> xarray.DataArray

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


.. py:function:: calc_Plumb_wave_activity_horizontal_flux(z_prime_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: str, lon_dim='lon', lat_dim='lat', omega=7.292e-05, g=9.8, R=6370000) -> xarray.Dataset

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


.. py:function:: calc_TN_wave_activity_horizontal_flux(z_prime_data: xarray.DataArray, u_climatology_data: xarray.DataArray, v_climatology_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: str, lon_dim: str = 'lon', lat_dim: str = 'lat', omega: float = 7.292e-05, g: float = 9.8, R: float = 6370000) -> xarray.DataArray

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
       - http://500hpa.cn/pyinmet/tnflux/
       - http://tytd.gx.cn/exchange/tnflux/
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

   Calculate a linear least-squares regression for spatial data of time.

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
   dim: :py:class:`str <str>`
       Dimension(s) over which to detrend. By default dimension is applied over the `time` dimension.

   Returns
   -------
   - :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::
       :py:func:`scipy.signal.detrend <scipy:scipy.signal.detrend>`.


.. py:function:: calc_ttestSpatialPattern_spatial(data_input1: xarray.DataArray, data_input2: xarray.DataArray, dim: str = 'time') -> xarray.Dataset

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

   Returns
   -------
   - **statistic**, **pvalue**: :py:class:`xarray.Dataset<xarray.Dataset>`.

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


.. py:function:: calc_seasonal_mean(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', extract_season=None, **kwargs) -> xarray.DataArray

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
   G. R. Ackermann. (1983). Means and Standard Deviations of Horizontal Wind Components.
   Website: https://doi.org/10.1175/1520-0450(1983)022%3C0959:MASDOH%3E2.0.CO;2


.. py:function:: populate_monmean2everymon(data_monthly: xarray.DataArray, data_climatology_monthly_data: xarray.DataArray = None, time_dim: str = 'time') -> xarray.DataArray

   Populate the data of each month using the monthly mean state of the `data_monthly` or given dataset.

   Parameters
   ----------
   - data_monthly: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   - data_climatology_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`, default `None`.
       The monthly climatology dataset. If it is `None`, the climatology is derived from `data_monthly`.
   - time_dim: :py:class:`str <str>`, default: `time`.
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
       If True, will print a progress bar of the download to standard error (stderr). Requires `tqdm` to be installed.
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


.. py:data:: __version__
   :value: '2024.11.0'


