easyclimate.core.diagnosis
==========================

.. py:module:: easyclimate.core.diagnosis

.. autoapi-nested-parse::

   Functions for Weather and climate variable diagnosis.



Functions
---------

.. autoapisummary::

   easyclimate.core.diagnosis.calc_brunt_vaisala_frequency_atm
   easyclimate.core.diagnosis.get_coriolis_parameter
   easyclimate.core.diagnosis.calc_potential_temperature
   easyclimate.core.diagnosis.calc_virtual_temperature
   easyclimate.core.diagnosis.calc_virtual_temperature_Hobbs2006
   easyclimate.core.diagnosis.calc_static_stability
   easyclimate.core.diagnosis.calc_dewpoint
   easyclimate.core.diagnosis.calc_mixing_ratio
   easyclimate.core.diagnosis.calc_vapor_pressure
   easyclimate.core.diagnosis.calc_saturation_vapor_pressure
   easyclimate.core.diagnosis.calc_saturation_mixing_ratio
   easyclimate.core.diagnosis.transfer_mixing_ratio_2_specific_humidity
   easyclimate.core.diagnosis.transfer_specific_humidity_2_mixing_ratio
   easyclimate.core.diagnosis.transfer_dewpoint_2_specific_humidity
   easyclimate.core.diagnosis.transfer_specific_humidity_2_dewpoint
   easyclimate.core.diagnosis.transfer_dewpoint_2_relative_humidity
   easyclimate.core.diagnosis.transfer_mixing_ratio_2_relative_humidity
   easyclimate.core.diagnosis.transfer_specific_humidity_2_relative_humidity


Module Contents
---------------

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


