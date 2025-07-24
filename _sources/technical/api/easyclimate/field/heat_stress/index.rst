easyclimate.field.heat_stress
=============================

.. py:module:: easyclimate.field.heat_stress


Submodules
----------

.. toctree::
   :maxdepth: 1

   /technical/api/easyclimate/field/heat_stress/humanindexmod_2020/index


Functions
---------

.. autoapisummary::

   easyclimate.field.heat_stress.calc_apparent_temperature
   easyclimate.field.heat_stress.calc_simplified_human_discomfort_index
   easyclimate.field.heat_stress.calc_simplified_human_discomfort_index_stull
   easyclimate.field.heat_stress.calc_swamp_cooler_temperatures
   easyclimate.field.heat_stress.calc_heat_thic_thip
   easyclimate.field.heat_stress.calc_simplified_wbgt_index
   easyclimate.field.heat_stress.calc_human_feels_temperature


Package Contents
----------------

.. py:function:: calc_apparent_temperature(temperature_data: xarray.DataArray, vapor_pressure_data: xarray.DataArray, wind_10m_data: xarray.DataArray, temperature_data_units: Literal['degC', 'degF', 'degK'], vapor_pressure_data_units: Literal['hPa', 'Pa', 'mbar'], wind_10m_data_units: Literal['m/s'] = 'm/s') -> xarray.DataArray

   Calculate apparent temperature.

   Parameters
   ------------
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The temperature(s).
   vapor_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The vapor pressure.
   wind_10m_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The 10-meter winds.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.
   vapor_pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `vapor_pressure_data` value. Optional values are ``Pa``, ``hPa``, ``mbar`` and so on.
   wind_10m_data_units: :py:class:`str <str>`, default ``m/s``.
       The unit corresponding to `wind_10m_data` value. default value is  ``m/s``.

   Returns
   ---------
   The apparent temperature (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - https://www.ncl.ucar.edu/Document/Functions/Heat_stress/heat_apptemp.shtml
   - https://github.com/jrbuzan/HumanIndexMod_2020
   - http://www.bom.gov.au/info/thermal_stress/#atapproximation
   - Buzan, J. R., Oleson, K., and Huber, M.: Implementation and comparison of a suite of heat stress metrics within the Community Land Model version 4.5, Geosci. Model Dev., 8, 151–170, https://doi.org/10.5194/gmd-8-151-2015, 2015.
   - Steadman, R.G., 1994: Norms of apparent temperature in Australia, Aust. Met. Mag., 43, 1-16.


.. py:function:: calc_simplified_human_discomfort_index(temperature_data: xarray.DataArray, vapor_pressure_data: xarray.DataArray, temperature_data_units: Literal['degC', 'degF', 'degK'], vapor_pressure_data_units: Literal['hPa', 'Pa', 'mbar']) -> xarray.DataArray

   Calculate a simplified human discomfort index.

   Parameters
   ------------
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The temperature(s).
   vapor_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The vapor pressure.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.
   vapor_pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `vapor_pressure_data` value. Optional values are ``Pa``, ``hPa``, ``mbar`` and so on.

   Returns
   ---------
   The simplified human discomfort index (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - https://www.ncl.ucar.edu/Document/Functions/Heat_stress/heat_discoi.shtml
   - https://github.com/jrbuzan/HumanIndexMod_2020
   - Buzan, J. R., Oleson, K., and Huber, M.: Implementation and comparison of a suite of heat stress metrics within the Community Land Model version 4.5, Geosci. Model Dev., 8, 151–170, https://doi.org/10.5194/gmd-8-151-2015, 2015.
   - Steadman, R.G., 1994: Norms of apparent temperature in Australia, Aust. Met. Mag., 43, 1-16.


.. py:function:: calc_simplified_human_discomfort_index_stull(temperature_2m_data: xarray.DataArray, stull_wet_bulb_temperature_data: xarray.DataArray, relative_humidity_data: xarray.DataArray, temperature_2m_data_units: Literal['degC', 'degF', 'degK'], stull_wet_bulb_temperature_data_units: Literal['degC', 'degF', 'degK'], relative_humidity_data_units: Literal['%', 'dimensionless']) -> xarray.DataArray

   Calculate the human discomfort index due to excessive heat and humidity using the Stull wet bulb temperature.

   Parameters
   ------------
   temperature_2m_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The 2m temperature(s).
   stull_wet_bulb_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The Stull wet bulb temperature (`wetbulb_stull <https://www.ncl.ucar.edu/Document/Functions/Contributed/wetbulb_stull.shtml>`__).
   relative_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The relative humidity.
   temperature_2m_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_2m_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.
   stull_wet_bulb_temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `stull_wet_bulb_temperature_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.
   relative_humidity_data_units: :py:class:`str <str>`.
       The unit corresponding to `vapor_pressure_data` value. Optional values are ``%``, ``dimensionless``.

   Returns
   ---------
   The human discomfort index due to excessive heat and humidity using the Stull wet bulb temperature (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - https://www.ncl.ucar.edu/Document/Functions/Heat_stress/heat_discoi_stull.shtml
   - https://github.com/jrbuzan/HumanIndexMod_2020
   - Buzan, J. R., Oleson, K., and Huber, M.: Implementation and comparison of a suite of heat stress metrics within the Community Land Model version 4.5, Geosci. Model Dev., 8, 151–170, https://doi.org/10.5194/gmd-8-151-2015, 2015.
   - Epstein, Y., and D.S. Moran (2006) Thermal comfort and the heat stress indices, Ind. Health, 44, 388-398 doi:https://doi.org/10.2486/indhealth.44.388.


.. py:function:: calc_swamp_cooler_temperatures(temperature_data: xarray.DataArray, wet_bulb_temperature_data: xarray.DataArray, temperature_data_units: Literal['degC', 'degF', 'degK'], wet_bulb_temperature_data_units: Literal['degC', 'degF', 'degK']) -> xarray.DataArray

   Calculate the swamp cooler temperatures at 65% amd 80% efficiency.

   Parameters
   ------------
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The temperature(s).
   wet_bulb_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The wet bulb temperature.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.
   wet_bulb_temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `wet_bulb_temperature_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.

   Returns
   ---------
   The swamp cooler temperatures at 65% amd 80% efficiency (:py:class:`xarray.Dataset<xarray.Dataset>`).

   Reference
   --------------
   - https://www.ncl.ucar.edu/Document/Functions/Heat_stress/heat_swamp_cooleff.shtml
   - https://github.com/jrbuzan/HumanIndexMod_2020
   - Buzan, J. R., Oleson, K., and Huber, M.: Implementation and comparison of a suite of heat stress metrics within the Community Land Model version 4.5, Geosci. Model Dev., 8, 151–170, https://doi.org/10.5194/gmd-8-151-2015, 2015.


.. py:function:: calc_heat_thic_thip(temperature_data: xarray.DataArray, wet_bulb_temperature_data: xarray.DataArray, temperature_data_units: Literal['degC', 'degF', 'degK'], wet_bulb_temperature_data_units: Literal['degC', 'degF', 'degK']) -> xarray.DataArray

   Calculate the thermal humidity comfort index (thic) and the thermal humidity physiology index (thip).

   Parameters
   ------------
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The temperature(s).
   wet_bulb_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The wet bulb temperature.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.
   wet_bulb_temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `wet_bulb_temperature_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.

   Returns
   ---------
   The thermal humidity comfort index (thic) and the thermal humidity physiology index (thip).

   Quantified estimates for Comfort (THIC) and Physiology (THIP)

   +------------+------------------------+
   |    THIC    |    Description         |
   +============+========================+
   |   75-78    |    alert               |
   +------------+------------------------+
   |   79-83    |    dangerous           |
   +------------+------------------------+
   |    84+     |    very dangerous      |
   +------------+------------------------+

   Reference
   --------------
   - https://www.ncl.ucar.edu/Document/Functions/Heat_stress/heat_thic_thip.shtml
   - https://github.com/jrbuzan/HumanIndexMod_2020
   - Buzan, J. R., Oleson, K., and Huber, M.: Implementation and comparison of a suite of heat stress metrics within the Community Land Model version 4.5, Geosci. Model Dev., 8, 151–170, https://doi.org/10.5194/gmd-8-151-2015, 2015.


.. py:function:: calc_simplified_wbgt_index(temperature_data: xarray.DataArray, vapor_pressure_data: xarray.DataArray, temperature_data_units: Literal['degC', 'degF', 'degK'], vapor_pressure_data_units: Literal['hPa', 'Pa', 'mbar']) -> xarray.DataArray

   Calculate Simplified WBGT index.

   Parameters
   ------------
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The temperature(s).
   vapor_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The vapor pressure.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.
   vapor_pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `vapor_pressure_data` value. Optional values are ``Pa``, ``hPa``, ``mbar`` and so on.

   Returns
   ---------
   The simplified WBGT index (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - https://www.ncl.ucar.edu/Document/Functions/Heat_stress/heat_wbgt_simplified.shtml
   - https://github.com/jrbuzan/HumanIndexMod_2020
   - Willett, K.M. and Sherwood, S. (2012), Exceedance of heat index thresholds for 15 regions under a warming climate using the wet-bulb globe temperature. Int. J. Climatol., 32: 161-177. https://doi.org/10.1002/joc.2257
   - Buzan, J. R., Oleson, K., and Huber, M.: Implementation and comparison of a suite of heat stress metrics within the Community Land Model version 4.5, Geosci. Model Dev., 8, 151–170, https://doi.org/10.5194/gmd-8-151-2015, 2015.


.. py:function:: calc_human_feels_temperature(temperature_data: xarray.DataArray, vapor_pressure_data: xarray.DataArray, temperature_data_units: Literal['degC', 'degF', 'degK'], vapor_pressure_data_units: Literal['hPa', 'Pa', 'mbar']) -> xarray.DataArray

   Calculate the 'feels-like' temperature for humans.

   Parameters
   ------------
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The temperature(s).
   vapor_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The vapor pressure.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.
   vapor_pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `vapor_pressure_data` value. Optional values are ``Pa``, ``hPa``, ``mbar`` and so on.

   Returns
   ---------
   The 'feels-like' temperature for humans (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - https://www.ncl.ucar.edu/Document/Functions/Heat_stress/heat_humidex.shtml
   - https://github.com/jrbuzan/HumanIndexMod_2020
   - Masterson, J., and F. Richardson, 1979: Humidex, a method of quantifying human discomfort due to excessive heat and humidity CLI 1-79, Environment Canada, Atmosheric Environment Servic website: https://publications.gc.ca/site/eng/9.865813/publication.html, https://publications.gc.ca/collections/collection_2018/eccc/En57-23-1-79-eng.pdf.
   - Buzan, J. R., Oleson, K., and Huber, M.: Implementation and comparison of a suite of heat stress metrics within the Community Land Model version 4.5, Geosci. Model Dev., 8, 151–170, https://doi.org/10.5194/gmd-8-151-2015, 2015.


