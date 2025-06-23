easyclimate.physics.transfer
============================

.. py:module:: easyclimate.physics.transfer

.. autoapi-nested-parse::

   Unit conversion for Meteorological variables



Functions
---------

.. autoapisummary::

   easyclimate.physics.transfer.transfer_mixing_ratio_2_specific_humidity
   easyclimate.physics.transfer.transfer_specific_humidity_2_mixing_ratio
   easyclimate.physics.transfer.transfer_dewpoint_2_specific_humidity
   easyclimate.physics.transfer.transfer_dewpoint_2_mixing_ratio
   easyclimate.physics.transfer.transfer_specific_humidity_2_dewpoint
   easyclimate.physics.transfer.transfer_dewpoint_2_relative_humidity
   easyclimate.physics.transfer.transfer_mixing_ratio_2_relative_humidity
   easyclimate.physics.transfer.transfer_specific_humidity_2_relative_humidity


Module Contents
---------------

.. py:function:: transfer_mixing_ratio_2_specific_humidity(mixing_ratio_data: xarray.DataArray) -> xarray.DataArray

   Calculate the specific humidity from mixing ratio.

   Parameters
   ----------
   mixing_ratio_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The mixing ratio of a gas.

   Returns
   -------
   The specific humidity (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless (e.g. :math:`\mathrm{kg/kg}`, :math:`\mathrm{g/g}`).

   .. seealso::
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.specific_humidity_from_mixing_ratio.html


.. py:function:: transfer_specific_humidity_2_mixing_ratio(specific_humidity_data: xarray.DataArray, specific_humidity_data_units: Literal['kg/kg', 'g/g', 'g/kg']) -> xarray.DataArray

   Calculate the mixing ratio from specific humidity.

   Parameters
   ----------
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The Specific humidity of air.
   specific_humidity_data_units: :py:class:`str <str>`.
       The unit corresponding to `specific_humidity` value. Optional values are :math:`\mathrm{kg/kg}`, :math:`\mathrm{g/g}`, :math:`\mathrm{g/kg}` and so on.

   Returns
   -------
   The mixing ratio (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless.

   .. seealso::
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.mixing_ratio_from_specific_humidity.html


.. py:function:: transfer_dewpoint_2_specific_humidity(dewpoint_data: xarray.DataArray, pressure_data: xarray.DataArray, dewpoint_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], pressure_data_units: Literal['hPa', 'Pa', 'mbar']) -> xarray.DataArray

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


.. py:function:: transfer_dewpoint_2_mixing_ratio(dewpoint_data: xarray.DataArray, pressure_data: xarray.DataArray, dewpoint_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], pressure_data_units: Literal['hPa', 'Pa', 'mbar'])

   Calculate the mixing ratio from the dewpoint temperature and pressure.

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
   The mixing ratio (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless.


.. py:function:: transfer_specific_humidity_2_dewpoint(specific_humidity_data: xarray.DataArray, pressure_data: xarray.DataArray, specific_humidity_data_units: Literal['kg/kg', 'g/g', 'g/kg'], pressure_data_units: Literal['hPa', 'Pa', 'mbar'], epsilon: float = 0.6219569100577033) -> xarray.DataArray

   Calculate the dewpoint from specific humidity and pressure.

   Parameters
   ----------
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The pressure data set.
   specific_humidity_data_units: :py:class:`str <str>`.
       The unit corresponding to `specific_humidity` value. Optional values are :math:`\mathrm{kg/kg}`, :math:`\mathrm{g/g}`, :math:`\mathrm{g/kg}` and so on.
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


.. py:function:: transfer_mixing_ratio_2_relative_humidity(pressure_data: xarray.DataArray, temperature_data: xarray.DataArray, mixing_ratio_data: xarray.DataArray, pressure_data_units: Literal['hPa', 'Pa', 'mbar'], temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], epsilon: float = 0.6219569100577033) -> xarray.DataArray

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


.. py:function:: transfer_specific_humidity_2_relative_humidity(pressure_data: xarray.DataArray, temperature_data: xarray.DataArray, specific_humidity_data: xarray.DataArray, pressure_data_units: Literal['hPa', 'Pa', 'mbar'], temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], specific_humidity_data_units: Literal['kg/kg', 'g/g', 'g/kg']) -> xarray.DataArray

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
       The unit corresponding to `specific_humidity` value. Optional values are :math:`\mathrm{kg/kg}`, :math:`\mathrm{g/g}`, :math:`\mathrm{g/kg}` and so on.

   Returns
   -------
   The relative humidity (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless.

   .. seealso::
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.relative_humidity_from_specific_humidity.html


