easyclimate.physics.moisture.vapor_pressure
===========================================

.. py:module:: easyclimate.physics.moisture.vapor_pressure

.. autoapi-nested-parse::

   Vapor Pressure



Functions
---------

.. autoapisummary::

   easyclimate.physics.moisture.vapor_pressure.calc_vapor_pressure
   easyclimate.physics.moisture.vapor_pressure.calc_saturation_vapor_pressure


Module Contents
---------------

.. py:function:: calc_vapor_pressure(pressure_data: xarray.DataArray, mixing_ratio_data: xarray.DataArray, pressure_data_units: Literal['hPa', 'Pa', 'mbar'] = None, epsilon: float = 0.6219569100577033) -> xarray.DataArray

   Calculate the vapor pressure.

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
   The water vapor (partial) pressure, units according to ``pressure_data_units``.
       :py:class:`xarray.DataArray<xarray.DataArray>`

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
   The saturation water vapor (partial) pressure ( :math:`\mathrm{hPa}` ).
       :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::
       - Bolton, D. (1980). The Computation of Equivalent Potential Temperature. Monthly Weather Review, 108(7), 1046-1053. https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_1046_tcoept_2_0_co_2.xml
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.saturation_vapor_pressure.html


