easyclimate.physics.temperature.equivalent_potential_temperature
================================================================

.. py:module:: easyclimate.physics.temperature.equivalent_potential_temperature

.. autoapi-nested-parse::

   Equivalent Potential Temperature



Functions
---------

.. autoapisummary::

   easyclimate.physics.temperature.equivalent_potential_temperature.calc_equivalent_potential_temperature
   easyclimate.physics.temperature.equivalent_potential_temperature.calc_equivalent_potential_temperature_davies_jones2009


Module Contents
---------------

.. py:function:: calc_equivalent_potential_temperature(pressure_data: xarray.DataArray, temperature_data: xarray.DataArray, dewpoint_data: xarray.DataArray, pressure_data_units: Literal['hPa', 'Pa'], temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], dewpoint_data_units: Literal['celsius', 'kelvin', 'fahrenheit']) -> xarray.DataArray

   Calculate equivalent potential temperature using Bolton (1980) approximation.


   Parameters
   ----------
   pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The pressure data set.
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric temperature.
   dewpoint_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The dew point temperature.
   pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
   dewpoint_data_units: :py:class:`str <str>`.
       The unit corresponding to `dewpoint_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.

   Returns
   -------
   Equivalent potential temperature ( :math:`\mathrm{K}` ).
       :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::
       - Bolton, D. (1980). The Computation of Equivalent Potential Temperature. Monthly Weather Review, 108(7), 1046-1053. https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_1046_tcoept_2_0_co_2.xml


.. py:function:: calc_equivalent_potential_temperature_davies_jones2009(pressure_data: xarray.DataArray, temperature_data: xarray.DataArray, dewpoint_data: xarray.DataArray, pressure_data_units: Literal['hPa', 'Pa'], temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], dewpoint_data_units: Literal['celsius', 'kelvin', 'fahrenheit']) -> xarray.DataArray

   Calculate equivalent potential temperature using Robert Davies-Jones (2009) approximation.


   Parameters
   ----------
   pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The pressure data set.
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric temperature.
   dewpoint_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The dew point temperature.
   pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
   dewpoint_data_units: :py:class:`str <str>`.
       The unit corresponding to `dewpoint_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.

   Returns
   -------
   Equivalent potential temperature ( :math:`\mathrm{K}` ).
       :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::
       - Davies-Jones, R. (2009). On Formulas for Equivalent Potential Temperature. Monthly Weather Review, 137(9), 3137-3148. https://doi.org/10.1175/2009MWR2774.1


