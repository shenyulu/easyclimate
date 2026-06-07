easyclimate.physics.moisture.lapse_rate
=======================================

.. py:module:: easyclimate.physics.moisture.lapse_rate

.. autoapi-nested-parse::

   Moist Adiabatic Lapse Rate



Functions
---------

.. autoapisummary::

   easyclimate.physics.moisture.lapse_rate.calc_moist_adiabatic_lapse_rate


Module Contents
---------------

.. py:function:: calc_moist_adiabatic_lapse_rate(pressure_data: xarray.DataArray, temperature_data: xarray.DataArray, pressure_data_units: Literal['hPa', 'Pa', 'mbar'], temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit']) -> xarray.DataArray

   Calculate moist adiabatic lapse rate.

   Parameters
   -----------
   pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The pressure data set.
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric temperature.
   pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.

   Returns:
   --------
   dtdp : :py:class:`xarray.DataArray<xarray.DataArray>` ( :math:`\mathrm{K/hPa}` ).
       Moist adiabatic lapse rate.


