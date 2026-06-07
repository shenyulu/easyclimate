easyclimate.physics.moisture.mix
================================

.. py:module:: easyclimate.physics.moisture.mix

.. autoapi-nested-parse::

   Mixing Ratio of a Gas



Functions
---------

.. autoapisummary::

   easyclimate.physics.moisture.mix.calc_mixing_ratio
   easyclimate.physics.moisture.mix.calc_saturation_mixing_ratio


Module Contents
---------------

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
   The mixing ratio ( :math:`\mathrm{g/g}` ).
       :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.mixing_ratio.html


.. py:function:: calc_saturation_mixing_ratio(total_pressure_data: xarray.DataArray, temperature_data: xarray.DataArray, temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], total_pressure_data_units: Literal['hPa', 'Pa', 'mbar']) -> xarray.DataArray

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
   The saturation mixing ratio ( :math:`\mathrm{g/g}` ).
       :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.saturation_mixing_ratio.html


