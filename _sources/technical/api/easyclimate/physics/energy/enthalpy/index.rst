easyclimate.physics.energy.enthalpy
===================================

.. py:module:: easyclimate.physics.energy.enthalpy

.. autoapi-nested-parse::

   Atmospheric Enthalpy



Functions
---------

.. autoapisummary::

   easyclimate.physics.energy.enthalpy.calc_enthalpy


Module Contents
---------------

.. py:function:: calc_enthalpy(temperature_data: xarray.DataArray, mixing_ratio_data: xarray.DataArray, temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], mixing_ratio_data_units: Literal['kg/kg', 'g/g', 'g/kg']) -> xarray.DataArray

   Calculate atmospheric enthalpy from temperature and humidity mixing ratio.

   Enthalpy is a thermodynamic quantity equivalent to the internal energy plus
   the energy the system exerts on its surroundings. The enthalpy is a constant
   pressure function. As such, it includes the work term for expansion
   against the atmosphere.

   .. math::

       T \cdot (1.01 + 0.00189 \cdot W) + 2.5 \cdot W

   where the unit of :math:`T` (atmospheric temperature) is ``degC``, and :math:`W` (mixing ratio) is ``g/kg``.

   Parameters
   -----------
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric temperature.
   mixing_ratio_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The mixing ratio of a gas.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
   mixing_ratio_data_units: :py:class:`str <str>`.
       The unit corresponding to ``mixing_ratio_data`` value. Optional values are :math:`\mathrm{kg/kg}`, :math:`\mathrm{g/g}`, :math:`\mathrm{g/kg}` and so on.

   Returns
   --------
   Atmospheric enthalpy ( :math:`\mathrm{kJ/kg}` ).
       :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::

       - https://www.ncl.ucar.edu/Document/Functions/Contributed/enthalpy.shtml


