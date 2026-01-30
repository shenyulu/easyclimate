easyclimate.physics.energy.latent_heat_water
============================================

.. py:module:: easyclimate.physics.energy.latent_heat_water

.. autoapi-nested-parse::

   Latent Heat Flux for Water



Functions
---------

.. autoapisummary::

   easyclimate.physics.energy.latent_heat_water.calc_latent_heat_water


Module Contents
---------------

.. py:function:: calc_latent_heat_water(temperature_data, temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], latent_heat_type: Literal['evaporation_condensation', 'melting_freezing', 'sublimation_deposition']) -> xarray.DataArray

   Estimate latent heat flux for water: evaporization (condensation), melting (freezing) or sublimation (deposition).

   .. tip::

       This function returns the latent heat of

       - evaporation/condensation
       - melting/freezing
       - sublimation/deposition

       for water. The latent heatis a function of temperature t. The formulas are polynomial approximations
       to the values in Table 92, p. 343 of the Smithsonian Meteorological Tables, Sixth Revised Edition,
       1963 by Roland List. The approximations were developed by Eric Smith at Colorado State University.

       - Source: Thomas W. Schlatter and Donald V. Baker: PROFS Program Office, NOAA Environmental Research Laboratories, Boulder, Colorado.

   Parameters
   -----------
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric temperature.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are ``celsius``, ``kelvin``, ``fahrenheit``.
   latent_heat_type: :py:class:`str <str>`.
       The type of latent heat to estimate. Optional values are ``evaporation_condensation``, ``melting_freezing``, ``sublimation_deposition``.

   Returns
   --------
   Latent heat flux for water ( :math:`\mathrm{J/kg}` ).
       :py:class:`xarray.DataArray <xarray.DataArray>`.

   .. seealso::

       - https://www.ncl.ucar.edu/Document/Functions/Contributed/latent_heat_water.shtml


