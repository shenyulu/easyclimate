easyclimate.physics.energy
==========================

.. py:module:: easyclimate.physics.energy


Submodules
----------

.. toctree::
   :maxdepth: 1

   /technical/api/easyclimate/physics/energy/angmom_atm/index
   /technical/api/easyclimate/physics/energy/enthalpy/index
   /technical/api/easyclimate/physics/energy/latent_heat_water/index


Functions
---------

.. autoapisummary::

   easyclimate.physics.energy.calc_enthalpy
   easyclimate.physics.energy.calc_latent_heat_water
   easyclimate.physics.energy.calc_relative_angular_momentum


Package Contents
----------------

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


.. py:function:: calc_relative_angular_momentum(zonal_wind_speed_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], lon_dim: str = 'lon', lat_dim: str = 'lat', weights=None)

   Calculate atmospheric relative angular momentum.

   Parameters
   -----------
   zonal_wind_speed_data : :py:class:`xarray.DataArray <xarray.DataArray>` ( :math:`\mathrm{m/s}` )
       Zonal wind component with the least similar dimensions ``(vertical_dim, lon_dim, lat_dim)``
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are ``hPa``, ``Pa``, ``mbar``.
   lon_dim: :py:class:`str <str>`, default: ``lon``.
       Longitude coordinate dimension name. By default extracting is applied over the ``lon`` dimension.
   lat_dim: :py:class:`str <str>`, default: ``lat``.
       Latitude coordinate dimension name. By default extracting is applied over the ``lat`` dimension.
   weights : :py:class:`xarray.DataArray <xarray.DataArray>`, optional
       Weights for each latitude, same dimension as lat.
       If None, computed as :math:`\cos(lat)*\mathrm{d}lat` with the values of ``latitude`` spacing.

   Returns
   --------
   aam : :py:class:`xarray.DataArray <xarray.DataArray>` ( :math:`\mathrm{kg} \cdot \mathrm{m^2/s}` )
       Atmospheric angular momentum.

   .. seealso::

       - https://www.ncl.ucar.edu/Document/Functions/Built-in/angmom_atm.shtml


