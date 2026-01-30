easyclimate.physics.energy.angmom_atm
=====================================

.. py:module:: easyclimate.physics.energy.angmom_atm

.. autoapi-nested-parse::

   Atmospheric Relative Angular Momentum



Functions
---------

.. autoapisummary::

   easyclimate.physics.energy.angmom_atm.calc_relative_angular_momentum


Module Contents
---------------

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


