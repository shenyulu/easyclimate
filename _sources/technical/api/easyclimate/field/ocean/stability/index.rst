easyclimate.field.ocean.stability
=================================

.. py:module:: easyclimate.field.ocean.stability

.. autoapi-nested-parse::

   The calculation of ocean instability.



Functions
---------

.. autoapisummary::

   easyclimate.field.ocean.stability.find_dims_axis
   easyclimate.field.ocean.stability.calc_N2_from_temp_salt
   easyclimate.field.ocean.stability.calc_potential_density_from_temp_salt


Module Contents
---------------

.. py:function:: find_dims_axis(data: xarray.DataArray, dim: str) -> int

   Find the index of `dim` in the xarray DataArray.

   Parameters
   ----------
   - data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   - dim : :py:class:`str <str>`
       Dimension(s) over which to find axis.

   Returns
   -------
   :py:class:`int <int>`.


.. py:function:: calc_N2_from_temp_salt(seawater_temperature_data: xarray.DataArray, seawater_practical_salinity_data: xarray.DataArray, time_dim: str | None, depth_dim: str = 'depth', lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.Dataset

   Calculate the frequency of seawater buoyancy.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{^\circ C}`)
       Vertical seawater temperature data.
   seawater_practical_salinity_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{PSU}`)
       Vertical seawater salinity data (practical salinity).
   time_dim: :py:class:`str <str>` or `None`, default: `time`.
       The time coordinate dimension name.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The frequency of seawater buoyancy (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::

       - http://www.teos-10.org/pubs/gsw/html/gsw_Nsquared.html
       - http://www.teos-10.org/pubs/gsw/html/gsw_contents.html


.. py:function:: calc_potential_density_from_temp_salt(seawater_temperature_data: xarray.DataArray, seawater_practical_salinity_data: xarray.DataArray, time_dim: str | None, depth_dim: str = 'depth', lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.Dataset

   Calculate the potential density of seawater.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{^\circ C}`)
       Vertical seawater temperature data.
   seawater_practical_salinity_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\mathrm{PSU}`)
       Vertical seawater salinity data (practical salinity).
   time_dim: :py:class:`str <str>` or `None`, default: `time`.
       The time coordinate dimension name.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The potential density of seawater (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::

       - http://www.teos-10.org/pubs/gsw/html/gsw_Nsquared.html
       - http://www.teos-10.org/pubs/gsw/html/gsw_contents.html


