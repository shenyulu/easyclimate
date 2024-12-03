easyclimate.field.teleconnection.index_eu
=========================================

.. py:module:: easyclimate.field.teleconnection.index_eu

.. autoapi-nested-parse::

   The Eurasian pattern (EU) pattern



Functions
---------

.. autoapisummary::

   easyclimate.field.teleconnection.index_eu.calc_index_EU_Wallace_Gutzler_1981


Module Contents
---------------

.. py:function:: calc_index_EU_Wallace_Gutzler_1981(z_monthly_data: xarray.DataArray, time_range: slice = slice(None, None), lon_dim: str = 'lon', lat_dim: str = 'lat', time_dim: str = 'time') -> xarray.DataArray

   The calculation of monthly mean Eurasian pattern (EU) index using Pointwise method following Wallace and Gutzler (1981):

   .. math::
       \mathrm{EU = - \frac{1}{4} Z^*(55^{\circ}N, 20^{\circ}E) + \frac{1}{2} Z^*(55^{\circ}N, 75^{\circ}E) - \frac{1}{4} Z^*(40^{\circ}N, 145^{\circ}E)}

   where :math:`Z^*` denotes monthly mean 500 hPa height anomaly.

   Parameters
   ----------
   z_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The monthly geopotential height. The 500hPa layer is recommended.
   time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
       The time range of seasonal cycle means to be calculated. The default value is the entire time range.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name.

   Returns
   -------
   The monthly mean EU index (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - `Wallace, J. M., & Gutzler, D. S. (1981). Teleconnections in the Geopotential Height Field during the Northern Hemisphere Winter. Monthly Weather Review, 109(4), 784-812. <https://journals.ametsoc.org/view/journals/mwre/109/4/1520-0493_1981_109_0784_titghf_2_0_co_2.xml>`__


