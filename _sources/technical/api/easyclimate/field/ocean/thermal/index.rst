easyclimate.field.ocean.thermal
===============================

.. py:module:: easyclimate.field.ocean.thermal

.. autoapi-nested-parse::

   The calculation of ocean thermocline variables.



Functions
---------

.. autoapisummary::

   easyclimate.field.ocean.thermal.calc_seawater_thermocline_depth
   easyclimate.field.ocean.thermal.calc_Dx_depth
   easyclimate.field.ocean.thermal.calc_D14_depth
   easyclimate.field.ocean.thermal.calc_D17_depth
   easyclimate.field.ocean.thermal.calc_D20_depth
   easyclimate.field.ocean.thermal.calc_D26_depth
   easyclimate.field.ocean.thermal.calc_D28_depth


Module Contents
---------------

.. py:function:: calc_seawater_thermocline_depth(seawater_temperature_data: xarray.DataArray, depth_dim: str = 'depth') -> xarray.DataArray

   Caculate thermocline depth of ocean temperature.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\mathrm{^\circ C}`).
       ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: calc_Dx_depth(seawater_temperature_data: xarray.DataArray, value: float, depth_dim: str = 'depth') -> xarray.DataArray

   Caculate `value` depth of ocean temperature.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\mathrm{^\circ C}`).
       ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
   value: float.
       The depth of ocean temperature to be calculated.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: calc_D14_depth(seawater_temperature_data: xarray.DataArray, value: float = 14, depth_dim: str = 'depth') -> xarray.DataArray

   Caculate 14m depth of ocean temperature.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\mathrm{^\circ C}`).
       ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
   value: float (:math:`\mathrm{m}`).
       The depth of ocean temperature to be calculated.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: calc_D17_depth(seawater_temperature_data: xarray.DataArray, value: float = 17, depth_dim: str = 'depth') -> xarray.DataArray

   Caculate 17m depth of ocean temperature.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\mathrm{^\circ C}`).
       ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
   value: float (:math:`\mathrm{m}`).
       The depth of ocean temperature to be calculated.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: calc_D20_depth(seawater_temperature_data: xarray.DataArray, value: float = 20, depth_dim: str = 'depth') -> xarray.DataArray

   Caculate 20m depth of ocean temperature.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\mathrm{^\circ C}`).
       ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
   value: float (:math:`\mathrm{m}`).
       The depth of ocean temperature to be calculated.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: calc_D26_depth(seawater_temperature_data: xarray.DataArray, value: float = 26, depth_dim: str = 'depth') -> xarray.DataArray

   Caculate 26m depth of ocean temperature.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\mathrm{^\circ C}`).
       ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
   value: float (:math:`\mathrm{m}`).
       The depth of ocean temperature to be calculated.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: calc_D28_depth(seawater_temperature_data: xarray.DataArray, value: float = 28, depth_dim: str = 'depth') -> xarray.DataArray

   Caculate 28m depth of ocean temperature.

   Parameters
   ----------
   seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\mathrm{^\circ C}`).
       ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
   value: float (:math:`\mathrm{m}`).
       The depth of ocean temperature to be calculated.
   depth_dim: :py:class:`str <str>`, default: `depth`.
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


