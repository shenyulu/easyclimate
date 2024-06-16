easyclimate.field.ocean.thermal
===============================

.. py:module:: easyclimate.field.ocean.thermal

.. autoapi-nested-parse::

   The calculation of ocean thermocline variables.



Functions
---------

.. autoapisummary::

   easyclimate.field.ocean.thermal.calc_gradient
   easyclimate.field.ocean.thermal.calc_seawater_thermocline_depth
   easyclimate.field.ocean.thermal.calc_Dx_depth
   easyclimate.field.ocean.thermal.calc_D14_depth
   easyclimate.field.ocean.thermal.calc_D17_depth
   easyclimate.field.ocean.thermal.calc_D20_depth
   easyclimate.field.ocean.thermal.calc_D26_depth
   easyclimate.field.ocean.thermal.calc_D28_depth


Module Contents
---------------

.. py:function:: calc_gradient(data_input: xarray.DataArray | xarray.Dataset, dim: str, varargs: int = 1, edge_order: int = 2) -> xarray.DataArray | xarray.Dataset

   Compute the gradient along the coordinate `dim` direction.

   The gradient is computed using **second order accurate central differences** in the interior points
   and either first or second order accurate one-sides (forward or backwards) differences at the boundaries.
   The returned gradient hence has the same shape as the input array.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
        The spatio-temporal data to be calculated.
   dim : :py:class:`str <str>`.
       Dimension(s) over which to apply gradient. By default gradient is applied over the `time` dimension.
   varargs: :py:class:`list <list>` of scalar or array, optional
       Spacing between f values. Default unitary spacing for all dimensions. Spacing can be specified using:

       1. Single scalar to specify a sample distance for all dimensions.
       2. N scalars to specify a constant sample distance for each dimension. i.e. :math:`\mathrm{d}x, \mathrm{d}y, \mathrm{d}z, ...`
       3. N arrays to specify the coordinates of the values along each dimension of F.
          The length of the array must match the size of the corresponding dimension.
       4. Any combination of N scalars/arrays with the meaning of 2. and 3.

   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.

   Returns
   -------
   The gradient along the coordinate `dim` direction (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`numpy.gradient <numpy:numpy.gradient>`


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


