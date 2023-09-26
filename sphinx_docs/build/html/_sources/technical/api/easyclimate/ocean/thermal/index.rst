:py:mod:`easyclimate.ocean.thermal`
===================================

.. py:module:: easyclimate.ocean.thermal

.. autoapi-nested-parse::

   Functions for calculation of ocean thermocline variables.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.ocean.thermal.calc_temp_thermocline_depth
   easyclimate.ocean.thermal.calc_Dx_depth
   easyclimate.ocean.thermal.calc_D14_depth
   easyclimate.ocean.thermal.calc_D17_depth
   easyclimate.ocean.thermal.calc_D20_depth
   easyclimate.ocean.thermal.calc_D26_depth
   easyclimate.ocean.thermal.calc_D28_depth



.. py:function:: calc_temp_thermocline_depth(temp_input, depth_dim='depth')

   Caculate thermocline depth of ocean temperature.

   Parameters
   ----------
   - temp_input: :py:class:`xarray.DataArray<xarray.DataArray>`
       ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
   - depth_dim: str
       `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: calc_Dx_depth(data_input, var, depth_dim='depth', mask=True)

       
       


.. py:function:: calc_D14_depth(data_input, var=14, depth_dim='depth')

       
       


.. py:function:: calc_D17_depth(data_input, var=17, depth_dim='depth')

       
       


.. py:function:: calc_D20_depth(data_input, var=20, depth_dim='depth')

       
       


.. py:function:: calc_D26_depth(data_input, var=26, depth_dim='depth')

       
       


.. py:function:: calc_D28_depth(data_input, var=28, depth_dim='depth')

       
       


