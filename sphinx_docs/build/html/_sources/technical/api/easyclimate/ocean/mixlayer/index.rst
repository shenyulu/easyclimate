:py:mod:`easyclimate.ocean.mixlayer`
====================================

.. py:module:: easyclimate.ocean.mixlayer

.. autoapi-nested-parse::

   Functions for calculation of ocean mixed layer variables.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.ocean.mixlayer.calc_MLD
   easyclimate.ocean.mixlayer.calc_MLD_temper_tendency
   easyclimate.ocean.mixlayer.calc_MLD_depth_weighted
   easyclimate.ocean.mixlayer.calc_temper_MLD
   easyclimate.ocean.mixlayer.calc_Horizontal_advection
   easyclimate.ocean.mixlayer.calc_Vertical_advection
   easyclimate.ocean.mixlayer.calc_Heat_flux



.. py:function:: calc_MLD()


.. py:function:: calc_MLD_temper_tendency(temper_anomaly_vertical, mld, depth_weight, depth_dim='depth', time_dim='month')

       
       


.. py:function:: calc_MLD_depth_weighted(temper_vertical, mld, depth_dim='depth')

       
       


.. py:function:: calc_temper_MLD(temper_vertical, mld, depth_dim='depth')

       
       


.. py:function:: calc_Horizontal_advection(u, v, temper_mld, depth_weight, lat_dim='lat', lon_dim='lon', depth_dim='depth', R=6370000)

       
       


.. py:function:: calc_Vertical_advection(w, temper_vertical, mld, depth_weight, depth_dim='depth')

       
       


.. py:function:: calc_Heat_flux(qnet_anomaly, mld, rho_0=1027, c_p=4007)

       
       


