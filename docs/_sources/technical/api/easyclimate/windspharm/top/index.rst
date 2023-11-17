:py:mod:`easyclimate.windspharm.top`
====================================

.. py:module:: easyclimate.windspharm.top

.. autoapi-nested-parse::

   easy cliamte top interface for the windspharm

   https://ajdawson.github.io/windspharm/latest/api/windspharm.xarray.html



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.windspharm.top.calc_wind_speed
   easyclimate.windspharm.top.calc_relative_vorticity_and_horizontal_divergence
   easyclimate.windspharm.top.calc_relative_vorticity
   easyclimate.windspharm.top.calc_divergence
   easyclimate.windspharm.top.calc_planetary_vorticity
   easyclimate.windspharm.top.calc_absolute_vorticity
   easyclimate.windspharm.top.calc_streamfunction_and_velocity_potential
   easyclimate.windspharm.top.calc_streamfunction
   easyclimate.windspharm.top.calc_velocity_potential
   easyclimate.windspharm.top.calc_helmholtz
   easyclimate.windspharm.top.calc_irrotational_component
   easyclimate.windspharm.top.calc_nondivergent_component
   easyclimate.windspharm.top.calc_rossby_wave_source
   easyclimate.windspharm.top.calc_gradient



.. py:function:: calc_wind_speed(u, v, R=6371200.0, legfunc='stored')

       
       


.. py:function:: calc_relative_vorticity_and_horizontal_divergence(u, v, truncation=None, R=6371200.0, legfunc='stored')

       
       


.. py:function:: calc_relative_vorticity(u, v, truncation=None, R=6371200.0, legfunc='stored')

       
       


.. py:function:: calc_divergence(u, v, truncation=None, R=6371200.0, legfunc='stored')

       
       


.. py:function:: calc_planetary_vorticity(u, v, omega=7.292115, R=6371200.0, legfunc='stored')

       
       


.. py:function:: calc_absolute_vorticity(u, v, truncation=None, omega=7.292115, R=6371200.0, legfunc='stored')

       


.. py:function:: calc_streamfunction_and_velocity_potential(u, v, truncation=None, R=6371200.0, legfunc='stored')

       
       


.. py:function:: calc_streamfunction(u, v, truncation=None, R=6371200.0, legfunc='stored')

       
       


.. py:function:: calc_velocity_potential(u, v, truncation=None, R=6371200.0, legfunc='stored')

       
       


.. py:function:: calc_helmholtz(u, v, truncation=None, R=6371200.0, legfunc='stored')

       
       


.. py:function:: calc_irrotational_component(u, v, truncation=None, R=6371200.0, legfunc='stored')

       
       


.. py:function:: calc_nondivergent_component(u, v, truncation=None, R=6371200.0, legfunc='stored')

       
       


.. py:function:: calc_rossby_wave_source(u, v, truncation=None, R=6371200.0, legfunc='stored')

       
       


.. py:function:: calc_gradient(data, truncation=None, R=6371200.0, legfunc='stored')

       
       


