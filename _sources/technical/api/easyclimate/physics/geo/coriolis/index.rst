easyclimate.physics.geo.coriolis
================================

.. py:module:: easyclimate.physics.geo.coriolis

.. autoapi-nested-parse::

   Coriolis Force



Functions
---------

.. autoapisummary::

   easyclimate.physics.geo.coriolis.get_coriolis_parameter


Module Contents
---------------

.. py:function:: get_coriolis_parameter(lat_data: xarray.DataArray | numpy.array, omega: float = 7.292e-05) -> xarray.DataArray | numpy.array

   Calculate the Coriolis parameter at each point.

   .. math::
       f = 2 \Omega \sin(\phi)

   Parameters
   ----------
   lat_data: :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`numpy.array <numpy.array>`.
       Latitude at each point.
   omega: :py:class:`float <float>`, default: `7.292e-5` ( :math:`\mathrm{rad/s}` ).
       The angular speed of the earth.

   Returns
   -------
   Corresponding Coriolis force at each point ( :math:`\mathrm{s^{-1}}` ).
       :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`numpy.array <numpy.array>`.

   Reference
   --------------
   - `Coriolis parameter - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Coriolis_parameter>`__

   .. seealso::
       - `coriolis_parameter — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.coriolis_parameter.html>`__
       - `coriolis_param - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/coriolis_param.shtml>`__

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_co_coeff.py


