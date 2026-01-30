easyclimate.physics.geo
=======================

.. py:module:: easyclimate.physics.geo


Submodules
----------

.. toctree::
   :maxdepth: 1

   /technical/api/easyclimate/physics/geo/coriolis/index
   /technical/api/easyclimate/physics/geo/linrood_latwgt/index


Functions
---------

.. autoapisummary::

   easyclimate.physics.geo.get_coriolis_parameter
   easyclimate.physics.geo.calc_lat_weight_lin_rood


Package Contents
----------------

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


.. py:function:: calc_lat_weight_lin_rood(nlat: int) -> Tuple[numpy.ndarray, numpy.ndarray]

   Calculate the latitudes and weights used by the Lin-Rood model.

   The Lin-Rood model requires a specific distribution of latitudes and corresponding weights
   for numerical integration on a spherical grid. This function generates these values based
   on the number of desired latitudes.

   Parameters
   ----------
   nlat : :py:class:`int <int>`
       Number of latitudes. Must be at least 2 to define a valid grid (from pole to pole).

   Returns
   -------
   Tuple[ndarray, ndarray]
       A tuple containing two numpy arrays:
       - lat : ndarray
           Array of latitudes in degrees, ranging from -90 (South Pole) to 90 (North Pole).
       - weight : ndarray
           Array of weights corresponding to each latitude, used for numerical integration.

   .. tip::

       The weights are computed such that they are suitable for use in the Lin-Rood semi-Lagrangian
       transport scheme. The latitudes are uniformly spaced between the poles.

   References
   ----------
   - Lin, S., & Rood, R. B. (1996). Multidimensional Flux-Form Semi-Lagrangian Transport Schemes. Monthly Weather Review, 124(9), 2046-2070. https://journals.ametsoc.org/view/journals/mwre/124/9/1520-0493_1996_124_2046_mffslt_2_0_co_2.xml
   - Lin, S.-J. and Rood, R.B. (1997), An explicit flux-form semi-lagrangian shallow-water model on the sphere. Q.J.R. Meteorol. Soc., 123: 2477-2498. https://doi.org/10.1002/qj.49712354416

   .. seealso::

       - https://www.ncl.ucar.edu/Document/Functions/Built-in/linrood_latwgt.shtml


