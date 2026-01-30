easyclimate.physics.geo.linrood_latwgt
======================================

.. py:module:: easyclimate.physics.geo.linrood_latwgt

.. autoapi-nested-parse::

   Latitudes and weights in the Lin-Rood model

   This module provides functions to calculate the latitudes and weights used in the Lin-Rood model,
   which is commonly employed in atmospheric modeling and semi-Lagrangian transport schemes.



Functions
---------

.. autoapisummary::

   easyclimate.physics.geo.linrood_latwgt.calc_lat_weight_lin_rood


Module Contents
---------------

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


