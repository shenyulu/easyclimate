:py:mod:`easyclimate.windspharm._common`
========================================

.. py:module:: easyclimate.windspharm._common

.. autoapi-nested-parse::

   Common functionality shared across interfaces.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.windspharm._common.get_apiorder
   easyclimate.windspharm._common.inspect_gridtype
   easyclimate.windspharm._common.to3d



.. py:function:: get_apiorder(ndim, latitude_dim, longitude_dim)

   Get the dimension ordering for a transpose to the required API
   dimension ordering.

   **Arguments:**

   *ndim*
       Total number of dimensions to consider.

   *latitude_dim*
       Index of the latitude dimension.

   *longitude_dim*
       Index of the longitude dimension.

   **Returns:**

   *apiorder*
       A list of indices corresponding to the order required to
       conform to the specified API order.

   *reorder*
       The inverse indices corresponding to *apiorder*.



.. py:function:: inspect_gridtype(latitudes)

   Determine a grid type by examining the points of a latitude
   dimension.

   Raises a ValueError if the grid type cannot be determined.

   **Argument:**

   *latitudes*
       An iterable of latitude point values.

   **Returns:**

   *gridtype*
       Either 'gaussian' for a Gaussian grid or 'regular' for an
       equally-spaced grid.



.. py:function:: to3d(array)


