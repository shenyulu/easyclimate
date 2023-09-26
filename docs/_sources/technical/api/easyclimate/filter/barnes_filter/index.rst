:py:mod:`easyclimate.filter.barnes_filter`
==========================================

.. py:module:: easyclimate.filter.barnes_filter

.. autoapi-nested-parse::

   https://github.com/LinOuyang/pybarnes



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   easyclimate.filter.barnes_filter.BarnesFilter



Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.filter.barnes_filter.field_grids



.. py:function:: field_grids(data, grids)

   For each point of the grids, its nearby region are chosen.

   Parameters
   ----------
   data : array
       An N-dimensional array.
   grids : int or 2-size tuple
       The total number of grid points of the x and y directions.

   Returns
   -------
   D : ndarray
       A ndarray with extra two dimensions. The last two dimensions are
       each points's nearby region.


.. py:class:: BarnesFilter(data_arr, lon=None, lat=None, radius_degree=10)

   The Barnes method performs grid point interpolation by selecting appropriate
   filtering parameters *c* and *g* to filter out shortwave noise in the original field,
   making the analysis results stable and smooth. In addition, it can form a bandpass filter
   to separate various sub weather scales that affect weather processes according to actual needs,
   achieving the purpose of scale separation.

   Reference:
   DOI : https://doi.org/10.1175/1520-0493(1980)108<1108:AOTFSM>2.0.CO;2

   .. py:method:: __convert_data(data)


   .. py:method:: __calculate_distance(lon, lat)


   .. py:method:: __lowpass(g=0.3, c=150000)


   .. py:method:: lowpass(g=0.3, c=150000)

      Selecting different parameters *g* and *c*
      will result in different filtering characteristics.

      Reference:
      DOI : https://doi.org/10.1175/1520-0493(1980)108<1108:AOTFSM>2.0.CO;2

      Parameters
      ----------
      g : float, generally between (0, 1]
          Constant parameter.
      c : int
          Constant parameter. When *c* takes a larger value, the filter function converges
          at a larger wavelength, and the response function slowly approaches the maximum value,
          which means that high-frequency fluctuations have been filtered out.

      Returns
      -------
      data_vars : array
          Data field after filtering out high-frequency fluctuations



   .. py:method:: bandpass(g1=0.3, c1=30000, g2=0.3, c2=150000)

      Select two different filtering schemes 1 and 2, and perform the filtering separately.
      And then perform the difference, that means *scheme1 - scheme2*.
      The mesoscale fluctuations are thus preserved.

      Parameters
      ----------
      g1 : float, generally between (0, 1]
          Constant parameter of scheme1.
      c1 : int
          Constant parameterof scheme1.
      g2 : float, generally between (0, 1]
          Constant parameter of scheme2.
      c2 : int
          Constant parameterof scheme2.

      Returns
      -------
      data_vars : array
          Mesoscale wave field filtered out from raw data



