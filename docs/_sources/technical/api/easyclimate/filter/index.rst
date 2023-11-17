:py:mod:`easyclimate.filter`
============================

.. py:module:: easyclimate.filter


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   barnes_filter/index.rst
   butter_filter/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   easyclimate.filter.BarnesFilter



Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.filter.field_grids
   easyclimate.filter.find_dims_axis
   easyclimate.filter.calc_butter_bandpass
   easyclimate.filter.calc_butter_lowpass
   easyclimate.filter.calc_butter_highpass



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



.. py:function:: find_dims_axis(data: xarray.DataArray, dim: str) -> int

   Find the index of `dim` in the xarray DataArray.

   Parameters
   ----------
   - data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   - dim : :py:class:`str<python.str>`
       Dimension(s) over which to find axis.

   Returns
   -------
   :py:class:`int<python.int>`.


.. py:function:: calc_butter_bandpass(data: xarray.DataArray, sampling_frequency: int, period: list, N=3, dim='time') -> xarray.DataArray

   Butterworth bandpass filter.

   Parameters
   ----------
   - data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The array of data to be filtered.
   - sampling_frequency: int.
       Data sampling frequency. If it is daily data with only one time level record per day, 
       then the parameter value is 1; If a day contains four time level records, the parameter value is 4.
   - period: list.
       The time period interval of the bandpass filter to be acquired. 
       If we are obtaining a 3-10 day bandpass filter, the value of this parameter is `[3, 10]`. 
       Note that the units of this parameter should be consistent with `sampling_frequency`.
   - N: int.
       The order of the filter. Default is 3.
   - dim: str.
       Dimension(s) over which to apply bandpass filter. By default gradient is applied over the `time` dimension.

   .. seealso::
       :py:func:`scipy.signal.butter <scipy:scipy.signal.butter>`, :py:func:`scipy.signal.sosfiltfilt <scipy:scipy.signal.sosfiltfilt>`


.. py:function:: calc_butter_lowpass(data: xarray.DataArray, sampling_frequency: int, period: int, N=3, dim='time') -> xarray.DataArray

   Butterworth lowpass filter.

   Parameters
   ----------
   - data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The array of data to be filtered.
   - sampling_frequency: int.
       Data sampling frequency. If it is daily data with only one time level record per day, 
       then the parameter value is 1; If a day contains four time level records, the parameter value is 4.
   - period: list.
       The low-pass filtering time period, above which the signal (low frequency signal) will pass. 
       If you are getting a 10-day low-pass filter, the value of this parameter is `10`. 
       Note that the units of this parameter should be consistent with `sampling_frequency`.
   - N: int.
       The order of the filter. Default is 3.
   - dim: str.
       Dimension(s) over which to apply lowpass filter. By default gradient is applied over the `time` dimension.

   .. seealso::
       :py:func:`scipy.signal.butter <scipy:scipy.signal.butter>`, :py:func:`scipy.signal.sosfiltfilt <scipy:scipy.signal.sosfiltfilt>`


.. py:function:: calc_butter_highpass(data: xarray.DataArray, sampling_frequency: int, period: int, N=3, dim='time') -> xarray.DataArray

   Butterworth highpass filter.

   Parameters
   ----------
   - data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The array of data to be filtered.
   - sampling_frequency: int.
       Data sampling frequency. If it is daily data with only one time level record per day, 
       then the parameter value is 1; If a day contains four time level records, the parameter value is 4.
   - period: list.
       The high-pass filtering time period below which the signal (high-frequency signal) will pass. 
       If you are obtaining a 10-day high-pass filter, the value of this parameter is `10`. 
       Note that the units of this parameter should be consistent with `sampling_frequency`.
   - N: int.
       The order of the filter. Default is 3.
   - dim: str.
       Dimension(s) over which to apply highpass filter. By default gradient is applied over the `time` dimension.

   .. seealso::
       :py:func:`scipy.signal.butter <scipy:scipy.signal.butter>`, :py:func:`scipy.signal.sosfiltfilt <scipy:scipy.signal.sosfiltfilt>`


