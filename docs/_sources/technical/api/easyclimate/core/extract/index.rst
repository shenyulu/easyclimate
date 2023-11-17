:py:mod:`easyclimate.core.extract`
==================================

.. py:module:: easyclimate.core.extract

.. autoapi-nested-parse::

   Obtain data within a specified time period



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.core.extract.get_specific_years_data
   easyclimate.core.extract.get_specific_months_data
   easyclimate.core.extract.get_specific_days_data
   easyclimate.core.extract.get_specific_hours_data
   easyclimate.core.extract.get_specific_minutes_data
   easyclimate.core.extract.get_specific_seconds_data
   easyclimate.core.extract.get_specific_microseconds_data
   easyclimate.core.extract.get_specific_nanoseconds_data
   easyclimate.core.extract.get_specific_dayofweek_data
   easyclimate.core.extract.get_yearmean_for_specific_months_data
   easyclimate.core.extract.get_year_exceed_index_upper_bound
   easyclimate.core.extract.get_year_exceed_index_lower_bound



.. py:function:: get_specific_years_data(data_input: easyclimate.core.yearstat.xr.DataArray, year_array: easyclimate.core.yearstat.np.array, dim='time') -> easyclimate.core.yearstat.xr.DataArray

   Slicing and extracting the part of the data containing the specified year based on an array of given integer years.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   year_array: numpy.array containing int data-type objects
       Year(s) to be extracted.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_months_data(data_input: easyclimate.core.yearstat.xr.DataArray, month_array: easyclimate.core.yearstat.np.array, dim='time') -> easyclimate.core.yearstat.xr.DataArray

   Slicing and extracting the part of the data containing the specified year based on an array of given integer months.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   month_array: numpy.array containing int data-type objects
       Month(s) to be extracted.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_days_data(data_input: easyclimate.core.yearstat.xr.DataArray, day_array: easyclimate.core.yearstat.np.array, dim='time') -> easyclimate.core.yearstat.xr.DataArray

   Slicing and extracting the part of the data containing the specified year based on an array of given integer days.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   day_array: numpy.array containing int data-type objects
       Days(s) to be extracted.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_hours_data(data_input: easyclimate.core.yearstat.xr.DataArray, hour_array: easyclimate.core.yearstat.np.array, dim='time') -> easyclimate.core.yearstat.xr.DataArray

   Slicing and extracting the part of the data containing the specified year based on an array of given integer hours.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   hour_array: numpy.array containing int data-type objects
       Hour(s) to be extracted.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_minutes_data(data_input: easyclimate.core.yearstat.xr.DataArray, minute_array: easyclimate.core.yearstat.np.array, dim='time') -> easyclimate.core.yearstat.xr.DataArray

   Slicing and extracting the part of the data containing the specified year based on an array of given integer minutes.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   minute_array: numpy.array containing int data-type objects
       Minute(s) to be extracted.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_seconds_data(data_input: easyclimate.core.yearstat.xr.DataArray, second_array: easyclimate.core.yearstat.np.array, dim='time') -> easyclimate.core.yearstat.xr.DataArray

   Slicing and extracting the part of the data containing the specified year based on an array of given integer seconds.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   second_array: numpy.array containing int data-type objects
       Second(s) to be extracted.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_microseconds_data(data_input: easyclimate.core.yearstat.xr.DataArray, microsecond_array: easyclimate.core.yearstat.np.array, dim='time') -> easyclimate.core.yearstat.xr.DataArray

   Slicing and extracting the part of the data containing the specified year based on an array of given integer microseconds.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   microsecond_array: numpy.array containing int data-type objects
       Microsecond(s) to be extracted.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_nanoseconds_data(data_input: easyclimate.core.yearstat.xr.DataArray, nanosecond_array: easyclimate.core.yearstat.np.array, dim='time') -> easyclimate.core.yearstat.xr.DataArray

   Slicing and extracting the part of the data containing the specified year based on an array of given integer nanoseconds.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   nanosecond_array: numpy.array containing int data-type objects
       Nanosecond(s) to be extracted.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_dayofweek_data(data_input: easyclimate.core.yearstat.xr.DataArray, dayofweek_array: easyclimate.core.yearstat.np.array, dim='time') -> easyclimate.core.yearstat.xr.DataArray

   Slicing and extracting the part of the data containing the specified year based on an array of given integer dayofweek.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   dayofweek_array: numpy.array containing int data-type objects
       The days of the week to be extracted. 

       The integer numbers correspond to the days of the week as follows.

   +-------------------+-------------------+
   | Day of the week   | Integer numbers   |
   +===================+===================+
   |      Monday       |         0         |
   +-------------------+-------------------+
   |      Tuesday      |         1         |
   +-------------------+-------------------+
   |      Wednesday    |         2         |
   +-------------------+-------------------+
   |      Thursday     |         3         |
   +-------------------+-------------------+
   |      Friday       |         4         |
   +-------------------+-------------------+
   |      Saturday     |         5         |
   +-------------------+-------------------+
   |      Sunday       |         6         |
   +-------------------+-------------------+

   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_yearmean_for_specific_months_data(data_input: easyclimate.core.yearstat.xr.DataArray, month_array: easyclimate.core.yearstat.np.array, dim='time', kwargs=None) -> easyclimate.core.yearstat.xr.DataArray

   Get the annual average of certain months.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   month_array: numpy.array containing int data-type objects
       Month(s) to be extracted.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data. 
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_year_exceed_index_upper_bound(data_input: easyclimate.core.yearstat.xr.DataArray, thresh: float, time_dim: str = 'time') -> easyclimate.core.yearstat.np.array

   Extract the years under the specified threshold (upper bound) in the annual average index (one-dimensional data with only a `time` dimension).

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The one-dimensional data with only a `time` dimension.
   thresh: :py:class:`float<python.float>`.
       The threshold value.
   time_dim: :py:class:`str<python.str>`.
       The time coordinate dimension name.

   Returns
   -------
   :py:class:`numpy.array <numpy:numpy.array>`.


.. py:function:: get_year_exceed_index_lower_bound(data_input: easyclimate.core.yearstat.xr.DataArray, thresh: float, time_dim: str = 'time') -> easyclimate.core.yearstat.np.array

   Extract the years under the specified threshold (lower bound) in the annual average index (one-dimensional data with only a `time` dimension).

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The one-dimensional data with only a `time` dimension.
   thresh: :py:class:`float<python.float>`.
       The threshold value.
   time_dim: :py:class:`str<python.str>`.
       The time coordinate dimension name.

   Returns
   -------
   :py:class:`numpy.array <numpy:numpy.array>`.


