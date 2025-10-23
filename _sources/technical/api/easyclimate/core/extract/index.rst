easyclimate.core.extract
========================

.. py:module:: easyclimate.core.extract

.. autoapi-nested-parse::

   Obtain data within a specified time period



Functions
---------

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
   easyclimate.core.extract.get_time_exceed_index_upper_bound
   easyclimate.core.extract.get_time_exceed_index_lower_bound


Module Contents
---------------

.. py:function:: get_specific_years_data(data_input: xarray.DataArray | xarray.Dataset, year_array: np.array(int) | List[int], dim: str = 'time') -> xarray.DataArray | xarray.Dataset

   Slicing and extracting the part of the data containing the specified year based on an array of given integer years.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   year_array: :py:class:`list[int]`
       Year(s) to be extracted.
   dim : :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_months_data(data_input: xarray.DataArray | xarray.Dataset, month_array: numpy.array, dim: str = 'time') -> xarray.DataArray | xarray.Dataset

   Slicing and extracting the part of the data containing the specified year based on an array of given integer months.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   month_array: :py:class:`list[int]`
       Month(s) to be extracted.
   dim : :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_basic_statistical_analysis.py


.. py:function:: get_specific_days_data(data_input: xarray.DataArray | xarray.Dataset, day_array: np.array(int) | List[int], dim: str = 'time') -> xarray.DataArray | xarray.Dataset

   Slicing and extracting the part of the data containing the specified year based on an array of given integer days.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   day_array: :py:class:`list[int]`
       Days(s) to be extracted.
   dim : :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_hours_data(data_input: xarray.DataArray | xarray.Dataset, hour_array: np.array(int) | List[int], dim: str = 'time') -> xarray.DataArray | xarray.Dataset

   Slicing and extracting the part of the data containing the specified year based on an array of given integer hours.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   hour_array: :py:class:`list[int]`
       Hour(s) to be extracted.
   dim : :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_minutes_data(data_input: xarray.DataArray | xarray.Dataset, minute_array: np.array(int) | List[int], dim: str = 'time') -> xarray.DataArray | xarray.Dataset

   Slicing and extracting the part of the data containing the specified year based on an array of given integer minutes.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   minute_array: :py:class:`list[int]`
       Minute(s) to be extracted.
   dim : :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_seconds_data(data_input: xarray.DataArray | xarray.Dataset, second_array: np.array(int) | List[int], dim: str = 'time') -> xarray.DataArray | xarray.Dataset

   Slicing and extracting the part of the data containing the specified year based on an array of given integer seconds.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   second_array: :py:class:`list[int]`
       Second(s) to be extracted.
   dim : :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_microseconds_data(data_input: xarray.DataArray | xarray.Dataset, microsecond_array: np.array(int) | List[int], dim: str = 'time') -> xarray.DataArray | xarray.Dataset

   Slicing and extracting the part of the data containing the specified year based on an array of given integer microseconds.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   microsecond_array: :py:class:`list[int]`
       Microsecond(s) to be extracted.
   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_nanoseconds_data(data_input: xarray.DataArray | xarray.Dataset, nanosecond_array: np.array(int) | List[int], dim: str = 'time') -> xarray.DataArray | xarray.Dataset

   Slicing and extracting the part of the data containing the specified year based on an array of given integer nanoseconds.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   nanosecond_array: :py:class:`list[int]`
       Nanosecond(s) to be extracted.
   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_specific_dayofweek_data(data_input: xarray.DataArray | xarray.Dataset, dayofweek_array: np.array(int) | List[int], dim: str = 'time') -> xarray.DataArray | xarray.Dataset

   Slicing and extracting the part of the data containing the specified year based on an array of given integer dayofweek.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   dayofweek_array: :py:class:`list[int]`
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

   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_yearmean_for_specific_months_data(data_input: xarray.DataArray | xarray.Dataset, month_array: np.array(int) | List[int], dim: str = 'time', **kwargs) -> xarray.DataArray | xarray.Dataset

   Get the annual average of certain months.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
   month_array: :py:class:`list[int]`
       Month(s) to be extracted.
   dim: :py:class:`str <str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_year_exceed_index_upper_bound(data_input: xarray.DataArray, thresh: float, time_dim: str = 'time') -> numpy.array

   Extract the years under the specified threshold (upper bound) in the annual average index (one-dimensional data with only a `time` dimension).

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The one-dimensional data with only a `time` dimension.
   thresh: :py:class:`float <float>`.
       The threshold value.
   time_dim: :py:class:`str <str>`.
       The time coordinate dimension name.

   Returns
   -------
   :py:class:`numpy.array <numpy:numpy.array>`.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_basic_statistical_analysis.py


.. py:function:: get_year_exceed_index_lower_bound(data_input: xarray.DataArray, thresh: float, time_dim: str = 'time') -> numpy.array

   Extract the years under the specified threshold (lower bound) in the annual average index (one-dimensional data with only a `time` dimension).

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The one-dimensional data with only a `time` dimension.
   thresh: :py:class:`float <float>`.
       The threshold value.
   time_dim: :py:class:`str <str>`.
       The time coordinate dimension name.

   Returns
   -------
   :py:class:`numpy.array <numpy:numpy.array>`.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_basic_statistical_analysis.py


.. py:function:: get_time_exceed_index_upper_bound(data_input: xarray.DataArray, thresh: float, time_dim: str = 'time') -> numpy.array

   Extract the time under the specified threshold (upper bound) in the annual average index (one-dimensional data with only a `time` dimension).

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The one-dimensional data with only a `time` dimension.
   thresh: :py:class:`float <float>`.
       The threshold value.
   time_dim: :py:class:`str <str>`.
       The time coordinate dimension name.

   Returns
   -------
   Time array.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_da_bbo.py


.. py:function:: get_time_exceed_index_lower_bound(data_input: xarray.DataArray, thresh: float, time_dim: str = 'time') -> numpy.array

   Extract the time under the specified threshold (lower bound) in the annual average index (one-dimensional data with only a `time` dimension).

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The one-dimensional data with only a `time` dimension.
   thresh: :py:class:`float <float>`.
       The threshold value.
   time_dim: :py:class:`str <str>`.
       The time coordinate dimension name.

   Returns
   -------
   Time array.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_da_bbo.py


