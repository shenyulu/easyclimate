easyclimate.filter.butter_filter
================================

.. py:module:: easyclimate.filter.butter_filter

.. autoapi-nested-parse::

   Butterworth bandpass filter



Functions
---------

.. autoapisummary::

   easyclimate.filter.butter_filter.calc_butter_bandpass
   easyclimate.filter.butter_filter.calc_butter_lowpass
   easyclimate.filter.butter_filter.calc_butter_highpass


Module Contents
---------------

.. py:function:: calc_butter_bandpass(data: xarray.DataArray | xarray.Dataset, sampling_frequency: int, period: list[int], N: int = 3, dim: str = 'time') -> xarray.DataArray

   Butterworth bandpass filter.

   Parameters
   ----------
   data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
       The array of data to be filtered.
   sampling_frequency: :py:class:`int <int>`.
       Data sampling frequency. If it is daily data with only one time level record per day,
       then the parameter value is 1; If a day contains four time level records, the parameter value is 4.
   period: :py:class:`list[int]`.
       The time period interval of the bandpass filter to be acquired.
       If we are obtaining a 3-10 day bandpass filter, the value of this parameter is `[3, 10]`.
       Note that the units of this parameter should be consistent with `sampling_frequency`.
   N: :py:class:`int <int>`.
       The order of the filter. Default is 3.
   dim: :py:class:`str <str>`..
       Dimension(s) over which to apply bandpass filter. By default gradient is applied over the `time` dimension.

   .. seealso::
       :py:func:`scipy.signal.butter <scipy:scipy.signal.butter>`, :py:func:`scipy.signal.sosfiltfilt <scipy:scipy.signal.sosfiltfilt>`


.. py:function:: calc_butter_lowpass(data: xarray.DataArray | xarray.Dataset, sampling_frequency: int, period: float, N: int = 3, dim: str = 'time') -> xarray.DataArray

   Butterworth lowpass filter.

   Parameters
   ----------
   data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
       The array of data to be filtered.
   sampling_frequency: :py:class:`int <int>`.
       Data sampling frequency. If it is daily data with only one time level record per day,
       then the parameter value is 1; If a day contains four time level records, the parameter value is 4.
   period: :py:class:`float <float>`.
       The low-pass filtering time period, above which the signal (low frequency signal) will pass.
       If you are getting a 10-day low-pass filter, the value of this parameter is `10`.
       Note that the units of this parameter should be consistent with `sampling_frequency`.
   N: :py:class:`int <int>`.
       The order of the filter. Default is 3.
   dim: :py:class:`str <str>`..
       Dimension(s) over which to apply lowpass filter. By default gradient is applied over the `time` dimension.

   .. seealso::
       :py:func:`scipy.signal.butter <scipy:scipy.signal.butter>`, :py:func:`scipy.signal.sosfiltfilt <scipy:scipy.signal.sosfiltfilt>`


.. py:function:: calc_butter_highpass(data: xarray.DataArray | xarray.Dataset, sampling_frequency: int, period: float, N: int = 3, dim: str = 'time') -> xarray.DataArray

   Butterworth highpass filter.

   Parameters
   ----------
   data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
       The array of data to be filtered.
   sampling_frequency: :py:class:`int <int>`.
       Data sampling frequency. If it is daily data with only one time level record per day,
       then the parameter value is 1; If a day contains four time level records, the parameter value is 4.
   period: :py:class:`float <float>`.
       The high-pass filtering time period below which the signal (high-frequency signal) will pass.
       If you are obtaining a 10-day high-pass filter, the value of this parameter is `10`.
       Note that the units of this parameter should be consistent with `sampling_frequency`.
   N: :py:class:`int <int>`.
       The order of the filter. Default is 3.
   dim: :py:class:`str <str>`..
       Dimension(s) over which to apply highpass filter. By default gradient is applied over the `time` dimension.

   .. seealso::
       :py:func:`scipy.signal.butter <scipy:scipy.signal.butter>`, :py:func:`scipy.signal.sosfiltfilt <scipy:scipy.signal.sosfiltfilt>`


