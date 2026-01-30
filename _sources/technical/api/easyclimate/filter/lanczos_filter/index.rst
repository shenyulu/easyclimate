easyclimate.filter.lanczos_filter
=================================

.. py:module:: easyclimate.filter.lanczos_filter

.. autoapi-nested-parse::

   Lanczos filter



Functions
---------

.. autoapisummary::

   easyclimate.filter.lanczos_filter.calc_lanczos_lowpass
   easyclimate.filter.lanczos_filter.calc_lanczos_highpass
   easyclimate.filter.lanczos_filter.calc_lanczos_bandpass


Module Contents
---------------

.. py:function:: calc_lanczos_lowpass(data: xarray.DataArray | xarray.Dataset, window_length: int, period: int, dim: str = 'time', method: Literal['rolling', 'convolve'] = 'rolling', drop_edge: bool = True) -> xarray.DataArray

   Lanczos lowpass filter

   data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
       The array of data to be filtered.
   window_length: :py:class:`int <int>`.
       Slide the size of the window.
   period: :py:class:`float <float>`.
       The low-pass filtering time period, above which the signal (low frequency signal) will pass.
       If you are getting a 10-day low-pass filter, the value of this parameter is `10`.
   dim: :py:class:`str <str>`.
       Dimension(s) over which to apply lowpass filter. By default gradient is applied over the `time` dimension.
   method: :py:class:`str <str>`, default: `rolling`.
       Filter method. Optional values are `rolling`, `convolve`.

   .. seealso::
       - https://github.com/liv0505/Lanczos-Filter/tree/master
       - https://scitools-iris.readthedocs.io/en/stable/generated/gallery/general/plot_SOI_filtering.html
       - `Duchon, C. E. (1979). Lanczos Filtering in One and Two Dimensions. Journal of Applied Meteorology and Climatology, 18(8), 1016-1022. <https://journals.ametsoc.org/view/journals/apme/18/8/1520-0450_1979_018_1016_lfioat_2_0_co_2.xml>`__


.. py:function:: calc_lanczos_highpass(data: xarray.DataArray | xarray.Dataset, window_length: int, period: int, dim: str = 'time', method: Literal['rolling', 'convolve'] = 'rolling', drop_edge: bool = True) -> xarray.DataArray

   Lanczos highpass filter

   data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
       The array of data to be filtered.
   window_length: :py:class:`int <int>`.
       Slide the size of the window.
   period: :py:class:`int <int>`.
       The high-pass filtering time period below which the signal (high-frequency signal) will pass.
       If you are obtaining a 10-day high-pass filter, the value of this parameter is `10`.
   dim: :py:class:`str <str>`.
       Dimension(s) over which to apply highpass filter. By default gradient is applied over the `time` dimension.
   method: :py:class:`str <str>`, default: `rolling`.
       Filter method. Optional values are `rolling`, `convolve`.

   .. seealso::
       - https://github.com/liv0505/Lanczos-Filter/tree/master
       - https://scitools-iris.readthedocs.io/en/stable/generated/gallery/general/plot_SOI_filtering.html
       - `Duchon, C. E. (1979). Lanczos Filtering in One and Two Dimensions. Journal of Applied Meteorology and Climatology, 18(8), 1016-1022. <https://journals.ametsoc.org/view/journals/apme/18/8/1520-0450_1979_018_1016_lfioat_2_0_co_2.xml>`__


.. py:function:: calc_lanczos_bandpass(data: xarray.DataArray | xarray.Dataset, window_length: int, period: list[int], dim: str = 'time', method: Literal['rolling', 'convolve', 'fft'] = 'rolling', drop_edge: bool = True) -> xarray.DataArray

   Lanczos bandpass filter

   data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
       The array of data to be filtered.
   window_length: :py:class:`int <int>`.
       Slide the size of the window.
   period: :py:class:`list[int]`.
       The time period interval of the bandpass filter to be acquired.
       If we are obtaining a 3-10 day bandpass filter, the value of this parameter is `[3, 10]`.
   dim: :py:class:`str <str>`.
       Dimension(s) over which to apply bandpass filter. By default gradient is applied over the `time` dimension.
   method: :py:class:`str <str>`, default: `rolling`.
       Filter method. Optional values are `rolling`, `convolve`.

   .. seealso::
       - https://github.com/liv0505/Lanczos-Filter/tree/master
       - https://scitools-iris.readthedocs.io/en/stable/generated/gallery/general/plot_SOI_filtering.html
       - `Duchon, C. E. (1979). Lanczos Filtering in One and Two Dimensions. Journal of Applied Meteorology and Climatology, 18(8), 1016-1022. <https://journals.ametsoc.org/view/journals/apme/18/8/1520-0450_1979_018_1016_lfioat_2_0_co_2.xml>`__


