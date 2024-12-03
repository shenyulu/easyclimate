easyclimate.filter
==================

.. py:module:: easyclimate.filter


Submodules
----------

.. toctree::
   :maxdepth: 1

   /technical/api/easyclimate/filter/butter_filter/index
   /technical/api/easyclimate/filter/smooth/index
   /technical/api/easyclimate/filter/wavelet/index
   /technical/api/easyclimate/filter/waveletFunctions/index


Functions
---------

.. autoapisummary::

   easyclimate.filter.calc_butter_bandpass
   easyclimate.filter.calc_butter_lowpass
   easyclimate.filter.calc_butter_highpass
   easyclimate.filter.calc_spatial_smooth_gaussian
   easyclimate.filter.calc_spatial_smooth_rectangular
   easyclimate.filter.calc_spatial_smooth_5or9_point
   easyclimate.filter.calc_forward_smooth
   easyclimate.filter.calc_reverse_smooth
   easyclimate.filter.calc_spatial_smooth_circular
   easyclimate.filter.calc_spatial_smooth_window
   easyclimate.filter.calc_timeseries_wavelet_transform
   easyclimate.filter.draw_global_wavelet_spectrum
   easyclimate.filter.draw_wavelet_transform


Package Contents
----------------

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


.. py:function:: calc_spatial_smooth_gaussian(data: xarray.DataArray | xarray.Dataset, n: int = 3, lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.DataArray | xarray.Dataset

   Filter with normal distribution of weights.

   Parameters
   ----------
   data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
       Some n-dimensional latitude-longitude dataset.
   n: :py:class:`int <int>`.
       Degree of filtering.

   Returns
   -------
   The filtered latitude-longitude dataset (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso:

       :py:func:`metpy.calc.smooth_gaussian <metpy:metpy.calc.smooth_gaussian>`


.. py:function:: calc_spatial_smooth_rectangular(data: xarray.DataArray | xarray.Dataset, rectangle_shapes: int | Sequence[int] = 3, times: int = 1, lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.DataArray | xarray.Dataset

   Filter with a rectangular window smoother.

   Parameters
   ----------
   data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
       Some n-dimensional latitude-longitude dataset.
   rectangle_shapes: :py:class:`int <int>` or :py:class:`Sequence[int]`, default `3`.
       Shape of rectangle along the trailing dimension(s) of the data.
   times: :py:class:`int <int>` or , default `1`.
       The number of times to apply the filter to the data.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The filtered latitude-longitude dataset (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso:

       :py:func:`metpy.calc.smooth_rectangular <metpy:metpy.calc.smooth_rectangular>`


.. py:function:: calc_spatial_smooth_5or9_point(data: xarray.DataArray | xarray.Dataset, n: {5, 9} = 5, times: int = 1, lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.DataArray | xarray.Dataset

   Filter with an 5-point or 9-point smoother.

   Parameters
   ----------
   data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
       Some n-dimensional latitude-longitude dataset.
   n: {5, 9}, default `5`.
       The number of points to use in smoothing, only valid inputs are 5 and 9.
   times: :py:class:`int <int>` or , default `1`.
       The number of times to apply the filter to the data.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The filtered latitude-longitude dataset (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso:

       :py:func:`metpy.calc.smooth_n_point <metpy:metpy.calc.smooth_n_point>`


.. py:function:: calc_forward_smooth(data: xarray.DataArray | xarray.Dataset, n: int = 5, S: float = 0.5, times: int = 1, normalize_weights=False, lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.DataArray | xarray.Dataset

   Filter with the forward smooth.

   Parameters
   ----------
   data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
       Some n-dimensional latitude-longitude dataset.
   n: {5, 9}, default `5`.
       The number of points to use in smoothing, only valid inputs are 5 and 9.
   S: :py:class:`float <float>`, default `0.5`.
       The smoothing factor.
   times: :py:class:`int <int>`, default `1`.
       The number of times to apply the filter to the data.
   normalize_weights: :py:class:`bool <bool>`, default `False`.
       If `True`, divide the values in window by the sum of all values in the window to obtain the normalized smoothing weights. If `False`, use supplied values directly as the weights.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The filtered latitude-longitude dataset (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. note::

       For any discrete variable, its value on the :math:`x`-axis can be described as :math:`F_i = F(x_i)`, where :math:`x_i = i \Delta x, i = 0, \pm 1, \pm 2`.

       Define a one-dimensional smoothing operator

       .. math::

           \overset{\sim}{F}_{i}^x = (1-S)F_i + \frac{S}{2} (F_{i+1} + F_{i-1})

       Here, :math:`S` is the spatial smoothing coefficient. Since the above formula only involves the values of three grid data points (:math:`F_i, F_{i+1}, F_{i-1}`), it is also referred to as a three-point smoothing.

       Correspondingly, for a two-dimensional variable :math:`F_{i,j}`, smoothing is performed in the :math:`x` direction as follows

       .. math::

           \overset{\sim}{F}_{i,j}^x = F_{i,j} + \frac{S}{2} (F_{i+1,j} + F_{i-1,j}-2F_{i,j})

       And, smoothing is done in the :math:`y` direction separately

       .. math::

           \overset{\sim}{F}_{i,j}^y = F_{i,j} + \frac{S}{2} (F_{i,j+1} + F_{i,j-1}-2F_{i,j})

       Taking the average of both results:

       .. math::

           \overset{\sim}{F}_{i,j}^{x,y} = \frac{\overset{\sim}{F}_{i,j}^x + \overset{\sim}{F}_{i,j}^y}{2} = F_{i,j} + \frac{S}{4}(F_{i+1,j}+F_{i,j+1}+F_{i-1,j}+F_{i,j-1}-4F_{i,j})

       At this point, the following "window" matrix :math:`\mathbf{W}` can be obtained

       .. math::

           \mathbf{W} =
           \left(
           \begin{array}{ccc}
               F_{i-1,j+1} &  F_{i,j+1} & F_{i+1,j+1} \\
               F_{i-1,j} &  F_{i,j} & F_{i+1,j} \\
               F_{i-1,j-1} &  F_{i,j-1} & F_{i+1,j-1}
           \end{array}
           \right)
           =
           \left(
           \begin{array}{ccc}
               0 &  \frac{S}{4} & 0 \\
               \frac{S}{4} &  (1-S) & \frac{S}{4} \\
               0 &  \frac{S}{4} & 0
           \end{array}
           \right)

       When :math:`S=0.5`, we have

       .. math::

           \mathbf{W} =
           \left(
           \begin{array}{ccc}
               0 &  0.125 & 0 \\
               0.125 &  0.5 & 0.125 \\
               0 &  0.125 & 0
           \end{array}
           \right)

       This is the "window" matrix used by the function when :math:`n=5`.

       For the variable :math:`F_{i,j}`, applying a three-point smoothing in the :math:`x`-axis direction followed by a three-point smoothing in the :math:`y`-axis direction results in the following nine-point smoothing scheme

       .. math::

           \overset{\sim}{F}_{i,j}^{x,y} = {F}_{i,j} + \frac{S}{2}(1-S)(F_{i+1,j}+F_{i,j+1}+F_{i-1,j}+F_{i,j-1}-4F_{i,j})+\frac{S^2}{4}(F_{i+1,j+1}+F_{i-1,j+1}+F_{i-1,j-1}+F_{i+1,j-1}-4F_{i,j})

       At this point, the corresponding "window" matrix :math:`\mathbf{W}` is given by

       .. math::

           \mathbf{W} =
           \left(
           \begin{array}{ccc}
               F_{i-1,j+1} &  F_{i,j+1} & F_{i+1,j+1} \\
               F_{i-1,j} &  F_{i,j} & F_{i+1,j} \\
               F_{i-1,j-1} &  F_{i,j-1} & F_{i+1,j-1}
           \end{array}
           \right)
           =
           \left(
           \begin{array}{ccc}
               \frac{S^2}{4} &  \frac{S}{2}(1-S) & \frac{S^2}{4} \\
               \frac{S}{2}(1-S) &  1-2S(1-S) & \frac{S}{2}(1-S) \\
               \frac{S^2}{4} &  \frac{S}{2}(1-S) & \frac{S^2}{4}
           \end{array}
           \right)

       When :math:`S=0.5`, we have

       .. math::

           \mathbf{W} =
           \left(
           \begin{array}{ccc}
               0.0625 &  0.125 & 0.0625 \\
               0.125 &  0.5 & 0.125 \\
               0.0625 &  0.125 & 0.0625
           \end{array}
           \right)

       This is the "window" matrix used by the function when :math:`n=9`.

   .. seealso::
       :py:func:`metpy.calc.smooth_n_point <metpy:metpy.calc.smooth_n_point>`


.. py:function:: calc_reverse_smooth(data: xarray.DataArray | xarray.Dataset, n: int = 5, S: float = -0.5, times: int = 1, normalize_weights=False, lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.DataArray | xarray.Dataset

   Filter with the reverse smooth.

   Parameters
   ----------
   data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
       Some n-dimensional latitude-longitude dataset.
   n: {5, 9}, default `5`.
       The number of points to use in smoothing, only valid inputs are 5 and 9.
   S: :py:class:`float <float>`, default `-0.5`.
       The smoothing factor.
   times: :py:class:`int <int>`, default `1`.
       The number of times to apply the filter to the data.
   normalize_weights: :py:class:`bool <bool>`, default `False`.
       If `True`, divide the values in window by the sum of all values in the window to obtain the normalized smoothing weights. If `False`, use supplied values directly as the weights.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The filtered latitude-longitude dataset (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. note::

       For any discrete variable, its value on the :math:`x`-axis can be described as :math:`F_i = F(x_i)`, where :math:`x_i = i \Delta x, i = 0, \pm 1, \pm 2`.

       Define a one-dimensional smoothing operator

       .. math::

           \overset{\sim}{F}_{i}^x = (1-S)F_i + \frac{S}{2} (F_{i+1} + F_{i-1})

       Here, :math:`S` is the spatial smoothing coefficient. Since the above formula only involves the values of three grid data points (:math:`F_i, F_{i+1}, F_{i-1}`), it is also referred to as a three-point smoothing.

       Correspondingly, for a two-dimensional variable :math:`F_{i,j}`, smoothing is performed in the :math:`x` direction as follows

       .. math::

           \overset{\sim}{F}_{i,j}^x = F_{i,j} + \frac{S}{2} (F_{i+1,j} + F_{i-1,j}-2F_{i,j})

       And, smoothing is done in the :math:`y` direction separately

       .. math::

           \overset{\sim}{F}_{i,j}^y = F_{i,j} + \frac{S}{2} (F_{i,j+1} + F_{i,j-1}-2F_{i,j})

       Taking the average of both results:

       .. math::

           \overset{\sim}{F}_{i,j}^{x,y} = \frac{\overset{\sim}{F}_{i,j}^x + \overset{\sim}{F}_{i,j}^y}{2} = F_{i,j} + \frac{S}{4}(F_{i+1,j}+F_{i,j+1}+F_{i-1,j}+F_{i,j-1}-4F_{i,j})

       At this point, the following "window" matrix :math:`\mathbf{W}` can be obtained

       .. math::

           \mathbf{W} =
           \left(
           \begin{array}{ccc}
               F_{i-1,j+1} &  F_{i,j+1} & F_{i+1,j+1} \\
               F_{i-1,j} &  F_{i,j} & F_{i+1,j} \\
               F_{i-1,j-1} &  F_{i,j-1} & F_{i+1,j-1}
           \end{array}
           \right)
           =
           \left(
           \begin{array}{ccc}
               0 &  \frac{S}{4} & 0 \\
               \frac{S}{4} &  (1-S) & \frac{S}{4} \\
               0 &  \frac{S}{4} & 0
           \end{array}
           \right)

       When :math:`S=0.5`, we have

       .. math::

           \mathbf{W} =
           \left(
           \begin{array}{ccc}
               0 &  0.125 & 0 \\
               0.125 &  0.5 & 0.125 \\
               0 &  0.125 & 0
           \end{array}
           \right)

       This is the "window" matrix used by the function when :math:`n=5`.

       For the variable :math:`F_{i,j}`, applying a three-point smoothing in the :math:`x`-axis direction followed by a three-point smoothing in the :math:`y`-axis direction results in the following nine-point smoothing scheme

       .. math::

           \overset{\sim}{F}_{i,j}^{x,y} = {F}_{i,j} + \frac{S}{2}(1-S)(F_{i+1,j}+F_{i,j+1}+F_{i-1,j}+F_{i,j-1}-4F_{i,j})+\frac{S^2}{4}(F_{i+1,j+1}+F_{i-1,j+1}+F_{i-1,j-1}+F_{i+1,j-1}-4F_{i,j})

       At this point, the corresponding "window" matrix :math:`\mathbf{W}` is given by

       .. math::

           \mathbf{W} =
           \left(
           \begin{array}{ccc}
               F_{i-1,j+1} &  F_{i,j+1} & F_{i+1,j+1} \\
               F_{i-1,j} &  F_{i,j} & F_{i+1,j} \\
               F_{i-1,j-1} &  F_{i,j-1} & F_{i+1,j-1}
           \end{array}
           \right)
           =
           \left(
           \begin{array}{ccc}
               \frac{S^2}{4} &  \frac{S}{2}(1-S) & \frac{S^2}{4} \\
               \frac{S}{2}(1-S) &  1-2S(1-S) & \frac{S}{2}(1-S) \\
               \frac{S^2}{4} &  \frac{S}{2}(1-S) & \frac{S^2}{4}
           \end{array}
           \right)

       When :math:`S=0.5`, we have

       .. math::

           \mathbf{W} =
           \left(
           \begin{array}{ccc}
               0.0625 &  0.125 & 0.0625 \\
               0.125 &  0.5 & 0.125 \\
               0.0625 &  0.125 & 0.0625
           \end{array}
           \right)

       This is the "window" matrix used by the function when :math:`n=9`.

   .. seealso::
       :py:func:`metpy.calc.smooth_n_point <metpy:metpy.calc.smooth_n_point>`


.. py:function:: calc_spatial_smooth_circular(data: xarray.DataArray | xarray.Dataset, radius: int, times: int = 1, lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.DataArray | xarray.Dataset

   Filter with a circular window smoother.

   Parameters
   ----------
   data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
       Some n-dimensional latitude-longitude dataset.
   radius: :py:class:`int <int>`.
       Radius of the circular smoothing window. The “diameter” of the circle (width of smoothing window) is :math:`2 \cdot \mathrm{radius} + 1` to provide a smoothing window with odd shape.
   times: :py:class:`int <int>`, default `1`.
       The number of times to apply the filter to the data.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The filtered latitude-longitude dataset (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso:

       :py:func:`metpy.calc.smooth_circular <metpy:metpy.calc.smooth_circular>`


.. py:function:: calc_spatial_smooth_window(data: xarray.DataArray | xarray.Dataset, window: numpy.ndarray | xarray.DataArray, times: int = 1, normalize_weights: bool = False, lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.DataArray | xarray.Dataset

   Filter with an arbitrary window smoother.

   Parameters
   ----------
   data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
       Some n-dimensional latitude-longitude dataset.
   window : :py:class:`numpy.ndarray <numpy.ndarray>` or :py:class:`xarray.DataArray<xarray.DataArray>`.
       Window to use in smoothing. Can have dimension less than or equal to :math:`N`. If dimension less than :math:`N`, the scalar grid will be smoothed along its trailing dimensions. Shape along each dimension must be odd.
   times: :py:class:`int <int>`, default `1`.
       The number of times to apply the filter to the data.
   normalize_weights: :py:class:`bool <bool>`, default `False`.
       If `True`, divide the values in window by the sum of all values in the window to obtain the normalized smoothing weights. If `False`, use supplied values directly as the weights.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The filtered latitude-longitude dataset (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso:

       :py:func:`metpy.calc.smooth_window <metpy:metpy.calc.smooth_window>`


.. py:function:: calc_timeseries_wavelet_transform(timeseries_data: xarray.DataArray, dt: float, time_dim: str = 'time', pad: float = 1, dj: float = 0.25, j1_n_div_dj: int = 7, s0=None, lag1: float = None, mother: str = 'morlet', mother_param=None, sigtest_wavelet: str = 'regular chi-square test', sigtest_global: str = 'time-average test', significance_level: float = 0.95) -> xarray.Dataset

   Wavelet transform parameters calculation.

   Parameters
   ----------
   timeseries_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Timeseries data.
   dt: :py:class:`float<float>`.
       Amount of time between each timeseries data value, i.e. the sampling time.
   time_dim: :py:class:`str <str>`.
       The time coordinate dimension name.
   pad: :py:class:`float<float>`, default: `1`.
       if set to 1 (default is 0), pad time series with zeroes to get N up to the next higher power of 2.
       This prevents wraparound from the end of the time series to the beginning, and also
       speeds up the FFT's used to do the wavelet transform. This will not eliminate all edge effects (see COI below).
   dj: :py:class:`float<float>`, default: `0.25`.
       The spacing between discrete scales. A smaller `dj` will give better scale resolution, but be slower to plot.
   j1_n_div_dj: :py:class:`int<int>`, default: `7`.
       Do `j1_n_div_dj` powers-of-two with `dj` sub-octaves each.

       .. math::
           j1 = \mathrm{j1\underline{ }n\underline{ }div\underline{ }dj} / dj

       This can adjust the size of the period.

   s0: :py:class:`float<float>`, default: :math:`2 \times \mathrm{d}t`.
       The smallest scale of the wavelet.
   lag1: :py:class:`float<float>`, default: `None`.
       Lag-1 autocorrelation for red noise background. The value is generated by :py:func:`statsmodels.api.tsa.acf <statsmodels.api.tsa.acf>`.
   mother: :py:class:`str<str>`, {'morlet', 'paul', 'dog'}, default: 'morlet'.
       The mother wavelet function.

       .. tab:: Morlet

           - Name: Morlet (:math:`\omega_0` = frequency)
           - :math:`\psi_0(\eta)`: :math:`\pi^{-1/4} e^{i \omega_{0} \eta} e^{-\eta^2/2}`.

       .. tab:: Paul

           - Name: Paul (:math:`m` = order)
           - :math:`\psi_0(\eta)`: :math:`\frac{2^m i^m m!}{\sqrt{\pi (2m)!}} (1-i \eta)^{-(m+1)}`.

       .. tab:: DOG

           - Name: DOG (:math:`m` = derivative)
           - :math:`\psi_0(\eta)`: :math:`\frac{(-1)^{m+1}}{\sqrt{\Gamma (m+\frac{1}{2})}} \frac{d^m}{d \eta^m} (e^{-\eta^2 /2})`.


   mother_param: :py:class:`float<float>`.
       The mother wavelet parameter.
           - For 'morlet' this is :math:`k_0` (wavenumber), default is 6.
           - For 'paul' this is :math:`m` (order), default is 4.
           - For 'dog' this is :math:`m` (m-th derivative), default is 2.

   sigtest_wavelet: {'regular chi-square test', 'time-average test', 'scale-average test'}, default: `'regular chi-square test'`.
       The type of significance test.

       1. Regular chi-square test
       i.e. Eqn (18) from Torrence & Compo.

       .. math::

           \frac{\left|W_n(s)\right|^2}{\sigma^2}\Longrightarrow\frac{1}{2} P_k\chi_2^2

       2. The "time-average" test, i.e. Eqn (23).

       .. math::

           \nu=2\sqrt{1+\left(\frac{n_a\delta t}{\gamma s}\right)^2}

       In this case, DOF should be set to NA, the number
       of local wavelet spectra that were averaged together.
       For the Global Wavelet Spectrum, this would be NA=N,
       where N is the number of points in your time series.

       3. The "scale-average" test, i.e. Eqns (25)-(28).

       .. math::

           \overline{P}=S_{\mathrm{avg}}\sum_{j=j_1}^{j_2}\frac{P_j}{S_j}, \ \mathrm{where} \ S_{\mathrm{avg}}=\left(\sum_{j=j_1}^{j_2}\frac1{s_j}\right)^{-1}, \frac{C_\delta S_\mathrm{avg}}{\delta j\delta t\sigma^2}\overline{W}_n^2\Rightarrow\overline{P}\frac{\chi_\nu^2}\nu, \nu=\frac{2n_aS_{\mathrm{avg}}}{S_{\mathrm{mid}}}\sqrt{1+\left(\frac{n_a\delta j}{\delta j_0}\right)^2}.

       In this case, DOF should be set to a
       two-element vector [S1,S2], which gives the scale
       range that was averaged together.
       e.g. if one scale-averaged scales between 2 and 8,
       then DOF=[2,8].

   sigtest_global: {'regular chi-square test', 'time-average test', 'scale-average test'}, default: `'time-average test'`.
       See also the description of `sigtest_wavelet`.

   significance_level: :py:class:`float<float>`, default: `0.95`.
       Significance level to use.

   Returns
   -------
   Timeseries wavelet transform result (:py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       - https://github.com/regeirk/pycwt, https://pycwt.readthedocs.io/en/latest/index.html
       - http://nicolasfauchereau.github.io/climatecode/posts/wavelet-analysis-in-python/
       - https://bbs.06climate.com/forum.php?mod=viewthread&tid=95016
       - https://blog.csdn.net/weixin_43304836/article/details/119752767

   Reference
   --------------
   - Torrence, C., & Compo, G. P. (1998). A Practical Guide to Wavelet Analysis. Bulletin of the American Meteorological Society, 79(1), 61-78. https://doi.org/10.1175/1520-0477(1998)079<0061:APGTWA>2.0.CO;2
   - Torrence, C., & Webster, P. J. (1999). Interdecadal Changes in the ENSO–Monsoon System. Journal of Climate, 12(8), 2679-2690. https://doi.org/10.1175/1520-0442(1999)012<2679:ICITEM>2.0.CO;2
   - Grinsted, A., Moore, J. C., and Jevrejeva, S.: Application of the cross wavelet transform and wavelet coherence to geophysical time series, Nonlin. Processes Geophys., 11, 561–566, https://doi.org/10.5194/npg-11-561-2004, 2004.


.. py:function:: draw_global_wavelet_spectrum(timeseries_wavelet_transform_result: xarray.Dataset, ax: matplotlib.axes.Axes = None, global_ws_kwargs: dict = {}, global_signif_kwargs: dict = {'ls': '--'})

   Draw global wavelet spectrum

   Parameters
   ----------
   timeseries_wavelet_transform_result: :py:class:`xarray.Dataset<xarray.Dataset>`.
       Timeseries wavelet transform result.
   ax: :py:class:`matplotlib.axes.Axes`
       The axes to which the boundary will be applied.
   **global_ws_kwargs, :py:class:`dict <dict>`, optional:
       Additional keyword arguments to :py:func:`xarray.DataArray.plot.line<xarray.DataArray.plot.line>` for ploting `global_ws`.
   **global_signif_kwargs, :py:class:`dict <dict>`, optional, default {'ls': '--'}:
       Additional keyword arguments to :py:func:`xarray.DataArray.plot.line<xarray.DataArray.plot.line>` for ploting `global_signif`.


.. py:function:: draw_wavelet_transform(timeseries_wavelet_transform_result: xarray.Dataset, ax: matplotlib.axes.Axes = None, power_kwargs: dict = {'levels': [0, 0.5, 1, 2, 4, 999], 'colors': ['white', 'bisque', 'orange', 'orangered', 'darkred']}, sig_kwargs: dict = {'levels': [-99, 1], 'colors': 'k'}, coi_kwargs: dict = {'color': 'k'}, fill_between_kwargs: dict = {'facecolor': 'none', 'edgecolor': '#00000040', 'hatch': 'x'})

   Draw wavelet transform

   Parameters
   ----------
   timeseries_wavelet_transform_result: :py:class:`xarray.Dataset<xarray.Dataset>`.
       Timeseries wavelet transform result.
   ax : :py:class:`matplotlib.axes.Axes`
       The axes to which the boundary will be applied.
   **power_kwargs, optional, :py:class:`dict <dict>`, default {'levels': [0, 0.5, 1, 2, 4, 999], 'colors': ['white', 'bisque', 'orange', 'orangered', 'darkred']}:
       Additional keyword arguments to :py:func:`xarray.DataArray.plot.contourf<xarray.DataArray.plot.contourf>` for ploting `power`.
   **sig_kwargs, optional, :py:class:`dict <dict>`, default {'levels': [-99, 1], 'colors': 'k'}:
       Additional keyword arguments to :py:func:`xarray.DataArray.plot.contourf<xarray.DataArray.plot.contourf>` for ploting `sig`.
   **coi_kwargs, optional, :py:class:`dict <dict>`, default {'color': 'k'}:
       Additional keyword arguments to :py:func:`xarray.DataArray.plot.contourf<xarray.DataArray.plot.contourf>` for ploting `coi`.
   **fill_between_kwargs, :py:class:`dict <dict>`, optional, default {'facecolor': 'none', 'edgecolor': '#00000040', 'hatch': 'x'}:
       Additional keyword arguments to :py:func:`matplotlib.pyplot.fill_between<matplotlib.pyplot.fill_between>`.


