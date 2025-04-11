easyclimate.filter
==================

.. py:module:: easyclimate.filter


Submodules
----------

.. toctree::
   :maxdepth: 1

   /technical/api/easyclimate/filter/barnes_filter/index
   /technical/api/easyclimate/filter/butter_filter/index
   /technical/api/easyclimate/filter/kf_filter/index
   /technical/api/easyclimate/filter/lanczos_filter/index
   /technical/api/easyclimate/filter/redfit/index
   /technical/api/easyclimate/filter/smooth/index
   /technical/api/easyclimate/filter/spatial_pcf/index
   /technical/api/easyclimate/filter/wavelet/index


Functions
---------

.. autoapisummary::

   easyclimate.filter.calc_barnes_lowpass
   easyclimate.filter.calc_barnes_bandpass
   easyclimate.filter.calc_butter_bandpass
   easyclimate.filter.calc_butter_lowpass
   easyclimate.filter.calc_butter_highpass
   easyclimate.filter.calc_lanczos_bandpass
   easyclimate.filter.calc_lanczos_lowpass
   easyclimate.filter.calc_lanczos_highpass
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
   easyclimate.filter.calc_redfit
   easyclimate.filter.calc_redfit_cross
   easyclimate.filter.kf_filter_wheeler_and_kiladis_1999
   easyclimate.filter.kf_filter_lf_wave
   easyclimate.filter.kf_filter_mjo_wave
   easyclimate.filter.kf_filter_er_wave
   easyclimate.filter.kf_filter_kelvin_wave
   easyclimate.filter.kf_filter_mt_wave
   easyclimate.filter.kf_filter_mrg_wave
   easyclimate.filter.kf_filter_td_wave
   easyclimate.filter.filter_2D_spatial_parabolic_cylinder_function


Package Contents
----------------

.. py:function:: calc_barnes_lowpass(data: xarray.DataArray, g: float = 0.3, c: int = 150000, lon_dim='lon', lat_dim='lat', radius_degree=8, print_progress=True) -> xarray.DataArray

   Selecting different parameters **g** and **c**
   will result in different filtering characteristics.

   .. warning::

       **NOT** support ``data`` contains ``np.nan``.


   Parameters
   ----------
   g : :py:class:`float <float>`, generally between (0, 1]
       Constant parameter.
   c : :py:class:`int <int>`
       Constant parameter. When *c* takes a larger value, the filter function converges
       at a larger wavelength, and the response function slowly approaches the maximum value,
       which means that high-frequency fluctuations have been filtered out.
   print_progress : :py:class:`bool <bool>`
       Whether to print the progress bar when executing computation.

   Returns
   -------
   data_vars : :py:class:`xarray.DataArray <xarray.DataArray>`
       Data field after filtering out high-frequency fluctuations

   .. seealso::
       - Maddox, R. A. (1980). An Objective Technique for Separating Macroscale and Mesoscale Features in Meteorological Data. Monthly Weather Review, 108(8), 1108-1121. https://journals.ametsoc.org/view/journals/mwre/108/8/1520-0493_1980_108_1108_aotfsm_2_0_co_2.xml
       - https://github.com/LinOuyang/pybarnes


.. py:function:: calc_barnes_bandpass(data: xarray.DataArray, g1: float = 0.3, g2: float = 0.3, c1: int = 30000, c2: int = 150000, r=1.2, lon_dim='lon', lat_dim='lat', radius_degree=8, print_progress=True) -> xarray.DataArray

   Select two different filtering schemes 1 and 2, and perform the filtering separately.
   And then perform the difference, that means **scheme1 - scheme2**.
   The mesoscale fluctuations are thus preserved.

   .. warning::

       **NOT** support ``data`` contains ``np.nan``.

   Parameters
   ----------
   g1 : :py:class:`float <float>`, generally between (0, 1]
       Constant parameter of scheme1.
   c1 : :py:class:`int <int>`
       Constant parameterof scheme1.
   g2 : :py:class:`float <float>`, generally between (0, 1]
       Constant parameter of scheme2.
   c2 : :py:class:`int <int>`
       Constant parameterof scheme2.
   r :  :py:class:`float <float>`
       The inverse of the maximum response differenc.
       It is prevented from being unduly large and very small difference fields are not greatly amplified.
   print_progress : :py:class:`bool <bool>`
       Whether to print the progress bar when executing computation.

   Returns
   -------
   data_vars : :py:class:`xarray.DataArray <xarray.DataArray>`
       Mesoscale wave field filtered out from raw data

   .. seealso::
       - Maddox, R. A. (1980). An Objective Technique for Separating Macroscale and Mesoscale Features in Meteorological Data. Monthly Weather Review, 108(8), 1108-1121. https://journals.ametsoc.org/view/journals/mwre/108/8/1520-0493_1980_108_1108_aotfsm_2_0_co_2.xml
       - https://github.com/LinOuyang/pybarnes


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


.. py:function:: calc_lanczos_bandpass(data: xarray.DataArray | xarray.Dataset, window_length: int, period: list[int], dim: str = 'time', method: Literal['rolling', 'convolve'] = 'rolling') -> xarray.DataArray

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


.. py:function:: calc_lanczos_lowpass(data: xarray.DataArray | xarray.Dataset, window_length: int, period: int, dim: str = 'time', method: Literal['rolling', 'convolve'] = 'rolling') -> xarray.DataArray

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


.. py:function:: calc_lanczos_highpass(data: xarray.DataArray | xarray.Dataset, window_length: int, period: int, dim: str = 'time', method: Literal['rolling', 'convolve'] = 'rolling') -> xarray.DataArray

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
       - https://github.com/ct6502/wavelets
       - https://github.com/regeirk/pycwt, https://pycwt.readthedocs.io/en/latest/index.html
       - http://nicolasfauchereau.github.io/climatecode/posts/wavelet-analysis-in-python/
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


.. py:function:: calc_redfit(data: xarray.DataArray, timearray: numpy.array = None, nsim: int = 1000, mctest: bool = False, rhopre: float = -99.0, ofac: float = 1.0, hifac: float = 1.0, n50: int = 1, iwin: Literal['rectangular', 'welch', 'hanning', 'triangular', 'blackmanharris'] = 'rectangular')

   Estimating red-noise spectra directly from unevenly spaced paleoclimatic time series.

   Parameters
   ----------
   data: :py:class:`xarray.DataArray<xarray.DataArray>`
       Input time series data
   timearray: :py:class:`numpy.array<numpy.array>`
       Time series data array
   nsim: :py:class:`int<int>`
       Number of Monte-Carlo simulations (1000-2000 should be o.k. in most cases)
   mctest: :py:class:`bool<bool>`
       Toggle calculation of false-alarm levels based on Monte-Carlo simulation,
       if set to `True` : perform Monte-Carlo test,
       if set to `False` : skip Monte-Carlo test (default).
   rhopre: :py:class:`float<float>`
       Prescibed value for :math:`\rho`; unused if < 0 (default = -99.0)
   ofac: :py:class:`float<float>`
       Oversampling factor for Lomb-Scargle Fourier transform (typical values: 2.0-4.0)
   hifac: :py:class:`float<float>`
       Max. frequency to analyze is set to hifac * <fNyq> (default = 1.0)
   n50: :py:class:`int<int>`
       Number of WOSA segments (with 50 % overlap)
   iwin: {"rectangular", "welch", "hanning", "triangular", "blackmanharris"}
       Window-type identifier used to suppress sidelobes in spectral analysis:
       ({"rectangular", "welch", "hanning", "triangular", "blackmanharris"}, optional)

   .. caution::
       Parameters `ofac`, `hifac`, `n50` and window type are identical to the SPECTRUM program
       (see Schulz and Stattegger, 1997 for further details).
       Except mctest, hifac and rhopre all parameters must be specified.

   Returns
   -------
   The red-noise spectra (:py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       - Schulz, M., & Mudelsee, M. (2002). REDFIT: estimating red-noise spectra directly from unevenly spaced paleoclimatic time series [Software]. Computers & Geosciences, 28(3), 421-426. https://doi.org/10.1016/S0098-3004(01)00044-9
       - https://www.marum.de/Prof.-Dr.-michael-schulz/Michael-Schulz-Software.html


.. py:function:: calc_redfit_cross(data_x: xarray.DataArray, data_y: xarray.DataArray, timearray_x: numpy.array = None, timearray_y: numpy.array = None, x_sign: bool = False, y_sign: bool = False, nsim: int = 1000, mctest: bool = True, mctest_phi: bool = True, rhopre_1: float = -999.0, rhopre_2: float = -999.0, ofac: float = 1.0, hifac: float = 1.0, n50: int = 1, alpha: float = 0.05, iwin: Literal['rectangular', 'welch', 'hanning', 'triangular', 'blackmanharris'] = 'rectangular')

   Estimating red-noise spectra directly from unevenly spaced paleoclimatic time series.

   Parameters
   ----------
   data_x::py:class:`xarray.DataArray<xarray.DataArray>`
       First input time series data
   data_y: :py:class:`xarray.DataArray<xarray.DataArray>`
       Second input time series data
   timearray_x: :py:class:`numpy.array<numpy.array>`
       First time series data array
   timearray_y: :py:class:`numpy.array<numpy.array>`
       Second time series data array
   x_sign: :py:class:`bool<bool>`
       Change the sign of the first time series:
       if `True`: The sign of the data is changed
       if `False`: The sign of the data is not changed (default)
   y_sign: :py:class:`bool<bool>`
       Change the sign of the second time series:
       if `True`: The sign of the data is changed
       if `False`: The sign of the data is not changed (default)
   nsim: :py:class:`int<int>`
       Number of Monte Carlo simulations (1000-2000 is recommended)
   mctest: :py:class:`bool<bool>`
       Estimate the significance of auto and coherency spectrum with Monte Carlo simulations
       if `True`: perform Monte Carlo simulations
       if `False`: do not perform Monte Carlo simulations
   mctest_phi: :py:class:`bool<bool>`
       Estimate Monte Carlo confidence interval for the phase spectrum
       if `True`: perform Monte Carlo simulations (mctest needs to be true as well)
       if `False`: do not perform Monte Carlo simulations
   rhopre_1: :py:class:`float<float>`
       Prescribed value for :math:`\rho` for the first time series, not used if :math:`\rho < 0` (default = -999.0).
   rhopre_2: :py:class:`float<float>`
       Prescribed value for :math:`\rho` for the second time series, not used if :math:`\rho< 0` (default = -999.0).
   ofac: :py:class:`float<float>`
       Oversampling factor for Lomb-Scargle Fourier transform (typical values: 2.0-4.0).
   hifac: :py:class:`float<float>`
       Max. frequency to analyze is set to hifac * <fNyq> (default = 1.0).
   n50: :py:class:`int<int>`
       Number of WOSA segments (with 50 % overlap)
   alpha: :py:class:`float<float>`
       Significance level (Note: only 0.01, 0.05 [default], or 0.1 are allowed).
   iwin: {"rectangular", "welch", "hanning", "triangular", "blackmanharris"}
       Window-type identifier used to suppress sidelobes in spectral analysis:
       ({"rectangular", "welch", "hanning", "triangular", "blackmanharris"}, optional).

   .. caution::
       Parameters ofac, hifac, n50 and window type are identical to the SPECTRUM program
       (see Schulz and Stattegger, 1997 for further details).
       Except mctest, hifac, rhopre(1) and rhopre(2) all parameters must be specified.

   .. seealso::
       - Schulz, M., & Mudelsee, M. (2002). REDFIT: estimating red-noise spectra directly from unevenly spaced paleoclimatic time series [Software]. Computers & Geosciences, 28(3), 421-426. https://doi.org/10.1016/S0098-3004(01)00044-9
       - https://www.marum.de/Prof.-Dr.-michael-schulz/Michael-Schulz-Software.html


.. py:function:: kf_filter_wheeler_and_kiladis_1999(input_data: xarray.DataArray, steps_per_day: int | float, tMin: float | Literal[-999], tMax: float | Literal[-999], kMin: int | float | Literal[-999], kMax: int | float | Literal[-999], hMin: float | Literal[-999], hMax: float | Literal[-999], waveName: Literal['kelvin', 'er', 'mrg', 'ig0', 'ig1', 'ig2', None] = None, time_dim='time', lon_dim='lon') -> xarray.DataArray

   Extract equatorial waves by filtering in the Wheeler-Kiladis wavenumber-frequency domain.

   The `kf_filter` function applies a space-time filter to isolate equatorial wave modes based on the Wheeler-Kiladis (WK99) methodology.
   It filters the input data in both the zonal (longitude) and temporal dimensions to retain specific wave modes, such as Kelvin,
   Equatorial Rossby (ER), Mixed Rossby-Gravity (MRG), or Inertia-Gravity (IG) waves.

   At each point, the data are space-time bandpass filtered following Wheeler and Kiladis (1999).
   The data are first detrended with dtrend and tapered with taper in time,
   and then they are filtered using 2-dimensional FFT.
   The filter bounds can simply be rectangular (as in Wheeler and Kiladis's MJO filter),
   or they can bounded by the dispersion curves of the shallow water equatorial waves.
   At this time, other filter shapes (e.g., the TD filter from Kiladis et al. 2006) are not supported.

   Parameters
   -----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The input data from which to remove the smooth daily annual cycle mean.

       .. attention::

               The input data should be periodic (global) in longitude. In addition,
               filtered anomalies near the temporal ends of the dataset should generally be ignored.
               The longer the periods filtered for, the more data should be ignored at the ends.

   steps_per_day : :py:class:`int <int>` or :py:class:`float <float>`
       Number of observations (time steps) per day. For daily data, usually set ``steps_per_day=1``.
   tMin : :py:class:`float <float>`
       Minimum period (in days) for the temporal filter.
   tMax : :py:class:`float <float>`
       Maximum period (in days) for the temporal filter.
   kMin : :py:class:`int <int>` or :py:class:`float <float>`
       Minimum zonal wavenumber for the filter. Positive values indicate eastward propagation,
       while negative values indicate westward propagation.
   kMax : :py:class:`int <int>` or :py:class:`float <float>`
       Maximum zonal wavenumber for the filter. Similar to `kMin`, positive values indicate eastward propagation,
       and negative values indicate westward propagation.
   hMin : :py:class:`float <float>`
       Minimum equivalent depth (in meters) for the dispersion curve filter.
   hMax : :py:class:`float <float>`
       Maximum equivalent depth (in meters) for the dispersion curve filter.
   waveName : str, optional
       Name of dispersion curve to use. Supported options include:
           - ``"kelvin"``: Kelvin waves
           - ``"er"``: Equatorial Rossby waves
           - ``"mrg"`` or ``"ig0"``: Mixed Rossby-Gravity waves (Inertia-Gravity waves [:math:`n=0`])
           - ``"ig1"``: Inertia-Gravity waves [:math:`n=1`]
           - ``"ig2"``: Inertia-Gravity waves [:math:`n=2`]
           - ``None``: Do NOT use dispersion curve.
       **Default**: ``None``.

   Returns
   --------
   :py:class:`xarray.DataArray<xarray.DataArray>`
       Filtered data with the same shape as `data_input`.

   Reference
   --------------
   - Wheeler, M., & Kiladis, G. N. (1999). Convectively Coupled Equatorial Waves: Analysis of Clouds and Temperature in the Wavenumber–Frequency Domain. Journal of the Atmospheric Sciences, 56(3), 374-399. https://journals.ametsoc.org/view/journals/atsc/56/3/1520-0469_1999_056_0374_ccewao_2.0.co_2.xml
   - Kiladis, G. N., Thorncroft, C. D., & Hall, N. M. J. (2006). Three-Dimensional Structure and Dynamics of African Easterly Waves. Part I: Observations. Journal of the Atmospheric Sciences, 63(9), 2212-2230. https://doi.org/10.1175/JAS3741.1
   - Hall, N. M. J., Kiladis, G. N., & Thorncroft, C. D. (2006). Three-Dimensional Structure and Dynamics of African Easterly Waves. Part II: Dynamical Modes. Journal of the Atmospheric Sciences, 63(9), 2231-2245. https://doi.org/10.1175/JAS3742.1
   - Thorncroft, C. D., Hall, N. M. J., & Kiladis, G. N. (2008). Three-Dimensional Structure and Dynamics of African Easterly Waves. Part III: Genesis. Journal of the Atmospheric Sciences, 65(11), 3596-3607. https://doi.org/10.1175/2008JAS2575.1

   .. seealso::

       - https://www.ncl.ucar.edu/Document/Functions/User_contributed/kf_filter.shtml
       - https://ncics.org/portfolio/monitor/mjo/
       - https://k3.cicsnc.org/carl/monitor


.. py:function:: kf_filter_lf_wave(input_data: xarray.DataArray, steps_per_day: int | float, tMin: float | Literal[-999] = 120, tMax: float | Literal[-999] = -999, kMin: int | float | Literal[-999] = -999, kMax: int | float | Literal[-999] = -999, hMin: float | Literal[-999] = -999, hMax: float | Literal[-999] = -999, waveName: Literal['kelvin', 'er', 'mrg', 'ig0', 'ig1', 'ig2', None] = None, time_dim='time', lon_dim='lon') -> xarray.DataArray

   Extract low-frequency waves by filtering in the Wheeler-Kiladis wavenumber-frequency domain. The maximum period is beyond 120 days.

   Parameters
   -----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The input data from which to remove the smooth daily annual cycle mean.

       .. attention::

               The input data should be periodic (global) in longitude. In addition,
               filtered anomalies near the temporal ends of the dataset should generally be ignored.
               The longer the periods filtered for, the more data should be ignored at the ends.

   steps_per_day : :py:class:`int <int>` or :py:class:`float <float>`
       Number of observations (time steps) per day. For daily data, usually set ``steps_per_day=1``.
   tMin : :py:class:`float <float>`
       Minimum period (in days) for the temporal filter.
   tMax : :py:class:`float <float>`
       Maximum period (in days) for the temporal filter.
   kMin : :py:class:`int <int>` or :py:class:`float <float>`
       Minimum zonal wavenumber for the filter. Positive values indicate eastward propagation,
       while negative values indicate westward propagation.
   kMax : :py:class:`int <int>` or :py:class:`float <float>`
       Maximum zonal wavenumber for the filter. Similar to `kMin`, positive values indicate eastward propagation,
       and negative values indicate westward propagation.
   hMin : :py:class:`float <float>`
       Minimum equivalent depth (in meters) for the dispersion curve filter.
   hMax : :py:class:`float <float>`
       Maximum equivalent depth (in meters) for the dispersion curve filter.
   waveName : str, optional
       Name of dispersion curve to use. Supported options include:
           - ``"kelvin"``: Kelvin waves
           - ``"er"``: Equatorial Rossby waves
           - ``"mrg"`` or ``"ig0"``: Mixed Rossby-Gravity waves (Inertia-Gravity waves [:math:`n=0`])
           - ``"ig1"``: Inertia-Gravity waves [:math:`n=1`]
           - ``"ig2"``: Inertia-Gravity waves [:math:`n=2`]
           - ``None``: Do NOT use dispersion curve.
       **Default**: ``None``.

   Returns
   --------
   :py:class:`xarray.DataArray<xarray.DataArray>`
       Filtered data with the same shape as `data_input`.


.. py:function:: kf_filter_mjo_wave(input_data: xarray.DataArray, steps_per_day: int | float, tMin: float | Literal[-999] = 30, tMax: float | Literal[-999] = 96, kMin: int | float | Literal[-999] = 1, kMax: int | float | Literal[-999] = 5, hMin: float | Literal[-999] = -999, hMax: float | Literal[-999] = -999, waveName: Literal['kelvin', 'er', 'mrg', 'ig0', 'ig1', 'ig2', None] = None, time_dim='time', lon_dim='lon') -> xarray.DataArray

   Extract Madden-Julian Oscillation (MJO) waves by filtering in the Wheeler-Kiladis wavenumber-frequency domain

   Parameters
   -----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The input data from which to remove the smooth daily annual cycle mean.

       .. attention::

               The input data should be periodic (global) in longitude. In addition,
               filtered anomalies near the temporal ends of the dataset should generally be ignored.
               The longer the periods filtered for, the more data should be ignored at the ends.

   steps_per_day : :py:class:`int <int>` or :py:class:`float <float>`
       Number of observations (time steps) per day. For daily data, usually set ``steps_per_day=1``.
   tMin : :py:class:`float <float>`
       Minimum period (in days) for the temporal filter.
   tMax : :py:class:`float <float>`
       Maximum period (in days) for the temporal filter.
   kMin : :py:class:`int <int>` or :py:class:`float <float>`
       Minimum zonal wavenumber for the filter. Positive values indicate eastward propagation,
       while negative values indicate westward propagation.
   kMax : :py:class:`int <int>` or :py:class:`float <float>`
       Maximum zonal wavenumber for the filter. Similar to `kMin`, positive values indicate eastward propagation,
       and negative values indicate westward propagation.
   hMin : :py:class:`float <float>`
       Minimum equivalent depth (in meters) for the dispersion curve filter.
   hMax : :py:class:`float <float>`
       Maximum equivalent depth (in meters) for the dispersion curve filter.
   waveName : str, optional
       Name of dispersion curve to use. Supported options include:
           - ``"kelvin"``: Kelvin waves
           - ``"er"``: Equatorial Rossby waves
           - ``"mrg"`` or ``"ig0"``: Mixed Rossby-Gravity waves (Inertia-Gravity waves [:math:`n=0`])
           - ``"ig1"``: Inertia-Gravity waves [:math:`n=1`]
           - ``"ig2"``: Inertia-Gravity waves [:math:`n=2`]
           - ``None``: Do NOT use dispersion curve.
       **Default**: ``None``.

   Returns
   --------
   :py:class:`xarray.DataArray<xarray.DataArray>`
       Filtered data with the same shape as `data_input`.

   Reference
   --------------
   - Kiladis, G. N., Straub, K. H., & Haertel, P. T. (2005). Zonal and Vertical Structure of the Madden–Julian Oscillation. Journal of the Atmospheric Sciences, 62(8), 2790-2809. https://doi.org/10.1175/JAS3520.1


.. py:function:: kf_filter_er_wave(input_data: xarray.DataArray, steps_per_day: int | float, tMin: float | Literal[-999] = 9.7, tMax: float | Literal[-999] = 48, kMin: int | float | Literal[-999] = -10, kMax: int | float | Literal[-999] = -1, hMin: float | Literal[-999] = 8, hMax: float | Literal[-999] = 90, waveName: Literal['kelvin', 'er', 'mrg', 'ig0', 'ig1', 'ig2', None] = 'er', time_dim='time', lon_dim='lon') -> xarray.DataArray

   Extract equatorial Rossby (ER) waves by filtering in the Wheeler-Kiladis wavenumber-frequency domain.

   Parameters
   -----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The input data from which to remove the smooth daily annual cycle mean.

       .. attention::

               The input data should be periodic (global) in longitude. In addition,
               filtered anomalies near the temporal ends of the dataset should generally be ignored.
               The longer the periods filtered for, the more data should be ignored at the ends.

   steps_per_day : :py:class:`int <int>` or :py:class:`float <float>`
       Number of observations (time steps) per day. For daily data, usually set ``steps_per_day=1``.
   tMin : :py:class:`float <float>`
       Minimum period (in days) for the temporal filter.
   tMax : :py:class:`float <float>`
       Maximum period (in days) for the temporal filter.
   kMin : :py:class:`int <int>` or :py:class:`float <float>`
       Minimum zonal wavenumber for the filter. Positive values indicate eastward propagation,
       while negative values indicate westward propagation.
   kMax : :py:class:`int <int>` or :py:class:`float <float>`
       Maximum zonal wavenumber for the filter. Similar to `kMin`, positive values indicate eastward propagation,
       and negative values indicate westward propagation.
   hMin : :py:class:`float <float>`
       Minimum equivalent depth (in meters) for the dispersion curve filter.
   hMax : :py:class:`float <float>`
       Maximum equivalent depth (in meters) for the dispersion curve filter.
   waveName : str, optional
       Name of dispersion curve to use. Supported options include:
           - ``"kelvin"``: Kelvin waves
           - ``"er"``: Equatorial Rossby waves
           - ``"mrg"`` or ``"ig0"``: Mixed Rossby-Gravity waves (Inertia-Gravity waves [:math:`n=0`])
           - ``"ig1"``: Inertia-Gravity waves [:math:`n=1`]
           - ``"ig2"``: Inertia-Gravity waves [:math:`n=2`]
           - ``None``: Do NOT use dispersion curve.
       **Default**: ``None``.

   Returns
   --------
   :py:class:`xarray.DataArray<xarray.DataArray>`
       Filtered data with the same shape as `data_input`.

   Reference
   --------------
   - Kiladis, G. N., M. C. Wheeler, P. T. Haertel, K. H. Straub, and P. E. Roundy (2009), Convectively coupled equatorial waves, Rev. Geophys., 47, RG2003, doi:https://doi.org/10.1029/2008RG000266.


.. py:function:: kf_filter_kelvin_wave(input_data: xarray.DataArray, steps_per_day: int | float, tMin: float | Literal[-999] = 2.5, tMax: float | Literal[-999] = 30, kMin: int | float | Literal[-999] = 1, kMax: int | float | Literal[-999] = 14, hMin: float | Literal[-999] = 8, hMax: float | Literal[-999] = 90, waveName: Literal['kelvin', 'er', 'mrg', 'ig0', 'ig1', 'ig2', None] = 'kelvin', time_dim='time', lon_dim='lon') -> xarray.DataArray

   Extract Kelvin waves by filtering in the Wheeler-Kiladis wavenumber-frequency domain.

   Parameters
   -----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The input data from which to remove the smooth daily annual cycle mean.

       .. attention::

               The input data should be periodic (global) in longitude. In addition,
               filtered anomalies near the temporal ends of the dataset should generally be ignored.
               The longer the periods filtered for, the more data should be ignored at the ends.

   steps_per_day : :py:class:`int <int>` or :py:class:`float <float>`
       Number of observations (time steps) per day. For daily data, usually set ``steps_per_day=1``.
   tMin : :py:class:`float <float>`
       Minimum period (in days) for the temporal filter.
   tMax : :py:class:`float <float>`
       Maximum period (in days) for the temporal filter.
   kMin : :py:class:`int <int>` or :py:class:`float <float>`
       Minimum zonal wavenumber for the filter. Positive values indicate eastward propagation,
       while negative values indicate westward propagation.
   kMax : :py:class:`int <int>` or :py:class:`float <float>`
       Maximum zonal wavenumber for the filter. Similar to `kMin`, positive values indicate eastward propagation,
       and negative values indicate westward propagation.
   hMin : :py:class:`float <float>`
       Minimum equivalent depth (in meters) for the dispersion curve filter.
   hMax : :py:class:`float <float>`
       Maximum equivalent depth (in meters) for the dispersion curve filter.
   waveName : str, optional
       Name of dispersion curve to use. Supported options include:
           - ``"kelvin"``: Kelvin waves
           - ``"er"``: Equatorial Rossby waves
           - ``"mrg"`` or ``"ig0"``: Mixed Rossby-Gravity waves (Inertia-Gravity waves [:math:`n=0`])
           - ``"ig1"``: Inertia-Gravity waves [:math:`n=1`]
           - ``"ig2"``: Inertia-Gravity waves [:math:`n=2`]
           - ``None``: Do NOT use dispersion curve.
       **Default**: ``None``.

   Returns
   --------
   :py:class:`xarray.DataArray<xarray.DataArray>`
       Filtered data with the same shape as `data_input`.

   Reference
   --------------
   - Straub, K. H., & Kiladis, G. N. (2002). Observations of a Convectively Coupled Kelvin Wave in the Eastern Pacific ITCZ. Journal of the Atmospheric Sciences, 59(1), 30-53. https://journals.ametsoc.org/view/journals/atsc/59/1/1520-0469_2002_059_0030_ooacck_2.0.co_2.xml


.. py:function:: kf_filter_mt_wave(input_data: xarray.DataArray, steps_per_day: int | float, tMin: float | Literal[-999] = 2.5, tMax: float | Literal[-999] = 10, kMin: int | float | Literal[-999] = -14, kMax: int | float | Literal[-999] = 0, hMin: float | Literal[-999] = -999.0, hMax: float | Literal[-999] = -999.0, waveName: Literal['kelvin', 'er', 'mrg', 'ig0', 'ig1', 'ig2', None] = None, time_dim='time', lon_dim='lon') -> xarray.DataArray

   Extract mixed Rossby-gravity (MRG)-tropical depression (TD) type waves by filtering in the Wheeler-Kiladis wavenumber-frequency domain

   Parameters
   -----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The input data from which to remove the smooth daily annual cycle mean.

       .. attention::

               The input data should be periodic (global) in longitude. In addition,
               filtered anomalies near the temporal ends of the dataset should generally be ignored.
               The longer the periods filtered for, the more data should be ignored at the ends.

   steps_per_day : :py:class:`int <int>` or :py:class:`float <float>`
       Number of observations (time steps) per day. For daily data, usually set ``steps_per_day=1``.
   tMin : :py:class:`float <float>`
       Minimum period (in days) for the temporal filter.
   tMax : :py:class:`float <float>`
       Maximum period (in days) for the temporal filter.
   kMin : :py:class:`int <int>` or :py:class:`float <float>`
       Minimum zonal wavenumber for the filter. Positive values indicate eastward propagation,
       while negative values indicate westward propagation.
   kMax : :py:class:`int <int>` or :py:class:`float <float>`
       Maximum zonal wavenumber for the filter. Similar to `kMin`, positive values indicate eastward propagation,
       and negative values indicate westward propagation.
   hMin : :py:class:`float <float>`
       Minimum equivalent depth (in meters) for the dispersion curve filter.
   hMax : :py:class:`float <float>`
       Maximum equivalent depth (in meters) for the dispersion curve filter.
   waveName : str, optional
       Name of dispersion curve to use. Supported options include:
           - ``"kelvin"``: Kelvin waves
           - ``"er"``: Equatorial Rossby waves
           - ``"mrg"`` or ``"ig0"``: Mixed Rossby-Gravity waves (Inertia-Gravity waves [:math:`n=0`])
           - ``"ig1"``: Inertia-Gravity waves [:math:`n=1`]
           - ``"ig2"``: Inertia-Gravity waves [:math:`n=2`]
           - ``None``: Do NOT use dispersion curve.
       **Default**: ``None``.

   Returns
   --------
   :py:class:`xarray.DataArray<xarray.DataArray>`
       Filtered data with the same shape as `data_input`.

   Reference
   --------------
   - Frank, W. M., & Roundy, P. E. (2006). The Role of Tropical Waves in Tropical Cyclogenesis. Monthly Weather Review, 134(9), 2397-2417. https://doi.org/10.1175/MWR3204.1


.. py:function:: kf_filter_mrg_wave(input_data: xarray.DataArray, steps_per_day: int | float, tMin: float | Literal[-999] = 3, tMax: float | Literal[-999] = 9.6, kMin: int | float | Literal[-999] = -10, kMax: int | float | Literal[-999] = -1, hMin: float | Literal[-999] = 8, hMax: float | Literal[-999] = 90, waveName: Literal['kelvin', 'er', 'mrg', 'ig0', 'ig1', 'ig2', None] = 'mrg', time_dim='time', lon_dim='lon') -> xarray.DataArray

   Extract mixed Rossby-gravity waves (MRG) by filtering in the Wheeler-Kiladis wavenumber-frequency domain.

   Parameters
   -----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The input data from which to remove the smooth daily annual cycle mean.

       .. attention::

               The input data should be periodic (global) in longitude. In addition,
               filtered anomalies near the temporal ends of the dataset should generally be ignored.
               The longer the periods filtered for, the more data should be ignored at the ends.

   steps_per_day : :py:class:`int <int>` or :py:class:`float <float>`
       Number of observations (time steps) per day. For daily data, usually set ``steps_per_day=1``.
   tMin : :py:class:`float <float>`
       Minimum period (in days) for the temporal filter.
   tMax : :py:class:`float <float>`
       Maximum period (in days) for the temporal filter.
   kMin : :py:class:`int <int>` or :py:class:`float <float>`
       Minimum zonal wavenumber for the filter. Positive values indicate eastward propagation,
       while negative values indicate westward propagation.
   kMax : :py:class:`int <int>` or :py:class:`float <float>`
       Maximum zonal wavenumber for the filter. Similar to `kMin`, positive values indicate eastward propagation,
       and negative values indicate westward propagation.
   hMin : :py:class:`float <float>`
       Minimum equivalent depth (in meters) for the dispersion curve filter.
   hMax : :py:class:`float <float>`
       Maximum equivalent depth (in meters) for the dispersion curve filter.
   waveName : str, optional
       Name of dispersion curve to use. Supported options include:
           - ``"kelvin"``: Kelvin waves
           - ``"er"``: Equatorial Rossby waves
           - ``"mrg"`` or ``"ig0"``: Mixed Rossby-Gravity waves (Inertia-Gravity waves [:math:`n=0`])
           - ``"ig1"``: Inertia-Gravity waves [:math:`n=1`]
           - ``"ig2"``: Inertia-Gravity waves [:math:`n=2`]
           - ``None``: Do NOT use dispersion curve.
       **Default**: ``None``.

   Returns
   --------
   :py:class:`xarray.DataArray<xarray.DataArray>`
       Filtered data with the same shape as `data_input`.


.. py:function:: kf_filter_td_wave(input_data: xarray.DataArray, steps_per_day: int | float, tMin: float | Literal[-999] = 2, tMax: float | Literal[-999] = 8.5, kMin: int | float | Literal[-999] = -15, kMax: int | float | Literal[-999] = -6, hMin: float | Literal[-999] = 90, hMax: float | Literal[-999] = -999, waveName: Literal['kelvin', 'er', 'mrg', 'ig0', 'ig1', 'ig2', None] = 'mrg', time_dim='time', lon_dim='lon') -> xarray.DataArray

   Extract tropical depression (TD) by filtering in the Wheeler-Kiladis wavenumber-frequency domain.

   Parameters
   -----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The input data from which to remove the smooth daily annual cycle mean.

       .. attention::

               The input data should be periodic (global) in longitude. In addition,
               filtered anomalies near the temporal ends of the dataset should generally be ignored.
               The longer the periods filtered for, the more data should be ignored at the ends.

   steps_per_day : :py:class:`int <int>` or :py:class:`float <float>`
       Number of observations (time steps) per day. For daily data, usually set ``steps_per_day=1``.
   tMin : :py:class:`float <float>`
       Minimum period (in days) for the temporal filter.
   tMax : :py:class:`float <float>`
       Maximum period (in days) for the temporal filter.
   kMin : :py:class:`int <int>` or :py:class:`float <float>`
       Minimum zonal wavenumber for the filter. Positive values indicate eastward propagation,
       while negative values indicate westward propagation.
   kMax : :py:class:`int <int>` or :py:class:`float <float>`
       Maximum zonal wavenumber for the filter. Similar to `kMin`, positive values indicate eastward propagation,
       and negative values indicate westward propagation.
   hMin : :py:class:`float <float>`
       Minimum equivalent depth (in meters) for the dispersion curve filter.
   hMax : :py:class:`float <float>`
       Maximum equivalent depth (in meters) for the dispersion curve filter.
   waveName : str, optional
       Name of dispersion curve to use. Supported options include:
           - ``"kelvin"``: Kelvin waves
           - ``"er"``: Equatorial Rossby waves
           - ``"mrg"`` or ``"ig0"``: Mixed Rossby-Gravity waves (Inertia-Gravity waves [:math:`n=0`])
           - ``"ig1"``: Inertia-Gravity waves [:math:`n=1`]
           - ``"ig2"``: Inertia-Gravity waves [:math:`n=2`]
           - ``None``: Do NOT use dispersion curve.
       **Default**: ``None``.

   Returns
   --------
   :py:class:`xarray.DataArray<xarray.DataArray>`
       Filtered data with the same shape as `data_input`.


.. py:function:: filter_2D_spatial_parabolic_cylinder_function(zonal_wind_speed_data: xarray.DataArray, meridional_wind_speed_data: xarray.DataArray, z_data: xarray.DataArray, period_min=3.0, period_max=30.0, wavenumber_min=2, wavenumber_max=20, trapping_scale_deg=6.0, lon_dim='lon', lat_dim='lat', time_dim: str = 'time')

   Perform space-time spectral analysis to extract equatorial wave components.

   This function filters atmospheric fields to isolate equatorial wave modes, including Kelvin waves,
   westward-moving mixed Rossby-gravity (WMRG) waves, and equatorial Rossby waves of the first and
   second kind (R1 and R2), within user-specified period and wavenumber ranges. It processes input
   data using a combination of detrending, windowing, Fourier transforms, and projection onto
   parabolic cylinder functions to decompose the fields into these wave types.

   Parameters
   ----------

   zonal_wind_speed_data: :class:`xarray.DataArray <xarray.DataArray>`
       The zonal (east-west) wind speed component.
   meridional_wind_speed_data: :class:`xarray.DataArray <xarray.DataArray>`
       The meridional (north-south) wind speed component.
   z_data: :class:`xarray.DataArray <xarray.DataArray>`
       The geopotential height.
   period_min: :class:`float`, optional (default=3.0)
       The minimum period (in days) of the waves to be extracted.
   period_max: :class:`float`, optional (default=30.0)
       The maximum period (in days) of the waves to be extracted.
   wavenumber_min: :class:`int`, optional (default=2)
       The minimum zonal wavenumber of the waves to be extracted.
   wavenumber_max: :class:`int`, optional (default=20)
       The maximum zonal wavenumber of the waves to be extracted.
   trapping_scale_deg: :class:`float`, optional (default=6.0)
       The meridional trapping scale in degrees, defining the width of the equatorial waveguide
       for the parabolic cylinder functions.
   lon_dim: :class:`str`, optional (default="lon")
       The name of the longitude dimension in the input data arrays.
   lat_dim: :class:`str`, optional (default="lat")
       The name of the latitude dimension in the input data arrays.
   time_dim: :class:`str`, optional (default="time")
       The name of the time dimension in the input data arrays.

   Returns
   -------

   :class:`xarray.Dataset <xarray.Dataset>`

   A dataset containing the filtered fields for each wave type:
       - ``'u'``: Zonal wind component (m/s)
       - ``'v'``: Meridional wind component (m/s)
       - ``'z'``: Geopotential height (m)
   Each variable includes an additional ``'wave_type'`` dimension with values:
   ``['kelvin', 'wmrg', 'r1', 'r2']``.

   Mathematical Explanation
   ------------------------

   The function implements a space-time spectral analysis method to decompose atmospheric fields into
   equatorial wave modes, based on the normal modes of the shallow water equations on an equatorial
   beta-plane. The meridional structure of these waves is represented by parabolic cylinder functions,
   which arise as solutions to the quantum harmonic oscillator problem adapted to the equatorial
   waveguide.

   Key Steps
   ~~~~~~~~~

   1. Detrending and Windowing:
       - The input fields are spatially detrended to remove large-scale trends.
       - A Tukey window (alpha=0.1) is applied along the time dimension to minimize spectral leakage.

   2. Variable Transformation:
       Following Yang et al. (2003), new variables :math:`q` and :math:`r` are defined:

       .. math::
           q = \frac{g}{c_e} z + u, \quad r = \frac{g}{c_e} z - u

       where:
           - :math:`g = 9.8 \text{m/s}^2` is gravitational acceleration.
           - :math:`c_e = 2 \beta L^2` is a characteristic speed.
           - :math:`\beta = 2.3 \times 10^{-11}, \text{m}^{-1}\text{s}^{-1}` is the beta parameter.
           - :math:`L = \frac{2\pi R_e \cdot \text{trapping_scale_deg}}{360}` is the trapping scale in meters (:math:`R_e = 6.371 \times 10^6, \text{m}` is Earth's radius).
           - :math:`z` is geopotential height, and :math:`u` is zonal wind speed.

   3. Fourier Transform:
       A 2D Fast Fourier Transform (FFT) is applied along the time and longitude dimensions to convert
       the fields into frequency-wavenumber space.

   4. Projection onto Parabolic Cylinder Functions:
       - The spectral coefficients are projected onto parabolic cylinder functions :math:`D_n(y)`, where
           :math:`y = \text{latitude} / \text{trapping_scale_deg}` is the scaled latitude.
       - The first four modes (:math:`n = 0, 1, 2, 3`) are computed with normalization factors:

       .. math::

           D_0(y) = e^{-y^2/4}, \quad D_1(y) = y e^{-y^2/4}, \quad D_2(y) = (y^2 - 1) e^{-y^2/4}, \quad D_3(y) = y (y^2 - 3) e^{-y^2/4}

       - These functions define the meridional structure of the waves.

   5. Wave Type Identification and Filtering:
       - Kelvin Wave: Uses :math:`n = 0` mode, eastward propagating within specified frequency (:math:`1/\text{period_max}` to :math:`1/\text{period_min}`) and wavenumber ranges.
       - WMRG Wave: Combines :math:`n = 1` for :math:`q` and :math:`n = 0` for meridional wind :math:`v`, westward propagating.
       - R1 Wave: Uses :math:`n = 2` for :math:`q`, :math:`n = 0` for :math:`r`, and :math:`n = 1` for :math:`v`, westward propagating.
       - R2 Wave: Uses :math:`n = 3` for :math:`q`, :math:`n = 1` for :math:`r`, and :math:`n = 2` for :math:`v`, westward propagating.
       - Frequencies and wavenumbers are selected based on ``period_min``, ``period_max``, ``wavenumber_min``, and ``wavenumber_max``.

   6. Reconstruction:
       - The filtered spectral coefficients are recombined with the parabolic cylinder functions and transformed back to physical space using an inverse 2D FFT.
       - The physical fields are reconstructed as:
           - :math:`u = (q + r) / 2`.
           - :math:`v` directly from its coefficients.
           - :math:`z = (q - r) \cdot c_e / (2g)`.

   References
   ----------
   - Gill, A.E. (1980), Some simple solutions for heat-induced tropical circulation. Q.J.R. Meteorol. Soc., 106: 447-462. https://doi.org/10.1002/qj.49710644905
   - Li, X.f., Cho, HR. Development and propagation of equatorial waves. Adv. Atmos. Sci. 14, 323–338 (1997). https://doi.org/10.1007/s00376-997-0053-6
   - Yang, G.-Y., Hoskins, B., & Slingo, J. (2003). Convectively coupled equatorial waves: A new methodology for identifying wave structures in observational data. *Journal of the Atmospheric Sciences*, 60(14), 1637-1654.
   - Knippertz, P., Gehne, M., Kiladis, G.N., Kikuchi, K., Rasheeda Satheesh, A., Roundy, P.E., et al. (2022) The intricacies of identifying equatorial waves. Quarterly Journal of the Royal Meteorological Society, 148(747), 2814–2852. Available from: https://doi.org/10.1002/qj.4338


