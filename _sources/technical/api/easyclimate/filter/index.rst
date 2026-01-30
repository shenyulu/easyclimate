easyclimate.filter
==================

.. py:module:: easyclimate.filter


Submodules
----------

.. toctree::
   :maxdepth: 1

   /technical/api/easyclimate/filter/barnes_filter/index
   /technical/api/easyclimate/filter/butter_filter/index
   /technical/api/easyclimate/filter/emd/index
   /technical/api/easyclimate/filter/gaussian_filter/index
   /technical/api/easyclimate/filter/kf_filter/index
   /technical/api/easyclimate/filter/lanczos_filter/index
   /technical/api/easyclimate/filter/redfit/index
   /technical/api/easyclimate/filter/smooth/index
   /technical/api/easyclimate/filter/spatial_pcf/index
   /technical/api/easyclimate/filter/spectrum/index
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
   easyclimate.filter.filter_emd
   easyclimate.filter.filter_eemd
   easyclimate.filter.calc_gaussian_filter
   easyclimate.filter.calc_time_spectrum
   easyclimate.filter.calc_mean_fourier_amplitude
   easyclimate.filter.filter_fourier_harmonic_analysis


Package Contents
----------------

.. py:function:: calc_barnes_lowpass(data: xarray.DataArray, g: float = 0.3, c: int = 150000, lon_dim='lon', lat_dim='lat', radius_degree=8, print_progress=True) -> xarray.DataArray

   Selecting different parameters **g** and **c**
   will result in different filtering characteristics.

   .. warning::

       **NOT** support ``data`` contains ``np.nan``.


   Parameters
   ----------
   data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The spatio-temporal data to be calculated.
   g : :py:class:`float <float>`, generally between (0, 1], default `0.3`.
       Constant parameter.
   c : :py:class:`int <int>`, default `150000`.
       Constant parameter. When *c* takes a larger value, the filter function converges
       at a larger wavelength, and the response function slowly approaches the maximum value,
       which means that high-frequency fluctuations have been filtered out.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   radius_degree : :py:class:`int <int>` or :py:class:`tuple <tuple>` (degree), default `8`.
       The radius of each point when caculating the distance of each other.

       It is recommended to set this with your schemes.
       For the constant ``c``, this parameter is recommended to be:

       for the ``c`` is ``[500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000]``
       ``radius_degree`` is recommended for ``[1, 1.5, 2, 3, 4, 5, 7, 8, 12]``
   print_progress : :py:class:`bool <bool>`
       Whether to print the progress bar when executing computation.

   Returns
   -------
   data_vars : :py:class:`xarray.DataArray <xarray.DataArray>`
       Data field after filtering out high-frequency fluctuations

   .. seealso::
       - Maddox, R. A. (1980). An Objective Technique for Separating Macroscale and Mesoscale Features in Meteorological Data. Monthly Weather Review, 108(8), 1108-1121. https://journals.ametsoc.org/view/journals/mwre/108/8/1520-0493_1980_108_1108_aotfsm_2_0_co_2.xml
       - https://github.com/LinOuyang/pybarnes

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_barnes_filter.py


.. py:function:: calc_barnes_bandpass(data: xarray.DataArray, g1: float = 0.3, g2: float = 0.3, c1: int = 30000, c2: int = 150000, r=1.2, lon_dim='lon', lat_dim='lat', radius_degree=8, print_progress=True) -> xarray.DataArray

   Select two different filtering schemes 1 and 2, and perform the filtering separately.
   And then perform the difference, that means **scheme1 - scheme2**.
   The mesoscale fluctuations are thus preserved.

   .. warning::

       **NOT** support ``data`` contains ``np.nan``.

   Parameters
   ----------
   data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The spatio-temporal data to be calculated.
   g1 : :py:class:`float <float>`, generally between (0, 1], default `0.3`.
       Constant parameter of scheme1.
   g2 : :py:class:`float <float>`, generally between (0, 1], default `0.3`.
       Constant parameter of scheme2.
   c1 : :py:class:`int <int>`, default `30000`.
       Constant parameterof scheme1.
   c2 : :py:class:`int <int>`, default `150000`.
       Constant parameterof scheme2.
   r :  :py:class:`float <float>`, default `1.2`.
       The inverse of the maximum response differenc.
       It is prevented from being unduly large and very small difference fields are not greatly amplified.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   radius_degree : :py:class:`int <int>` or :py:class:`tuple <tuple>` (degree), default `8`.
       The radius of each point when caculating the distance of each other.

       It is recommended to set this with your schemes.
       For the constant ``c``, this parameter is recommended to be:

       for the ``c`` is ``[500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000]``,
       ``radius_degree`` is recommended for ``[1, 1.5, 2, 3, 4, 5, 7, 8, 12]``
   print_progress : :py:class:`bool <bool>`
       Whether to print the progress bar when executing computation.

   Returns
   -------
   data_vars : :py:class:`xarray.DataArray <xarray.DataArray>`
       Mesoscale wave field filtered out from raw data

   .. seealso::
       - Maddox, R. A. (1980). An Objective Technique for Separating Macroscale and Mesoscale Features in Meteorological Data. Monthly Weather Review, 108(8), 1108-1121. https://journals.ametsoc.org/view/journals/mwre/108/8/1520-0493_1980_108_1108_aotfsm_2_0_co_2.xml
       - https://github.com/LinOuyang/pybarnes

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_barnes_filter.py


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

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_wavelet.py


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

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_wavelet.py


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

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_wavelet.py


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

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_redfit.py


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

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_kf_filter.py


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

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_kf_filter.py


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

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_kf_filter.py


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

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_kf_filter.py


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

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_kf_filter.py


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

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_kf_filter.py


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

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_kf_filter.py


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

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_kf_filter.py


.. py:function:: filter_2D_spatial_parabolic_cylinder_function(zonal_wind_speed_data: xarray.DataArray, meridional_wind_speed_data: xarray.DataArray, z_data: xarray.DataArray, period_min=3.0, period_max=30.0, wavenumber_min=2, wavenumber_max=20, trapping_scale_deg=6.0, lon_dim='lon', lat_dim='lat', time_dim: str = 'time', complex_dtype=np.complex128, real_dtype=np.float64)

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

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_spatial_pcf.py


.. py:function:: filter_emd(input_data: xarray.DataArray, time_step: Literal['ns', 'us', 'ms', 's', 'm', 'h', 'D', 'W', 'M', 'Y'], time_array=None, time_dim: str = 'time', spline_kind: Literal['cubic', 'akima', 'pchip', 'cubic_hermite', 'slinear', 'quadratic', 'linear'] = 'cubic', nbsym: int = 2, max_iteration: int = 1000, energy_ratio_thr: float = 0.2, std_thr: float = 0.2, svar_thr: float = 0.001, total_power_thr: float = 0.005, range_thr: float = 0.001, extrema_detection: Literal['simple', 'parabol'] = 'simple', max_imf: int = -1, dtype=np.float64)

   Empirical Mode Decomposition

   Method of decomposing signal into Intrinsic Mode Functions (IMFs) based on algorithm presented in Huang et al (1998).

   Algorithm was validated with Rilling et al (2003). Matlab’s version from 3.2007.

   Threshold which control the goodness of the decomposition:

   - ``std_thr``: Test for the proto-IMF how variance changes between siftings.
   - ``svar_thr``: Test for the proto-IMF how energy changes between siftings.
   - ``total_power_thr``: Test for the whole decomp how much of energy is solved.
   - ``range_thr``: Test for the whole decomp whether the difference is tiny.

   Parameters
   ----------
   input_data: :py:class:`xarray.DataArray <xarray.DataArray>`
       Input signal data to be decomposed.
   time_step: :py:class:`str <str>`
       Time step unit for datetime conversion (e.g., 's', 'ms', 'us', 'ns').
   time_array: :py:class:`array-like <array-like>`, optional
       Custom time array for the input signal. If None, uses input_data's time dimension.
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name.
   spline_kind: ``Literal["cubic", "akima", "pchip", "cubic_hermite", "slinear", "quadratic", "linear"]``, default: "cubic"
       Type of spline used for envelope interpolation.
   nbsym: :py:class:`int <int>`, default: 2
       Number of points to add at signal boundaries for mirroring.
   max_iteration: :py:class:`int <int>`, default: 1000
       Maximum number of iterations per single sifting in EMD.
   energy_ratio_thr: :py:class:`float <float>`, default: 0.2
       Threshold value on energy ratio per IMF check.
   std_thr: :py:class:`float <float>`, default: 0.2
       Threshold value on standard deviation per IMF check.
   svar_thr: :py:class:`float <float>`, default: 0.001
       Threshold value on scaled variance per IMF check.
   total_power_thr: :py:class:`float <float>`, default: 0.005
       Threshold value on total power per EMD decomposition.
   range_thr: :py:class:`float <float>`, default: 0.001
       Threshold for amplitude range (after scaling) per EMD decomposition.
   extrema_detection: ``Literal["simple", "parabol"]``, default: "simple"
       Method used to finding extrema.
   max_imf: :py:class:`int <int>`, default: -1
       IMF number to which decomposition should be performed. Negative value means all.
   dtype: :py:class:`numpy.dtype <numpy.dtype>`, default: np.float64
       Data type used for calculations.


   Reference
   --------------
   - Huang Norden E., Shen Zheng, Long Steven R., Wu Manli C., Shih Hsing H., Zheng Quanan, Yen Nai-Chyuan, Tung Chi Chao and Liu Henry H. 1998 The empirical mode decomposition and the Hilbert spectrum for nonlinear and non-stationary time series analysisProc. R. Soc. Lond. A.454903–995 http://doi.org/10.1098/rspa.1998.0193
   - Gabriel Rilling, Patrick Flandrin, Paulo Gonçalves. On empirical mode decomposition and its algorithms. IEEE-EURASIP Workshop on Nonlinear Signal and Image Processing NSIP-03, Jun 2003, Grado, Italy. https://inria.hal.science/inria-00570628v1
   - Colominas, M. A., Schlotthauer, G., and Torres, M. E. (2014). Improved complete ensemble EMD: A suitable tool for biomedical signal processing. Biomedical Signal Processing and Control, 14, 19-29. https://doi.org/10.1016/j.bspc.2014.06.009
   - Hsuan, R. (2014). Ensemble Empirical Mode Decomposition Parameters Optimization for Spectral Distance Measurement in Hyperspectral Remote Sensing Data. Remote Sens. 6(3), 2069-2083. http://doi.org/10.3390/rs6032069. http://www.mdpi.com/2072-4292/6/3/2069
   - Kim, D., and HS. Uh (2009). EMD: A Package for Empirical Mode Decomposition and Hilbert Spectrum. https://journal.r-project.org/archive/2009-1/RJournal_2009-1_Kim+Oh.pdf
   - Lambert et al. Empirical Mode Decomposition. https://www.clear.rice.edu/elec301/Projects02/empiricalMode/
   - Meta Trader. Introduction to the Empirical Mode Decomposition Method. https://www.mql5.com/en/articles/439
   - Salisbury, J.I. and Wimbush, M. (2002). Using modern time series analysis techniques to predict ENSO events from the SOI time series. Nonlinear Processes in Geophysics, 9, 341-345. http://www.nonlin-processes-geophys.net/9/341/2002/npg-9-341-2002.pdf
   - Torres, M. E., Colominas, M. A., Schlotthauer, G., & Flandrin, P. (2011). A complete ensemble empirical mode decomposition with adaptive noise. ICASSP, 4144-4147. http://doi.org/10.1109/ICASSP.2011.5947265.
   - Wang, T., M. Zhang, Q. Yu, and H. Zhang (2012). Comparing the applications of EMD and EEMD on time-frequency analysis of seismic signal. J. Appl. Geophys., 83, 29-34. http://doi.org/10.1016/j.jappgeo.2012.05.002.
   - Wu, Z., & Huang, N. E. (2009). Ensemble empirical mode decomposition: a noise-assisted data analysis method. Advances in Adaptive Data Analysis, 01(01), 1-41. https://doi.org/10.1142/S1793536909000047
   - Wu, Z, et al (2015). Fast multidimensional ensemble empirical mode decomposition for the analysis of big spatio-temporal datasets. Philos Trans A Math Phys Eng Sci, 374(2065), 20150197. http://doi.org/10.1098/rsta.2015.0197. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4792406/
   - Wu, Y. and Shen, BW (2016). An Evaluation of the Parallel Ensemble Empirical Mode Decomposition Method in Revealing the Role of Downscaling Processes Associated with African Easterly Waves in Tropical Cyclone Genesis. http://doi.org/10.1175/JTECH-D-15-0257.1. http://journals.ametsoc.org/doi/abs/10.1175/JTECH-D-15-0257.1

   .. seealso::

       - https://pyemd.readthedocs.io/
       - https://www.ncl.ucar.edu/Applications/eemd.shtml

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_emd.py


.. py:function:: filter_eemd(input_data: xarray.DataArray, time_step, time_array=None, time_dim: str = 'time', noise_seed: None | int = None, trials: int = 100, noise_width: float = 0.05, parallel: bool = False, processes: None | int = None, separate_trends: bool = False, spline_kind: Literal['cubic', 'akima', 'pchip', 'cubic_hermite', 'slinear', 'quadratic', 'linear'] = 'cubic', nbsym: int = 2, max_iteration: int = 1000, energy_ratio_thr: float = 0.2, std_thr: float = 0.2, svar_thr: float = 0.001, total_power_thr: float = 0.005, range_thr: float = 0.001, extrema_detection: Literal['simple', 'parabol'] = 'parabol', dtype=np.float64)

   Ensemble Empirical Mode Decomposition (EEMD)

   Ensemble empirical mode decomposition (EEMD) is noise-assisted technique (Wu & Huang, 2009),
   which is meant to be more robust than simple Empirical Mode Decomposition (EMD).
   The robustness is checked by performing many decompositions on signals slightly perturbed from their initial position.
   In the grand average over all IMF results the noise will cancel each other out and the result is pure decomposition.

   Parameters
   ----------
   input_data: :py:class:`xarray.DataArray <xarray.DataArray>`
       Input signal data to be decomposed.
   time_step: :py:class:`str <str>`
       Time step unit for datetime conversion (e.g., 's', 'ms', 'us', 'ns').
   time_array: :py:class:`array-like <array-like>`, optional
       Custom time array for the input signal. If None, uses input_data's time dimension.
   time_dim: :py:class:`str <str>`, default: "time"
       The time coordinate dimension name.
   noise_seed: :py:class:`int <int>` or None, default: None
       Set seed for noise generation.

       .. warning::

           Given the nature of EEMD, each time you decompose a signal you will obtain a different set of components. That’s the expected consequence of adding noise which is going to be random. To make the decomposition reproducible, one needs to set a seed for the random number generator used in EEMD.

   trials: :py:class:`int <int>`, default: 100
       Number of trials or EMD performance with added noise.
   noise_width: :py:class:`float <float>`, default: 0.05
       Standard deviation of Gaussian noise (:math:`\hat\sigma`).
       It's relative to absolute amplitude of the signal, i.e.
       :math:`\hat\sigma = \sigma\cdot|\max(S)-\min(S)|`, where
       :math:`\sigma` is noise_width.
   parallel: :py:class:`bool <bool>`, default: False
       Flag whether to use multiprocessing in EEMD execution. Since each EMD(s+noise) is independent this should improve execution speed considerably. Note that it’s disabled by default because it’s the most common problem when EEMD takes too long time to finish. If you set the flag to True, make also sure to set processes to some reasonable value.
   processes: :py:class:`int <int>` or None, default: None
       Number of processes harness when executing in parallel mode. The value should be between 1 and max that depends on your hardware. If None, uses all available cores.
   separate_trends: :py:class:`bool <bool>`, default: False
       Flag whether to isolate trends from each EMD decomposition into a separate component. If ``True``, the resulting EEMD will contain ensemble only from IMFs and the mean residue will be stacked as the last element.
   spline_kind: ``Literal["cubic", "akima", "pchip", "cubic_hermite", "slinear", "quadratic", "linear"]``, default: "cubic"
       Type of spline used for envelope interpolation.
   nbsym: :py:class:`int <int>`, default: 2
       Number of points to add at signal boundaries for mirroring.
   max_iteration: :py:class:`int <int>`, default: 1000
       Maximum number of iterations per single sifting in EMD.
   energy_ratio_thr: :py:class:`float <float>`, default: 0.2
       Threshold value on energy ratio per IMF check.
   std_thr: :py:class:`float <float>`, default: 0.2
       Threshold value on standard deviation per IMF check.
   svar_thr: :py:class:`float <float>`, default: 0.001
       Threshold value on scaled variance per IMF check.
   total_power_thr: :py:class:`float <float>`, default: 0.005
       Threshold value on total power per EMD decomposition.
   range_thr: :py:class:`float <float>`, default: 0.001
       Threshold for amplitude range (after scaling) per EMD decomposition.
   extrema_detection: ``Literal["simple", "parabol"]``, default: "parabol"
       Method used to finding extrema.
   dtype: :py:class:`numpy.dtype <numpy.dtype>`, default: np.float64
       Data type used for calculations.

   Returns
   -------
   :py:class:`xarray.Dataset <xarray.Dataset>`
       Dataset containing the input data and the decomposed eIMFs (Ensemble IMFs).

   Reference
   --------------
   Wu, Z., & Huang, N. E. (2009). Ensemble empirical mode decomposition: a noise-assisted data analysis method. Advances in Adaptive Data Analysis, 01(01), 1-41. https://doi.org/10.1142/S1793536909000047

   .. seealso::
       - :func:`filter_emd` : Standard Empirical Mode Decomposition
       - https://pyemd.readthedocs.io/
       - https://www.ncl.ucar.edu/Applications/eemd.shtml

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_emd.py


.. py:function:: calc_gaussian_filter(data: xarray.DataArray, window_length: float, sigma: float = None, dim: str = 'time', keep_attrs: bool = False) -> xarray.DataArray

   Apply a Gaussian filter to data along a specified dimension.

   Parameters
   ----------
   da : :py:class:`xarray.DataArray<xarray.DataArray>`
       Input data array
   window_length: :py:class:`int <int>`, optional
       The window width.
   sigma : :py:class:`float <float>`, optional
       Standard deviation for Gaussian kernel. If None, calculated as :math:`\mathrm{window_length} / \sqrt{8 log2}`.
   dim : :py:class:`str <str>`, optional
       Dimension along which to filter (default: 'time')
   keep_attrs : :py:class:`bool <bool>`, optional
       Whether to preserve attributes (default: False)

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`
       Smoothed data with the same dimensions as input.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_ao_index.py


.. py:function:: calc_time_spectrum(data: xarray.DataArray, time_dim: str = 'time', inv: bool = False) -> tuple[xarray.DataArray, xarray.DataArray]

   Calculate the time spectrum of a multi-dimensional :py:class:`xarray.DataArray <xarray.DataArray>` along the time dimension.

   Parameters
   ----------
   data : :py:class:`xarray.DataArray <xarray.DataArray>`
       Input data with at least a time dimension.
   time_dim : :py:class:`str <str>`, optional
       Name of the time dimension in the DataArray, by default ``'time'``
   inv : :py:class:`bool <bool>`, optional
       If True, use forward FFT instead of inverse, by default ``False``

   Returns
   -------
   :py:class:`easyclimate.DataNode <easyclimate.core.datanode.DataNode>`, containing:

   - Amplitude spectrum (same dimensions as input but with ``freq`` and ``period`` dimension)
   - Frequency or period values

   .. tip::

       The function automatically handles even/odd length time series and removes the
       zero-frequency component. The amplitude is calculated as the power spectrum (:math:`|\mathrm{fft}|^2`).


.. py:function:: calc_mean_fourier_amplitude(data: xarray.DataArray, time_dim: str = 'time', lower: float = None, upper: float = None) -> xarray.DataArray

   Calculate mean Fourier amplitude between specified **period** bounds.

   Parameters
   ----------
   data : :py:class:`xarray.DataArray <xarray.DataArray>`
       Input data with at least a time dimension. Can have additional dimensions.
   time_dim : :py:class:`str <str>`, optional
       Name of the time coordinate in the DataArray, by default ``'time'``
   lower : :py:class:`float <float>`
       Lower bound of **period** range
   upper : :py:class:`float <float>`
       Upper bound of **period** range

   Returns
   -------
   :py:class:`xarray.DataArray <xarray.DataArray>`
       Mean amplitude within specified **period** range,
       with same dimensions as input but without time dimension.

   .. tip::

       - The bounds should be in the same units as returned by :py:func:`easyclimate.filter.spectrum.calc_time_spectrum <easyclimate.filter.spectrum.calc_time_spectrum>`
       - Uses :py:func:`easyclimate.filter.spectrum.calc_time_spectrum <easyclimate.filter.spectrum.calc_time_spectrum>` internally
       - The amplitude is scaled by the variance of the input data


.. py:function:: filter_fourier_harmonic_analysis(da: xarray.DataArray, time_dim: str = 'time', period_bounds: tuple = (None, None), filter_type: Literal['highpass', 'lowpass', 'bandpass'] = 'bandpass', sampling_interval: float = 1.0, apply_window: bool = True) -> xarray.DataArray

   Apply Fourier harmonic analysis to filter an dataset along a time dimension.

   Parameters
   ----------
   da : :py:class:`xarray.DataArray <xarray.DataArray>`
       Input data array with a time dimension (e.g., z200 with dims [time, lat, lon]).
   time_dim : :py:class:`str <str>`, optional
       Name of the time dimension (default: 'time').
   period_bounds : :py:class:`tuple <tuple>`, optional
       Period range for filtering in units of sampling_interval (e.g., years if sampling_interval=1).
       Format: ``(min_period, max_period)``. Use None for unbounded limits.

       - High-pass: ``(None, max_period)`` to retain periods < max_period.
       - Low-pass: ``(min_period, None)`` to retain periods > min_period.
       - Bandpass: ``(min_period, max_period)`` to retain min_period < periods < max_period.

   filter_type : :py:class:`str <str>`, optional
       Type of filter: ``'highpass'``, ``'lowpass'``, or ``'bandpass'`` (default: ``'bandpass'``).
   sampling_interval : :py:class:`float <float>`, optional
       Sampling interval of the time dimension (default: 1.0, e.g., 1 year for annual data).
   apply_window : :py:class:`bool <bool>`, optional
       Apply a Hann window to reduce boundary effects (default: True).

   Returns
   -------
   :py:class:`xarray.DataArray <xarray.DataArray>`
       Filtered data array with same dimensions and coordinates as input.

   Examples
   --------
   >>> # Create example data
   >>> ds = xr.DataArray(
   ...    np.random.randn(56, 90, 180),
   ...    dims=['time', 'lat', 'lon'],
   ...    coords={
   ...        'time': np.arange(1948, 2004),
   ...        'lat': np.linspace(-90, 90, 90),
   ...        'lon': np.linspace(0, 360, 180, endpoint=False)
   ...    },
   ...    name='z200'
   ... )
   >>> # Apply low-pass filter to retain periods > 8 years
   >>> ds_filtered = filter_fourier_harmonic_analysis(
   ...    da=ds,
   ...    time_dim='time',
   ...    period_bounds=(8.0, None),
   ...    filter_type='lowpass',
   ...    sampling_interval=1.0,
   ...    apply_window=True
   ... )



