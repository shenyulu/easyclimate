easyclimate.filter.spectrum
===========================

.. py:module:: easyclimate.filter.spectrum

.. autoapi-nested-parse::

   Spatio-temporal Spectrum Analysis



Functions
---------

.. autoapisummary::

   easyclimate.filter.spectrum.calc_time_spectrum
   easyclimate.filter.spectrum.calc_mean_fourier_amplitude
   easyclimate.filter.spectrum.filter_fourier_harmonic_analysis


Module Contents
---------------

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



