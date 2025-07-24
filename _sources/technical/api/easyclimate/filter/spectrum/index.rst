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


