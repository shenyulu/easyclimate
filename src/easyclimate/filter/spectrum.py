"""
Spatio-temporal Spectrum Analysis
"""

from __future__ import annotations
import numpy as np
import xarray as xr
from scipy.fft import fft, ifft, fftfreq
from scipy.signal.windows import hann
from typing import Literal
from ..core.datanode import DataNode

__all__ = [
    "calc_time_spectrum",
    "calc_mean_fourier_amplitude",
    "filter_fourier_harmonic_analysis",
]


def calc_time_spectrum(
    data: xr.DataArray, time_dim: str = "time", inv: bool = False
) -> tuple[xr.DataArray, xr.DataArray]:
    """
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
        zero-frequency component. The amplitude is calculated as the power spectrum (:math:`|\\mathrm{fft}|^2`).
    """
    # Get time values from the DataArray
    ntime = data[time_dim].shape[0]
    time = np.arange(ntime)
    dt = float(time[1] - time[0])  # time step

    def _spectrum_1d(data_1d: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Helper function to compute spectrum for 1D array"""
        # Compute FFT
        if inv:
            datafft = np.fft.fft(data_1d)
        else:
            datafft = np.fft.ifft(data_1d)

        # Calculate power spectrum
        amp_val = np.abs(datafft * np.conj(datafft))

        # Determine frequency cutoff (handle even/odd lengths)
        n = len(time)
        fcut = n // 2 if n % 2 == 0 else (n - 1) // 2

        # Get frequencies (excluding zero frequency)
        freq = np.fft.fftfreq(n, d=dt)[1:fcut]
        # convert to period
        # period = 1.0 / freq
        return amp_val[1:fcut], freq

    # Apply along time dimension using xarray's apply_ufunc
    amp, freq = xr.apply_ufunc(
        _spectrum_1d,
        data,
        input_core_dims=[[time_dim]],
        output_core_dims=[["freq"], ["freq"]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[float, float],
    )

    # Get the frequency values (same for all positions)
    freq_values = freq.isel({d: 0 for d in freq.dims if d != "freq"}).values
    period_values = 1.0 / freq_values

    # Create new DataArray with proper coordinates
    node = DataNode(name="root")
    node["freq"] = xr.DataArray(
        freq_values, dims=["freq"], coords={"freq": freq_values}, name="frequency"
    )
    node["spectrum_freq"] = amp.assign_coords({"freq": freq_values})
    node["period"] = xr.DataArray(
        period_values, dims=["period"], coords={"period": period_values}, name="period"
    )
    node["spectrum_period"] = amp.rename({"freq": "period"}).assign_coords(
        {"period": period_values}
    )

    return node


def calc_mean_fourier_amplitude(
    data: xr.DataArray, time_dim: str = "time", lower: float = None, upper: float = None
) -> xr.DataArray:
    """
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
    """
    # Calculate spectrum (using wavelength by default)
    result = calc_time_spectrum(data, time_dim=time_dim)
    spec = result["spectrum_period"]
    freq = result["period"]

    # Create selection mask
    select = np.logical_and(freq >= lower, freq <= upper)

    # Calculate mean amplitude in selected range
    with np.errstate(invalid="ignore"):  # Ignore NaN warnings
        # Sum over selected frequencies
        selected_sum = spec.where(select, drop=True).sum(dim="period")
        # Total sum over all frequencies
        total_sum = spec.sum(dim="period")
        # Variance along time dimension
        data_var = data.var(dim=time_dim)

        # Calculate scaled amplitude
        amp_scaled = (selected_sum / total_sum) * np.sqrt(data_var)

    # Preserve attributes and name
    amp_scaled.attrs.update(data.attrs)
    amp_scaled.name = data.name or "mean_amplitude"

    return amp_scaled


def filter_fourier_harmonic_analysis(
    da: xr.DataArray,
    time_dim: str = "time",
    period_bounds: tuple = (None, None),
    filter_type: Literal["highpass", "lowpass", "bandpass"] = "bandpass",
    sampling_interval: float = 1.0,
    apply_window: bool = True,
) -> xr.DataArray:
    """
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

    """
    # Input validation
    if time_dim not in da.dims:
        raise ValueError(
            f"Time dimension '{time_dim}' not found in DataArray dims: {da.dims}"
        )
    if filter_type not in ["highpass", "lowpass", "bandpass"]:
        raise ValueError(
            f"Invalid filter_type: {filter_type}. Choose 'highpass', 'lowpass', or 'bandpass'."
        )

    min_period, max_period = period_bounds
    if filter_type == "highpass" and max_period is None:
        raise ValueError("High-pass filter requires max_period to be specified.")
    if filter_type == "lowpass" and min_period is None:
        raise ValueError("Low-pass filter requires min_period to be specified.")
    if filter_type == "bandpass" and (min_period is None or max_period is None):
        raise ValueError("Bandpass filter requires both min_period and max_period.")
    if min_period is not None and max_period is not None and min_period >= max_period:
        raise ValueError(
            f"min_period ({min_period}) must be less than max_period ({max_period})."
        )

    # Check for NaN values
    if np.any(np.isnan(da.values)):
        raise ValueError(
            "Input DataArray contains NaN values. Please handle missing data first."
        )

    # Get time dimension size
    n = len(da[time_dim])
    time_axis = da.dims.index(time_dim)

    # Apply window function if requested
    if apply_window:
        window = hann(n)
        # Reshape window to broadcast along time axis
        window_shape = [-1 if i == time_axis else 1 for i in range(da.ndim)]
        da_windowed = da * window.reshape(window_shape)
    else:
        da_windowed = da

    # Compute FFT
    z_fft = fft(da_windowed.values, axis=time_axis)

    # Compute frequencies and periods
    freq = fftfreq(n, d=sampling_interval)
    periods = np.where(freq != 0, 1.0 / np.abs(freq), np.inf)

    # Create frequency mask based on filter type
    if filter_type == "highpass":
        mask = (periods < max_period) | (
            freq == 0
        )  # Retain periods < max_period and zero frequency
    elif filter_type == "lowpass":
        mask = (periods > min_period) | (
            freq == 0
        )  # Retain periods > min_period and zero frequency
    else:  # bandpass
        mask = ((periods > min_period) & (periods < max_period)) | (freq == 0)

    # Reshape mask to broadcast along time axis
    mask_shape = [n if i == time_axis else 1 for i in range(da.ndim)]
    mask = mask.reshape(mask_shape)

    # Apply mask to FFT coefficients
    z_fft_filtered = z_fft.copy()
    z_fft_filtered = z_fft * mask

    # Inverse FFT
    z_filtered = ifft(z_fft_filtered, axis=time_axis).real

    # Create output DataArray
    da_filtered = xr.DataArray(
        z_filtered, dims=da.dims, coords=da.coords, name=da.name, attrs=da.attrs
    )

    # Add filter metadata to attributes
    da_filtered.attrs["filter_type"] = filter_type
    da_filtered.attrs["period_bounds"] = period_bounds
    da_filtered.attrs["sampling_interval"] = sampling_interval
    da_filtered.attrs["window_applied"] = apply_window

    return da_filtered
