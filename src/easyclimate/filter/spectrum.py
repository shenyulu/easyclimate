"""
Spatio-temporal Spectrum Analysis
"""

from __future__ import annotations
import numpy as np
import xarray as xr
from ..core.datanode import DataNode

__all__ = ["calc_time_spectrum", "calc_mean_fourier_amplitude"]


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
