"""
Butterworth bandpass filter
"""

from __future__ import annotations
import numpy as np
import xarray as xr
import scipy.signal as signal
from ..core.utility import find_dims_axis, generate_dataset_dispatcher

__all__ = ["calc_butter_bandpass", "calc_butter_lowpass", "calc_butter_highpass"]


@generate_dataset_dispatcher
def calc_butter_bandpass(
    data: xr.DataArray | xr.Dataset,
    sampling_frequency: int,
    period: list[int],
    N: int = 3,
    dim: str = "time",
) -> xr.DataArray:
    """Butterworth bandpass filter.

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
    """
    time_axis = find_dims_axis(data, dim=dim)
    Wn_value = 2 * (1 / np.array(period)) * sampling_frequency
    sos = signal.butter(N=N, Wn=np.sort(Wn_value), btype="bandpass", output="sos")
    filter_data = signal.sosfiltfilt(sos, data.data, axis=time_axis)
    return data.copy(data=filter_data, deep=True)


@generate_dataset_dispatcher
def calc_butter_lowpass(
    data: xr.DataArray | xr.Dataset,
    sampling_frequency: int,
    period: float,
    N: int = 3,
    dim: str = "time",
) -> xr.DataArray:
    """Butterworth lowpass filter.

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
    """
    time_axis = find_dims_axis(data, dim=dim)
    Wn_value = 2 * (1 / np.array(period)) * sampling_frequency
    sos = signal.butter(N=N, Wn=Wn_value, btype="lowpass", output="sos")
    filter_data = signal.sosfiltfilt(sos, data.data, axis=time_axis)
    return data.copy(data=filter_data, deep=True)


@generate_dataset_dispatcher
def calc_butter_highpass(
    data: xr.DataArray | xr.Dataset,
    sampling_frequency: int,
    period: float,
    N: int = 3,
    dim: str = "time",
) -> xr.DataArray:
    """Butterworth highpass filter.

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
    """
    time_axis = find_dims_axis(data, dim=dim)
    Wn_value = 2 * (1 / np.array(period)) * sampling_frequency
    sos = signal.butter(N=N, Wn=Wn_value, btype="highpass", output="sos")
    filter_data = signal.sosfiltfilt(sos, data.data, axis=time_axis)
    return data.copy(data=filter_data, deep=True)
