"""
Lanczos filter
"""

from __future__ import annotations
import numpy as np
import xarray as xr
import scipy.signal as signal
from ..core.utility import find_dims_axis, generate_dataset_dispatcher
from typing import Literal

__all__ = ["calc_lanczos_bandpass", "calc_lanczos_lowpass", "calc_lanczos_highpass"]


def lanczos_lowpass_weights(window, cutoff):
    """Calculate weights for a low pass Lanczos filter.

    Args:

    window: int
        The length of the filter window.

    cutoff: float
        The cutoff frequency in inverse time steps.

    """
    order = ((window - 1) // 2) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1.0, n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2.0 * np.pi * cutoff * k) / (np.pi * k)
    w[n - 1 : 0 : -1] = firstfactor * sigma
    w[n + 1 : -1] = firstfactor * sigma
    low_pass = w[1:-1]
    return xr.DataArray(low_pass, dims=["window"])


def lanczos_highpass_weights(window, cutoff):
    """Calculate weights for a high pass Lanczos filter.

    Args:
    window: int
        The length of the filter window.
    cutoff: float
        The cutoff frequency in inverse time steps.

    .. seealso::
        https://wyhtsai.github.io/pyaos-wks/docs/A3_lanczos_filter.html
    """
    # # Get low-pass filter weights
    # low_pass = lanczos_lowpass_weights(window, cutoff)
    # low_pass = low_pass.data
    # # Create a unit impulse
    # impulse = np.zeros_like(low_pass)
    # impulse[len(impulse) // 2] = 1
    # # Calculate high-pass filter weights
    # high_pass = impulse - low_pass
    # return xr.DataArray(high_pass, dims=["window"])
    order = ((window - 1) // 2) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 1 - 2 * cutoff
    k = np.arange(1.0, n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2.0 * np.pi * cutoff * k) / (np.pi * k)
    w[n - 1 : 0 : -1] = -firstfactor * sigma
    w[n + 1 : -1] = -firstfactor * sigma
    high_pass = w[1:-1]
    return xr.DataArray(high_pass, dims=["window"])


def apply_pass_multidim(data, fw, dim):
    """
    Apply filter along the time dimension for multi-dimensional data.

    Args:
        data: xarray.DataArray
            Input data (can be multi-dimensional).
        fw: xarray.DataArray
            Filter weights.

    Returns:
        xarray.DataArray:
            Filtered data along the time dimension.
    """

    # Define a 1D filtering function
    def lowpass_1d(signal, weights):
        return np.convolve(signal, weights, mode="same")

    # Use apply_ufunc to apply the 1D filter along the `dim` dimension
    return xr.apply_ufunc(
        lowpass_1d,
        data,
        fw,
        input_core_dims=[[dim], ["window"]],  # Specify dimensions to apply filter
        output_core_dims=[[dim]],
        vectorize=True,  # Apply along each slice independently
        dask="parallelized",  # Enable parallel computation if using Dask
        dask_gufunc_kwargs={"allow_rechunk": True},
        output_dtypes=[data.dtype],  # Ensure correct output type
    )


def calc_lanczos_lowpass(
    data: xr.DataArray | xr.Dataset,
    window_length: int,
    period: int,
    dim: str = "time",
    method: Literal["rolling", "convolve"] = "rolling",
) -> xr.DataArray:
    """
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
    """
    fw = lanczos_lowpass_weights(window_length, 1.0 / period)

    if method == "rolling":
        lowpass_hf = (
            data.rolling({dim: len(fw)}, center=True).construct("window").dot(fw)
        )
    elif method == "convolve":
        lowpass_hf = apply_pass_multidim(data, fw, dim)

    lowpass_hf.attrs = {}
    return lowpass_hf


def calc_lanczos_highpass(
    data: xr.DataArray | xr.Dataset,
    window_length: int,
    period: int,
    dim: str = "time",
    method: Literal["rolling", "convolve"] = "rolling",
) -> xr.DataArray:
    """
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
    """
    fw = lanczos_highpass_weights(window_length, 1.0 / period)

    if method == "rolling":
        highpass_hf = (
            data.rolling({dim: len(fw)}, center=True).construct("window").dot(fw)
        )
    elif method == "convolve":
        highpass_hf = apply_pass_multidim(data, fw, dim)

    highpass_hf.attrs = {}
    return highpass_hf


def calc_lanczos_bandpass(
    data: xr.DataArray | xr.Dataset,
    window_length: int,
    period: list[int],
    dim: str = "time",
    method: Literal["rolling", "convolve"] = "rolling",
) -> xr.DataArray:
    """
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
    """
    period = np.array(period)
    period = np.sort(period)

    fw_0 = lanczos_highpass_weights(window_length, 1.0 / period[0])
    fw_1 = lanczos_highpass_weights(window_length, 1.0 / period[1])

    if method == "rolling":
        lowpass_hf0 = (
            data.rolling({dim: len(fw_0)}, center=True).construct("window").dot(fw_0)
        )
        lowpass_hf1 = (
            data.rolling({dim: len(fw_1)}, center=True).construct("window").dot(fw_1)
        )
        bandpass = lowpass_hf0 - lowpass_hf1
    elif method == "convolve":
        lowpass_hf0 = apply_pass_multidim(data, fw_0, dim)
        lowpass_hf1 = apply_pass_multidim(data, fw_1, dim)
        bandpass = lowpass_hf0 - lowpass_hf1

    bandpass.attrs = {}
    return bandpass
