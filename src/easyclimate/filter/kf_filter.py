"""
Extract equatorial waves by filtering in the Wheeler-Kiladis wavenumber-frequency domain

.. seealso::
    - Kiladis, G. N., M. C. Wheeler, P. T. Haertel, K. H. Straub, and P. E. Roundy (2009), Convectively coupled equatorial waves, Rev. Geophys., 47, RG2003, doi: https://doi.org/10.1029/2008RG000266.
    - https://k3.cicsnc.org/carl/monitor
"""

import numpy as np
import xarray as xr
from scipy import signal
from typing import Literal, Optional

__all__ = [
    "kf_filter_wheeler_and_kiladis_1999",
    "kf_filter_lf_wave",
    "kf_filter_mjo_wave",
    "kf_filter_er_wave",
    "kf_filter_kelvin_wave",
    "kf_filter_mt_wave",
    "kf_filter_mrg_wave",
    "kf_filter_td_wave",
]


def kf_filter_wheeler_and_kiladis_1999(
    input_data: xr.DataArray,
    steps_per_day: int | float,
    tMin: float | Literal[-999],
    tMax: float | Literal[-999],
    kMin: int | float | Literal[-999],
    kMax: int | float | Literal[-999],
    hMin: float | Literal[-999],
    hMax: float | Literal[-999],
    waveName: Literal["kelvin", "er", "mrg", "ig0", "ig1", "ig2", None] = None,
    time_dim="time",
    lon_dim="lon",
) -> xr.DataArray:
    """
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
    """
    mis = -999

    def _kf_filter(input):
        timeDim, lonDim = input.shape

        # Reshape data with lon-dim reversed and transposed [time,lon] => [lon,time]
        originalData = np.flip(input, axis=1).T

        # Detrend along time-axis (axis=1)
        detrendData = signal.detrend(originalData, axis=1)

        # Create Tukey taper window and apply
        taper = signal.windows.tukey(originalData.shape[1], 0.05, True)
        taperData = detrendData * taper

        # Perform 2-D Fourier Transform
        fftData = np.fft.rfft2(taperData)
        kDim = lonDim
        freqDim = np.round(fftData.shape[1])

        # Find the indeces for the period cut-offs
        jMin = int(np.round((timeDim * 1.0 / (tMax * steps_per_day)), 0))
        jMax = int(np.round((timeDim * 1.0 / (tMin * steps_per_day)), 0))
        jMax = np.min((jMax, freqDim))

        # Find the indices for the wavenumber cut-offs
        # This is more complicated because east and west are separate
        if kMin < 0:
            iMin = np.round((kDim + kMin), 3)
            iMin = np.max((iMin, (kDim / 2)))
        else:
            iMin = np.round(kMin, 3)
            iMin = np.min((iMin, (kDim / 2)))

        if kMax < 0:
            iMax = np.round((kDim + kMax), 3)
            iMax = np.max((iMax, (kDim / 2)))
        else:
            iMax = np.round(kMax, 3)
            iMax = np.min((iMax, (kDim / 2)))

        # set the appropriate coefficients to zero
        iMin = int(iMin)
        iMax = int(iMax)
        jMin = int(jMin)
        jMax = int(jMax)
        if jMin > 0:
            fftData[:, : jMin - 1] = 0
        if jMax < (freqDim - 1):
            fftData[:, jMax + 1 :] = 0

        if iMin < iMax:
            # Set things outside the range to zero, this is more normal
            if iMin > 0:
                fftData[: iMin - 1, :] = 0
            if iMax < (kDim - 1):
                fftData[iMax + 1 :, :] = 0
        else:
            # Set things inside the range to zero, this should be somewhat unusual
            fftData[iMax + 1 : iMin - 1, :] = 0

        # Find constants
        beta = 2.28e-11
        if hMin != -999:
            cMin = float(9.8 * float(hMin)) ** 0.5
        else:
            cMin = hMin
        if hMax != -999:
            cMax = float(9.8 * float(hMax)) ** 0.5
        else:
            cMax = hMax
        c = np.array([cMin, cMax])
        spc = 24 * 3600.0 / (2 * np.pi * steps_per_day)  # seconds per cycle

        # Now set things to zero that are outside the Kelvin dispersion
        for i in range(0, kDim):
            # Find Non-Dimensional WaveNumber (k)
            if i > (kDim / 2):
                # k is negative
                k = (i - kDim) * 1.0 / (6.37e6)  # adjusting for circumfrence of earth
            else:
                # k is positive
                k = i * 1.0 / (6.37e6)  # adjusting for circumfrence of earth

            # Find Frequency
            freq = np.array([0, freqDim * (1.0 / spc)])  # waveName='None'
            jMinWave = 0
            jMaxWave = freqDim

            if waveName == "kelvin":
                freq = k * c
            if waveName == "er":
                freq = -beta * k / (k**2 + 3.0 * beta / c)
            if waveName == "ig1":
                freq = (3 * beta * c + k**2 * c**2) ** 0.5
            if waveName == "ig2":
                freq = (5 * beta * c + k**2 * c**2) ** 0.5
            if waveName == "mrg" or waveName == "ig0":
                if k == 0:
                    freq = (beta * c) ** 0.5
                else:
                    if k > 0:
                        freq = k * c * (0.5 + 0.5 * (1 + 4 * beta / (k**2 * c)) ** 0.5)
                    else:
                        freq = k * c * (0.5 - 0.5 * (1 + 4 * beta / (k**2 * c)) ** 0.5)

            # Get Min/Max Wave
            if hMin == mis:
                jMinWave = 0
            else:
                jMinWave = int(np.floor(freq[0] * spc * timeDim))

            if hMax == mis:
                jMaxWave = freqDim
            else:
                jMaxWave = int(np.ceil(freq[1] * spc * timeDim))

            jMaxWave = np.max([jMaxWave, 0])
            jMinWave = np.min([jMinWave, freqDim])

            # set the appropriate coefficients to zero
            i = int(i)
            jMinWave = int(jMinWave)
            jMaxWave = int(jMaxWave)
            if jMinWave > 0:
                fftData[i, : jMinWave - 1] = 0
            if jMaxWave < (freqDim - 1):
                fftData[i, jMaxWave + 1 :] = 0

        # perform the inverse transform to reconstruct the data
        returnedData = np.fft.irfft2(fftData)

        # Reshape data from [lon,time] to [time,lon]
        outData = np.flip(returnedData, axis=0).T

        return outData

    result = xr.apply_ufunc(
        _kf_filter,
        input_data,
        input_core_dims=[[time_dim, lon_dim]],
        output_core_dims=[[time_dim, lon_dim]],
        vectorize=True,
    )
    result.attrs = {}
    result.attrs["wavenumber"] = [kMin, kMax]
    result.attrs["period"] = [tMin, tMax]
    result.attrs["depth"] = f"[{hMin}, {hMax}] m"
    result.attrs["waveName"] = str(waveName)
    return result


def kf_filter_lf_wave(
    input_data: xr.DataArray,
    steps_per_day: int | float,
    tMin: float | Literal[-999] = 120,
    tMax: float | Literal[-999] = -999,
    kMin: int | float | Literal[-999] = -999,
    kMax: int | float | Literal[-999] = -999,
    hMin: float | Literal[-999] = -999,
    hMax: float | Literal[-999] = -999,
    waveName: Literal["kelvin", "er", "mrg", "ig0", "ig1", "ig2", None] = None,
    time_dim="time",
    lon_dim="lon",
) -> xr.DataArray:
    """
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
    """
    result = kf_filter_wheeler_and_kiladis_1999(
        input_data,
        steps_per_day,
        tMin,
        tMax,
        kMin,
        kMax,
        hMin,
        hMax,
        waveName,
        time_dim,
        lon_dim,
    )
    return result


def kf_filter_mjo_wave(
    input_data: xr.DataArray,
    steps_per_day: int | float,
    tMin: float | Literal[-999] = 30,
    tMax: float | Literal[-999] = 96,
    kMin: int | float | Literal[-999] = 1,
    kMax: int | float | Literal[-999] = 5,
    hMin: float | Literal[-999] = -999,
    hMax: float | Literal[-999] = -999,
    waveName: Literal["kelvin", "er", "mrg", "ig0", "ig1", "ig2", None] = None,
    time_dim="time",
    lon_dim="lon",
) -> xr.DataArray:
    """
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
    """
    result = kf_filter_wheeler_and_kiladis_1999(
        input_data,
        steps_per_day,
        tMin,
        tMax,
        kMin,
        kMax,
        hMin,
        hMax,
        waveName,
        time_dim,
        lon_dim,
    )
    return result


def kf_filter_er_wave(
    input_data: xr.DataArray,
    steps_per_day: int | float,
    tMin: float | Literal[-999] = 9.7,
    tMax: float | Literal[-999] = 48,
    kMin: int | float | Literal[-999] = -10,
    kMax: int | float | Literal[-999] = -1,
    hMin: float | Literal[-999] = 8,
    hMax: float | Literal[-999] = 90,
    waveName: Literal["kelvin", "er", "mrg", "ig0", "ig1", "ig2", None] = "er",
    time_dim="time",
    lon_dim="lon",
) -> xr.DataArray:
    """
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
    """
    result = kf_filter_wheeler_and_kiladis_1999(
        input_data,
        steps_per_day,
        tMin,
        tMax,
        kMin,
        kMax,
        hMin,
        hMax,
        waveName,
        time_dim,
        lon_dim,
    )
    return result


def kf_filter_kelvin_wave(
    input_data: xr.DataArray,
    steps_per_day: int | float,
    tMin: float | Literal[-999] = 2.5,
    tMax: float | Literal[-999] = 30,
    kMin: int | float | Literal[-999] = 1,
    kMax: int | float | Literal[-999] = 14,
    hMin: float | Literal[-999] = 8,
    hMax: float | Literal[-999] = 90,
    waveName: Literal["kelvin", "er", "mrg", "ig0", "ig1", "ig2", None] = "kelvin",
    time_dim="time",
    lon_dim="lon",
) -> xr.DataArray:
    """
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
    """
    result = kf_filter_wheeler_and_kiladis_1999(
        input_data,
        steps_per_day,
        tMin,
        tMax,
        kMin,
        kMax,
        hMin,
        hMax,
        waveName,
        time_dim,
        lon_dim,
    )
    return result


def kf_filter_mt_wave(
    input_data: xr.DataArray,
    steps_per_day: int | float,
    tMin: float | Literal[-999] = 2.5,
    tMax: float | Literal[-999] = 10,
    kMin: int | float | Literal[-999] = -14,
    kMax: int | float | Literal[-999] = 0,
    hMin: float | Literal[-999] = -999.0,
    hMax: float | Literal[-999] = -999.0,
    waveName: Literal["kelvin", "er", "mrg", "ig0", "ig1", "ig2", None] = None,
    time_dim="time",
    lon_dim="lon",
) -> xr.DataArray:
    """
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
    """
    result = kf_filter_wheeler_and_kiladis_1999(
        input_data,
        steps_per_day,
        tMin,
        tMax,
        kMin,
        kMax,
        hMin,
        hMax,
        waveName,
        time_dim,
        lon_dim,
    )
    return result


def kf_filter_mrg_wave(
    input_data: xr.DataArray,
    steps_per_day: int | float,
    tMin: float | Literal[-999] = 3,
    tMax: float | Literal[-999] = 9.6,
    kMin: int | float | Literal[-999] = -10,
    kMax: int | float | Literal[-999] = -1,
    hMin: float | Literal[-999] = 8,
    hMax: float | Literal[-999] = 90,
    waveName: Literal["kelvin", "er", "mrg", "ig0", "ig1", "ig2", None] = "mrg",
    time_dim="time",
    lon_dim="lon",
) -> xr.DataArray:
    """
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
    """
    result = kf_filter_wheeler_and_kiladis_1999(
        input_data,
        steps_per_day,
        tMin,
        tMax,
        kMin,
        kMax,
        hMin,
        hMax,
        waveName,
        time_dim,
        lon_dim,
    )
    return result


def kf_filter_td_wave(
    input_data: xr.DataArray,
    steps_per_day: int | float,
    tMin: float | Literal[-999] = 2,
    tMax: float | Literal[-999] = 8.5,
    kMin: int | float | Literal[-999] = -15,
    kMax: int | float | Literal[-999] = -6,
    hMin: float | Literal[-999] = 90,
    hMax: float | Literal[-999] = -999,
    waveName: Literal["kelvin", "er", "mrg", "ig0", "ig1", "ig2", None] = "mrg",
    time_dim="time",
    lon_dim="lon",
) -> xr.DataArray:
    """
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
    """
    result = kf_filter_wheeler_and_kiladis_1999(
        input_data,
        steps_per_day,
        tMin,
        tMax,
        kMin,
        kMax,
        hMin,
        hMax,
        waveName,
        time_dim,
        lon_dim,
    )
    return result
