"""
Wheeler-Kiladis Space-Time Spectra

This module provides functions for analyzing and visualizing Wheeler-Kiladis space-time spectra,
including signal processing, symmetric/asymmetric decomposition, spectral analysis, and plotting
of equatorial wave dispersion relationships.

.. seealso::

    - https://github.com/mmaiergerber/wk_spectra
    - Wheeler, M., & Kiladis, G. N. (1999). Convectively Coupled Equatorial Waves: Analysis of Clouds and Temperature in the Wavenumber–Frequency Domain. Journal of the Atmospheric Sciences, 56(3), 374-399. https://journals.ametsoc.org/view/journals/atsc/56/3/1520-0469_1999_056_0374_ccewao_2.0.co_2.xml
    - Kiladis, G. N., M. C. Wheeler, P. T. Haertel, K. H. Straub, and P. E. Roundy (2009), Convectively coupled equatorial waves, Rev. Geophys., 47, RG2003, doi:https://doi.org/10.1029/2008RG000266
    - Wheeler, M. C., & Nguyen, H. (2015). TROPICAL METEOROLOGY AND CLIMATE | Equatorial Waves. In Encyclopedia of Atmospheric Sciences (pp. 102–112). Elsevier. https://doi.org/10.1016/B978-0-12-382225-3.00414-X
    - Yoshikazu Hayashi, A Generalized Method of Resolving Disturbances into Progressive and Retrogressive Waves by Space Fourier and Time Cross-Spectral Analyses, Journal of the Meteorological Society of Japan. Ser. II, 1971, Volume 49, Issue 2, Pages 125-128, Released on J-STAGE May 27, 2008, Online ISSN 2186-9057, Print ISSN 0026-1165, https://doi.org/10.2151/jmsj1965.49.2_125, https://www.jstage.jst.go.jp/article/jmsj1965/49/2/49_2_125/_article/-char/en
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from ...backend import wk_analysis, matsuno_plot

__all__ = [
    "remove_dominant_signals",
    "decompose_symasym",
    "calc_spectral_coefficients",
    "draw_wk_anti_analysis",
    "draw_wk_sym_analysis",
]


def remove_dominant_signals(
    data: xr.DataArray,
    spd: float,
    nDayWin: float,
    nDaySkip: float,
    time_dim: str = "time",
    lon_dim: str = "lon",
    lat_dim: str = "lat",
) -> xr.DataArray:
    """
    Removes the dominant signals by removing the long term linear trend (conserving the mean) and
    by eliminating the annual cycle by removing all time periods less than a corresponding critical frequency.

    Parameters
    ----------
    data : :py:class:`xarray.DataArray<xarray.DataArray>`
        Input data array with time, lat, lon dimensions

    .. caution:: `data` must be **daily** data, and the length of time should be longer than one year (i.e., >365 day).

    spd : :py:class:`float <float>`
        Samples per day (frequency of observations)
    nDayWin : :py:class:`float <float>`
        Number of days in each analysis window
    nDaySkip : :py:class:`float <float>`
        Number of days to skip between windows
    time_dim : :py:class:`str <str>`, optional
        Name of time dimension, default 'time'
    lon_dim : :py:class:`str <str>`, optional
        Name of longitude dimension, default 'lon'
    lat_dim : :py:class:`str <str>`, optional
        Name of latitude dimension, default 'lat'

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`
        Data with dominant signals removed

    Example
    -------
    >>> data = xr.DataArray(np.random.rand(730, 10, 20),
    ...                    dims=['time', 'lat', 'lon'])
    >>> cleaned = remove_dominant_signals(data, spd=1, nDayWin=90, nDaySkip=30)

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_wk_spectra.py
    """
    data = data.transpose(time_dim, lat_dim, lon_dim)
    result_array = wk_analysis.remove_dominant_signals(
        data.data, spd, nDayWin, nDaySkip
    )
    return data.copy(deep=True, data=result_array)


def decompose_symasym(da, lat_dim="lat"):
    """
    Decompose data into symmetric and asymmetric parts about the equator.

    The symmetric part is stored in the Southern Hemisphere and the asymmetric
    part in the Northern Hemisphere.

    Parameters
    ----------
    da : :py:class:`xarray.DataArray<xarray.DataArray>`
        Input array with latitude dimension
    lat_dim : :py:class:`str <str>`, optional
        Name of latitude dimension, default 'lat'

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`
        Array with decomposed components (symmetric in SH, asymmetric in NH)

    Example
    -------
    >>> data = xr.DataArray(np.random.rand(10, 20), dims=['lat', 'lon'])
    >>> decomposed = decompose_symasym(data)

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_wk_spectra.py
    """

    def sym_asym_1d(array_lat):
        """Decompose a 1D array into symmetric and asymmetric parts about the equator."""
        nlat = len(array_lat)
        N = nlat // 2
        SymAsym = np.copy(array_lat)
        for i in range(N):
            SymAsym[i] = 0.5 * (
                array_lat[nlat - 1 - i] + array_lat[i]
            )  # Symmetric in South
            SymAsym[nlat - 1 - i] = 0.5 * (
                array_lat[nlat - 1 - i] - array_lat[i]
            )  # Asymmetric in North
        return SymAsym

    if lat_dim not in da.dims:
        raise ValueError(f"Input DataArray must have {lat_dim} dimension")
    return xr.apply_ufunc(
        sym_asym_1d,
        da,
        input_core_dims=[[lat_dim]],
        output_core_dims=[[lat_dim]],
        vectorize=True,
    )


def calc_spectral_coefficients(
    data: xr.DataArray,
    spd: float,
    nDayWin: float,
    nDaySkip: float,
    time_dim: str = "time",
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    max_freq: float = 0.5,
    max_wn: float = 15,
):
    """
    Calculate Wheeler-Kiladis spectral coefficients.

    Parameters
    ----------
    data : :py:class:`xarray.DataArray<xarray.DataArray>`
        Input data array.

    .. caution:: `data` must be **daily** data, and the length of time should be longer than one year (i.e., >365 day).

    spd : :py:class:`float <float>`
        Samples per day
    nDayWin : :py:class:`float <float>`
        Number of days in analysis window
    nDaySkip : :py:class:`float <float>`
        Number of days to skip between windows
    time_dim : :py:class:`str <str>`, optional
        Time dimension name, default 'time'
    lon_dim : :py:class:`str <str>`, optional
        Longitude dimension name, default 'lon'
    lat_dim : :py:class:`str <str>`, optional
        Latitude dimension name, default 'lat'
    max_freq : :py:class:`float <float>`, optional
        Maximum frequency to return (CPD), default 0.5
    max_wn : :py:class:`float <float>`, optional
        Maximum wavenumber to return, default 15

    Returns
    -------
    :py:class:`xarray.Dataset<xarray.Dataset>`
        Dataset containing:

        - ``psumanti_r``: Antisymmetric power spectrum (background removed)
        - ``psumsym_r``: Symmetric power spectrum (background removed)

    Example
    -------
    >>> data = xr.DataArray(np.random.rand(730, 10, 20),
    ...                    dims=['time', 'lat', 'lon'])
    >>> spectra = calc_spectral_coefficients(data, spd=1, nDayWin=90, nDaySkip=30)

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_wk_spectra.py
    """
    data = data.transpose(time_dim, lat_dim, lon_dim)

    ntim, nlat, nlon, nDayTot, nSampWin, nSampSkip = wk_analysis.sampling_vars(
        data.data, spd, nDayWin, nDaySkip
    )
    wavefft, freqfft, peeAS = wk_analysis.spectral_coefficients(
        data.data, spd, nDayWin, nDaySkip
    )

    # Calculate the power spectrum
    power = (np.abs(peeAS)) ** 2  # power[window,freq,lat,lon]

    psumanti, psumsym = wk_analysis.separate_power(
        power, nlat, nSampWin, wavefft, freqfft
    )

    # Now derive the background spectrum (red noise)
    # Sum power over all latitude
    # Apply smoothing to the spectrum. This smoothing DOES include wavenumber zero.
    psumb = wk_analysis.derive_background(power, nlat, nSampWin, wavefft, freqfft)

    # Cropping the output
    indwave = np.where(np.logical_and(wavefft >= -max_wn, wavefft <= max_wn))[0]
    indfreq = np.where(np.logical_and(freqfft > 0, freqfft <= max_freq))[0]

    wavefft = wavefft[indwave]
    freqfft = freqfft[indfreq]

    psumanti = psumanti[indfreq, :]
    psumanti = psumanti[:, indwave]

    psumsym = psumsym[indfreq, :]
    psumsym = psumsym[:, indwave]

    psumb = psumb[indfreq, :]
    psumb = psumb[:, indwave]

    psumanti_r = psumanti / psumb
    psumsym_r = psumsym / psumb

    #
    indwave = np.where(
        np.logical_and(wavefft >= -max_wn - 0.1, wavefft <= max_wn + 0.1)
    )[0]
    indfreq = np.where(np.logical_and(freqfft > 0, freqfft <= max_freq + 0.1))[0]

    waveplt = wavefft[indwave]
    freqplt = freqfft[indfreq]

    #
    powerplt = psumanti_r[indfreq, :]
    powerplt = powerplt[:, indwave]

    psumanti_r_xarray = xr.DataArray(
        powerplt, dims=("freq", "k"), coords={"freq": freqplt, "k": waveplt}
    )

    powerplt = psumsym_r[indfreq, :]
    powerplt = powerplt[:, indwave]

    psumsym_r_xarray = xr.DataArray(
        powerplt, dims=("freq", "k"), coords={"freq": freqplt, "k": waveplt}
    )

    result = xr.Dataset()
    result["psumanti_r"] = psumanti_r_xarray
    result["psumsym_r"] = psumsym_r_xarray

    return result


def draw_wk_anti_analysis(
    max_freq: float = 0.5,
    max_wn: float = 15,
    ax=None,
    add_xylabel: bool = True,
    add_central_line: bool = True,
    add_westward_and_eastward: bool = True,
    auto_determine_xyrange: bool = True,
    freq_lines: bool = True,
    matsuno_modes_labels: bool = True,
    cpd_lines_levels: list = [3, 6, 30],
    matsuno_lines: bool = True,
    he: list = [12, 25, 50],
    meridional_modes: list = [1],
):
    """
    Plot antisymmetric Wheeler-Kiladis analysis with Matsuno dispersion curves.

    Parameters
    ----------
    max_freq : :py:class:`float <float>`, optional
        Maximum frequency to plot (CPD), default 0.5
    max_wn : :py:class:`float <float>`, optional
        Maximum wavenumber to plot, default 15
    ax : matplotlib.axes.Axes, optional
        Axes to plot on, creates new if None
    add_xylabel : :py:class:`bool <bool>`, optional
        Add x/y labels, default True
    add_central_line : :py:class:`bool <bool>`, optional
        Add central vertical line, default True
    add_westward_and_eastward : :py:class:`bool <bool>`, optional
        Add eastward/westward labels, default True
    auto_determine_xyrange : :py:class:`bool <bool>`, optional
        Auto-set axis ranges, default True
    freq_lines : :py:class:`bool <bool>`, optional
        Add frequency lines, default True
    matsuno_modes_labels : :py:class:`bool <bool>`, optional
        Add Matsuno mode labels, default True
    cpd_lines_levels : list, optional
        Periods (days) for frequency lines, default [3, 6, 30]
    matsuno_lines : :py:class:`bool <bool>`, optional
        Plot Matsuno dispersion curves, default True
    he : list, optional
        Equivalent depths for Matsuno curves, default [12, 25, 50]
    meridional_modes : list, optional
        Meridional mode numbers, default [1]

    Example
    -------
    >>> fig, ax = plt.subplots()
    >>> draw_wk_anti_analysis(ax=ax)

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_wk_spectra.py
    """
    matsuno_modes = matsuno_plot.matsuno_modes_wk(
        he=he, n=meridional_modes, max_wn=max_wn
    )

    if ax is None:
        ax = plt.gca()

    if add_xylabel:
        ax.set_xlabel("Zonal Wavenumber")
        ax.set_ylabel("Frequency (CPD)")

    if add_central_line:
        ax.axvline(x=0, color="k", linestyle="--")
    if add_westward_and_eastward:
        ax.text(max_wn - 2 * 0.25 * max_wn, -0.08, "Eastward")
        ax.text(-max_wn + 0.25 * max_wn, -0.08, "Westward")
    if auto_determine_xyrange:
        ax.set_xlim((-max_wn, max_wn))
        ax.set_ylim((0.02, max_freq))

    if freq_lines:
        for d in cpd_lines_levels:
            if (1.0 / d) <= max_freq:
                ax.axhline(y=1.0 / d, color="k", linestyle="--")
                ax.text(
                    -max_wn + 0.2,
                    (1.0 / d + 0.01),
                    str(d) + " days",
                    bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": "none"},
                )

    if matsuno_lines:
        for key in matsuno_modes:
            ax.plot(
                matsuno_modes[key]["MRG(he={}m)".format(key)], color="k", linestyle="--"
            )
            ax.plot(
                matsuno_modes[key]["EIG(n=0,he={}m)".format(key)],
                color="k",
                linestyle="--",
            )

    if matsuno_modes_labels:
        key = list(matsuno_modes.keys())[len(list(matsuno_modes.keys())) // 2]
        wn = matsuno_modes[key].index.values

        # Print EIG(n=0) Label
        i = int((len(wn) / 2) + 0.1 * (len(wn) / 2))
        (i,) = np.where(wn == wn[i])[0]
        ax.text(
            wn[i] - 1,
            matsuno_modes[key]["EIG(n=0,he={}m)".format(key)].iloc[i],
            "EIG (n=0)",
            bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": "none"},
        )

        # Print MRG Label
        i = int(0.7 * (len(wn) / 2))
        (i,) = np.where(wn == wn[i])[0]
        ax.text(
            wn[i] - 1,
            matsuno_modes[key]["MRG(he={}m)".format(key)].iloc[i],
            "MRG",
            bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": "none"},
        )


def draw_wk_sym_analysis(
    max_freq: float = 0.5,
    max_wn: float = 15,
    ax=None,
    add_xylabel: bool = True,
    add_central_line: bool = True,
    add_westward_and_eastward: bool = True,
    auto_determine_xyrange: bool = True,
    freq_lines: bool = True,
    matsuno_modes_labels: bool = True,
    cpd_lines_levels: list = [3, 6, 30],
    matsuno_lines: bool = True,
    he: list = [12, 25, 50],
    meridional_modes: list = [1],
):
    """
    Plot symmetric Wheeler-Kiladis analysis with Matsuno dispersion curves.

    Parameters
    ----------
    max_freq : :py:class:`float <float>`, optional
        Maximum frequency to plot (CPD), default 0.5
    max_wn : :py:class:`float <float>`, optional
        Maximum wavenumber to plot, default 15
    ax : matplotlib.axes.Axes, optional
        Axes to plot on, creates new if None
    add_xylabel : :py:class:`bool <bool>`, optional
        Add x/y labels, default True
    add_central_line : :py:class:`bool <bool>`, optional
        Add central vertical line, default True
    add_westward_and_eastward : :py:class:`bool <bool>`, optional
        Add eastward/westward labels, default True
    auto_determine_xyrange : :py:class:`bool <bool>`, optional
        Auto-set axis ranges, default True
    freq_lines : :py:class:`bool <bool>`, optional
        Add frequency lines, default True
    matsuno_modes_labels : :py:class:`bool <bool>`, optional
        Add Matsuno mode labels, default True
    cpd_lines_levels : list, optional
        Periods (days) for frequency lines, default [3, 6, 30]
    matsuno_lines : :py:class:`bool <bool>`, optional
        Plot Matsuno dispersion curves, default True
    he : list, optional
        Equivalent depths for Matsuno curves, default [12, 25, 50]
    meridional_modes : list, optional
        Meridional mode numbers, default [1]

    Example
    -------
    >>> fig, ax = plt.subplots()
    >>> draw_wk_sym_analysis(ax=ax)

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_wk_spectra.py
    """
    matsuno_modes = matsuno_plot.matsuno_modes_wk(
        he=he, n=meridional_modes, max_wn=max_wn
    )

    if ax is None:
        ax = plt.gca()

    if add_xylabel:
        ax.set_xlabel("Zonal Wavenumber")
        ax.set_ylabel("Frequency (CPD)")

    if add_central_line:
        ax.axvline(x=0, color="k", linestyle="--")
    if add_westward_and_eastward:
        ax.text(max_wn - 2 * 0.25 * max_wn, -0.08, "Eastward")
        ax.text(-max_wn + 0.25 * max_wn, -0.08, "Westward")
    if auto_determine_xyrange:
        ax.set_xlim((-max_wn, max_wn))
        ax.set_ylim((0.02, max_freq))

    if freq_lines:
        for d in cpd_lines_levels:
            if (1.0 / d) <= max_freq:
                ax.axhline(y=1.0 / d, color="k", linestyle="--")
                ax.text(
                    -max_wn + 0.2,
                    (1.0 / d + 0.01),
                    str(d) + " days",
                    bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": "none"},
                )

    if matsuno_lines:
        for key in matsuno_modes:
            ax.plot(
                matsuno_modes[key]["Kelvin(he={}m)".format(key)],
                color="k",
                linestyle="--",
            )
            ax.plot(
                matsuno_modes[key]["ER(n=1,he={}m)".format(key)],
                color="k",
                linestyle="--",
            )
            ax.plot(
                matsuno_modes[key]["EIG(n=1,he={}m)".format(key)],
                color="k",
                linestyle="--",
            )
            ax.plot(
                matsuno_modes[key]["WIG(n=1,he={}m)".format(key)],
                color="k",
                linestyle="--",
            )

    if matsuno_modes_labels:
        key = list(matsuno_modes.keys())[len(list(matsuno_modes.keys())) // 2]
        wn = matsuno_modes[key].index.values

        # Print Kelvin Label
        i = int((len(wn) / 2) + 0.3 * (len(wn) / 2))
        (i,) = np.where(wn == wn[i])[0]
        ax.text(
            wn[i] - 1,
            float(matsuno_modes[key]["Kelvin(he={}m)".format(key)].iloc[i]),
            "Kelvin",
            bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": "none"},
        )

        # Print ER Label
        i = int(0.7 * (len(wn) / 2))
        i = np.where(wn == wn[i])[0]
        y = matsuno_modes[key]["ER(n=1,he={}m)".format(key)].iloc[i]
        y = float(y.iloc[0]) + 0.01
        ax.text(
            wn[i] - 1,
            y,
            "ER",
            bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": "none"},
        )

        key2 = list(matsuno_modes.keys())[0]
        wn2 = matsuno_modes[key].index.values

        # Print EIG Label
        i = int((len(wn2) / 2) + 0.3 * (len(wn2) / 2))
        (i,) = np.where(wn2 == wn2[i])[0]
        ax.text(
            wn2[i] - 1,
            float(matsuno_modes[key2]["EIG(n=1,he={}m)".format(key2)].iloc[i]),
            "EIG",
            bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": "none"},
        )

        # Print WIG Label
        i = int(0.55 * (len(wn2) / 2))
        (i,) = np.where(wn2 == wn2[i])[0]
        ax.text(
            wn2[i] - 1,
            float(matsuno_modes[key2]["WIG(n=1,he={}m)".format(key2)].iloc[i]),
            "WIG",
            bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": "none"},
        )
