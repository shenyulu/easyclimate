"""
Wheeler-Kiladis Space-Time Spectra
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from easyclimate_backend.wk_spectra import wk_analysis, matsuno_plot

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
    data = data.transpose(time_dim, lat_dim, lon_dim)
    result_array = wk_analysis.remove_dominant_signals(
        data.data, spd, nDayWin, nDaySkip
    )
    return data.copy(deep=True, data=result_array)


def decompose_symasym(da, lat_dim="lat"):
    """
    Decompose an xarray DataArray into symmetric and asymmetric parts about the equator
    along the 'lat' dimension. Supports arrays with 'lat' and optional dimensions like
    'lon', 'time', 'level', etc.

    Parameters
    ----------
    da : xarray.DataArray
        Input array with 'lat' dimension

    Returns
    -------
    xarray.DataArray
        Decomposed array with symmetric part in Southern Hemisphere and asymmetric
        part in Northern Hemisphere
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
    cpd_lines_levels=[3, 6, 30],
    matsuno_lines: bool = True,
    he=[12, 25, 50],
    meridional_modes=[1],
):
    """ """
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
    cpd_lines_levels=[3, 6, 30],
    matsuno_lines: bool = True,
    he=[12, 25, 50],
    meridional_modes=[1],
):
    """ """
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
