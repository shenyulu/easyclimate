"""
Functions for Taylor Diagrams plotting.
"""

from __future__ import annotations

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import xarray as xr

__all__ = [
    "calc_TaylorDiagrams_metadata",
    "draw_TaylorDiagrams_base",
    "draw_TaylorDiagrams_metadata",
    "calc_TaylorDiagrams_values",
]


def calc_correlation_coefficient(f: xr.DataArray, r: xr.DataArray) -> xr.DataArray:
    """Calculate the correlation coefficient.

    .. math::
        R = \\frac{\\frac{1}{N} \\sum_{n=1}^N (f_n -\\bar{f})(r_n - \\bar{r})}{\\sigma_f \\sigma_r}

    Parameters
    ----------
    - f: :py:class:`xarray.DataArray<xarray.DataArray>`, required.
        A spatial array of models to be compared.
    - r: :py:class:`xarray.DataArray<xarray.DataArray>`, required.
        A spatial array for model reference comparisons (observations).

    .. attention::
        `f` and `r` must have the same two-dimensional spatial dimension.

    Returns
    -------
    - numpy array
        Pattern correlation coefficient array.

    """
    numerator = (f - f.mean()) * (r - r.mean())
    denominator = f.std() * r.std()
    return (numerator.mean() / denominator).data


def calc_standard_deviation(f: xr.DataArray, ddof=0) -> xr.DataArray:
    """Calculate the standard deviation.

    .. math::
        STD = \\frac{1}{N} \ \\sum_{n=1}^{N} \\left ( x_n - \\bar{x} \ \\right ) ^2

    Parameters
    ----------
    - f: :py:class:`xarray.DataArray<xarray.DataArray>`, required.
        A spatial array of models to be compared.
    - ddof: :py:class:`int <int>`, optional.
        "Delta Degrees of Freedom": the divisor used in the calculation is `N - ddof`, where `N` represents the number of elements.

    .. attention::
        `f` and `r` must have the same two-dimensional spatial dimension.

    Returns
    -------
    - numpy array
        Pattern standard deviation array.

    """
    return f.std(ddof=ddof).data


def calc_centeredRMS(
    f: xr.DataArray,
    r: xr.DataArray,
) -> xr.DataArray:
    """Calculate the center root mean square error.

    .. math::
        E'= \\left \\lbrace \\frac{1}{N} \\sum_{n=1}^N \\left[ (f_n - \\bar{f}) - (r_n - \\bar{r}) \\right] ^2  \\right \\rbrace ^{1/2}

    Parameters
    ----------
    f: :py:class:`xarray.DataArray<xarray.DataArray>`, required.
        A spatial array of models to be compared.
    r: :py:class:`xarray.DataArray<xarray.DataArray>`, required.
        A spatial array for model reference comparisons (observations).

    .. attention::
        `f` and `r` must have the same two-dimensional spatial dimension.

    Returns
    -------
    - centerRMS: :py:class:`numpy.array <numpy.array>`
        Pattern center root mean square error array.

    """
    tmp = ((f - f.mean()) - (r - r.mean())) ** 2
    tmp = tmp.mean()
    return np.sqrt(tmp).data


def calc_Taylor_skill_score(
    r: xr.DataArray, sigma_f: xr.DataArray, sigma_r: xr.DataArray, r0: float = 0.999
) -> xr.DataArray:
    """Calculate Taylor skill score (TSS).

    .. math::
        TSS = \\frac{4 * (1+r)^4}{(SDR + \\frac{1}{1+SDR})^2 (1+r_0)^4}

    Parameters
    ----------
    r: :py:class:`float <float>`, required.
        The correlation coefficient value.
    sigma_f: :py:class:`float <float>`, required.
        The standard deviation value of the model.
    sigma_r: :py:class:`float <float>`, required.
        The standard deviation value of the observation.
    r0: :py:class:`float <float>`, optional.
        Maximum correlation obtainable.
    """
    SDR = sigma_f / sigma_r
    numerator = 4 * (1 + r) ** 4
    denominator = (SDR + 1 / SDR) ** 2 * (1 + r0) ** 4
    return numerator / denominator


def calc_TaylorDiagrams_values(
    f: xr.DataArray,
    r: xr.DataArray,
    model_name: str,
    weighted: bool = False,
    lat_dim: str = "lat",
    normalized: bool = True,
    r0: float = 0.999,
) -> pd.DataFrame:
    """Calculate the center root mean square error.

    where :math:`N` is the number of points in spatial pattern.

    Parameters
    ----------
    f : :py:class:`xarray.DataArray<xarray.DataArray>`, required.
        A spatial array of models to be compared.
    r : :py:class:`xarray.DataArray<xarray.DataArray>`, required.
        A spatial array for model reference comparisons (observations).
    model_name: :py:class:`str <str>`, required.
        The name of the model.
    weighted: :py:class:`bool <bool>`, default `False`.
        Whether to weight the data by latitude or not? The default value is `False`.
    lat_dim: :py:class:`str <str>`, default `lat`.
        The name of `latitude` coordinate name.
    normalized: :py:class:`bool <bool>`, default `True`, optional.
        Whether the standard deviations is normalized, that is standard deviation of the model divided by that of the observations.
    r0 : :py:class:`float <float>`, optional.
        Maximum correlation obtainable.

    .. attention::
        `f` and `r` must have the same two-dimensional spatial dimension.

    Returns
    -------
    :py:class:`pandas.DataFrame <pandas.DataFrame>`.

    Reference
    --------------
    Taylor, K. E. (2001), Summarizing multiple aspects of model performance in a single diagram, J. Geophys. Res., 106(D7), 7183-7192, doi:`10.1029/2000JD900719 <https://doi.org/10.1029/2000JD900719>`__.

    """
    if weighted == True:
        # Spatial weighting
        weights_f = np.cos(np.deg2rad(f[lat_dim]))
        weights_r = np.cos(np.deg2rad(r[lat_dim]))
        f = f * weights_f
        r = r * weights_r
    elif weighted == False:
        pass

    correlation_coefficient = calc_correlation_coefficient(f, r)
    standard_deviation_f = calc_standard_deviation(f)
    standard_deviation_r = calc_standard_deviation(r)
    centered_RMS = calc_centeredRMS(f, r)
    TSS_value = calc_Taylor_skill_score(
        correlation_coefficient, standard_deviation_f, standard_deviation_r, r0=r0
    )

    correlation_coefficient = float(correlation_coefficient)
    standard_deviation_f = float(standard_deviation_f)
    standard_deviation_r = float(standard_deviation_r)
    centered_RMS = float(centered_RMS)
    TSS_value = float(TSS_value)

    if normalized == False:
        taylor_diagrams_data = {
            "item": [model_name],
            "std": [standard_deviation_f],
            "cc": [correlation_coefficient],
            "centeredRMS": [centered_RMS],
            "TSS": [TSS_value],
        }
    elif normalized == True:
        taylor_diagrams_data = {
            "item": [model_name],
            "std": [standard_deviation_f / standard_deviation_r],
            "cc": [correlation_coefficient],
            "centeredRMS": [centered_RMS / standard_deviation_r],
            "TSS": [TSS_value],
        }
    return pd.DataFrame(taylor_diagrams_data)


def calc_TaylorDiagrams_metadata(
    f: xr.DataArray,
    r: xr.DataArray,
    models_name: list[str] = [],
    weighted: bool = False,
    lat_dim: str = "lat",
    normalized: bool = True,
):
    """
    Calculating Taylor diagram metadata

    Parameters
    ----------
    f: :py:class:`xarray.DataArray<xarray.DataArray>`, required.
        A spatial array of models to be compared.
    r: :py:class:`xarray.DataArray<xarray.DataArray>`, required.
        A spatial array for model reference comparisons (observations).
    model_name: :py:class:`list[str]`, required.
        The `list` of the models' name.
    weighted: :py:class:`bool <bool>`, default `False`.
        Whether to weight the data by latitude or not? The default value is `False`.
    lat_dim: :py:class:`str <str>`, default `lat`.
        The name of `latitude` coordinate name.
    normalized: :py:class:`bool <bool>`, default `True`, optional.
        Whether the standard deviations is normalized, that is standard deviation of the model divided by that of the observations.

    Returns
    --------------
    :py:class:`pandas.DataFrame <pandas.DataFrame>`.

    Examples
    ---------------

    .. code:: python

        >>> import xarray as xr
        >>> import pandas as pd
        >>> import numpy as np
        >>> import easyclimate as ecl
        >>> da_a = xr.DataArray(
        ...:     np.array([[1, 2, 3], [0.1, 0.2, 0.3], [3.2, 0.6, 1.8]]),
        ...:     dims = ("lat", "time"),
        ...:     coords= {'lat': np.array([-30, 0, 30]),
        ...:              'time': pd.date_range("2000-01-01", freq="D", periods=3)
        ...:              }
        ...:)
        >>> da_a
        <xarray.DataArray (lat: 3, time: 3)>
        array([[1. , 2. , 3. ],
            [0.1, 0.2, 0.3],
            [3.2, 0.6, 1.8]])
        Coordinates:
        * lat      (lat) int32 -30 0 30
        * time     (time) datetime64[ns] 2000-01-01 2000-01-02 2000-01-03
        >>>  da_b = xr.DataArray(
        ...:     np.array([[0.2, 0.4, 0.6], [15, 10, 5], [3.2, 0.6, 1.8]]),
        ...:     dims = ("lat", "time"),
        ...:     coords= {'lat': np.array([-30, 0, 30]),
        ...:              'time': pd.date_range("2000-01-01", freq="D", periods=3)
        ...:              }
        ...:)
        >>>  da_b
        <xarray.DataArray (lat: 3, time: 3)>
        array([[ 0.2,  0.4,  0.6],
            [15. , 10. ,  5. ],
            [ 3.2,  0.6,  1.8]])
        Coordinates:
        * lat      (lat) int32 -30 0 30
        * time     (time) datetime64[ns] 2000-01-01 2000-01-02 2000-01-03
        >>>  da_obs = (da_a + da_b) / 1.85
        >>>  da_obs
        <xarray.DataArray (lat: 3, time: 3)>
        array([[0.64864865, 1.2972973 , 1.94594595],
            [8.16216216, 5.51351351, 2.86486486],
            [3.45945946, 0.64864865, 1.94594595]])
        Coordinates:
        * lat      (lat) int32 -30 0 30
        * time     (time) datetime64[ns] 2000-01-01 2000-01-02 2000-01-03
        >>>  ecl.calc_TaylorDiagrams_metadata(
        ...:     f = [da_a, da_b],
        ...:     r = [da_obs, da_obs],
        ...:     models_name = ['f1', 'f2'],
        ...:     weighted = True,
        ...:     normalized = True,
        ...:)
        item       std                   cc  centeredRMS       TSS
        0  Obs  1.000000                  1.0     0.000000  1.002003
        1   f1  0.404621  -0.4293981636461462     1.229311  0.003210
        2   f2  2.056470    0.984086060161888     1.087006  0.600409

    """
    # Check that the list lengths of `f`, `r` and `models_name` are aligned.
    length_models = len(f)
    length_referance = len(r)
    length_models_name = len(models_name)
    if length_models == length_referance == length_models_name:
        pass
    else:
        raise ValueError(
            f"Parameters `f`, `r`, and `models_name` should have same length! But the length of `f`, `r`, and `models_name` are {length_models}, {length_referance}, {length_models_name}, respectively."
        )

    dataframe_all = []
    # Obs
    dataframe_part = calc_TaylorDiagrams_values(
        r[0],
        r[0],
        model_name="Obs",
        weighted=weighted,
        lat_dim=lat_dim,
        normalized=normalized,
    )
    dataframe_all.append(dataframe_part)
    # Models
    for intern in np.arange(length_models):
        f_intern = f[intern]
        r_intern = r[intern]
        models_name_intern = models_name[intern]
        dataframe_part = calc_TaylorDiagrams_values(
            f_intern,
            r_intern,
            model_name=models_name_intern,
            weighted=weighted,
            lat_dim=lat_dim,
            normalized=normalized,
        )
        dataframe_all.append(dataframe_part)

    return pd.concat(objs=dataframe_all, ignore_index=True)


def draw_TaylorDiagrams_base(
    ax: matplotlib.axes.Axes = None,
    half_circle: bool = False,
    normalized: bool = True,
    std_min: float = 0,
    std_max: float = 2,
    std_interval: float = 0.25,
    arc_label: str = "Correlation",
    arc_label_pad: float = 0.2,
    arc_label_kwargs: dict = {"fontsize": 12},
    arc_ticker_kwargs: dict = {"lw": 0.8, "c": "black"},
    arc_tickerlabel_kwargs: dict = {"labelsize": 12, "pad": 8},
    arc_ticker_length: float = 0.02,
    arc_minorticker_length: float = 0.01,
    x_label: str = "Std (Normalized)",
    x_label_pad: float = 0.25,
    x_label_kwargs: dict = {"fontsize": 12},
    x_ticker_length: float = 0.02,
    x_tickerlabel_kwargs: dict = {"fontsize": 12},
    x_ticker_kwargs: dict = {"lw": 0.8, "c": "black"},
    y_ticker_kwargs: dict = {"lw": 0.8, "c": "black"},
) -> matplotlib.collections.Collection:
    """
    Drawing Taylor Graphics Basic Framework

    Parameters
    ----------
    ax: :py:class:`matplotlib.axes.Axes <matplotlib.axes.Axes>`, optional.
        Axes on which to plot. By default, use the current axes, i.e. `ax = plt.gca()`.
    half_circle: :py:class:`bool <bool>`, default `False`, optional.
        Whether to draw the `'half-circle'` version of the Taylor diagram.
    normalized: :py:class:`bool <bool>`, default `True`, optional.
        Whether the standard deviations is normalized, that is standard deviation of the model divided by that of the observations.
        This parameter mainly affects the label `x=1` on the `x` axis, if normalized to True, it is rewritten to `REF`.
    std_min: :py:class:`float <float>`, default `0.0`, optional.
        Minimum value of x-axis (standard deviation) on Taylor diagram.

        .. note:: The value of `std_min` shoud be 0 in the `'half-circle'` version of the Taylor diagram.

    std_max: :py:class:`float <float>`, default `2.0`, optional.
        Maximum value of x-axis (standard deviation) on Taylor diagram.
    std_interval: :py:class:`float <float>`, default `0.25`, optional.
        The interval between the ticker on the x-axis (standard deviation) between the minimum and maximum values on the Taylor diagram.
    arc_label: :py:class:`str <str>`, default `'Correlation'`, optional.
        Label on Taylor chart arc, default value is `'Correlation'`.
    arc_label_pad: :py:class:`float <float>`, default `0.2`, optional.
        The offset of the title from the top of the arc, based on x-axis based coordinate system.
    arc_label_kwargs: :py:class:`dict <dict>`, default `{'fontsize': 12}`, optional.
        Additional keyword arguments passed on to labels on arcs, according to other miscellaneous parameters in`matplotlib.axes.Axes.text`.
    arc_ticker_kwargs: :py:class:`dict <dict>`, default `{'lw': 0.8, 'c': 'black'}`, optional.
        Additional keyword arguments passed on to tickers on arcs, according to other miscellaneous parameters in`matplotlib.axes.Axes.plot`.
    arc_tickerlabel_kwargs: :py:class:`dict <dict>`, default `{'labelsize': 12, 'pad': 8}`, optional.
        Additional keyword arguments passed on to tickers on arcs, according to other miscellaneous parameters in`matplotlib.axes.Axes.tick_params`.
    arc_ticker_length: :py:class:`float <float>`, default `0.02`, optional.
        Ticker length on arc.
    arc_minorticker_length: :py:class:`float <float>`, default `0.01`, optional.
        Minor ticker length on arc.
    x_label: :py:class:`str <str>`, default `'Std (Normalized)'`, optional.
        Label on Taylor chart x axis, default value is `'Std (Normalized)'`.
    x_label_pad: :py:class:`float <float>`, default `0.25`, optional.
        The offset of the title from the top of the x-axis, based on x-axis based coordinate system.
    x_label_kwargs: :py:class:`dict <dict>`, default `{'fontsize': 12}`, optional.
        Additional keyword arguments passed on to labels on x-axis, according to other miscellaneous parameters in`matplotlib.axes.Axes.text`.
    x_ticker_length: :py:class:`float <float>`, default `0.02`, optional.
        Ticker length on x-axis
    x_tickerlabel_kwargs: :py:class:`dict <dict>`, default `{'fontsize': 12}`, optional.
        Additional keyword arguments passed on to tickers' labels on x-axis, according to other miscellaneous parameters in`matplotlib.axes.Axes.text`.
    x_ticker_kwargs: :py:class:`dict <dict>`, default `{'lw': 0.8, 'c': 'black'}`, optional.
        Additional keyword arguments passed on to tickers on x-axis, according to other miscellaneous parameters in`matplotlib.axes.Axes.plot`.
    y_ticker_kwargs: :py:class:`dict <dict>`, default `{'lw': 0.8, 'c': 'black'}`, optional.
        Additional keyword arguments passed on to tickers on y-axis, according to other miscellaneous parameters in`matplotlib.axes.Axes.plot`.

    Returns
    -------
    :py:class:`matplotlib.collections.Collection <matplotlib.collections.Collection>`.
    """

    # Get Axes
    if ax == None:
        ax = plt.gca()
    else:
        pass

    # Check if `projection` is `polar`
    if ax.name == "polar":
        pass
    else:
        raise TypeError(
            "The projection type of the Axes should be `polar`, consider to specify the parameter `projection` as `polar`. E.g. `fig, ax = plt.subplots(subplot_kw = {'projection': 'polar'})`."
        )

    if half_circle == True:
        if std_min != 0:
            raise ValueError(
                "The value of `std_min` shoud be 0 in the `'half-circle'` version of the Taylor diagram."
            )
    elif half_circle == False:
        if std_min < 0:
            raise ValueError("`std_min` should be greater than 0!")
    else:
        raise ValueError("`half_circle` should be Boolean type.")

    if std_max < 0:
        raise ValueError("`std_max` should be greater than 0!")

    # The 'half-circle' version of the Taylor Diagram
    if half_circle == True:
        # Major ticker of R value
        rad_list = np.array(
            [
                -1.0,
                -0.99,
                -0.95,
                -0.9,
                -0.85,
                -0.8,
                -0.7,
                -0.6,
                -0.4,
                -0.2,
                0,
                0.2,
                0.4,
                0.6,
                0.7,
                0.8,
                0.85,
                0.9,
                0.95,
                0.99,
                1,
            ]
        )
        # Minor ticker of R value
        minor_rad_list = np.array(
            [
                -1.0,
                -0.99,
                -0.98,
                -0.97,
                -0.96,
                -0.95,
                -0.94,
                -0.93,
                -0.92,
                -0.91,
                -0.9,
                -0.89,
                -0.88,
                -0.87,
                -0.86,
                -0.85,
                -0.8,
                -0.75,
                -0.7,
                -0.65,
                -0.6,
                -0.5,
                -0.4,
                -0.3,
                -0.2,
                -0.1,
                0,
                0.1,
                0.2,
                0.3,
                0.4,
                0.5,
                0.6,
                0.65,
                0.7,
                0.75,
                0.8,
                0.85,
                0.86,
                0.87,
                0.88,
                0.89,
                0.9,
                0.91,
                0.92,
                0.93,
                0.94,
                0.95,
                0.96,
                0.97,
                0.98,
                0.99,
                1,
            ]
        )
        # polar plot theta range
        ax.set_thetalim(thetamin=0, thetamax=180)
        # x_label
        ax.text(
            np.deg2rad(270),
            x_label_pad,
            s=x_label,
            ha="center",
            va="top",
            **x_label_kwargs,
        )
        ax.text(
            np.deg2rad(45),
            std_max + arc_label_pad,
            s=arc_label,
            ha="center",
            va="bottom",
            rotation=-45,
            **arc_label_kwargs,
        )

        # Draw x, y label
        if normalized == True:
            for i in np.arange(std_min, std_max, std_interval):
                if i == 1:
                    # The first coordinate of `text` is the angle (radian system) and the second is the distance
                    ax.text(
                        0,
                        i,
                        s="\n" + "REF",
                        ha="center",
                        va="top",
                        **x_tickerlabel_kwargs,
                    )
                else:
                    # The first coordinate of `text` is the angle (radian system) and the second is the distance
                    ax.text(
                        0,
                        i,
                        s="\n" + str(i),
                        ha="center",
                        va="top",
                        **x_tickerlabel_kwargs,
                    )
            for i in np.trim_zeros(np.arange(std_min, std_max, std_interval)) * (-1):
                if i == -1:
                    # The first coordinate of `text` is the angle (radian system) and the second is the distance
                    ax.text(
                        np.deg2rad(180),
                        -i,
                        s="\n" + "REF",
                        ha="center",
                        va="top",
                        **x_tickerlabel_kwargs,
                    )
                else:
                    # The first coordinate of `text` is the angle (radian system) and the second is the distance
                    ax.text(
                        np.deg2rad(180),
                        -i,
                        s="\n" + str(i),
                        ha="center",
                        va="top",
                        **x_tickerlabel_kwargs,
                    )
        elif normalized == False:
            for i in np.arange(std_min, std_max, std_interval):
                # The first coordinate of `text` is the angle (radian system) and the second is the distance
                ax.text(
                    0, i, s="\n" + str(i), ha="center", va="top", **x_tickerlabel_kwargs
                )

            for i in np.trim_zeros(np.arange(std_min, std_max, std_interval)) * (-1):
                # The first coordinate of `text` is the angle (radian system) and the second is the distance
                ax.text(
                    np.deg2rad(180),
                    -i,
                    s="\n" + str(i),
                    ha="center",
                    va="top",
                    **x_tickerlabel_kwargs,
                )

        # draw line `x=0` on polar plot
        # https://stackoverflow.com/questions/53730598/how-to-draw-a-curved-line-arc-in-a-polar-plot-with-matplotlib
        rad = np.deg2rad(90)
        ax.plot([rad, rad], [0, std_max], color="grey", linewidth=0.8)

        for value in np.arange(std_interval, std_max * 2, std_interval):
            rect = patches.Circle(
                (1, 0),
                value,
                transform=ax.transData._b,
                facecolor=(0, 0, 0, 0),
                edgecolor="grey",
                linestyle="-",
                linewidth=0.8,
            )
            ax.add_patch(rect)
    elif half_circle == False:
        # Major ticker of R value
        rad_list = [0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.99, 1]
        # Minor ticker of R value
        minor_rad_list = [
            0,
            0.1,
            0.2,
            0.3,
            0.4,
            0.5,
            0.6,
            0.65,
            0.7,
            0.75,
            0.8,
            0.85,
            0.86,
            0.87,
            0.88,
            0.89,
            0.9,
            0.91,
            0.92,
            0.93,
            0.94,
            0.95,
            0.96,
            0.97,
            0.98,
            0.99,
            1,
        ]

        ax.set_thetalim(thetamin=0, thetamax=90)
        ax.text(
            2 * np.pi - np.arctan(x_label_pad / (std_max / 2)),
            np.sqrt((std_max / 2) ** 2 + x_label_pad**2),
            s=x_label,
            ha="center",
            va="top",
            **x_label_kwargs,
        )
        ax.text(
            np.deg2rad(45),
            std_max + arc_label_pad,
            s=arc_label,
            ha="center",
            va="bottom",
            rotation=-45,
            **arc_label_kwargs,
        )

        if normalized == True:
            for i in np.arange(std_min, std_max, std_interval):
                if i == 1:
                    # The first coordinate of `text` is the angle (radian system) and the second is the distance
                    ax.text(
                        0,
                        i,
                        s="\n" + "REF",
                        ha="center",
                        va="top",
                        **x_tickerlabel_kwargs,
                    )
                else:
                    # The first coordinate of `text` is the angle (radian system) and the second is the distance
                    ax.text(
                        0,
                        i,
                        s="\n" + str(i),
                        ha="center",
                        va="top",
                        **x_tickerlabel_kwargs,
                    )

                # The first coordinate of `text` is the angle (radian system) and the second is the distance
                ax.text(
                    np.pi / 2,
                    i,
                    s=str(i) + "  ",
                    ha="right",
                    va="center",
                    **x_tickerlabel_kwargs,
                )

            # Make a circle with `ref` as the center and a multiple of `std_interval` as the radius
            for value in np.arange(std_interval, std_max, std_interval):
                rect = patches.Circle(
                    (1, 0),
                    value,
                    transform=ax.transData._b,
                    facecolor=(0, 0, 0, 0),
                    edgecolor="grey",
                    linestyle="-",
                    linewidth=0.8,
                )
                ax.add_patch(rect)
        elif normalized == False:
            for i in np.arange(std_min, std_max, std_interval):
                # The first coordinate of `text` is the angle (radian system) and the second is the distance
                ax.text(
                    0, i, s="\n" + str(i), ha="center", va="top", **x_tickerlabel_kwargs
                )

                # The first coordinate of `text` is the angle (radian system) and the second is the distance
                ax.text(
                    np.pi / 2,
                    i,
                    s=str(i) + "  ",
                    ha="right",
                    va="center",
                    **x_tickerlabel_kwargs,
                )

    # Remove default setting of polar plot
    ax.set_rgrids([])
    ax.grid(False)

    # arc label
    ax.set_rlim(std_min, std_max)
    angle_list = np.rad2deg(np.arccos(rad_list))
    angle_list_rad = np.arccos(rad_list)
    angle_minor_list = np.arccos(minor_rad_list)
    ax.set_thetagrids(angle_list, rad_list)
    ax.tick_params(axis="x", **arc_tickerlabel_kwargs)

    # arc ticker
    tick = [ax.get_rmax(), ax.get_rmax() * (1 - arc_ticker_length)]
    tick_minor = [ax.get_rmax(), ax.get_rmax() * (1 - arc_minorticker_length)]
    for t in angle_list_rad:
        # The first coordinate is the angle (angle system), the second is the distance
        ax.plot([t, t], tick, **arc_ticker_kwargs)
    for t in angle_minor_list:
        # The first coordinate is the angle (angle system), the second is the distance
        ax.plot([t, t], tick_minor, **arc_ticker_kwargs)

    # x ticker
    if half_circle == True:
        # Positive x-axis
        x_list = np.arange(std_min, std_max, std_interval)
        for t in x_list:
            if t == 0:
                rad = np.deg2rad(90)
                ax.plot([rad, rad], [0, x_ticker_length], **x_ticker_kwargs)
            else:
                ax.plot(
                    [0, np.arctan(x_ticker_length / t)],
                    [t, np.sqrt(t**2 + x_ticker_length**2)],
                    **x_ticker_kwargs,
                )
        # Negative x-axis
        x_list = np.arange(std_min, std_max, std_interval)
        for t in x_list:
            if t == 0:
                rad = np.deg2rad(90)
                ax.plot([rad, rad], [0, x_ticker_length], **x_ticker_kwargs)
            else:
                ax.plot(
                    [np.pi, np.pi - np.arctan(x_ticker_length / t)],
                    [t, np.sqrt(t**2 + x_ticker_length**2)],
                    **x_ticker_kwargs,
                )
    elif half_circle == False:
        # x-axis
        x_list = np.arange(std_min, std_max, std_interval)
        for t in x_list:
            if t == 0:
                rad = np.deg2rad(90)
                ax.plot([rad, rad], [0, x_ticker_length], **x_ticker_kwargs)
            else:
                ax.plot(
                    [0, np.arctan(x_ticker_length / t)],
                    [t, np.sqrt(t**2 + x_ticker_length**2)],
                    **x_ticker_kwargs,
                )
        # y-axis
        x_list = np.arange(std_min, std_max, std_interval)
        for t in x_list:
            if t == 0:
                rad = np.deg2rad(90)

                ax.plot([rad, rad], [0, x_ticker_length], **y_ticker_kwargs)
            else:
                ax.plot(
                    [np.pi / 2, np.arctan(t / x_ticker_length)],
                    [t, np.sqrt(t**2 + x_ticker_length**2)],
                    **y_ticker_kwargs,
                )

    # Make a circle with 0 as the center and 2 times `std_interval` as the radius
    for value in np.arange(std_interval * 2, std_max, std_interval * 2):
        rect = patches.Circle(
            (0, 0),
            value,
            transform=ax.transData._b,
            facecolor=(0, 0, 0, 0),
            edgecolor="grey",
            linestyle=(5, (10, 3)),
            linewidth=0.8,
        )
        ax.add_patch(rect)

    # Draw Rays
    for rad_num in rad_list:
        ax.plot(
            [0, np.arccos(rad_num)], [0, std_max], lw=0.7, color="grey", linestyle="--"
        )


def draw_TaylorDiagrams_metadata(
    taylordiagrams_metadata: pd.DataFrame,
    marker_list: list,
    color_list: list,
    label_list: list,
    legend_list: list,
    ax: matplotlib.axes.Axes = None,
    normalized: bool = True,
    cc: str = "cc",
    std: str = "std",
    point_label_xoffset: float = 0,
    point_label_yoffset: float = 0.05,
    point_kwargs: dict = {"alpha": 1, "markersize": 6.5},
    point_label_kwargs: dict = {"fontsize": 14},
) -> matplotlib.collections.Collection:
    """
    Draw points to Taylor Graphics Basic Framework according to Taylor diagram metadata.

    Parameters
    ----------
    taylordiagrams_metadata: :py:class:`pandas.DataFrame <pandas.DataFrame>`, required.
        Taylor diagram metadata generated by the function `calc_TaylorDiagrams_metadata`.
    marker_list: :py:class:`list <list>`, required.
        The list of markers. The order of `marker` in `marker_list` is determined by the order in `taylordiagrams_metadata`.
        See `matplotlib.markers` for full description of possible arguments.
    color_list: :py:class:`list <list>`, required.
        The list of colors. The order of `color` in `color_list` is determined by the order in `taylordiagrams_metadata`.
    label_list: :py:class:`list <list>`, required.
        The list of data point labels (marked next to plotted points).
        The order of label in `label_list` is determined by the order in `taylordiagrams_metadata`.
    legend_list: :py:class:`list <list>`, required.
        The list of legend label.
        The order of label in `legend_list` is determined by the order in `taylordiagrams_metadata`.
    ax: :py:class:`matplotlib.axes.Axes <matplotlib.axes.Axes>`, optional.
        Axes on which to plot. By default, use the current axes, i.e. `ax = plt.gca()`.
    normalized: :py:class:`bool <bool>`, default `True`, optional.
        Whether the standard deviations is normalized, that is standard deviation of the model divided by that of the observations.
    cc: :py:class:`str <str>`, default `'cc'`, optional.
        The name of correlation coefficient in `taylordiagrams_metadata`.
    std: :py:class:`str <str>`, default `'std'`, optional.
        The name of standard deviation in `taylordiagrams_metadata`.
    point_label_xoffset: :py:class:`float <float>`, optional.
        The offset of the labels from the points, based on x-axis based coordinate system.
    point_label_yoffset: :py:class:`float <float>`, optional.
        The offset of the labels from the points, based on y-axis based coordinate system.
    point_kwargs: :py:class:`dict <dict>`, optional.
        Additional keyword arguments passed on to data points, according to other miscellaneous parameters in`matplotlib.axes.Axes.plot`.
    point_label_kwargs: :py:class:`dict <dict>`, optional.
        Additional keyword arguments passed on to the labels of data points, according to other miscellaneous parameters in`matplotlib.axes.Axes.text`.

    Returns
    -------
    :py:class:`matplotlib.collections.Collection <matplotlib.collections.Collection>`.
    """
    # Get Axes
    if ax == None:
        ax = plt.gca()
    else:
        pass

    # Check if `projection` is `polar`
    if ax.name == "polar":
        pass
    else:
        raise TypeError(
            "The projection type of the Axes should be `polar`, consider to specify the parameter `projection` as `polar`. E.g. `fig, ax = plt.subplots(subplot_kw = {'projection': 'polar'})`."
        )

    # Number of rows
    datalist_row = taylordiagrams_metadata.shape[0]

    plot_list = []

    for intern in np.arange(datalist_row):
        taylor_diagrams_values_single_list = taylordiagrams_metadata[
            intern : intern + 1
        ]

        if isinstance(point_label_xoffset, list):
            point_label_xoffset_item = point_label_xoffset[intern]
        elif isinstance(point_label_xoffset, (float, int)):
            point_label_xoffset_item = point_label_xoffset
        else:
            raise ValueError(
                "`point_label_xoffset` shuold be `list`, `int` or `float`."
            )

        if isinstance(point_label_yoffset, list):
            point_label_yoffset_item = point_label_yoffset[intern]
        elif isinstance(point_label_yoffset, (float, int)):
            point_label_yoffset_item = point_label_yoffset
        else:
            raise ValueError(
                "`point_label_yoffset` shuold be `list`, `int` or `float`."
            )

        if (
            taylor_diagrams_values_single_list["item"].to_list()[0] == "Obs"
            and normalized == True
        ):
            x = 0
            y = 1
            (point_tmp,) = ax.plot(
                x,
                y,
                marker_list[intern],
                color=color_list[intern],
                label=legend_list[intern],
                **point_kwargs,
            )
            ax.text(
                x + point_label_xoffset_item,
                y + point_label_yoffset_item,
                s=label_list[intern],
                **point_label_kwargs,
            )
            plot_list.append(point_tmp)
        else:
            x = taylor_diagrams_values_single_list[cc][intern]
            y = taylor_diagrams_values_single_list[std][intern]
            (point_tmp,) = ax.plot(
                float(np.arccos(x)),
                y,
                marker_list[intern],
                color=color_list[intern],
                label=legend_list[intern],
                **point_kwargs,
            )
            ax.text(
                np.arccos(x) + point_label_xoffset_item,
                y + point_label_yoffset_item,
                s=label_list[intern],
                **point_label_kwargs,
            )
            plot_list.append(point_tmp)

    return plot_list
