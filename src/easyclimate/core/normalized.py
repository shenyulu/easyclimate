"""
Normalized Data
"""

import xarray as xr
from typing import Literal
from xarray.groupers import SeasonResampler

__all__ = [
    "timeseries_normalize_zscore",
    "timeseries_normalize_minmax",
    "timeseries_normalize_robust",
    "timeseries_normalize_mean",
    "calc_precip_anomaly_percentage",
]


def timeseries_normalize_zscore(
    da: xr.DataArray,
    dim: str = "time",
    time_range: slice = slice(None, None),
    ddof: int = 1,
) -> xr.DataArray:
    """
    Perform Z-Score standardization on an xarray time series.

    This function standardizes the input data by transforming it to have a mean of 0 and a standard deviation of 1,
    using the formula:

    .. math::

        z = \\frac{x - \\mu}{\\sigma}

    where :math:`\\mu` is the mean and :math:`\\sigma` is the standard deviation.

    Parameters
    ----------
    da: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The input time series data to be standardized.
    dim: :py:class:`str <str>`, default: `time`.
        The dimension along which to compute the mean and standard deviation. By default, standardization is applied over the `time` dimension.
    time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
        The time range of ``da`` to be normalized. The default value is the entire time range.
    ddof: :py:class:`int <int>`, default: `1`.
        Delta degrees of freedom for standard deviation calculation. The divisor used in calculations is :math:`N - \\mathrm{ddof}`, where :math:`N` is the number of elements.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
        The standardized data with mean 0 and standard deviation 1 along the specified dimension.

    .. note ::

        - **Applicable Scenarios**: Suitable for data that is approximately normally distributed or when comparing variables with different units in machine learning or statistical analysis.
        - **Advantages**: Retains the relative distribution characteristics of the data, widely used in algorithms requiring standardized inputs.
        - **Disadvantages**: Sensitive to outliers, which can skew the mean and standard deviation.

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_corr_reg.py
    """
    da_ = da.sel({dim: time_range})
    mean = da_.mean(dim=dim)
    std = da_.std(dim=dim, ddof=ddof)
    return (da - mean) / std


def timeseries_normalize_minmax(
    da: xr.DataArray,
    dim: str = "time",
    time_range: slice = slice(None, None),
    feature_range: tuple[float, float] = (0, 1),
) -> xr.DataArray:
    """
    Perform Min-Max standardization on an xarray time series.

    This function linearly scales the data to a specified range (default :math:`[0, 1]`) using the formula:

    .. math::

        x' = \\frac{x - x_{\\min}}{x_{\\max} - x_{\\min}} \\cdot (b - a) + a

    where :math:`(a, b)` are the target range bounds.

    Parameters
    ----------
    da: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The input time series data to be standardized.
    dim: :py:class:`str <str>`, default: `time`.
        The dimension along which to compute the minimum and maximum values. By default, standardization is applied over the `time` dimension.
    time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
        The time range of ``da`` to be normalized. The default value is the entire time range.
    feature_range: :py:class:`tuple[float, float] <tuple>`, default: ``(0, 1)``.
        The target range for scaling the data, specified as (min, max).

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
        The standardized data scaled to the specified range.

    .. note ::

        - **Applicable Scenarios**: Ideal for neural network inputs or when data needs to be constrained to a fixed range.
        - **Advantages**: Simple and intuitive, preserves relative relationships in the data.
        - **Disadvantages**: Sensitive to outliers, as the range depends on the minimum and maximum values.
    """
    da_ = da.sel({dim: time_range})
    min_val = da_.min(dim=dim)
    max_val = da_.max(dim=dim)
    a, b = feature_range
    return (da - min_val) / (max_val - min_val) * (b - a) + a


def timeseries_normalize_robust(
    da: xr.DataArray,
    dim: str = "time",
    time_range: slice = slice(None, None),
    q_low: float = 0.25,
    q_high: float = 0.75,
) -> xr.DataArray:
    """
    Perform Robust standardization on an xarray time series.

    This function standardizes the data using the median and interquartile range (IQR), with the formula:

    .. math::

        x' = \\frac{x - \\text{median}}{\\text{IQR}},

    where :math:`\\mathrm{IQR = Q_3 - Q_1}`.

    Parameters
    ----------
    da: :py:class:`xarray.DataArray <xarray.DataArray>`.
        The input time series data to be standardized.
    dim: :py:class:`str <str>`, default: `time`.
        The dimension along which to compute the median and IQR. By default, standardization is applied over the `time` dimension.
    time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
        The time range of ``da`` to be normalized. The default value is the entire time range.
    q_low: :py:class:`float <float>`, default: `0.25`.
        The lower quantile for IQR calculation (:math:`Q_1`).
    q_high: :py:class:`float <float>`, default: `0.75`.
        The upper quantile for IQR calculation (:math:`Q_3`).

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
        The standardized data based on median and IQR.

    .. note ::

        - **Applicable Scenarios**: Suitable for data with many outliers or non-normal distributions.
        - **Advantages**: Robust to outliers, providing a more stable standardization for skewed data.
        - **Disadvantages**: May lose some distribution information compared to Z-Score standardization.
    """
    da_ = da.sel({dim: time_range})
    median = da_.median(dim=dim)
    q75 = da_.quantile(q_high, dim=dim)
    q25 = da_.quantile(q_low, dim=dim)
    iqr = q75 - q25
    return (da - median) / iqr


def timeseries_normalize_mean(
    da: xr.DataArray,
    dim: str = "time",
    time_range: slice = slice(None, None),
) -> xr.DataArray:
    """
    Perform Mean normalization on an xarray time series.

    This function centers the data around zero and scales it by the range, using the formula:

    .. math::

        x' = \\frac{x - \\mu}{x_{\\max} - x_{\\min}},

    where :math:`\\mu` is the mean.

    Parameters
    ----------
    da: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The input time series data to be standardized.
    dim: :py:class:`str <str>`, default: `time`.
        The dimension along which to compute the mean and range. By default, standardization is applied over the `time` dimension.
    time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
        The time range of ``da`` to be normalized. The default value is the entire time range.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
        The normalized data centered around zero.

    .. note ::

        - **Applicable Scenarios**: Useful when centering data is needed without enforcing a standard deviation of 1.
        - **Advantages**: Simple, partially preserves data distribution characteristics.
        - **Disadvantages**: Sensitive to outliers, and the scaling range is not fixed.
    """
    da_ = da.sel({dim: time_range})
    mean = da_.mean(dim=dim)
    min_val = da_.min(dim=dim)
    max_val = da_.max(dim=dim)
    return (da - mean) / (max_val - min_val)


def calc_precip_anomaly_percentage(
    precip_data: xr.DataArray,
    freq: Literal["monthly", "seasonly", "yearly"] = "monthly",
    time_range: slice = slice(None, None),
    time_dim: str = "time",
) -> xr.DataArray:
    """
    Calculate Precipitation Anomaly Percentage (PAP, 降水距平百分率)

    .. math::

        P _{a} = \\frac {P - \\bar {P}} {\\bar {P}} \\times 100 \\%

    Where, :math:`P_a` is PAP, :math:`P` is the rainfall of a certain period, :math:`\\bar{P}= \\frac {1} {n}\\sum _{i=1}^{n} P_{ i }` is the long-term average rainfall of the period, :math:`n` is :math:`1` to :math:`n` years, :math:`i=1,2,\\cdots,n`.

    Parameters
    ----------
    precip_data : :py:class:`xarray.DataArray <xarray.DataArray>`.
        Precipitation data, recommended units: ``mm/month`` or ``mm/day`` (converted to monthly cumulative).
        Dimensions must include time, e.g., (time, lat, lon)

    .. caution::

        ``precip_data`` should be applied to **monthly means** precipitation.

    time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
        The time range of baseline climatology period, e.g., ``slice('1991-01', '2020-12')``. The default value is the entire time range.

    time_dim : str, default ``"time"``
        Time dimension name

    freq : {"monthly", "seasonly", "yearly", or custom seasons}.
        Time grouping method, options:

        - "monthly": Calculate climatology by month
        - "seasonly": Calculate by meteorological seasons (DJF, MAM, JJA, SON)
        - "yearly": Calculate climatology by year
        - custom seasons: Calculate climatology by custom seasons, e.g., ``JJAS``, ``OND``.

        .. note::

            This method is based on :py:class:`xarray.groupers.SeasonResampler <xarray.groupers.SeasonResampler>` for resampling, and you can customize the seasons by referring to the `examples <https://docs.xarray.dev/en/stable/generated/xarray.groupers.SeasonResampler.html#xarray.groupers.SeasonResampler>`__.

    Returns
    ----------
    pap : :py:class:`xarray.DataArray <xarray.DataArray>` (%)
        Precipitation Anomaly Percentage (PAP).

    Reference
    --------------
    - Zhai, P., Zhang, X., Wan, H., & Pan, X. (2005). Trends in total precipitation and frequency of daily precipitation extremes over China. Journal of Climate, 18(7), 1096–1108. https://doi.org/10.1175/JCLI-3318.1
    - Zou, Y., Wu, H., Lin, X., & Wang, Y. (2019). A quantitative method for the assessment of annual state of climate (气候年景定量化评价方法). Acta Meteorologica Sinica (in Chinese), 77(6), 1124–1133. https://doi.org/10.11676/qxxb2019.067
    - GB/T 20481-2017, Classification of meteorological drought (气象干旱等级, in Chinese) https://std.samr.gov.cn/gb/search/gbDetailed?id=71F772D81C2DD3A7E05397BE0A0AB82A.
    - Yang Shao-E and Wu Bing-fang, "Calculation of monthly precipitation anomaly percentage using web-serviced remote sensing data," 2010 2nd International Conference on Advanced Computer Control, Shenyang, China, 2010, pp. 621-625, doi: http://doi.org/10.1109/ICACC.2010.5486796.
    - Ma, Y., Zhao, L., Wang, J.-S., & Yu, T. (2021). Increasing difference of China summer precipitation statistics between percentage anomaly and probability distribution methods due to tropical warming. Earth and Space Science, 8, e2021EA001777. https://doi.org/10.1029/2021EA001777
    - Wang, Y., Wang, S., Luo, F., & Wang, H. (2022). Strengthened impacts of Indian Ocean Dipole on the Yangtze precipitation contribute to the extreme rainfall of 2020 Meiyu season. Journal of Geophysical Research: Atmospheres, 127, e2022JD037028. https://doi.org/10.1029/2022JD037028
    """

    # Select baseline climatology period
    if isinstance(time_range, (tuple, list)):
        time_range = slice(time_range[0], time_range[1])

    # Group by month or season
    if freq == "seasonly":
        precip_seasonal_data = precip_data.resample({time_dim: "QS-DEC"}).mean()
        pr_clim = (
            precip_seasonal_data.sel({time_dim: time_range})
            .groupby("time.month")
            .mean()
        )
        tmp = precip_seasonal_data.groupby("time.month") - pr_clim
        pap = tmp.groupby("time.month") / pr_clim * 100
    elif freq == "monthly":
        pr_clim = precip_data.sel({time_dim: time_range}).groupby("time.month").mean()
        tmp = precip_data.groupby("time.month") - pr_clim
        pap = tmp.groupby("time.month") / pr_clim * 100
    elif freq == "yearly":
        precip_yearly_data = precip_data.resample({time_dim: "YS"}).mean()
        pr_clim = (
            precip_yearly_data.sel({time_dim: time_range}).groupby("time.month").mean()
        )
        tmp = precip_yearly_data.groupby("time.month") - pr_clim
        pap = tmp.groupby("time.month") / pr_clim * 100
    else:
        precip_sr_data = precip_data.resample({time_dim: SeasonResampler([freq])}).mean(
            dim=time_dim
        )
        pr_sr_clim = (
            precip_sr_data.sel({time_dim: time_range}).groupby("time.month").mean()
        )
        tmp = precip_sr_data.groupby("time.month") - pr_sr_clim
        pap = tmp.groupby("time.month") / pr_sr_clim * 100

    pap.attrs["long_name"] = "Precipitation Anomaly Percentage"
    pap.attrs["units"] = "%"
    pap.attrs["time_range"] = f"{time_range.start}-{time_range.stop}"
    pap.attrs["method"] = (
        f"Grouped by {freq}, computed as (P - climatology)/climatology * 100"
    )

    return pap
