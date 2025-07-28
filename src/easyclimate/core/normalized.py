"""
Normalized Data
"""

import xarray as xr

__all__ = ["normalize_zscore", "normalize_minmax", "normalize_robust", "normalize_mean"]


def normalize_zscore(
    da: xr.DataArray,
    dim: str = "time",
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
    """
    mean = da.mean(dim=dim)
    std = da.std(dim=dim, ddof=ddof)
    return (da - mean) / std


def normalize_minmax(
    da: xr.DataArray,
    dim: str = "time",
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
    min_val = da.min(dim=dim)
    max_val = da.max(dim=dim)
    a, b = feature_range
    return (da - min_val) / (max_val - min_val) * (b - a) + a


def normalize_robust(
    da: xr.DataArray,
    dim: str = "time",
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
    da: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The input time series data to be standardized.
    dim: :py:class:`str <str>`, default: `time`.
        The dimension along which to compute the median and IQR. By default, standardization is applied over the `time` dimension.
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
    median = da.median(dim=dim)
    q75 = da.quantile(q_high, dim=dim)
    q25 = da.quantile(q_low, dim=dim)
    iqr = q75 - q25
    return (da - median) / iqr


def normalize_mean(
    da: xr.DataArray,
    dim: str = "time",
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

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
        The normalized data centered around zero.

    .. note ::

        - **Applicable Scenarios**: Useful when centering data is needed without enforcing a standard deviation of 1.
        - **Advantages**: Simple, partially preserves data distribution characteristics.
        - **Disadvantages**: Sensitive to outliers, and the scaling range is not fixed.
    """
    mean = da.mean(dim=dim)
    min_val = da.min(dim=dim)
    max_val = da.max(dim=dim)
    return (da - mean) / (max_val - min_val)
