"""
Mann-Kendall Test

The Mann-Kendall Trend Test (sometimes called the MK test) is used to analyze time series data for consistently increasing or decreasing trends (monotonic trends). It is a non-parametric test, which means it works for all distributions (i.e. data doesn't have to meet the assumption of normality), but data should have no serial correlation. If the data has a serial correlation, it could affect in significant level (p-value). It could lead to misinterpretation. To overcome this problem, researchers proposed several modified Mann-Kendall tests (Hamed and Rao Modified MK Test, Yue and Wang Modified MK Test, Modified MK test using Pre-Whitening method, etc.). Seasonal Mann-Kendall test also developed to remove the effect of seasonality.

.. seealso::
    - https://github.com/mmhs013/pymannkendall
    - Hussain et al., (2019). pyMannKendall: a python package for non parametric Mann Kendall family of trend tests.. Journal of Open Source Software, 4(39), 1556, https://doi.org/10.21105/joss.01556
"""

import numpy as np
import xarray as xr
import pymannkendall as mk

from .utility import generate_datatree_dispatcher
from datatree import DataTree

__all__ = [
    "original_test",
    "hamed_rao_modification_test",
    "yue_wang_modification_test",
    "pre_whitening_modification_test",
    "trend_free_pre_whitening_modification_test",
    "seasonal_test",
    "regional_test",
    "correlated_seasonal_test",
    "sens_slope",
    "seasonal_sens_slope",
]


@generate_datatree_dispatcher
def original_test(
    data_input: xr.DataArray | xr.Dataset,
    dim: str,
    alpha: float = 0.05,
    returns_type: str = "dataset_returns",
) -> xr.Dataset | DataTree:
    """
    Original Mann-Kendall test is a nonparametric test, which does not consider serial correlation or seasonal effects.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The spatio-temporal data to be calculated.
    dim: :py:class:`str <str>`
        Dimension(s) over which to apply MK test.
    alpha: :py:class:`float <float>`, default 0.05.
        significance level (0.05 is the default).
    returns_type: :py:class:`str <str>`, default `dataset_returns`.
        Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

    Returns
    -------
    Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

    - **trend**: tells the trend (increasing, decreasing or no trend)
    - **h**: True (if trend is present) or False (if the trend is absence)
    - **p**: p-value of the significance test
    - **z**: normalized test statistics
    - **Tau**: Kendall Tau
    - **s**: Mann-Kendal's score
    - **var_s**: Variance S
    - **slope**: Theil-Sen estimator/slope
    - **intercept**: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.
    """

    def _top_test(data_input, dim, alpha):
        def mk_test(data):
            mk_test_result = mk.original_test(data, alpha=alpha)

            trend = mk_test_result.trend
            if trend == "increasing":
                trend = 1
            elif trend == "decreasing":
                trend = -1
            elif trend == "no trend":
                trend = 0
            else:
                raise ValueError("Error `trend` type.")

            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims=[["parameter"]],
            output_dtypes=["float64"],
            dask="parallelized",
            vectorize=True,
            dask_gufunc_kwargs={
                "output_sizes": {"parameter": 1},
                "allow_rechunk": True,
            },
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars={
                "trend": mk_test_dataarray[..., 0],
                "h": mk_test_dataarray[..., 1],
                "p": mk_test_dataarray[..., 2],
                "z": mk_test_dataarray[..., 3],
                "Tau": mk_test_dataarray[..., 4],
                "s": mk_test_dataarray[..., 5],
                "var_s": mk_test_dataarray[..., 6],
                "slope": mk_test_dataarray[..., 7],
                "intercept": mk_test_dataarray[..., 8],
            }
        )

        mk_test_dataset["trend"].attrs[
            "Description"
        ] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset["h"].attrs[
            "Description"
        ] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset["p"].attrs["Description"] = "p-value of the significance test"
        mk_test_dataset["z"].attrs["Description"] = "normalized test statistics"
        mk_test_dataset["Tau"].attrs["Description"] = "Kendall Tau"
        mk_test_dataset["s"].attrs["Description"] = "Mann-Kendal's score"
        mk_test_dataset["var_s"].attrs["Description"] = "Variance S"
        mk_test_dataset["slope"].attrs["Description"] = "Theil-Sen estimator/slope"
        mk_test_dataset["intercept"].attrs[
            "Description"
        ] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs["Description"] = "original_test"
        return mk_test_dataset

    return _top_test(data_input, dim, alpha)


@generate_datatree_dispatcher
def hamed_rao_modification_test(
    data_input: xr.DataArray | xr.Dataset,
    dim: str,
    alpha: float = 0.05,
    lag: int = None,
    returns_type: str = "dataset_returns",
) -> xr.Dataset | DataTree:
    """
    This modified MK test proposed by Hamed and Rao (1998) to address serial autocorrelation issues. They suggested a variance correction approach to improve trend analysis. User can consider first n significant lag by insert lag number in this function. By default, it considered all significant lags.

    .. seealso::
        Hamed, K. H., & Rao, A. R. (1998). A modified Mann-Kendall trend test for autocorrelated data. Journal of hydrology, 204(1-4), 182-196. doi: http://doi.org/10.1016/S0022-1694(97)00125-X

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The spatio-temporal data to be calculated.
    dim: :py:class:`str <str>`
        Dimension(s) over which to apply MK test.
    alpha: :py:class:`float <float>`, default 0.05.
        significance level (0.05 is the default).
    lag: :py:class:`int <int>`.
        No. of First Significant Lags
    returns_type: :py:class:`str <str>`, default `dataset_returns`.
        Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

    Returns
    -------
    Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

    - **trend**: tells the trend (increasing, decreasing or no trend)
    - **h**: True (if trend is present) or False (if the trend is absence)
    - **p**: p-value of the significance test
    - **z**: normalized test statistics
    - **Tau**: Kendall Tau
    - **s**: Mann-Kendal's score
    - **var_s**: Variance S
    - **slope**: Theil-Sen estimator/slope
    - **intercept**: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.
    """

    def _top_test(data_input, dim, alpha, lag=lag):
        def mk_test(data):
            mk_test_result = mk.hamed_rao_modification_test(data, alpha=alpha, lag=lag)

            trend = mk_test_result.trend
            if trend == "increasing":
                trend = 1
            elif trend == "decreasing":
                trend = -1
            elif trend == "no trend":
                trend = 0
            else:
                raise ValueError("Error `trend` type.")

            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims=[["parameter"]],
            output_dtypes=["float64"],
            dask="parallelized",
            vectorize=True,
            dask_gufunc_kwargs={
                "output_sizes": {"parameter": 1},
                "allow_rechunk": True,
            },
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars={
                "trend": mk_test_dataarray[..., 0],
                "h": mk_test_dataarray[..., 1],
                "p": mk_test_dataarray[..., 2],
                "z": mk_test_dataarray[..., 3],
                "Tau": mk_test_dataarray[..., 4],
                "s": mk_test_dataarray[..., 5],
                "var_s": mk_test_dataarray[..., 6],
                "slope": mk_test_dataarray[..., 7],
                "intercept": mk_test_dataarray[..., 8],
            }
        )

        mk_test_dataset["trend"].attrs[
            "Description"
        ] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset["h"].attrs[
            "Description"
        ] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset["p"].attrs["Description"] = "p-value of the significance test"
        mk_test_dataset["z"].attrs["Description"] = "normalized test statistics"
        mk_test_dataset["Tau"].attrs["Description"] = "Kendall Tau"
        mk_test_dataset["s"].attrs["Description"] = "Mann-Kendal's score"
        mk_test_dataset["var_s"].attrs["Description"] = "Variance S"
        mk_test_dataset["slope"].attrs["Description"] = "Theil-Sen estimator/slope"
        mk_test_dataset["intercept"].attrs[
            "Description"
        ] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs["Description"] = "hamed_rao_modification_test"
        return mk_test_dataset

    return _top_test(data_input, dim, alpha, lag=lag)


@generate_datatree_dispatcher
def yue_wang_modification_test(
    data_input: xr.DataArray | xr.Dataset,
    dim: str,
    alpha: float = 0.05,
    lag: int = None,
    returns_type: str = "dataset_returns",
) -> xr.Dataset | DataTree:
    """
    This is also a variance correction method for considered serial autocorrelation proposed by Yue, S., & Wang, C. Y. (2004). User can also set their desired significant n lags for the calculation.

    .. seealso::
        Yue, S., & Wang, C. (2004). The Mann-Kendall test modified by effective sample size to detect trend in serially correlated hydrological series. Water resources management, 18(3), 201-218. doi: http://doi.org/10.1023/B:WARM.0000043140.61082.60

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The spatio-temporal data to be calculated.
    dim: :py:class:`str <str>`
        Dimension(s) over which to apply MK test.
    alpha: :py:class:`float <float>`, default 0.05.
        significance level (0.05 is the default).
    lag: :py:class:`int <int>`.
        No. of First Significant Lags
    returns_type: :py:class:`str <str>`, default `dataset_returns`.
        Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

    Returns
    -------
    Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

    - **trend**: tells the trend (increasing, decreasing or no trend)
    - **h**: True (if trend is present) or False (if the trend is absence)
    - **p**: p-value of the significance test
    - **z**: normalized test statistics
    - **Tau**: Kendall Tau
    - **s**: Mann-Kendal's score
    - **var_s**: Variance S
    - **slope**: Theil-Sen estimator/slope
    - **intercept**: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.
    """

    def _top_test(data_input, dim, alpha, lag=lag):
        def mk_test(data):
            mk_test_result = mk.yue_wang_modification_test(data, alpha=alpha, lag=lag)

            trend = mk_test_result.trend
            if trend == "increasing":
                trend = 1
            elif trend == "decreasing":
                trend = -1
            elif trend == "no trend":
                trend = 0
            else:
                raise ValueError("Error `trend` type.")

            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims=[["parameter"]],
            output_dtypes=["float64"],
            dask="parallelized",
            vectorize=True,
            dask_gufunc_kwargs={
                "output_sizes": {"parameter": 1},
                "allow_rechunk": True,
            },
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars={
                "trend": mk_test_dataarray[..., 0],
                "h": mk_test_dataarray[..., 1],
                "p": mk_test_dataarray[..., 2],
                "z": mk_test_dataarray[..., 3],
                "Tau": mk_test_dataarray[..., 4],
                "s": mk_test_dataarray[..., 5],
                "var_s": mk_test_dataarray[..., 6],
                "slope": mk_test_dataarray[..., 7],
                "intercept": mk_test_dataarray[..., 8],
            }
        )

        mk_test_dataset["trend"].attrs[
            "Description"
        ] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset["h"].attrs[
            "Description"
        ] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset["p"].attrs["Description"] = "p-value of the significance test"
        mk_test_dataset["z"].attrs["Description"] = "normalized test statistics"
        mk_test_dataset["Tau"].attrs["Description"] = "Kendall Tau"
        mk_test_dataset["s"].attrs["Description"] = "Mann-Kendal's score"
        mk_test_dataset["var_s"].attrs["Description"] = "Variance S"
        mk_test_dataset["slope"].attrs["Description"] = "Theil-Sen estimator/slope"
        mk_test_dataset["intercept"].attrs[
            "Description"
        ] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs["Description"] = "yue_wang_modification_test"
        return mk_test_dataset

    return _top_test(data_input, dim, alpha, lag=lag)


@generate_datatree_dispatcher
def pre_whitening_modification_test(
    data_input: xr.DataArray | xr.Dataset,
    dim: str,
    alpha: float = 0.05,
    returns_type: str = "dataset_returns",
) -> xr.Dataset | DataTree:
    """
    This test suggested by Yue and Wang (2002) to using Pre-Whitening the time series before the application of trend test.

    .. seealso::
        Yue, S., & Wang, C. Y. (2002). Applicability of prewhitening to eliminate the influence of serial correlation on the Mann-Kendall test. Water resources research, 38(6), 4-1. doi: http://doi.org/10.1029/2001WR000861

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The spatio-temporal data to be calculated.
    dim: :py:class:`str <str>`
        Dimension(s) over which to apply MK test.
    alpha: :py:class:`float <float>`, default 0.05.
        significance level (0.05 is the default).
    returns_type: :py:class:`str <str>`, default `dataset_returns`.
        Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

    Returns
    -------
    Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

    - **trend**: tells the trend (increasing, decreasing or no trend)
    - **h**: True (if trend is present) or False (if the trend is absence)
    - **p**: p-value of the significance test
    - **z**: normalized test statistics
    - **Tau**: Kendall Tau
    - **s**: Mann-Kendal's score
    - **var_s**: Variance S
    - **slope**: Theil-Sen estimator/slope
    - **intercept**: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.
    """

    def _top_test(data_input, dim, alpha):
        def mk_test(data):
            mk_test_result = mk.pre_whitening_modification_test(data, alpha=alpha)

            trend = mk_test_result.trend
            if trend == "increasing":
                trend = 1
            elif trend == "decreasing":
                trend = -1
            elif trend == "no trend":
                trend = 0
            else:
                raise ValueError("Error `trend` type.")

            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims=[["parameter"]],
            output_dtypes=["float64"],
            dask="parallelized",
            vectorize=True,
            dask_gufunc_kwargs={
                "output_sizes": {"parameter": 1},
                "allow_rechunk": True,
            },
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars={
                "trend": mk_test_dataarray[..., 0],
                "h": mk_test_dataarray[..., 1],
                "p": mk_test_dataarray[..., 2],
                "z": mk_test_dataarray[..., 3],
                "Tau": mk_test_dataarray[..., 4],
                "s": mk_test_dataarray[..., 5],
                "var_s": mk_test_dataarray[..., 6],
                "slope": mk_test_dataarray[..., 7],
                "intercept": mk_test_dataarray[..., 8],
            }
        )

        mk_test_dataset["trend"].attrs[
            "Description"
        ] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset["h"].attrs[
            "Description"
        ] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset["p"].attrs["Description"] = "p-value of the significance test"
        mk_test_dataset["z"].attrs["Description"] = "normalized test statistics"
        mk_test_dataset["Tau"].attrs["Description"] = "Kendall Tau"
        mk_test_dataset["s"].attrs["Description"] = "Mann-Kendal's score"
        mk_test_dataset["var_s"].attrs["Description"] = "Variance S"
        mk_test_dataset["slope"].attrs["Description"] = "Theil-Sen estimator/slope"
        mk_test_dataset["intercept"].attrs[
            "Description"
        ] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs["Description"] = "pre_whitening_modification_test"
        return mk_test_dataset

    return _top_test(data_input, dim, alpha)


@generate_datatree_dispatcher
def trend_free_pre_whitening_modification_test(
    data_input: xr.DataArray | xr.Dataset,
    dim: str,
    alpha: float = 0.05,
    returns_type: str = "dataset_returns",
) -> xr.Dataset | DataTree:
    """
    This test also proposed by Yue and Wang (2002) to remove trend component and then Pre-Whitening the time series before application of trend test.

    .. seealso::
        Yue, S., & Wang, C. Y. (2002). Applicability of prewhitening to eliminate the influence of serial correlation on the Mann-Kendall test. Water resources research, 38(6), 4-1. doi: http://doi.org/10.1029/2001WR000861

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The spatio-temporal data to be calculated.
    dim: :py:class:`str <str>`
        Dimension(s) over which to apply MK test.
    alpha: :py:class:`float <float>`, default 0.05.
        significance level (0.05 is the default).
    returns_type: :py:class:`str <str>`, default `dataset_returns`.
        Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

    Returns
    -------
    Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

    - **trend**: tells the trend (increasing, decreasing or no trend)
    - **h**: True (if trend is present) or False (if the trend is absence)
    - **p**: p-value of the significance test
    - **z**: normalized test statistics
    - **Tau**: Kendall Tau
    - **s**: Mann-Kendal's score
    - **var_s**: Variance S
    - **slope**: Theil-Sen estimator/slope
    - **intercept**: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.
    """

    def _top_test(data_input, dim, alpha):
        def mk_test(data):
            mk_test_result = mk.trend_free_pre_whitening_modification_test(
                data, alpha=alpha
            )

            trend = mk_test_result.trend
            if trend == "increasing":
                trend = 1
            elif trend == "decreasing":
                trend = -1
            elif trend == "no trend":
                trend = 0
            else:
                raise ValueError("Error `trend` type.")

            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims=[["parameter"]],
            output_dtypes=["float64"],
            dask="parallelized",
            vectorize=True,
            dask_gufunc_kwargs={
                "output_sizes": {"parameter": 1},
                "allow_rechunk": True,
            },
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars={
                "trend": mk_test_dataarray[..., 0],
                "h": mk_test_dataarray[..., 1],
                "p": mk_test_dataarray[..., 2],
                "z": mk_test_dataarray[..., 3],
                "Tau": mk_test_dataarray[..., 4],
                "s": mk_test_dataarray[..., 5],
                "var_s": mk_test_dataarray[..., 6],
                "slope": mk_test_dataarray[..., 7],
                "intercept": mk_test_dataarray[..., 8],
            }
        )

        mk_test_dataset["trend"].attrs[
            "Description"
        ] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset["h"].attrs[
            "Description"
        ] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset["p"].attrs["Description"] = "p-value of the significance test"
        mk_test_dataset["z"].attrs["Description"] = "normalized test statistics"
        mk_test_dataset["Tau"].attrs["Description"] = "Kendall Tau"
        mk_test_dataset["s"].attrs["Description"] = "Mann-Kendal's score"
        mk_test_dataset["var_s"].attrs["Description"] = "Variance S"
        mk_test_dataset["slope"].attrs["Description"] = "Theil-Sen estimator/slope"
        mk_test_dataset["intercept"].attrs[
            "Description"
        ] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs["Description"] = (
            "trend_free_pre_whitening_modification_test"
        )
        return mk_test_dataset

    return _top_test(data_input, dim, alpha)


@generate_datatree_dispatcher
def seasonal_test(
    data_input: xr.DataArray | xr.Dataset,
    dim: str,
    alpha: float = 0.05,
    period: int = 12,
    returns_type: str = "dataset_returns",
) -> xr.Dataset | DataTree:
    """
    For seasonal time series data, Hirsch, R.M., Slack, J.R. and Smith, R.A. (1982) proposed this test to calculate the seasonal trend.

    .. seealso::
        Hirsch, R. M., Slack, J. R., & Smith, R. A. (1982). Techniques of trend analysis for monthly water quality data. Water resources research, 18(1), 107-121. doi: http://doi.org/10.1029/WR018i001p00107

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The spatio-temporal data to be calculated.
    dim: :py:class:`str <str>`
        Dimension(s) over which to apply MK test.
    alpha: :py:class:`float <float>`, default 0.05.
        Significance level (0.05 is the default).
    period: :py:class:`int <int>`, default 12.
        Seasonal cycle. For monthly data it is 12, weekly data it is 52.
    returns_type: :py:class:`str <str>`, default `dataset_returns`.
        Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

    Returns
    -------
    Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

    - **trend**: tells the trend (increasing, decreasing or no trend)
    - **h**: True (if trend is present) or False (if the trend is absence)
    - **p**: p-value of the significance test
    - **z**: normalized test statistics
    - **Tau**: Kendall Tau
    - **s**: Mann-Kendal's score
    - **var_s**: Variance S
    - **slope**: Theil-Sen estimator/slope
    - **intercept**: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.
    """

    def _top_test(data_input, dim, alpha, period):
        def mk_test(data):
            mk_test_result = mk.seasonal_test(data, alpha=alpha, period=period)

            trend = mk_test_result.trend
            if trend == "increasing":
                trend = 1
            elif trend == "decreasing":
                trend = -1
            elif trend == "no trend":
                trend = 0
            else:
                raise ValueError("Error `trend` type.")

            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims=[["parameter"]],
            output_dtypes=["float64"],
            dask="parallelized",
            vectorize=True,
            dask_gufunc_kwargs={
                "output_sizes": {"parameter": 1},
                "allow_rechunk": True,
            },
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars={
                "trend": mk_test_dataarray[..., 0],
                "h": mk_test_dataarray[..., 1],
                "p": mk_test_dataarray[..., 2],
                "z": mk_test_dataarray[..., 3],
                "Tau": mk_test_dataarray[..., 4],
                "s": mk_test_dataarray[..., 5],
                "var_s": mk_test_dataarray[..., 6],
                "slope": mk_test_dataarray[..., 7],
                "intercept": mk_test_dataarray[..., 8],
            }
        )

        mk_test_dataset["trend"].attrs[
            "Description"
        ] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset["h"].attrs[
            "Description"
        ] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset["p"].attrs["Description"] = "p-value of the significance test"
        mk_test_dataset["z"].attrs["Description"] = "normalized test statistics"
        mk_test_dataset["Tau"].attrs["Description"] = "Kendall Tau"
        mk_test_dataset["s"].attrs["Description"] = "Mann-Kendal's score"
        mk_test_dataset["var_s"].attrs["Description"] = "Variance S"
        mk_test_dataset["slope"].attrs["Description"] = "Theil-Sen estimator/slope"
        mk_test_dataset["intercept"].attrs[
            "Description"
        ] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs["Description"] = "seasonal_test"
        return mk_test_dataset

    return _top_test(data_input, dim, alpha, period)


@generate_datatree_dispatcher
def regional_test(
    data_input: xr.DataArray | xr.Dataset,
    dim: str,
    alpha: float = 0.05,
    returns_type: str = "dataset_returns",
) -> xr.Dataset | DataTree:
    """
    Based on Hirsch (1982) proposed seasonal mk test, Helsel, D.R. and Frans, L.M., (2006) suggest regional mk test to calculate the overall trend in a regional scale.

    .. seealso::
        Hirsch, R. M., Slack, J. R., & Smith, R. A. (1982). Techniques of trend analysis for monthly water quality data. Water resources research, 18(1), 107-121. doi: http://doi.org/10.1029/WR018i001p00107

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The spatio-temporal data to be calculated.
    dim: :py:class:`str <str>`
        Dimension(s) over which to apply MK test.
    alpha: :py:class:`float <float>`, default 0.05.
        significance level (0.05 is the default).
    returns_type: :py:class:`str <str>`, default `dataset_returns`.
        Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

    Returns
    -------
    Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

    - **trend**: tells the trend (increasing, decreasing or no trend)
    - **h**: True (if trend is present) or False (if the trend is absence)
    - **p**: p-value of the significance test
    - **z**: normalized test statistics
    - **Tau**: Kendall Tau
    - **s**: Mann-Kendal's score
    - **var_s**: Variance S
    - **slope**: Theil-Sen estimator/slope
    - **intercept**: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.
    """

    def _top_test(data_input, dim, alpha):
        def mk_test(data):
            mk_test_result = mk.regional_test(data, alpha=alpha)

            trend = mk_test_result.trend
            if trend == "increasing":
                trend = 1
            elif trend == "decreasing":
                trend = -1
            elif trend == "no trend":
                trend = 0
            else:
                raise ValueError("Error `trend` type.")

            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims=[["parameter"]],
            output_dtypes=["float64"],
            dask="parallelized",
            vectorize=True,
            dask_gufunc_kwargs={
                "output_sizes": {"parameter": 1},
                "allow_rechunk": True,
            },
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars={
                "trend": mk_test_dataarray[..., 0],
                "h": mk_test_dataarray[..., 1],
                "p": mk_test_dataarray[..., 2],
                "z": mk_test_dataarray[..., 3],
                "Tau": mk_test_dataarray[..., 4],
                "s": mk_test_dataarray[..., 5],
                "var_s": mk_test_dataarray[..., 6],
                "slope": mk_test_dataarray[..., 7],
                "intercept": mk_test_dataarray[..., 8],
            }
        )

        mk_test_dataset["trend"].attrs[
            "Description"
        ] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset["h"].attrs[
            "Description"
        ] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset["p"].attrs["Description"] = "p-value of the significance test"
        mk_test_dataset["z"].attrs["Description"] = "normalized test statistics"
        mk_test_dataset["Tau"].attrs["Description"] = "Kendall Tau"
        mk_test_dataset["s"].attrs["Description"] = "Mann-Kendal's score"
        mk_test_dataset["var_s"].attrs["Description"] = "Variance S"
        mk_test_dataset["slope"].attrs["Description"] = "Theil-Sen estimator/slope"
        mk_test_dataset["intercept"].attrs[
            "Description"
        ] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs["Description"] = "regional_test"
        return mk_test_dataset

    return _top_test(data_input, dim, alpha)


@generate_datatree_dispatcher
def correlated_seasonal_test(
    data_input: xr.DataArray | xr.Dataset,
    dim: str,
    alpha: float = 0.05,
    period: int = 12,
    returns_type: str = "dataset_returns",
) -> xr.Dataset | DataTree:
    """
    This method proposed by Hipel (1994) used, when time series significantly correlated with the preceding one or more months/seasons.

    .. seealso::
        Hipel, K. W., & McLeod, A. I. (1994). Time series modelling of water resources and environmental systems (Vol. 45). Elsevier.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The spatio-temporal data to be calculated.
    dim: :py:class:`str <str>`
        Dimension(s) over which to apply MK test.
    alpha: :py:class:`float <float>`, default 0.05.
        significance level (0.05 is the default).
    period: :py:class:`int <int>`, default 12.
        Seasonal cycle. For monthly data it is 12, weekly data it is 52.
    returns_type: :py:class:`str <str>`, default `dataset_returns`.
        Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

    Returns
    -------
    Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

    - **trend**: tells the trend (increasing, decreasing or no trend)
    - **h**: True (if trend is present) or False (if the trend is absence)
    - **p**: p-value of the significance test
    - **z**: normalized test statistics
    - **Tau**: Kendall Tau
    - **s**: Mann-Kendal's score
    - **var_s**: Variance S
    - **slope**: Theil-Sen estimator/slope
    - **intercept**: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.
    """

    def _top_test(data_input, dim, alpha, period):
        def mk_test(data):
            mk_test_result = mk.correlated_seasonal_test(
                data, alpha=alpha, period=period
            )

            trend = mk_test_result.trend
            if trend == "increasing":
                trend = 1
            elif trend == "decreasing":
                trend = -1
            elif trend == "no trend":
                trend = 0
            else:
                raise ValueError("Error `trend` type.")

            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims=[["parameter"]],
            output_dtypes=["float64"],
            dask="parallelized",
            vectorize=True,
            dask_gufunc_kwargs={
                "output_sizes": {"parameter": 1},
                "allow_rechunk": True,
            },
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars={
                "trend": mk_test_dataarray[..., 0],
                "h": mk_test_dataarray[..., 1],
                "p": mk_test_dataarray[..., 2],
                "z": mk_test_dataarray[..., 3],
                "Tau": mk_test_dataarray[..., 4],
                "s": mk_test_dataarray[..., 5],
                "var_s": mk_test_dataarray[..., 6],
                "slope": mk_test_dataarray[..., 7],
                "intercept": mk_test_dataarray[..., 8],
            }
        )

        mk_test_dataset["trend"].attrs[
            "Description"
        ] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset["h"].attrs[
            "Description"
        ] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset["p"].attrs["Description"] = "p-value of the significance test"
        mk_test_dataset["z"].attrs["Description"] = "normalized test statistics"
        mk_test_dataset["Tau"].attrs["Description"] = "Kendall Tau"
        mk_test_dataset["s"].attrs["Description"] = "Mann-Kendal's score"
        mk_test_dataset["var_s"].attrs["Description"] = "Variance S"
        mk_test_dataset["slope"].attrs["Description"] = "Theil-Sen estimator/slope"
        mk_test_dataset["intercept"].attrs[
            "Description"
        ] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs["Description"] = "correlated_seasonal_test"
        return mk_test_dataset

    return _top_test(data_input, dim, alpha, period)


@generate_datatree_dispatcher
def sens_slope(
    data_input: xr.DataArray | xr.Dataset,
    dim: str,
    returns_type: str = "dataset_returns",
) -> xr.Dataset | DataTree:
    """
    This method proposed by Theil (1950) and Sen (1968) to estimate the magnitude of the monotonic trend. Intercept is calculate using Conover, W.J. (1980) method.

    .. seealso::
        - Theil, H. (1950). A rank-invariant method of linear and polynominal regression analysis (parts 1-3). In Ned. Akad. Wetensch. Proc. Ser. A (Vol. 53, pp. 1397-1412).
        - Sen, P. K. (1968). Estimates of the regression coefficient based on Kendall's tau. Journal of the American statistical association, 63(324), 1379-1389. doi: http://doi.org/10.1080/01621459.1968.10480934

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The spatio-temporal data to be calculated.
    dim: :py:class:`str <str>`
        Dimension(s) over which to apply MK test.
    alpha: :py:class:`float <float>`, default 0.05.
        significance level (0.05 is the default).
    returns_type: :py:class:`str <str>`, default `dataset_returns`.
        Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

    Returns
    -------
    Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

    - **slope**: Theil-Sen estimator/slope
    - **intercept**: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.
    """

    def _top_test(data_input, dim):
        def mk_test(data):
            mk_test_result = mk.sens_slope(data)

            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims=[["parameter"]],
            output_dtypes=["float64"],
            dask="parallelized",
            vectorize=True,
            dask_gufunc_kwargs={
                "output_sizes": {"parameter": 1},
                "allow_rechunk": True,
            },
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars={
                "slope": mk_test_dataarray[..., 0],
                "intercept": mk_test_dataarray[..., 1],
            }
        )

        mk_test_dataset["slope"].attrs["Description"] = "Theil-Sen estimator/slope"
        mk_test_dataset["intercept"].attrs[
            "Description"
        ] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs["Description"] = "sens_slope"
        return mk_test_dataset

    return _top_test(data_input, dim)


@generate_datatree_dispatcher
def seasonal_sens_slope(
    data_input: xr.DataArray | xr.Dataset,
    dim: str,
    period: int = 12,
    returns_type: str = "dataset_returns",
) -> xr.Dataset | DataTree:
    """
    This method proposed by Hipel (1994) to estimate the magnitude of the monotonic trend, when data has seasonal effects. Intercept is calculate using Conover, W.J. (1980) method.

    .. seealso::
        Hipel, K. W., & McLeod, A. I. (1994). Time series modelling of water resources and environmental systems (Vol. 45). Elsevier.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The spatio-temporal data to be calculated.
    dim: :py:class:`str <str>`
        Dimension(s) over which to apply MK test.
    period: :py:class:`int <int>`, default 12.
        Seasonal cycle. For monthly data it is 12, weekly data it is 52.
    returns_type: :py:class:`str <str>`, default `dataset_returns`.
        Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

    Returns
    -------
    Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

    - **trend**: tells the trend (increasing, decreasing or no trend)
    - **h**: True (if trend is present) or False (if the trend is absence)
    - **p**: p-value of the significance test
    - **z**: normalized test statistics
    - **Tau**: Kendall Tau
    - **s**: Mann-Kendal's score
    - **var_s**: Variance S
    - **slope**: Theil-Sen estimator/slope
    - **intercept**: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.
    """

    def _top_test(data_input, dim, period):
        def mk_test(data):
            mk_test_result = mk.seasonal_sens_slope(data, period=period)

            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims=[["parameter"]],
            output_dtypes=["float64"],
            dask="parallelized",
            vectorize=True,
            dask_gufunc_kwargs={
                "output_sizes": {"parameter": 1},
                "allow_rechunk": True,
            },
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars={
                "slope": mk_test_dataarray[..., 0],
                "intercept": mk_test_dataarray[..., 1],
            }
        )

        mk_test_dataset["slope"].attrs["Description"] = "Theil-Sen estimator/slope"
        mk_test_dataset["intercept"].attrs[
            "Description"
        ] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs["Description"] = "seasonal_sens_slope"
        return mk_test_dataset

    return _top_test(data_input, dim, period)
