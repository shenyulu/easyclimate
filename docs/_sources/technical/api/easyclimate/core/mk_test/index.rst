:py:mod:`easyclimate.core.mk_test`
==================================

.. py:module:: easyclimate.core.mk_test

.. autoapi-nested-parse::

   Mann-Kendall Test

   The Mann-Kendall Trend Test (sometimes called the MK test) is used to analyze time series data for consistently increasing or decreasing trends (monotonic trends). It is a non-parametric test, which means it works for all distributions (i.e. data doesn't have to meet the assumption of normality), but data should have no serial correlation. If the data has a serial correlation, it could affect in significant level (p-value). It could lead to misinterpretation. To overcome this problem, researchers proposed several modified Mann-Kendall tests (Hamed and Rao Modified MK Test, Yue and Wang Modified MK Test, Modified MK test using Pre-Whitening method, etc.). Seasonal Mann-Kendall test also developed to remove the effect of seasonality.

   .. seealso::
       - https://github.com/mmhs013/pymannkendall
       - Hussain et al., (2019). pyMannKendall: a python package for non parametric Mann Kendall family of trend tests.. Journal of Open Source Software, 4(39), 1556, https://doi.org/10.21105/joss.01556



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.core.mk_test.original_test
   easyclimate.core.mk_test.hamed_rao_modification_test
   easyclimate.core.mk_test.yue_wang_modification_test
   easyclimate.core.mk_test.pre_whitening_modification_test
   easyclimate.core.mk_test.trend_free_pre_whitening_modification_test
   easyclimate.core.mk_test.seasonal_test
   easyclimate.core.mk_test.regional_test
   easyclimate.core.mk_test.correlated_seasonal_test
   easyclimate.core.mk_test.partial_test
   easyclimate.core.mk_test.sens_slope
   easyclimate.core.mk_test.seasonal_sens_slope



.. py:function:: original_test(data_input: xr.DataArray | xr.Dataset, dim: str, alpha: float = 0.05, returns_type: str = 'dataset_returns') -> xr.Dataset | DataTree

   Original Mann-Kendall test is a nonparametric test, which does not consider serial correlation or seasonal effects.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
   dim: str
       Dimension(s) over which to apply MK test.
   alpha: float, default 0.05.
       significance level (0.05 is the default).
   returns_type: str, default `dataset_returns`.
       Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

   Returns
   -------
   Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

   - trend: tells the trend (increasing, decreasing or no trend)
   - h: True (if trend is present) or False (if the trend is absence)
   - p: p-value of the significance test
   - z: normalized test statistics
   - Tau: Kendall Tau
   - s: Mann-Kendal's score
   - var_s: Variance S
   - slope: Theil-Sen estimator/slope
   - intercept: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.


.. py:function:: hamed_rao_modification_test(data_input: xr.DataArray | xr.Dataset, dim: str, alpha: float = 0.05, lag: int = None, returns_type: str = 'dataset_returns') -> xr.Dataset | DataTree

   This modified MK test proposed by Hamed and Rao (1998) to address serial autocorrelation issues. They suggested a variance correction approach to improve trend analysis. User can consider first n significant lag by insert lag number in this function. By default, it considered all significant lags.

   .. seealso::
       Hamed, K. H., & Rao, A. R. (1998). A modified Mann-Kendall trend test for autocorrelated data. Journal of hydrology, 204(1-4), 182-196. doi: http://doi.org/10.1016/S0022-1694(97)00125-X

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
   dim: str
       Dimension(s) over which to apply MK test.
   alpha: float, default 0.05.
       significance level (0.05 is the default).
   lag: int.
       No. of First Significant Lags
   returns_type: str, default `dataset_returns`.
       Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

   Returns
   -------
   Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

   - trend: tells the trend (increasing, decreasing or no trend)
   - h: True (if trend is present) or False (if the trend is absence)
   - p: p-value of the significance test
   - z: normalized test statistics
   - Tau: Kendall Tau
   - s: Mann-Kendal's score
   - var_s: Variance S
   - slope: Theil-Sen estimator/slope
   - intercept: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.


.. py:function:: yue_wang_modification_test(data_input: xr.DataArray | xr.Dataset, dim: str, alpha: float = 0.05, lag: int = None, returns_type: str = 'dataset_returns') -> xr.Dataset | DataTree

   This is also a variance correction method for considered serial autocorrelation proposed by Yue, S., & Wang, C. Y. (2004). User can also set their desired significant n lags for the calculation.

   .. seealso::
       Yue, S., & Wang, C. (2004). The Mann-Kendall test modified by effective sample size to detect trend in serially correlated hydrological series. Water resources management, 18(3), 201-218. doi: http://doi.org/10.1023/B:WARM.0000043140.61082.60

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
   dim: str
       Dimension(s) over which to apply MK test.
   alpha: float, default 0.05.
       significance level (0.05 is the default).
   lag: int.
       No. of First Significant Lags
   returns_type: str, default `dataset_returns`.
       Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

   Returns
   -------
   Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

   - trend: tells the trend (increasing, decreasing or no trend)
   - h: True (if trend is present) or False (if the trend is absence)
   - p: p-value of the significance test
   - z: normalized test statistics
   - Tau: Kendall Tau
   - s: Mann-Kendal's score
   - var_s: Variance S
   - slope: Theil-Sen estimator/slope
   - intercept: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.


.. py:function:: pre_whitening_modification_test(data_input: xr.DataArray | xr.Dataset, dim: str, alpha: float = 0.05, returns_type: str = 'dataset_returns') -> xr.Dataset | DataTree

   This test suggested by Yue and Wang (2002) to using Pre-Whitening the time series before the application of trend test.

   .. seealso::
       Yue, S., & Wang, C. Y. (2002). Applicability of prewhitening to eliminate the influence of serial correlation on the Mann-Kendall test. Water resources research, 38(6), 4-1. doi: http://doi.org/10.1029/2001WR000861

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
   dim: str
       Dimension(s) over which to apply MK test.
   alpha: float, default 0.05.
       significance level (0.05 is the default).
   returns_type: str, default `dataset_returns`.
       Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

   Returns
   -------
   Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

   - trend: tells the trend (increasing, decreasing or no trend)
   - h: True (if trend is present) or False (if the trend is absence)
   - p: p-value of the significance test
   - z: normalized test statistics
   - Tau: Kendall Tau
   - s: Mann-Kendal's score
   - var_s: Variance S
   - slope: Theil-Sen estimator/slope
   - intercept: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.


.. py:function:: trend_free_pre_whitening_modification_test(data_input: xr.DataArray | xr.Dataset, dim: str, alpha: float = 0.05, returns_type: str = 'dataset_returns') -> xr.Dataset | DataTree

   This test also proposed by Yue and Wang (2002) to remove trend component and then Pre-Whitening the time series before application of trend test.

   .. seealso::
       Yue, S., & Wang, C. Y. (2002). Applicability of prewhitening to eliminate the influence of serial correlation on the Mann-Kendall test. Water resources research, 38(6), 4-1. doi: http://doi.org/10.1029/2001WR000861

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
   dim: str
       Dimension(s) over which to apply MK test.
   alpha: float, default 0.05.
       significance level (0.05 is the default).
   returns_type: str, default `dataset_returns`.
       Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

   Returns
   -------
   Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

   - trend: tells the trend (increasing, decreasing or no trend)
   - h: True (if trend is present) or False (if the trend is absence)
   - p: p-value of the significance test
   - z: normalized test statistics
   - Tau: Kendall Tau
   - s: Mann-Kendal's score
   - var_s: Variance S
   - slope: Theil-Sen estimator/slope
   - intercept: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.


.. py:function:: seasonal_test(data_input: xr.DataArray | xr.Dataset, dim: str, alpha: float = 0.05, period: int = 12, returns_type: str = 'dataset_returns') -> xr.Dataset | DataTree

   For seasonal time series data, Hirsch, R.M., Slack, J.R. and Smith, R.A. (1982) proposed this test to calculate the seasonal trend.

   .. seealso::
       Hirsch, R. M., Slack, J. R., & Smith, R. A. (1982). Techniques of trend analysis for monthly water quality data. Water resources research, 18(1), 107-121. doi: http://doi.org/10.1029/WR018i001p00107

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
   dim: str
       Dimension(s) over which to apply MK test.
   alpha: float, default 0.05.
       Significance level (0.05 is the default).
   period: int, default 12.
       Seasonal cycle. For monthly data it is 12, weekly data it is 52.
   returns_type: str, default `dataset_returns`.
       Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

   Returns
   -------
   Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

   - trend: tells the trend (increasing, decreasing or no trend)
   - h: True (if trend is present) or False (if the trend is absence)
   - p: p-value of the significance test
   - z: normalized test statistics
   - Tau: Kendall Tau
   - s: Mann-Kendal's score
   - var_s: Variance S
   - slope: Theil-Sen estimator/slope
   - intercept: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.


.. py:function:: regional_test(data_input: xr.DataArray | xr.Dataset, dim: str, alpha: float = 0.05, returns_type: str = 'dataset_returns') -> xr.Dataset | DataTree

   Based on Hirsch (1982) proposed seasonal mk test, Helsel, D.R. and Frans, L.M., (2006) suggest regional mk test to calculate the overall trend in a regional scale.

   .. seealso::
       Hirsch, R. M., Slack, J. R., & Smith, R. A. (1982). Techniques of trend analysis for monthly water quality data. Water resources research, 18(1), 107-121. doi: http://doi.org/10.1029/WR018i001p00107

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
   dim: str
       Dimension(s) over which to apply MK test.
   alpha: float, default 0.05.
       significance level (0.05 is the default).
   returns_type: str, default `dataset_returns`.
       Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

   Returns
   -------
   Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

   - trend: tells the trend (increasing, decreasing or no trend)
   - h: True (if trend is present) or False (if the trend is absence)
   - p: p-value of the significance test
   - z: normalized test statistics
   - Tau: Kendall Tau
   - s: Mann-Kendal's score
   - var_s: Variance S
   - slope: Theil-Sen estimator/slope
   - intercept: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.


.. py:function:: correlated_seasonal_test(data_input: xr.DataArray | xr.Dataset, dim: str, alpha: float = 0.05, period: int = 12, returns_type: str = 'dataset_returns') -> xr.Dataset | DataTree

   This method proposed by Hipel (1994) used, when time series significantly correlated with the preceding one or more months/seasons.

   .. seealso::
       Hipel, K. W., & McLeod, A. I. (1994). Time series modelling of water resources and environmental systems (Vol. 45). Elsevier.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
   dim: str
       Dimension(s) over which to apply MK test.
   alpha: float, default 0.05.
       significance level (0.05 is the default).
   returns_type: str, default `dataset_returns`.
       Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

   Returns
   -------
   Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

   - trend: tells the trend (increasing, decreasing or no trend)
   - h: True (if trend is present) or False (if the trend is absence)
   - p: p-value of the significance test
   - z: normalized test statistics
   - Tau: Kendall Tau
   - s: Mann-Kendal's score
   - var_s: Variance S
   - slope: Theil-Sen estimator/slope
   - intercept: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.


.. py:function:: partial_test(data_input: xr.DataArray | xr.Dataset, dim: str, alpha: float = 0.05, returns_type: str = 'dataset_returns') -> xr.Dataset | DataTree

   In a real event, many factors are affecting the main studied response parameter, which can bias the trend results. To overcome this problem, Libiseller (2002) proposed this partial mk test. It required two parameters as input, where, one is response parameter and other is an independent parameter.

   .. seealso::
       Libiseller, C., & Grimvall, A. (2002). Performance of partial Mann-Kendall tests for trend detection in the presence of covariates. Environmetrics: The official journal of the International Environmetrics Society, 13(1), 71-84. doi: http://doi.org/10.1002/env.507

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
   dim: str
       Dimension(s) over which to apply MK test.
   alpha: float, default 0.05.
       significance level (0.05 is the default).
   returns_type: str, default `dataset_returns`.
       Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

   Returns
   -------
   Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

   - trend: tells the trend (increasing, decreasing or no trend)
   - h: True (if trend is present) or False (if the trend is absence)
   - p: p-value of the significance test
   - z: normalized test statistics
   - Tau: Kendall Tau
   - s: Mann-Kendal's score
   - var_s: Variance S
   - slope: Theil-Sen estimator/slope
   - intercept: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.


.. py:function:: sens_slope(data_input: xr.DataArray | xr.Dataset, dim: str, alpha: float = 0.05, returns_type: str = 'dataset_returns') -> xr.Dataset | DataTree

   This method proposed by Theil (1950) and Sen (1968) to estimate the magnitude of the monotonic trend. Intercept is calculate using Conover, W.J. (1980) method.

   .. seealso::
       - Theil, H. (1950). A rank-invariant method of linear and polynominal regression analysis (parts 1-3). In Ned. Akad. Wetensch. Proc. Ser. A (Vol. 53, pp. 1397-1412).
       - Sen, P. K. (1968). Estimates of the regression coefficient based on Kendall's tau. Journal of the American statistical association, 63(324), 1379-1389. doi: http://doi.org/10.1080/01621459.1968.10480934

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
   dim: str
       Dimension(s) over which to apply MK test.
   alpha: float, default 0.05.
       significance level (0.05 is the default).
   returns_type: str, default `dataset_returns`.
       Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

   Returns
   -------
   Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

   - trend: tells the trend (increasing, decreasing or no trend)
   - h: True (if trend is present) or False (if the trend is absence)
   - p: p-value of the significance test
   - z: normalized test statistics
   - Tau: Kendall Tau
   - s: Mann-Kendal's score
   - var_s: Variance S
   - slope: Theil-Sen estimator/slope
   - intercept: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.


.. py:function:: seasonal_sens_slope(data_input: xr.DataArray | xr.Dataset, dim: str, alpha: float = 0.05, period: int = 12, returns_type: str = 'dataset_returns') -> xr.Dataset | DataTree

   This method proposed by Hipel (1994) to estimate the magnitude of the monotonic trend, when data has seasonal effects. Intercept is calculate using Conover, W.J. (1980) method.

   .. seealso::
       Hipel, K. W., & McLeod, A. I. (1994). Time series modelling of water resources and environmental systems (Vol. 45). Elsevier.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
   dim: str
       Dimension(s) over which to apply MK test.
   alpha: float, default 0.05.
       significance level (0.05 is the default).
   returns_type: str, default `dataset_returns`.
       Returns the sorting type of the value. Sort by variable (`dataset_returns`) or by return value (`dataset_vars`).

   Returns
   -------
   Mann-Kendall Test results (:py:class:`xarray.Dataset <xarray.Dataset>` or :py:class:`DataTree <DataTree>`).

   - trend: tells the trend (increasing, decreasing or no trend)
   - h: True (if trend is present) or False (if the trend is absence)
   - p: p-value of the significance test
   - z: normalized test statistics
   - Tau: Kendall Tau
   - s: Mann-Kendal's score
   - var_s: Variance S
   - slope: Theil-Sen estimator/slope
   - intercept: intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step.


