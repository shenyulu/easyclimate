easyclimate.core.stat
=====================

.. py:module:: easyclimate.core.stat

.. autoapi-nested-parse::

   Basic statistical analysis of weather and climate variables



Functions
---------

.. autoapisummary::

   easyclimate.core.stat.calc_linregress_spatial
   easyclimate.core.stat.calc_detrend_spatial
   easyclimate.core.stat.calc_corr_spatial
   easyclimate.core.stat.calc_leadlag_corr_spatial
   easyclimate.core.stat.calc_multiple_linear_regression_spatial
   easyclimate.core.stat.calc_ttestSpatialPattern_spatial
   easyclimate.core.stat.calc_levenetestSpatialPattern_spatial
   easyclimate.core.stat.calc_levenetestSpatialPattern_spatial
   easyclimate.core.stat.calc_skewness_spatial
   easyclimate.core.stat.calc_kurtosis_spatial
   easyclimate.core.stat.calc_theilslopes_spatial
   easyclimate.core.stat.calc_lead_lag_correlation_coefficients
   easyclimate.core.stat.calc_timeseries_correlations


Module Contents
---------------

.. py:function:: calc_linregress_spatial(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', x: numpy.array = None, alternative: str = 'two-sided', returns_type: {'dataset_returns', 'dataset_vars'} = 'dataset_returns') -> xarray.Dataset | xarray.DataTree

   Calculate a linear least-squares regression (**trend**) for spatial data of time.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>` to be regression.
   dim: :py:class:`str <str>`, default `time`.
       Dimension(s) over which to apply linregress. By default linregress is applied over the `time` dimension.
   x: :py:class:`numpy.array <numpy.array>`
       Independent variable. If None, use `np.arange(len(data_input['time'].shape[0]))` instead.
   returns_type: :py:class:`str <str>`, default `'dataset_returns'`.
       Return data type.

   Returns
   -------
   result : ``LinregressResult`` Dataset
       The return Dataset have following data_var:

       **slope**: :py:class:`float <float>`
           Slope of the regression line.
       **intercept**: :py:class:`float <float>`
           Intercept of the regression line.
       **rvalue**: :py:class:`float <float>`
           The Pearson correlation coefficient. The square of ``rvalue``
           is equal to the coefficient of determination.
       **pvalue**: :py:class:`float <float>`
           The p-value for a hypothesis test whose null hypothesis is
           that the slope is zero, using Wald Test with t-distribution of
           the test statistic. See `alternative` above for alternative
           hypotheses.
       **stderr**: :py:class:`float <float>`
           Standard error of the estimated slope (gradient), under the
           assumption of residual normality.
       **intercept_stderr**: :py:class:`float <float>`
           Standard error of the estimated intercept, under the assumption
           of residual normality.

   .. seealso::
       :py:func:`scipy.stats.linregress <scipy:scipy.stats.linregress>`.


.. py:function:: calc_detrend_spatial(data_input: xarray.DataArray | xarray.Dataset, time_dim: str = 'time') -> xarray.DataArray | xarray.DataTree

   Remove linear trend along axis from data.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`
        The spatio-temporal data of :py:class:`xarray.DataArray<xarray.DataArray>` to be detrended.
   time_dim: :py:class:`str <str>`
       Dimension(s) over which to detrend. By default dimension is applied over the `time` dimension.

   Returns
   -------
   - :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::
       :py:func:`scipy.signal.detrend <scipy:scipy.signal.detrend>`.


.. py:function:: calc_corr_spatial(data_input: xarray.DataArray, x: xarray.DataArray | numpy.ndarray, time_dim: str = 'time', method: Literal['scipy', 'xarray'] = 'xarray') -> xarray.Dataset

   Calculate Pearson correlation coefficients and corresponding p-values between spatial data
   and a time series using ``scipy.stats.pearsonr``.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       Input spatial data with dimensions ``(time, ...)``.

       .. note::

           NaN values are automatically skipped in calculations.

   x : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`numpy.ndarray<numpy.ndarray>`
       Time series data with dimension ``(time,)``. Must have the same length as data_input's time dimension.

       .. note::

           NaN values are automatically skipped in calculations.

   time_dim: :py:class:`str <str>`
       Dimension(s) over which to detrend. By default dimension is applied over the `time` dimension.
   method : {'scipy', 'xarray'}, optional
       Method used to compute correlations:

       - 'scipy': Uses :py:func:`scipy.stats.pearsonr<scipy:scipy.stats.pearsonr>` for direct calculation
       - 'xarray': Uses xarray's built-in correlation with t-test conversion (faster)

       Default is 'xarray'.

   Returns
   -------
   corr & pvalue (:py:class:`xarray.Dataset<xarray.Dataset>`)

   corr : :py:class:`xarray.DataArray<xarray.DataArray>`
       Pearson correlation coefficients with dimensions.
       Values range from -1 to 1 where:

       - 1 means perfect positive correlation
       - -1 means perfect negative correlation
       - 0 means no correlation

   pvalue : :py:class:`xarray.DataArray<xarray.DataArray>`
       Two-tailed p-values with dimensions.
       Small p-values (<0.05) indicate statistically significant correlations.

   Examples
   --------
   >>> data_input = xr.DataArray(np.random.rand(10, 3, 4),
   ...                           dims=['time', 'lat', 'lon'],
   ...                           coords={'time': pd.date_range('2000-01-01', periods=10)})
   >>> x = xr.DataArray(np.random.rand(10), dims=['time'])
   >>> corr_dataset = ecl.calc_corr_spatial(data_input, x)

   .. seealso::
       :py:func:`scipy.stats.pearsonr<scipy:scipy.stats.pearsonr>`:
       The underlying correlation function used for calculations.


.. py:function:: calc_leadlag_corr_spatial(data_input: xarray.DataArray, x: xarray.DataArray | numpy.ndarray, leadlag_array: numpy.array | List[int], time_dim: str = 'time', method: Literal['scipy', 'xarray'] = 'xarray')

   Calculate Pearson correlation coefficients and corresponding p-values between spatial data
   and a time series with specified lead or lag shifts, using ``scipy.stats.pearsonr`` or xarray methods.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       Input spatial data with dimensions ``(time, ...)`` representing spatial fields over time.

       .. note::
           NaN values are automatically skipped in calculations.

   x : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`numpy.ndarray<numpy.ndarray>`
       Time series data with dimension ``(time,)``. Must have the same length as ``data_input``'s time dimension.

       .. note::
           NaN values are automatically skipped in calculations.

   leadlag_array : :py:class:`numpy.ndarray<numpy.ndarray>` or :py:class:`List[int]<list>`
       Array or list of integers specifying the lead or lag shifts (in time steps) to apply to the time series `x`
       relative to `data_input`.

       - **Positive values** indicate a **lag**: the time series `x` is shifted forward in time (e.g., a value of +2 means `x` is delayed by 2 time steps relative to `data_input`).
       - **Negative values** indicate a **lead**: the time series `x` is shifted backward in time (e.g., a value of -2 means `x` leads `data_input` by 2 time steps).
       - A value of **0** means no shift (synchronous correlation).

       Example: If ``leadlag_array = [-2, 0, 2]``, correlations are computed for :math:`x` leading by 2 time steps, no shift, and lagging by 2 time steps, respectively.

   time_dim : :py:class:`str<str>`
       Name of the time dimension in `data_input` and `x`. Default is `"time"`.

   method : {'scipy', 'xarray'}, optional
       Method used to compute correlations:

       - `'scipy'`: Uses :py:func:`scipy.stats.pearsonr<scipy:scipy.stats.pearsonr>` for direct calculation, which may be more precise but slower.
       - `'xarray'`: Uses xarray's built-in correlation function with t-test conversion, which is typically faster.

       Default is `'xarray'`.

   Returns
   -------
   result : :py:class:`xarray.Dataset<xarray.Dataset>`
       Dataset containing two variables:

       - **corr** : :py:class:`xarray.DataArray<xarray.DataArray>`
           Pearson correlation coefficients with dimensions ``(leadlag, ...)``.
           Values range from -1 to 1, where:
           - 1 indicates a perfect positive correlation.
           - -1 indicates a perfect negative correlation.
           - 0 indicates no correlation.

       - **pvalue** : :py:class:`xarray.DataArray<xarray.DataArray>`
           Two-tailed p-values with dimensions ``(leadlag, ...)``.
           Small p-values (<0.05) indicate statistically significant correlations.

   Notes
   -----
   - The function iterates over each lead/lag value in `leadlag_array`, computes the correlation between the shifted `x` and `data_input`, and concatenates results along a new `leadlag` dimension.
   - Shifting `x` may introduce NaN values at the edges of the time series, which are handled automatically during correlation calculations.
   - Ensure `data_input` and `x` have compatible time dimensions to avoid errors.

   Examples
   --------
   >>> import xarray as xr
   >>> import numpy as np
   >>> data = xr.DataArray(np.random.rand(100, 10, 10), dims=["time", "lat", "lon"])
   >>> ts = xr.DataArray(np.random.rand(100), dims=["time"])
   >>> leadlag = [-2, 0, 2]
   >>> result = calc_leadlag_corr_spatial(data, ts, leadlag, time_dim="time", method="xarray")
   >>> print(result)
   Processing leadlag: 2 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
   <xarray.Dataset> Size: 7kB
   Dimensions:    (leadlag: 3, lat: 10, lon: 10)
   Coordinates:
   * leadlag    (leadlag) int64 24B -2 0 2
   Dimensions without coordinates: lat, lon
   Data variables:
       reg_coeff  (leadlag, lat, lon) float64 2kB 0.006322 0.002647 ... -0.02781
       corr       (leadlag, lat, lon) float64 2kB 0.02141 0.00894 ... -0.09169
       pvalue     (leadlag, lat, lon) float64 2kB 0.8326 0.9297 ... 0.3053 0.3643


.. py:function:: calc_multiple_linear_regression_spatial(y_data: xarray.DataArray, x_datas: List[xarray.DataArray], dim='time') -> xarray.Dataset

   Apply multiple linear regression to dataset across spatial dimensions.

   .. math::

       y = a_1 x_1 + a_2 x_2 + \cdots

   Parameters
   -----------
   y_data : :py:class:`xarray.DataArray<xarray.DataArray>`
       Dependent variable with dimensions, each with dimensions ``(time,)``.
   x_datas : :py:class:`list <list>` of :py:class:`xarray.DataArray<xarray.DataArray>`
       List of independent variables, each with dimensions ``(time,)``.
   dim : :py:class:`str <str>`, optional
       Time dimension name (default: ``'time'``)

   Returns
   --------
   :py:class:`xarray.Dataset <xarray.Dataset>`
       :py:class:`xarray.Dataset <xarray.Dataset>` containing regression results with:

       - slopes: slope coefficients for each predictor ``(coef, lat, lon)``
       - intercept: intercept values ``(lat, lon)``
       - r_squared: coefficient of determination ``(lat, lon)``
       - slopes_p: p-values for slope coefficients ``(coef, lat, lon)``
       - intercept_p: p-values for intercept ``(lat, lon)``

   Raises
   -------
   ValueError
       If the time coordinates of input variables don't match.


.. py:function:: calc_ttestSpatialPattern_spatial(data_input1: xarray.DataArray, data_input2: xarray.DataArray, dim: str = 'time') -> xarray.Dataset

   Calculate the T-test for the means of two independent sptial samples along with other axis (i.e. 'time') of scores.

   Parameters
   ----------
   data_input1: :py:class:`xarray.DataArray<xarray.DataArray>`
        The first spatio-temporal data of xarray DataArray to be calculated.
   data_input2: :py:class:`xarray.DataArray<xarray.DataArray>`
        The second spatio-temporal data of xarray DataArray to be calculated.

   .. note::
       - The order of `data_input1` and `data_input2` has no effect on the calculation result.
       - The non-time dimensions of the two data sets must be exactly the same, and the dimensionality values must be arranged in the same order (ascending or descending).

   dim: :py:class:`str <str>`
       Dimension(s) over which to apply the test. By default the test is applied over the `time` dimension.

   Returns
   -------
   - **statistic**, **pvalue**: :py:class:`xarray.Dataset<xarray.Dataset>`.

   .. seealso::
       :py:func:`scipy.stats.ttest_ind <scipy:scipy.stats.ttest_ind>`.


.. py:function:: calc_levenetestSpatialPattern_spatial(data_input1: xarray.DataArray, data_input2: xarray.DataArray, dim: str = 'time', center: {'mean', 'median', 'trimmed'} = 'median', proportiontocut: float = 0.05) -> xarray.Dataset

   Perform Levene test for equal variances of two independent sptial samples along with other axis (i.e. 'time') of scores.

   The Levene test tests the null hypothesis that all input samples are from populations with equal variances.
   Levene's test is an alternative to Bartlett's test in the case where there are significant deviations from normality.

   Parameters
   ----------
   data_input1: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The first spatio-temporal data of xarray DataArray to be calculated.
   data_input2: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The second spatio-temporal data of xarray DataArray to be calculated.

   .. note::
       - The order of `data_input1` and `data_input2` has no effect on the calculation result.
       - The non-time dimensions of the two data sets must be exactly the same, and the dimensionality values must be arranged in the same order (ascending or descending).

   dim: :py:class:`str <str>`.
       Dimension(s) over which to apply the test. By default the test is applied over the `time` dimension.
   center: {'mean', 'median', 'trimmed'}, default `'median'`.
       Which function of the data to use in the test.

       .. note::

           Three variations of Levene’s test are possible. The possibilities and their recommended usages are:

           - median: Recommended for skewed (non-normal) distributions.
           - mean: Recommended for symmetric, moderate-tailed distributions.
           - trimmed: Recommended for heavy-tailed distributions.

           The test version using the mean was proposed in the original article of Levene (Levene, H., 1960) while the median and trimmed mean have been studied by Brown and Forsythe (Brown, M. B. and Forsythe, A. B., 1974), sometimes also referred to as Brown-Forsythe test.


   proportiontocut: :py:class:`float <float>`, default `0.05`.
       When center is `'trimmed'`, this gives the proportion of data points to cut from each end (See :py:func:`scipy.stats.trim_mean <scipy:scipy.stats.trim_mean>`).

   Returns
   -------
   - **statistic**, **pvalue**: :py:class:`xarray.Dataset<xarray.Dataset>`.

   Reference
   --------------
   - Levene, H. (1960). In Contributions to Probability and Statistics: Essays in Honor of Harold Hotelling, I. Olkin et al. eds., Stanford University Press, pp. 278-292.
   - Morton B. Brown & Alan B. Forsythe (1974) Robust Tests for the Equality of Variances, Journal of the American Statistical Association, 69:346, 364-367, DOI: https://doi.org/10.1080/01621459.1974.10482955

   .. seealso::
       :py:func:`scipy.stats.levene <scipy:scipy.stats.levene>`.


.. py:function:: calc_levenetestSpatialPattern_spatial(data_input1: xarray.DataArray, data_input2: xarray.DataArray, dim: str = 'time', center: {'mean', 'median', 'trimmed'} = 'median', proportiontocut: float = 0.05) -> xarray.Dataset

   Perform Levene test for equal variances of two independent sptial samples along with other axis (i.e. 'time') of scores.

   The Levene test tests the null hypothesis that all input samples are from populations with equal variances.
   Levene's test is an alternative to Bartlett's test in the case where there are significant deviations from normality.

   Parameters
   ----------
   data_input1: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The first spatio-temporal data of xarray DataArray to be calculated.
   data_input2: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The second spatio-temporal data of xarray DataArray to be calculated.

   .. note::
       - The order of `data_input1` and `data_input2` has no effect on the calculation result.
       - The non-time dimensions of the two data sets must be exactly the same, and the dimensionality values must be arranged in the same order (ascending or descending).

   dim: :py:class:`str <str>`.
       Dimension(s) over which to apply the test. By default the test is applied over the `time` dimension.
   center: {'mean', 'median', 'trimmed'}, default `'median'`.
       Which function of the data to use in the test.

       .. note::

           Three variations of Levene’s test are possible. The possibilities and their recommended usages are:

           - median: Recommended for skewed (non-normal) distributions.
           - mean: Recommended for symmetric, moderate-tailed distributions.
           - trimmed: Recommended for heavy-tailed distributions.

           The test version using the mean was proposed in the original article of Levene (Levene, H., 1960) while the median and trimmed mean have been studied by Brown and Forsythe (Brown, M. B. and Forsythe, A. B., 1974), sometimes also referred to as Brown-Forsythe test.


   proportiontocut: :py:class:`float <float>`, default `0.05`.
       When center is `'trimmed'`, this gives the proportion of data points to cut from each end (See :py:func:`scipy.stats.trim_mean <scipy:scipy.stats.trim_mean>`).

   Returns
   -------
   - **statistic**, **pvalue**: :py:class:`xarray.Dataset<xarray.Dataset>`.

   Reference
   --------------
   - Levene, H. (1960). In Contributions to Probability and Statistics: Essays in Honor of Harold Hotelling, I. Olkin et al. eds., Stanford University Press, pp. 278-292.
   - Morton B. Brown & Alan B. Forsythe (1974) Robust Tests for the Equality of Variances, Journal of the American Statistical Association, 69:346, 364-367, DOI: https://doi.org/10.1080/01621459.1974.10482955

   .. seealso::
       :py:func:`scipy.stats.levene <scipy:scipy.stats.levene>`.


.. py:function:: calc_skewness_spatial(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time') -> xarray.Dataset | xarray.DataTree

   Calculate the skewness of the spatial field on the time axis and its significance test.

   The :math:`k` th statistical moment about the mean is given by

   .. math::
       m_k = \sum_{i=1}^{N} \frac{(x_i-\bar{x})^k}{N}

   where :math:`x_i` is the :math:`i` th observation, :math:`\bar{x}` the mean and :math:`N` the number of observations.

   One definition of the coefficient of skewness is

   .. math::
       a_3 = \frac{m_3}{(m_2)^{3/2}}

   Skewness is a measure of the asymmetry of a distribution and is zero for a normal distribution. If the longer wing of a distribution
   occurs for values of :math:`x` higher than the mean, that distribution is said to have positive skewness. If thelonger wing occurs for
   values of :math:`x` lower than the mean, the distribution is said to have negative skewness.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`
        The spatio-temporal data of xarray DataArray to be calculated.
   dim: :py:class:`str <str>`
       Dimension(s) over which to apply skewness. By default skewness is applied over the `time` dimension.

   Returns
   -------
   - **skewness**, **pvalue**: :py:class:`xarray.Dataset<xarray.Dataset>`.

   Reference
   --------------
   White, G. H. (1980). Skewness, Kurtosis and Extreme Values of
   Northern Hemisphere Geopotential Heights, Monthly Weather Review, 108(9), 1446-1455.
   Website: https://journals.ametsoc.org/view/journals/mwre/108/9/1520-0493_1980_108_1446_skaevo_2_0_co_2.xml

   .. seealso::
       :py:func:`scipy.stats.skew <scipy:scipy.stats.skew>`, :py:func:`scipy.stats.normaltest <scipy:scipy.stats.normaltest>`.


.. py:function:: calc_kurtosis_spatial(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time') -> xarray.DataArray | xarray.DataTree

   Calculate the kurtosis of the spatial field on the time axis and its significance test.

   The :math:`k` th statistical moment about the mean is given by

   .. math::
       m_k = \sum_{i=1}^{N} \frac{(x_i-\bar{x})^k}{N}

   where :math:`x_i` is the :math:`i` th observation, :math:`\bar{x}` the mean and :math:`N` the number of observations.

   The coefficient of kurtosis is defined by

   .. math::
       a_4 = \frac{m_4}{(m_2)^{2}}

   The kurtosis of a normal distribution is 3. If a distribution has a large central region which is flatter than a normal distribution
   with the same mean and variance, it has a kurtosis of less than 3. If the distribution has a central maximum more peaked and with
   longer wings than the equivalent normal distribution, its kurtosis is higher than 3 (Brooks and Carruthers, 1954).
   Extreme departures from the mean will cause very high values of kurtosis. Consequently, high kurtosis has been used as
   an indicator of bad data (Craddock and Flood, 1969). For the same reason, high values of kurtosis can be a result of one or two
   extreme events in a period of several years.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`
        The spatio-temporal data of xarray DataArray to be calculated.
   dim: :py:class:`str <str>`
       Dimension(s) over which to apply kurtosis. By default kurtosis is applied over the `time` dimension.

   Returns
   -------
   - kurtosis: :py:class:`xarray.DataArray<xarray.DataArray>`.

   Reference
   --------------
   White, G. H. (1980). Skewness, Kurtosis and Extreme Values of
   Northern Hemisphere Geopotential Heights, Monthly Weather Review, 108(9), 1446-1455.
   Website: https://journals.ametsoc.org/view/journals/mwre/108/9/1520-0493_1980_108_1446_skaevo_2_0_co_2.xml

   Køie, M., Brooks, C.E., & Carruthers, N. (1954). Handbook of Statistical Methods in Meteorology. Oikos, 4, 202.

   Craddock, J.M. and Flood, C.R. (1969), Eigenvectors for representing the 500 mb geopotential
   surface over the Northern Hemisphere. Q.J.R. Meteorol. Soc., 95: 576-593.
   doi: https://doi.org/10.1002/qj.49709540510

   .. seealso::
       :py:func:`scipy.stats.kurtosis <scipy:scipy.stats.kurtosis>`.


.. py:function:: calc_theilslopes_spatial(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', x=None, alpha: float = 0.95, method: {'joint', 'separate'} = 'separate', returns_type: {'dataset_returns', 'dataset_vars'} = 'dataset_returns') -> xarray.Dataset | xarray.DataTree

   Computes the Theil-Sen estimator.

   Theilslopes implements a method for robust linear regression. It computes the slope as the median of all slopes between paired values.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>` to be regression.
   dim: :py:class:`str <str>`, default `time`.
       Dimension(s) over which to apply linregress. By default linregress is applied over the `time` dimension.
   x: :py:class:`numpy.array <numpy.array>`
       Independent variable. If None, use `np.arange(len(data_input['time'].shape[0]))` instead.
   alpha: :py:class:`float <float>`, default 0.95.
       Confidence degree between 0 and 1. Default is 95% confidence. Note that alpha is symmetric around 0.5, i.e. both 0.1 and 0.9 are interpreted as "find the 90% confidence interval".
   method: {'joint', 'separate'}, default `'separate'`.
       Method to be used for computing estimate for intercept. Following methods are supported,

       - *joint*: Uses `np.median(y - slope * x)` as intercept.
       - *separate*: Uses `np.median(y) - slope * np.median(x)` as intercept.

   returns_type: :py:class:`str <str>`, default `'dataset_returns'`.
       Return data type.

   Returns
   -------
   result : ``TheilslopesResult`` Dataset
       The return Dataset have following data_var:

       **slope**: :py:class:`float <float>`
           Theil slope.
       **intercept**: :py:class:`float <float>`
           Intercept of the Theil line.
       **low_slope**: :py:class:`float <float>`
           Lower bound of the confidence interval on `slope`.
       **high_slope**: :py:class:`float <float>`
           Upper bound of the confidence interval on `slope`.

   .. seealso::
       :py:func:`scipy.stats.theilslopes <scipy:scipy.stats.theilslopes>`.


.. py:function:: calc_lead_lag_correlation_coefficients(pcs: dict, pairs: List[tuple], max_lag: int) -> xarray.Dataset

   Compute lead-lag correlation coefficients for specified pairs of indexes.

   This function calculates the cross-correlation between pairs of time series (e.g., MJO/BSISO principal components PC1 vs. PC2)
   to determine their lead-lag relationships. The correlation coefficients are computed for a range
   of lags, and the maximum correlation and corresponding lag are stored as attributes in the output
   dataset.

   Parameters
   ----------
   pcs : py:class:`dict <dict>`
       Dictionary mapping PC names (e.g., 'PC1', 'PC2') to :py:class:`xarray.DataArray<xarray.DataArray>` objects, where each
       DataArray represents a principal component time series with a ``'time'`` dimension.
   pairs : py:class:`list <list>` of tuples
       List of tuples, where each tuple contains (pair_name, pc_name1, pc_name2). For example,
       [('PC1_vs_PC2', 'PC1', 'PC2'), ('PC3_vs_PC4', 'PC3', 'PC4')].
   max_lag : py:class:`int <int>`
       Maximum lag (in time steps) to consider for the cross-correlation. The function computes
       correlations for lags from ``-max_lag`` to ``+max_lag``.

   Returns
   -------
   corr_da : :py:class:`xarray.Dataset<xarray.Dataset>`
       Dataset containing correlation coefficients for each pair, with a 'lag' dimension.
       Each variable (e.g., 'PC1_vs_PC2') has attributes ``'max_correlation' (float)`` and
       ``'lag_at_max_correlation' (int)`` indicating the maximum correlation and the lag at which it occurs.

   Notes
   -----

   - Positive lags indicate that the first PC leads the second; negative lags indicate the opposite.
   - The input PCs should have no missing values. Use ``fillna`` or ``interpolate_na`` if needed.
   - The correlation coefficients are normalized to range between :math:`-1` and :math:`1`.


.. py:function:: calc_timeseries_correlations(da: dict[str, xarray.DataArray] | list[xarray.DataArray], dim: str = 'time') -> xarray.DataArray

   Calculate the correlation matrix between multiple DataArray time series.

   This function calculates pairwise correlations between time series in the input DataArrays
   using the specified correlation method along the given dimension. The output is a symmetric
   correlation matrix stored as an xarray DataArray with dimensions (var1, var2).

   Parameters
   ----------
   da : :py:class:`dict[str, xarray.DataArray<xarray.DataArray>]` or :py:class:`list[xarray.DataArray<xarray.DataArray>]`.
       A dictionary with names as keys and DataArrays as values, or a list of DataArrays.
       Each DataArray must contain the specified dimension.
   dim : :py:class:`str <str>`, default: `time`.
       The dimension along which to compute correlations. All DataArrays must have this dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.
       A DataArray containing the correlation matrix with dimensions (var1, var2).
       Coordinates are set to the names of the input time series.

   Examples
   --------
   >>> time = pd.date_range('2020-01-01', '2020-12-31', freq='D')
   >>> da1 = xr.DataArray(np.random.randn(len(time)), dims='time', coords={'time': time}, name='series1')
   >>> da2 = xr.DataArray(da1 * 0.5 + np.random.randn(len(time)) * 0.5, dims='time', coords={'time': time}, name='series2')
   >>> data = {'series1': da1, 'series2': da2}
   >>> corr_matrix = calc_timeseries_correlations(data, method='pearson')
   >>> print(corr_matrix)


