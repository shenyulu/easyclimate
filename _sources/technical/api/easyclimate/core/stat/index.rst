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
   easyclimate.core.stat.calc_ttestSpatialPattern_spatial
   easyclimate.core.stat.calc_levenetestSpatialPattern_spatial
   easyclimate.core.stat.calc_levenetestSpatialPattern_spatial
   easyclimate.core.stat.calc_skewness_spatial
   easyclimate.core.stat.calc_kurtosis_spatial
   easyclimate.core.stat.calc_theilslopes_spatial


Module Contents
---------------

.. py:function:: calc_linregress_spatial(data_input: xarray.DataArray | xarray.Dataset, dim: str = 'time', x: numpy.array = None, alternative: str = 'two-sided', returns_type: {'dataset_returns', 'dataset_vars'} = 'dataset_returns') -> xarray.Dataset | xarray.DataTree

   Calculate a linear least-squares regression for spatial data of time.

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
   dim: :py:class:`str <str>`
       Dimension(s) over which to detrend. By default dimension is applied over the `time` dimension.

   Returns
   -------
   - :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::
       :py:func:`scipy.signal.detrend <scipy:scipy.signal.detrend>`.


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


