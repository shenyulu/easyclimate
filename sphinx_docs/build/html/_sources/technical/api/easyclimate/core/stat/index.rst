:py:mod:`easyclimate.core.stat`
===============================

.. py:module:: easyclimate.core.stat


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.core.stat.calc_linregress_spatial
   easyclimate.core.stat.calc_detrend_data
   easyclimate.core.stat.calc_ttestSpatialPattern_spatial
   easyclimate.core.stat.calc_skewness_spatial
   easyclimate.core.stat.calc_kurtosis_spatial



.. py:function:: calc_linregress_spatial(data_input, dim='time', x=None, alternative='two-sided', returns_type='dataset_returns', engine='scipy_linregress')

   Calculate a linear least-squares regression for spatial data of time.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>` to be regression.
   dim : str
       Dimension(s) over which to apply linregress. By default linregress is applied over the `time` dimension.
   x : numpy.array
   returns_type: str

   Returns
   -------
   result : ``LinregressResult`` Dataset
       The return Dataset have following data_var:

       slope : float
           Slope of the regression line.
       intercept : float
           Intercept of the regression line.
       rvalue : float
           The Pearson correlation coefficient. The square of ``rvalue``
           is equal to the coefficient of determination.
       pvalue : float
           The p-value for a hypothesis test whose null hypothesis is
           that the slope is zero, using Wald Test with t-distribution of
           the test statistic. See `alternative` above for alternative
           hypotheses.
       stderr : float
           Standard error of the estimated slope (gradient), under the
           assumption of residual normality.
       intercept_stderr : float
           Standard error of the estimated intercept, under the assumption
           of residual normality.

   .. seealso::
       :py:func:`scipy.stats.linregress <scipy:scipy.stats.linregress>`.


.. py:function:: calc_detrend_data(data_input, time_dim='time')

   Remove linear trend along axis from data.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        The spatio-temporal data of :py:class:`xarray.DataArray<xarray.DataArray>` to be detrended.
   dim : `str`
       Dimension(s) over which to detrend. By default dimension is applied over the `time` dimension.

   Returns
   -------
   - :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::
       :py:func:`scipy.signal.detrend <scipy:scipy.signal.detrend>`.


.. py:function:: calc_ttestSpatialPattern_spatial(data_input1, data_input2, dim='time')

   Calculate the T-test for the means of two independent sptial samples along with other axis (i.e. 'time') of scores.

   Parameters
   ----------
   data_input1 : :py:class:`xarray.DataArray<xarray.DataArray>`
        The first spatio-temporal data of xarray DataArray to be calculated.
   data_input2 : :py:class:`xarray.DataArray<xarray.DataArray>`
        The second spatio-temporal data of xarray DataArray to be calculated.

   .. note::
       The order of `data_input1` and `data_input2` has no effect on the calculation result.

   dim : `str`
       Dimension(s) over which to apply skewness. By default skewness is applied over the `time` dimension.

   Returns
   -------
   - statistic, p: :py:class:`xarray.Dataset<xarray.Dataset>`.

   .. seealso::
       :py:func:`scipy.stats.ttest_ind <scipy:scipy.stats.ttest_ind>`.


.. py:function:: calc_skewness_spatial(data_input, dim='time')

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
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        The spatio-temporal data of xarray DataArray to be calculated.
   dim : str
       Dimension(s) over which to apply skewness. By default skewness is applied over the `time` dimension.

   Returns
   -------
   - skewness, pvalue: :py:class:`xarray.Dataset<xarray.Dataset>`.

   Reference
   --------------
   White, G. H. (1980). Skewness, Kurtosis and Extreme Values of 
   Northern Hemisphere Geopotential Heights, Monthly Weather Review, 108(9), 1446-1455. 
   Website: https://journals.ametsoc.org/view/journals/mwre/108/9/1520-0493_1980_108_1446_skaevo_2_0_co_2.xml

   .. seealso::
       :py:func:`scipy.stats.skew <scipy:scipy.stats.skew>`, :py:func:`scipy.stats.normaltest <scipy:scipy.stats.normaltest>`.


.. py:function:: calc_kurtosis_spatial(data_input, dim='time')

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
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        The spatio-temporal data of xarray DataArray to be calculated.
   dim : str
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


