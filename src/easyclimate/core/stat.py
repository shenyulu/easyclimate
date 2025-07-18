"""
Basic statistical analysis of weather and climate variables
"""

import numpy as np
import xarray as xr
import warnings

from scipy import signal, stats
from scipy.stats import pearsonr
from scipy.signal import correlate
from .utility import generate_datanode_dispatcher, find_dims_axis, validate_dataarrays

from xarray import DataTree
from .datanode import DataNode
from typing import List

__all__ = [
    "calc_linregress_spatial",
    "calc_detrend_spatial",
    "calc_corr_spatial",
    "calc_multiple_linear_regression_spatial",
    "calc_ttestSpatialPattern_spatial",
    "calc_levenetestSpatialPattern_spatial",
    "calc_skewness_spatial",
    "calc_kurtosis_spatial",
    "calc_theilslopes_spatial",
    "calc_lead_lag_correlation_coefficients",
]


@generate_datanode_dispatcher
def calc_linregress_spatial(
    data_input: xr.DataArray | xr.Dataset,
    dim: str = "time",
    x: np.array = None,
    alternative: str = "two-sided",
    returns_type: {"dataset_returns", "dataset_vars"} = "dataset_returns",
) -> xr.Dataset | DataTree:
    """
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
    """

    def _calc_linregress_spatial_scipy_linregress(data_input, dim, x, alternative):

        if data_input.chunks is not None:
            # Dask routine
            data_input = data_input.chunk({dim: -1})
        else:
            pass

        # y shape
        n = data_input[dim].shape[0]

        if x is None:
            # Regression parameter x
            x_data = np.arange(0, data_input[dim].shape[0])
            x = xr.DataArray(x_data, dims=dim, coords={dim: data_input[dim].data})

        x_shape = x.shape[0]
        if x_shape != n:
            raise ValueError(
                "`data_input` array size along dimension `dim` should be the same as the `x` array size, but data_input[dim]: "
                + str(n)
                + "; x: "
                + str(x_shape)
                + "."
            )

        if isinstance(x, np.ndarray):
            warnings.warn(
                f"Assuming that the coordinate value of '{dim}' in `data_input` and `x` is same. Ignoring."
            )
            x = xr.DataArray(x, dims=dim, coords={dim: data_input[dim].data})
        if isinstance(x, xr.DataArray):
            if x.dims[0] != data_input[dim].dims[0]:
                raise ValueError(
                    "The coordinate name of `data_input` array along dimension `dim` should be the same as the `x`."
                )
            if (x[dim].data == data_input[dim].data).all() == False:
                raise ValueError(
                    f"Coordinate value of '{dim}' in `data_input` and `x` is not same. If you are sure that the `x` dimension '{dim}' is the same as the value of the dimension in `data_input`, pass in the numpy array corresponding to `x`, e.g. `x.data`."
                )
        else:
            raise ValueError(
                "Unsupported input type. Expected numpy.ndarray or xarray.DataArray."
            )

        # `scipy.stats.linregress` only supports numpy arrays
        x = x.data

        # scipy function scipy.stats.linregress calculate regression parameter
        def linregress_scipy(data):
            LinregressResult = stats.linregress(x, data, alternative)
            slope = LinregressResult.slope
            intercept = LinregressResult.intercept
            rvalue = LinregressResult.rvalue
            pvalue = LinregressResult.pvalue
            stderr = LinregressResult.stderr
            intercept_stderr = LinregressResult.intercept_stderr
            return np.array(
                [slope, intercept, rvalue, pvalue, stderr, intercept_stderr]
            )

        # Use xarray apply_ufunc to create DataArray
        LinregressResult_dataarray = xr.apply_ufunc(
            linregress_scipy,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims=[["parameter"]],
            output_dtypes=["float64"],
            dask="parallelized",
            vectorize=True,
            dask_gufunc_kwargs={
                "output_sizes": {"parameter": 6},
                "allow_rechunk": True,
            },
        )

        # Transform DataArray to Dataset
        return xr.Dataset(
            data_vars={
                "slope": LinregressResult_dataarray[..., 0],
                "intercept": LinregressResult_dataarray[..., 1],
                "rvalue": LinregressResult_dataarray[..., 2],
                "pvalue": LinregressResult_dataarray[..., 3],
                "stderr": LinregressResult_dataarray[..., 4],
                "intercept_stderr": LinregressResult_dataarray[..., 5],
            }
        )

    return _calc_linregress_spatial_scipy_linregress(data_input, dim, x, alternative)


@generate_datanode_dispatcher
def calc_detrend_spatial(
    data_input: xr.DataArray | xr.Dataset, time_dim: str = "time"
) -> xr.DataArray | DataTree:
    """
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
    """

    # Because `scipy.signal.detrend` cannot detrend `np.nan`,
    # so we need to get the data mask first, then assign `np.nan` to 1 for calculation,
    # and then remove the mask area again.
    mask_bool = np.isnan(data_input).mean(dim=time_dim)
    mask_float = mask_bool + 0.0

    detrenddata_withoutmask = data_input.fillna(1).reduce(signal.detrend, dim=time_dim)
    result = detrenddata_withoutmask.where(mask_float < 0.5)
    return result


def calc_corr_spatial(
    data_input: xr.DataArray, x: xr.DataArray | np.ndarray, time_dim: str = "time"
) -> xr.Dataset:
    """
    Calculate Pearson correlation coefficients and corresponding p-values between spatial data
    and a time series using ``scipy.stats.pearsonr``.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        Input spatial data with dimensions ``(time, ...)``.

        NaN values are automatically skipped in calculations.
    x : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`numpy.ndarray<numpy.ndarray>`
        Time series data with dimension ``(time,)``. Must have the same length as data_input's time dimension.
        NaN values are automatically skipped in calculations.
    time_dim: :py:class:`str <str>`
        Dimension(s) over which to detrend. By default dimension is applied over the `time` dimension.

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

    Raises
    ------
    ValueError
        If the time dimension length of ``data_input`` and ``x`` don't match.

    Notes
    -----
    - The calculation automatically skips pairs where either value is NaN
    - Requires at least 2 valid observations for each grid point to compute correlation
    - Uses :py:func:`scipy.stats.pearsonr<scipy:scipy.stats.pearsonr>` for calculations
    - Maintains input coordinates for the output arrays

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
    """
    # 检查时间维度是否一致
    # Check whether the time dimensions are consistent
    if len(data_input[time_dim]) != len(x):
        raise ValueError(
            f"The time dimension is not consistent! data_input.time Length:{len(data_input[time_dim])}, the length of `x`: {len(x)}"
        )

    # 定义处理每个网格点的函数
    # Define a function that handles each grid point
    def _pearsonr_wrapper(x1, x2):
        # x1是空间数据的时间序列，x2是输入的时间序列
        # 跳过NaN值
        # x1 is the spatial data time series, x2 is the input time series
        # Skip NaN values
        mask = ~np.isnan(x1) & ~np.isnan(x2)
        if np.sum(mask) < 2:  # 至少需要2个有效点 You need at least 2 valid points
            return np.nan, np.nan
        return pearsonr(x1[mask], x2[mask])

    # 使用apply_ufunc计算
    # # Calculate using apply_ufunc
    result = xr.apply_ufunc(
        _pearsonr_wrapper,
        data_input,
        x,
        input_core_dims=[
            [time_dim],
            [time_dim],
        ],  # 每个输入的核心维度 Core dimension for each input
        output_core_dims=[[], []],  # 输出没有核心维度 The output has no core dimensions
        vectorize=True,  # 自动循环处理非核心维度 Automatic loop handling of non-core dimensions
        output_dtypes=[float, float],
        dask=(
            "parallelized" if data_input.chunks else False
        ),  # 支持dask数组 Support dask array
    )

    # 分离相关性和p值
    # Separate correlation and p-value
    corr = result[0].rename("correlation")
    pvalue = result[1].rename("pvalue")

    result_dataset = xr.Dataset()
    result_dataset["corr"] = corr
    result_dataset["pvalue"] = pvalue

    return result_dataset


def calc_multiple_linear_regression_spatial(
    y_data: xr.DataArray, x_datas: List[xr.DataArray], dim="time"
) -> xr.Dataset:
    """
    Apply multiple linear regression to dataset across spatial dimensions.

    .. math::

        y = a_1 x_1 + a_2 x_2 + \\cdots

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
    """

    def _multiple_linear_regression(y, *x_vars):
        """
        Perform multiple linear regression and compute relevant statistics.

        Parameters:
        -----------
        y : numpy.ndarray
            Dependent variable array with shape (time,)
        *x_vars : tuple of numpy.ndarray
            Independent variable arrays, each with shape (time,)

        Returns:
        --------
        tuple
            A tuple containing:

            - slopes (numpy.ndarray): Array of slope coefficients
            - intercept (float): Intercept value
            - r_squared (float): Coefficient of determination (R²)
            - slopes_p (numpy.ndarray): p-values for each slope coefficient
            - intercept_p (float): p-value for the intercept

        Notes:
        ------
        This function uses ordinary least squares (OLS) to estimate the regression coefficients.
        The R-squared value is calculated as 1 - (SS_res / SS_tot).
        Standard errors and p-values are computed using the residual mean squared error.
        """
        n_obs = len(y)
        n_vars = len(x_vars)

        # Stack the independent variables into a matrix (n_time, n_vars)
        x_matrix = np.column_stack(x_vars)

        # Add constant term (intercept)
        x_matrix = np.column_stack([x_matrix, np.ones(x_matrix.shape[0])])

        # Perform linear regression
        coefficients, _, rank, _ = np.linalg.lstsq(x_matrix, y, rcond=None)

        # Calculation of predicted values and residuals
        y_pred = np.dot(x_matrix, coefficients)
        residuals = y - y_pred

        # Compute the residual sum of squares
        ss_res = np.sum(residuals**2)

        # Calculate the total sum of squares
        y_mean = np.mean(y)
        ss_tot = np.sum((y - y_mean) ** 2)

        # Calculating R-squared
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else np.nan

        # Standard error of calculated coefficients
        if n_obs > rank:
            dof = n_obs - rank
            mse = ss_res / dof
        else:
            dof = 0
            mse = np.nan

        try:
            cov = np.linalg.inv(x_matrix.T @ x_matrix) * mse
            se = np.sqrt(np.diag(cov))

            # Calculate t-statistics and p-values
            t = coefficients / se
            p_values = 2 * stats.t.sf(np.abs(t), dof)
        except:
            se = np.full_like(coefficients, np.nan)
            p_values = np.full_like(coefficients, np.nan)

        # Separate slopes and intercepts
        slopes = coefficients[:-1]
        intercept = coefficients[-1]

        slopes_p = p_values[:-1]
        intercept_p = p_values[-1]

        return slopes, intercept, r_squared, slopes_p, intercept_p

    # Ensure all data has the same time dimension
    for x in x_datas:
        if not np.array_equal(y_data[dim].values, x[dim].values):
            warnings.warn("All variables must have the same time coordinates")

    # Applying regression in the spatial dimension using apply_ufunc
    result = xr.apply_ufunc(
        _multiple_linear_regression,
        y_data,
        *x_datas,
        input_core_dims=[[dim]] + [[dim]] * len(x_datas),
        output_core_dims=[
            ["coef"],
            [],
            [],
            ["coef"],
            [],
        ],  # Slope, Intercept, R-squared, Slope p-value, Intercept p-value
        exclude_dims=set((dim,)),
        vectorize=True,
    )

    # Collation results
    slopes = result[0].rename("slopes")
    intercept = result[1].rename("intercept")
    r_squared = result[2].rename("r_squared")
    slopes_p = result[3].rename("slopes_p")
    intercept_p = result[4].rename("intercept_p")

    # Creating a Dataset
    ds = xr.Dataset(
        {
            "slopes": slopes.assign_coords(coef=np.arange(len(x_datas))),
            "intercept": intercept,
            "r_squared": r_squared,
            "slopes_p": slopes_p.assign_coords(coef=np.arange(len(x_datas))),
            "intercept_p": intercept_p,
        }
    )

    return ds


def calc_ttestSpatialPattern_spatial(
    data_input1: xr.DataArray, data_input2: xr.DataArray, dim: str = "time"
) -> xr.Dataset:
    """
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
    """

    new_dims = tuple(dim for dim in data_input1.dims if dim != dim)
    validate_dataarrays([data_input1, data_input2], dims=new_dims, time_dims=dim)

    # scipy function scipy.stats.ttest_ind calculate the T-test for the means of two independent samples of scores.
    def _ttest_ind_scipy(data1, data2):
        statistic, pvalue = stats.ttest_ind(data1, data2)
        return np.array([statistic, pvalue])

    # Use xarray apply_ufunc to create DataArray
    ttest_ind_dataarray = xr.apply_ufunc(
        _ttest_ind_scipy,
        data_input1,
        data_input2,
        input_core_dims=[[dim], [dim]],
        output_core_dims=[["parameter"]],
        output_dtypes=["float64"],
        dask="parallelized",
        vectorize=True,
        dask_gufunc_kwargs={"output_sizes": {"parameter": 2}, "allow_rechunk": True},
        exclude_dims=set((dim,)),  # allow change size
    )

    return xr.Dataset(
        data_vars={
            "statistic": ttest_ind_dataarray[..., 0],
            "pvalue": ttest_ind_dataarray[..., 1],
        }
    )


def calc_levenetestSpatialPattern_spatial(
    data_input1: xr.DataArray,
    data_input2: xr.DataArray,
    dim: str = "time",
    center: {"mean", "median", "trimmed"} = "median",
    proportiontocut: float = 0.05,
) -> xr.Dataset:
    """
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
    """

    new_dims = tuple(dim for dim in data_input1.dims if dim != dim)
    validate_dataarrays([data_input1, data_input2], dims=new_dims, time_dims=dim)

    # scipy function scipy.stats.levene calculate the F-test for the means of two independent samples of scores.
    def _levenetest_ind_scipy(data1, data2):
        statistic, pvalue = stats.levene(
            data1, data2, center=center, proportiontocut=proportiontocut
        )
        return np.array([statistic, pvalue])

    # Use xarray apply_ufunc to create DataArray
    levenetest_ind_dataarray = xr.apply_ufunc(
        _levenetest_ind_scipy,
        data_input1,
        data_input2,
        input_core_dims=[[dim], [dim]],
        output_core_dims=[["parameter"]],
        output_dtypes=["float64"],
        dask="parallelized",
        vectorize=True,
        dask_gufunc_kwargs={"output_sizes": {"parameter": 2}, "allow_rechunk": True},
        exclude_dims=set((dim,)),  # allow change size
    )

    return xr.Dataset(
        data_vars={
            "statistic": levenetest_ind_dataarray[..., 0],
            "pvalue": levenetest_ind_dataarray[..., 1],
        }
    )


def calc_levenetestSpatialPattern_spatial(
    data_input1: xr.DataArray,
    data_input2: xr.DataArray,
    dim: str = "time",
    center: {"mean", "median", "trimmed"} = "median",
    proportiontocut: float = 0.05,
) -> xr.Dataset:
    """
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
    """

    new_dims = tuple(dim for dim in data_input1.dims if dim != dim)
    validate_dataarrays([data_input1, data_input2], dims=new_dims, time_dims=dim)

    # scipy function scipy.stats.levene calculate the F-test for the means of two independent samples of scores.
    def _levenetest_ind_scipy(data1, data2):
        statistic, pvalue = stats.levene(
            data1, data2, center=center, proportiontocut=proportiontocut
        )
        return np.array([statistic, pvalue])

    # Use xarray apply_ufunc to create DataArray
    levenetest_ind_dataarray = xr.apply_ufunc(
        _levenetest_ind_scipy,
        data_input1,
        data_input2,
        input_core_dims=[[dim], [dim]],
        output_core_dims=[["parameter"]],
        output_dtypes=["float64"],
        dask="parallelized",
        vectorize=True,
        dask_gufunc_kwargs={"output_sizes": {"parameter": 2}, "allow_rechunk": True},
        exclude_dims=set((dim,)),  # allow change size
    )

    return xr.Dataset(
        data_vars={
            "statistic": levenetest_ind_dataarray[..., 0],
            "pvalue": levenetest_ind_dataarray[..., 1],
        }
    )


@generate_datanode_dispatcher
def calc_skewness_spatial(
    data_input: xr.DataArray | xr.Dataset,
    dim: str = "time",
) -> xr.Dataset | DataTree:
    """
    Calculate the skewness of the spatial field on the time axis and its significance test.

    The :math:`k` th statistical moment about the mean is given by

    .. math::
        m_k = \\sum_{i=1}^{N} \\frac{(x_i-\\bar{x})^k}{N}

    where :math:`x_i` is the :math:`i` th observation, :math:`\\bar{x}` the mean and :math:`N` the number of observations.

    One definition of the coefficient of skewness is

    .. math::
        a_3 = \\frac{m_3}{(m_2)^{3/2}}

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
    """
    # Find the index of `dim` in the xarray DataArray for `time`.
    time_dim_index = find_dims_axis(data_input, dim=dim)

    # Calculate skewness
    def _calc_skew_core(data_input, data_all):
        time_length = data_input[dim].shape[0]
        m_3 = (((data_input - data_all) ** 3) / time_length).sum(dim=dim)
        m_2 = (((data_input - data_all) ** 2) / time_length).sum(dim=dim)
        return m_3 / (m_2 ** (3 / 2))

    data_all = data_input.mean(dim=dim)
    skewness = _calc_skew_core(data_input, data_all)

    # Significance test
    k2, p_numpy = stats.normaltest(
        data_input, axis=time_dim_index, nan_policy="propagate"
    )
    format_coordniate = data_all
    p = format_coordniate.copy(data=p_numpy, deep=True)

    # Merge multiple :py:class:`xarray.DataArray<xarray.DataArray>`s into one `xarray.Dataset`.
    dateset = xr.Dataset(data_vars={"skewness": skewness, "pvalue": p})

    return dateset


@generate_datanode_dispatcher
def calc_kurtosis_spatial(
    data_input: xr.DataArray | xr.Dataset,
    dim: str = "time",
) -> xr.DataArray | DataTree:
    """
    Calculate the kurtosis of the spatial field on the time axis and its significance test.

    The :math:`k` th statistical moment about the mean is given by

    .. math::
        m_k = \\sum_{i=1}^{N} \\frac{(x_i-\\bar{x})^k}{N}

    where :math:`x_i` is the :math:`i` th observation, :math:`\\bar{x}` the mean and :math:`N` the number of observations.

    The coefficient of kurtosis is defined by

    .. math::
        a_4 = \\frac{m_4}{(m_2)^{2}}

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
    """

    # Calculate kurtosis
    def _calc_kurt_core(data_input, data_all):
        time_length = data_input[dim].shape[0]
        m_4 = (((data_input - data_all) ** 4) / time_length).sum(dim=dim)
        m_2 = (((data_input - data_all) ** 2) / time_length).sum(dim=dim)
        return m_4 / (m_2**2)

    data_all = data_input.mean(dim=dim)
    kurtosis = _calc_kurt_core(data_input, data_all)

    return kurtosis


@generate_datanode_dispatcher
def calc_theilslopes_spatial(
    data_input: xr.DataArray | xr.Dataset,
    dim: str = "time",
    x=None,
    alpha: float = 0.95,
    method: {"joint", "separate"} = "separate",
    returns_type: {"dataset_returns", "dataset_vars"} = "dataset_returns",
) -> xr.Dataset | DataTree:
    """
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
    """

    def _calc_theilslopes_spatial_scipy_theilslopes(data_input, dim, x, alpha, method):

        if data_input.chunks is not None:
            # Dask routine
            data_input = data_input.chunk({dim: -1})
        else:
            pass

        # y shape
        n = data_input[dim].shape[0]

        if x is None:
            # Regression parameter x
            x_data = np.arange(0, data_input[dim].shape[0])
            x = xr.DataArray(x_data, dims=dim, coords={dim: data_input[dim].data})

        x_shape = x.shape[0]
        if x_shape != n:
            raise ValueError(
                "`data_input` array size along dimension `dim` should be the same as the `x` array size, but data_input[dim]: "
                + str(n)
                + "; x: "
                + str(x_shape)
                + "."
            )

        if isinstance(x, np.ndarray):
            warnings.warn(
                f"Assuming that the coordinate value of '{dim}' in `data_input` and `x` is same. Ignoring."
            )
            x = xr.DataArray(x, dims=dim, coords={dim: data_input[dim].data})
        if isinstance(x, xr.DataArray):
            if x.dims[0] != data_input[dim].dims[0]:
                raise ValueError(
                    "The coordinate name of `data_input` array along dimension `dim` should be the same as the `x`."
                )
            if (x[dim].data == data_input[dim].data).all() == False:
                raise ValueError(
                    f"Coordinate value of '{dim}' in `data_input` and `x` is not same. If you are sure that the `x` dimension '{dim}' is the same as the value of the dimension in `data_input`, pass in the numpy array corresponding to `x`, e.g. `x.data`."
                )
        else:
            raise ValueError(
                "Unsupported input type. Expected numpy.ndarray or xarray.DataArray."
            )

        # `scipy.stats.linregress` only supports numpy arrays
        x = x.data

        # scipy function scipy.stats.linregress calculate regression parameter
        def theilslopes_scipy(data):
            TheilslopesResult = stats.theilslopes(data, x, alpha=alpha, method=method)
            slope = TheilslopesResult[0]
            intercept = TheilslopesResult[1]
            low_slope = TheilslopesResult[2]
            high_slope = TheilslopesResult[3]
            return np.array([slope, intercept, low_slope, high_slope])

        # Use xarray apply_ufunc to create DataArray
        TheilslopesResult_dataarray = xr.apply_ufunc(
            theilslopes_scipy,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims=[["parameter"]],
            output_dtypes=["float64"],
            dask="parallelized",
            vectorize=True,
            dask_gufunc_kwargs={
                "output_sizes": {"parameter": 4},
                "allow_rechunk": True,
            },
        )

        # Transform DataArray to Dataset
        return xr.Dataset(
            data_vars={
                "slope": TheilslopesResult_dataarray[..., 0],
                "intercept": TheilslopesResult_dataarray[..., 1],
                "low_slope": TheilslopesResult_dataarray[..., 2],
                "high_slope": TheilslopesResult_dataarray[..., 3],
            }
        )

    return _calc_theilslopes_spatial_scipy_theilslopes(
        data_input, dim, x, alpha, method
    )


def calc_lead_lag_correlation_coefficients(
    pcs: dict, pairs: List[tuple], max_lag: int
) -> xr.Dataset:
    """
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
    """

    def cross_correlation(da1, da2, max_lag):
        """
        Compute cross-correlation between two xarray DataArrays for a range of lags.

        Parameters
        ----------
        da1, da2 : xarray.DataArray
            Time series to correlate, with a 'time' dimension.
        max_lag : int
            Maximum lag to consider (in time steps).

        Returns
        -------
        lags : numpy.ndarray
            Array of lags from -max_lag to +max_lag.
        corr : numpy.ndarray
            Array of correlation coefficients corresponding to the lags.
        """
        x = da1.values - da1.mean().values
        y = da2.values - da2.mean().values
        corr = correlate(x, y, mode="full", method="auto")
        corr = corr / (np.std(x) * np.std(y) * len(x))
        lags = np.arange(-len(x) + 1, len(x))
        mask = (lags >= -max_lag) & (lags <= max_lag)
        return lags[mask], corr[mask]

    # Compute correlations for each pair
    results = {}
    for pair_name, pc_name1, pc_name2 in pairs:
        lags, corr = cross_correlation(pcs[pc_name1], pcs[pc_name2], max_lag)
        results[pair_name] = (lags, corr)

    # Create xarray Dataset
    corr_da = xr.Dataset(
        {
            pair_name: xr.DataArray(
                corr, dims=["lag"], coords={"lag": lags}, name=pair_name
            )
            for pair_name, (lags, corr) in results.items()
        }
    )

    # Add maximum correlation and lag as attributes
    for pair_name, (lags, corr) in results.items():
        max_corr_idx = np.argmax(np.abs(corr))
        max_lag = lags[max_corr_idx]
        max_corr = corr[max_corr_idx]
        corr_da[pair_name].attrs["max_correlation"] = float(max_corr)
        corr_da[pair_name].attrs["lag_at_max_correlation"] = int(max_lag)

    return corr_da
