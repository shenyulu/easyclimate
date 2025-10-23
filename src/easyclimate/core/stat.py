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
from .normalized import normalize_zscore
from rich.progress import Progress, BarColumn, TimeRemainingColumn
from xarray import DataTree
from .datanode import DataNode
from typing import List, Literal

__all__ = [
    "calc_linregress_spatial",
    "calc_detrend_spatial",
    "calc_corr_spatial",
    "calc_leadlag_corr_spatial",
    "calc_multiple_linear_regression_spatial",
    "calc_ttestSpatialPattern_spatial",
    "calc_windmask_ttestSpatialPattern_spatial",
    "calc_levenetestSpatialPattern_spatial",
    "calc_skewness_spatial",
    "calc_kurtosis_spatial",
    "calc_theilslopes_spatial",
    "calc_lead_lag_correlation_coefficients",
    "calc_timeseries_correlations",
    "calc_non_centered_corr",
    "calc_pattern_corr",
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

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_basic_statistical_analysis.py
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
    data_input: xr.DataArray,
    x: xr.DataArray | np.ndarray,
    time_dim: str = "time",
    method: Literal["scipy", "xarray"] = "xarray",
) -> xr.Dataset:
    """
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
    reg_coeff, corr & pvalue (:py:class:`xarray.Dataset<xarray.Dataset>`)

    reg_coeff: :py:class:`xarray.DataArray<xarray.DataArray>`
        Regression coefficient, in units of ``data_input`` per standard deviation of the index.

    corr : :py:class:`xarray.DataArray<xarray.DataArray>`
        Pearson correlation coefficients with dimensions.
        Values range from -1 to 1 where:

        - 1: perfect positive correlation
        - -1: perfect negative correlation
        - 0: no correlation

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

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_corr_reg.py
    """
    # Check whether the time dimensions are consistent
    if len(data_input[time_dim]) != len(x):
        raise ValueError(
            f"The time dimension is not consistent! data_input.time Length:{len(data_input[time_dim])}, the length of `x`: {len(x)}"
        )

    # Covariance (regression coefficient)
    reg_coeff = xr.cov(data_input, normalize_zscore(x, dim=time_dim), dim=time_dim)

    if method == "scipy":
        # Define a function that handles each grid point
        def _pearsonr_wrapper(x1, x2):
            # x1 is the spatial data time series, x2 is the input time series
            # Skip NaN values
            mask = ~np.isnan(x1) & ~np.isnan(x2)
            if np.sum(mask) < 2:  # You need at least 2 valid points
                return np.nan, np.nan
            return pearsonr(x1[mask], x2[mask])

        # # Calculate using apply_ufunc
        result = xr.apply_ufunc(
            _pearsonr_wrapper,
            data_input,
            x,
            input_core_dims=[
                [time_dim],
                [time_dim],
            ],  # Core dimension for each input
            output_core_dims=[[], []],  # The output has no core dimensions
            vectorize=True,  # Automatic loop handling of non-core dimensions
            output_dtypes=[float, float],
            dask=("parallelized" if data_input.chunks else False),  # Support dask array
        )

        # Separate correlation and p-value
        corr = result[0].rename("correlation")
        pvalue = result[1].rename("pvalue")

    elif method == "xarray":
        # Calculate the correlation coefficient
        corr = xr.corr(data_input, x, dim=time_dim)

        # Calculate the t-statistic
        N = len(data_input[time_dim])
        t_stats = (corr * np.sqrt(N - 2)) / np.sqrt(1 - corr**2)
        # Degree of freedom
        df = N - 2

        def _t_test_func(t):
            if np.isnan(t).any():
                return np.nan
            else:
                return 2 * (1 - stats.t.cdf(np.abs(t), df))

        pvalue = xr.apply_ufunc(
            _t_test_func,
            t_stats,
            vectorize=True,
            dask="parallelized",
            dask_gufunc_kwargs={
                "allow_rechunk": True,
            },
        )

    else:
        raise ValueError("The parameter of `method` should be `xarray` or `scipy.`")

    result = xr.Dataset()
    result["reg_coeff"] = reg_coeff
    result["corr"] = corr
    result["pvalue"] = pvalue
    return result


def calc_leadlag_corr_spatial(
    data_input: xr.DataArray,
    x: xr.DataArray | np.ndarray,
    leadlag_array: np.array | List[int],
    time_dim: str = "time",
    method: Literal["scipy", "xarray"] = "xarray",
):
    """
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
    """
    leadlag_list = []

    with Progress(
        "[progress.description]{task.description}",
        BarColumn(),
        "[progress.percentage]{task.percentage:>3.0f}%",
        TimeRemainingColumn(),
    ) as progress:
        task = progress.add_task(
            "[cyan]Calculating correlations...", total=len(leadlag_array)
        )

        for leadlag_item in leadlag_array:
            tmp = calc_corr_spatial(
                data_input,
                x=x.shift({time_dim: leadlag_item}),
                time_dim=time_dim,
                method=method,
            )
            tmp = tmp.assign_coords({"leadlag": leadlag_item}).expand_dims("leadlag")
            leadlag_list.append(tmp)
            progress.update(
                task, advance=1, description=f"[cyan]Processing leadlag: {leadlag_item}"
            )

        result = xr.concat(leadlag_list, dim="leadlag")
        # Make sure the progress bar refreshes to 100%.
        progress.update(task, completed=len(leadlag_array), refresh=True)

    return result


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

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_multi_linear_reg.py
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
    data_input1: xr.DataArray,
    data_input2: xr.DataArray,
    dim: str = "time",
    equal_var: bool = True,
    alternative: Literal["two-sided", "less", "greater"] = "two-sided",
    method: Literal["scipy", "xarray"] = "xarray",
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

    equal_var: :py:class:`bool <bool>`
        If True (default), perform a standard independent 2 sample test that assumes equal population variances (see https://en.wikipedia.org/wiki/T-test#Independent_two-sample_t-test).
        If False, perform Welch’s t-test, which does not assume equal population variance (see https://en.wikipedia.org/wiki/Welch%27s_t-test).

    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis.
        The following options are available (default is 'two-sided'):

        - 'two-sided': the means of the distributions underlying the samples are unequal.
        - 'less': the mean of the distribution underlying the first sample is less than the mean of the distribution underlying the second sample.
        - 'greater': the mean of the distribution underlying the first sample is greater than the mean of the distribution underlying the second sample.

    method : {'scipy', 'xarray'}, optional
        Method used to compute correlations:

        - 'scipy': Uses :py:func:`scipy.stats.ttest_ind<scipy:scipy.stats.ttest_ind>` for direct calculation
        - 'xarray': Uses xarray's built-in method to calculate (faster)

        Default is 'xarray'.

    Returns
    -------
    - **statistic**, **pvalue**: :py:class:`xarray.Dataset<xarray.Dataset>`.

    .. seealso::
        :py:func:`scipy.stats.ttest_ind <scipy:scipy.stats.ttest_ind>`.
    """

    new_dims = tuple(dim for dim in data_input1.dims if dim != dim)
    validate_dataarrays([data_input1, data_input2], dims=new_dims, time_dims=dim)

    if method == "scipy":
        # scipy function scipy.stats.ttest_ind calculate the T-test for the means of two independent samples of scores.
        def _ttest_ind_scipy(data1, data2):
            # Check if there are any NaN values in data1 or data2.
            if np.isnan(data1).any() or np.isnan(data2).any():
                return np.array([np.nan, np.nan])
            statistic, pvalue = stats.ttest_ind(
                data1, data2, equal_var=equal_var, alternative=alternative
            )
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
            dask_gufunc_kwargs={
                "output_sizes": {"parameter": 2},
                "allow_rechunk": True,
            },
            exclude_dims=set((dim,)),  # allow change size
        )

        return xr.Dataset(
            data_vars={
                "statistic": ttest_ind_dataarray[..., 0],
                "pvalue": ttest_ind_dataarray[..., 1],
            }
        )

    elif method == "xarray":
        # Compute mask for points with no NaNs along the test dimension (nan_policy='propagate')
        has_nan1 = data_input1.isnull().any(dim=dim)
        has_nan2 = data_input2.isnull().any(dim=dim)
        mask = ~(has_nan1 | has_nan2)

        # Compute means and sample variances along the dimension (skipna=True, but masked points will be NaN)
        mean1 = data_input1.mean(dim=dim).where(mask)
        mean2 = data_input2.mean(dim=dim).where(mask)
        var1 = data_input1.var(dim=dim, ddof=1).where(mask)
        var2 = data_input2.var(dim=dim, ddof=1).where(mask)

        # Sample sizes (scalars; for valid points, full lengths apply as no NaNs)
        n1 = data_input1.sizes[dim]
        n2 = data_input2.sizes[dim]

        # Compute t-statistic and df based on equal_var
        diff = mean1 - mean2
        if equal_var:
            # Pooled variance for equal_var=True
            sp2 = ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2)
            denom = np.sqrt(sp2 * (1 / n1 + 1 / n2))
            df = n1 + n2 - 2
        else:
            # Welch's t-test (equal_var=False)
            denom = np.sqrt(var1 / n1 + var2 / n2)
            var_term = var1 / n1 + var2 / n2
            df = var_term**2 / (
                (var1 / n1) ** 2 / (n1 - 1) + (var2 / n2) ** 2 / (n2 - 1)
            )

        statistic = diff / denom

        # Compute p-value based on alternative
        def _pvalue_func(x, df_val):
            t_abs = np.abs(x)
            if alternative == "two-sided":
                return 2 * stats.t.sf(t_abs, df_val)
            elif alternative == "less":
                return stats.t.cdf(x, df_val)
            elif alternative == "greater":
                return stats.t.sf(x, df_val)
            else:
                raise ValueError(
                    f"Invalid alternative: {alternative}. Must be 'two-sided', 'less', or 'greater'."
                )

        pvalue = xr.apply_ufunc(
            _pvalue_func,
            statistic,
            df,
            input_core_dims=[[], []],
            output_core_dims=[[]],
            output_dtypes=["float64"],
            dask="parallelized",
        )

        return xr.Dataset(
            data_vars={
                "statistic": statistic,
                "pvalue": pvalue,
            }
        )


def calc_windmask_ttestSpatialPattern_spatial(
    data_input1: xr.Dataset,
    data_input2: xr.Dataset,
    dim: str = "time",
    u_dim: str = "u",
    v_dim: str = "v",
    mask_method: Literal["or", "and"] = "or",
    thresh: float = 0.05,
    equal_var: bool = True,
    alternative: Literal["two-sided", "less", "greater"] = "two-sided",
    method: Literal["scipy", "xarray"] = "xarray",
):
    """
    Generate a significance mask for T-tests on the means of two independent spatial zonal (u) and meridional (v) wind samples,
    aggregated over the specified dimension (default 'time').

    Parameters
    ----------
    data_input1 : :py:class:`xarray.Dataset`
         The first spatio-temporal data of xarray Dataset to be calculated. It is necessary to include the zonal wind component (u_dim) and the meridional wind component (v_dim).
    data_input2 : :py:class:`xarray.Dataset`
         The second spatio-temporal data of xarray Dataset to be calculated. It is necessary to include the zonal wind component (u_dim) and the meridional wind component (v_dim).

    .. note::
        - The order of `data_input1` and `data_input2` has no effect on the calculation result.
        - The non-time dimensions of the two data sets must be exactly the same, and the dimensionality values must be arranged in the same order (ascending or descending).

    dim : :py:class:`str`, default: `time`
        Dimension(s) over which to apply the test. By default the test is applied over the `time` dimension.
    u_dim : :py:class:`str`, default: `u`
        Variable name for the u velocity (zonal wind, in x direction).
    v_dim : :py:class:`str`, default: `v`
        Variable name for the v velocity (meridional wind, in y direction).
    mask_method : Literal["or", "and"], default: "or"
        Method to combine the significance masks for u and v components:

        - "or": A grid point is considered significant if either the u or v component is significant (p <= thresh).
        - "and": A grid point is considered significant if both the u and v components are significant (p <= thresh).

    thresh : :py:class:`float`, default: 0.05
        The significance level (alpha) for the p-value threshold used to create the mask.
    equal_var : :py:class:`bool`, default: True
        If True (default), perform a standard independent 2 sample test that assumes equal population variances (see https://en.wikipedia.org/wiki/T-test#Independent_two-sample_t-test).
        If False, perform Welch’s t-test, which does not assume equal population variance (see https://en.wikipedia.org/wiki/Welch%27s_t-test).

    alternative : {'two-sided', 'less', 'greater'}, optional=
        Defines the alternative hypothesis.
        The following options are available (default is 'two-sided'):

        - 'two-sided': the means of the distributions underlying the samples are unequal.
        - 'less': the mean of the distribution underlying the first sample is less than the mean of the distribution underlying the second sample.
        - 'greater': the mean of the distribution underlying the first sample is greater than the mean of the distribution underlying the second sample.

    method : {'scipy', 'xarray'}, optional
        Method used to compute t-tests:

        - 'scipy': Uses :py:func:`scipy.stats.ttest_ind<scipy:scipy.stats.ttest_ind>` for direct calculation
        - 'xarray': Uses xarray's built-in method to calculate (faster)

        Default is 'xarray'.

    Returns
    -------
    masked_pvalue : :py:class:`xarray.DataArray`
        A boolean mask indicating significant regions (True where p <= thresh, combined via mask_method for u and v).

    .. seealso::
        :py:func:`scipy.stats.ttest_ind <scipy:scipy.stats.ttest_ind>`.
    """

    def _get_wind_mask(u_pvalue, v_pvalue, p_thresh=thresh):
        if mask_method == "or":
            return (u_pvalue <= p_thresh) | (v_pvalue <= p_thresh)
        elif mask_method == "and":
            return (u_pvalue <= p_thresh) & (v_pvalue <= p_thresh)
        else:
            raise ValueError(
                f"Invalid mask_method: {mask_method}. Must be 'or' or 'and'."
            )

    uwnd_pvalue = calc_ttestSpatialPattern_spatial(
        data_input1[u_dim],
        data_input2[u_dim],
        dim=dim,
        equal_var=equal_var,
        alternative=alternative,
        method=method,
    ).pvalue

    vwnd_pvalue = calc_ttestSpatialPattern_spatial(
        data_input1[v_dim],
        data_input2[v_dim],
        dim=dim,
        equal_var=equal_var,
        alternative=alternative,
        method=method,
    ).pvalue

    uvwnd_pvalue = _get_wind_mask(uwnd_pvalue, vwnd_pvalue)
    return uvwnd_pvalue


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

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_basic_statistical_analysis.py
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

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_basic_statistical_analysis.py
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


def calc_timeseries_correlations(
    da: dict[str, xr.DataArray] | list[xr.DataArray],
    dim: str = "time",
) -> xr.DataArray:
    """
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
    """
    # Handle input format
    if isinstance(da, dict):
        names = list(da.keys())
        arrays = list(da.values())
    else:
        names = [f"var_{i}" for i in range(len(da))]
        arrays = da

    # Validate inputs
    if not arrays:
        raise ValueError("data_arrays cannot be empty")
    if not all(isinstance(da, xr.DataArray) for da in arrays):
        raise TypeError("All inputs must be xarray.DataArray")
    if not all(dim in da.dims for da in arrays):
        raise ValueError(f"All DataArrays must contain the '{dim}' dimension")

    # Initialize correlation matrix
    n = len(arrays)
    corr_matrix = np.zeros((n, n))

    # Compute pairwise correlations
    for i in range(n):
        for j in range(i, n):  # Compute upper triangle (including diagonal)
            corr = xr.corr(arrays[i], arrays[j], dim=dim)
            corr_matrix[i, j] = corr
            if i != j:  # Fill lower triangle (symmetric matrix)
                corr_matrix[j, i] = corr

    # Create DataArray output
    corr_da = xr.DataArray(
        corr_matrix,
        dims=("var1", "var2"),
        coords={"var1": names, "var2": names},
        name="correlation",
    )

    return corr_da


def calc_non_centered_corr(data_input1, data_input2, dim=None):
    """
    Compute the non-centered (uncentered) correlation coefficient between two xarray DataArrays.
    This is equivalent to the cosine similarity, calculated as the sum of the product of the two arrays
    divided by the product of their L2 norms (Euclidean norms), without subtracting the means.

    The formula is:

    .. math::

        r = \\frac{\\sum (x \\cdot y)}{\\sqrt{\\sum x^2} \\cdot \\sqrt{\\sum y^2}}

    Parameters
    ----------
    data_input1 : :py:class:`xarray.DataArray`
        The first input data array to be correlated.
    data_input2 : :py:class:`xarray.DataArray`
        The second input data array to be correlated.

    .. note::
        - Both inputs must be xarray DataArray objects.
        - The arrays must have compatible shapes: if `dim` is specified, it must be a shared dimension;
          if `dim` is None, all dimensions are flattened into a single vector for computation.
        - The result is set to 0 where the denominator (product of norms) is zero to avoid division by zero.

    dim : :py:class:`str` or None, optional
        Dimension(s) over which to compute the correlation. If None (default), the arrays are flattened
        across all dimensions into a single vector before computation. If a string, the correlation is
        computed along the specified dimension, preserving other dimensions.

    Returns
    -------
    corr : :py:class:`xarray.DataArray`
        The non-centered correlation coefficient, with the same dimensions as the input arrays
        (or broadcasted appropriately).

    .. seealso::
        :py:func:`scipy.spatial.distance.cosine <scipy:scipy.spatial.distance.cosine>`
        (for the related cosine distance metric).

    Examples
    --------
    Compute correlation along a dimension:

    >>> import xarray as xr
    >>> import numpy as np
    >>> import easyclimate as ecl
    >>> da1 = xr.DataArray(np.array([[1, 2], [3, 4]]), dims=['x', 'y'])
    >>> da2 = xr.DataArray(np.array([[2, 3], [4, 5]]), dims=['x', 'y'])
    >>> corr = ecl.calc_non_centered_corr(da1, da2, dim='y')
    >>> print(corr)
    <xarray.DataArray (x: 2)> Size: 16B
    array([0.99227788, 0.99951208])

    Flatten and compute scalar correlation:

    >>> corr_flat = calc_non_centered_corr(da1, da2)
    >>> print(corr_flat)
    Dimensions without coordinates: x
    <xarray.DataArray ()> Size: 8B
    array(0.99380799)
    """
    # Ensure inputs are DataArray
    if not isinstance(data_input1, xr.DataArray) or not isinstance(
        data_input2, xr.DataArray
    ):
        raise TypeError("Inputs must be xarray.DataArray objects")

    # If no dimension is specified, flatten the data
    if dim is None:
        da1_flat = data_input1.stack(flat_dim=data_input1.dims)
        da2_flat = data_input2.stack(flat_dim=data_input2.dims)
        numerator = (da1_flat * da2_flat).sum()
        denominator = np.sqrt((da1_flat**2).sum() * (da2_flat**2).sum())
    else:
        # Compute along the specified dimension
        numerator = (data_input1 * data_input2).sum(dim=dim)
        denominator = np.sqrt(
            (data_input1**2).sum(dim=dim) * (data_input2**2).sum(dim=dim)
        )

    # Avoid division by zero error
    return xr.where(denominator != 0, numerator / denominator, 0)


def calc_pattern_corr(
    data_input1: xr.DataArray, data_input2: xr.DataArray, time_dim: str = "time"
):
    """
    Compute the pattern correlation (non-centered) between two xarray DataArrays over their common
    spatial dimensions. This is useful for comparing spatial patterns, such as in climate data.

    It uses the non-centered correlation formula:

    .. math::

        r = \\frac{\\sum (x \\cdot y)}{\\sqrt{\\sum x^2} \\cdot \\sqrt{\\sum y^2}}

    where the summation is over the stacked spatial (pattern) dimensions.

    The spatial pattern dimensions are automatically detected as the intersection of the input dimensions,
    excluding 'time' (if present). Both inputs are stacked along these pattern dimensions into a temporary
    'pattern' dimension, and the non-centered correlation is computed along it.

    - If both inputs lack 'time', the result is a scalar.
    - If one input has 'time' and the other does not, the result preserves the 'time' dimension from the timed input.
    - Broadcasting occurs automatically for compatible shapes.

    Parameters
    ----------
    data_input1 : :py:class:`xarray.DataArray`
        The first input data array (e.g., spatial pattern or time series of patterns).
    data_input2 : :py:class:`xarray.DataArray`
        The second input data array (must have compatible spatial dimensions).
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.

    Returns
    -------
    corr : :py:class:`xarray.DataArray` or scalar
        The pattern correlation coefficient(s). Dimensions match the non-spatial dimensions of the inputs
        (e.g., 'time' if present in one input).

    .. note::
        - Assumes inputs have compatible shapes and the only differing dimension is 'time'.
        - Equivalent to cosine similarity over the spatial pattern.
        - For zero-norm cases, correlation is set to 0.

    .. seealso::
        - `pattern_cor -NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/pattern_cor.shtml>`__
        - :py:func:`calc_non_centered_corr`

    Examples
    --------
    Scalar correlation between two spatial patterns:

    >>> import xarray as xr
    >>> import numpy as np
    >>> import easyclimate as ecl
    >>> # Create a random number generator with a fixed seed.
    >>> rng = np.random.default_rng(42)
    >>> pat1 = xr.DataArray(rng.random((2, 3)), dims=['lat', 'lon'])
    >>> pat2 = xr.DataArray(rng.random((2, 3)), dims=['lat', 'lon'])
    >>> pcc = ecl.calc_pattern_corr(pat1, pat2)
    >>> print(pcc)
    <xarray.DataArray 'pcc' ()> Size: 8B
    array(0.85730639)
    Attributes:
        long_name:  Pattern Correlation Coefficient (non-centered)
        units:      dimensionless

    Time series correlation (one with time):

    >>> # Create a random number generator with a fixed seed.
    >>> rng = np.random.default_rng(42)
    >>> time = xr.DataArray(np.arange(4), dims=['time'])
    >>> timed_pat = xr.DataArray(rng.random((4, 2, 3)), dims=['time', 'lat', 'lon'])
    >>> pcc_time = ecl.calc_pattern_corr(timed_pat, pat2)
    >>> print(pcc_time)
    <xarray.DataArray 'pcc' (time: 4)> Size: 32B
    array([0.85730639, 1.        , 0.78188174, 0.88162673])
    Dimensions without coordinates: time
    Attributes:
        long_name:  Pattern Correlation Coefficient (non-centered)
        units:      dimensionless
    """
    # Ensure inputs are DataArray
    if not isinstance(data_input1, xr.DataArray) or not isinstance(
        data_input2, xr.DataArray
    ):
        raise TypeError("Inputs must be xarray.DataArray objects")

    # Detect common spatial dimensions (exclude 'time')
    dims1 = set(data_input1.dims)
    dims2 = set(data_input2.dims)
    common_dims = dims1.intersection(dims2)
    pattern_dims_set = common_dims - {time_dim}
    if not pattern_dims_set:
        raise ValueError("No common spatial (non-time) dimensions found between inputs")
    pattern_dims = sorted(pattern_dims_set)  # Sorted for consistent stacking

    # Stack both inputs along pattern dimensions
    data1_stacked = data_input1.stack(pattern=tuple(pattern_dims))
    data2_stacked = data_input2.stack(pattern=tuple(pattern_dims))

    # Compute correlation along the pattern dimension
    pcc = calc_non_centered_corr(data1_stacked, data2_stacked, dim="pattern")

    # Add necessary name and attributes
    pcc.name = "pcc"
    pcc.attrs = {
        "long_name": "Pattern Correlation Coefficient (non-centered)",
        "units": "dimensionless",
    }

    return pcc
