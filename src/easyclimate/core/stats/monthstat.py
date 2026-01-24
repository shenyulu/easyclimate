"""
This module calculates statistical values over timesteps of the same month
"""

from __future__ import annotations
import xarray as xr

__all__ = [
    "calc_monthly_mean",
    "calc_monthly_sum",
    "calc_monthly_std",
    "calc_monthly_var",
    "calc_monthly_max",
    "calc_monthly_min",
]


def calc_monthly_mean(data_input: xr.DataArray, dim: str = "time", **kwargs):
    """
    Calculate monthly mean.

    For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same month it is:

    .. math::
        o(t, x) = \\mathrm{mean} \\left \\lbrace i(t', x), t_1 < t' \\leqslant t_n \\right\\rbrace

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

        .. note::

            The recommended frequence of the `data_input` is daily.

    dim: :py:class:`str <str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.

    Examples
    --------
    >>> import xarray as xr
    >>> import numpy as np
    >>> import pandas as pd
    >>> import easyclimate as ecl
    >>> # Create sample data with daily frequency
    >>> time_index = pd.date_range('2020-01-01', '2020-03-31', freq='D')
    >>> rng = np.random.default_rng(42)
    >>> data = rng.random((len(time_index), 3, 3))
    >>> da = xr.DataArray(data, dims=['time', 'x', 'y'], coords={'time': time_index})
    >>> # Calculate monthly mean
    >>> monthly_mean = ecl.calc_monthly_mean(da)
    >>> print(monthly_mean)
    <xarray.DataArray (time: 3, x: 3, y: 3)> Size: 216B
    array([[[0.50635175, 0.45767908, 0.50271707],
            [0.52091523, 0.44830133, 0.44293946],
            [0.47944591, 0.53314083, 0.48073062]],
        [[0.53431127, 0.48259521, 0.47464862],
            [0.41070456, 0.51619935, 0.4872374 ],
            [0.6009132 , 0.43963445, 0.6028882 ]],
        [[0.48363606, 0.60589154, 0.42622008],
            [0.47519641, 0.46989711, 0.45327877],
            [0.44193025, 0.45050389, 0.641573  ]]])
    Coordinates:
    * time     (time) datetime64[ns] 24B 2020-01-01 2020-02-01 2020-03-01
    Dimensions without coordinates: x, y

    .. seealso::
        :py:func:`numpy.mean <numpy:numpy.mean>`, :py:func:`dask.array.mean <dask:dask.array.mean>`,
        :py:meth:`xarray.DataArray.mean <xarray:xarray.DataArray.mean>`,
        :py:meth:`xarray.core.groupby.DataArrayGroupBy.mean <xarray:xarray.core.groupby.DataArrayGroupBy.mean>`.
    """
    ds_mon = data_input.resample({dim: "MS"}).mean(**kwargs)
    return ds_mon


def calc_monthly_sum(data_input: xr.DataArray, dim: str = "time", **kwargs):
    """
    Calculate monthly sum.

    For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same month it is:

    .. math::
        o(t, x) = \\mathrm{sum} \\left \\lbrace i(t', x), t_1 < t' \\leqslant t_n \\right\\rbrace

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

        .. note::

            The recommended frequence of the `data_input` is daily.

    dim: :py:class:`str <str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.

    Examples
    --------
    >>> import xarray as xr
    >>> import numpy as np
    >>> import pandas as pd
    >>> import easyclimate as ecl
    >>> # Create sample data with daily frequency
    >>> time_index = pd.date_range('2020-01-01', '2020-03-31', freq='D')
    >>> rng = np.random.default_rng(42)
    >>> data = rng.random((len(time_index), 3, 3))
    >>> da = xr.DataArray(data, dims=['time', 'x', 'y'], coords={'time': time_index})
    >>> # Calculate monthly sum
    >>> monthly_sum = ecl.calc_monthly_sum(da)
    >>> print(monthly_sum)
    <xarray.DataArray (time: 3, x: 3, y: 3)> Size: 216B
    array([[[15.69690421, 14.18805154, 15.58422905],
            [16.14837203, 13.89734131, 13.7311233 ],
            [14.86282316, 16.52736588, 14.90264935]],
        [[15.49502686, 13.99526102, 13.76481005],
            [11.91043229, 14.96978113, 14.1298847 ],
            [17.42648281, 12.7493991 , 17.48375789]],
        [[14.99271788, 18.78263759, 13.21282259],
            [14.73108859, 14.56681052, 14.05164186],
            [13.69983774, 13.96562059, 19.88876291]]])
    Coordinates:
    * time     (time) datetime64[ns] 24B 2020-01-01 2020-02-01 2020-03-01
    Dimensions without coordinates: x, y

    .. seealso::
        :py:func:`numpy.sum <numpy:numpy.sum>`, :py:func:`dask.array.sum <dask:dask.array.sum>`,
        :py:meth:`xarray.DataArray.sum <xarray:xarray.DataArray.sum>`,
        :py:meth:`xarray.core.groupby.DataArrayGroupBy.sum <xarray:xarray.core.groupby.DataArrayGroupBy.sum>`.
    """
    ds_mon = data_input.resample({dim: "MS"}).sum(**kwargs)
    return ds_mon


def calc_monthly_std(data_input: xr.DataArray, dim: str = "time", **kwargs):
    """
    Calculate monthly standard deviation.

    For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same month it is:

    .. math::
        o(t, x) = \\mathrm{std} \\left \\lbrace i(t', x), t_1 < t' \\leqslant t_n \\right\\rbrace

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

        .. note::

            The recommended frequence of the `data_input` is daily.

    dim: :py:class:`str <str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.

    Examples
    --------
    >>> import xarray as xr
    >>> import numpy as np
    >>> import pandas as pd
    >>> import easyclimate as ecl
    >>> # Create sample data with daily frequency
    >>> time_index = pd.date_range('2020-01-01', '2020-03-31', freq='D')
    >>> rng = np.random.default_rng(42)
    >>> data = rng.random((len(time_index), 3, 3))
    >>> da = xr.DataArray(data, dims=['time', 'x', 'y'], coords={'time': time_index})
    >>> # Calculate monthly std
    >>> monthly_std = ecl.calc_monthly_std(da)
    >>> print(monthly_std)
    <xarray.DataArray (time: 3, x: 3, y: 3)> Size: 216B
    array([[[0.30528844, 0.29905447, 0.26472868],
            [0.23456056, 0.30879525, 0.29333846],
            [0.26139562, 0.306974  , 0.27987361]],
        [[0.30196879, 0.24783961, 0.26078164],
            [0.2708643 , 0.3012602 , 0.29801453],
            [0.24816804, 0.33863555, 0.25623523]],
        [[0.29399346, 0.31595077, 0.30336434],
            [0.31117807, 0.3130123 , 0.28909393],
            [0.27104435, 0.26864038, 0.22912052]]])
    Coordinates:
    * time     (time) datetime64[ns] 24B 2020-01-01 2020-02-01 2020-03-01
    Dimensions without coordinates: x, y

    .. seealso::
        :py:func:`numpy.std <numpy:numpy.std>`, :py:func:`dask.array.std <dask:dask.array.std>`,
        :py:meth:`xarray.DataArray.std <xarray:xarray.DataArray.std>`,
        :py:meth:`xarray.core.groupby.DataArrayGroupBy.std <xarray:xarray.core.groupby.DataArrayGroupBy.std>`.
    """
    ds_mon = data_input.resample({dim: "MS"}).std(**kwargs)
    return ds_mon


def calc_monthly_var(data_input: xr.DataArray, dim: str = "time", **kwargs):
    """
    Calculate monthly variance.

    For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same month it is:

    .. math::
        o(t, x) = \\mathrm{var} \\left \\lbrace i(t', x), t_1 < t' \\leqslant t_n \\right\\rbrace

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

        .. note::

            The recommended frequence of the `data_input` is daily.

    dim: :py:class:`str <str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.

    Examples
    --------
    >>> import xarray as xr
    >>> import numpy as np
    >>> import pandas as pd
    >>> import easyclimate as ecl
    >>> # Create sample data with daily frequency
    >>> time_index = pd.date_range('2020-01-01', '2020-03-31', freq='D')
    >>> rng = np.random.default_rng(42)
    >>> data = rng.random((len(time_index), 3, 3))
    >>> da = xr.DataArray(data, dims=['time', 'x', 'y'], coords={'time': time_index})
    >>> # Calculate monthly var
    >>> monthly_var = ecl.calc_monthly_var(da)
    >>> print(monthly_var)
    <xarray.DataArray (time: 3, x: 3, y: 3)> Size: 216B
    array([[[0.09320103, 0.08943358, 0.07008127],
            [0.05501865, 0.09535451, 0.08604745],
            [0.06832767, 0.09423304, 0.07832924]],
        [[0.09118515, 0.06142447, 0.06800706],
            [0.07336747, 0.09075771, 0.08881266],
            [0.06158738, 0.11467404, 0.0656565 ]],
        [[0.08643216, 0.09982489, 0.09202992],
            [0.09683179, 0.0979767 , 0.0835753 ],
            [0.07346504, 0.07216766, 0.05249621]]])
    Coordinates:
    * time     (time) datetime64[ns] 24B 2020-01-01 2020-02-01 2020-03-01
    Dimensions without coordinates: x, y

    .. seealso::
        :py:func:`numpy.var <numpy:numpy.var>`, :py:func:`dask.array.var <dask:dask.array.var>`,
        :py:meth:`xarray.DataArray.var <xarray:xarray.DataArray.var>`,
        :py:meth:`xarray.core.groupby.DataArrayGroupBy.var <xarray:xarray.core.groupby.DataArrayGroupBy.var>`.
    """
    ds_mon = data_input.resample({dim: "MS"}).var(**kwargs)
    return ds_mon


def calc_monthly_max(data_input: xr.DataArray, dim: str = "time", **kwargs):
    """
    Calculate monthly maximum.

    For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same month it is:

    .. math::
        o(t, x) = \\mathrm{max} \\left \\lbrace i(t', x), t_1 < t' \\leqslant t_n \\right\\rbrace

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

        .. note::

            The recommended frequence of the `data_input` is daily.

    dim: :py:class:`str <str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.

    Examples
    --------
    >>> import xarray as xr
    >>> import numpy as np
    >>> import pandas as pd
    >>> import easyclimate as ecl
    >>> # Create sample data with daily frequency
    >>> time_index = pd.date_range('2020-01-01', '2020-03-31', freq='D')
    >>> rng = np.random.default_rng(42)
    >>> data = rng.random((len(time_index), 3, 3))
    >>> da = xr.DataArray(data, dims=['time', 'x', 'y'], coords={'time': time_index})
    >>> # Calculate monthly max
    >>> monthly_max = ecl.calc_monthly_max(da)
    >>> print(monthly_max)
    <xarray.DataArray (time: 3, x: 3, y: 3)> Size: 216B
    <xarray.DataArray (time: 3, x: 3, y: 3)> Size: 216B
    array([[[0.96189766, 0.95855921, 0.93604357],
            [0.97182643, 0.97069802, 0.97562235],
            [0.99237556, 0.96623191, 0.91601185]],
        [[0.95119466, 0.88414571, 0.85053368],
            [0.94602445, 0.99910473, 0.99546447],
            [0.98002718, 0.98663154, 0.99703466]],
        [[0.989133  , 0.99874337, 0.99308458],
            [0.99032166, 0.98180595, 0.92746046],
            [0.99758004, 0.91879368, 0.99470175]]])
    Coordinates:
    * time     (time) datetime64[ns] 24B 2020-01-01 2020-02-01 2020-03-01
    Dimensions without coordinates: x, y

    .. seealso::
        :py:func:`numpy.max <numpy:numpy.max>`, :py:func:`dask.array.max <dask:dask.array.max>`,
        :py:meth:`xarray.DataArray.max <xarray:xarray.DataArray.max>`,
        :py:meth:`xarray.core.groupby.DataArrayGroupBy.max <xarray:xarray.core.groupby.DataArrayGroupBy.max>`.
    """
    ds_mon = data_input.resample({dim: "MS"}).max(**kwargs)
    return ds_mon


def calc_monthly_min(data_input: xr.DataArray, dim: str = "time", **kwargs):
    """
    Calculate monthly minimum.

    For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same month it is:

    .. math::
        o(t, x) = \\mathrm{min} \\left \\lbrace i(t', x), t_1 < t' \\leqslant t_n \\right\\rbrace

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

        .. note::

            The recommended frequence of the `data_input` is daily.

    dim: :py:class:`str <str>`
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data.
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.

    Examples
    --------
    >>> import xarray as xr
    >>> import numpy as np
    >>> import pandas as pd
    >>> import easyclimate as ecl
    >>> # Create sample data with daily frequency
    >>> time_index = pd.date_range('2020-01-01', '2020-03-31', freq='D')
    >>> rng = np.random.default_rng(42)
    >>> data = rng.random((len(time_index), 3, 3))
    >>> da = xr.DataArray(data, dims=['time', 'x', 'y'], coords={'time': time_index})
    >>> # Calculate monthly min
    >>> monthly_min = ecl.calc_monthly_min(da)
    >>> print(monthly_min)
    <xarray.DataArray (time: 3, x: 3, y: 3)> Size: 216B
    array([[[0.02280387, 0.02485949, 0.05338193],
            [0.02271207, 0.01783678, 0.02161208],
            [0.00736227, 0.04161417, 0.02114802]],
        [[0.0289995 , 0.04347506, 0.0401513 ],
            [0.01072764, 0.09172101, 0.01468284],
            [0.07205915, 0.01230269, 0.00542983]],
        [[0.0040076 , 0.0165798 , 0.00166071],
            [0.04896371, 0.01903415, 0.00123306],
            [0.04737402, 0.00450012, 0.22825288]]])
    Coordinates:
    * time     (time) datetime64[ns] 24B 2020-01-01 2020-02-01 2020-03-01
    Dimensions without coordinates: x, y

    .. seealso::
        :py:func:`numpy.min <numpy:numpy.min>`, :py:func:`dask.array.min <dask:dask.array.min>`,
        :py:meth:`xarray.DataArray.min <xarray:xarray.DataArray.min>`,
        :py:meth:`xarray.core.groupby.DataArrayGroupBy.min <xarray:xarray.core.groupby.DataArrayGroupBy.min>`.
    """
    ds_mon = data_input.resample({dim: "MS"}).min(**kwargs)
    return ds_mon
