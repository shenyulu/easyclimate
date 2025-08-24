"""
Interpolate from pressure layer to altitude layers
"""

from __future__ import annotations

import numpy as np
import xarray as xr

__all__ = ["interp1d_vertical_pressure2altitude"]


def interp1d_vertical_pressure2altitude(
    z_data: xr.DataArray,
    variable_data: xr.DataArray,
    target_heights: np.array,
    vertical_input_dim: str,
    vertical_output_dim: str,
    kind: str = "cubic",
    bounds_error=None,
    fill_value="extrapolate",
    assume_sorted: bool = False,
) -> xr.DataArray:
    """
    Interpolating variables from the pressure levels (e.g., sigma vertical coordinate) to the altitude levels by 1-D function.

    z_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The vertical atmospheric geopotential height on each pressure level.
    variable_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        variable to be interpolated on model levels.

        .. note::
            The shape of `z_data` `z_data` and `variable_data` shuold be the same.

    target_heights: :py:class:`numpy.array<numpy.array>`.
        Altitude interpolation sampling range. e.g., `np.linspace(0, 15000, 100)`.

        .. note::
            The unit of `target_heights` shuold be same as `z_data`, e.g., the unit of `target_heights` is `meter`, and the unit of `z_data` shuold be `meter`.

    vertical_input_dim: :py:class:`str <str>`.
        The name of the standard pressure levels dimension, often assigned the value `'plev'` (Pa) or `'lev'` (hPa).

        .. note::
            The `vertical_input_dim` `z_data` and `variable_data` shuold be the same.

    vertical_output_dim: :py:class:`str <str>`.
        The name of the altitude dimension, often assigned the value `'altitude'`.
    kind: :py:class:`str <str>` or :py:class:`int <int>`, default: `'cubic'`.
        Specifies the kind of interpolation as a string or as an integer specifying the order of the spline interpolator to use.
        The string has to be one of 'linear', 'nearest', 'nearest-up', 'zero', 'slinear', 'quadratic', 'cubic', 'previous',
        or 'next'. 'zero', 'slinear', 'quadratic' and 'cubic' refer to a spline interpolation of zeroth, first, second or third order;
        'previous' and 'next' simply return the previous or next value of the point;
        'nearest-up' and 'nearest' differ when interpolating half-integers (e.g. 0.5, 1.5) in that 'nearest-up' rounds up and 'nearest' rounds down.
    bounds_error: :py:class:`bool <bool>`, default: `None`.
        If True, a ValueError is raised any time interpolation is attempted on a value outside of the range of `pressure_data` (where extrapolation is necessary).
        If False, out of bounds values are assigned fill_value. By default, an error is raised unless `fill_value="extrapolate"`.
    fill_value: array-like or (array-like, array_like) or "extrapolate", default: `extrapolate`.
        - If a ndarray (or float), this value will be used to fill in for requested points outside of the data range. If not provided, then the default is NaN. The array-like must broadcast properly to the dimensions of the non-interpolation axes.
        - If a two-element tuple, then the first element is used as a fill value for x_new < x[0] and the second element is used for x_new > x[-1]. Anything that is not a 2-element tuple (e.g., list or ndarray, regardless of shape) is taken to be a single array-like argument meant to be used for both bounds as below, above = fill_value, fill_value. Using a two-element tuple or ndarray requires bounds_error=False.
        - If "extrapolate", then points outside the data range will be extrapolated.
    assume_sortedbool: :py:class:`bool <bool>`, default: `False`.
        If `False`, values of x can be in any order and they are sorted first. If `True`, `pressure_data` has to be an array of monotonically increasing values.

    .. seealso::
        - `scipy.interpolate.interp1d <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html#scipy-interpolate-interp1d>`__

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_interp.py
    """
    from scipy import interpolate

    def _interp_height(z_data, variable_data, target_heights=target_heights):
        object_interp = interpolate.interp1d(
            z_data,
            variable_data,
            kind=kind,
            bounds_error=bounds_error,
            fill_value=fill_value,
            assume_sorted=assume_sorted,
        )
        result_interp = object_interp(target_heights)
        return result_interp

    interped_data = xr.apply_ufunc(
        _interp_height,
        z_data,
        variable_data,
        input_core_dims=[[vertical_input_dim], [vertical_input_dim]],
        output_core_dims=[[vertical_output_dim]],
        output_dtypes=["float64"],
        dask="parallelized",
        vectorize=True,
        dask_gufunc_kwargs={
            "output_sizes": {vertical_output_dim: len(target_heights)},
            "allow_rechunk": True,
        },
    )
    interped_data = interped_data.assign_coords({vertical_output_dim: target_heights})
    return interped_data
