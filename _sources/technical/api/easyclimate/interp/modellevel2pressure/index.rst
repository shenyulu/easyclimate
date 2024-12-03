easyclimate.interp.modellevel2pressure
======================================

.. py:module:: easyclimate.interp.modellevel2pressure

.. autoapi-nested-parse::

   Interpolate from model level layers to pressure layer



Functions
---------

.. autoapisummary::

   easyclimate.interp.modellevel2pressure.interp1d_vertical_model2pressure


Module Contents
---------------

.. py:function:: interp1d_vertical_model2pressure(pressure_data: xarray.DataArray, variable_data: xarray.DataArray, vertical_input_dim: str, vertical_output_dim: str, vertical_output_level: list[int | float], kind: str = 'linear', bounds_error=None, fill_value=np.nan, assume_sorted: bool = False) -> xarray.DataArray

   Interpolating variables from the model levels (e.g., sigma vertical coordinate) to the standard pressure levels by 1-D function.

   pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The pressure on each model level.
   variable_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       variable to be interpolated on model levels.
   vertical_input_dim: :py:class:`str <str>`.
       The name of the model levels dimension.
   vertical_output_dim: :py:class:`str <str>`.
       The name of the standard pressure levels dimension, often assigned the value `'plev'` (Pa) or `'lev'` (hPa).
   vertical_output_level: :py:class:`list[int | float]`.
       Customized standard pressure levels to be interpolated, e.g., `[5000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 60000, 70000, 85000, 92500, 100000]`.

       .. note::
           The unit of `vertical_output_level` shuold be same as `pressure_data`, e.g., the unit of `pressure_data` is `Pa`, and the unit of `vertical_output_level` shuold be `Pa`.

   kind: :py:class:`str <str>` or :py:class:`int <int>`, default: `'linear'`.
       Specifies the kind of interpolation as a string or as an integer specifying the order of the spline interpolator to use.
       The string has to be one of 'linear', 'nearest', 'nearest-up', 'zero', 'slinear', 'quadratic', 'cubic', 'previous',
       or 'next'. 'zero', 'slinear', 'quadratic' and 'cubic' refer to a spline interpolation of zeroth, first, second or third order;
       'previous' and 'next' simply return the previous or next value of the point;
       'nearest-up' and 'nearest' differ when interpolating half-integers (e.g. 0.5, 1.5) in that 'nearest-up' rounds up and 'nearest' rounds down.
   bounds_error: :py:class:`bool <bool>`, default: `None`.
       If True, a ValueError is raised any time interpolation is attempted on a value outside of the range of `pressure_data` (where extrapolation is necessary).
       If False, out of bounds values are assigned fill_value. By default, an error is raised unless `fill_value="extrapolate"`.
   fill_value: array-like or (array-like, array_like) or "extrapolate", default: `np.nan`.
       - If a ndarray (or float), this value will be used to fill in for requested points outside of the data range. If not provided, then the default is NaN. The array-like must broadcast properly to the dimensions of the non-interpolation axes.
       - If a two-element tuple, then the first element is used as a fill value for x_new < x[0] and the second element is used for x_new > x[-1]. Anything that is not a 2-element tuple (e.g., list or ndarray, regardless of shape) is taken to be a single array-like argument meant to be used for both bounds as below, above = fill_value, fill_value. Using a two-element tuple or ndarray requires bounds_error=False.
       - If "extrapolate", then points outside the data range will be extrapolated.
   assume_sortedbool: :py:class:`bool <bool>`, default: `False`.
       If `False`, values of x can be in any order and they are sorted first. If `True`, `pressure_data` has to be an array of monotonically increasing values.

   .. seealso::
       - `scipy.interpolate.interp1d <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html#scipy-interpolate-interp1d>`__


