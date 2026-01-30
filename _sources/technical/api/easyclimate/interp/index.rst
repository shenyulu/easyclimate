easyclimate.interp
==================

.. py:module:: easyclimate.interp


Submodules
----------

.. toctree::
   :maxdepth: 1

   /technical/api/easyclimate/interp/barnes/index
   /technical/api/easyclimate/interp/mesh2mesh/index
   /technical/api/easyclimate/interp/mesh2point/index
   /technical/api/easyclimate/interp/modellevel2pressure/index
   /technical/api/easyclimate/interp/pressure2altitude/index
   /technical/api/easyclimate/interp/vinth2p/index


Functions
---------

.. autoapisummary::

   easyclimate.interp.interp_mesh2mesh
   easyclimate.interp.interp_spatial_barnes
   easyclimate.interp.interp_spatial_barnesS2
   easyclimate.interp.interp1d_vertical_model2pressure
   easyclimate.interp.interp1d_vertical_pressure2altitude
   easyclimate.interp.interp1d_vertical_pressure2altitude_linear_rs
   easyclimate.interp.interp_mesh2point
   easyclimate.interp.interp_mesh2point_withtime
   easyclimate.interp.interp_vinth2p_dp
   easyclimate.interp.interp_vinth2p_ecmwf
   easyclimate.interp.interp_vintp2p_ecmwf


Package Contents
----------------

.. py:function:: interp_mesh2mesh(data_input: xarray.DataArray | xarray.Dataset, target_grid: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon', lat_dim: str = 'lat', method: Literal['linear', 'nearest', 'cubic', 'conservative'] = 'linear')

   Regridding regular or lat-lon grid data.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   target_grid: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       Target grid to be regridding.

       :py:class:`xarray.DataArray<xarray.DataArray>` version sample

       .. code:: python

           target_grid = xr.DataArray(
               dims=('lat', 'lon'),
               coords={'lat': np.arange(-89, 89, 3) + 1 / 1.0, 'lon': np.arange(-180, 180, 3) + 1 / 1.0}
           )

       :py:class:`xarray.Dataset<xarray.Dataset>` version sample

       .. code:: python

           target_grid = xr.Dataset()
           target_grid['lat'] = np.arange(-89, 89, 3) + 1 / 1.0
           target_grid['lon'] = np.arange(-180, 180, 3) + 1 / 1.0

   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   method: :py:class:`str <str>`, default: `linear`.
       The methods of regridding.

       - `linear`: linear, bilinear, or higher dimensional linear interpolation.
       - `nearest`: nearest-neighbor regridding.
       - `cubic`: cubic spline regridding.
       - `conservative`: conservative regridding.

   Reference
   --------------
   - https://github.com/EXCITED-CO2/xarray-regrid

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_interp.py


.. py:function:: interp_spatial_barnes(data: pandas.DataFrame, var_name: str, point: list[int], grid_x: float, grid_y: float, resolution: float, sigma: float, lon_dim: str = 'lon', lat_dim: str = 'lat', method: Literal['optimized_convolution', 'convolution', 'radius', 'naive'] = 'optimized_convolution', num_iter: int = 4, max_dist: float = 3.5, min_weight: float = 0.001) -> xarray.DataArray

   Computes the Barnes interpolation for observation values `var_name` taken at sample
   points `data` using Gaussian weights for the width parameter `sigma`.
   The underlying grid embedded in a Euclidean space is given with start point
   `point`, regular x-direction length `grid_x` (degree), regular y-direction length `grid_y` (degree),
   and resolution `resolution`.

   Barnes interpolation is a method that is widely used in geospatial sciences like meteorology
   to remodel data values recorded at irregularly distributed points into a representative
   analytical field. It is defined as

   .. math::
       f(\boldsymbol{x})=\frac{\sum_{k=1}^N f_k\cdot w_k(\boldsymbol{x})}{\sum_{k=1}^N w_k(\boldsymbol{x})}

   with Gaussian weights

   .. math::
       w_k(\boldsymbol{x})=\text{e}^{-\frac{1}{2\sigma^2}\left|x-\boldsymbol{x}_k\right|^2}

   Naive computation of Barnes interpolation leads to an algorithmic complexity of :math:`O(N \times W \times H)`,
   where :math:`N` is the number of sample points and :math:`W \times H` the size of the underlying grid.

   For sufficiently large :math:`n` (in general in the range from 3 to 6) a good approximation of
   Barnes interpolation with a reduced complexity :math:`O(N + W \times H)` can be obtained by the convolutional expression

   .. math::
       f(\boldsymbol{x})\approx \frac{ (\sum_{k=1}^{N}f_k\cdot\delta_{\boldsymbol{x}_k}) *  ( r_n^{*n[x]}(x)\cdot r_n^{*n[y]}(y) )   }{ ( \sum_{k=1}^{N} \delta_{\boldsymbol{x}_k}  ) *  (  r_{n}^{*n[x]}(x)\cdot r_{n}^{*n[y]}(y)  )   }

   where :math:`\delta` is the Dirac impulse function and :math:`r(.)` an elementary rectangular function of a specific length that depends on :math:`\sigma` and :math:`n`.

   data : :py:class:`pandas.DataFrame<pandas.DataFrame>`
       Longitude and latitude grid point discrete data. There should be a similar structure as follows

       +------------+------------+-----------+
       |    lon     |    lat     |    qff    |
       +============+============+===========+
       |   -3.73    |   56.33    |   995.1   |
       +------------+------------+-----------+
       |    2.64    |   47.05    |  1012.5   |
       +------------+------------+-----------+
       |    ...     |   ...      |   ...     |
       +------------+------------+-----------+

       .. note::
           Data points should contain longitude (`lon`), latitude (`lat`) and data variables (the above data variable name is `qff`).

   var_name: :py:class:`str <str>`
       The name of the data variable. This should match the one in the parameter `data`.
   lat_dim: :py:class:`str <str>`.
       Latitude dimension name. This should match the one in the parameter `data`. By default is `lat`.
   lon_dim: :py:class:`str <str>`.
       Longitude dimension name. This should match the one in the parameter `data`. By default is `lon`.
   point : :py:class:`numpy.ndarray <numpy.ndarray>`.
       A 1-dimensional array of size 2 containing the coordinates of the
       start point of the grid to be used.
   grid_x : :py:class:`int <int>`.
       Length in degrees in the x-direction of the interpolated rectangular grid.
   grid_y : :py:class:`int <int>`.
       Length in degrees in the y-direction of the interpolated rectangular grid.
   resolution: :py:class:`float <float>`
       Grid resolution. The distance between regular grid points is the reciprocal of the value. Common values: 4.0, 8.0, 16.0, 32.0, 64.0.
   sigma : :py:class:`float <float>`
       The Gaussian width parameter to be used. Common values: 0.25, 0.5, 1.0, 2.0, 4.0.
   method : {'optimized_convolution', 'convolution', 'radius', 'naive'}, default: 'optimized_convolution'.
       Designates the Barnes interpolation method to be used. The possible
       implementations that can be chosen are 'naive' for the straightforward
       implementation (algorithm A from paper), 'radius' to consider only sample
       points within a specific radius of influence, both with an algorithmic
       complexity of :math:`O(N \times W \times H)`.
       The choice 'convolution' implements algorithm B specified in the paper
       and 'optimized_convolution' is its optimization by appending tail values
       to the rectangular kernel. The latter two algorithms reduce the complexity
       down to :math:`O(N + W \times H)`.
       The default is 'optimized_convolution'.
   num_iter : :py:class:`int <int>`, optional
       The number of performed self-convolutions of the underlying rect-kernel.
       Applies only if method is 'optimized_convolution' or 'convolution'.
       The default is 4. Applies only to Convol interpolations: one of 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50.
   max_dist : float, optional
       The maximum distance between a grid point and the next sample point for which
       the Barnes interpolation is still calculated. Specified in sigma distances.
       Applies only if method is 'optimized_convolution' or 'convolution'.
       The default is 3.5, i.e. the maximum distance is 3.5 * sigma.
   min_weight : :py:class:`float <float>`, optional
       Choose radius of influence such that Gaussian weight of considered sample
       points is greater than `min_weight`.
       Applies only if method is 'radius'. Recommended values are 0.001 and less.
       The default is 0.001, which corresponds to a radius of 3.717 * sigma.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::
       - https://github.com/MeteoSwiss/fast-barnes-py
       - Zürcher, B. K.: Fast approximate Barnes interpolation: illustrated by Python-Numba implementation fast-barnes-py v1.0, Geosci. Model Dev., 16, 1697–1711, https://doi.org/10.5194/gmd-16-1697-2023, 2023.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_interp.py


.. py:function:: interp_spatial_barnesS2(data: pandas.DataFrame, var_name: str, point: list[int], grid_x: float, grid_y: float, resolution: float, sigma: float, lon_dim: str = 'lon', lat_dim: str = 'lat', method: Literal['optimized_convolution_S2', 'naive_S2'] = 'optimized_convolution_S2', num_iter: int = 4, max_dist: float = 3.5, resample: bool = True) -> xarray.DataArray

   Computes the Barnes interpolation for observation values `var_name` taken at sample
   points `data` using Gaussian weights for the width parameter `sigma`.

   The underlying grid embedded on the unit sphere :math:`S^2` and thus inherits the
   spherical distance measure (taken in degrees). The grid is given by the start
   point `point`, regular x-direction length `grid_x` (degree), regular y-direction length `grid_y` (degree),
   and resolution `resolution`.

   Parameters
   ----------
   data : :py:class:`pandas.DataFrame<pandas.DataFrame>`
       Longitude and latitude grid point discrete data. There should be a similar structure as follows

       +------------+------------+-----------+
       |    lon     |    lat     |    qff    |
       +============+============+===========+
       |   -3.73    |   56.33    |   995.1   |
       +------------+------------+-----------+
       |    2.64    |   47.05    |  1012.5   |
       +------------+------------+-----------+
       |    ...     |   ...      |   ...     |
       +------------+------------+-----------+

       .. note::
           Data points should contain longitude (`lon`), latitude (`lat`) and data variables (the above data variable name is `qff`).

   var_name: :py:class:`str <str>`
       The name of the data variable. This should match the one in the parameter `data`.
   lat_dim: :py:class:`str <str>`.
       Latitude dimension name. This should match the one in the parameter `data`. By default is `lat`.
   lon_dim: :py:class:`str <str>`.
       Longitude dimension name. This should match the one in the parameter `data`. By default is `lon`.
   point : :py:class:`numpy.ndarray <numpy.ndarray>`
       A 1-dimensional array of size 2 containing the coordinates of the
       start point of the grid to be used.
   grid_x : :py:class:`int <int>`.
       Length in degrees in the x-direction of the interpolated rectangular grid.
   grid_y : :py:class:`int <int>`.
       Length in degrees in the y-direction of the interpolated rectangular grid.
   resolution: :py:class:`float <float>`
       Grid resolution. The distance between regular grid points is the reciprocal of the value. Common values: 4.0, 8.0, 16.0, 32.0, 64.0.
   sigma : :py:class:`float <float>`
       The Gaussian width parameter to be used. Common values: 0.25, 0.5, 1.0, 2.0, 4.0.
   method : {'optimized_convolution_S2', 'naive_S2'}, default: 'optimized_convolution_S2'.
       Designates the Barnes interpolation method to be used. The possible
       implementations that can be chosen are 'naive_S2' for the straightforward
       implementation (algorithm A from the paper) with an algorithmic complexity
       of :math:`O(N \times W \times H)`.
       The choice 'optimized_convolution_S2' implements the optimized algorithm B
       specified in the paper by appending tail values to the rectangular kernel.
       The latter algorithm has a reduced complexity of :math:`O(N + W \times H)`.
       The default is 'optimized_convolution_S2'.
   num_iter : :py:class:`int <int>`, optional, default: 4.
       The number of performed self-convolutions of the underlying rect-kernel.
       Applies only if method is 'optimized_convolution_S2'.
       The default is 4. Applies only to Convol interpolations: one of 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50.
   max_dist : float, optional
       The maximum distance between a grid point and the next sample point for which
       the Barnes interpolation is still calculated. Specified in sigma distances.
       The default is 3.5, i.e. the maximum distance is 3.5 * sigma.
   resample : :py:class:`bool <bool>`, optional, default: `True`.
       Specifies whether to resample Lambert grid field to lonlat grid.
       Applies only if method is 'optimized_convolution_S2'.
       The default is True.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::
       - https://github.com/MeteoSwiss/fast-barnes-py
       - Zürcher, B. K.: Fast approximate Barnes interpolation: illustrated by Python-Numba implementation fast-barnes-py v1.0, Geosci. Model Dev., 16, 1697–1711, https://doi.org/10.5194/gmd-16-1697-2023, 2023.


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

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_interp.py


.. py:function:: interp1d_vertical_pressure2altitude(z_data: xarray.DataArray, variable_data: xarray.DataArray, target_heights: numpy.array, vertical_input_dim: str, vertical_output_dim: str, kind: str = 'cubic', bounds_error=None, fill_value='extrapolate', assume_sorted: bool = False) -> xarray.DataArray

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


.. py:function:: interp1d_vertical_pressure2altitude_linear_rs(z_data: xarray.DataArray, variable_data: xarray.DataArray, target_heights: numpy.array, vertical_input_dim: str, vertical_output_dim: str) -> xarray.DataArray

   Interpolating variables from the pressure levels (e.g., sigma vertical coordinate) to the altitude levels using `easyclimate-rust`.

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

   .. seealso::
       - `scipy.interpolate.interp1d <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html#scipy-interpolate-interp1d>`__

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_interp.py


.. py:function:: interp_mesh2point(data_input: xarray.DataArray, df: pandas.DataFrame, lon_dim_mesh: str = 'lon', lat_dim_mesh: str = 'lat', lon_dim_df: str = 'lon', lat_dim_df: str = 'lat', method: Literal['linear', 'nearest', 'slinear', 'cubic', 'quintic', 'pchip'] = 'linear')

   Interpolate values from a regular grid/mesh to specific point locations.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray <xarray.DataArray>`.
       Input grid data with latitude and longitude dimensions.
   df : :py:class:`pandas.DataFrame <pandas.DataFrame>`.
       DataFrame containing point locations with latitude and longitude columns.
   lon_dim_mesh : :py:class:`str <str>`, optional
       Name of the longitude dimension in the input grid, by default ``"lon"``.
   lat_dim_mesh : :py:class:`str <str>`, optional
       Name of the latitude dimension in the input grid, by default ``"lat"``.
   lon_dim_df : :py:class:`str <str>`, optional
       Name of the longitude column in the DataFrame, by default ``"lon"``.
   lat_dim_df : :py:class:`str <str>`, optional
       Name of the latitude column in the DataFrame, by default ``"lat"``.
   method : :py:class:`str <str>`, optional
       Interpolation method to use. Options are:

       - ``"linear"``: bilinear interpolation (default)
       - ``"nearest"``: nearest neighbor
       - ``"slinear"``: spline linear
       - ``"cubic"``: spline cubic
       - ``"quintic"``: spline quintic
       - ``"pchip"``: piecewise cubic Hermite interpolating polynomial

   Returns
   -------
   :py:class:`pandas.DataFrame <pandas.DataFrame>`.
       Original DataFrame with an additional column ``"interpolated_value"`` containing
       the interpolated values. Points outside the grid range will have NaN values.

   Raises
   ------
   ValueError
       If no points in the DataFrame fall within the grid's spatial extent.

   Examples
   --------
   >>> import xarray as xr
   >>> import pandas as pd
   >>> # Create sample grid data
   >>> lats = np.linspace(-90, 90, 181)
   >>> lons = np.linspace(-180, 180, 361)
   >>> data = xr.DataArray(np.random.rand(181, 361), dims=['lat', 'lon'],
   ...                     coords={'lat': lats, 'lon': lons})
   >>> # Create sample points
   >>> points = pd.DataFrame({'lat': [45.5, 30.2], 'lon': [-120.3, 150.7]})
   >>> # Interpolate
   >>> result = interp_mesh2point(data, points)

   .. seealso::

       - https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html
       - https://www.ncl.ucar.edu/Applications/station.shtml

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_interp_mesh2point.py


.. py:function:: interp_mesh2point_withtime(data_input: xarray.DataArray, stations_df: pandas.DataFrame, lon_dim_df: str, lat_dim_df: str, station_dim_df: str, lon_dim_mesh: str = 'lon', lat_dim_mesh: str = 'lat', time_dim_mesh: str = 'time', method: Literal['linear', 'nearest', 'slinear', 'cubic', 'quintic', 'pchip'] = 'linear') -> xarray.DataArray

   Interpolate gridded data to specific point locations with time series preservation.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray <xarray.DataArray>`.
       Input grid data with time, latitude and longitude dimensions.
   stations_df : :py:class:`pandas.DataFrame <pandas.DataFrame>`.
       DataFrame containing station locations with ID, latitude and longitude.
   lon_dim_df : :py:class:`str <str>`.
       Name of the longitude column in the DataFrame.
   lat_dim_df : :py:class:`str <str>`.
       Name of the latitude column in the DataFrame.
   station_dim_df : :py:class:`str <str>`.
       Name of the station ID column in the DataFrame.
   lon_dim_mesh : :py:class:`str <str>`, default: "lon"
       Name of the longitude dimension in the input grid.
   lat_dim_mesh : :py:class:`str <str>`, default: "lat"
       Name of the latitude dimension in the input grid.
   time_dim_mesh : :py:class:`str <str>`, default: "time"
       Name of the time dimension in the input grid.
   method : :py:class:`str <str>`, optional
       Interpolation method:

       - "linear": bilinear interpolation (default)
       - "nearest": nearest neighbor
       - "slinear": spline linear
       - "cubic": spline cubic
       - "quintic": spline quintic
       - "pchip": piecewise cubic Hermite interpolating polynomial

   Returns
   -------
   :py:class:`xarray.Dataset <xarray.Dataset>`.
       :py:class:`xarray.Dataset <xarray.Dataset>` with dimensions ``(station, time)`` containing:

       - Interpolated values
       - Station coordinates


   Examples
   --------
   >>> import xarray as xr
   >>> import pandas as pd
   >>> import numpy as np
   >>> # Create sample grid
   >>> times = pd.date_range("2020-01-01", periods=5)
   >>> lats = np.linspace(-45, -10, 100)
   >>> lons = np.linspace(110, 156, 120)
   >>> data = xr.DataArray(
   ...     np.random.rand(5, 100, 120),
   ...     dims=["time", "lat", "lon"],
   ...     coords={"time": times, "lat": lats, "lon": lons}
   ... )
   >>> # Create stations
   >>> stations = pd.DataFrame({
   ...     "station_id_col": [1001, 1002],
   ...     "lat_col": [-15.5, -20.3],
   ...     "lon_col": [125.5, 130.2]
   ... })
   >>> # Interpolate
   >>> result = interp_mesh2point_withtime(
   ...     data, stations,
   ...     lon_dim_df="lon_col", lat_dim_df="lat_col", station_dim_df="station_id_col"
   ... )

   .. seealso::

       :py:class:`scipy.interpolate.RegularGridInterpolator <scipy.interpolate.RegularGridInterpolator>`.


.. py:function:: interp_vinth2p_dp(data_input: xarray.DataArray, surface_pressure_data: xarray.DataArray, surface_pressure_data_units: Literal['hPa', 'Pa', 'mbar'], hybrid_A_coefficients: xarray.DataArray, hybrid_B_coefficients: xarray.DataArray, vertical_output_level: list[int | float], vertical_input_dim: str, vertical_output_dim: str, vertical_output_dim_units: str, interp_method: Literal['linear', 'log', 'loglog'] = 'linear', extrapolation: bool = False, lon_dim: str = 'lon', lat_dim: str = 'lat', time_dim: str = 'time', p0_hPa=1000.0) -> xarray.DataArray

   Interpolate atmospheric data from Community Atmosphere Model (CAM) hybrid sigma-pressure coordinates to pressure levels.

   This function performs vertical interpolation of atmospheric data (typically temperature)
   from hybrid sigma-pressure coordinates to specified pressure levels using the NCL's
   ``vinth2p`` algorithm implemented in Fortran.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray <xarray.DataArray>`.
       Input 3D temperature field on hybrid levels with dimensions, e.g., (time, lev, lat, lon).
   surface_pressure_data : :py:class:`xarray.DataArray <xarray.DataArray>`.
       Surface pressure field with dimensions, e.g., (time, lat, lon).
   surface_pressure_data_units : Literal["hPa", "Pa", "mbar"]
       Units of the surface pressure data.
   hybrid_A_coefficients : :py:class:`xarray.DataArray <xarray.DataArray>`.
       Hybrid A coefficients (pressure term) for the model levels.
   hybrid_B_coefficients : :py:class:`xarray.DataArray <xarray.DataArray>`.
       Hybrid B coefficients (sigma term) for the model levels.
   vertical_output_level : list[:py:class:`int <in>`t | :py:class:`float <float>`]
       List of target pressure levels for interpolation.
   vertical_input_dim : :py:class:`str <str>`.
       Name of the vertical dimension in the input data.
   vertical_output_dim : :py:class:`str <str>`.
       Name to use for the vertical dimension in the output data.
   vertical_output_dim_units : :py:class:`str <str>`.
       Units for the output pressure levels (must be convertible to hPa).
   interp_method : Literal["linear", "log", "loglog"], optional
       Interpolation method:

       - ``"linear"``: Linear interpolation
       - ``"log"``: Logarithmic interpolation
       - ``"loglog"``: Log-log interpolation

       Default is ``"linear"``.
   extrapolation : :py:class:`bool <bool>`., optional
       Whether to extrapolate below the lowest model level when needed.
       Default is ``False``.
   lon_dim : :py:class:`str <str>`., optional
       Name of the longitude dimension. Default is ``"lon"``.
   lat_dim : :py:class:`str <str>`., optional
       Name of the latitude dimension. Default is ``"lat"``.
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name. Default is ``"time"``.
   p0_hPa : :py:class:`float <float>`., optional
       Reference pressure in hPa for hybrid level calculation. Default is ``1000.0`` hPa.

   Returns
   -------
   :py:class:`xarray.DataArray <xarray.DataArray>`
       Interpolated data on pressure levels with dimensions, e.g., (time, plev, lat, lon),
       where plev corresponds to vertical_output_level.

   .. seealso::

       - https://www.ncl.ucar.edu/Document/Functions/Built-in/vinth2p.shtml

   Examples
   --------
   >>> # Interpolate temperature to pressure levels
   >>> interp_vinth2p_dp(
   ...     data_input=temp_data,
   ...     surface_pressure_data=psfc_data,
   ...     surface_pressure_data_units="Pa",
   ...     hybrid_A_coefficients=hyam,
   ...     hybrid_B_coefficients=hybm,
   ...     vertical_output_level=[1000, 850, 700, 500, 300],
   ...     vertical_input_dim="lev",
   ...     vertical_output_dim="plev",
   ...     vertical_output_dim_units="hPa",
   ...     interp_method="log"
   ... )


.. py:function:: interp_vinth2p_ecmwf(data_input: xarray.DataArray, surface_pressure_data: xarray.DataArray, surface_pressure_data_units: Literal['hPa', 'Pa', 'mbar'], hybrid_A_coefficients: xarray.DataArray, hybrid_B_coefficients: xarray.DataArray, vertical_output_level: list[int | float], vertical_input_dim: str, vertical_output_dim: str, vertical_output_dim_units: str, variable_flag: Literal['T', 'Z', 'other'], temperature_bottom_data: Optional[xarray.DataArray] = None, surface_geopotential_data: Optional[xarray.DataArray] = None, interp_method: Literal['linear', 'log', 'loglog'] = 'linear', extrapolation: bool = True, lon_dim: str = 'lon', lat_dim: str = 'lat', time_dim: str = 'time', p0_hPa: float = 1000.0) -> xarray.DataArray

   Interpolate atmospheric data from Community Atmosphere Model (CAM) hybrid sigma-pressure
   coordinates to pressure levels using ECMWF extrapolation methods.

   This function performs vertical interpolation of atmospheric data from hybrid sigma-pressure
   coordinates to specified pressure levels using the NCL's ``vinth2p_ecmwf`` algorithm implemented
   in Fortran. It supports ECMWF-specific extrapolation for temperature ('T') and geopotential
   height ('Z') below the lowest hybrid level.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray <xarray.DataArray>`
       Input 3D field (e.g., temperature or geopotential) on hybrid levels with dimensions,
       e.g., (time, lev, lat, lon).
   surface_pressure_data : :py:class:`xarray.DataArray <xarray.DataArray>`
       Surface pressure field with dimensions, e.g., (time, lat, lon).
   surface_pressure_data_units : ``Literal["hPa", "Pa", "mbar"]``
       Units of the surface pressure data.
   hybrid_A_coefficients : :py:class:`xarray.DataArray <xarray.DataArray>`
       Hybrid A coefficients (pressure term) for the model levels.
   hybrid_B_coefficients : :py:class:`xarray.DataArray <xarray.DataArray>`
       Hybrid B coefficients (sigma term) for the model levels.
   vertical_output_level : list[:py:class:`int <int>` | :py:class:`float <float>`]
       List of target pressure levels for interpolation (in specified units).
   vertical_input_dim : :py:class:`str <str>`
       Name of the vertical dimension in the input data.
   vertical_output_dim : :py:class:`str <str>`
       Name to use for the vertical dimension in the output data.
   vertical_output_dim_units : :py:class:`str <str>`
       Units for the output pressure levels (must be convertible to hPa).
   variable_flag : ``Literal["T", "Z", "other"]``
       Indicates the type of variable being interpolated:
       - "T": Temperature (uses ECMWF extrapolation if enabled)
       - "Z": Geopotential height (uses ECMWF extrapolation if enabled)
       - "other": Any other variable (uses lowest level value for extrapolation)
   temperature_bottom_data : Optional[:py:class:`xarray.DataArray <xarray.DataArray>`], optional
       Temperature at the lowest model level (required for 'Z' extrapolation).
       Dimensions, e.g., (time, lat, lon). Default is None.
   surface_geopotential_data : Optional[:py:class:`xarray.DataArray <xarray.DataArray>`], optional
       Surface geopotential (required for 'T' or 'Z' extrapolation).
       Dimensions, e.g., (time, lat, lon). Default is None.
   interp_method : ``Literal["linear", "log", "loglog"]``, optional
       Interpolation method:

       - ``"linear"``: Linear interpolation
       - ``"log"``: Logarithmic interpolation
       - ``"loglog"``: Log-log interpolation

       Default is ``"linear"``.
   extrapolation : :py:class:`bool <bool>`, optional
       Whether to extrapolate below the lowest model level when needed.
       Default is False.
   lon_dim : :py:class:`str <str>`, optional
       Name of the longitude dimension. Default is "lon".
   lat_dim : :py:class:`str <str>`, optional
       Name of the latitude dimension. Default is "lat".
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name. Default is "time".
   p0_hPa : :py:class:`float <float>`, optional
       Reference pressure in **hPa** for hybrid level calculation. Default is ``1000.0 hPa``.

   Returns
   -------
   :py:class:`xarray.DataArray <xarray.DataArray>`
       Interpolated data on pressure levels with dimensions, e.g., (time, plev, lat, lon),
       where plev corresponds to ``vertical_output_level``.

   Notes
   -----
   - The hybrid level pressure is calculated as: :math:`P = A \cdot p_0 + B \cdot \mathrm{psfc}`
   - Output pressure levels are converted to hPa internally for calculations.
   - Missing values are converted to NaN in the output.
   - ECMWF extrapolation is applied only when ``extrapolation=True`` and variable_flag is 'T' or 'Z'.
   - For 'T' or 'Z', tbot and phis must be provided when ``extrapolation=True``.

   .. seealso::

       - https://www.ncl.ucar.edu/Document/Functions/Built-in/vinth2p_ecmwf.shtml

   Examples
   --------
   >>> # Interpolate temperature to pressure levels with ECMWF extrapolation
   >>> interp_vinth2p_ecmwf(
   ...     data_input=temp_data,
   ...     surface_pressure_data=psfc_data,
   ...     surface_pressure_data_units="Pa",
   ...     hybrid_A_coefficients=hyam,
   ...     hybrid_B_coefficients=hybm,
   ...     vertical_output_level=[1000, 850, 700, 500, 300],
   ...     vertical_input_dim="lev",
   ...     vertical_output_dim="plev",
   ...     vertical_output_dim_units="hPa",
   ...     variable_flag="T",
   ...     temperature_bottom_data=tbot_data,
   ...     surface_geopotential_data=phis_data,
   ...     interp_method="log",
   ...     extrapolation=True
   ... )


.. py:function:: interp_vintp2p_ecmwf(data_input: xarray.DataArray, pressure_data: xarray.DataArray, pressure_data_units: Literal['hPa', 'Pa', 'mbar'], surface_pressure_data: xarray.DataArray, surface_pressure_data_units: Literal['hPa', 'Pa', 'mbar'], vertical_output_level: list[int | float], vertical_input_dim: str, vertical_output_dim: str, vertical_output_dim_units: str, variable_flag: Literal['T', 'Z', 'other'], temperature_bottom_data: Optional[xarray.DataArray] = None, surface_geopotential_data: Optional[xarray.DataArray] = None, interp_method: Literal['linear', 'log', 'loglog'] = 'linear', extrapolation: bool = False, lon_dim: str = 'lon', lat_dim: str = 'lat', time_dim: str = 'time') -> xarray.DataArray

   Interpolates data at multidimensional pressure levels to constant pressure coordinates and uses an ECMWF formulation to extrapolate values below ground.

   This function performs vertical interpolation of atmospheric data from input pressure levels
   to specified output pressure levels using the NCL's `vintp2p_ecmwf` algorithm implemented
   in Fortran. It supports ECMWF-specific extrapolation for temperature ('T') and geopotential
   height ('Z') below the lowest pressure level.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray <xarray.DataArray>`
       Input 3D field (e.g., temperature or geopotential) on pressure levels with dimensions,
       e.g., (time, lev, lat, lon).
   pressure_data : :py:class:`xarray.DataArray <xarray.DataArray>`
       3D pressure field corresponding to the input data levels, with dimensions,
       e.g., (time, lev, lat, lon).
   pressure_data_units : Literal["hPa", "Pa", "mbar"]
       Units of the pressure data.
   surface_pressure_data : :py:class:`xarray.DataArray <xarray.DataArray>`
       Surface pressure field with dimensions, e.g., (time, lat, lon).
   surface_pressure_data_units : Literal["hPa", "Pa", "mbar"]
       Units of the surface pressure data.
   vertical_output_level : list[:py:class:`int <int>` | :py:class:`float <float>`]
       List of target pressure levels for interpolation (in specified units).
   vertical_input_dim : :py:class:`str <str>`
       Name of the vertical dimension in the input data.
   vertical_output_dim : :py:class:`str <str>`
       Name to use for the vertical dimension in the output data.
   vertical_output_dim_units : :py:class:`str <str>`
       Units for the output pressure levels (must be convertible to hPa).
   variable_flag : ``Literal["T", "Z", "other"]``
       Indicates the type of variable being interpolated:
       - "T": Temperature (uses ECMWF extrapolation if enabled)
       - "Z": Geopotential height (uses ECMWF extrapolation if enabled)
       - "other": Any other variable (uses lowest level value for extrapolation)
   temperature_bottom_data : Optional[:py:class:`xarray.DataArray <xarray.DataArray>`], optional
       Temperature at the lowest model level (required for 'Z' extrapolation).
       Dimensions, e.g., (time, lat, lon). Default is None.
   surface_geopotential_data : Optional[:py:class:`xarray.DataArray <xarray.DataArray>`], optional
       Surface geopotential (required for 'T' or 'Z' extrapolation).
       Dimensions, e.g., (time, lat, lon). Default is None.
   interp_method : ``Literal["linear", "log", "loglog"]``, optional
       Interpolation method:

       - "linear": Linear interpolation
       - "log": Logarithmic interpolation
       - "loglog": Log-log interpolation

       Default is "linear".
   extrapolation : :py:class:`bool <bool>`, optional
       Whether to extrapolate below the lowest pressure level when needed.
       Default is False.
   lon_dim : :py:class:`str <str>`, optional
       Name of the longitude dimension. Default is "lon".
   lat_dim : :py:class:`str <str>`, optional
       Name of the latitude dimension. Default is "lat".
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name. Default is "time".

   Returns
   -------
   :py:class:`xarray.DataArray <xarray.DataArray>`
       Interpolated data on pressure levels with dimensions, e.g., (time, plev, lat, lon),
       where plev corresponds to vertical_output_level.

   Notes
   -----
   - The input pressure levels are provided directly via ``pressure_data``.
   - Output pressure levels and surface pressure are converted to hPa internally for calculations.
   - Missing values are converted to NaN in the output.
   - ECMWF extrapolation is applied only when ``extrapolation=True`` and variable_flag is 'T' or 'Z'.
   - For 'T' or 'Z', ``temperature_bottom_data`` and ``surface_geopotential_data`` must be provided when ``extrapolation=True``.

   .. seealso::

       - https://www.ncl.ucar.edu/Document/Functions/Built-in/vintp2p_ecmwf.shtml

   Examples
   --------
   >>> # Interpolate temperature to pressure levels with ECMWF extrapolation
   >>> interp_vintp2p_ecmwf(
   ...     data=temp_data,
   ...     pressure_data=pres_data,
   ...     pressure_data_units="Pa",
   ...     surface_pressure_data=psfc_data,
   ...     surface_pressure_data_units="Pa",
   ...     vertical_output_level=[1000, 850, 700, 500, 300],
   ...     vertical_input_dim="lev",
   ...     vertical_output_dim="plev",
   ...     vertical_output_dim_units="hPa",
   ...     variable_flag="T",
   ...     temperature_bottom_data=tbot_data,
   ...     surface_geopotential_data=phis_data,
   ...     interp_method="log",
   ...     extrapolation=True
   ... )


