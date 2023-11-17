:py:mod:`easyclimate.interp.point2mesh`
=======================================

.. py:module:: easyclimate.interp.point2mesh

.. autoapi-nested-parse::

   Functions for interpolate from points to grid.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.interp.point2mesh.interp_point2mesh
   easyclimate.interp.point2mesh.interp_point2mesh_S2



.. py:function:: interp_point2mesh(data, var_name, lon_dim_name='lon', lat_dim_name='lat', point=[-9.0, 47.0], grid_x=12, grid_y=12, resolution=32.0, sigma=1.0, method='optimized_convolution', num_iter=4, min_weight=0.001)

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
       w_k(\boldsymbol{x})=\text{e}^{-\frac{1}{2\sigma^2}\left\|x-\boldsymbol{x}_k\right\|^2}

   Naive computation of Barnes interpolation leads to an algorithmic complexity of O(N x W x H), 
   where N is the number of sample points and W x H the size of the underlying grid.

   For sufficiently large n (in general in the range from 3 to 6) a good approximation of 
   Barnes interpolation with a reduced complexity O(N + W x H) can be obtained by the convolutional expression

   .. math::
       f(\boldsymbol{x})\approx \frac{ (\sum_{k=1}^{N}f_k\cdot\delta_{\boldsymbol{x}_k}) *  ( r_n^{*n[x]}(x)\cdot r_n^{*n[y]}(y) )   }{ ( \sum_{k=1}^{N} \delta_{\boldsymbol{x}_k}  ) *  (  r_{n}^{*n[x]}(x)\cdot r_{n}^{*n[y]}(y)  )   }

   where :math:`\delta` is the Dirac impulse function and :math:`r(.)` an elementary rectangular function of a specific length that depends on :math:`\sigma` and :math:`n`.

   - data : :py:class:`pandas.DataFrame<pandas.DataFrame>`
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

   - var_name: :py:class:`str<python.str>`
       The name of the data variable. This should match the one in the parameter `data`.
   - lat_dim_name: :py:class:`str<python.str>`.
       Latitude dimension name. This should match the one in the parameter `data`. By default is `lat`.
   - lon_dim_name: :py:class:`str<python.str>`.
       Longitude dimension name. This should match the one in the parameter `data`. By default is `lon`.
   - point : numpy ndarray
       A 1-dimensional array of size 2 containing the coordinates of the
       start point of the grid to be used.
   - grid_x : :py:class:`int<python.int>`.
       Length in degrees in the x-direction of the interpolated rectangular grid.
   - grid_y : :py:class:`int<python.int>`.
       Length in degrees in the y-direction of the interpolated rectangular grid.
   - resolution: float
       Grid resolution. The distance between regular grid points is the reciprocal of the value.
   - sigma : float
       The Gaussian width parameter to be used.
   - method : {'optimized_convolution', 'convolution', 'radius', 'naive'}
       Designates the Barnes interpolation method to be used. The possible
       implementations that can be chosen are 'naive' for the straightforward
       implementation (algorithm A from paper), 'radius' to consider only sample
       points within a specific radius of influence, both with an algorithmic
       complexity of O(N x W x H).
       The choice 'convolution' implements algorithm B specified in the paper
       and 'optimized_convolution' is its optimization by appending tail values
       to the rectangular kernel. The latter two algorithms reduce the complexity
       down to O(N + W x H).
       The default is 'optimized_convolution'.
   - num_iter : int, optional
       The number of performed self-convolutions of the underlying rect-kernel.
       Applies only if method is 'optimized_convolution' or 'convolution'.
       The default is 4.
   - min_weight : float, optional
       Choose radius of influence such that Gaussian weight of considered sample
       points is greater than `min_weight`.
       Applies only if method is 'radius'. Recommended values are 0.001 and less.
       The default is 0.001, which corresponds to a radius of 3.717 * sigma.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::   
       - `fast-barnes-py <https://github.com/MeteoSwiss/fast-barnes-py>`__


.. py:function:: interp_point2mesh_S2(data, var_name, lon_dim_name='lon', lat_dim_name='lat', point=[-9.0, 47.0], grid_x=12, grid_y=12, resolution=32.0, sigma=1.0, method='optimized_convolution_S2', num_iter=4, resample=True)

   Computes the Barnes interpolation for observation values `var_name` taken at sample
   points `data` using Gaussian weights for the width parameter `sigma`.

   The underlying grid embedded on the unit sphere S^2 and thus inherits the
   spherical distance measure (taken in degrees). The grid is given by the start
   point `point`, regular x-direction length `grid_x` (degree), regular y-direction length `grid_y` (degree),
   and resolution `resolution`.

   Parameters
   ----------
   - data : :py:class:`pandas.DataFrame<pandas.DataFrame>`
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

   - var_name: :py:class:`str<python.str>`
       The name of the data variable. This should match the one in the parameter `data`.
   - lat_dim_name: :py:class:`str<python.str>`.
       Latitude dimension name. This should match the one in the parameter `data`. By default is `lat`.
   - lon_dim_name: :py:class:`str<python.str>`.
       Longitude dimension name. This should match the one in the parameter `data`. By default is `lon`.
   - point : numpy ndarray
       A 1-dimensional array of size 2 containing the coordinates of the
       start point of the grid to be used.
   - grid_x : :py:class:`int<python.int>`.
       Length in degrees in the x-direction of the interpolated rectangular grid.
   - grid_y : :py:class:`int<python.int>`.
       Length in degrees in the y-direction of the interpolated rectangular grid.
   - resolution: float
       Grid resolution. The distance between regular grid points is the reciprocal of the value.
   - sigma : float
       The Gaussian width parameter to be used.
   - method : {'optimized_convolution_S2', 'naive_S2'}
       Designates the Barnes interpolation method to be used. The possible
       implementations that can be chosen are 'naive_S2' for the straightforward
       implementation (algorithm A from the paper) with an algorithmic complexity
       of O(N x W x H).
       The choice 'optimized_convolution_S2' implements the optimized algorithm B
       specified in the paper by appending tail values to the rectangular kernel.
       The latter algorithm has a reduced complexity of O(N + W x H).
       The default is 'optimized_convolution_S2'.
   - num_iter : int, optional
       The number of performed self-convolutions of the underlying rect-kernel.
       Applies only if method is 'optimized_convolution_S2'.
       The default is 4.
   - resample : bool, optional
       Specifies whether to resample Lambert grid field to lonlat grid.
       Applies only if method is 'optimized_convolution_S2'.
       The default is True.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::   
       - `fast-barnes-py <https://github.com/MeteoSwiss/fast-barnes-py>`__


