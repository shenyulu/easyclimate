easyclimate.core.diff
=====================

.. py:module:: easyclimate.core.diff

.. autoapi-nested-parse::

   The calculation of geographic finite difference



Functions
---------

.. autoapisummary::

   easyclimate.core.diff.calc_gradient
   easyclimate.core.diff.calc_dx_gradient
   easyclimate.core.diff.calc_dlon_radian_gradient
   easyclimate.core.diff.calc_dlon_degree_gradient
   easyclimate.core.diff.calc_dy_gradient
   easyclimate.core.diff.calc_dlat_radian_gradient
   easyclimate.core.diff.calc_dlat_degree_gradient
   easyclimate.core.diff.calc_dx_laplacian
   easyclimate.core.diff.calc_dy_laplacian
   easyclimate.core.diff.calc_dxdy_mixed_derivatives
   easyclimate.core.diff.calc_p_gradient
   easyclimate.core.diff.calc_time_gradient
   easyclimate.core.diff.calc_delta_pressure
   easyclimate.core.diff.calc_p_integral
   easyclimate.core.diff.calc_top2surface_integral
   easyclimate.core.diff.calc_dxdy_laplacian
   easyclimate.core.diff.calc_divergence
   easyclimate.core.diff.calc_vorticity
   easyclimate.core.diff.calc_geostrophic_wind
   easyclimate.core.diff.calc_geostrophic_wind_vorticity
   easyclimate.core.diff.calc_horizontal_water_flux
   easyclimate.core.diff.calc_vertical_water_flux
   easyclimate.core.diff.calc_water_flux_top2surface_integral
   easyclimate.core.diff.calc_divergence_watervaporflux
   easyclimate.core.diff.calc_divergence_watervaporflux_top2surface_integral
   easyclimate.core.diff.calc_u_advection
   easyclimate.core.diff.calc_v_advection
   easyclimate.core.diff.calc_p_advection
   easyclimate.core.diff.calc_shear_stretch_deform


Module Contents
---------------

.. py:function:: calc_gradient(data_input: xarray.DataArray | xarray.Dataset, dim: str, varargs: int = 1, edge_order: int = 2) -> xarray.DataArray | xarray.Dataset

   Compute the gradient along the coordinate `dim` direction.

   The gradient is computed using **second order accurate central differences** in the interior points
   and either first or second order accurate one-sides (forward or backwards) differences at the boundaries.
   The returned gradient hence has the same shape as the input array.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
        The spatio-temporal data to be calculated.
   dim : :py:class:`str <str>`.
       Dimension(s) over which to apply gradient. By default gradient is applied over the `time` dimension.
   varargs: :py:class:`list <list>` of scalar or array, optional
       Spacing between f values. Default unitary spacing for all dimensions. Spacing can be specified using:

       1. Single scalar to specify a sample distance for all dimensions.
       2. N scalars to specify a constant sample distance for each dimension. i.e. :math:`\mathrm{d}x, \mathrm{d}y, \mathrm{d}z, ...`
       3. N arrays to specify the coordinates of the values along each dimension of F.
          The length of the array must match the size of the corresponding dimension.
       4. Any combination of N scalars/arrays with the meaning of 2. and 3.

   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 2.

   Returns
   -------
   The gradient along the coordinate `dim` direction (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`numpy.gradient <numpy:numpy.gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_dx_gradient(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon', lat_dim: str = 'lat', min_dx: float = 1.0, edge_order: int = 2, R: float = 6371200.0) -> xarray.DataArray | xarray.Dataset

   Calculate the zonal gradient of the input data in physical units (meters).

   This function computes the partial derivative :math:`\partial F / \partial x`,
   where :math:`x` is the eastward distance along a parallel (in meters). It is
   the full physical zonal gradient on a sphere, given by:

   .. math::

       \frac{\partial F}{\partial x} = \frac{1}{R \cos \varphi} \cdot \frac{\partial F}{\partial \lambda}

   where :math:`R` is the Earth's radius, :math:`\varphi` is latitude, and
   :math:`\lambda` is longitude in radians. This is essential for dynamical
   calculations like advection or wave propagation in atmospheric/oceanic models.

   The computation uses finite differences along the longitude dimension:
   :math:`\partial F / \partial x = (\partial F / \partial i) / (\partial x / \partial i)`,
   where :math:`i` is the grid index. Longitude is assumed in degrees and
   converted to radians; latitude is broadcasted for the cosine factor.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated. If a Dataset, the gradient is applied to all
       data variables.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default, the gradient is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. Used to compute the cosine factor; broadcasted if necessary.
   min_dx: :py:class:`float <float>`, default: `1.0` (:math:`\mathrm{m}`).
       The minimum acceptable value of `dx` (zonal spacing in meters), below which the output
       is set to NaN to avoid numerical instabilities from very small grid spacings.
       Set to a negative value to disable this check.
   edge_order: {1, 2}, optional
       Order of the finite difference used at the boundaries. 1 uses first-order accurate
       one-sided differences; 2 uses second-order accurate one-sided differences. Default: 2.
   R: :py:class:`float <float>`, default: `6370000` (:math:`\mathrm{m}`).
       Radius of the Earth in meters (approximate mean radius). Can be adjusted for specific models.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The zonal gradient :math:`\partial F / \partial x`, with the same shape and coordinates
       as the input. Units are those of the input data divided by meters (e.g., if F is in K,
       output is K/m). Invalid regions (dx < min_dx) are NaN.

   .. seealso::
       - :py:func:`calc_gradient <calc_gradient>`
       - :py:func:`calc_dlon_radian_gradient <calc_dlon_radian_gradient>`
       - :py:func:`calc_dy_gradient <calc_dy_gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_dlon_radian_gradient(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon', edge_order: int = 2) -> xarray.DataArray | xarray.Dataset

   Calculate the gradient of the input data with respect to longitude in radians.

   This function computes the partial derivative :math:`\partial F / \partial \lambda`,
   where :math:`\lambda` is the longitude in radians. It is useful for spherical coordinate
   calculations, such as in wave activity flux (WAF) formulations, where angular gradients
   must be in radians for consistency with trigonometric functions and the Earth's radius.

   The computation uses finite differences:

   .. math::

       \frac{\partial F}{\partial \lambda} = \frac{\partial F}{\partial i } / \frac{\partial \lambda}{\partial i},

   where :math:`i` is the grid index along the longitude dimension. Longitude values are
   assumed to be in degrees initially and converted to radians for the denominator.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated. If a Dataset, the gradient is applied to all data variables.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default, the gradient is applied over the `lon` dimension.
   edge_order: {1, 2}, optional
       Order of the finite difference used at the boundaries. 1 uses first-order accurate
       one-sided differences; 2 uses second-order accurate one-sided differences. Default: 2.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The gradient :math:`\partial F / \partial \lambda` (in radians), with the same shape
       and coordinates as the input. Units are inherited from the input data divided by radians
       (e.g., if F is in K, output is K/rad).


.. py:function:: calc_dlon_degree_gradient(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon', edge_order: int = 2) -> xarray.DataArray | xarray.Dataset

   Calculate the gradient of the input data with respect to longitude in degrees.

   This function computes the partial derivative :math:`\partial F / \partial \lambda`,
   where :math:`\lambda` is the longitude in degrees. It is suitable for general-purpose
   gradient calculations where degree units are preferred for interpretability.

   The computation uses finite differences:

   .. math::

       \frac{\partial F}{\partial \lambda} = \frac{\partial F}{\partial i} / \frac{\partial \lambda}{\partial i},

   where :math:`i` is the grid index along the longitude dimension. Longitude values remain
   in degrees for the denominator.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated. If a Dataset, the gradient is applied to all
       data variables.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default, the gradient is applied over the `lon` dimension.
   edge_order: {1, 2}, optional
       Order of the finite difference used at the boundaries. 1 uses first-order accurate
       one-sided differences; 2 uses second-order accurate one-sided differences. Default: 2.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The gradient :math:`\partial F / \partial \lambda` (in degrees), with the same shape
       and coordinates as the input. Units are inherited from the input data divided by degrees
       (e.g., if F is in K, output is K/deg).

   .. seealso::
       - :py:func:`calc_gradient <calc_gradient>`
       - :py:func:`calc_dlon_radian_gradient <calc_dlon_radian_gradient>`
       - :py:func:`calc_dlat_degree_gradient <calc_dlat_degree_gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_dy_gradient(data_input: xarray.DataArray | xarray.Dataset, lat_dim: str = 'lat', min_dy: float = 1.0, edge_order: int = 2, R: float = 6371200.0) -> xarray.DataArray | xarray.Dataset

   Calculate the meridional gradient of the input data in physical units (meters).

   This function computes the partial derivative :math:`\partial F / \partial y`,
   where :math:`y` is the northward distance along a meridian (in meters). It is
   the full physical meridional gradient on a sphere, given by:

   .. math::

       \frac{\partial F}{\partial y} = \frac{1}{R} \cdot \frac{\partial F}{\partial \varphi}

   where :math:`R` is the Earth's radius and :math:`\varphi` is latitude in radians.
   This is essential for dynamical calculations like advection or vorticity in models.

   The computation uses finite differences along the latitude dimension:
   :math:`\partial F / \partial y = (\partial F / \partial j) / (\partial y / \partial j)`,
   where :math:`j` is the grid index. Latitude is assumed in degrees and
   converted to radians.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated. If a Dataset, the gradient is applied to all
       data variables.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default, the gradient is applied over the `lat` dimension.
   min_dy: :py:class:`float <float>`, default: `1.0`.
       The minimum acceptable value of `dy` (meridional spacing in meters), below which the output
       is set to NaN to avoid numerical instabilities from very small grid spacings.
       Set to a negative value to disable this check. Unit: meters.
   edge_order: {1, 2}, optional
       Order of the finite difference used at the boundaries. 1 uses first-order accurate
       one-sided differences; 2 uses second-order accurate one-sided differences. Default: 2.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth in meters (approximate mean radius). Can be adjusted for specific models.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The meridional gradient :math:`\partial F / \partial y`, with the same shape and coordinates
       as the input. Units are those of the input data divided by meters (e.g., if F is in K,
       output is K/m). Invalid regions (dy < min_dy) are NaN.

   .. seealso::
       - :py:func:`calc_gradient <calc_gradient>`
       - :py:func:`calc_dlat_radian_gradient <calc_dlat_radian_gradient>`
       - :py:func:`calc_dx_gradient <calc_dx_gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_dlat_radian_gradient(data_input: xarray.DataArray | xarray.Dataset, lat_dim: str = 'lat', edge_order: int = 2) -> xarray.DataArray | xarray.Dataset

   Calculate the gradient of the input data with respect to latitude in radians.

   This function computes the partial derivative :math:`\partial F / \partial \phi`,
   where :math:`\phi` is the latitude in radians. It is useful for spherical coordinate
   calculations, such as in wave activity flux (WAF) formulations.

   The computation uses finite differences:

   .. math::

       \frac{\partial F}{\partial \phi} = \frac{\partial F}{\partial j} / \frac{\partial \phi}{\partial j},

   where :math:`j` is the grid index along the latitude dimension. Latitude values are
   assumed to be in degrees initially and converted to radians for the denominator.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated. If a Dataset, the gradient is applied to all
       data variables.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default, the gradient is applied over the `lat` dimension.
   edge_order: {1, 2}, optional
       Order of the finite difference used at the boundaries. 1 uses first-order accurate
       one-sided differences; 2 uses second-order accurate one-sided differences. Default: 2.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The gradient :math:`\partial F / \partial \phi` (in radians), with the same shape
       and coordinates as the input. Units are inherited from the input data divided by radians
       (e.g., if F is in K, output is K/rad).

   .. seealso::
       - :py:func:`calc_gradient <calc_gradient>`
       - :py:func:`calc_dlat_degree_gradient <calc_dlat_degree_gradient>`
       - :py:func:`calc_dlon_radian_gradient <calc_dlon_radian_gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_dlat_degree_gradient(data_input: xarray.DataArray | xarray.Dataset, lat_dim: str = 'lat', edge_order: int = 2) -> xarray.DataArray | xarray.Dataset

   Calculate the gradient of the input data with respect to latitude in degrees.

   This function computes the partial derivative :math:`\partial F / \partial \phi`,
   where :math:`\phi` is the latitude in degrees. It is suitable for general-purpose
   gradient calculations where degree units are preferred.

   The computation uses finite differences:

   .. math::

       \frac{\partial F}{\partial \phi} = \frac{\partial F}{\partial j} / \frac{\partial \phi}{\partial j},

   where :math:`j` is the grid index along the latitude dimension. Latitude values remain
   in degrees for the denominator.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated. If a Dataset, the gradient is applied to all
       data variables.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default, the gradient is applied over the `lat` dimension.
   edge_order: {1, 2}, optional
       Order of the finite difference used at the boundaries. 1 uses first-order accurate
       one-sided differences; 2 uses second-order accurate one-sided differences. Default: 2.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The gradient :math:`\partial F / \partial \phi` (in degrees), with the same shape
       and coordinates as the input. Units are inherited from the input data divided by degrees
       (e.g., if F is in K, output is K/deg).

   .. seealso::
       - :py:func:`calc_gradient <calc_gradient>`
       - :py:func:`calc_dlat_radian_gradient <calc_dlat_radian_gradient>`
       - :py:func:`calc_dlon_degree_gradient <calc_dlon_degree_gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_dx_laplacian(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon', lat_dim: str = 'lat', min_dx2: float = 1000000000.0, edge_order: int = 2, R: float = 6371200.0) -> xarray.DataArray | xarray.Dataset

   Calculation of the second-order partial derivative term (Laplace term) along longitude.

   .. math::
       \frac{\partial^2 F}{\partial x^2} = \frac{1}{(R \cos\varphi)^2} \cdot \frac{\partial^2 F}{\partial \lambda^2}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   min_dx2: :py:class:`float <float>`, default: `1e9`.
       The minimum acceptable value of :math:`(\mathrm{d}x)^2`, below which parts will set `nan` to avoid large computational errors.
       The unit is m. You can set it to a negative value in order to remove this benefit.
   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The second-order partial derivative term (Laplace term) along longitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_dy_laplacian(data_input: xarray.DataArray | xarray.Dataset, lat_dim: str = 'lat', min_dy2: float = 1.0, edge_order: int = 2, R: float = 6371200.0) -> xarray.DataArray | xarray.Dataset

   Calculation of the second-order partial derivative term (Laplace term) along latitude.

   .. math::
       \frac{\partial^2 F}{\partial y^2} = \frac{1}{R^2} \cdot \frac{\partial^2 F}{\partial \varphi^2}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   min_dy2: :py:class:`float <float>`, default: `1.0`.
       The minimum acceptable value of :math:`(\mathrm{d}y)^2`, below which parts will set `nan` to avoid large computational errors.
       The unit is m. You can set it to a negative value in order to remove this benefit.
   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The second-order partial derivative term (Laplace term) along latitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_dxdy_mixed_derivatives(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon', lat_dim: str = 'lat', min_dxdy: float = 10000000000.0, edge_order: int = 2, R: float = 6371200.0) -> xarray.DataArray | xarray.Dataset

   Calculation of second-order mixed partial derivative terms along longitude and latitude.

   .. math::
       \frac{\partial^2 F}{\partial x \partial y} = \frac{1}{R^2 \cos\varphi} \cdot \frac{\partial^2 F}{\partial \lambda \partial \varphi}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   min_dxdy: :py:class:`float <float>`, default: `1e10`.
       The minimum acceptable value of :math:`\mathrm{d}x\mathrm{d}y`, below which parts will set `nan` to avoid large computational errors.
       The unit is m. You can set it to a negative value in order to remove this benefit.
   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The second-order mixed partial derivative terms along longitude and latitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_p_gradient(data_input: xarray.DataArray, vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar']) -> xarray.DataArray

   Calculate the gradient along the barometric pressure direction in the p-coordinate system.

   .. math::
       \frac{\partial F}{\partial p}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.

   Returns
   -------
   The gradient along the barometric pressure direction in the p-coordinate system (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>`


.. py:function:: calc_time_gradient(data_input: xarray.DataArray, time_units: str, time_dim: str = 'time') -> xarray.DataArray

   Calculate the gradient along the time direction.

   .. math::
       \frac{\partial F}{\partial t}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   time_units: :py:class:`str <str>`.
       The unit corresponding to the time dimension value. Optional values are `seconds`, `months`, `years` and so on.
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name.

   Returns
   -------
   The gradient along the time direction (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. caution:: The units for partial derivative of `time` are :math:`\mathrm{s^{-1}}`.

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>`


.. py:function:: calc_delta_pressure(data_input: xarray.DataArray, surface_pressure_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], surface_pressure_data_units: Literal['hPa', 'Pa', 'mbar']) -> xarray.DataArray

   Calculates the pressure layer thickness (delta pressure) of a constant
   pressure level coordinate system.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The spatio-temporal data to be calculated.
   surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Mean surface sea level pressure.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   surface_pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.

   Returns
   -------
   The pressure layer thickness (delta pressure) of a constant pressure level coordinate system (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       - :py:func:`geocat.comp.meteorology.delta_pressure <geocat-comp:geocat.comp.meteorology.delta_pressure>`
       - `dpres_plevel - NCL <https://www.ncl.ucar.edu/Document/Functions/Built-in/dpres_plevel.shtml>`__

   Examples
   --------
   The results in :py:func:`geocat.comp.meteorology.delta_pressure <geocat-comp:geocat.comp.meteorology.delta_pressure>`:

   >>> from geocat.comp.meteorology import delta_pressure
   >>> dp = delta_pressure(
   ...     pressure_lev= np.array([1000.,925.,850.,700.,600.,500., 400.,300.,250.,200.,150.,100., 70.,50.,30.,20.,10.]),
   ...     surface_pressure = np.array([1013]),
   ... )
   >>> print(dp)
   [[ 50.5  75.  112.5 125.  100.  100.  100.   75.   50.   50.   50.   40.
      25.   20.   15.   10.    5. ]]

   For comparison, the results in :py:func:`easyclimate.calc_delta_pressure <calc_delta_pressure>`:

   >>> temp_sample = xr.DataArray(
   ...     np.array([[292.,285.,283.,277.,270.,260., 250.,235.,225.,215.,207.,207., 213.,220.,225.,228.,230.]]),
   ...     dims = ("lat", "plev"),
   ...     coords = {"plev": np.array([1000.,925.,850.,700.,600.,500., 400.,300.,250.,200.,150.,100., 70.,50.,30.,20.,10.]),
   ...             "lat": np.array([0])}
   ... )
   >>> dp = ecl.calc_delta_pressure(
   ...     data_input = temp_sample,
   ...     surface_pressure_data = xr.DataArray([1013], dims = "lat"),
   ...     vertical_dim = "plev",
   ...     surface_pressure_data_units = "Pa",
   ...     vertical_dim_units = "Pa",
   ... ).transpose("lat", "plev")
   >>> print(dp)
   <xarray.DataArray 'plev' (lat: 1, plev: 17)> Size: 136B
   array([[ 50.5,  75. , 112.5, 125. , 100. , 100. , 100. ,  75. ,  50. ,
           50. ,  50. ,  40. ,  25. ,  20. ,  15. ,  10. ,   5. ]])
   Coordinates:
   * lat      (lat) int64 8B 0
   * plev     (plev) float64 136B 1e+03 925.0 850.0 700.0 ... 50.0 30.0 20.0 10.0


.. py:function:: calc_p_integral(data_input: xarray.DataArray, vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], normalize: bool = False) -> xarray.DataArray

   Calculate the vertical integral along the barometric pressure direction in the p-coordinate system.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The spatio-temporal data to be calculated.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   normalize: :py:class:`bool<bool>`, default: `True`.
       Whether or not the integral results are averaged over the entire layer.

   Returns
   -------
   The vertical integral along the barometric pressure direction in the p-coordinate system (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. attention::
       This method ignores the effect of topography, so it applies to altitudes **above 900hPa** and is **NOT applicable to the Tibetan Plateau region**.
       For a fully accurate vertical integration, please use the :py:func:`calc_top2surface_integral <calc_top2surface_integral>` function to calculate,
       but the speed of the calculation is slightly slowed down.


.. py:function:: calc_top2surface_integral(data_input: xarray.DataArray, surface_pressure_data: xarray.DataArray, vertical_dim: str, surface_pressure_data_units: Literal['hPa', 'Pa', 'mbar'], vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], method: Literal['Boer1982', 'Trenberth1991', 'vibeta-ncl'] = 'vibeta-ncl', normalize: bool = False) -> xarray.DataArray

   Calculate the vertical integral in the p-coordinate system from the ground to the zenith along the barometric pressure direction.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The spatio-temporal data to be calculated.
   surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Surface level pressure.

   .. warning::

       This parameter only accepts local pressure (slp) and must **NOT** be substituted
       with mean sea level pressure (msl). There is a fundamental difference in their physical meaning and numerical
       characteristics—local pressure reflects atmospheric pressure at the actual elevation of local area,
       while mean sea level pressure is a theoretical value adjusted to sea level.

   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   surface_pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   method: {"Boer1982", "Trenberth1991", "vibeta-ncl"}, default: `vibeta-ncl`.
       vertical integration method. Optional values are ``Boer1982``, ``Trenberth1991`` and ``vibeta-ncl``.

       .. note::
           The trapezoidal rule of integration is exactly equivalent to

           .. math::
               I = \sum_{j=1,2J-1,2} (\beta M)_j \Delta p_j,

           where Kevin E. Trenberth (1991) define

           .. math::
               \beta_j = \left\lbrace
               \begin{array}{ll}
               1, & \mathrm{if} \ p_{j-1} < p_s,\\
               0, & \mathrm{if} \ p_{j+1} > p_s ,\\
               \frac{p_s - p_{j+1}}{p_{j-1} - p_{j+1}}, & \mathrm{if}  \ p_{j-1} > p_s > p_{j+1}.
               \end{array}
               \right.

           While G. J. Boer (1982) define :math:`\beta = 0, 1` only.

   normalize: :py:class:`bool<bool>`, default: `True`.
       Whether or not the integral results are averaged over the entire layer.

   Returns
   -------
   The vertical integral in the p-coordinate system from the ground to the zenith along the barometric pressure direction. (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - `Boer, G. J., 1982: Diagnostic Equations in Isobaric Coordinates. Mon. Wea. Rev., 110, 1801–1820, <https://doi.org/10.1175/1520-0493(1982)110%3C1801:DEIIC%3E2.0.CO;2>`__
   - `Trenberth, K. E., 1991: Climate Diagnostics from Global Analyses: Conservation of Mass in ECMWF Analyses. J. Climate, 4, 707–722, <https://doi.org/10.1175/1520-0442(1991)004%3C0707:CDFGAC%3E2.0.CO;2>`__

   .. seealso::
       - `vibeta - NCL <https://www.ncl.ucar.edu/Document/Functions/Built-in/vibeta.shtml>`__
       - `dpres_plevel - NCL <https://www.ncl.ucar.edu/Document/Functions/Built-in/dpres_plevel.shtml>`__


.. py:function:: calc_dxdy_laplacian(data_input: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', R: float = 6371200.0, spherical_coord: bool = True) -> xarray.DataArray

   Calculate the horizontal Laplace term.

   rectangular coordinates

   .. math::
       \nabla^2 F = \frac{\partial^2 F}{\partial x^2} + \frac{\partial^2 F}{\partial y^2}

   Spherical coordinates

   .. math::
       \nabla^2 F = \frac{\partial^2 F}{\partial x^2} + \frac{\partial^2 F}{\partial y^2} - \frac{1}{R} \frac{\partial F}{\partial y} \tan \varphi

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The spatio-temporal data to be calculated.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.
   spherical_coord: :py:class:`bool <bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates.

   Returns
   -------
   The horizontal Laplace term. (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_divergence(u_data: xarray.DataArray, v_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', R: float = 6371200.0, spherical_coord=True, cyclic_boundary: bool = False, method: Literal['easyclimate', 'uv2dv_cfd-ncl'] = 'uv2dv_cfd-ncl') -> xarray.DataArray

   Calculate the horizontal divergence term.

   rectangular coordinates

   .. math::
       \mathrm{D} = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y}

   Spherical coordinates

   .. math::
       \mathrm{D} = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} - \frac{v}{R} \tan \varphi

   Parameters
   ----------
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   R: :py:class:`float <float>`, default: `6371200.0`.
       Radius of the Earth.
   spherical_coord: :py:class:`bool<bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates. The parameter is applicable only when ``method = easyclimate``.
   cyclic_boundary: :py:class:`bool <bool>`, default: `False`.
       If True, assume cyclic (periodic) boundaries in longitude. The parameter is applicable only when ``method = ddvfidf-ncl``.
   method: {"easyclimate", "ddvfidf-ncl"}, default: `ddvfidf-ncl`.
       The method to calculate horizontal divergence term. Optional values are ``easyclimate`` and ``ddvfidf-ncl``.

   Returns
   -------
   The horizontal divergence term. (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::

       - https://www.ncl.ucar.edu/Document/Functions/Built-in/uv2dv_cfd.shtml
       - Howard B. Bluestein. (1992). Synoptic-Dynamic Meteorology in Midlatitudes: Principles of Kinematics and Dynamics, Vol. 1. p113-114

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_vorticity(u_data: xarray.DataArray, v_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', R: float = 6371200.0, spherical_coord: bool = True, cyclic_boundary: bool = False, method: Literal['easyclimate', 'uv2vr_cfd-ncl'] = 'uv2vr_cfd-ncl') -> xarray.DataArray

   Calculate the horizontal relative vorticity term.

   rectangular coordinates

   .. math::
       \zeta = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}

   Spherical coordinates

   .. math::
       \zeta = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y} + \frac{u}{R} \tan \varphi

   Parameters
   ----------
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.
   spherical_coord: :py:class:`bool<bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates.
   cyclic_boundary: :py:class:`bool <bool>`, default: `False`.
       If True, assume cyclic (periodic) boundaries in longitude. The parameter is applicable only when ``method = uv2vr_cfd-ncl``.
   method: {"easyclimate", "uv2vr_cfd-ncl"}, default: `uv2vr_cfd-ncl`.
       The method to calculate horizontal divergence term. Optional values are ``easyclimate`` and ``uv2vr_cfd-ncl``.

   Returns
   -------
   The horizontal relative vorticity term. (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::

       - https://www.ncl.ucar.edu/Document/Functions/Built-in/uv2vr_cfd.shtml
       - Howard B. Bluestein. (1992). Synoptic-Dynamic Meteorology in Midlatitudes: Principles of Kinematics and Dynamics, Vol. 1. p113-114

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_geostrophic_wind(z_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', omega: float = 7.292e-05, g: float = 9.8, R: float = 6371200.0) -> xarray.DataArray

   Calculate the geostrophic wind.

   .. math::
       u_g = - \frac{g}{f} \frac{\partial H}{\partial y}

   .. math::
       v_g = \frac{g}{f} \frac{\partial H}{\partial x}

   Parameters
   ----------
   z_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric geopotential height.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   omega: :py:class:`float <float>`, default: `7.292e-5`.
       The angular speed of the earth.
   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The geostrophic wind term. (:py:class:`xarray.DataArray<xarray.DataArray>`).
       - ug
       - vg

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_geostrophic_wind_vorticity(z_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', spherical_coord: bool = True, omega: float = 7.292e-05, g: float = 9.8, R: float = 6371200.0, cyclic_boundary: bool = False, method: Literal['easyclimate', 'uv2vr_cfd-ncl'] = 'uv2vr_cfd-ncl') -> xarray.DataArray

   Calculate the geostrophic vorticity.

   Rectangular coordinates

   .. math::
       \zeta_g = \frac{\partial v_g}{\partial x} - \frac{\partial u_g}{\partial y}

   Spherical coordinates

   .. math::
       \zeta_g = \frac{\partial v_g}{\partial x} - \frac{\partial u_g}{\partial y} + \frac{u_g}{R} \tan \varphi

   Parameters
   ----------
   z_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric geopotential height.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   spherical_coord: :py:class:`bool<bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates.
   omega: :py:class:`float <float>`, default: `7.292e-5`.
       The angular speed of the earth.
   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.
   cyclic_boundary: :py:class:`bool <bool>`, default: `False`.
       If True, assume cyclic (periodic) boundaries in longitude. The parameter is applicable only when ``method = uv2vr_cfd-ncl``.
   method: {"easyclimate", "uv2vr_cfd-ncl"}, default: `uv2vr_cfd-ncl`.
       The method to calculate horizontal divergence term. Optional values are ``easyclimate`` and ``uv2vr_cfd-ncl``.

   Returns
   -------
   The geostrophic vorticity term. (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_horizontal_water_flux(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, g: float = 9.8) -> xarray.Dataset

   Calculate horizontal water vapor flux at each vertical level.

   .. math::
       \frac{1}{g} q \mathbf{V} = \frac{1}{g} (u q\ \mathbf{i} + vq\ \mathbf{j})

   Parameters
   ----------
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.

   Returns
   -------
   The water vapor flux. (:py:class:`xarray.Dataset<xarray.Dataset>`).

   - :math:`qu`: zonal water vapor flux.
   - :math:`qv`: meridional water vapor flux.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_vertical_water_flux(specific_humidity_data: xarray.DataArray, omega_data: xarray.DataArray, g: float = 9.8) -> xarray.DataArray

   Calculate vertical water vapor flux.

   .. math::
       -\omega \frac{q}{g}

   Parameters
   ----------
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   omega_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The vertical velocity data (:math:`\frac{\mathrm{d} p}{\mathrm{d} t}`).
   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.

   Returns
   -------
   The vertical water flux. (:py:class:`xarray.DataArray <xarray.DataArray>`).

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_water_flux_top2surface_integral(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, surface_pressure_data: xarray.DataArray, surface_pressure_data_units: Literal['hPa', 'Pa', 'mbar'], specific_humidity_data_units: Literal['kg/kg', 'g/kg', 'g/g'], vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], method: Literal['Boer1982', 'Trenberth1991', 'vibeta-ncl'] = 'vibeta-ncl', g: float = 9.8) -> xarray.DataArray

   Calculate the water vapor flux across the vertical level.

   .. math::

       \frac{1}{g} \int_0^{p_s} (q\mathbf{v}),dp

   Parameters
   ----------
   specific_humidity: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Mean surface sea level pressure.
   surface_pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
   specific_humidity_data_units: :py:class:`str <str>`.
       The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   method: :py:class:`str <str>`, default: `'Trenberth1991'`.
       Vertical integration method. Optional values are `Boer1982`, `'Trenberth1991'`.
   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.

   Returns
   -------
   The water vapor flux. (:py:class:`xarray.Dataset<xarray.Dataset>`, :math:`\mathrm{kg \cdot m^-1 \cdot s^-1 }`).

   - :math:`qu`: zonal water vapor flux.
   - :math:`qv`: meridional water vapor flux.

   .. seealso::
       :py:func:`calc_top2surface_integral <calc_top2surface_integral>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_divergence_watervaporflux(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, specific_humidity_data_units: Literal['kg/kg', 'g/kg', 'g/g'], spherical_coord: bool = True, cyclic_boundary: bool = False, method: Literal['easyclimate', 'uv2dv_cfd-ncl'] = 'uv2dv_cfd-ncl', lon_dim: str = 'lon', lat_dim: str = 'lat', g: float = 9.8, R: float = 6371200.0) -> xarray.DataArray

   Calculate water vapor flux divergence at each vertical level.

   .. math::
       \nabla \left( \frac{1}{g} q \mathbf{V} \right) = \frac{1}{g} \nabla \cdot \left( q \mathbf{V} \right)


   Parameters
   ----------
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   specific_humidity_data_units: :py:class:`str <str>`.
       The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
   spherical_coord: :py:class:`bool<bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates. The parameter is applicable only when ``method = easyclimate``.
   cyclic_boundary: :py:class:`bool <bool>`, default: `False`.
       If True, assume cyclic (periodic) boundaries in longitude. The parameter is applicable only when ``method = ddvfidf-ncl``.
   method: {"easyclimate", "ddvfidf-ncl"}, default: `ddvfidf-ncl`.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The water vapor flux divergence. (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_divergence_watervaporflux_top2surface_integral(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, surface_pressure_data: xarray.DataArray, vertical_dim: str, specific_humidity_data_units: Literal['kg/kg', 'g/kg', 'g/g'], surface_pressure_data_units: Literal['hPa', 'Pa', 'mbar'], vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], spherical_coord: bool = True, cyclic_boundary: bool = False, lon_dim: str = 'lon', lat_dim: str = 'lat', integral_method: Literal['Boer1982', 'Trenberth1991', 'vibeta-ncl'] = 'vibeta-ncl', div_method: Literal['easyclimate', 'uv2dv_cfd-ncl'] = 'uv2dv_cfd-ncl', g: float = 9.8, R: float = 6371200.0) -> xarray.DataArray

   Calculate water vapor flux divergence across the vertical level.

   .. math::

       \nabla \cdot \frac{1}{g} \int_0^{p_s} (q\mathbf{v}),dp

   Parameters
   ----------
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Mean surface sea level pressure.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   specific_humidity_data_units: :py:class:`str <str>`.
       The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
   surface_pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   spherical_coord: :py:class:`bool<bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates. The parameter is applicable only when ``method = easyclimate``.
   cyclic_boundary: :py:class:`bool <bool>`, default: `False`.
       If True, assume cyclic (periodic) boundaries in longitude. The parameter is applicable only when ``method = ddvfidf-ncl``.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   integral_method: {"Boer1982", "Trenberth1991", "vibeta-ncl"}, default: `vibeta-ncl`.
       The vertical integration method. Optional values are ``Boer1982``, ``Trenberth1991`` and ``vibeta-ncl``.

       .. note::
           The trapezoidal rule of integration is exactly equivalent to

           .. math::
               I = \sum_{j=1,2J-1,2} (\beta M)_j \Delta p_j,

           where Kevin E. Trenberth (1991) define

           .. math::
               \beta_j = \left\lbrace
               \begin{array}{ll}
               1, & \mathrm{if} \ p_{j-1} < p_s,\\
               0, & \mathrm{if} \ p_{j+1} > p_s ,\\
               \frac{p_s - p_{j+1}}{p_{j-1} - p_{j+1}}, & \mathrm{if}  \ p_{j-1} > p_s > p_{j+1}.
               \end{array}
               \right.

           While G. J. Boer (1982) define :math:`\beta = 0, 1` only.

   div_method: {"easyclimate", "ddvfidf-ncl"}, default: `ddvfidf-ncl`.
       The method to calculate horizontal divergence term. Optional values are ``easyclimate`` and ``ddvfidf-ncl``.

   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The water vapor flux divergence. (:py:class:`xarray.DataArray<xarray.DataArray>`, :math:`\mathrm{kg \cdot m^-2 \cdot s^-1 }`).

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_u_advection(u_data: xarray.DataArray, temper_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', min_dx: float = 1.0, edge_order: int = 2, R: float = 6371200.0) -> xarray.DataArray

   Calculate zonal temperature advection at each vertical level.

   .. math::
       -u \frac{\partial T}{\partial x}

   Parameters
   ----------
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   min_dx: :py:class:`float <float>`, default: `1.0`.
       The minimum acceptable value of `dx`, below which parts will set `nan` to avoid large computational errors.
       The unit is m. You can set it to a negative value in order to remove this benefit.
   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The zonal temperature advection. (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_v_advection(v_data: xarray.DataArray, temper_data: xarray.DataArray, lat_dim: str = 'lat', min_dy: float = 1.0, edge_order: int = 2, R: float = 6371200.0) -> xarray.DataArray

   Calculate meridional temperature advection at each vertical level.

   .. math::
       -v \frac{\partial T}{\partial y}

   Parameters
   ----------
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The meridional temperature advection. (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_geographic_finite_difference.py


.. py:function:: calc_p_advection(omega_data: xarray.DataArray, temper_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar']) -> xarray.DataArray

   Calculate vertical temperature transport at each vertical level.

   .. math::
       -\omega \frac{\partial T}{\partial p}

   Parameters
   ----------
   omega: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The vertical velocity data (:math:`\frac{\mathrm{d} p}{\mathrm{d} t}`).
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.

   Returns
   -------
   The vertical temperature transport. (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_shear_stretch_deform(u_data: xarray.DataArray, v_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', edge_order: {1, 2} = 2, R: float = 6371200.0)

   - https://www.ncl.ucar.edu/Document/Functions/Contributed/shear_stretch_deform_cfd.shtml
   - Spensberger, C., & Spengler, T. (2014). A New Look at Deformation as a Diagnostic for Large-Scale Flow. Journal of the Atmospheric Sciences, 71(11), 4221-4234. https://doi.org/10.1175/JAS-D-14-0108.1


