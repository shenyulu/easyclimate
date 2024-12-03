easyclimate.core.diff
=====================

.. py:module:: easyclimate.core.diff

.. autoapi-nested-parse::

   The calculation of geographic finite difference



Functions
---------

.. autoapisummary::

   easyclimate.core.diff.calc_gradient
   easyclimate.core.diff.calc_lon_gradient
   easyclimate.core.diff.calc_lat_gradient
   easyclimate.core.diff.calc_lon_laplacian
   easyclimate.core.diff.calc_lat_laplacian
   easyclimate.core.diff.calc_lon_lat_mixed_derivatives
   easyclimate.core.diff.calc_p_gradient
   easyclimate.core.diff.calc_time_gradient
   easyclimate.core.diff.calc_delta_pressure
   easyclimate.core.diff.calc_p_integral
   easyclimate.core.diff.calc_top2surface_integral
   easyclimate.core.diff.calc_laplacian
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
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.

   Returns
   -------
   The gradient along the coordinate `dim` direction (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`numpy.gradient <numpy:numpy.gradient>`


.. py:function:: calc_lon_gradient(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon', lat_dim: str = 'lat', min_dx: float = 1.0, edge_order: int = 2, R: float = 6370000) -> xarray.DataArray | xarray.Dataset

   Calculate the gradient along the longitude.

   .. math::
       \frac{\partial F}{\partial x} = \frac{1}{R \cos\varphi} \cdot \frac{\partial F}{\partial \lambda}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   min_dx: :py:class:`float <float>`, default: `1.0`.
       The minimum acceptable value of `dx`, below which parts will set `nan` to avoid large computational errors.
       The unit is m. You can set it to a negative value in order to remove this benefit.
   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 2.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The gradient along the longitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>`


.. py:function:: calc_lat_gradient(data_input: xarray.DataArray | xarray.Dataset, lat_dim: str = 'lat', min_dy: float = 1.0, edge_order: int = 2, R: float = 6370000) -> xarray.DataArray | xarray.Dataset

   Calculate the gradient along the latitude.

   .. math::
       \frac{\partial F}{\partial y} = \frac{1}{R} \cdot \frac{\partial F}{\partial \varphi}

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   min_dy: :py:class:`float <float>`, default: `1.0`.
       The minimum acceptable value of `dy`, below which parts will set `nan` to avoid large computational errors.
       The unit is m. You can set it to a negative value in order to remove this benefit.
   edge_order: {1, 2}, optional
       Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The gradient along the latitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       :py:func:`calc_gradient <calc_gradient>`


.. py:function:: calc_lon_laplacian(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon', lat_dim: str = 'lat', min_dx2: float = 1000000000.0, edge_order: int = 2, R: float = 6370000) -> xarray.DataArray | xarray.Dataset

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


.. py:function:: calc_lat_laplacian(data_input: xarray.DataArray | xarray.Dataset, lat_dim: str = 'lat', min_dy2: float = 1.0, edge_order: int = 2, R: float = 6370000) -> xarray.DataArray | xarray.Dataset

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


.. py:function:: calc_lon_lat_mixed_derivatives(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon', lat_dim: str = 'lat', min_dxdy: float = 10000000000.0, edge_order: int = 2, R: float = 6370000) -> xarray.DataArray | xarray.Dataset

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


.. py:function:: calc_p_gradient(data_input: xarray.DataArray, vertical_dim: str, vertical_dim_units: str) -> xarray.DataArray

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


.. py:function:: calc_delta_pressure(data_input: xarray.DataArray, surface_pressure_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: str, surface_pressure_data_units: str) -> xarray.DataArray

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


.. py:function:: calc_p_integral(data_input: xarray.DataArray, vertical_dim: str, normalize: bool = True) -> xarray.DataArray

   Calculate the vertical integral along the barometric pressure direction in the p-coordinate system.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The spatio-temporal data to be calculated.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   normalize: :py:class:`bool<bool>`, default: `True`.
       Whether or not the integral results are averaged over the entire layer.

   Returns
   -------
   The vertical integral along the barometric pressure direction in the p-coordinate system (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

   .. attention::
       This method ignores the effect of topography, so it applies to altitudes **above 900hPa** and is **NOT applicable to the Tibetan Plateau region**.
       For a fully accurate vertical integration, please use the :py:func:`calc_top2surface_integral <calc_top2surface_integral>` function to calculate,
       but the speed of the calculation is slightly slowed down.


.. py:function:: calc_top2surface_integral(data_input: xarray.DataArray, surface_pressure_data: xarray.DataArray, vertical_dim: str, surface_pressure_data_units: str, vertical_dim_units: str, method: str = 'Trenberth-vibeta', normalize: bool = True) -> xarray.DataArray

   Calculate the vertical integral in the p-coordinate system from the ground to the zenith along the barometric pressure direction.

   Parameters
   ----------
   data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The spatio-temporal data to be calculated.
   surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Mean surface sea level pressure.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   surface_pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   method: :py:class:`str <str>`, default: `'Trenberth-vibeta'`.
       vertical integration method. Optional values are `Boer-vibeta`, `'Trenberth-vibeta'`.

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


.. py:function:: calc_laplacian(data_input: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', R: float = 6370000, spherical_coord: bool = True) -> xarray.DataArray

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


.. py:function:: calc_divergence(u_data: xarray.DataArray, v_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', R: float = 6370000, spherical_coord=True) -> xarray.DataArray

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
   R: :py:class:`float <float>`, default: `6370000`.
       Radius of the Earth.
   spherical_coord: :py:class:`bool<bool>`, default: `True`.
       Whether or not to compute the horizontal Laplace term in spherical coordinates.

   Returns
   -------
   The horizontal divergence term. (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_vorticity(u_data: xarray.DataArray, v_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', R: float = 6370000, spherical_coord: bool = True) -> xarray.DataArray

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

   Returns
   -------
   The horizontal relative vorticity term. (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_geostrophic_wind(z_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', omega: float = 7.292e-05, g: float = 9.8, R: float = 6370000) -> xarray.DataArray

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


.. py:function:: calc_geostrophic_wind_vorticity(z_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', spherical_coord: bool = True, omega: float = 7.292e-05, g: float = 9.8, R: float = 6370000) -> xarray.DataArray

   Calculate the geostrophic vorticity.

   rectangular coordinates

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


.. py:function:: calc_water_flux_top2surface_integral(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, surface_pressure_data: xarray.DataArray, surface_pressure_data_units: str, vertical_dim: str, vertical_dim_units: str, method: str = 'Trenberth-vibeta', g: float = 9.8) -> xarray.DataArray

   Calculate the water vapor flux across the vertical level.

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
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   method: :py:class:`str <str>`, default: `'Trenberth-vibeta'`.
       Vertical integration method. Optional values are `Boer-vibeta`, `'Trenberth-vibeta'`.
   g: :py:class:`float <float>`, default: `9.8`.
       The acceleration of gravity.

   Returns
   -------
   The water vapor flux. (:py:class:`xarray.Dataset<xarray.Dataset>`).

   - :math:`qu`: zonal water vapor flux.
   - :math:`qv`: meridional water vapor flux.

   .. seealso::
       :py:func:`calc_top2surface_integral <calc_top2surface_integral>`


.. py:function:: calc_divergence_watervaporflux(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, specific_humidity_data_units: str, spherical_coord: bool = True, lon_dim: str = 'lon', lat_dim: str = 'lat', g: float = 9.8, R: float = 6370000) -> xarray.DataArray

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
       Whether or not to compute the horizontal Laplace term in spherical coordinates.
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


.. py:function:: calc_divergence_watervaporflux_top2surface_integral(specific_humidity_data: xarray.DataArray, u_data: xarray.DataArray, v_data: xarray.DataArray, surface_pressure_data: xarray.DataArray, vertical_dim: str, specific_humidity_data_units: str, surface_pressure_data_units: str, vertical_dim_units: str, spherical_coord: bool = True, lon_dim: str = 'lon', lat_dim: str = 'lat', method: str = 'Trenberth-vibeta', g: float = 9.8, R: float = 6370000) -> xarray.DataArray

   Calculate water vapor flux divergence across the vertical level.

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
       Whether or not to compute the horizontal Laplace term in spherical coordinates.
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


.. py:function:: calc_u_advection(u_data: xarray.DataArray, temper_data: xarray.DataArray, lon_dim: str = 'lon', lat_dim: str = 'lat', min_dx: float = 1.0, edge_order: int = 2, R: float = 6370000) -> xarray.DataArray

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


.. py:function:: calc_v_advection(v_data: xarray.DataArray, temper_data: xarray.DataArray, lat_dim: str = 'lat', min_dy: float = 1.0, edge_order: int = 2, R: float = 6370000) -> xarray.DataArray

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


.. py:function:: calc_p_advection(omega_data: xarray.DataArray, temper_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: str) -> xarray.DataArray

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


