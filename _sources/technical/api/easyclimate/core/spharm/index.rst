easyclimate.core.spharm
=======================

.. py:module:: easyclimate.core.spharm

.. autoapi-nested-parse::

   Easy climate top interface for the Pyspharm

   This is the top layer of packaging for the Pyspharm package.



Functions
---------

.. autoapisummary::

   easyclimate.core.spharm.calc_gaussian_latitudes
   easyclimate.core.spharm.calc_geodesic_points
   easyclimate.core.spharm.calc_spherical_harmonic_coefficients
   easyclimate.core.spharm.calc_legendre_functions
   easyclimate.core.spharm.transfer_grid2spectral_transform
   easyclimate.core.spharm.transfer_spectral_transform2grid


Module Contents
---------------

.. py:function:: calc_gaussian_latitudes(nlat: int)

   Calculate the gaussian latitudes (in degrees) and quadrature weights.

   Parameters
   ----------
   nlat: :py:class:`int <int>`.
       Number of gaussian latitudes desired.

   Returns
   -------
   The gaussian latitudes (in degrees north) and gaussian quadrature weights (:py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::

       - https://www.ncl.ucar.edu/Document/Functions/Built-in/gaus.shtml


.. py:function:: calc_geodesic_points(m: int)

   Calculate the lat/lon values of the points on the surface of the sphere
   corresponding to a twenty-sided (icosahedral) geodesic.

   Parameters
   ----------
   m: :py:class:`int <int>`.
       The number of points on the edge of a single geodesic triangle.
       There are :math:`10(m-1)^2+2` total geodesic points, including the poles.

   Returns
   -------
   The latitudes and longitudes of the geodesic points (in degrees).
   These points are nearly evenly distributed on the surface of the sphere. (:py:class:`xarray.Dataset<xarray.Dataset>`).


.. py:function:: calc_spherical_harmonic_coefficients(ntrunc: int)

   Calculate indices of zonal wavenumber (indxm) and degree (indxn) for complex spherical harmonic coefficients.

   Parameters
   ----------
   ntrunc: :py:class:`int <int>`.
       The spherical harmonic triangular truncation limit, i.e, truncation wavenumber (e.g., ``T42``).

   Returns
   -------
   The latitudes and longitudes of the geodesic points (in degrees).
   These points are nearly evenly distributed on the surface of the sphere. (:py:class:`xarray.Dataset<xarray.Dataset>`).


.. py:function:: calc_legendre_functions(lat: float, ntrunc: int) -> xarray.DataArray

   Calculate associated legendre functions for triangular truncation T(ntrunc), at a given latitude.

   Parameters
   ----------
   lat: :py:class:`float <float>`.
       The latitude (in degrees) to compute the associate legendre functions.
   ntrunc: :py:class:`int <int>`.
       The spherical harmonic triangular truncation limit, i.e, truncation wavenumber (e.g., ``T42``).

   Returns
   -------
   :math:`(\mathrm{ntrunc} + 1) (\mathrm{ntrunc} + 2) /2` associated legendre functions at latitude ``lat``.


.. py:function:: transfer_grid2spectral_transform(grid_data: xarray.DataArray, grid_data_type: Literal['regular', 'gaussian'], lon_dim: str = 'lon', lat_dim: str = 'lat', ntrunc: int = None) -> xarray.DataArray

   Transform grid data to spectral representation (spherical harmonic analysis).

   Parameters
   ------------
   grid_data: :py:class:`xarray.DataArray <xarray.DataArray>`.
       Input grid data, must contain longitude and latitude dimensions
   grid_data_type:
       Type of grid ('regular' or 'gaussian')
   lon_dim: :py:class:`str <str>`.
       Name of longitude dimension, default is ``'lon'``.
   lat_dim: :py:class:`str <str>`.
       Name of latitude dimension, default is ``'lat'``.
   ntrunc: :py:class:`int <int>`.
       Spectral truncation wavenumber, defaults to ``nlat-1``.

   Returns
   ------------

   :py:class:`xarray.DataArray <xarray.DataArray>` containing complex spherical harmonic coefficients with triangular spectral dimension


.. py:function:: transfer_spectral_transform2grid(spec_data: xarray.DataArray, nlon: int, nlat: int, grid_data_type: Literal['regular', 'gaussian'], spec_dim: str = 'spec_dim', lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.DataArray

   Transform spectral data back to grid space representation (spherical harmonic synthesis)

   Parameters
   ------------
   spec_data: :py:class:`xarray.DataArray <xarray.DataArray>`.
       Input spectral coefficient data, must contain the spectral dimension.
   nlon: :py:class:`int <int>`.
       Number of longitude points in output grid.
   nlat: :py:class:`int <int>`.
       Number of latitude points in output grid.
   grid_data_type:
       Type of output grid (``'regular'`` or ``'gaussian'``).
   spec_dim:
       Name of spectral dimension, default is ``'spec_dim'``.
   lon_dim: :py:class:`str <str>`.
       Name for output longitude dimension, default is ``'lon'``.
   lat_dim: :py:class:`str <str>`.
       Name for output latitude dimension, default is ``'lat'``.

   Returns
   ------------

   :py:class:`xarray.DataArray <xarray.DataArray>` in grid space representation with ``(lat, lon)`` dimensions


