:py:mod:`easyclimate.windspharm.top`
====================================

.. py:module:: easyclimate.windspharm.top

.. autoapi-nested-parse::

   easy cliamte top interface for the windspharm

   This is the top layer of packaging for the windspharm package.

   .. seealso::
       - `Spherical harmonic vector wind computations (xarray interface) <https://ajdawson.github.io/windspharm/latest/api/windspharm.xarray.html>`__
       - https://github.com/ajdawson/windspharm
       - Dawson, A. (2016). Windspharm: A High-Level Library for Global Wind Field Computations Using Spherical Harmonics. Journal of Open Research Software, 4(1), e31.DOI: https://doi.org/10.5334/jors.129



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.windspharm.top.calc_wind_speed
   easyclimate.windspharm.top.calc_relative_vorticity_and_horizontal_divergence
   easyclimate.windspharm.top.calc_relative_vorticity
   easyclimate.windspharm.top.calc_divergence
   easyclimate.windspharm.top.calc_planetary_vorticity
   easyclimate.windspharm.top.calc_absolute_vorticity
   easyclimate.windspharm.top.calc_streamfunction_and_velocity_potential
   easyclimate.windspharm.top.calc_streamfunction
   easyclimate.windspharm.top.calc_velocity_potential
   easyclimate.windspharm.top.calc_helmholtz
   easyclimate.windspharm.top.calc_irrotational_component
   easyclimate.windspharm.top.calc_nondivergent_component
   easyclimate.windspharm.top.calc_rossby_wave_source
   easyclimate.windspharm.top.calc_gradient
   easyclimate.windspharm.top._format_lat_lon_coordinate



.. py:function:: calc_wind_speed(u_data: xarray.DataArray, v_data: xarray.DataArray, R: float = 6371200.0, legfunc: str = 'stored', lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.DataArray

   Calculate the wind speed (magnitude of vector wind).

   Parameters
   ----------
   u_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal component of the wind.
   v_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional component of vector wind.
   R: :py:class:`float`.
       The radius in metres of the sphere used in the spherical
       harmonic computations. Default is 6371200 m, the approximate
       mean spherical Earth radius.
   legfunc: :py:class:`str`, 'stored' (default) or 'computed'.
       If 'stored', associated legendre
       functions are precomputed and stored when the class instance is
       created. This uses :math:`O(\mathrm{nlat}^3)` memory, but speeds up the spectral
       transforms. If 'computed', associated legendre functions are
       computed on the fly when transforms are requested. This uses
       :math:`O(\mathrm{nlat}^2)` memory, but slows down the spectral transforms a bit.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.    

   Returns
   -------
   The wind speed (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_relative_vorticity_and_horizontal_divergence(u_data: xarray.DataArray, v_data: xarray.DataArray, truncation: int = None, R: float = 6371200.0, legfunc: str = 'stored', lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.Dataset

   Calculate relative vorticity and horizontal divergence.

   Parameters
   ----------
   u_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal component of the wind.
   v_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional component of vector wind.
   truncation: :py:class:`int`.
       Truncation limit (triangular truncation) for the spherical harmonic computation.
   R: :py:class:`float`.
       The radius in metres of the sphere used in the spherical
       harmonic computations. Default is 6371200 m, the approximate
       mean spherical Earth radius.
   legfunc: :py:class:`str`, 'stored' (default) or 'computed'.
       If 'stored', associated legendre
       functions are precomputed and stored when the class instance is
       created. This uses :math:`O(\mathrm{nlat}^3)` memory, but speeds up the spectral
       transforms. If 'computed', associated legendre functions are
       computed on the fly when transforms are requested. This uses
       :math:`O(\mathrm{nlat}^2)` memory, but slows down the spectral transforms a bit.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   Relative vorticity and horizontal divergence (:py:class:`xarray.Dataset<xarray.Dataset>`).


.. py:function:: calc_relative_vorticity(u_data: xarray.DataArray, v_data: xarray.DataArray, truncation: int = None, R: float = 6371200.0, legfunc: str = 'stored', lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.DataArray

   Calculate relative vorticity.

   Parameters
   ----------
   u_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal component of the wind.
   v_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional component of vector wind.
   truncation: :py:class:`int`.
       Truncation limit (triangular truncation) for the spherical harmonic computation.
   R: :py:class:`float`.
       The radius in metres of the sphere used in the spherical
       harmonic computations. Default is 6371200 m, the approximate
       mean spherical Earth radius.
   legfunc: :py:class:`str`, 'stored' (default) or 'computed'.
       If 'stored', associated legendre
       functions are precomputed and stored when the class instance is
       created. This uses :math:`O(\mathrm{nlat}^3)` memory, but speeds up the spectral
       transforms. If 'computed', associated legendre functions are
       computed on the fly when transforms are requested. This uses
       :math:`O(\mathrm{nlat}^2)` memory, but slows down the spectral transforms a bit.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   Relative vorticity (:py:class:`xarray.DataArray<xarray.DataArray>`).    


.. py:function:: calc_divergence(u_data: xarray.DataArray, v_data: xarray.DataArray, truncation: int = None, R: float = 6371200.0, legfunc: str = 'stored', lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.DataArray

   Calculate horizontal divergence.

   Parameters
   ----------
   u_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal component of the wind.
   v_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional component of vector wind.
   truncation: :py:class:`int`.
       Truncation limit (triangular truncation) for the spherical harmonic computation.
   R: :py:class:`float`.
       The radius in metres of the sphere used in the spherical
       harmonic computations. Default is 6371200 m, the approximate
       mean spherical Earth radius.
   legfunc: :py:class:`str`, 'stored' (default) or 'computed'.
       If 'stored', associated legendre
       functions are precomputed and stored when the class instance is
       created. This uses :math:`O(\mathrm{nlat}^3)` memory, but speeds up the spectral
       transforms. If 'computed', associated legendre functions are
       computed on the fly when transforms are requested. This uses
       :math:`O(\mathrm{nlat}^2)` memory, but slows down the spectral transforms a bit.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   Horizontal divergence (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_planetary_vorticity(u_data: xarray.DataArray, v_data: xarray.DataArray, omega: float = 7.292115, R: float = 6371200.0, legfunc: str = 'stored', lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.DataArray

   Calculate planetary vorticity (Coriolis parameter).

   Parameters
   ----------
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal component of the wind.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional component of vector wind.
   omega: :py:class:`float`.
       Earth's angular velocity. The default value if not specified is :math:`7.292    imes 10^{-5} \mathrm{s^{-1}}`.
   R: :py:class:`float`.
       The radius in metres of the sphere used in the spherical
       harmonic computations. Default is 6371200 m, the approximate
       mean spherical Earth radius.
   legfunc: :py:class:`str`, 'stored' (default) or 'computed'.
       If 'stored', associated legendre
       functions are precomputed and stored when the class instance is
       created. This uses :math:`O(\mathrm{nlat}^3)` memory, but speeds up the spectral
       transforms. If 'computed', associated legendre functions are
       computed on the fly when transforms are requested. This uses
       :math:`O(\mathrm{nlat}^2)` memory, but slows down the spectral transforms a bit.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   Planetary vorticity (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_absolute_vorticity(u_data: xarray.DataArray, v_data: xarray.DataArray, truncation: int = None, omega: float = 7.292115, R: float = 6371200.0, legfunc: str = 'stored', lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.DataArray

   Calculate absolute vorticity (sum of relative and planetary vorticity).

   Parameters
   ----------
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal component of the wind.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional component of vector wind.
   truncation: :py:class:`int`.
       Truncation limit (triangular truncation) for the spherical harmonic computation.
   omega: :py:class:`float`.
       Earth's angular velocity. The default value if not specified is :math:`7.292    imes 10^{-5} \mathrm{s^{-1}}`.
   R: :py:class:`float`.
       The radius in metres of the sphere used in the spherical
       harmonic computations. Default is 6371200 m, the approximate
       mean spherical Earth radius.
   legfunc: :py:class:`str`, 'stored' (default) or 'computed'.
       If 'stored', associated legendre
       functions are precomputed and stored when the class instance is
       created. This uses :math:`O(\mathrm{nlat}^3)` memory, but speeds up the spectral
       transforms. If 'computed', associated legendre functions are
       computed on the fly when transforms are requested. This uses
       :math:`O(\mathrm{nlat}^2)` memory, but slows down the spectral transforms a bit.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   Absolute vorticity (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_streamfunction_and_velocity_potential(u_data: xarray.DataArray, v_data: xarray.DataArray, truncation: int = None, R: float = 6371200.0, legfunc: str = 'stored', lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.Dataset

   Calculate stream function and velocity potential.

   Parameters
   ----------
   u_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal component of the wind.
   v_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional component of vector wind.    
   truncation: :py:class:`int`.
       Truncation limit (triangular truncation) for the spherical harmonic computation.
   R: :py:class:`float`.
       The radius in metres of the sphere used in the spherical
       harmonic computations. Default is 6371200 m, the approximate
       mean spherical Earth radius.
   legfunc: :py:class:`str`, 'stored' (default) or 'computed'.
       If 'stored', associated legendre
       functions are precomputed and stored when the class instance is
       created. This uses :math:`O(\mathrm{nlat}^3)` memory, but speeds up the spectral
       transforms. If 'computed', associated legendre functions are
       computed on the fly when transforms are requested. This uses
       :math:`O(\mathrm{nlat}^2)` memory, but slows down the spectral transforms a bit.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   Stream function and velocity potential (:py:class:`xarray.Dataset<xarray.Dataset>`).


.. py:function:: calc_streamfunction(u_data: xarray.DataArray, v_data: xarray.DataArray, truncation: int = None, R: float = 6371200.0, legfunc: str = 'stored', lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.DataArray

   Calculate stream function.

   Parameters
   ----------
   u_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal component of the wind.
   v_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional component of vector wind.    
   truncation: :py:class:`int`.
       Truncation limit (triangular truncation) for the spherical harmonic computation.
   R: :py:class:`float`.
       The radius in metres of the sphere used in the spherical
       harmonic computations. Default is 6371200 m, the approximate
       mean spherical Earth radius.
   legfunc: :py:class:`str`, 'stored' (default) or 'computed'.
       If 'stored', associated legendre
       functions are precomputed and stored when the class instance is
       created. This uses :math:`O(\mathrm{nlat}^3)` memory, but speeds up the spectral
       transforms. If 'computed', associated legendre functions are
       computed on the fly when transforms are requested. This uses
       :math:`O(\mathrm{nlat}^2)` memory, but slows down the spectral transforms a bit.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   stream function (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_velocity_potential(u_data: xarray.DataArray, v_data: xarray.DataArray, truncation: int = None, R: float = 6371200.0, legfunc: str = 'stored', lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.DataArray

   Calculate velocity potential.

   Parameters
   ----------
   u_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal component of the wind.
   v_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional component of vector wind.    
   truncation: :py:class:`int`.
       Truncation limit (triangular truncation) for the spherical harmonic computation.
   R: :py:class:`float`.
       The radius in metres of the sphere used in the spherical
       harmonic computations. Default is 6371200 m, the approximate
       mean spherical Earth radius.
   legfunc: :py:class:`str`, 'stored' (default) or 'computed'.
       If 'stored', associated legendre
       functions are precomputed and stored when the class instance is
       created. This uses :math:`O(\mathrm{nlat}^3)` memory, but speeds up the spectral
       transforms. If 'computed', associated legendre functions are
       computed on the fly when transforms are requested. This uses
       :math:`O(\mathrm{nlat}^2)` memory, but slows down the spectral transforms a bit.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   Velocity potential (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_helmholtz(u_data: xarray.DataArray, v_data: xarray.DataArray, truncation: int = None, R: float = 6371200.0, legfunc: str = 'stored', lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.Dataset

   Calculate irrotational and non-divergent components of the vector wind.

   Parameters
   ----------
   u_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal component of the wind.
   v_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional component of vector wind.    
   truncation: :py:class:`int`.
       Truncation limit (triangular truncation) for the spherical harmonic computation.
   R: :py:class:`float`.
       The radius in metres of the sphere used in the spherical
       harmonic computations. Default is 6371200 m, the approximate
       mean spherical Earth radius.
   legfunc: :py:class:`str`, 'stored' (default) or 'computed'.
       If 'stored', associated legendre
       functions are precomputed and stored when the class instance is
       created. This uses :math:`O(\mathrm{nlat}^3)` memory, but speeds up the spectral
       transforms. If 'computed', associated legendre functions are
       computed on the fly when transforms are requested. This uses
       :math:`O(\mathrm{nlat}^2)` memory, but slows down the spectral transforms a bit.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   Irrotational and non-divergent components of the vector wind (:py:class:`xarray.Dataset<xarray.Dataset>`).


.. py:function:: calc_irrotational_component(u_data: xarray.DataArray, v_data: xarray.DataArray, truncation: int = None, R: float = 6371200.0, legfunc: str = 'stored', lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.Dataset

   Calculate irrotational (divergent) component of the vector wind.

   Parameters
   ----------
   u_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal component of the wind.
   v_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional component of vector wind.    
   truncation: :py:class:`int`.
       Truncation limit (triangular truncation) for the spherical harmonic computation.
   R: :py:class:`float`.
       The radius in metres of the sphere used in the spherical
       harmonic computations. Default is 6371200 m, the approximate
       mean spherical Earth radius.
   legfunc: :py:class:`str`, 'stored' (default) or 'computed'.
       If 'stored', associated legendre
       functions are precomputed and stored when the class instance is
       created. This uses :math:`O(\mathrm{nlat}^3)` memory, but speeds up the spectral
       transforms. If 'computed', associated legendre functions are
       computed on the fly when transforms are requested. This uses
       :math:`O(\mathrm{nlat}^2)` memory, but slows down the spectral transforms a bit.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   Irrotational (divergent) component of the vector wind (:py:class:`xarray.Dataset<xarray.Dataset>`).


.. py:function:: calc_nondivergent_component(u_data: xarray.DataArray, v_data: xarray.DataArray, truncation: int = None, R: float = 6371200.0, legfunc: str = 'stored', lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.Dataset

   Calculate non-divergent (rotational) component of the vector wind.

   Parameters
   ----------
   u_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal component of the wind.
   v_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional component of vector wind.    
   truncation: :py:class:`int`.
       Truncation limit (triangular truncation) for the spherical harmonic computation.
   R: :py:class:`float`.
       The radius in metres of the sphere used in the spherical
       harmonic computations. Default is 6371200 m, the approximate
       mean spherical Earth radius.
   legfunc: :py:class:`str`, 'stored' (default) or 'computed'.
       If 'stored', associated legendre
       functions are precomputed and stored when the class instance is
       created. This uses :math:`O(\mathrm{nlat}^3)` memory, but speeds up the spectral
       transforms. If 'computed', associated legendre functions are
       computed on the fly when transforms are requested. This uses
       :math:`O(\mathrm{nlat}^2)` memory, but slows down the spectral transforms a bit.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   Non-divergent (rotational) component of the vector wind (:py:class:`xarray.Dataset<xarray.Dataset>`).


.. py:function:: calc_rossby_wave_source(u_data: xarray.DataArray, v_data: xarray.DataArray, truncation: int = None, R: float = 6371200.0, legfunc: str = 'stored', lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.DataArray

   Calculate Rossby wave sources (RWS).

   .. math::
       RWS=-\nabla \cdot \left({v}_{x}\zeta \right)=-\left(\zeta \nabla \cdot {v}_{x}+{v}_{x}\cdot \nabla \zeta \right)

   with :math:`\zeta` being the absolute vorticity.

   Parameters
   ----------
   u_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal component of the wind.
   v_data : :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional component of vector wind.
   truncation: :py:class:`int`.
       Truncation limit (triangular truncation) for the spherical harmonic computation.
   R: :py:class:`float`.
       The radius in metres of the sphere used in the spherical
       harmonic computations. Default is 6371200 m, the approximate
       mean spherical Earth radius.
   legfunc: :py:class:`str`, 'stored' (default) or 'computed'.
       If 'stored', associated legendre
       functions are precomputed and stored when the class instance is
       created. This uses :math:`O(\mathrm{nlat}^3)` memory, but speeds up the spectral
       transforms. If 'computed', associated legendre functions are
       computed on the fly when transforms are requested. This uses
       :math:`O(\mathrm{nlat}^2)` memory, but slows down the spectral transforms a bit.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   Rossby wave sources (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - Sardeshmukh, P. D., & Hoskins, B. J. (1988). The Generation of Global Rotational Flow by Steady Idealized Tropical Divergence. Journal of Atmospheric Sciences, 45(7), 1228-1251. https://doi.org/10.1175/1520-0469(1988)045<1228:TGOGRF>2.0.CO;2
   - James IN (1994) Low frequency variability of the circulation. Introduction to Circulating Atmospheres. Cambridge University Press, Cambridge, UK, pp 255–301
   - Trenberth, K. E., Branstator, G. W., Karoly, D., Kumar, A., Lau, N.-C., and Ropelewski, C. (1998), Progress during TOGA in understanding and modeling global teleconnections associated with tropical sea surface temperatures, J. Geophys. Res., 103(C7), 14291–14324, doi: https://doi.org/10.1029/97JC01444.
   - Nie, Y., Zhang, Y., Yang, X.-Q., & Ren, H.-L. (2019). Winter and summer Rossby wave sources in the CMIP5 models. Earth and Space Science, 6, 1831–1846. https://doi.org/10.1029/2019EA000674
   - Fuentes-Franco, R., Koenigk, T., Docquier, D. et al. Exploring the influence of the North Pacific Rossby wave sources on the variability of summer atmospheric circulation and precipitation over the Northern Hemisphere. Clim Dyn 59, 2025–2039 (2022). https://doi.org/10.1007/s00382-022-06194-4


.. py:function:: calc_gradient(data_input: xarray.DataArray, truncation: int = None, R: float = 6371200.0, legfunc: str = 'stored', lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.Dataset

   Computes the vector gradient of a scalar field on the sphere.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       The spatio-temporal data to be calculated.
   truncation: :py:class:`int`.
       Truncation limit (triangular truncation) for the spherical harmonic computation.
   R: :py:class:`float`.
       The radius in metres of the sphere used in the spherical
       harmonic computations. Default is 6371200 m, the approximate
       mean spherical Earth radius.
   legfunc: :py:class:`str`, 'stored' (default) or 'computed'.
       If 'stored', associated legendre
       functions are precomputed and stored when the class instance is
       created. This uses :math:`O(\mathrm{nlat}^3)` memory, but speeds up the spectral
       transforms. If 'computed', associated legendre functions are
       computed on the fly when transforms are requested. This uses
       :math:`O(\mathrm{nlat}^2)` memory, but slows down the spectral transforms a bit.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The zonal and meridional components of the vector gradient respectively (:py:class:`xarray.Dataset<xarray.Dataset>`).


.. py:function:: _format_lat_lon_coordinate(data_input: xarray.DataArray, lat_dim: str, lon_dim: str) -> xarray.DataArray

   Add attrs to lat lon coordinate


