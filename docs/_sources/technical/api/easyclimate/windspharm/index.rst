:py:mod:`easyclimate.windspharm`
================================

.. py:module:: easyclimate.windspharm

.. autoapi-nested-parse::

   Spherical harmonic vector wind analysis.



Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   _common/index.rst
   standard/index.rst
   tools/index.rst
   top/index.rst
   xarray/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   easyclimate.windspharm.VectorWind
   easyclimate.windspharm.VectorWind



Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.windspharm.calc_wind_speed
   easyclimate.windspharm.calc_relative_vorticity_and_horizontal_divergence
   easyclimate.windspharm.calc_relative_vorticity
   easyclimate.windspharm.calc_divergence
   easyclimate.windspharm.calc_planetary_vorticity
   easyclimate.windspharm.calc_absolute_vorticity
   easyclimate.windspharm.calc_streamfunction_and_velocity_potential
   easyclimate.windspharm.calc_streamfunction
   easyclimate.windspharm.calc_velocity_potential
   easyclimate.windspharm.calc_helmholtz
   easyclimate.windspharm.calc_irrotational_component
   easyclimate.windspharm.calc_nondivergent_component
   easyclimate.windspharm.calc_rossby_wave_source
   easyclimate.windspharm.calc_gradient
   easyclimate.windspharm._format_lat_lon_coordinate
   easyclimate.windspharm.get_apiorder
   easyclimate.windspharm.inspect_gridtype
   easyclimate.windspharm.to3d
   easyclimate.windspharm._reverse
   easyclimate.windspharm._find_coord_and_dim
   easyclimate.windspharm._find_latitude_coordinate
   easyclimate.windspharm._find_longitude_coordinate



.. py:class:: VectorWind(u, v, rsphere=6371200.0, legfunc='stored')

   Bases: :py:obj:`object`

   Vector wind computations (`xarray` interface).

   .. py:method:: _metadata(var, name, **attributes)


   .. py:method:: u()

      Zonal component of vector wind.

      **Returns:**

      *u*
          The zonal component of the wind.

      **Example:**

      Get the zonal component of the vector wind::

          u = w.u()



   .. py:method:: v()

      Meridional component of vector wind.

      **Returns:**

      *v*
          The meridional component of the wind.

      **Example:**

      Get the meridional component of the vector wind::

          v = w.v()



   .. py:method:: magnitude()

      Wind speed (magnitude of vector wind).

      **Returns:**

      *speed*
          The wind speed.

      **Example:**

      Get the magnitude of the vector wind::

          spd = w.magnitude()



   .. py:method:: vrtdiv(truncation=None)

      Relative vorticity and horizontal divergence.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *vrt*, *div*
          The relative vorticity and divergence respectively.

      **See also:**

      `~VectorWind.vorticity`, `~VectorWind.divergence`.

      **Examples:**

      Compute the relative vorticity and divergence::

          vrt, div = w.vrtdiv()

      Compute the relative vorticity and divergence and apply spectral
      truncation at triangular T13::

          vrtT13, divT13 = w.vrtdiv(truncation=13)



   .. py:method:: vorticity(truncation=None)

      Relative vorticity.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *vrt*
          The relative vorticity.

      **See also:**

      `~VectorWind.vrtdiv`, `~VectorWind.absolutevorticity`.

      **Examples:**

      Compute the relative vorticity::

          vrt = w.vorticity()

      Compute the relative vorticity and apply spectral truncation at
      triangular T13::

          vrtT13 = w.vorticity(truncation=13)



   .. py:method:: divergence(truncation=None)

      Horizontal divergence.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *div*
          The divergence.

      **See also:**

      `~VectorWind.vrtdiv`.

      **Examples:**

      Compute the divergence::

          div = w.divergence()

      Compute the divergence and apply spectral truncation at
      triangular T13::

          divT13 = w.divergence(truncation=13)



   .. py:method:: planetaryvorticity(omega=None)

      Planetary vorticity (Coriolis parameter).

      **Optional argument:**

      *omega*
          Earth's angular velocity. The default value if not specified
          is 7.292x10**-5 s**-1.

      **Returns:**

      *pvorticity*
          The planetary vorticity.

      **See also:**

      `~VectorWind.absolutevorticity`.

      **Example:**

      Compute planetary vorticity using default values::

          pvrt = w.planetaryvorticity()

      Override the default value for Earth's angular velocity::

          pvrt = w.planetaryvorticity(omega=7.2921150)



   .. py:method:: absolutevorticity(omega=None, truncation=None)

      Absolute vorticity (sum of relative and planetary vorticity).

      **Optional arguments:**

      *omega*
          Earth's angular velocity. The default value if not specified
          is 7.292x10**-5 s**-1.

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *avorticity*
          The absolute (relative + planetary) vorticity.

      **See also:**

      `~VectorWind.vorticity`, `~VectorWind.planetaryvorticity`.

      **Examples:**

      Compute absolute vorticity::

          avrt = w.absolutevorticity()

      Compute absolute vorticity and apply spectral truncation at
      triangular T13, also override the default value for Earth's
      angular velocity::

          avrt = w.absolutevorticity(omega=7.2921150, truncation=13)



   .. py:method:: sfvp(truncation=None)

      Streamfunction and velocity potential.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *sf*, *vp*
          The streamfunction and velocity potential respectively.

      **See also:**

      `~VectorWind.streamfunction`, `~VectorWind.velocitypotential`.

      **Examples:**

      Compute streamfunction and velocity potential::

          sf, vp = w.sfvp()

      Compute streamfunction and velocity potential and apply spectral
      truncation at triangular T13::

          sfT13, vpT13 = w.sfvp(truncation=13)



   .. py:method:: streamfunction(truncation=None)

      Streamfunction.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *sf*
          The streamfunction.

      **See also:**

      `~VectorWind.sfvp`.

      **Examples:**

      Compute streamfunction::

          sf = w.streamfunction()

      Compute streamfunction and apply spectral truncation at
      triangular T13::

          sfT13 = w.streamfunction(truncation=13)



   .. py:method:: velocitypotential(truncation=None)

      Velocity potential.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *vp*
          The velocity potential.

      **See also:**

      `~VectorWind.sfvp`.

      **Examples:**

      Compute velocity potential::

          vp = w.velocity potential()

      Compute velocity potential and apply spectral truncation at
      triangular T13::

          vpT13 = w.velocity potential(truncation=13)



   .. py:method:: helmholtz(truncation=None)

      Irrotational and non-divergent components of the vector wind.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *uchi*, *vchi*, *upsi*, *vpsi*
          Zonal and meridional components of irrotational and
          non-divergent wind components respectively.

      **See also:**

      `~VectorWind.irrotationalcomponent`,
      `~VectorWind.nondivergentcomponent`.

      **Examples:**

      Compute the irrotational and non-divergent components of the
      vector wind::

          uchi, vchi, upsi, vpsi = w.helmholtz()

      Compute the irrotational and non-divergent components of the
      vector wind and apply spectral truncation at triangular T13::

          uchiT13, vchiT13, upsiT13, vpsiT13 = w.helmholtz(truncation=13)



   .. py:method:: irrotationalcomponent(truncation=None)

      Irrotational (divergent) component of the vector wind.

      .. note::

         If both the irrotational and non-divergent components are
         required then `~VectorWind.helmholtz` should be used instead.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *uchi*, *vchi*
          The zonal and meridional components of the irrotational wind
          respectively.

      **See also:**

      `~VectorWind.helmholtz`.

      **Examples:**

      Compute the irrotational component of the vector wind::

          uchi, vchi = w.irrotationalcomponent()

      Compute the irrotational component of the vector wind and apply
      spectral truncation at triangular T13::

          uchiT13, vchiT13 = w.irrotationalcomponent(truncation=13)



   .. py:method:: nondivergentcomponent(truncation=None)

      Non-divergent (rotational) component of the vector wind.

      .. note::

         If both the non-divergent and irrotational components are
         required then `~VectorWind.helmholtz` should be used instead.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *upsi*, *vpsi*
          The zonal and meridional components of the non-divergent
          wind respectively.

      **See also:**

      `~VectorWind.helmholtz`.

      **Examples:**

      Compute the non-divergent component of the vector wind::

          upsi, vpsi = w.nondivergentcomponent()

      Compute the non-divergent component of the vector wind and apply
      spectral truncation at triangular T13::

          upsiT13, vpsiT13 = w.nondivergentcomponent(truncation=13)



   .. py:method:: gradient(chi, truncation=None)

      Computes the vector gradient of a scalar field on the sphere.

      **Argument:**

      *chi*
          A scalar field. It must be a `~xarray.DataArray` with the
          same latitude and longitude dimensions as the vector wind
          components that initialized the `VectorWind` instance.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *uchi*, *vchi*
          The zonal and meridional components of the vector gradient
          respectively.

      **Examples:**

      Compute the vector gradient of absolute vorticity::

          avrt = w.absolutevorticity()
          avrt_zonal, avrt_meridional = w.gradient(avrt)

      Compute the vector gradient of absolute vorticity and apply
      spectral truncation at triangular T13::

          avrt = w.absolutevorticity()
          avrt_zonalT13, avrt_meridionalT13 = w.gradient(avrt, truncation=13)



   .. py:method:: truncate(field, truncation=None)

      Apply spectral truncation to a scalar field.

      This is useful to represent other fields in a way consistent
      with the output of other `VectorWind` methods.

      **Argument:**

      *field*
          A scalar field. It must be a `~xarray.DataArray` with the
          same latitude and longitude dimensions as the vector wind
          components that initialized the `VectorWind` instance.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation. If not specified it will default to
          *nlats - 1* where *nlats* is the number of latitudes.

      **Returns:**

      *truncated_field*
          The field with spectral truncation applied.

      **Examples:**

      Truncate a scalar field to the computational resolution of the
      `VectorWind` instance::

          scalar_field_truncated = w.truncate(scalar_field)

      Truncate a scalar field to T21::

          scalar_field_T21 = w.truncate(scalar_field, truncation=21)




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


.. py:function:: get_apiorder(ndim, latitude_dim, longitude_dim)

   Get the dimension ordering for a transpose to the required API
   dimension ordering.

   **Arguments:**

   *ndim*
       Total number of dimensions to consider.

   *latitude_dim*
       Index of the latitude dimension.

   *longitude_dim*
       Index of the longitude dimension.

   **Returns:**

   *apiorder*
       A list of indices corresponding to the order required to
       conform to the specified API order.

   *reorder*
       The inverse indices corresponding to *apiorder*.



.. py:function:: inspect_gridtype(latitudes)

   Determine a grid type by examining the points of a latitude
   dimension.

   Raises a ValueError if the grid type cannot be determined.

   **Argument:**

   *latitudes*
       An iterable of latitude point values.

   **Returns:**

   *gridtype*
       Either 'gaussian' for a Gaussian grid or 'regular' for an
       equally-spaced grid.



.. py:function:: to3d(array)


.. py:class:: VectorWind(u, v, rsphere=6371200.0, legfunc='stored')

   Bases: :py:obj:`object`

   Vector wind computations (`xarray` interface).

   .. py:method:: _metadata(var, name, **attributes)


   .. py:method:: u()

      Zonal component of vector wind.

      **Returns:**

      *u*
          The zonal component of the wind.

      **Example:**

      Get the zonal component of the vector wind::

          u = w.u()



   .. py:method:: v()

      Meridional component of vector wind.

      **Returns:**

      *v*
          The meridional component of the wind.

      **Example:**

      Get the meridional component of the vector wind::

          v = w.v()



   .. py:method:: magnitude()

      Wind speed (magnitude of vector wind).

      **Returns:**

      *speed*
          The wind speed.

      **Example:**

      Get the magnitude of the vector wind::

          spd = w.magnitude()



   .. py:method:: vrtdiv(truncation=None)

      Relative vorticity and horizontal divergence.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *vrt*, *div*
          The relative vorticity and divergence respectively.

      **See also:**

      `~VectorWind.vorticity`, `~VectorWind.divergence`.

      **Examples:**

      Compute the relative vorticity and divergence::

          vrt, div = w.vrtdiv()

      Compute the relative vorticity and divergence and apply spectral
      truncation at triangular T13::

          vrtT13, divT13 = w.vrtdiv(truncation=13)



   .. py:method:: vorticity(truncation=None)

      Relative vorticity.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *vrt*
          The relative vorticity.

      **See also:**

      `~VectorWind.vrtdiv`, `~VectorWind.absolutevorticity`.

      **Examples:**

      Compute the relative vorticity::

          vrt = w.vorticity()

      Compute the relative vorticity and apply spectral truncation at
      triangular T13::

          vrtT13 = w.vorticity(truncation=13)



   .. py:method:: divergence(truncation=None)

      Horizontal divergence.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *div*
          The divergence.

      **See also:**

      `~VectorWind.vrtdiv`.

      **Examples:**

      Compute the divergence::

          div = w.divergence()

      Compute the divergence and apply spectral truncation at
      triangular T13::

          divT13 = w.divergence(truncation=13)



   .. py:method:: planetaryvorticity(omega=None)

      Planetary vorticity (Coriolis parameter).

      **Optional argument:**

      *omega*
          Earth's angular velocity. The default value if not specified
          is 7.292x10**-5 s**-1.

      **Returns:**

      *pvorticity*
          The planetary vorticity.

      **See also:**

      `~VectorWind.absolutevorticity`.

      **Example:**

      Compute planetary vorticity using default values::

          pvrt = w.planetaryvorticity()

      Override the default value for Earth's angular velocity::

          pvrt = w.planetaryvorticity(omega=7.2921150)



   .. py:method:: absolutevorticity(omega=None, truncation=None)

      Absolute vorticity (sum of relative and planetary vorticity).

      **Optional arguments:**

      *omega*
          Earth's angular velocity. The default value if not specified
          is 7.292x10**-5 s**-1.

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *avorticity*
          The absolute (relative + planetary) vorticity.

      **See also:**

      `~VectorWind.vorticity`, `~VectorWind.planetaryvorticity`.

      **Examples:**

      Compute absolute vorticity::

          avrt = w.absolutevorticity()

      Compute absolute vorticity and apply spectral truncation at
      triangular T13, also override the default value for Earth's
      angular velocity::

          avrt = w.absolutevorticity(omega=7.2921150, truncation=13)



   .. py:method:: sfvp(truncation=None)

      Streamfunction and velocity potential.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *sf*, *vp*
          The streamfunction and velocity potential respectively.

      **See also:**

      `~VectorWind.streamfunction`, `~VectorWind.velocitypotential`.

      **Examples:**

      Compute streamfunction and velocity potential::

          sf, vp = w.sfvp()

      Compute streamfunction and velocity potential and apply spectral
      truncation at triangular T13::

          sfT13, vpT13 = w.sfvp(truncation=13)



   .. py:method:: streamfunction(truncation=None)

      Streamfunction.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *sf*
          The streamfunction.

      **See also:**

      `~VectorWind.sfvp`.

      **Examples:**

      Compute streamfunction::

          sf = w.streamfunction()

      Compute streamfunction and apply spectral truncation at
      triangular T13::

          sfT13 = w.streamfunction(truncation=13)



   .. py:method:: velocitypotential(truncation=None)

      Velocity potential.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *vp*
          The velocity potential.

      **See also:**

      `~VectorWind.sfvp`.

      **Examples:**

      Compute velocity potential::

          vp = w.velocity potential()

      Compute velocity potential and apply spectral truncation at
      triangular T13::

          vpT13 = w.velocity potential(truncation=13)



   .. py:method:: helmholtz(truncation=None)

      Irrotational and non-divergent components of the vector wind.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *uchi*, *vchi*, *upsi*, *vpsi*
          Zonal and meridional components of irrotational and
          non-divergent wind components respectively.

      **See also:**

      `~VectorWind.irrotationalcomponent`,
      `~VectorWind.nondivergentcomponent`.

      **Examples:**

      Compute the irrotational and non-divergent components of the
      vector wind::

          uchi, vchi, upsi, vpsi = w.helmholtz()

      Compute the irrotational and non-divergent components of the
      vector wind and apply spectral truncation at triangular T13::

          uchiT13, vchiT13, upsiT13, vpsiT13 = w.helmholtz(truncation=13)



   .. py:method:: irrotationalcomponent(truncation=None)

      Irrotational (divergent) component of the vector wind.

      .. note::

         If both the irrotational and non-divergent components are
         required then `~VectorWind.helmholtz` should be used instead.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *uchi*, *vchi*
          The zonal and meridional components of the irrotational wind
          respectively.

      **See also:**

      `~VectorWind.helmholtz`.

      **Examples:**

      Compute the irrotational component of the vector wind::

          uchi, vchi = w.irrotationalcomponent()

      Compute the irrotational component of the vector wind and apply
      spectral truncation at triangular T13::

          uchiT13, vchiT13 = w.irrotationalcomponent(truncation=13)



   .. py:method:: nondivergentcomponent(truncation=None)

      Non-divergent (rotational) component of the vector wind.

      .. note::

         If both the non-divergent and irrotational components are
         required then `~VectorWind.helmholtz` should be used instead.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *upsi*, *vpsi*
          The zonal and meridional components of the non-divergent
          wind respectively.

      **See also:**

      `~VectorWind.helmholtz`.

      **Examples:**

      Compute the non-divergent component of the vector wind::

          upsi, vpsi = w.nondivergentcomponent()

      Compute the non-divergent component of the vector wind and apply
      spectral truncation at triangular T13::

          upsiT13, vpsiT13 = w.nondivergentcomponent(truncation=13)



   .. py:method:: gradient(chi, truncation=None)

      Computes the vector gradient of a scalar field on the sphere.

      **Argument:**

      *chi*
          A scalar field. It must be a `~xarray.DataArray` with the
          same latitude and longitude dimensions as the vector wind
          components that initialized the `VectorWind` instance.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation.

      **Returns:**

      *uchi*, *vchi*
          The zonal and meridional components of the vector gradient
          respectively.

      **Examples:**

      Compute the vector gradient of absolute vorticity::

          avrt = w.absolutevorticity()
          avrt_zonal, avrt_meridional = w.gradient(avrt)

      Compute the vector gradient of absolute vorticity and apply
      spectral truncation at triangular T13::

          avrt = w.absolutevorticity()
          avrt_zonalT13, avrt_meridionalT13 = w.gradient(avrt, truncation=13)



   .. py:method:: truncate(field, truncation=None)

      Apply spectral truncation to a scalar field.

      This is useful to represent other fields in a way consistent
      with the output of other `VectorWind` methods.

      **Argument:**

      *field*
          A scalar field. It must be a `~xarray.DataArray` with the
          same latitude and longitude dimensions as the vector wind
          components that initialized the `VectorWind` instance.

      **Optional argument:**

      *truncation*
          Truncation limit (triangular truncation) for the spherical
          harmonic computation. If not specified it will default to
          *nlats - 1* where *nlats* is the number of latitudes.

      **Returns:**

      *truncated_field*
          The field with spectral truncation applied.

      **Examples:**

      Truncate a scalar field to the computational resolution of the
      `VectorWind` instance::

          scalar_field_truncated = w.truncate(scalar_field)

      Truncate a scalar field to T21::

          scalar_field_T21 = w.truncate(scalar_field, truncation=21)




.. py:function:: _reverse(array, dim)

   Reverse an `xarray.DataArray` along a given dimension.


.. py:function:: _find_coord_and_dim(array, predicate, name)

   Find a dimension coordinate in an `xarray.DataArray` that satisfies
   a predicate function.



.. py:function:: _find_latitude_coordinate(array)

   Find a latitude dimension coordinate in an `xarray.DataArray`.


.. py:function:: _find_longitude_coordinate(array)

   Find a longitude dimension coordinate in an `xarray.DataArray`.


