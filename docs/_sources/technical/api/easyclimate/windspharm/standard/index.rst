:py:mod:`easyclimate.windspharm.standard`
=========================================

.. py:module:: easyclimate.windspharm.standard

.. autoapi-nested-parse::

   Spherical harmonic vector wind computations.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   easyclimate.windspharm.standard.VectorWind




.. py:class:: VectorWind(u, v, gridtype='regular', rsphere=6371200.0, legfunc='stored')

   Bases: :py:obj:`object`

   Vector Wind computations (standard `numpy` interface).

   .. py:method:: magnitude()

      Wind speed (magnitude of vector wind).

      **Returns:**

      *speed*
          The wind speed.

      **Example:**

      Magnitude of the vector wind::

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
          A scalar field. Its shape must be either (nlat, nlon) or
          (nlat, nlon, nfields) where nlat and nlon are the same
          as those for the vector wind components that initialized the
          `VectorWind` instance.

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
          A scalar field. Its shape must be either (nlat, nlon) or
          (nlat, nlon, nfields) where nlat and nlon are the same
          as those for the vector wind components that initialized the
          `VectorWind` instance.

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




