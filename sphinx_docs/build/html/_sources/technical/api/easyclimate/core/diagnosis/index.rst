:py:mod:`easyclimate.core.diagnosis`
====================================

.. py:module:: easyclimate.core.diagnosis

.. autoapi-nested-parse::

   Functions for transient eddy.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.core.diagnosis.calc_brunt_vaisala_frequency_atm
   easyclimate.core.diagnosis.get_coriolis_parameter
   easyclimate.core.diagnosis.get_potential_temperature
   easyclimate.core.diagnosis.calc_static_stability



.. py:function:: calc_brunt_vaisala_frequency_atm(potential_temperature_data: xarray.DataArray, z_data: xarray.DataArray, vertical_dim: str, g=9.8) -> xarray.DataArray

   Calculation of the Brunt-väisälä frequency for the vertical atmosphere.

   .. math::
       N = \left( \frac{g}{\theta} \frac{\mathrm{d}\theta}{\mathrm{d}z} \right)^\frac{1}{2}

   Parameters
   ----------
   potential_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Vertical atmospheric potential temperature.
   z_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Vertical atmospheric geopotential height.

   .. attention:: The unit of `z_data` should be **meters**, NOT :math:`\mathrm{m^2 \cdot s^2}` which is the unit used in the representation of potential energy.

   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   g: :py:class:`float<python.float>`, default: `9.8`.
       The acceleration of gravity.

   Returns
   -------
   Brunt-väisälä frequency (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - `Brunt-väisälä frequency - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Brunt-v%C3%A4is%C3%A4l%C3%A4_frequency>`__

   .. seealso::
       - `brunt_vaisala_frequency — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.brunt_vaisala_frequency.html>`__
       - `brunt_vaisala_atm - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/brunt_vaisala_atm.shtml>`__
       


.. py:function:: get_coriolis_parameter(lat_data, omega=7.292e-05) -> xarray.DataArray

   Calculate the Coriolis parameter at each point.

   .. math::
       f = 2 \Omega \sin(\phi)

   Parameters
   ----------
   lat_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Latitude at each point.
   omega: :py:class:`float<python.float>`, default: `7.292e-5`.
       The angular speed of the earth.

   Returns
   -------
   Corresponding Coriolis force at each point (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - `Coriolis parameter - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Coriolis_parameter>`__

   .. seealso::
       - `coriolis_parameter — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.coriolis_parameter.html>`__
       - `coriolis_param - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/coriolis_param.shtml>`__


.. py:function:: get_potential_temperature(temper_data: xarray.DataArray, vertical_dim: str, vertical_units: str, kappa=287 / 1005.7) -> xarray.DataArray

   Calculate the potential temperature.

   Uses the Poisson equation to calculation the potential temperature given pressure and temperature.

   .. math::
       \theta = T \left( \frac{p_0}{p} \right) ^\kappa

   Parameters
   ----------
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   kappa: :py:class:`float<python.float>`, default: `287/1005.7`.
       Poisson constant :math:`\kappa`.

       .. note::
           `Poisson constant - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Poisson_constant>`__

   Returns
   -------
   Potential temperature corresponding to the temperature and pressure (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - `Potential temperature - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Potential_temperature>`__
   - `Potential-temperature.pdf <http://weatherclimatelab.mit.edu/wp-content/uploads/2018/02/Potential-temperature.pdf>`__
   - `大气位温、相当位温、饱和相当位温、静力稳定度 <https://renqlsysu.github.io/2019/10/23/potential_temperature/>`__

   .. seealso::
       - `potential_temperature — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.potential_temperature.html>`__
       - `pot_temp - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/pot_temp.shtml>`__


.. py:function:: calc_static_stability(temper_data: xarray.DataArray, vertical_dim: str, vertical_units: str) -> xarray.DataArray

   Calculate the static stability within a vertical profile.

   .. math::
       \sigma = - T \frac{\partial \ln \theta}{\partial p}

   Parameters
   ----------
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.

   Returns
   -------
   Static stability (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Howard B. Bluestein. (1992). Synoptic-Dynamic Meteorology in Midlatitudes: Principles of Kinematics and Dynamics, Vol. 1

   .. seealso::
       - `static_stability - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/static_stability.shtml>`__
       - `static_stability — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.static_stability.html>`__
       - `Static stability parameters · Issue #2535 · Unidata/MetPy <https://github.com/Unidata/MetPy/issues/2535>`__


