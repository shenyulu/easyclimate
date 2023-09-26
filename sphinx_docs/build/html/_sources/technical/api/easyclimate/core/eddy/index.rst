:py:mod:`easyclimate.core.eddy`
===============================

.. py:module:: easyclimate.core.eddy

.. autoapi-nested-parse::

   Functions for transient eddy.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.core.eddy.calc_eady_growth_rate
   easyclimate.core.eddy.calc_apparent_heat_source
   easyclimate.core.eddy.calc_total_diabatic_heating
   easyclimate.core.eddy.calc_apparent_moisture_sink
   easyclimate.core.eddy.calc_Plumb_wave_activity_horizontal_flux
   easyclimate.core.eddy.calc_TN_wave_activity_horizontal_flux
   easyclimate.core.eddy.calc_TN_wave_activity_3D_flux
   easyclimate.core.eddy.calc_EP_horizontal_flux



.. py:function:: calc_eady_growth_rate(u_daily_data: xarray.DataArray, z_daily_data: xarray.DataArray, temper_daily_data: xarray.DataArray, vertical_dim: str, vertical_units: str, lat_dim='lat', g=9.8) -> xarray.Dataset

   Calculate the maximum Eady growth rate.

   .. math::
       \sigma = 0.3098 \frac{f}{N} \frac{\mathrm{d} U}{\mathrm{d} z}

   .. caution::
       Eady growth rate (EGR) is a non-linear quantity. Hence, `calc_eady_growth_rate` should **NOT** be **directly applied to monthly means** variables. 
       If a monthly climatology of EGR is desired, the EGR values at the high frequency temporal time steps should be calculated; 
       then, use calculate monthly mean.

   Parameters
   ----------
   u_daily_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind daily data.
   z_daily_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Daily atmospheric geopotential height.
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Daily air temperature.
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   g: :py:class:`float<python.float>`, default: `9.8`.
       The acceleration of gravity. 

   Returns
   -------
   The maximum Eady growth rate. (:py:class:`xarray.Dataset<xarray.Dataset>`).

   - `eady_growth_rate`: The maximum Eady growth rate.
   - `dudz`: :math:`\frac{\mathrm{d} U}{\mathrm{d} z}`
   - `brunt_vaisala_frequency`: Brunt-väisälä frequency.

   .. seealso::
       - `eady_growth_rate -NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/eady_growth_rate.shtml>`__
       - `瞬变涡旋诊断量 <https://renqlsysu.github.io/2020/02/16/wave_activity_flux/>`__


.. py:function:: calc_apparent_heat_source(u_data: xarray.DataArray, v_data: xarray.DataArray, omega_data: xarray.DataArray, temper_data: xarray.DataArray, vertical_dim: str, vertical_units: str, time_units: str, lon_dim='lon', lat_dim='lat', time_dim='time', c_p=1005.7) -> xarray.DataArray

   Calculate the apparent heat source.

   .. math::
       Q_1 = C_p \frac{T}{\theta} \left( \frac{\partial \theta}{\partial t} + u \frac{\partial \theta}{\partial x} + v \frac{\partial \theta}{\partial y} + \omega \frac{\partial \theta}{\partial p} \right)

   Parameters
   ----------
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   omega_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The vertical velocity data (:math:`\frac{\mathrm{d} p}{\mathrm{d} t}`).
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   time_units: :py:class:`str<python.str>`.
       The unit corresponding to the time dimension value. Optional values are `seconds`, `months`, `years` and so on.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   time_dim: :py:class:`str<python.str>`.
       The time coordinate dimension name.
   c_p: :py:class:`float<python.float>`, default: `1005.7`.
       The specific heat at constant pressure of dry air.

       .. note::
           `specific heat capacity - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Specific_heat_capacity>`__

   Returns
   -------
   The apparent heat source (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::
       - `Yanai, M., & Tomita, T. (1998). Seasonal and Interannual Variability of Atmospheric Heat Sources and Moisture Sinks as Determined from NCEP–NCAR Reanalysis, Journal of Climate, 11(3), 463-482. <https://journals.ametsoc.org/view/journals/clim/11/3/1520-0442_1998_011_0463_saivoa_2.0.co_2.xml>`__
       - `Ling, J., & Zhang, C. (2013). Diabatic Heating Profiles in Recent Global Reanalyses, Journal of Climate, 26(10), 3307-3325. <https://doi.org/10.1175/JCLI-D-12-00384.1>`__


.. py:function:: calc_total_diabatic_heating(u_data: xarray.DataArray, v_data: xarray.DataArray, omega_data: xarray.DataArray, temper_data: xarray.DataArray, vertical_dim: str, vertical_units: str, time_units: str, lat_dim='lat', lon_dim='lon', time_dim='time', c_p=1005.7) -> xarray.DataArray

   Calculate the total diabatic heating.

   Calculated in exactly the same way as for the apparent heat source.

   Parameters
   ----------
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   omega_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The vertical velocity data (:math:`\frac{\mathrm{d} p}{\mathrm{d} t}`).
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   time_units: :py:class:`str<python.str>`.
       The unit corresponding to the time dimension value. Optional values are `seconds`, `months`, `years` and so on.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   time_dim: :py:class:`str<python.str>`, default: `time`.
       The time coordinate dimension name.
   c_p: :py:class:`float<python.float>`, default: `1005.7` (:math:`\mathrm{J \cdot kg^{-1} \cdot K^{-1}}`).
       The specific heat at constant pressure of dry air.

       .. note::
           `specific heat capacity - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Specific_heat_capacity>`__

   Returns
   -------
   The total diabatic heating (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::
       :py:func:`calc_apparent_heat_source <calc_apparent_heat_source>`


.. py:function:: calc_apparent_moisture_sink(u_data: xarray.DataArray, v_data: xarray.DataArray, omega_data: xarray.DataArray, specific_humidity_data: xarray.DataArray, vertical_dim: str, vertical_units: str, time_units: str, specific_humidity_units: str, lon_dim='lon', lat_dim='lat', time_dim='time', latent_heat_of_condensation=2501000.0)

   Calculate the apparent moisture sink.

   .. math::
       Q_2 = -L \left( \frac{\partial q}{\partial t} + u \frac{\partial q}{\partial x} + v \frac{\partial q}{\partial y} + \omega \frac{\partial q}{\partial p}  \right)

   Parameters
   ----------
   u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The zonal wind data.
   v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The meridional wind data.
   omega_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The vertical velocity data (:math:`\frac{\mathrm{d} p}{\mathrm{d} t}`).
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   time_units: :py:class:`str<python.str>`.
       The unit corresponding to the time dimension value. Optional values are `seconds`, `months`, `years` and so on.
   specific_humidity_units: :py:class:`str<python.str>`.
       The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   time_dim: :py:class:`str<python.str>`, default: `time`.
       The time coordinate dimension name.
   latent_heat_of_condensation: :py:class:`float<python.float>`, default: `2.5008e6` (:math:`\mathrm{J \cdot kg^{-1}}`).
       Latent heat of condensation of water at 0°C.

       .. note::
           - `latent heat - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Latent_heat>`__
           - `Latent heat - Wikipedia <https://en.wikipedia.org/wiki/Latent_heat>`__

   Returns
   -------
   The apparent moisture sink (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::
       - `Yanai, M., & Tomita, T. (1998). Seasonal and Interannual Variability of Atmospheric Heat Sources and Moisture Sinks as Determined from NCEP–NCAR Reanalysis, Journal of Climate, 11(3), 463-482. <https://journals.ametsoc.org/view/journals/clim/11/3/1520-0442_1998_011_0463_saivoa_2.0.co_2.xml>`__
       - `HAO Lisheng, MA Ning, HE Liye. Circulation anomalies characteritics of the abnormal drought and high temperature event in the middle and lower reaches of the Yangtze River in summer of 2022[J]. Arid Meteorology, 2022, 40(5): 721-732 <https://doi.org/10.11755/j.issn.1006-7639(2022)-05-0721>`__


.. py:function:: calc_Plumb_wave_activity_horizontal_flux(z_prime_data: xarray.DataArray, vertical_dim: str, vertical_units: str, lon_dim='lon', lat_dim='lat', omega=7.292e-05, g=9.8, R=6370000)

   Calculate Plumb wave activity horizontal flux.

   Parameters
   ----------
   z_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The anormaly of atmospheric geopotential height.
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   omega: :py:class:`float<python.float>`, default: `7.292e-5`.
       The angular speed of the earth.
   g: :py:class:`float<python.float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float<python.float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The Plumb wave activity horizontal flux (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::
       - `Plumb, R. A., 1985: On the Three-Dimensional Propagation of Stationary Waves. J. Atmos. Sci., 42, 217–229 <https://journals.ametsoc.org/view/journals/atsc/42/3/1520-0469_1985_042_0217_ottdpo_2_0_co_2.xml>`__


.. py:function:: calc_TN_wave_activity_horizontal_flux(z_prime_data, u_climatology_data, v_climatology_data, vertical_dim, vertical_dim_units, lon_dim='lon', lat_dim='lat', omega=7.292e-05, g=9.8, R=6370000)

   Calculate TN wave activity horizontal flux.

   .. math::
       \mathbf{W_h} = \frac{p\cos\varphi}{2\lvert \mathbf{U_c} \rvert}\begin{pmatrix}
                             \frac{U_c}{R^2 \cos^2 \varphi} \left[ \left( \frac{\partial \psi'}{\partial \lambda} \right)^2 - \psi'\frac{\partial^2 \psi'}{\partial \lambda^2} \right] + \frac{V_c}{R^2 \cos \varphi} \left[ \frac{\partial \psi'}{\partial \lambda} \frac{\partial \psi'}{\partial \varphi} - \psi' \frac{\partial^2 \psi'}{\partial \lambda \partial \varphi} \right] \\
                             \frac{U_c}{R^2 \cos \varphi} \left[ \frac{\partial \psi'}{\partial \lambda} \frac{\partial \psi'}{\partial \varphi} - \psi' \frac{\partial^2 \psi'}{\partial \lambda \partial \varphi} \right] + \frac{V_c}{R^2} \left[ \left( \frac{\partial \psi'}{\partial \varphi} \right)^2 - \psi'\frac{\partial^2 \psi'}{\partial \varphi^2} \right] \\
                              \end{pmatrix}    

   Parameters
   ----------
   z_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The anormaly of atmospheric geopotential height.
   u_climatology_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The climatology of zonal wind data.
   v_climatology_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The climatology of meridional wind data.
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   omega: :py:class:`float<python.float>`, default: `7.292e-5`.
       The angular speed of the earth.
   g: :py:class:`float<python.float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float<python.float>`, default: `6370000`.
       Radius of the Earth.

   Returns
   -------
   The TN wave activity horizontal flux (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::
       - http://www.atmos.rcast.u-tokyo.ac.jp/nishii/programs/index.html
       - http://500hpa.cn/pyinmet/tnflux/
       - http://tytd.gx.cn/exchange/tnflux/
       - https://github.com/laishenggx/T-N_Wave-Activity-Flux


.. py:function:: calc_TN_wave_activity_3D_flux(z_prime_data, u_climatology_data, v_climatology_data, temper_data, vertical_dim, vertical_units, z_data=None, lon_dim='lon', lat_dim='lat', omega=7.292e-05, g=9.8, R=6370000, scale_height=8000, kappa=287 / 1005.7, method='practical_height')

   Calculate TN wave activity 3D flux.

   Parameters
   ----------
   z_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The anormaly of atmospheric geopotential height.
   u_climatology_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The climatology of zonal wind data.
   v_climatology_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The climatology of meridional wind data.
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   z_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric geopotential height.
   vertical_dim: :py:class:`str<python.str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str<python.str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   omega: :py:class:`float<python.float>`, default: `7.292e-5`.
       The angular speed of the earth.
   g: :py:class:`float<python.float>`, default: `9.8`.
       The acceleration of gravity.
   R: :py:class:`float<python.float>`, default: `6370000`.
       Radius of the Earth.
   scale_height: :py:class:`float<python.float>`, default: `8000`.
       Scale height.
   kappa: :py:class:`float<python.float>`, default: `287/1005.7`.
       Poisson constant :math:`\kappa`.

       .. note::
           `Poisson constant - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Poisson_constant>`__

   method: :py:class:`str<python.str>`, default: `'practical_height'`.
       The calculation method of :math:`\mathrm{d}z`. Optional values are `'practical_height'`, `'scale_height'`.

   Returns
   -------
   The TN wave activity 3D flux (:py:class:`xarray.DataArray<xarray.DataArray>`).


.. py:function:: calc_EP_horizontal_flux(u_prime_data, v_prime_data, time_dim='time', lat_dim='lat')

   Calculate horizontal Eliassen–Palm Flux.

   Parameters
   ----------
   u_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The anormaly of zonal wind data.
   v_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The anormaly of meridional wind data.
   time_dim: :py:class:`str<python.str>`, default: `time`.
       The time coordinate dimension name.        
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   The Eliassen–Palm Flux (:py:class:`xarray.DataArray<xarray.DataArray>`).

   .. seealso::
       - https://www.ncl.ucar.edu/Applications/EPflux.shtml
       - https://renqlsysu.github.io/2020/02/16/wave_activity_flux/


