easyclimate.field.boundary_layer
================================

.. py:module:: easyclimate.field.boundary_layer


Submodules
----------

.. toctree::
   :maxdepth: 1

   /technical/api/easyclimate/field/boundary_layer/aerobulk/index


Functions
---------

.. autoapisummary::

   easyclimate.field.boundary_layer.calc_turbulent_fluxes_without_skin_correction
   easyclimate.field.boundary_layer.calc_turbulent_fluxes_skin_correction


Package Contents
----------------

.. py:function:: calc_turbulent_fluxes_without_skin_correction(sst_data: xarray.DataArray, sst_data_units: Literal['degC', 'degK', 'degF'], absolute_temperature_data: xarray.DataArray, absolute_temperature_data_units: Literal['degC', 'degK', 'degF'], specific_humidity_data: xarray.DataArray, specific_humidity_data_units: Literal['g/g', 'g/kg', 'kg/kg'], zonal_wind_speed_data: xarray.DataArray, meridional_wind_speed_data: xarray.DataArray, mean_sea_level_pressure_data: xarray.DataArray, mean_sea_level_pressure_data_units: Literal['hPa', 'Pa', 'mbar'], zonal_wind_speed_data_units: Literal['m/s'] = 'm/s', meridional_wind_speed_data_units: Literal['m/s'] = 'm/s', algorithm: Literal['coare3p0', 'coare3p6', 'ecmwf', 'ncar', 'andreas'] = 'coare3p0', height_for_temperature_specific_humidity: float = 2, height_for_wind: float = 10, iteration: int = 8, check_data_valid=True) -> xarray.Dataset

   Aerobulk without skin correction.

   Parameters
   ----------
   sst_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Bulk sea surface temperature.
   sst_data_units: Literal["degC", "degK", "degF"]
       The units of ``sst_data``.
   absolute_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Absolute air temperature at height ``height_for_temperature_specific_humidity``.
   absolute_temperature_data_units: Literal["degC", "degK", "degF"]
       The units of ``absolute_temperature_data``.
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       air humidity at ``height_for_temperature_specific_humidity``, given as specific humidity.
   specific_humidity_data_units: Literal["g/g", "g/kg", "kg/kg"]
       The units of ``specific_humidity_data``.
   zonal_wind_speed_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       zonal wind speed at ``height_for_wind``.
   zonal_wind_speed_data_units: Literal["m/s"]
       The units of ``zonal_wind_speed_data``.
   meridional_wind_speed_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       meridional wind speed at ``height_for_wind``.
   meridional_wind_speed_data_units: Literal["m/s"]
       The units of ``meridional_wind_speed_data``.
   mean_sea_level_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`, optional
       mean sea-level pressure. ~101000 Pa, by default 101000.0.
   mean_sea_level_pressure_data_units: Literal["hPa", "Pa", "mbar"]
       The units of ``mean_sea_level_pressure_data``.
   algorithm: Literal["coare3p0", "coare3p6", "ecmwf", "ncar", "andreas"], default ``coare3p0``.
       Algorithm, can be one of: ``"coare3p0"``, ``"coare3p6"``, ``"ecmwf"``, ``"ncar"``, ``"andreas"``.
   height_for_temperature_specific_humidity: float
       height (:math:`\mathrm{m}`) for temperature and specific humidity of air.
   height_for_wind: float
       height (:math:`\mathrm{m}`) for wind (10m = traditional anemometric height).
   iteration: int
       Number of iteration steps used in the algorithm.

   Returns
   -------
   :py:class:`xarray.Dataset<xarray.Dataset>`.

   * ql: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Latent heat flux (:math:`\mathrm{W/m^2}`).
   * qh: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Sensible heat flux (:math:`\mathrm{W/m^2}`).
   * taux: :py:class:`xarray.DataArray<xarray.DataArray>`.
       zonal wind stress (:math:`\mathrm{N/m^2}`).
   * tauy: :py:class:`xarray.DataArray<xarray.DataArray>`.
       meridional wind stress (:math:`\mathrm{N/m^2}`).
   * evap: :py:class:`xarray.DataArray<xarray.DataArray>`.
       evaporation (:math:`\mathrm{mm/s}`) aka (:math:`\mathrm{kg/m^2/s}`) (usually :math:`< 0`, as ocean loses water).

   .. seealso::

       - https://github.com/brodeau/aerobulk
       - https://github.com/xgcm/aerobulk-python
       - https://ams.confex.com/ams/103ANNUAL/meetingapp.cgi/Session/63444

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_aerobulk.py


.. py:function:: calc_turbulent_fluxes_skin_correction(sst_data: xarray.DataArray, sst_data_units: Literal['degC', 'degK', 'degF'], absolute_temperature_data: xarray.DataArray, absolute_temperature_data_units: Literal['degC', 'degK', 'degF'], specific_humidity_data: xarray.DataArray, specific_humidity_data_units: Literal['g/g', 'g/kg', 'kg/kg'], zonal_wind_speed_data: xarray.DataArray, meridional_wind_speed_data: xarray.DataArray, mean_sea_level_pressure_data: xarray.DataArray, mean_sea_level_pressure_data_units: Literal['hPa', 'Pa', 'mbar'], downwelling_shortwave_radiation: xarray.DataArray, downwelling_shortwave_radiation_units: Literal['W/m^2'], downwelling_longwave_radiation: xarray.DataArray, downwelling_longwave_radiation_units: Literal['W/m^2'], zonal_wind_speed_data_units: Literal['m/s'] = 'm/s', meridional_wind_speed_data_units: Literal['m/s'] = 'm/s', algorithm: Literal['coare3p0', 'coare3p6', 'ecmwf'] = 'coare3p0', height_for_temperature_specific_humidity: float = 2, height_for_wind: float = 10, iteration: int = 8, check_data_valid=True) -> xarray.Dataset

   Aerobulk with skin correction.

   Parameters
   ----------
   sst_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Bulk sea surface temperature.
   sst_data_units: Literal["degC", "degK", "degF"]
       The units of ``sst_data``.
   absolute_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Absolute air temperature at height ``height_for_temperature_specific_humidity``.
   absolute_temperature_data_units: Literal["degC", "degK", "degF"]
       The units of ``absolute_temperature_data``.
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       air humidity at ``height_for_temperature_specific_humidity``, given as specific humidity.
   specific_humidity_data_units: Literal["g/g", "g/kg", "kg/kg"]
       The units of ``specific_humidity_data``.
   zonal_wind_speed_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       zonal wind speed at ``height_for_wind``.
   zonal_wind_speed_data_units: Literal["m/s"]
       The units of ``zonal_wind_speed_data``.
   meridional_wind_speed_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       meridional wind speed at ``height_for_wind``.
   meridional_wind_speed_data_units: Literal["m/s"]
       The units of ``meridional_wind_speed_data``.
   mean_sea_level_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`, optional
       mean sea-level pressure. ~101000 Pa, by default 101000.0.
   mean_sea_level_pressure_data_units: Literal["hPa", "Pa", "mbar"]
       The units of ``mean_sea_level_pressure_data``.
   downwelling_shortwave_radiation: :py:class:`xarray.DataArray<xarray.DataArray>`.
       downwelling shortwave radiation at the surface (>0).
   downwelling_shortwave_radiation_units: Literal["W/m^2"]
       The units of ``downwelling_shortwave_radiation``.
   downwelling_longwave_radiation: :py:class:`xarray.DataArray<xarray.DataArray>`.
       downwelling longwave radiation at the surface (>0).
   downwelling_longwave_radiation_units: Literal["W/m^2"]
       The units of ``downwelling_longwave_radiation``.
   algorithm: Literal["coare3p0", "coare3p6", "ecmwf"], default ``coare3p0``.
       Algorithm, can be one of: ``"coare3p0"``, ``"coare3p6"``, ``"ecmwf"``.
   height_for_temperature_specific_humidity: float
       height (:math:`\mathrm{m}`) for temperature and specific humidity of air.
   height_for_wind: float
       height (:math:`\mathrm{m}`) for wind (10m = traditional anemometric height).
   iteration: int
       Number of iteration steps used in the algorithm.

   Returns
   -------
   :py:class:`xarray.Dataset<xarray.Dataset>`.

   * ql: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Latent heat flux (:math:`\mathrm{W/m^2}`).
   * qh: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Sensible heat flux (:math:`\mathrm{W/m^2}`).
   * taux: :py:class:`xarray.DataArray<xarray.DataArray>`.
       zonal wind stress (:math:`\mathrm{N/m^2}`).
   * tauy: :py:class:`xarray.DataArray<xarray.DataArray>`.
       meridional wind stress (:math:`\mathrm{N/m^2}`).
   * t_s: :py:class:`xarray.DataArray<xarray.DataArray>`.
       skin temperature (:math:`\mathrm{K}`).
   * evap: :py:class:`xarray.DataArray<xarray.DataArray>`.
       evaporation (:math:`\mathrm{mm/s}`) aka (:math:`\mathrm{kg/m^2/s}`) (usually :math:`< 0`, as ocean loses water).

   .. seealso::

       - https://github.com/brodeau/aerobulk
       - https://github.com/xgcm/aerobulk-python
       - https://ams.confex.com/ams/103ANNUAL/meetingapp.cgi/Session/63444

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_aerobulk.py


