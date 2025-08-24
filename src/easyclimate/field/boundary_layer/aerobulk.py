"""
Estimate turbulent air-sea fluxes


**AeroBulk** is a FORTRAN90-based library and suite of tools (including a C++ interface) that feature *state of the art* parameterizations to estimate turbulent air-sea fluxes by means of the traditional **aerodynamic bulk formulae**.

.. raw:: html

    <p align="center">
        <img width="300" src="https://raw.githubusercontent.com/brodeau/aerobulk/refs/heads/master/doc/figs/logo_gray_300.svg">
    </p>

These turbulent fluxes, namely, wind stress, evaporation (latent heat flux) and sensible heat flux,
are estimated using the sea surface temperature (bulk or skin), and the near-surface atmospheric surface state: wind speed,
air temperature and humidity. If the *cool-skin/warm-layer* schemes need to be called to estimate the skin temperature,
surface downwelling shortwave and longwave radiative fluxes are required.

.. seealso::

    - https://github.com/brodeau/aerobulk
    - https://github.com/xgcm/aerobulk-python
    - https://ams.confex.com/ams/103ANNUAL/meetingapp.cgi/Session/63444

Bulk formula and their parameterizations
==================================================
**AeroBulk** relies on the following traditional set of bulk formula to estimate turbulent fluxes at the air-sea interface:

.. math::

    \\vec{\\tau} = \\rho\\  C_D \\  \\vec{U}_z  \\  U_B

.. math::

    E   = \\rho \\ C_E \\     \\big[    q_s   - q_z \\big]  \\  U_B

.. math::

    Q_L = -L_v \\  E

.. math::

    Q_H = \\rho \\ C_H \\ C_P \\ \\big[ \\theta_z - T_s \\big] \\  U_B

where :math:`\\rho` is the density of air. The :math:`z` subscript relates to the reference height above the air-sea interface (generally **z=10m**).
:math:`\vec{U}_z` is the wind speed vector at the reference height.
:math:`U_B` is the bulk scalar wind speed at the reference height (very close to :math:`|\\vec{U}_z|` in most cases).
:math:`\\theta_z` and :math:`q_z` are the potential temperature and specific humidity of air at the reference height, respectively.
:math:`T_s` and :math:`q_s` are the absolute (= potential) temperature and specific humidity of air immediately at the air-sea interface (*z=0*),
respectively; if the *cool-skin/warm-layer* scheme is used, these two are deduced from the skin temperature,
otherwise they are deduced from the bulk SST (default).

Any decent level of accuracy from this set of formula can only be achieved through a consistent estimate of the 3 bulk transfer coefficients:
:math:`C_D`, :math:`C_E`, and :math:`C_H` (drag, evaporation, and sensible heat coefficients).
In **AeroBulk**, these bulk coefficients can be estimated thanks to a collection of *bulk parameterizations* a.k.a *bulk algorithms*,
which relate the value of these coefficients to the near-surface atmospheric stability, the wind speed, and (ideally) the roughness of the sea surface.

The following figure provides a schematic overview on the way turbulent fluxes are computed in AeroBulk:

.. raw:: html

    <p align="center">
        <img width="800" src="https://raw.githubusercontent.com/brodeau/aerobulk/refs/heads/master/doc/figs/fig_bulk_model_f2p.svg">
    </p>

Currently, in AeroBulk, 5 bulk parameterizations are available to compute :math:`C_D`, :math:`C_E`, and :math:`C_H` used in the bulk formula:

*   COARE v3.0 (`Fairall et al., 2003 <https://journals.ametsoc.org/view/journals/clim/16/4/1520-0442_2003_016_0571_bpoasf_2.0.co_2.xml>`__)
*   COARE v3.6 (`Edson et al., 2013 <http://dx.doi.org/10.1175/jpo-d-12-0173.1>`__ + Chris Fairall, *private communication*, 2016)
*   ECMWF (`IFS (Cy40) documentation <https://software.ecmwf.int/wiki/display/IFS/CY40R1+Official+IFS+Documentation>`__)
*   ANDREAS (`Andreas et al., 2015 <https://dx.doi.org/10.1002/qj.2424>`__)
*   NCAR (Large & Yeager 2004, `2009 <http://dx.doi.org/10.1007/s00382-008-0441-3>`__)

In the COARE and ECMWF algorithms, a cool-skin/warm layer scheme is included and can be activated if the input sea-surface temperature to be used is the bulk SST (usually measured a few tenths of meters below the surface). Activation of these cool-skin/warm layer schemes requires the surface downwelling shortwave and longwave radiative flux components to be provided. Other parameterizations, such as NCAR, are meant to be used with the bulk SST, and do not feature a cool-skin/warm layer scheme.

Beside bulk algorithms, AeroBulk also provides a collection of functions (module `mod_phymbl.f90`) to accurately estimate relevant atmospheric parameters such as: density of air, different expressions of the humidity of air, viscosity of air, specific humidity at saturation, *Obukhov* length, bulk *Richardson* number, wind gustiness, etc...

The focus in AeroBulk is readability, efficiency, and portability towards modern ocean & atmosphere GCMs (Fortran 90, set of modules and a library).

.. raw:: html

    <p align="center">
        <img width="800" src="https://raw.githubusercontent.com/brodeau/aerobulk/refs/heads/master/doc/figs/Comparaison_Psi.svg">
    </p>

*Fig 1/ Comparison of the stability correction profiles :math:`\\Psi(\\zeta)` as used in 4 different bulk algorithms.*

.. raw:: html

    <p align="center">
        <img width="800" src="https://raw.githubusercontent.com/brodeau/aerobulk/refs/heads/master/doc/figs/Comparaison_CxN10.svg">
    </p>

*Fig 2/ Comparison of the neutral drag (thick lines) and evaporation coefficients (thinner lines) as a function of the neutral wind speed at 10m.*
"""

from ...backend import aeronoskin, aeroskin
import numpy as np
import xarray as xr
import warnings
from typing import Literal
from ...core.utility import (
    transfer_data_temperature_units,
    transfer_data_multiple_units,
    transfer_data_difference_units,
)

__all__ = [
    "calc_turbulent_fluxes_without_skin_correction",
    "calc_turbulent_fluxes_skin_correction",
]

VALID_VALUE_RANGES = {
    "sst": np.array([270, 320]),
    "t_zt": np.array([180, 330]),
    "hum_zt": np.array([0, 0.08]),
    "u_zu": np.array([-50, 50]),
    "v_zu": np.array([-50, 50]),
    "slp": np.array([80000, 110000]),
    "rad_sw": np.array([0, 1500]),
    "rad_lw": np.array([0, 750]),
}


def check_value_range(var_name: str, var_value: xr.DataArray):
    """
    检查变量值是否在有效范围内。

    参数:
        var_name (str): 变量名，例如 "sst"
        var_value (float/int): 变量的值

    异常:
        KeyError: 如果变量名不在 VALID_VALUE_RANGES 中
        ValueError: 如果变量值不在有效范围内
    """
    # 检查变量名是否在字典中
    if var_name not in VALID_VALUE_RANGES:
        warnings.warn(f"The varable {var_name} don't have valid range.", UserWarning)
        return

    # 获取变量的有效范围
    range_min, range_max = VALID_VALUE_RANGES[var_name]
    try:
        var_value = var_value.compute()
    except:
        var_value = var_value

    # 检查值是否在范围内
    if not (range_min <= var_value <= range_max):
        warnings.warn(
            f"The varable {var_name} of values {var_value} is beyond [{range_min}, {range_max}]。",
            UserWarning,
        )


def calc_turbulent_fluxes_without_skin_correction(
    sst_data: xr.DataArray,
    sst_data_units: Literal["degC", "degK", "degF"],
    absolute_temperature_data: xr.DataArray,
    absolute_temperature_data_units: Literal["degC", "degK", "degF"],
    specific_humidity_data: xr.DataArray,
    specific_humidity_data_units: Literal["g/g", "g/kg", "kg/kg"],
    zonal_wind_speed_data: xr.DataArray,
    meridional_wind_speed_data: xr.DataArray,
    mean_sea_level_pressure_data: xr.DataArray,
    mean_sea_level_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    zonal_wind_speed_data_units: Literal["m/s"] = "m/s",
    meridional_wind_speed_data_units: Literal["m/s"] = "m/s",
    algorithm: Literal["coare3p0", "coare3p6", "ecmwf", "ncar", "andreas"] = "coare3p0",
    height_for_temperature_specific_humidity: float = 2,
    height_for_wind: float = 10,
    iteration: int = 8,
    check_data_valid=True,
) -> xr.Dataset:
    """
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
        height (:math:`\\mathrm{m}`) for temperature and specific humidity of air.
    height_for_wind: float
        height (:math:`\\mathrm{m}`) for wind (10m = traditional anemometric height).
    iteration: int
        Number of iteration steps used in the algorithm.

    Returns
    -------
    :py:class:`xarray.Dataset<xarray.Dataset>`.

    * ql: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Latent heat flux (:math:`\\mathrm{W/m^2}`).
    * qh: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Sensible heat flux (:math:`\\mathrm{W/m^2}`).
    * taux: :py:class:`xarray.DataArray<xarray.DataArray>`.
        zonal wind stress (:math:`\\mathrm{N/m^2}`).
    * tauy: :py:class:`xarray.DataArray<xarray.DataArray>`.
        meridional wind stress (:math:`\\mathrm{N/m^2}`).
    * evap: :py:class:`xarray.DataArray<xarray.DataArray>`.
        evaporation (:math:`\\mathrm{mm/s}`) aka (:math:`\\mathrm{kg/m^2/s}`) (usually :math:`< 0`, as ocean loses water).

    .. seealso::

        - https://github.com/brodeau/aerobulk
        - https://github.com/xgcm/aerobulk-python
        - https://ams.confex.com/ams/103ANNUAL/meetingapp.cgi/Session/63444

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_aerobulk.py
    """
    sst_data = transfer_data_temperature_units(sst_data, sst_data_units, "degK")
    absolute_temperature_data = transfer_data_temperature_units(
        absolute_temperature_data, absolute_temperature_data_units, "degK"
    )
    specific_humidity_data = transfer_data_multiple_units(
        specific_humidity_data, specific_humidity_data_units, "kg/kg"
    )
    zonal_wind_speed_data = transfer_data_multiple_units(
        zonal_wind_speed_data, zonal_wind_speed_data_units, "m/s"
    )
    meridional_wind_speed_data = transfer_data_multiple_units(
        meridional_wind_speed_data, meridional_wind_speed_data_units, "m/s"
    )
    mean_sea_level_pressure_data = transfer_data_multiple_units(
        mean_sea_level_pressure_data, mean_sea_level_pressure_data_units, "Pa"
    )

    if check_data_valid:
        check_value_range("sst", sst_data.max().data)
        check_value_range("sst", sst_data.min().data)
        check_value_range("t_zt", absolute_temperature_data.max().data)
        check_value_range("t_zt", absolute_temperature_data.min().data)
        check_value_range("hum_zt", specific_humidity_data.max().data)
        check_value_range("hum_zt", specific_humidity_data.min().data)
        check_value_range("u_zu", zonal_wind_speed_data.max().data)
        check_value_range("u_zu", zonal_wind_speed_data.min().data)
        check_value_range("v_zu", meridional_wind_speed_data.max().data)
        check_value_range("v_zu", meridional_wind_speed_data.min().data)
        check_value_range("slp", mean_sea_level_pressure_data.max().data)
        check_value_range("slp", mean_sea_level_pressure_data.min().data)

    def _noskin(sst, t_zt, hum_zt, u_zu, v_zu, slp):
        arrays = np.array([sst, t_zt, hum_zt, u_zu, v_zu, slp])
        if np.isnan(arrays).any():
            return np.array([np.nan, np.nan, np.nan, np.nan, np.nan])

        sst = arrays[0]
        t_zt = arrays[1]
        hum_zt = arrays[2]
        u_zu = arrays[3]
        v_zu = arrays[4]
        slp = arrays[5]

        args_shrunk = (
            sst.reshape(1, 1, 1),
            t_zt.reshape(1, 1, 1),
            hum_zt.reshape(1, 1, 1),
            u_zu.reshape(1, 1, 1),
            v_zu.reshape(1, 1, 1),
            slp.reshape(1, 1, 1),
        )

        ql, qh, taux, tauy, evap = (
            aeronoskin.mod_aerobulk_wrapper_noskin.aerobulk_model_noskin(
                algorithm,
                height_for_temperature_specific_humidity,
                height_for_wind,
                *args_shrunk,
                iteration,
            )
        )

        return np.array(
            [ql[0][0][0], qh[0][0][0], taux[0][0][0], tauy[0][0][0], evap[0][0][0]]
        )

    tmp_result = xr.apply_ufunc(
        _noskin,
        sst_data,
        absolute_temperature_data,
        specific_humidity_data,
        zonal_wind_speed_data,
        meridional_wind_speed_data,
        mean_sea_level_pressure_data,
        output_core_dims=[["parameter"]],
        output_dtypes=["float32"],
        vectorize=True,
        dask="parallelized",
        dask_gufunc_kwargs={"output_sizes": {"parameter": 5}, "allow_rechunk": True},
    )

    result1 = tmp_result[..., 0]
    result2 = tmp_result[..., 1]
    result3 = tmp_result[..., 2]
    result4 = tmp_result[..., 3]
    result5 = tmp_result[..., 4]

    result1.attrs["long_name"] = "Latent heat flux"
    result1.attrs["units"] = "W/m^2"
    result2.attrs["long_name"] = "Sensible heat flux"
    result2.attrs["units"] = "W/m^2"
    result3.attrs["long_name"] = "Zonal wind stress"
    result3.attrs["units"] = "W/m^2"
    result4.attrs["long_name"] = "Meridional wind stress"
    result4.attrs["units"] = "W/m^2"
    result5.attrs["long_name"] = "Evaporation (usually <0, as ocean loses water!)"
    result5.attrs["units"] = "[mm/s] aka [kg/m^2/s]"

    result = xr.Dataset(
        {
            "ql": result1,
            "qh": result2,
            "taux": result3,
            "tauy": result4,
            "evap": result5,
        }
    )
    result.attrs["method"] = "aerobulk without skin correction"
    result.attrs["algorithm"] = algorithm
    result.attrs["height_for_temperature_specific_humidity"] = (
        height_for_temperature_specific_humidity
    )
    result.attrs["height_for_wind"] = height_for_wind
    result.attrs["iteration"] = iteration
    return result


def calc_turbulent_fluxes_skin_correction(
    sst_data: xr.DataArray,
    sst_data_units: Literal["degC", "degK", "degF"],
    absolute_temperature_data: xr.DataArray,
    absolute_temperature_data_units: Literal["degC", "degK", "degF"],
    specific_humidity_data: xr.DataArray,
    specific_humidity_data_units: Literal["g/g", "g/kg", "kg/kg"],
    zonal_wind_speed_data: xr.DataArray,
    meridional_wind_speed_data: xr.DataArray,
    mean_sea_level_pressure_data: xr.DataArray,
    mean_sea_level_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    downwelling_shortwave_radiation: xr.DataArray,
    downwelling_shortwave_radiation_units: Literal["W/m^2"],
    downwelling_longwave_radiation: xr.DataArray,
    downwelling_longwave_radiation_units: Literal["W/m^2"],
    zonal_wind_speed_data_units: Literal["m/s"] = "m/s",
    meridional_wind_speed_data_units: Literal["m/s"] = "m/s",
    algorithm: Literal["coare3p0", "coare3p6", "ecmwf"] = "coare3p0",
    height_for_temperature_specific_humidity: float = 2,
    height_for_wind: float = 10,
    iteration: int = 8,
    check_data_valid=True,
) -> xr.Dataset:
    """
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
        height (:math:`\\mathrm{m}`) for temperature and specific humidity of air.
    height_for_wind: float
        height (:math:`\\mathrm{m}`) for wind (10m = traditional anemometric height).
    iteration: int
        Number of iteration steps used in the algorithm.

    Returns
    -------
    :py:class:`xarray.Dataset<xarray.Dataset>`.

    * ql: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Latent heat flux (:math:`\\mathrm{W/m^2}`).
    * qh: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Sensible heat flux (:math:`\\mathrm{W/m^2}`).
    * taux: :py:class:`xarray.DataArray<xarray.DataArray>`.
        zonal wind stress (:math:`\\mathrm{N/m^2}`).
    * tauy: :py:class:`xarray.DataArray<xarray.DataArray>`.
        meridional wind stress (:math:`\\mathrm{N/m^2}`).
    * t_s: :py:class:`xarray.DataArray<xarray.DataArray>`.
        skin temperature (:math:`\\mathrm{K}`).
    * evap: :py:class:`xarray.DataArray<xarray.DataArray>`.
        evaporation (:math:`\\mathrm{mm/s}`) aka (:math:`\\mathrm{kg/m^2/s}`) (usually :math:`< 0`, as ocean loses water).

    .. seealso::

        - https://github.com/brodeau/aerobulk
        - https://github.com/xgcm/aerobulk-python
        - https://ams.confex.com/ams/103ANNUAL/meetingapp.cgi/Session/63444

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_aerobulk.py
    """
    sst_data = transfer_data_temperature_units(sst_data, sst_data_units, "degK")
    absolute_temperature_data = transfer_data_temperature_units(
        absolute_temperature_data, absolute_temperature_data_units, "degK"
    )
    specific_humidity_data = transfer_data_multiple_units(
        specific_humidity_data, specific_humidity_data_units, "kg/kg"
    )
    zonal_wind_speed_data = transfer_data_multiple_units(
        zonal_wind_speed_data, zonal_wind_speed_data_units, "m/s"
    )
    meridional_wind_speed_data = transfer_data_multiple_units(
        meridional_wind_speed_data, meridional_wind_speed_data_units, "m/s"
    )
    mean_sea_level_pressure_data = transfer_data_multiple_units(
        mean_sea_level_pressure_data, mean_sea_level_pressure_data_units, "Pa"
    )
    downwelling_shortwave_radiation = transfer_data_multiple_units(
        downwelling_shortwave_radiation, downwelling_shortwave_radiation_units, "W/m^2"
    )
    downwelling_longwave_radiation = transfer_data_multiple_units(
        downwelling_longwave_radiation, downwelling_longwave_radiation_units, "W/m^2"
    )

    if check_data_valid:
        check_value_range("sst", sst_data.max().data)
        check_value_range("sst", sst_data.min().data)
        check_value_range("t_zt", absolute_temperature_data.max().data)
        check_value_range("t_zt", absolute_temperature_data.min().data)
        check_value_range("hum_zt", specific_humidity_data.max().data)
        check_value_range("hum_zt", specific_humidity_data.min().data)
        check_value_range("u_zu", zonal_wind_speed_data.max().data)
        check_value_range("u_zu", zonal_wind_speed_data.min().data)
        check_value_range("v_zu", meridional_wind_speed_data.max().data)
        check_value_range("v_zu", meridional_wind_speed_data.min().data)
        check_value_range("slp", mean_sea_level_pressure_data.max().data)
        check_value_range("slp", mean_sea_level_pressure_data.min().data)
        check_value_range("rad_sw", downwelling_shortwave_radiation.max().data)
        check_value_range("rad_sw", downwelling_shortwave_radiation.min().data)
        check_value_range("rad_lw", downwelling_longwave_radiation.max().data)
        check_value_range("rad_lw", downwelling_longwave_radiation.min().data)

    def _skin(sst, t_zt, hum_zt, u_zu, v_zu, slp, rad_sw, rad_lw):
        arrays = np.array([sst, t_zt, hum_zt, u_zu, v_zu, slp, rad_sw, rad_lw])
        if np.isnan(arrays).any():
            return np.array([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])

        sst = arrays[0]
        t_zt = arrays[1]
        hum_zt = arrays[2]
        u_zu = arrays[3]
        v_zu = arrays[4]
        slp = arrays[5]
        rad_sw = arrays[6]
        rad_lw = arrays[7]

        args_shrunk = (
            sst.reshape(1, 1, 1),
            t_zt.reshape(1, 1, 1),
            hum_zt.reshape(1, 1, 1),
            u_zu.reshape(1, 1, 1),
            v_zu.reshape(1, 1, 1),
            slp.reshape(1, 1, 1),
            rad_sw.reshape(1, 1, 1),
            rad_lw.reshape(1, 1, 1),
        )

        ql, qh, taux, tauy, t_s, evap = (
            aeroskin.mod_aerobulk_wrapper_skin.aerobulk_model_skin(
                algorithm,
                height_for_temperature_specific_humidity,
                height_for_wind,
                *args_shrunk,
                iteration,
            )
        )

        return np.array(
            [
                ql[0][0][0],
                qh[0][0][0],
                taux[0][0][0],
                tauy[0][0][0],
                t_s[0][0][0],
                evap[0][0][0],
            ]
        )

    tmp_result = xr.apply_ufunc(
        _skin,
        sst_data,
        absolute_temperature_data,
        specific_humidity_data,
        zonal_wind_speed_data,
        meridional_wind_speed_data,
        mean_sea_level_pressure_data,
        downwelling_shortwave_radiation,
        downwelling_longwave_radiation,
        output_core_dims=[["parameter"]],
        output_dtypes=["float64"],
        vectorize=True,
        dask="parallelized",
        dask_gufunc_kwargs={"output_sizes": {"parameter": 6}, "allow_rechunk": True},
    )

    result1 = tmp_result[..., 0]
    result2 = tmp_result[..., 1]
    result3 = tmp_result[..., 2]
    result4 = tmp_result[..., 3]
    result5 = tmp_result[..., 4]
    result6 = tmp_result[..., 5]

    result1.attrs["long_name"] = "Latent heat flux"
    result1.attrs["units"] = "W/m^2"
    result2.attrs["long_name"] = "Sensible heat flux"
    result2.attrs["units"] = "W/m^2"
    result3.attrs["long_name"] = "Zonal wind stress"
    result3.attrs["units"] = "W/m^2"
    result4.attrs["long_name"] = "Meridional wind stress"
    result4.attrs["units"] = "W/m^2"
    result5.attrs["long_name"] = "skin temperature"
    result5.attrs["units"] = "K"
    result6.attrs["long_name"] = "Evaporation (usually <0, as ocean loses water!)"
    result6.attrs["units"] = "[mm/s] aka [kg/m^2/s]"

    result = xr.Dataset(
        {
            "ql": result1,
            "qh": result2,
            "taux": result3,
            "tauy": result4,
            "t_s": result5,
            "evap": result6,
        }
    )
    result.attrs["method"] = "aerobulk with skin correction"
    result.attrs["algorithm"] = algorithm
    result.attrs["height_for_temperature_specific_humidity"] = (
        height_for_temperature_specific_humidity
    )
    result.attrs["height_for_wind"] = height_for_wind
    result.attrs["iteration"] = iteration
    return result
