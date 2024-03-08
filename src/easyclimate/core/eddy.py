"""
The calculation of transient eddy
"""

from __future__ import annotations
import xarray as xr
import numpy as np
from .diagnosis import (
    get_coriolis_parameter,
    calc_potential_temperature,
    calc_brunt_vaisala_frequency_atm,
)
from .diff import (
    calc_gradient,
    calc_lon_gradient,
    calc_lat_gradient,
    calc_lon_laplacian,
    calc_lon_lat_mixed_derivatives,
)
from .utility import (
    transfer_data_units,
    transfer_units_coeff,
    transfer_deg2rad,
    transfer_dFdp2dFdz,
)

__all__ = [
    "calc_eady_growth_rate",
    "calc_apparent_heat_source",
    "calc_total_diabatic_heating",
    "calc_apparent_moisture_sink",
    "calc_Plumb_wave_activity_horizontal_flux",
    "calc_TN_wave_activity_horizontal_flux",
    "calc_EP_horizontal_flux",
]


def calc_eady_growth_rate(
    u_daily_data: xr.DataArray,
    z_daily_data: xr.DataArray,
    temper_daily_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: str,
    lat_dim="lat",
    g=9.8,
) -> xr.Dataset:
    """
    Calculate the maximum Eady growth rate.

    .. math::
        \\sigma = 0.3098 \\frac{f}{N} \\frac{\mathrm{d} U}{\mathrm{d} z}

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
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    g: :py:class:`float <float>`, default: `9.8`.
        The acceleration of gravity.

    Returns
    -------
    The maximum Eady growth rate. (:py:class:`xarray.Dataset<xarray.Dataset>`).

    - `eady_growth_rate`: The maximum Eady growth rate.
    - `dudz`: :math:`\\frac{\\mathrm{d} U}{\\mathrm{d} z}`
    - `brunt_vaisala_frequency`: Brunt-väisälä frequency.

    .. seealso::
        - `eady_growth_rate -NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/eady_growth_rate.shtml>`__
        - `瞬变涡旋诊断量 <https://renqlsysu.github.io/2020/02/16/wave_activity_flux/>`__
    """
    dim_tuple = u_daily_data.dims
    f = get_coriolis_parameter(u_daily_data[lat_dim])
    dp = 1
    dudp = calc_gradient(u_daily_data, dim=vertical_dim) / dp
    dzdp = calc_gradient(z_daily_data, dim=vertical_dim) / dp
    dudz = dudp / dzdp
    pt = calc_potential_temperature(
        temper_daily_data,
        vertical_dim=vertical_dim,
        vertical_dim_units=vertical_dim_units,
    )
    brunt_vaisala_atm = calc_brunt_vaisala_frequency_atm(
        pt, z_daily_data, vertical_dim=vertical_dim
    )
    eady_growth_rate = 0.3098 * g * np.abs(f) * np.abs(dudz) / brunt_vaisala_atm
    eady_growth_rate = eady_growth_rate.transpose(*dim_tuple)
    result_dataset = xr.Dataset(
        data_vars={
            "eady_growth_rate": eady_growth_rate,
            "dudz": dudz,
            "brunt_vaisala_frequency": brunt_vaisala_atm,
        }
    )
    return result_dataset


def calc_apparent_heat_source(
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    omega_data: xr.DataArray,
    temper_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: str,
    time_units: str,
    lon_dim="lon",
    lat_dim="lat",
    time_dim="time",
    c_p=1005.7,
) -> xr.DataArray:
    """
    Calculate the apparent heat source.

    .. math::
        Q_1 = C_p \\frac{T}{\\theta} \\left( \\frac{\\partial \\theta}{\\partial t} + u \\frac{\\partial \\theta}{\\partial x} + v \\frac{\\partial \\theta}{\\partial y} + \\omega \\frac{\\partial \\theta}{\\partial p} \\right)

    Parameters
    ----------
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    omega_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The vertical velocity data (:math:`\\frac{\\mathrm{d} p}{\\mathrm{d} t}`).
    temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Air temperature.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    time_units: :py:class:`str <str>`.
        The unit corresponding to the time dimension value. Optional values are `seconds`, `months`, `years` and so on.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    time_dim: :py:class:`str <str>`.
        The time coordinate dimension name.
    c_p: :py:class:`float <float>`, default: `1005.7`.
        The specific heat at constant pressure of dry air.

        .. note::
            `specific heat capacity - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Specific_heat_capacity>`__

    Returns
    -------
    The apparent heat source (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. seealso::
        - `Yanai, M., & Tomita, T. (1998). Seasonal and Interannual Variability of Atmospheric Heat Sources and Moisture Sinks as Determined from NCEP–NCAR Reanalysis, Journal of Climate, 11(3), 463-482. <https://journals.ametsoc.org/view/journals/clim/11/3/1520-0442_1998_011_0463_saivoa_2.0.co_2.xml>`__
        - `Ling, J., & Zhang, C. (2013). Diabatic Heating Profiles in Recent Global Reanalyses, Journal of Climate, 26(10), 3307-3325. <https://doi.org/10.1175/JCLI-D-12-00384.1>`__
    """
    # Convert time units to seconds
    dt = transfer_units_coeff(time_units, "seconds")
    # Convert the pressure unit to Pascal
    dp_base = transfer_units_coeff(vertical_dim_units, "Pa")

    pt = calc_potential_temperature(
        temper_data, vertical_dim=vertical_dim, vertical_dim_units=vertical_dim_units
    )
    dtheta_dt = calc_gradient(pt, dim=time_dim) / dt
    dtheta_dx = calc_lon_gradient(pt, lon_dim=lon_dim, lat_dim=lat_dim)
    dtheta_dy = calc_lat_gradient(pt, lat_dim=lat_dim)

    dp = calc_gradient(pt[vertical_dim], dim=vertical_dim) * dp_base
    domega_dp = calc_gradient(pt, dim=vertical_dim) / dp

    Q1 = (
        c_p
        * (temper_data / pt)
        * (dtheta_dt + u_data * dtheta_dx + v_data * dtheta_dy + omega_data * domega_dp)
    )
    Q1.name = "apparent_heat_source"
    return Q1


def calc_total_diabatic_heating(
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    omega_data: xr.DataArray,
    temper_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: str,
    time_units: str,
    lat_dim="lat",
    lon_dim="lon",
    time_dim="time",
    c_p=1005.7,
) -> xr.DataArray:
    """
    Calculate the total diabatic heating.

    Calculated in exactly the same way as for the apparent heat source.

    Parameters
    ----------
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    omega_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The vertical velocity data (:math:`\\frac{\\mathrm{d} p}{\\mathrm{d} t}`).
    temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Air temperature.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    time_units: :py:class:`str <str>`.
        The unit corresponding to the time dimension value. Optional values are `seconds`, `months`, `years` and so on.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.
    c_p: :py:class:`float <float>`, default: `1005.7` (:math:`\\mathrm{J \\cdot kg^{-1} \\cdot K^{-1}}`).
        The specific heat at constant pressure of dry air.

        .. note::
            `specific heat capacity - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Specific_heat_capacity>`__

    Returns
    -------
    The total diabatic heating (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. seealso::
        :py:func:`calc_apparent_heat_source <calc_apparent_heat_source>`
    """
    Q1 = calc_apparent_heat_source(
        u_data=u_data,
        v_data=v_data,
        omega_data=omega_data,
        temper_data=temper_data,
        vertical_dim=vertical_dim,
        vertical_dim_units=vertical_dim_units,
        time_units=time_units,
        lat_dim=lat_dim,
        lon_dim=lon_dim,
        time_dim=time_dim,
        c_p=c_p,
    )
    Q1.name = "total_diabatic_heating"
    return Q1


def calc_apparent_moisture_sink(
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    omega_data: xr.DataArray,
    specific_humidity_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: str,
    time_units: str,
    specific_humidity_units: str,
    lon_dim="lon",
    lat_dim="lat",
    time_dim="time",
    latent_heat_of_condensation=2.501e6,
) -> xr.DataArray:
    """
    Calculate the apparent moisture sink.

    .. math::
        Q_2 = -L \\left( \\frac{\\partial q}{\\partial t} + u \\frac{\\partial q}{\\partial x} + v \\frac{\\partial q}{\\partial y} + \\omega \\frac{\\partial q}{\\partial p}  \\right)

    Parameters
    ----------
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    omega_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The vertical velocity data (:math:`\\frac{\\mathrm{d} p}{\\mathrm{d} t}`).
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    time_units: :py:class:`str <str>`.
        The unit corresponding to the time dimension value. Optional values are `seconds`, `months`, `years` and so on.
    specific_humidity_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.
    latent_heat_of_condensation: :py:class:`float <float>`, default: `2.5008e6` (:math:`\\mathrm{J \\cdot kg^{-1}}`).
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
    """
    specific_humidity_data = transfer_data_units(
        specific_humidity_data, specific_humidity_units, "kg/kg"
    )

    # Convert time units to seconds
    dt = transfer_units_coeff(time_units, "seconds")
    # Convert the pressure unit to Pascal
    dp_base = transfer_units_coeff(vertical_dim_units, "Pa")

    dqs_dt = calc_gradient(specific_humidity_data, dim=time_dim) / dt
    dqs_dx = calc_lon_gradient(specific_humidity_data, lon_dim=lon_dim, lat_dim=lat_dim)
    dqs_dy = calc_lat_gradient(specific_humidity_data, lat_dim=lat_dim)

    dp = calc_gradient(specific_humidity_data[vertical_dim], dim=vertical_dim) * dp_base
    domega_dp = calc_gradient(specific_humidity_data, dim=vertical_dim) / dp

    Q2 = (
        (-1)
        * latent_heat_of_condensation
        * (dqs_dt + u_data * dqs_dx + v_data * dqs_dy + omega_data * domega_dp)
    )
    Q2.name = "apparent_moisture_sink"
    return Q2


def calc_Plumb_wave_activity_horizontal_flux(
    z_prime_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: str,
    lon_dim="lon",
    lat_dim="lat",
    omega=7.292e-5,
    g=9.8,
    R=6370000,
) -> xr.Dataset:
    """
    Calculate Plumb wave activity horizontal flux.

    Parameters
    ----------
    z_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The anormaly of atmospheric geopotential height.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
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
    The Plumb wave activity horizontal flux (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. seealso::
        - `Plumb, R. A., 1985: On the Three-Dimensional Propagation of Stationary Waves. J. Atmos. Sci., 42, 217–229 <https://journals.ametsoc.org/view/journals/atsc/42/3/1520-0469_1985_042_0217_ottdpo_2_0_co_2.xml>`__
    """
    coordinate_sample_data = z_prime_data

    lat_array = coordinate_sample_data["lat"].astype("float64")
    coslat = np.cos(transfer_deg2rad(lat_array))

    f = get_coriolis_parameter(lat_array, omega=omega)
    psi_p = z_prime_data * g / f

    p_lev = transfer_data_units(
        coordinate_sample_data[vertical_dim], vertical_dim_units, "Pa"
    )
    p = p_lev / 1e5

    dpsi_dlambda = calc_gradient(psi_p, dim=lon_dim)
    dpsi_dphi = calc_gradient(psi_p, dim=lat_dim)

    d2psi_dlambda2 = calc_gradient(dpsi_dlambda, dim=lon_dim)
    d2psi_dlambdadphi = calc_gradient(dpsi_dphi, dim=lon_dim)

    term_xu = dpsi_dlambda**2 - (psi_p * d2psi_dlambda2)
    term_xv = dpsi_dlambda * dpsi_dphi - (psi_p * d2psi_dlambdadphi)

    fx = p * ((1 / (2 * R**2 * coslat)) * term_xu)
    fy = p * ((1 / (2 * R**2)) * term_xv)

    result = xr.Dataset()
    result["psi_p"] = psi_p
    result["fx"] = fx
    result["fy"] = fy
    return result


def calc_TN_wave_activity_horizontal_flux(
    z_prime_data: xr.DataArray,
    u_climatology_data: xr.DataArray,
    v_climatology_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: str,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    omega: float = 7.292e-5,
    g: float = 9.8,
    R: float = 6370000,
) -> xr.DataArray:
    """
    Calculate TN wave activity horizontal flux.

    .. math::
        \\mathbf{W_h} = \\frac{p\\cos\\varphi}{2\\lvert \\mathbf{U_c} \\rvert}\\begin{pmatrix}
                              \\frac{U_c}{R^2 \\cos^2 \\varphi} \\left[ \\left( \\frac{\\partial \\psi'}{\\partial \\lambda} \\right)^2 - \\psi'\\frac{\\partial^2 \\psi'}{\\partial \\lambda^2} \\right] + \\frac{V_c}{R^2 \\cos \\varphi} \\left[ \\frac{\\partial \\psi'}{\\partial \\lambda} \\frac{\\partial \\psi'}{\\partial \\varphi} - \\psi' \\frac{\\partial^2 \\psi'}{\\partial \\lambda \\partial \\varphi} \\right] \\\\
                              \\frac{U_c}{R^2 \\cos \\varphi} \\left[ \\frac{\\partial \\psi'}{\\partial \\lambda} \\frac{\\partial \\psi'}{\\partial \\varphi} - \\psi' \\frac{\\partial^2 \\psi'}{\\partial \\lambda \\partial \\varphi} \\right] + \\frac{V_c}{R^2} \\left[ \\left( \\frac{\\partial \\psi'}{\\partial \\varphi} \\right)^2 - \\psi'\\frac{\\partial^2 \\psi'}{\\partial \\varphi^2} \\right] \\\\
                               \\end{pmatrix}

    Parameters
    ----------
    z_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The anormaly of atmospheric geopotential height.
    u_climatology_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The climatology of zonal wind data.
    v_climatology_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The climatology of meridional wind data.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
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
    The TN wave activity horizontal flux (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. seealso::
        - http://www.atmos.rcast.u-tokyo.ac.jp/nishii/programs/index.html
        - http://500hpa.cn/pyinmet/tnflux/
        - http://tytd.gx.cn/exchange/tnflux/
        - https://github.com/laishenggx/T-N_Wave-Activity-Flux
    """
    coordinate_sample_data = z_prime_data
    u_c = u_climatology_data
    v_c = v_climatology_data

    lat_array = coordinate_sample_data["lat"].astype("float64")
    coslat = np.cos(transfer_deg2rad(lat_array))

    f = get_coriolis_parameter(lat_array, omega=omega)
    psi_p = z_prime_data * g / f

    p_lev = transfer_data_units(
        coordinate_sample_data[vertical_dim], vertical_dim_units, "Pa"
    )
    p = p_lev / 1e5

    dpsi_dlambda = calc_gradient(psi_p, dim=lon_dim)
    dpsi_dphi = calc_gradient(psi_p, dim=lat_dim)

    d2psi_dlambda2 = calc_gradient(dpsi_dlambda, dim=lon_dim)
    d2psi_dphi2 = calc_gradient(dpsi_dphi, dim=lat_dim)
    d2psi_dlambdadphi = calc_gradient(dpsi_dphi, dim=lon_dim)

    term_xu = dpsi_dlambda**2 - (psi_p * d2psi_dlambda2)
    term_xv = dpsi_dlambda * dpsi_dphi - (psi_p * d2psi_dlambdadphi)
    term_yv = dpsi_dphi**2 - (psi_p * d2psi_dphi2)

    magU = np.sqrt(u_c**2 + v_c**2)
    coeff = p / (2 * magU)

    fx = coeff * ((u_c / (R**2 * coslat) * term_xu) + (v_c / (R**2) * term_xv))
    fy = coeff * ((u_c / (R**2) * term_xv) + (v_c * coslat / (R**2) * term_yv))

    result = xr.Dataset()
    result["psi_p"] = psi_p
    result["fx"] = fx
    result["fy"] = fy
    return result


# def calc_TN_wave_activity_3D_flux(
#     z_prime_data: xr.DataArray,
#     u_climatology_data: xr.DataArray,
#     v_climatology_data: xr.DataArray,
#     temper_data: xr.DataArray,
#     vertical_dim: str,
#     vertical_dim_units: str,
#     z_data: xr.DataArray | None = None,
#     lon_dim: str = 'lon',
#     lat_dim: str = 'lat',
#     omega: float = 7.292e-5,
#     g: float = 9.8,
#     R: float = 6370000,
#     scale_height: float = 8000.,
#     kappa: float = 287/1005.7,
#     method: str ='practical_height',
# ) -> xr.Dataset:
#     """
#     Calculate TN wave activity 3D flux.

#     Parameters
#     ----------
#     z_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
#         The anormaly of atmospheric geopotential height.
#     u_climatology_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
#         The climatology of zonal wind data.
#     v_climatology_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
#         The climatology of meridional wind data.
#     temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
#         Air temperature.
#     z_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
#         Atmospheric geopotential height.
#     vertical_dim: :py:class:`str <str>`.
#         Vertical coordinate dimension name.
#     vertical_dim_units: :py:class:`str <str>`.
#         The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
#     lon_dim: :py:class:`str <str>`, default: `lon`.
#         Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
#     lat_dim: :py:class:`str <str>`, default: `lat`.
#         Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
#     omega: :py:class:`float <float>`, default: `7.292e-5`.
#         The angular speed of the earth.
#     g: :py:class:`float <float>`, default: `9.8`.
#         The acceleration of gravity.
#     R: :py:class:`float <float>`, default: `6370000`.
#         Radius of the Earth.
#     scale_height: :py:class:`float <float>`, default: `8000`.
#         Scale height.
#     kappa: :py:class:`float <float>`, default: `287/1005.7`.
#         Poisson constant :math:`\\kappa`.

#         .. note::
#             `Poisson constant - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Poisson_constant>`__

#     method: :py:class:`str <str>`, default: `'practical_height'`.
#         The calculation method of :math:`\\mathrm{d}z`. Optional values are `'practical_height'`, `'scale_height'`.

#     Returns
#     -------
#     The TN wave activity 3D flux (:py:class:`xarray.DataArray<xarray.DataArray>`).
#     """

#     coordinate_sample_data = z_prime_data
#     u_c = u_climatology_data
#     v_c = v_climatology_data

#     lat_array = coordinate_sample_data['lat'].astype('float64')
#     coslat = np.cos(transfer_deg2rad(lat_array))

#     f = get_coriolis_parameter(lat_array, omega = omega)
#     psi_p = z_prime_data *g /f

#     p_lev = transfer_data_units(coordinate_sample_data[vertical_dim], vertical_dim_units, 'Pa')
#     p = p_lev /1e5

#     if method == 'scale_height':
#         dz = -scale_height *np.log(p)
#     elif method == 'practical_height':
#         dz = calc_gradient(z_data, dim = vertical_dim)

#     dpsi_dlambda = calc_gradient(psi_p, dim = lon_dim)
#     dpsi_dphi = calc_gradient(psi_p, dim = lat_dim)
#     dpsi_dz = calc_gradient(psi_p, dim = vertical_dim) /dz

#     d2psi_dlambda2 = calc_gradient(dpsi_dlambda, dim = lon_dim)
#     d2psi_dphi2 = calc_gradient(dpsi_dphi, dim = lat_dim)
#     d2psi_dlambdadphi = calc_gradient(dpsi_dphi, dim = lon_dim)
#     d2psi_dlambdadz = calc_gradient(dpsi_dlambda, dim = vertical_dim) /dz
#     d2psi_dphidz = calc_gradient(dpsi_dphi, dim = vertical_dim) /dz

#     term_xu = dpsi_dlambda**2 - (psi_p *d2psi_dlambda2)
#     term_xv = dpsi_dlambda *dpsi_dphi - (psi_p *d2psi_dlambdadphi)
#     term_yv = dpsi_dphi**2 - (psi_p *d2psi_dphi2)

#     term_zu = dpsi_dlambda *dpsi_dz - psi_p *d2psi_dlambdadz
#     term_zv = dpsi_dphi *dpsi_dz - psi_p *d2psi_dphidz

#     magU = np.sqrt(u_c**2 + v_c**2)
#     coeff = p /(2 *magU)

#     potential_temperature_data = calc_potential_temperature(temper_data, vertical_dim = vertical_dim, vertical_units = vertical_units, kappa = kappa)
#     N = calc_brunt_vaisala_frequency_atm(potential_temperature_data, z_data, vertical_dim = vertical_dim, g = g)

#     fx = coeff *( (u_c /(R**2 *coslat) *term_xu) + (v_c /(R**2) *term_xv) )
#     fy = coeff *( (u_c /(R**2) *term_xv) + (v_c *coslat /(R**2) *term_yv) )
#     fz = coeff * (f**2 /N**2)( u_c *term_zu + v_c *coslat /R *term_zv )

#     result = xr.Dataset()
#     result['psi_p'] = psi_p
#     result['fx'] = fx
#     result['fy'] = fy
#     result['fz'] = fz
#     return result


def calc_EP_horizontal_flux(
    u_prime_data: xr.DataArray,
    v_prime_data: xr.DataArray,
    time_dim: str = "time",
    lat_dim: str = "lat",
) -> xr.Dataset:
    """
    Calculate horizontal Eliassen–Palm Flux.

    Parameters
    ----------
    u_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The anormaly of zonal wind data.
    v_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The anormaly of meridional wind data.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

    Returns
    -------
    The Eliassen–Palm Flux (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. seealso::
        - https://www.ncl.ucar.edu/Applications/EPflux.shtml
        - https://renqlsysu.github.io/2020/02/16/wave_activity_flux/
    """
    lat_array = u_prime_data[lat_dim].astype("float64")
    coslat = np.cos(transfer_deg2rad(lat_array))

    fx = 0.5 * (
        (v_prime_data**2).mean(dim=time_dim) - (u_prime_data**2).mean(dim=time_dim)
    )
    fy = -(u_prime_data * v_prime_data).mean(dim=time_dim)

    result = xr.Dataset()
    result["fx"] = fx * coslat
    result["fy"] = fy * coslat
    return result
