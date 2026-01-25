"""
The calculation of transient eddy
"""

from __future__ import annotations
import xarray as xr
import numpy as np
from ..physics.geo import get_coriolis_parameter
from ..physics.temperature import calc_potential_temperature_vertical
from ..physics.convection.stability import calc_brunt_vaisala_frequency_atm
from .diff import (
    calc_gradient,
    calc_dx_gradient,
    calc_dy_gradient,
    calc_dlon_radian_gradient,
    calc_dlat_radian_gradient,
)
from .utility import transfer_deg2rad
from .units import (
    transfer_data_multiple_units,
    transfer_units_coeff,
)
from typing import Literal

__all__ = [
    "calc_eady_growth_rate",
    "calc_apparent_heat_source",
    "calc_total_diabatic_heating",
    "calc_apparent_moisture_sink",
    "calc_Plumb_wave_activity_horizontal_flux",
    "calc_TN_wave_activity_horizontal_flux",
    "calc_EP_horizontal_flux",
    "calc_monthly_rossby_wave_source",
]


def calc_eady_growth_rate(
    u_daily_data: xr.DataArray,
    z_daily_data: xr.DataArray,
    temper_daily_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
    lat_dim: str = "lat",
) -> xr.Dataset:
    """
    Calculate the maximum Eady growth rate.

    .. math::
        \\sigma = 0.3098 \\frac{f}{N} \\frac{\\mathrm{d} U}{\\mathrm{d} z}

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

    .. attention:: The unit of `z_daily_data` should be **meters**, NOT :math:`\\mathrm{m^2 \\cdot s^2}` which is the unit used in the representation of potential energy.

    temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Daily air temperature.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

    Returns
    -------
    The maximum Eady growth rate (:py:class:`xarray.Dataset<xarray.Dataset>`).

    - `eady_growth_rate`: The maximum Eady growth rate.
    - `dudz`: :math:`\\frac{\\mathrm{d} U}{\\mathrm{d} z}`
    - `brunt_vaisala_frequency`: Brunt-väisälä frequency.

    .. seealso::
        - Eady, E. T. (1949). Long Waves and Cyclone Waves. Tellus, 1(3), 33–52. https://doi.org/10.3402/tellusa.v1i3.8507, https://www.tandfonline.com/doi/abs/10.3402/tellusa.v1i3.8507
        - Lindzen, R. S. , & Farrell, B. (1980). A Simple Approximate Result for the Maximum Growth Rate of Baroclinic Instabilities. Journal of Atmospheric Sciences, 37(7), 1648-1654. https://journals.ametsoc.org/view/journals/atsc/37/7/1520-0469_1980_037_1648_asarft_2_0_co_2.xml
        - Simmonds, I., and E.-P. Lim (2009), Biases in the calculation of Southern Hemisphere mean baroclinic eddy growth rate, Geophys. Res. Lett., 36, L01707, https://doi.org/10.1029/2008GL036320.
        - Sloyan, B. M., and T. J. O'Kane (2015), Drivers of decadal variability in the Tasman Sea, J. Geophys. Res. Oceans, 120, 3193–3210, https://doi.org/10.1002/2014JC010550.
        - `eady_growth_rate -NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/eady_growth_rate.shtml>`__
        - `瞬变涡旋诊断量 <https://renqlsysu.github.io/2020/02/16/wave_activity_flux/>`__

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_egr.py
    """
    dim_tuple = u_daily_data.dims
    f = get_coriolis_parameter(u_daily_data[lat_dim])
    dp = 1
    dudp = calc_gradient(u_daily_data, dim=vertical_dim) / dp
    dzdp = calc_gradient(z_daily_data, dim=vertical_dim) / dp
    dudz = dudp / dzdp
    pt = calc_potential_temperature_vertical(
        temper_daily_data,
        vertical_dim=vertical_dim,
        vertical_dim_units=vertical_dim_units,
    )
    brunt_vaisala_atm = calc_brunt_vaisala_frequency_atm(
        pt, z_daily_data, vertical_dim=vertical_dim
    )
    eady_growth_rate = 0.3098 * np.abs(f) * np.abs(dudz) / brunt_vaisala_atm
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
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
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

    pt = calc_potential_temperature_vertical(
        temper_data, vertical_dim=vertical_dim, vertical_dim_units=vertical_dim_units
    )
    dtheta_dt = calc_gradient(pt, dim=time_dim) / dt
    dtheta_dx = calc_dx_gradient(pt, lon_dim=lon_dim, lat_dim=lat_dim)
    dtheta_dy = calc_dy_gradient(pt, lat_dim=lat_dim)

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
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
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
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
    time_units: str,
    specific_humidity_data_units: str,
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
    specific_humidity_data_units: :py:class:`str <str>`.
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
    specific_humidity_data = transfer_data_multiple_units(
        specific_humidity_data, specific_humidity_data_units, "kg/kg"
    )

    # Convert time units to seconds
    dt = transfer_units_coeff(time_units, "seconds")
    # Convert the pressure unit to Pascal
    dp_base = transfer_units_coeff(vertical_dim_units, "Pa")

    dqs_dt = calc_gradient(specific_humidity_data, dim=time_dim) / dt
    dqs_dx = calc_dx_gradient(specific_humidity_data, lon_dim=lon_dim, lat_dim=lat_dim)
    dqs_dy = calc_dy_gradient(specific_humidity_data, lat_dim=lat_dim)

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
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
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

    dim_tuple = z_prime_data.dims

    lat_array = coordinate_sample_data["lat"].astype("float64")
    coslat = np.cos(transfer_deg2rad(lat_array))

    f = get_coriolis_parameter(lat_array, omega=omega)
    psi_p = z_prime_data * g / f

    p_lev = transfer_data_multiple_units(
        coordinate_sample_data[vertical_dim], vertical_dim_units, "Pa"
    )
    p = p_lev / 1e5

    dpsi_dlambda = calc_dlon_radian_gradient(psi_p, lon_dim=lon_dim)
    dpsi_dphi = calc_dlat_radian_gradient(psi_p, lat_dim=lat_dim)

    d2psi_dlambda2 = calc_dlon_radian_gradient(dpsi_dlambda, lon_dim=lon_dim)
    d2psi_dlambdadphi = calc_dlon_radian_gradient(dpsi_dphi, lon_dim=lon_dim)

    term_xu = dpsi_dlambda**2 - (psi_p * d2psi_dlambda2)
    term_xv = dpsi_dlambda * dpsi_dphi - (psi_p * d2psi_dlambdadphi)

    fx = p * ((1 / (2 * R**2 * coslat)) * term_xu)
    fy = p * ((1 / (2 * R**2)) * term_xv)

    result = xr.Dataset()

    psi_p = psi_p.transpose(*dim_tuple)
    psi_p.attrs = {
        "long_name": "Perturbation stream function",
        "units": "m^2/s",
        "description": "Geostrophic stream function perturbation calculated from geopotential height anomaly",
    }
    result["psi_p"] = psi_p

    fx = fx.transpose(*dim_tuple)
    fx.attrs = {
        "long_name": "Zonal component of Plumb wave activity horizontal flux",
        "units": "m^2/s^2",
        "description": "Eastward component of Plumb wave activity flux",
        "standard_name": "tn_wave_activity_flux_x",
    }
    result["fx"] = fx

    fy = fy.transpose(*dim_tuple)
    fy.attrs = {
        "long_name": "Meridional component of Plumb wave activity horizontal flux",
        "units": "m^2/s^2",
        "description": "Northward component of Plumb wave activity flux",
        "standard_name": "Plumb_wave_activity_flux_y",
    }
    result["fy"] = fy

    result.attrs = {
        "title": "Plumb Wave Activity Flux",
        "description": "Horizontal components of quasi-geostrophic wave activity flux following Plumb (1985)",
        "reference": "Plumb, R. A., 1985: On the Three-Dimensional Propagation of Stationary Waves. J. Atmos. Sci., 42, 217–229",
        "created_with": "easyclimate: calc_Plumb_wave_activity_horizontal_flux function",
    }
    return result.astype("float32")


def calc_TN_wave_activity_horizontal_flux(
    z_prime_data: xr.DataArray | None,
    u_climatology_data: xr.DataArray,
    v_climatology_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
    psi_prime_data: xr.DataArray | None = None,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    time_dim: str = "time",
    omega: float = 7.292e-5,
    g: float = 9.8,
    R: float = 6371200.0,
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

    .. attention:: The unit of `z_prime_data` should be **meters**, NOT :math:`\\mathrm{m^2 \\cdot s^2}` which is the unit used in the representation of potential energy.

    u_climatology_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The climatology of zonal wind data.
    v_climatology_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The climatology of meridional wind data.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    psi_prime_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Perturbation stream function. Geostrophic stream function perturbation calculated from geopotential height anomaly.
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

    Reference
    --------------
    - Takaya, K., & Nakamura, H. (2001). A Formulation of a Phase-Independent Wave-Activity Flux for Stationary and Migratory Quasigeostrophic Eddies on a Zonally Varying Basic Flow. Journal of the Atmospheric Sciences, 58(6), 608-627. https://journals.ametsoc.org/view/journals/atsc/58/6/1520-0469_2001_058_0608_afoapi_2.0.co_2.xml

    .. seealso::
        - http://www.atmos.rcast.u-tokyo.ac.jp/nishii/programs/index.html
        - https://github.com/laishenggx/T-N_Wave-Activity-Flux
    """
    # Check data validity
    if z_prime_data is None and psi_prime_data is None:
        raise ValueError(
            "`z_prime_data` and `psi_prime_data` should not both be None at the same time."
        )

    # Select the data to be used (`z_prime_data` is preferred)
    vertical_level_data = z_prime_data if z_prime_data is not None else psi_prime_data

    # Verify data type
    if not isinstance(vertical_level_data, xr.DataArray):
        raise TypeError(
            f"The input data must be of the xarray.DataArray type. Current type:{type(vertical_level_data)}"
        )

    # Extract month from `vertical_level_data` and select corresponding climatology data
    if time_dim in vertical_level_data.dims and time_dim in u_climatology_data.dims:
        # If `vertical_level_data` has multiple time steps, extract month for each
        months = vertical_level_data[time_dim].dt.month

        # Select corresponding climatology months using groupby
        u_c = u_climatology_data.sel(
            {time_dim: u_climatology_data[time_dim].dt.month.isin(months)}
        )
        v_c = v_climatology_data.sel(
            {time_dim: v_climatology_data[time_dim].dt.month.isin(months)}
        )

        # For multiple time steps, we need to align the data
        # This assumes u_climatology_data and v_climatology_data have the same time coordinates
        if len(vertical_level_data[time_dim]) > 1:
            # Create a mapping between month and climatology data
            u_c_monthly = u_climatology_data.groupby(
                u_climatology_data[time_dim].dt.month
            )
            v_c_monthly = v_climatology_data.groupby(
                v_climatology_data[time_dim].dt.month
            )

            # Select appropriate months for each time step in `vertical_level_data`
            u_c_list = []
            v_c_list = []

            for t in vertical_level_data[time_dim]:
                month = t.dt.month.values
                month = int(month)
                u_c_month = u_c_monthly[month].squeeze(drop=True)
                v_c_month = v_c_monthly[month].squeeze(drop=True)
                u_c_list.append(u_c_month)
                v_c_list.append(v_c_month)

            u_c = xr.concat(u_c_list, dim=time_dim)
            v_c = xr.concat(v_c_list, dim=time_dim)
            u_c[time_dim] = vertical_level_data[time_dim]
            v_c[time_dim] = vertical_level_data[time_dim]
    elif time_dim in u_climatology_data.dims:
        # If `vertical_level_data` has no time dimension, use the first month as default
        # or raise an error depending on your use case
        u_c = u_climatology_data.isel({time_dim: 0})
        v_c = v_climatology_data.isel({time_dim: 0})
    else:
        u_c = u_climatology_data
        v_c = v_climatology_data

    coordinate_sample_data = vertical_level_data
    dim_tuple = vertical_level_data.dims

    lat_array = coordinate_sample_data["lat"].astype("float64")
    coslat = np.cos(transfer_deg2rad(lat_array))

    if psi_prime_data is None:
        f = get_coriolis_parameter(lat_array, omega=omega)
        psi_p = z_prime_data * g / f
    else:
        psi_p = psi_prime_data

    p_lev = transfer_data_multiple_units(
        coordinate_sample_data[vertical_dim], vertical_dim_units, "Pa"
    )
    p = p_lev / 1e5

    dpsi_dlambda = calc_dlon_radian_gradient(psi_p, lon_dim=lon_dim)
    dpsi_dphi = calc_dlat_radian_gradient(psi_p, lat_dim=lat_dim)

    d2psi_dlambda2 = calc_dlon_radian_gradient(dpsi_dlambda, lon_dim=lon_dim)
    d2psi_dphi2 = calc_dlat_radian_gradient(dpsi_dphi, lat_dim=lat_dim)
    d2psi_dlambdadphi = calc_dlon_radian_gradient(dpsi_dphi, lon_dim=lon_dim)

    term_xu = dpsi_dlambda**2 - (psi_p * d2psi_dlambda2)
    term_xv = dpsi_dlambda * dpsi_dphi - (psi_p * d2psi_dlambdadphi)
    term_yv = dpsi_dphi**2 - (psi_p * d2psi_dphi2)

    magU = np.sqrt(u_c**2 + v_c**2)
    coeff = p / (2 * magU)

    fx = coeff * ((u_c / (R**2 * coslat) * term_xu) + (v_c / (R**2) * term_xv))
    fy = coeff * ((u_c / (R**2) * term_xv) + (v_c * coslat / (R**2) * term_yv))

    result = xr.Dataset()

    psi_p = psi_p.transpose(*dim_tuple)
    psi_p.attrs = {
        "long_name": "Perturbation stream function",
        "units": "m^2/s",
        "description": "Geostrophic stream function perturbation calculated from geopotential height anomaly",
    }
    result["psi_p"] = psi_p

    fx = fx.transpose(*dim_tuple)
    fx.attrs = {
        "long_name": "Zonal component of TN wave activity horizontal flux",
        "units": "m^2/s^2",
        "description": "Eastward component of Takaya and Nakamura wave activity flux",
        "standard_name": "tn_wave_activity_flux_x",
    }
    result["fx"] = fx

    fy = fy.transpose(*dim_tuple)
    fy.attrs = {
        "long_name": "Meridional component of TN wave activity horizontal flux",
        "units": "m^2/s^2",
        "description": "Northward component of Takaya and Nakamura wave activity flux",
        "standard_name": "tn_wave_activity_flux_y",
    }
    result["fy"] = fy

    result.attrs = {
        "title": "Takaya and Nakamura Wave Activity Flux",
        "description": "Horizontal components of quasi-geostrophic wave activity flux following Takaya and Nakamura (2001)",
        "reference": "Takaya, K. and H. Nakamura, 2001: A Formulation of a Phase-Independent Wave-Activity Flux for Stationary and Migratory Quasigeostrophic Eddies on a Zonally Varying Basic Flow. J. Atmos. Sci., 58, 608-627.",
        "created_with": "easyclimate: calc_TN_wave_activity_horizontal_flux function",
    }
    return result.astype("float32")


def calc_TN_wave_activity_vertical_flux(
    z_prime_data: xr.DataArray,
    u_climatology_data: xr.DataArray,
    v_climatology_data: xr.DataArray,
    temper_climatology_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    omega: float = 7.292e-5,
    g: float = 9.8,
    R: float = 6371200.0,
    scale_height: float = 8000.0,
    kappa: float = 287 / 1005.7,
    p_0=1000.0,
) -> xr.Dataset:
    """ """

    coordinate_sample_data = z_prime_data
    u_c = u_climatology_data
    v_c = v_climatology_data
    t_c = temper_climatology_data

    lat_array = coordinate_sample_data["lat"].astype("float64")
    coslat = np.cos(transfer_deg2rad(lat_array))

    f = get_coriolis_parameter(lat_array, omega=omega)
    psi_p = z_prime_data * g / f

    p_lev = transfer_data_multiple_units(
        coordinate_sample_data[vertical_dim], vertical_dim_units, "hPa"
    )  # hPa
    scale_z = -scale_height * np.log(p_lev / p_0)
    dz = calc_gradient(scale_z, dim=vertical_dim)  # m, hPa

    dpsi_dlambda = calc_dlon_radian_gradient(psi_p, lon_dim=lon_dim)
    dpsi_dphi = calc_dlat_radian_gradient(psi_p, lat_dim=lat_dim)
    dpsi_dz = calc_gradient(psi_p, dim=vertical_dim) / dz

    d2psi_dlambdadz = calc_gradient(dpsi_dlambda, dim=vertical_dim) / dz
    d2psi_dphidz = calc_gradient(dpsi_dphi, dim=vertical_dim) / dz

    term_zu = dpsi_dlambda * dpsi_dz - psi_p * d2psi_dlambdadz
    term_zv = dpsi_dphi * dpsi_dz - psi_p * d2psi_dphidz

    magU = np.sqrt(u_c**2 + v_c**2)
    coeff = p_lev / (2 * magU)

    potential_temperature_data = calc_potential_temperature_vertical(
        t_c, vertical_dim=vertical_dim, vertical_units=vertical_dim_units, kappa=kappa
    )
    N = calc_brunt_vaisala_frequency_atm(
        potential_temperature_data, scale_z, vertical_dim=vertical_dim, g=g
    )

    fz = coeff * (f**2 / N**2)(u_c * term_zu + v_c * coslat / R * term_zv)

    result = xr.Dataset()
    # result['psi_p'] = psi_p
    result["fz"] = fz
    return result


# def calc_TN_wave_activity_3D_flux(
#     z_prime_data: xr.DataArray,
#     u_climatology_data: xr.DataArray,
#     v_climatology_data: xr.DataArray,
#     temper_data: xr.DataArray,
#     vertical_dim: str,
#     vertical_dim_units: Literal["hPa", "Pa", "mbar"],
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

#     p_lev = transfer_data_multiple_units(coordinate_sample_data[vertical_dim], vertical_dim_units, 'Pa')
#     p = p_lev /1e5

#     if method == 'scale_height':
#         dz = -scale_height *np.log(p)
#     elif method == 'practical_height':
#         dz = calc_gradient(z_data, dim = vertical_dim)

#     dpsi_dlambda = calc_dlon_radian_gradient(psi_p, lon_dim = lon_dim)
#     dpsi_dphi = calc_dlat_radian_gradient(psi_p, lat_dim = lat_dim)
#     dpsi_dz = calc_gradient(psi_p, dim = vertical_dim) /dz

#     d2psi_dlambda2 = calc_dlon_radian_gradient(dpsi_dlambda, lon_dim = lon_dim)
#     d2psi_dphi2 = calc_dlat_radian_gradient(dpsi_dphi, lat_dim = lat_dim)
#     d2psi_dlambdadphi = calc_dlon_radian_gradient(dpsi_dphi, lon_dim = lon_dim)
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


def calc_monthly_rossby_wave_source(
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    u_climatology_data: xr.DataArray,
    v_climatology_data: xr.DataArray,
    lat_dim: str = "lat",
    omega: float = 7.292e-05,
    R: float = 6371200.0,
) -> xr.DataArray:
    """
    Calculate the Rossby wave source following Sardeshmukh and Hoskins (1988).

    The Rossby wave source (RWS) represents the forcing term for stationary Rossby waves
    and is widely used to diagnose atmospheric teleconnections. It is decomposed into
    five terms representing different physical mechanisms:

    .. math::

        S' = -\\nabla \\cdot (\\mathbf{v}_\\chi \\zeta)' = \\text{term1a} + \\text{term1b} + \\text{term2a} + \\text{term2b} + \\text{term3} + \\text{term4} + \\text{term5}

    where:

    - **term1a**: :math:`-\\bar{\\zeta} \\nabla \\cdot \\mathbf{v}'_\\chi`
      (stretching of climatological vorticity by divergent wind anomalies)
    - **Term1b**: :math:`-f \\nabla \\cdot \\mathbf{v'_\\chi}`
      (Stretching of planetary vorticity by anomalous divergence)
    - **term2a**: :math:`-\\mathbf{v}'_\\chi \\cdot \\nabla\\bar{\\zeta}`
      (advection of climatological vorticity by divergent wind anomalies)
    - **term2b**: :math:`-\\beta v'_\\chi`
      (planetary vorticity advection, where :math:`\\beta = \\partial f/\\partial y`)
    - **term3**: :math:`-\\zeta' \\nabla \\cdot \\bar{\\mathbf{v}}_\\chi`
      (stretching of vorticity anomalies by climatological divergent wind)
    - **term4**: :math:`-\\bar{\\mathbf{v}}_\\chi \\cdot \\nabla\\zeta'`
      (advection of vorticity anomalies by climatological divergent wind)
    - **term5**: :math:`-\\nabla \\cdot (\\mathbf{v}'_\\chi \\zeta')`
      (nonlinear term: divergence of transient vorticity flux)

    where :math:`\\mathbf{v}_\\chi = (u_\\chi, v_\\chi)` is the divergent (irrotational)
    wind component, :math:`\\zeta` is relative vorticity, overbar denotes climatology,
    and prime denotes anomaly.

    .. note::
        - Anomalies are computed as deviations from the monthly climatology.
        - The planetary term (term6) represents the beta effect, which is important for Rossby wave propagation, especially in the tropics and subtropics.
        - Terms 1, 2, and 6 typically dominate the Rossby wave source.
        - Positive RWS indicates a cyclonic wave source; negative indicates anticyclonic.

    Parameters
    ----------
    u_data : :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{m/s}`)
        Zonal wind component. Must contain a 'time' dimension with monthly or
        higher frequency data. Expected dimensions: ``(time, lat, lon)`` or ``(time, level, lat, lon)``.
    v_data : :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{m/s}`)
        Meridional wind component. Must have the same dimensions as `u_data`.
    u_climatology_data : :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{m/s}`)
        Monthly climatology of zonal wind. Expected dimensions: ``(time, lat, lon)``
        or ``(time, level, lat, lon)``, where time/month ranges from 1 to 12.
    v_climatology_data : :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{m/s}`)
        Monthly climatology of meridional wind. Must have the same dimensions
        as `u_climatology_data`.

    Returns
    -------
    result : :py:class:`xarray.Dataset<xarray.Dataset>`
        Dataset containing the Rossby wave source and its five decomposed terms:

        - **RWS**: Total Rossby wave source (sum of all terms)
        - **term1a**: Vorticity stretching by anomalous divergence
        - **term1b**: Stretching of planetary vorticity by anomalous divergence
        - **term2a**: Vorticity advection by anomalous divergent wind
        - **term2b**: Planetary vorticity advection (beta effect)
        - **term3**: Anomalous vorticity stretching by climatological divergence
        - **term4**: Anomalous vorticity advection by climatological divergent wind
        - **term5**: Nonlinear transient eddy term


        All variables have units of :math:`\\mathrm{s^{-2}}` and retain the spatial/temporal dimensions
        of the input data.

    .. seealso::

        - Sardeshmukh, P. D., & Hoskins, B. J. (1988). The generation of global rotational flow by steady idealized tropical divergence. Journal of the Atmospheric Sciences, 45(7), 1228-1251. https://journals.ametsoc.org/view/journals/atsc/45/7/1520-0469_1988_045_1228_tgogrf_2_0_co_2.xml, https://doi.org/10.1175/1520-0469(1988)045<1228:TGOGRF>2.0.CO;2

    Examples
    --------
    >>> # Calculate monthly climatology
    >>> u_clim = u.groupby('time.month').mean('time')
    >>> v_clim = v.groupby('time.month').mean('time')
    >>>
    >>> # Compute Rossby wave source
    >>> rws_result = calc_rossby_wave_source(u, v, u_clim, v_clim)
    >>>
    >>> # Visualize the dominant terms
    >>> rws_result['RWS'].sel(time='2015-12').plot()
    >>> (rws_result['term1'] + rws_result['term2']).sel(time='2015-12').plot()
    """
    from .windspharm import (
        calc_relative_vorticity,
        calc_irrotational_component,
        calc_gradient,
    )
    from ..physics.geo import get_coriolis_parameter
    from .utility import transfer_deg2rad

    # planetary vorticity
    lat = v_data[lat_dim]

    # Real Value
    zeta = calc_relative_vorticity(u_data, v_data)
    ir_result = calc_irrotational_component(u_data, v_data)
    uchi = ir_result.uchi
    vchi = ir_result.vchi

    # Climate
    zeta_climate = calc_relative_vorticity(u_climatology_data, v_climatology_data)
    ir_climate_result = calc_irrotational_component(
        u_climatology_data, v_climatology_data
    )
    uchi_climate = ir_climate_result.uchi
    vchi_climate = ir_climate_result.vchi

    # Prime
    zeta_prime = zeta.groupby("time.month") - zeta_climate.groupby("time.month").mean()
    uchi_prime = uchi.groupby("time.month") - uchi_climate.groupby("time.month").mean()
    vchi_prime = vchi.groupby("time.month") - vchi_climate.groupby("time.month").mean()

    # Term1a: \bar{\zeta} \nabla \cdot \bm{v'_\xi}
    uchi_prime_grdx = calc_gradient(uchi_prime)["zonal_gradient"]
    vchi_prime_grdy = calc_gradient(vchi_prime)["meridional_gradient"]
    term1a = -(
        zeta_climate.groupby("time.month").mean()
        * ((uchi_prime_grdx + vchi_prime_grdy).groupby("time.month"))
    )

    # Term1b: -f \nabla \cdot \bm{v'_\chi}
    f = get_coriolis_parameter(lat, omega=omega)
    term1b = -f * (uchi_prime_grdx + vchi_prime_grdy)

    # Term2a: \bm{v'_\xi} \cdot \nabla\bar{\zeta}
    zeta_climate_result = calc_gradient(zeta_climate)
    zeta_climate_grdx = zeta_climate_result["zonal_gradient"]
    zeta_climate_grdy = zeta_climate_result["meridional_gradient"]
    term2a = -(
        uchi_prime.groupby("time.month")
        * zeta_climate_grdx.groupby("time.month").mean()
        + vchi_prime.groupby("time.month")
        * zeta_climate_grdy.groupby("time.month").mean()
    )

    # Term2b: -\beta v'_\chi
    # Calculate beta = df/dy
    # beta = (2 * Omega * cos(phi)) / a
    lat_rad = transfer_deg2rad(lat)
    beta = (2 * omega * np.cos(lat_rad)) / R
    # beta * v'_chi
    term2b = -(beta * vchi_prime)

    # Term3: \zeta' \nabla \cdot \bar{v'_\xi}
    uchi_climate_grdx = calc_gradient(uchi_climate)["zonal_gradient"]
    vchi_climate_grdy = calc_gradient(vchi_climate)["meridional_gradient"]
    term3 = -(
        zeta_prime.groupby("time.month")
        * ((uchi_climate_grdx + vchi_climate_grdy).groupby("time.month").mean())
    )

    # Term4: \bm{\bar{v}_\xi} \cdot \nabla \zeta'
    zeta_prime_result = calc_gradient(zeta_prime)
    zeta_prime_grdx = zeta_prime_result["zonal_gradient"]
    zeta_prime_grdy = zeta_prime_result["meridional_gradient"]
    term4 = -(
        uchi_climate.groupby("time.month").mean()
        * zeta_prime_grdx.groupby("time.month")
        + vchi_climate.groupby("time.month").mean()
        * zeta_prime_grdy.groupby("time.month")
    )

    # Term5: \nabla \cdot (\bm{v}'_\chi \zeta')
    # u'_\chi * \zeta'
    uchi_prime_zeta_prime = uchi_prime * zeta_prime
    # v'_\chi * \zeta'
    vchi_prime_zeta_prime = vchi_prime * zeta_prime
    uchi_prime_zeta_prime_grdx = calc_gradient(uchi_prime_zeta_prime)["zonal_gradient"]
    vchi_prime_zeta_prime_grdy = calc_gradient(vchi_prime_zeta_prime)[
        "meridional_gradient"
    ]
    term5 = -(uchi_prime_zeta_prime_grdx + vchi_prime_zeta_prime_grdy)

    # Compute total RWS
    rws = term1a + term1b + term2a + term2b + term3 + term4 + term5

    # Obtain the dimension order
    dim_array = u_data.dims

    # Create Data Variables
    result = xr.Dataset()

    # Add the data variables for each item
    result["term1a"] = term1a.transpose(*dim_array)
    result["term1b"] = term1b.transpose(*dim_array)
    result["term2a"] = term2a.transpose(*dim_array)
    result["term2b"] = term2b.transpose(*dim_array)
    result["term3"] = term3.transpose(*dim_array)
    result["term4"] = term4.transpose(*dim_array)
    result["term5"] = term5.transpose(*dim_array)
    result["RWS"] = rws.transpose(*dim_array)

    # Set attributes for each variable
    result["term1a"].attrs = {
        "long_name": "RWS Term 1a: Vorticity stretching",
        "description": "Stretching of climatological vorticity by anomalous divergent wind",
        "formula": "-zeta_bar * div(v_chi_prime)",
        "units": "s-2",
        "standard_name": "rossby_wave_source_term1",
    }

    result["term1b"].attrs = {
        "long_name": "RWS Term 1b: Planetary vorticity stretching",
        "formula": "-f * div(v_chi_prime)",
        "units": "s-2",
        "description": "Stretching of planetary vorticity by anomalous divergence",
    }

    result["term2a"].attrs = {
        "long_name": "RWS Term 2a: Vorticity advection",
        "description": "Advection of climatological vorticity gradient by anomalous divergent wind",
        "formula": "-v_chi_prime · grad(zeta_bar)",
        "units": "s-2",
        "standard_name": "rossby_wave_source_term2",
    }

    result["term2b"].attrs = {
        "long_name": "RWS Term 2b: Planetary vorticity advection",
        "description": "Planetary vorticity advection (beta effect): meridional advection of planetary vorticity by anomalous divergent wind",
        "formula": "-beta * v_chi_prime, where beta = df/dy",
        "units": "s-2",
        "standard_name": "rossby_wave_source_planetary",
        "note": "Important for Rossby wave propagation, especially in tropics and subtropics",
    }

    result["term3"].attrs = {
        "long_name": "RWS Term 3: Anomalous vorticity stretching",
        "description": "Stretching of vorticity anomalies by climatological divergent wind",
        "formula": "-zeta_prime * div(v_chi_bar)",
        "units": "s-2",
        "standard_name": "rossby_wave_source_term3",
    }

    result["term4"].attrs = {
        "long_name": "RWS Term 4: Anomalous vorticity advection",
        "description": "Advection of vorticity anomalies by climatological divergent wind",
        "formula": "-v_chi_bar · grad(zeta_prime)",
        "units": "s-2",
        "standard_name": "rossby_wave_source_term4",
    }

    result["term5"].attrs = {
        "long_name": "RWS Term 5: Nonlinear transient term",
        "description": "Divergence of transient vorticity flux (nonlinear eddy term)",
        "formula": "-div(v_chi_prime * zeta_prime)",
        "units": "s-2",
        "standard_name": "rossby_wave_source_term5",
    }

    result["RWS"].attrs = {
        "long_name": "Total Rossby Wave Source",
        "description": "Total Rossby wave source including both relative and planetary vorticity effects",
        "formula": "-div(v_chi * eta)' where eta = zeta + f",
        "units": "s-2",
        "standard_name": "rossby_wave_source_total",
        "references": "Sardeshmukh and Hoskins (1988, JAS)",
        "positive": "cyclonic wave source",
        "negative": "anticyclonic wave source",
    }

    # Add global attributes
    result.attrs.update(
        {
            "title": "Rossby Wave Source Decomposition with Planetary Term",
            "method": "Six-term decomposition (5 relative + 1 planetary) following Sardeshmukh and Hoskins (1988)",
            "anomaly_method": "Deviation from monthly climatology",
            "convention": "Negative divergence convention (sources are positive for cyclonic forcing)",
            "references": "Sardeshmukh, P. D., & Hoskins, B. J. (1988). JAS, 45(7), 1228-1251.",
            "note": "Terms 1, 2, and 6 (planetary) typically dominate the Rossby wave source",
            "earth_radius": f"{R} m",
            "earth_rotation_rate": f"{omega} rad/s",
            "created_with": "easyclimate: calc_rossby_wave_source function",
        }
    )

    return result
