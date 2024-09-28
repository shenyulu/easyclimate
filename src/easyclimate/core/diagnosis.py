"""
Functions for Weather and climate variable diagnosis.
"""

from __future__ import annotations
import xarray as xr
import numpy as np
from .diff import calc_gradient, calc_p_gradient
from .utility import (
    transfer_deg2rad,
    transfer_units_coeff,
    transfer_data_multiple_units,
    transfer_data_difference_units,
)
from typing import Literal
import warnings

__all__ = [
    "calc_brunt_vaisala_frequency_atm",
    "get_coriolis_parameter",
    "calc_potential_temperature",
    "calc_virtual_temperature_Hobbs2006",
    "calc_virtual_temperature",
    "calc_static_stability",
    "calc_dewpoint",
    "calc_mixing_ratio",
    "calc_vapor_pressure",
    "calc_saturation_vapor_pressure",
    "calc_saturation_mixing_ratio",
    "transfer_mixing_ratio_2_specific_humidity",
    "transfer_specific_humidity_2_mixing_ratio",
    "transfer_dewpoint_2_specific_humidity",
    "transfer_specific_humidity_2_dewpoint",
    "transfer_dewpoint_2_relative_humidity",
    "transfer_mixing_ratio_2_relative_humidity",
    "transfer_specific_humidity_2_relative_humidity",
]


def calc_brunt_vaisala_frequency_atm(
    potential_temperature_data: xr.DataArray,
    z_data: xr.DataArray,
    vertical_dim: str,
    g: float = 9.8,
) -> xr.DataArray:
    """
    Calculation of the Brunt-väisälä frequency for the vertical atmosphere.

    .. math::
        N = \\left( \\frac{g}{\\theta} \\frac{\\mathrm{d}\\theta}{\\mathrm{d}z} \\right)^\\frac{1}{2}

    Parameters
    ----------
    potential_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Vertical atmospheric potential temperature.
    z_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Vertical atmospheric geopotential height.

    .. attention:: The unit of `z_data` should be **meters**, NOT :math:`\\mathrm{m^2 \\cdot s^2}` which is the unit used in the representation of potential energy.

    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    g: :py:class:`float <float>`, default: `9.8`.
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

    """
    dp = 1.0
    dtheta_dp = calc_gradient(potential_temperature_data, dim=vertical_dim) / dp
    dz_dp = calc_gradient(z_data, dim=vertical_dim) / dp
    N = np.sqrt((g / potential_temperature_data) * (dtheta_dp / dz_dp))

    try:
        z_data_units = z_data.attrs["units"]
        N.attrs["units"] = f"{z_data_units}^(1/2)"
    except:
        N.attrs["units"] = f"[L]^(1/2)"
        warnings.warn(
            "The variable of `z_data` do not have `units` attribution, so the attribution of units in it is assigned to the symbol for dimension!"
        )
    N.name = "brunt_vaisala_frequency_atm"
    return N


def get_coriolis_parameter(
    lat_data: xr.DataArray | np.array, omega: float = 7.292e-5
) -> xr.DataArray | np.array:
    """
    Calculate the Coriolis parameter at each point.

    .. math::
        f = 2 \\Omega \\sin(\\phi)

    Parameters
    ----------
    lat_data: :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`numpy.array <numpy.array>`.
        Latitude at each point.
    omega: :py:class:`float <float>`, default: `7.292e-5`.
        The angular speed of the earth.

    Returns
    -------
    Corresponding Coriolis force at each point (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`numpy.array <numpy.array>`).

    Reference
    --------------
    - `Coriolis parameter - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Coriolis_parameter>`__

    .. seealso::
        - `coriolis_parameter — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.coriolis_parameter.html>`__
        - `coriolis_param - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/coriolis_param.shtml>`__
    """
    lat_data = lat_data.astype("float64")
    return_data = 2 * omega * np.sin(transfer_deg2rad(lat_data))

    if isinstance(return_data, xr.DataArray):
        return_data.attrs["units"] = "s^-1"
        return_data.name = "coriolis_parameter"
    elif isinstance(return_data, np.ndarray):
        warnings.warn(
            "The unit for the output in `get_coriolis_parameter` is " + "s^-1"
        )
    else:
        raise TypeError("`lat_data` shuold be `xr.DataArray` or `np.array`.")
    return return_data


def calc_potential_temperature(
    temper_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: str,
    kappa: float = 287 / 1005.7,
) -> xr.DataArray:
    """
    Calculate the potential temperature.

    Uses the Poisson equation to calculation the potential temperature given pressure and temperature.

    .. math::
        \\theta = T \\left( \\frac{p_0}{p} \\right) ^\\kappa

    Parameters
    ----------
    temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Air temperature.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    kappa: :py:class:`float <float>`, default: `287/1005.7`.
        Poisson constant :math:`\\kappa`.

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
    """
    p_base = transfer_units_coeff(vertical_dim_units, "Pa")
    if p_base == 1.0:
        P_0 = 1000e2
    elif p_base == 100.0:
        P_0 = 1000
    else:
        raise ValueError("`vertical_dim_units` be `Pa`, `hPa`, `mbar`.")

    return_data = temper_data * (P_0 / temper_data[vertical_dim]) ** (kappa)

    try:
        return_data.attrs["units"] = temper_data.attrs["units"]
    except:
        return_data.attrs["units"] = "[T]"
        warnings.warn(
            "The variable of `temper_data` do not have `units` attribution, so the attribution of units in it is assigned to the symbol for dimension!"
        )
    return_data.name = "potential_temperature"
    return return_data


def calc_virtual_temperature(
    temper_data: xr.DataArray,
    specific_humidity_data: xr.DataArray,
    specific_humidity_units: Literal["kg/kg", "g/g", "g/kg"],
    epsilon: float = 0.608,
) -> xr.DataArray:
    """
    Calculate virtual temperature.

    The virtual temperature (:math:`T_v`) is the temperature at which dry air would have the same density as the moist air, at a given pressure.
    In other words, two air samples with the same virtual temperature have the same density, regardless of their actual temperature or relative humidity.
    The virtual temperature is always greater than  the absolute air temperature.

    .. math::
        T_v = T(1+ \\epsilon q)

    where :math:`\\epsilon = 0.608` when the mixing ratio (specific humidity) :math:`q` is expressed in :math:`\\mathrm{g \\cdot g^{-1}}`.

    Parameters
    ----------
    temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Air temperature.
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    specific_humidity_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/g`, `g/kg` and so on.
    epsilon: :py:class:`float <float>`.
        A constant.

    Reference
    --------------
    - Doswell , C. A., and E. N. Rasmussen, 1994: The Effect of Neglecting the Virtual Temperature Correction on CAPE Calculations. Wea. Forecasting, 9, 625–629, https://doi.org/10.1175/1520-0434(1994)009<0625:TEONTV>2.0.CO;2.
    - https://en.wikipedia.org/wiki/Virtual_temperature
    - https://glossary.ametsoc.org/wiki/Virtual_temperature

    .. seealso::
        - `virtual_temperature — MetPy <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.virtual_temperature.html>`__
        - `temp_virtual - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/temp_virtual.shtml>`__
    """
    specific_humidity_data = transfer_data_multiple_units(
        specific_humidity_data, specific_humidity_units, "g/g"
    )

    T_v = (1 + epsilon * specific_humidity_data) * temper_data

    try:
        T_v.attrs["units"] = temper_data.attrs["units"]
    except:
        T_v.attrs["units"] = "[T]"
        warnings.warn(
            "The variable of `temper_data` do not have `units` attribution, so the attribution of units in it is assigned to the symbol for dimension!"
        )
    T_v.name = "virtual_temperature"
    return T_v


def calc_virtual_temperature_Hobbs2006(
    temper_data: xr.DataArray,
    specific_humidity_data: xr.DataArray,
    specific_humidity_units: Literal["kg/kg", "g/g", "g/kg"],
    epsilon: float = 0.6219569100577033,
) -> xr.DataArray:
    """
    Calculate virtual temperature.

    The virtual temperature (:math:`T_v`) is the temperature at which dry air would have the same density as the moist air, at a given pressure.
    In other words, two air samples with the same virtual temperature have the same density, regardless of their actual temperature or relative humidity.
    The virtual temperature is always greater than  the absolute air temperature.

    This calculation must be given an air parcel's temperature and mixing ratio. The implementation uses the formula outlined in [Hobbs2006] pg.67 & 80.

    .. math::
        T_v = T \\frac{\\text{q} + \\epsilon}{\\epsilon\\,(1 + \\text{q})}

    where :math:`\\epsilon \\approx 0.622` when the mixing ratio (specific humidity) :math:`q` is expressed in :math:`\\mathrm{g \ g^{-1}}`.

    Parameters
    ----------
    temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Air temperature.
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    specific_humidity_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/g`, `g/kg` and so on.
    epsilon: :py:class:`float <float>`.
        The molecular weight ratio, which is molecular weight of the constituent gas to that assumed for air. Defaults to the ratio for water vapor to dry air. (:math:`\\epsilon \\approx 0.622`)

    Reference
    --------------
    - Hobbs, P. V., and J. M. Wallace, 2006: Atmospheric Science: An Introductory Survey. 2nd ed. Academic Press, 504 pp. https://www.sciencedirect.com/book/9780127329512/atmospheric-science
    - Doswell , C. A., and E. N. Rasmussen, 1994: The Effect of Neglecting the Virtual Temperature Correction on CAPE Calculations. Wea. Forecasting, 9, 625–629, https://doi.org/10.1175/1520-0434(1994)009<0625:TEONTV>2.0.CO;2.
    - https://en.wikipedia.org/wiki/Virtual_temperature
    - https://glossary.ametsoc.org/wiki/Virtual_temperature

    .. seealso::
        - `virtual_temperature — MetPy <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.virtual_temperature.html>`__
        - `temp_virtual - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/temp_virtual.shtml>`__
    """
    specific_humidity_data = transfer_data_multiple_units(
        specific_humidity_data, specific_humidity_units, "g/g"
    )

    T_v = temper_data * (
        (specific_humidity_data + epsilon) / (epsilon * (1 + specific_humidity_data))
    )

    try:
        T_v.attrs["units"] = temper_data.attrs["units"]
    except:
        T_v.attrs["units"] = "[T]"
        warnings.warn(
            "The variable of `temper_data` do not have `units` attribution, so the attribution of units in it is assigned to the symbol for dimension!"
        )
    T_v.name = "virtual_temperature"
    return T_v


def calc_static_stability(
    temper_data: xr.DataArray, vertical_dim: str, vertical_dim_units: str
) -> xr.DataArray:
    """
    Calculate the static stability within a vertical profile.

    .. math::
        \\sigma = - T \\frac{\\partial \\ln \\theta}{\\partial p}

    Parameters
    ----------
    temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Air temperature.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
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
    """
    theta = calc_potential_temperature(
        temper_data, vertical_dim=vertical_dim, vertical_dim_units=vertical_dim_units
    )
    ln_theta = np.log(theta)
    part = calc_p_gradient(
        ln_theta, vertical_dim=vertical_dim, vertical_dim_units=vertical_dim_units
    )
    return_data = -temper_data * part

    try:
        temper_data_units = temper_data.attrs["units"]
        return_data.attrs["units"] = f"{temper_data_units}^2 {vertical_dim_units}^-1"
    except:
        return_data.attrs["units"] = f"[T]^2 {vertical_dim_units}^-1"
        warnings.warn(
            "The variable of `temper_data` do not have `units` attribution, so the attribution of units in it is assigned to the symbol for dimension!"
        )
    return_data.name = "static_stability"
    return return_data


def calc_dewpoint(
    vapor_pressure_data: xr.DataArray,
    vapor_pressure_data_units: Literal["hPa", "Pa"],
) -> xr.DataArray:
    """
    Calculate the ambient dewpoint given the vapor pressure.

    Parameters
    ----------
    vapor_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Water vapor partial pressure.
    total_pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `total_pressure_data` value. Optional values are `hPa`, `Pa`.

    Returns
    -------
    The dew point (:py:class:`xarray.DataArray<xarray.DataArray>`), degrees Celsius.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.dewpoint.html
    """
    vapor_pressure_data = transfer_data_multiple_units(
        vapor_pressure_data, vapor_pressure_data_units, "hPa"
    )

    val = np.log(vapor_pressure_data / 6.112)
    dewpoint = 243.5 * val / (17.67 - val)
    dewpoint.attrs["units"] = "celsius"
    dewpoint.name = "dewpoint"
    return dewpoint


def calc_mixing_ratio(
    partial_pressure_data: xr.DataArray,
    total_pressure_data: xr.DataArray,
    molecular_weight_ratio: float = 0.6219569100577033,
) -> xr.DataArray:
    """
    Calculate the mixing ratio of a gas.

    This calculates mixing ratio given its partial pressure and the total pressure of the air.
    There are no required units for the input arrays, other than that they have the same units.

    Parameters
    ----------
    partial_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Partial pressure of the constituent gas.
    total_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Total air pressure.
    molecular_weight_ratio : :py:class:`float <float>`, optional.
        The ratio of the molecular weight of the constituent gas to that assumed for air.
        Defaults to the ratio for water vapor to dry air (:math:`\\epsilon\\approx0.622`).

    .. note::
        The units of `partial_pressure_data` and `total_pressure_data` should be the same.

    Returns
    -------
    The mixing ratio (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless (e.g. Kg/Kg or g/g).

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.mixing_ratio.html
    """
    return_data = (
        molecular_weight_ratio
        * partial_pressure_data
        / (total_pressure_data - partial_pressure_data)
    )
    return_data.attrs["units"] = "1"
    return_data.name = "mixing_ratio"
    return return_data


def calc_vapor_pressure(
    pressure_data: xr.DataArray,
    mixing_ratio_data: xr.DataArray,
    pressure_data_units: Literal["hPa", "Pa"] = None,
    epsilon: float = 0.6219569100577033,
) -> xr.DataArray:
    """
    Parameters
    ----------
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    mixing_ratio_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The mixing ratio of a gas.
    epsilon: :py:class:`float <float>`.
        The molecular weight ratio, which is molecular weight of the constituent gas to that assumed for air.
        Defaults to the ratio for water vapor to dry air. (:math:`\\epsilon \\approx 0.622`)
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.

    Returns
    -------
    The water vapor (partial) pressure (:py:class:`xarray.DataArray<xarray.DataArray>`), units according to `pressure_data_units`.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.vapor_pressure.html
    """
    return_data = pressure_data * mixing_ratio_data / (mixing_ratio_data + epsilon)

    if pressure_data_units is None:
        pressure_data_units = pressure_data.attrs["units"]
        return_data.attrs["units"] = f"{pressure_data_units}"
    else:
        return_data.attrs["units"] = f"{pressure_data_units}"
    return_data.name = "vapor_pressure"
    return return_data


def calc_saturation_vapor_pressure(
    temperature_data: xr.DataArray,
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
) -> xr.DataArray:
    """
    Calculate the saturation water vapor (partial) pressure.

    Parameters
    ----------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.

    Returns
    -------
    The mixing ratio (:py:class:`xarray.DataArray<xarray.DataArray>`), hPa.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.saturation_vapor_pressure.html
    """
    temperature_data = transfer_data_difference_units(
        input_data=temperature_data,
        input_units=temperature_data_units,
        output_units="celsius",
    )
    return_data = 6.112 * np.exp(17.67 * temperature_data / (temperature_data + 243.5))
    return_data.attrs["units"] = "hPa"
    return_data.name = "saturation_vapor_pressure"
    return return_data


def calc_saturation_mixing_ratio(
    total_pressure_data: xr.DataArray,
    temperature_data: xr.DataArray,
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    total_pressure_data_units: Literal["hPa", "Pa"],
) -> xr.DataArray:
    """
    Calculate the saturation mixing ratio of water vapor.

    This calculation is given total atmospheric pressure and air temperature.

    Parameters
    ----------
    total_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Total atmospheric pressure.
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    total_pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `total_pressure_data` value. Optional values are `hPa`, `Pa`.

    Returns
    -------
    The saturation mixing ratio (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.saturation_mixing_ratio.html
    """
    total_pressure_data = transfer_data_multiple_units(
        total_pressure_data, total_pressure_data_units, "hPa"
    )

    partial_pressure_data = calc_saturation_vapor_pressure(
        temperature_data=temperature_data, temperature_data_units=temperature_data_units
    )
    return_data = calc_mixing_ratio(
        partial_pressure_data=partial_pressure_data,
        total_pressure_data=total_pressure_data,
    )
    return_data.attrs["units"] = "1"
    return_data.name = "saturation_mixing_ratio"
    return return_data


def transfer_mixing_ratio_2_specific_humidity(
    mixing_ratio_data: xr.DataArray,
) -> xr.DataArray:
    """
    Calculate the specific humidity from mixing ratio.

    Parameters
    ----------
    mixing_ratio_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The mixing ratio of a gas.

    Returns
    -------
    The specific humidity (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless (e.g. Kg/Kg or g/g).

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.specific_humidity_from_mixing_ratio.html
    """
    return_data = mixing_ratio_data / (1 + mixing_ratio_data)
    return_data.attrs["units"] = "g/g"
    return_data.name = "specific_humidity"
    return return_data


def transfer_specific_humidity_2_mixing_ratio(
    specific_humidity_data: xr.DataArray,
    specific_humidity_units: Literal["kg/kg", "g/g", "g/kg"],
) -> xr.DataArray:
    """
    Calculate the mixing ratio from specific humidity.

    Parameters
    ----------
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The Specific humidity of air.
    specific_humidity_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/g`, `g/kg` and so on.

    Returns
    -------
    The mixing ratio (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.mixing_ratio_from_specific_humidity.html
    """
    specific_humidity_data = transfer_data_multiple_units(
        specific_humidity_data, specific_humidity_units, "g/g"
    )

    return_data = specific_humidity_data / (1 - specific_humidity_data)
    return_data.attrs["units"] = "1"
    return_data.name = "mixing_ratio"
    return return_data


def transfer_dewpoint_2_specific_humidity(
    dewpoint_data: xr.DataArray,
    pressure_data: xr.DataArray,
    dewpoint_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    pressure_data_units: Literal["hPa", "Pa"],
) -> xr.DataArray:
    """
    Calculate the specific humidity from the dewpoint temperature and pressure.

    Parameters
    ----------
    dewpoint_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The dewpoint temperature.
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    dewpoint_data_units: :py:class:`str <str>`.
        The unit corresponding to `dewpoint_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.

    Returns
    -------
    The specific humidity (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.specific_humidity_from_dewpoint.html
    """
    mixing_ratio = calc_saturation_mixing_ratio(
        total_pressure_data=pressure_data,
        temperature_data=dewpoint_data,
        temperature_data_units=dewpoint_data_units,
        total_pressure_data_units=pressure_data_units,
    )
    return_data = transfer_mixing_ratio_2_specific_humidity(mixing_ratio)
    return_data.attrs["units"] = "g/g"
    return_data.name = "specific_humidity"
    return return_data


def transfer_specific_humidity_2_dewpoint(
    specific_humidity_data: xr.DataArray,
    pressure_data: xr.DataArray,
    specific_humidity_units: Literal["kg/kg", "g/g", "g/kg"],
    pressure_data_units: Literal["hPa", "Pa"],
    epsilon: float = 0.6219569100577033,
) -> xr.DataArray:
    """
    Calculate the dewpoint from specific humidity and pressure.

    Parameters
    ----------
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    specific_humidity_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/g`, `g/kg` and so on.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
    epsilon: :py:class:`float <float>`.
        The molecular weight ratio, which is molecular weight of the constituent gas to that assumed for air.
        Defaults to the ratio for water vapor to dry air. (:math:`\\epsilon \\approx 0.622`)

    Returns
    -------
    The dewpoint (:py:class:`xarray.DataArray<xarray.DataArray>`), degrees Celsius.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.dewpoint_from_specific_humidity.html
    """
    specific_humidity_data = transfer_data_multiple_units(
        specific_humidity_data, specific_humidity_units, "g/g"
    )
    w = transfer_specific_humidity_2_mixing_ratio(
        specific_humidity_data=specific_humidity_data,
        specific_humidity_units=specific_humidity_units,
    )
    e = pressure_data * w / (epsilon + w)
    return_data = calc_dewpoint(
        vapor_pressure_data=e, vapor_pressure_data_units=pressure_data_units
    )
    return_data.attrs["units"] = "Celsius"
    return_data.name = "dewpoint"
    return return_data


def transfer_dewpoint_2_relative_humidity(
    temperature_data: xr.DataArray,
    dewpoint_data: xr.DataArray,
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    dewpoint_data_units: Literal["celsius", "kelvin", "fahrenheit"],
) -> xr.DataArray:
    """
    Calculate the relative humidity from dewpoint.

    Uses temperature and dewpoint to calculate relative humidity as the ratio of vapor pressure to saturation vapor pressures.

    Parameters
    ----------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    dewpoint_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The dewpoint temperature.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    dewpoint_data_units: :py:class:`str <str>`.
        The unit corresponding to `dewpoint_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.

    Returns
    -------
    The relative humidity (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.relative_humidity_from_dewpoint.html
    """
    temperature_data = transfer_data_difference_units(
        input_data=temperature_data,
        input_units=temperature_data_units,
        output_units="celsius",
    )
    dewpoint_data = transfer_data_difference_units(
        input_data=dewpoint_data,
        input_units=dewpoint_data_units,
        output_units="celsius",
    )

    e = calc_saturation_vapor_pressure(
        temperature_data=dewpoint_data, temperature_data_units=dewpoint_data_units
    )
    e_s = calc_saturation_vapor_pressure(
        temperature_data=temperature_data, temperature_data_units=temperature_data_units
    )
    return_data = e / e_s
    return_data.attrs["units"] = "1"
    return_data.name = "relative_humidity"
    return return_data


def transfer_mixing_ratio_2_relative_humidity(
    pressure_data: xr.DataArray,
    temperature_data: xr.DataArray,
    mixing_ratio_data: xr.DataArray,
    pressure_data_units: Literal["hPa", "Pa"],
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    epsilon: float = 0.6219569100577033,
) -> xr.DataArray:
    """
    Calculate the relative humidity from mixing ratio, temperature, and pressure.

    Parameters
    ----------
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    mixing_ratio_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The mixing ratio of a gas.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    epsilon: :py:class:`float <float>`.
        The molecular weight ratio, which is molecular weight of the constituent gas to that assumed for air.
        Defaults to the ratio for water vapor to dry air. (:math:`\\epsilon \\approx 0.622`)

    Returns
    -------
    The relative humidity (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.relative_humidity_from_mixing_ratio.html
    """
    w_s = calc_saturation_mixing_ratio(
        total_pressure_data=pressure_data,
        temperature_data=temperature_data,
        total_pressure_data_units=pressure_data_units,
        temperature_data_units=temperature_data_units,
    )
    return_data = (
        mixing_ratio_data / (epsilon + mixing_ratio_data) * (epsilon + w_s) / w_s
    )
    return_data.attrs["units"] = "1"
    return_data.name = "relative_humidity"
    return return_data


def transfer_specific_humidity_2_relative_humidity(
    pressure_data: xr.DataArray,
    temperature_data: xr.DataArray,
    specific_humidity_data: xr.DataArray,
    pressure_data_units: Literal["hPa", "Pa"],
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    specific_humidity_units: Literal["kg/kg", "g/g", "g/kg"],
) -> xr.DataArray:
    """
    Calculate the relative humidity from specific humidity, temperature, and pressure.

    Parameters
    ----------
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    specific_humidity_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/g`, `g/kg` and so on.

    Returns
    -------
    The relative humidity (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.relative_humidity_from_specific_humidity.html
    """
    mixing_ratio_data = transfer_specific_humidity_2_mixing_ratio(
        specific_humidity_data=specific_humidity_data,
        specific_humidity_units=specific_humidity_units,
    )
    return_data = transfer_mixing_ratio_2_relative_humidity(
        pressure_data=pressure_data,
        temperature_data=temperature_data,
        mixing_ratio_data=mixing_ratio_data,
        pressure_data_units=pressure_data_units,
        temperature_data_units=temperature_data_units,
    )
    return_data.attrs["units"] = "1"
    return_data.name = "relative_humidity"
    return return_data
