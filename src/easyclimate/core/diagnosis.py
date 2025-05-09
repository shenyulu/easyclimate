"""
Functions for Weather and climate variable diagnosis.

.. seealso::
    Stull, R. (2011): Meteorology for Scientists and Engineers. 3rd ed. Discount Textbooks, 924 pp. [Available online at https://www.eoas.ubc.ca/books/Practical_Meteorology/, https://www.eoas.ubc.ca/courses/atsc201/MSE3.html]
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
    transfer_data_temperature_units,
)
from easyclimate_backend.wet_bulb import _wet_bulb_temperature
from typing import Literal
import warnings

__all__ = [
    "calc_brunt_vaisala_frequency_atm",
    "get_coriolis_parameter",
    "calc_potential_temperature",
    "calc_potential_temperature_vertical",
    "calc_equivalent_potential_temperature",
    "calc_equivalent_potential_temperature_davies_jones2009",
    "calc_virtual_temperature_Hobbs2006",
    "calc_virtual_temperature",
    "calc_static_stability",
    "calc_dewpoint",
    "calc_mixing_ratio",
    "calc_vapor_pressure",
    "calc_saturation_vapor_pressure",
    "calc_saturation_mixing_ratio",
    "calc_wet_bulb_potential_temperature_iteration",
    "calc_wet_bulb_potential_temperature_davies_jones2008",
    "calc_wet_bulb_temperature_stull2011",
    "calc_wet_bulb_temperature_sadeghi2013",
    "calc_lifting_condensation_level_bolton1980",
    "calc_lifting_condensation_level_Bohren_Albrecht2023",
    "calc_moist_adiabatic_lapse_rate",
    "transfer_mixing_ratio_2_specific_humidity",
    "transfer_specific_humidity_2_mixing_ratio",
    "transfer_dewpoint_2_specific_humidity",
    "transfer_dewpoint_2_mixing_ratio",
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

    # clean other attrs
    N.attrs = dict()
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
    pressure_data: xr.DataArray,
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
    kappa: float = 287 / 1005.7,
) -> xr.DataArray:
    """
    Calculate the potential temperature for **dry air**.

    Uses the Poisson equation to calculation the potential temperature given pressure and temperature.

    .. math::
        \\theta = T \\left( \\frac{p_0}{p} \\right) ^\\kappa

    Parameters
    ----------
    temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Air temperature.
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
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
    - Bolton, D. (1980). The Computation of Equivalent Potential Temperature. Monthly Weather Review, 108(7), 1046-1053. https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_1046_tcoept_2_0_co_2.xml

    .. seealso::
        - `potential_temperature — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.potential_temperature.html>`__
        - `pot_temp - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/pot_temp.shtml>`__
    """
    pressure_data = transfer_data_multiple_units(
        pressure_data, pressure_data_units, "hPa"
    )
    pt = temper_data * (1000.0 / pressure_data) ** (kappa)
    return pt


def calc_potential_temperature_vertical(
    temper_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: str,
    kappa: float = 287 / 1005.7,
) -> xr.DataArray:
    """
    Calculate the potential temperature for vertical variables.

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

    # clean other attrs
    return_data.attrs = dict()
    try:
        return_data.attrs["units"] = temper_data.attrs["units"]
    except:
        return_data.attrs["units"] = "[T]"
        warnings.warn(
            "The variable of `temper_data` do not have `units` attribution, so the attribution of units in it is assigned to the symbol for dimension!"
        )
    return_data.name = "potential_temperature"
    return return_data


def calc_equivalent_potential_temperature(
    pressure_data: xr.DataArray,
    temperature_data: xr.DataArray,
    dewpoint_data: xr.DataArray,
    pressure_data_units: Literal["hPa", "Pa"],
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    dewpoint_data_units: Literal["celsius", "kelvin", "fahrenheit"],
) -> xr.DataArray:
    """
    Calculate equivalent potential temperature using Bolton (1980) approximation.


    Parameters
    ----------
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    dewpoint_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The dewpoint temperature.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    dewpoint_data_units: :py:class:`str <str>`.
        The unit corresponding to `dewpoint_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.

    Returns
    -------
    Equivalent potential temperature (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. seealso::
        - Bolton, D. (1980). The Computation of Equivalent Potential Temperature. Monthly Weather Review, 108(7), 1046-1053. https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_1046_tcoept_2_0_co_2.xml
    """
    pressure_data = transfer_data_multiple_units(
        pressure_data, pressure_data_units, "hPa"
    )
    temperature_data = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "K"
    )
    dewpoint_data = transfer_data_temperature_units(
        dewpoint_data, dewpoint_data_units, "K"
    )

    w = transfer_dewpoint_2_mixing_ratio(
        dewpoint_data=dewpoint_data,
        pressure_data=pressure_data,
        dewpoint_data_units="K",
        pressure_data_units="hPa",
    )

    e = calc_saturation_vapor_pressure(
        temperature_data=dewpoint_data, temperature_data_units="K"
    )

    theta_d = calc_potential_temperature(
        temper_data=temperature_data,
        pressure_data=pressure_data - e,
        pressure_data_units="hPa",
    )

    t_l = (
        1.0
        / (
            1.0 / (dewpoint_data - 56.0)
            + np.log(temperature_data / dewpoint_data) / 800.0
        )
        + 56.0
    )

    # see (24) in Bolton (1980)
    theta_dl = theta_d * ((temperature_data / t_l) ** (0.28 * w))

    # see (39) in Bolton (1980)
    theta_e = theta_dl * np.exp((3036.0 / t_l - 1.78) * w * (1 + 0.448 * w))

    # clean other attrs
    theta_e.attrs = dict()
    theta_e.name = "theta_e"
    theta_e.attrs["units"] = "K"
    return theta_e


def calc_equivalent_potential_temperature_davies_jones2009(
    pressure_data: xr.DataArray,
    temperature_data: xr.DataArray,
    dewpoint_data: xr.DataArray,
    pressure_data_units: Literal["hPa", "Pa"],
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    dewpoint_data_units: Literal["celsius", "kelvin", "fahrenheit"],
) -> xr.DataArray:
    """
    Calculate equivalent potential temperature using Robert Davies-Jones (2009) approximation.


    Parameters
    ----------
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    dewpoint_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The dewpoint temperature.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    dewpoint_data_units: :py:class:`str <str>`.
        The unit corresponding to `dewpoint_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.

    Returns
    -------
    Equivalent potential temperature (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. seealso::
        - Davies-Jones, R. (2009). On Formulas for Equivalent Potential Temperature. Monthly Weather Review, 137(9), 3137-3148. https://doi.org/10.1175/2009MWR2774.1
    """
    L0 = 2.56313 * 10**6
    L1 = 1754.0
    K2 = 1.137 * 10**6
    c_pd = 1005.7

    pressure_data = transfer_data_multiple_units(
        pressure_data, pressure_data_units, "hPa"
    )
    temperature_data = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "K"
    )
    dewpoint_data = transfer_data_temperature_units(
        dewpoint_data, dewpoint_data_units, "K"
    )

    w = transfer_dewpoint_2_mixing_ratio(
        dewpoint_data=dewpoint_data,
        pressure_data=pressure_data,
        dewpoint_data_units="K",
        pressure_data_units="hPa",
    )

    e = calc_saturation_vapor_pressure(
        temperature_data=dewpoint_data, temperature_data_units="K"
    )

    theta_d = calc_potential_temperature(
        temper_data=temperature_data,
        pressure_data=pressure_data - e,
        pressure_data_units="hPa",
    )

    t_l = (
        1.0
        / (
            1.0 / (dewpoint_data - 56.0)
            + np.log(temperature_data / dewpoint_data) / 800.0
        )
        + 56.0
    )

    # see (24) in Bolton (1980)
    theta_dl = theta_d * ((temperature_data / t_l) ** (0.28 * w))

    # see (6.5) in Robert Davies-Jones (2009)
    theta_e = theta_dl * np.exp((L0 - L1 * (t_l - 273.15) + K2 * w) * w / c_pd / t_l)

    # clean other attrs
    theta_e.attrs = dict()
    theta_e.name = "theta_e"
    theta_e.attrs["units"] = "K"
    return theta_e


def calc_virtual_temperature(
    temper_data: xr.DataArray,
    specific_humidity_data: xr.DataArray,
    specific_humidity_data_units: Literal["kg/kg", "g/g", "g/kg"],
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
    specific_humidity_data_units: :py:class:`str <str>`.
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
        specific_humidity_data, specific_humidity_data_units, "g/g"
    )

    T_v = (1 + epsilon * specific_humidity_data) * temper_data

    # clean other attrs
    T_v.attrs = dict()

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
    specific_humidity_data_units: Literal["kg/kg", "g/g", "g/kg"],
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
    specific_humidity_data_units: :py:class:`str <str>`.
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
        specific_humidity_data, specific_humidity_data_units, "g/g"
    )

    T_v = temper_data * (
        (specific_humidity_data + epsilon) / (epsilon * (1 + specific_humidity_data))
    )

    # clean other attrs
    T_v.attrs = dict()

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
    theta = calc_potential_temperature_vertical(
        temper_data, vertical_dim=vertical_dim, vertical_dim_units=vertical_dim_units
    )
    ln_theta = np.log(theta)
    part = calc_p_gradient(
        ln_theta, vertical_dim=vertical_dim, vertical_dim_units=vertical_dim_units
    )
    return_data = -temper_data * part

    # clean other attrs
    return_data.attrs = dict()
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
    vapor_pressure_data_units: Literal["hPa", "Pa", "mbar"],
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
    # clean other attrs
    dewpoint.attrs = dict()
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
    # clean other attrs
    return_data.attrs = dict()

    return_data.attrs["units"] = "dimensionless"
    return_data.attrs["standard_name"] = "humidity_mixing_ratio"
    return_data.name = "mixing_ratio"
    return return_data


def calc_vapor_pressure(
    pressure_data: xr.DataArray,
    mixing_ratio_data: xr.DataArray,
    pressure_data_units: Literal["hPa", "Pa", "mbar"] = None,
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

    # clean other attrs
    return_data.attrs = dict()
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
        - Bolton, D. (1980). The Computation of Equivalent Potential Temperature. Monthly Weather Review, 108(7), 1046-1053. https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_1046_tcoept_2_0_co_2.xml
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.saturation_vapor_pressure.html
    """
    temperature_data = transfer_data_temperature_units(
        input_data=temperature_data,
        input_units=temperature_data_units,
        output_units="celsius",
    )
    return_data = 6.112 * np.exp(17.67 * temperature_data / (temperature_data + 243.5))

    # clean other attrs
    return_data.attrs = dict()
    return_data.attrs["units"] = "hPa"
    return_data.name = "saturation_vapor_pressure"
    return return_data


def calc_saturation_mixing_ratio(
    total_pressure_data: xr.DataArray,
    temperature_data: xr.DataArray,
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    total_pressure_data_units: Literal["hPa", "Pa", "mbar"],
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

    # clean other attrs
    return_data.attrs = dict()
    return_data.attrs["units"] = "dimensionless"
    return_data.name = "saturation_mixing_ratio"
    return return_data


def calc_wet_bulb_potential_temperature_iteration(
    temperature_data: xr.DataArray,
    relative_humidity_data: xr.DataArray,
    pressure_data: xr.DataArray,
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    relative_humidity_data_units: Literal["%", "dimensionless"],
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
    A: float = 0.662 * 10 ** (-3),
    tolerance: float = 0.01,
    max_iter: int = 100,
) -> xr.DataArray:
    """
    Calculate wet-bulb potential temperature using iteration.

    The iterative formula

    .. math::

        e = e_{tw} - AP(t-t_{w})

    - :math:`e` is the water vapor pressure
    - :math:`e_{tw}` is the saturation water vapor pressure over a pure flat ice surface at wet-bulb temperature :math:`t_w` (when the wet-bulb thermometer is frozen, this becomes the saturation vapor pressure over a pure flat ice surface)
    - :math:`A` is the psychrometer constant
    - :math:`P` is the sea-level pressure
    - :math:`t` is the dry-bulb temperature
    - :math:`t_w` is the wet-bulb temperature

    Parameters
    ----------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    relative_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The relative humidity.
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    relative_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `vapor_pressure_data` value. Optional values are ``%``, ``dimensionless``.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
    A: :py:class:`float <float>`.
        Psychrometer coefficients.

        +-----------------------------------------+---------------------------------+-------------------------------+
        | Psychrometer Type and Ventilation Rate  | Wet Bulb Unfrozen (10^-3/°C^-1) | Wet Bulb Frozen (10^-3/°C^-1) |
        +=========================================+=================================+===============================+
        | Ventilated Psychrometer (2.5 m/s)       | 0.662                           | 0.584                         |
        +-----------------------------------------+---------------------------------+-------------------------------+
        | Spherical Psychrometer (0.4 m/s)        | 0.857                           | 0.756                         |
        +-----------------------------------------+---------------------------------+-------------------------------+
        | Cylindrical Psychrometer (0.4 m/s)      | 0.815                           | 0.719                         |
        +-----------------------------------------+---------------------------------+-------------------------------+
        | Chinese Spherical Psychrometer (0.8 m/s)| 0.7949                          | 0.7949                        |
        +-----------------------------------------+---------------------------------+-------------------------------+

    tolerance: :py:class:`float <float>`.
        Minimum acceptable deviation of the iterated value from the true value.
    max_iter: :py:class:`int <float>`.
        Maximum number of iterations.


    Returns
    ---------------
    tw: :py:class:`xarray.DataArray<xarray.DataArray>` ( :math:`\\mathrm{degC}` )
        Wet-bulb temperature

    Examples
    ----------------

    .. code:: python

        >>> import xarray as xr
        >>> import numpy as np

        # Create sample data
        >>> temp = xr.DataArray(np.array([20, 25, 30]), dims=['point'])
        >>> rh = xr.DataArray(np.array([50, 60, 70]), dims=['point'])
        >>> pressure = xr.DataArray(np.array([1000, 950, 900]), dims=['point'])

        # Calculate wet-bulb potential temperature
        >>> theta_w = calc_wet_bulb_potential_temperature_iteration(
        ...     temperature_data=temp,
        ...     relative_humidity_data=rh,
        ...     pressure_data=pressure,
        ...     temperature_data_units="celsius",
        ...     relative_humidity_data_units="%",
        ...     pressure_data_units="hPa"
        ... )

        # Example with 2D data
        >>> temp_2d = xr.DataArray(np.random.rand(10, 10) * 30, dims=['lat', 'lon'])
        >>> rh_2d = xr.DataArray(np.random.rand(10, 10) * 100, dims=['lat', 'lon'])
        >>> pres_2d = xr.DataArray(np.random.rand(10, 10) * 200 + 800, dims=['lat', 'lon'])
        >>> theta_w_2d = calc_wet_bulb_potential_temperature_iteration(
        ...     temp_2d, rh_2d, pres_2d, "celsius", "%", "hPa"
        ... )

    .. seealso::
        - Fan, J. (1987). Determination of the Psychrometer Coefficient A of the WMO Reference Psychrometer by Comparison with a Standard Gravimetric Hygrometer. Journal of Atmospheric and Oceanic Technology, 4(1), 239-244. https://journals.ametsoc.org/view/journals/atot/4/1/1520-0426_1987_004_0239_dotpco_2_0_co_2.xml
        - Wang Haijun. (2011). Two Wet-Bulb Temperature Estimation Methods and Error Analysis. Meteorological Monthly (Chinese), 37(4): 497-502. website: http://qxqk.nmc.cn/html/2011/4/20110415.html
        - Cheng Zhi, Wu Biwen, Zhu Baolin, et al, (2011). Wet-Bulb Temperature Looping Iterative Scheme and Its Application. Meteorological Monthly (Chinese), 37(1): 112-115. website: http://qxqk.nmc.cn/html/2011/1/20110115.html
    """
    pressure_data = transfer_data_multiple_units(
        pressure_data, pressure_data_units, "hPa"
    )
    temperature_data = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "degC"
    )
    relative_humidity_data = transfer_data_multiple_units(
        relative_humidity_data, relative_humidity_data_units, "%"
    )

    def _wet_bulb_iteration(t_dry, rh, p):
        t_w = _wet_bulb_temperature.wet_bulb_temperature(
            t_dry, rh, p, tolerance, A, max_iter
        )
        return np.array([t_w])

    result = xr.apply_ufunc(
        _wet_bulb_iteration,
        temperature_data,
        relative_humidity_data,
        pressure_data,
        vectorize=True,
        dask="parallelized",
        dask_gufunc_kwargs={
            "allow_rechunk": True,
        },
    )

    # clean other attrs
    result.attrs = dict()
    result.name = "tw"
    result.attrs["standard_name"] = "wet_bulb_temperature"
    result.attrs["units"] = "degC"

    return result


def calc_wet_bulb_potential_temperature_davies_jones2008(
    pressure_data: xr.DataArray,
    temperature_data: xr.DataArray,
    dewpoint_data: xr.DataArray,
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    dewpoint_data_units: Literal["celsius", "kelvin", "fahrenheit"],
) -> xr.DataArray:
    """
    Calculate wet-bulb potential temperature using Robert Davies-Jones (2008) approximation.

    Parameters
    ----------
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    dewpoint_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The dewpoint temperature.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    dewpoint_data_units: :py:class:`str <str>`.
        The unit corresponding to `dewpoint_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.

    Returns
    -------
    tw: :py:class:`xarray.DataArray<xarray.DataArray>` ( :math:`\\mathrm{K}` )
        Wet-bulb temperature

    .. seealso::

        - Davies-Jones, R. (2008). An Efficient and Accurate Method for Computing the Wet-Bulb Temperature along Pseudoadiabats. Monthly Weather Review, 136(7), 2764-2785. https://doi.org/10.1175/2007MWR2224.1
        - Knox, J. A., Nevius, D. S., & Knox, P. N. (2017). Two Simple and Accurate Approximations for Wet-Bulb Temperature in Moist Conditions, with Forecasting Applications. Bulletin of the American Meteorological Society, 98(9), 1897-1906. https://doi.org/10.1175/BAMS-D-16-0246.1
    """
    pressure_data = transfer_data_multiple_units(
        pressure_data, pressure_data_units, "hPa"
    )
    temperature_data = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "K"
    )
    dewpoint_data = transfer_data_temperature_units(
        dewpoint_data, dewpoint_data_units, "K"
    )

    # Calculate the equivalent potential temperature
    theta_e = calc_equivalent_potential_temperature(
        pressure_data=pressure_data,
        temperature_data=temperature_data,
        dewpoint_data=dewpoint_data,
        pressure_data_units="hPa",
        temperature_data_units="K",
        dewpoint_data_units="K",
    )  # units: K

    # Create the resulting DataArray, keeping the same coordinates and properties
    theta_w = xr.full_like(theta_e, np.nan)

    # Calculate x value
    x = theta_e / 273.15

    # see (3.8) in Davies-Jones, R. (2008)
    # Create mask condition
    mask = theta_e >= 173.15

    # The points that satisfy the conditions are calculated
    x_masked = x.where(mask)
    x2 = x_masked * x_masked
    x3 = x2 * x_masked
    x4 = x3 * x_masked

    a = 7.101574 - 20.68208 * x_masked + 16.11182 * x2 + 2.574631 * x3 - 5.205688 * x4
    b = 1 - 3.552497 * x_masked + 3.781782 * x2 - 0.6899655 * x3 - 0.5929340 * x4

    # Assign the result to theta_w
    theta_w = xr.where(mask, theta_e - np.exp(a / b), theta_e)

    # clean other attrs
    theta_w.attrs = dict()
    theta_w.name = "tw"
    theta_w.attrs["standard_name"] = "wet_bulb_temperature"
    theta_w.attrs["units"] = "K"
    return theta_w


def calc_wet_bulb_temperature_stull2011(
    temperature_data: xr.DataArray,
    relative_humidity_data: xr.DataArray,
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    relative_humidity_data_units: Literal["%", "dimensionless"],
) -> xr.DataArray:
    """
    Calculate wet-bulb temperature using Stull (2011) empirical formula.

    .. math::
        T_{w} =T\\operatorname{atan}[0.151977(\\mathrm{RH}\\%+8.313659)^{1/2}]+\\operatorname{atan}(T+\\mathrm{RH}\%)-\\operatorname{atan}(\\mathrm{RH}\\%-1.676331)
        +0.00391838(\\mathrm{RH}\\%)^{3/2}\\operatorname{atan}(0.023101\\mathrm{RH}\\%)-4.686035.

    .. tip::

        This methodology was not valid for ambient conditions with low values of :math:`T_a` (dry-bulb temperature; i.e., <10°C),
        and/or with low values of RH  (5% < RH < 10%).
        The Stull methodology was also only valid at sea level.

    Parameters
    ----------------------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    relative_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The relative humidity.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    relative_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `vapor_pressure_data` value. Optional values are ``%``, ``dimensionless``.

    Returns
    ----------------------
    tw: :py:class:`xarray.DataArray<xarray.DataArray>` ( :math:`\\mathrm{K}` )
        Wet-bulb temperature

    .. seealso::
        - Stull, R. (2011). Wet-Bulb Temperature from Relative Humidity and Air Temperature. Journal of Applied Meteorology and Climatology, 50(11), 2267-2269. https://doi.org/10.1175/JAMC-D-11-0143.1
        - Stull, R. (2011): Meteorology for Scientists and Engineers. 3rd ed. Discount Textbooks, 924 pp. [Available online at https://www.eoas.ubc.ca/books/Practical_Meteorology/, https://www.eoas.ubc.ca/courses/atsc201/MSE3.html]
        - Knox, J. A., Nevius, D. S., & Knox, P. N. (2017). Two Simple and Accurate Approximations for Wet-Bulb Temperature in Moist Conditions, with Forecasting Applications. Bulletin of the American Meteorological Society, 98(9), 1897-1906. https://doi.org/10.1175/BAMS-D-16-0246.1
    """
    # Convert to Celsius
    t_c = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "degC"
    )
    # Convert to `%`
    relative_humidity_data = transfer_data_multiple_units(
        relative_humidity_data, relative_humidity_data_units, "%"
    )

    tw_c = (
        t_c * np.arctan(0.151977 * np.sqrt(relative_humidity_data + 8.313659))
        + np.arctan(t_c + relative_humidity_data)
        - np.arctan(relative_humidity_data - 1.676331)
        + 0.00391838
        * (relative_humidity_data**1.5)
        * np.arctan(0.023101 * relative_humidity_data)
        - 4.686035
    )

    # Convert back to Kelvin
    tw = transfer_data_temperature_units(tw_c, "degC", "K")

    # clean other attrs
    tw.attrs = dict()
    tw.name = "tw"
    tw.attrs["standard_name"] = "wet_bulb_temperature"
    tw.attrs["units"] = "K"
    return tw


def calc_wet_bulb_temperature_sadeghi2013(
    temperature_data: xr.DataArray,
    height_data: xr.DataArray,
    relative_humidity_data: xr.DataArray,
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    height_data_units: Literal["m", "km"],
    relative_humidity_data_units: Literal["%", "dimensionless"],
) -> xr.DataArray:
    """
    Calculate wet-bulb temperature using Sadeghi et. al (2011) empirical formula.

    Parameters
    ----------------------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    height_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The elevation.
    relative_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The relative humidity.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    height_data_units: :py:class:`str <str>`.
        The unit corresponding to `height_data` value. Optional values are `m`, `km`.
    relative_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `vapor_pressure_data` value. Optional values are ``%``, ``dimensionless``.

    Returns
    ----------------------
    tw: :py:class:`xarray.DataArray<xarray.DataArray>` ( :math:`\\mathrm{degC}` )
        Wet-bulb temperature

    .. seealso::
        - Sadeghi, S., Peters, T. R., Cobos, D. R., Loescher, H. W., & Campbell, C. S. (2013). Direct Calculation of Thermodynamic Wet-Bulb Temperature as a Function of Pressure and Elevation. Journal of Atmospheric and Oceanic Technology, 30(8), 1757-1765. https://doi.org/10.1175/JTECH-D-12-00191.1
    """
    gamma = 0.4 * 0.001  # units: g/(kg K) => g/(g K)
    a = 0.611

    # T_a: dry-bulb temperature (K)
    temperature_data = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "degC"
    )
    # H: elevation (m).
    height_data = transfer_data_multiple_units(height_data, height_data_units, "m")
    relative_humidity_data = transfer_data_multiple_units(
        relative_humidity_data, relative_humidity_data_units, "dimensionless"
    )

    # saturation vapor pressure
    e_s_Ta = calc_saturation_vapor_pressure(
        temperature_data=temperature_data, temperature_data_units="degC"
    )  # units: hPa
    # vapor_pressure
    e_a_Ta = e_s_Ta * relative_humidity_data
    e_a = transfer_data_multiple_units(e_a_Ta, "hPa", "kPa")

    # see (11)
    lambda_value = 0.0014 * np.exp(0.027 * temperature_data)
    # see (12)
    zeta = (
        (-3 * 10 ** (-7) * temperature_data**3)
        - (10 ** (-5)) * temperature_data**2
        + (2 * 10 ** (-5)) * temperature_data
        + 4.44 * 10 ** (-2)
    )
    # see (3c)
    p_a = 101.3 * np.exp(-height_data / 8200.0)  # meter for `height_data`
    # see (9c)
    phi = zeta + gamma * p_a
    # see (9a)
    psi = a - gamma * p_a * temperature_data - e_a

    # see (8b)
    delta = phi**2 - 4 * lambda_value * psi

    # see (8a)
    tw = (-phi + np.sqrt(delta)) / (2 * lambda_value)

    # clean other attrs
    tw.attrs = dict()
    tw.name = "tw"
    tw.attrs["standard_name"] = "wet_bulb_temperature"
    tw.attrs["units"] = "degC"
    return tw


def calc_lifting_condensation_level_bolton1980(
    temperature_data: xr.DataArray,
    relative_humidity_data: xr.DataArray,
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    relative_humidity_data_units: Literal["%", "dimensionless"],
) -> xr.DataArray:
    """
    Calculate lifting condensation level using Bolton (1980) approximation.

    .. math::

        T_L = \\frac{1}{\\frac{1}{T_K - 55} - \\frac{\\ln (U/100)}{2840}} + 55


    Parameters
    ----------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    relative_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The relative humidity.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    relative_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `vapor_pressure_data` value. Optional values are ``%``, ``dimensionless``.

    Returns
    -------
    The lifting condensation level (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. seealso::
        - Bolton, D. (1980). The Computation of Equivalent Potential Temperature. Monthly Weather Review, 108(7), 1046-1053. https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_1046_tcoept_2_0_co_2.xml
    """
    t_k = transfer_data_temperature_units(temperature_data, temperature_data_units, "K")
    U = transfer_data_multiple_units(
        relative_humidity_data, relative_humidity_data_units, "%"
    )

    # see (22) in Bolton, D. (1980)
    part = (1 / t_k - 55) - (np.log(U / 100) / 2840)
    t_l = 1 / part + 55

    # clean other attrs
    t_l.attrs = dict()
    t_l.attrs["standard_name"] = "atmosphere_lifting_condensation_level"
    t_l.attrs["units"] = "K"
    return t_l


def calc_lifting_condensation_level_Bohren_Albrecht2023(
    pressure_data: xr.DataArray,
    temperature_data: xr.DataArray,
    dewpoint_data: xr.DataArray,
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    dewpoint_data_units: Literal["celsius", "kelvin", "fahrenheit"],
) -> xr.Dataset:
    """
    Calculate lifting condensation level using Bohren & Albrecht (2023) approximation.

    According to  formulation (6.32) in Bohren & Albrecht (2023),

    .. math::

        T_{LCL}=\\frac{1-AT_d}{1/T_d + B\\ln(T/T_d)-A}

    Where

    .. math::

        A=-\\left(\\frac{c_{pv}-c_{pw}}{l_r}-\\frac{c_{pd}}{\\epsilon l_v}\\right),\\quad B=\\frac{c_{pd}}{\\epsilon l_v}

    and

    .. math::

        l_v=l_{vr} - (c_{pv}-c_{pw})T_r, \\quad l_v(T_d)=l_r+(c_{pd} - c_{pw}) T_d


    Parameters
    ----------
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    dewpoint_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The dewpoint temperature.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    dewpoint_data_units: :py:class:`str <str>`.
        The unit corresponding to `dewpoint_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.

    Returns
    -------
    The lifting condensation level (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. seealso::
        Bohren, C. F., and B. A. Albrecht, 2023: Atmospheric Thermodynamics Second Edition. Oxford University Press, 579 pp. Website: http://gen.lib.rus.ec/book/index.php?md5=AA3B25841BE3AEBA2628EF9961F58C52
    """
    # Constants
    AIR_R_d = 287.0  # Dry air gas constant [J/K/kg]
    AIR_Cp_d = 1004.0  # Dry air isobaric specific heat [J/K/kg]
    WATER_Cp_v = 1885.0  # Water vapor isobaric specific heat [J/K/kg]
    WATER_Cp_l = 4186.0  # Liquid water isobaric specific heat [J/K/kg]
    WATER_Lv_0c = 2.50084e6  # Latent heat of vaporization at 0°C [J/kg]

    # Derived constants
    AIR_KAPPA_d = AIR_R_d / AIR_Cp_d  # Rd/Cp of dry air

    p = transfer_data_multiple_units(pressure_data, pressure_data_units, "hPa")
    t = transfer_data_temperature_units(temperature_data, temperature_data_units, "K")
    td = transfer_data_temperature_units(dewpoint_data, dewpoint_data_units, "K")

    lr = WATER_Lv_0c - (WATER_Cp_v - WATER_Cp_l) * 273.15
    a = -(WATER_Cp_v - WATER_Cp_l) / lr * AIR_Cp_d / (AIR_KAPPA_d * lr)
    b = AIR_Cp_d / (AIR_KAPPA_d * lr)

    t_lcl = (1 - a * td) / (1 / td + b * np.log(t / td) - a)
    p_lcl = p * (t_lcl / t) ** (1.0 / AIR_KAPPA_d)

    result = xr.Dataset()
    result["p_lcl"] = p_lcl
    result["t_lcl"] = t_lcl
    result["p_lcl"].attrs["units"] = "hPa"
    result["t_lcl"].attrs["units"] = "K"
    return result


def calc_moist_adiabatic_lapse_rate(
    pressure_data: xr.DataArray,
    temperature_data: xr.DataArray,
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
) -> xr.DataArray:
    """
    Calculate moist adiabatic lapse rate.

    Parameters
    -----------
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.

    Returns:
    --------
    dtdp : :py:class:`xarray.DataArray<xarray.DataArray>`.
        Moist adiabatic lapse rate [K/hPa]
    """
    # Constants
    AIR_R_d = 287.0  # Dry air gas constant [J/K/kg]
    AIR_Cp_d = 1004.0  # Dry air isobaric specific heat [J/K/kg]
    WATER_R_v = 461.50  # Water vapor gas constant [J/K/kg]
    WATER_Lv_0c = 2.50084e6  # Latent heat of vaporization at 0°C [J/kg]

    p = transfer_data_multiple_units(pressure_data, pressure_data_units, "hPa")
    t = transfer_data_temperature_units(temperature_data, temperature_data_units, "K")

    # Saturated air parcel
    w_s = transfer_dewpoint_2_mixing_ratio(
        pressure_data=pressure_data,
        temperature_data=temperature_data,
        pressure_data_units="hPa",
        dewpoint_data_units="K",
    )  # dimensionless
    numerator = AIR_R_d * t + WATER_Lv_0c * w_s
    denominator = AIR_Cp_d + (WATER_Lv_0c**2 * w_s / (WATER_R_v * t**2))

    result = (1.0 / p) * (numerator / denominator)

    # clean other attrs
    result.attrs = dict()
    result.attrs["units"] = "K/hPa"
    return result


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

    # clean other attrs
    return_data.attrs = dict()
    return_data.attrs["units"] = "g/g"
    return_data.name = "specific_humidity"
    return return_data


def transfer_specific_humidity_2_mixing_ratio(
    specific_humidity_data: xr.DataArray,
    specific_humidity_data_units: Literal["kg/kg", "g/g", "g/kg"],
) -> xr.DataArray:
    """
    Calculate the mixing ratio from specific humidity.

    Parameters
    ----------
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The Specific humidity of air.
    specific_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/g`, `g/kg` and so on.

    Returns
    -------
    The mixing ratio (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.mixing_ratio_from_specific_humidity.html
    """
    specific_humidity_data = transfer_data_multiple_units(
        specific_humidity_data, specific_humidity_data_units, "g/g"
    )

    return_data = specific_humidity_data / (1 - specific_humidity_data)

    # clean other attrs
    return_data.attrs = dict()
    return_data.attrs["units"] = "dimensionless"
    return_data.name = "mixing_ratio"
    return return_data


def transfer_dewpoint_2_specific_humidity(
    dewpoint_data: xr.DataArray,
    pressure_data: xr.DataArray,
    dewpoint_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
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
    # clean other attrs
    return_data.attrs = dict()
    return_data.attrs["units"] = "g/g"
    return_data.name = "specific_humidity"
    return return_data


def transfer_dewpoint_2_mixing_ratio(
    dewpoint_data: xr.DataArray,
    pressure_data: xr.DataArray,
    dewpoint_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
):
    """
    Calculate the mixing ratio from the dewpoint temperature and pressure.

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
    The mixing ratio (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless.
    """
    q = transfer_dewpoint_2_specific_humidity(
        dewpoint_data=dewpoint_data,
        pressure_data=pressure_data,
        dewpoint_data_units=dewpoint_data_units,
        pressure_data_units=pressure_data_units,
    )

    mixing_ratio = transfer_specific_humidity_2_mixing_ratio(
        specific_humidity_data=q, specific_humidity_data_units="g/g"
    )

    return mixing_ratio


def transfer_specific_humidity_2_dewpoint(
    specific_humidity_data: xr.DataArray,
    pressure_data: xr.DataArray,
    specific_humidity_data_units: Literal["kg/kg", "g/g", "g/kg"],
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
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
    specific_humidity_data_units: :py:class:`str <str>`.
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
        specific_humidity_data, specific_humidity_data_units, "g/g"
    )
    w = transfer_specific_humidity_2_mixing_ratio(
        specific_humidity_data=specific_humidity_data,
        specific_humidity_data_units=specific_humidity_data_units,
    )
    e = pressure_data * w / (epsilon + w)
    return_data = calc_dewpoint(
        vapor_pressure_data=e, vapor_pressure_data_units=pressure_data_units
    )

    # clean other attrs
    return_data.attrs = dict()
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
    temperature_data = transfer_data_temperature_units(
        input_data=temperature_data,
        input_units=temperature_data_units,
        output_units="celsius",
    )
    dewpoint_data = transfer_data_temperature_units(
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

    # clean other attrs
    return_data.attrs = dict()
    return_data.attrs["units"] = "dimensionless"
    return_data.name = "relative_humidity"
    return return_data


def transfer_mixing_ratio_2_relative_humidity(
    pressure_data: xr.DataArray,
    temperature_data: xr.DataArray,
    mixing_ratio_data: xr.DataArray,
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
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

    # clean other attrs
    return_data.attrs = dict()
    return_data.attrs["units"] = "dimensionless"
    return_data.name = "relative_humidity"
    return return_data


def transfer_specific_humidity_2_relative_humidity(
    pressure_data: xr.DataArray,
    temperature_data: xr.DataArray,
    specific_humidity_data: xr.DataArray,
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    specific_humidity_data_units: Literal["kg/kg", "g/g", "g/kg"],
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
    specific_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/g`, `g/kg` and so on.

    Returns
    -------
    The relative humidity (:py:class:`xarray.DataArray<xarray.DataArray>`), dimensionless.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.relative_humidity_from_specific_humidity.html
    """
    mixing_ratio_data = transfer_specific_humidity_2_mixing_ratio(
        specific_humidity_data=specific_humidity_data,
        specific_humidity_data_units=specific_humidity_data_units,
    )
    return_data = transfer_mixing_ratio_2_relative_humidity(
        pressure_data=pressure_data,
        temperature_data=temperature_data,
        mixing_ratio_data=mixing_ratio_data,
        pressure_data_units=pressure_data_units,
        temperature_data_units=temperature_data_units,
    )

    # clean other attrs
    return_data.attrs = dict()
    return_data.attrs["units"] = "dimensionless"
    return_data.name = "relative_humidity"
    return return_data
