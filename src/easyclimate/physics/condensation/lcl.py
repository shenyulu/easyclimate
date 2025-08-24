"""
Lifting Condensation Level (LCL)
"""

from __future__ import annotations
import numpy as np
import xarray as xr
from typing import Literal
from ...core.utility import (
    transfer_data_multiple_units,
    transfer_data_temperature_units,
)


__all__ = [
    "calc_lifting_condensation_level_bolton1980",
    "calc_lifting_condensation_level_Bohren_Albrecht2023",
]


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
        The unit corresponding to `temperature_data` value. Optional values are ``celsius``, ``kelvin``, ``fahrenheit``.
    relative_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `vapor_pressure_data` value. Optional values are ``%``, ``dimensionless``.

    Returns
    -------
    The lifting condensation level ( :math:`\\mathrm{K}` ).
        :py:class:`xarray.DataArray<xarray.DataArray>`.

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

    According to formulation (6.32) in Bohren & Albrecht (2023),

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
    The lifting condensation level (:py:class:`xarray.Dataset<xarray.Dataset>`)

    - p_lcl: lifting condensation level pressure ( :math:`\\mathrm{hPa}` ).
    - t_lcl: lifting condensation level temperature ( :math:`\\mathrm{K}` ).

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
