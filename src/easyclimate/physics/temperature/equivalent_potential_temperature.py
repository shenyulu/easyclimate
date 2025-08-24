"""
Equivalent Potential Temperature
"""

from __future__ import annotations
import numpy as np
import xarray as xr
from typing import Literal
from ...core.utility import (
    transfer_data_multiple_units,
    transfer_data_temperature_units,
)
from ..transfer import transfer_dewpoint_2_mixing_ratio
from ..moisture.vapor_pressure import calc_saturation_vapor_pressure
from .potential_temperature import calc_potential_temperature

__all__ = [
    "calc_equivalent_potential_temperature",
    "calc_equivalent_potential_temperature_davies_jones2009",
]


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
        The dew point temperature.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    dewpoint_data_units: :py:class:`str <str>`.
        The unit corresponding to `dewpoint_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.

    Returns
    -------
    Equivalent potential temperature ( :math:`\\mathrm{K}` ).
        :py:class:`xarray.DataArray<xarray.DataArray>`.

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
        The dew point temperature.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    dewpoint_data_units: :py:class:`str <str>`.
        The unit corresponding to `dewpoint_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.

    Returns
    -------
    Equivalent potential temperature ( :math:`\\mathrm{K}` ).
        :py:class:`xarray.DataArray<xarray.DataArray>`.

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
