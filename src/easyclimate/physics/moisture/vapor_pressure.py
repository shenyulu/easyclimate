"""
Vapor Pressure
"""

from __future__ import annotations
import numpy as np
import xarray as xr
import warnings
from typing import Literal
from ...core.utility import transfer_data_temperature_units

__all__ = ["calc_vapor_pressure", "calc_saturation_vapor_pressure"]


def calc_vapor_pressure(
    pressure_data: xr.DataArray,
    mixing_ratio_data: xr.DataArray,
    pressure_data_units: Literal["hPa", "Pa", "mbar"] = None,
    epsilon: float = 0.6219569100577033,
) -> xr.DataArray:
    """
    Calculate the vapor pressure.

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
    The water vapor (partial) pressure, units according to ``pressure_data_units``.
        :py:class:`xarray.DataArray<xarray.DataArray>`

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
    The saturation water vapor (partial) pressure ( :math:`\\mathrm{hPa}` ).
        :py:class:`xarray.DataArray<xarray.DataArray>`.

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
