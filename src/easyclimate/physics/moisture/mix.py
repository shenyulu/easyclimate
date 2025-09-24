"""
Mixing Ratio of a Gas
"""

from __future__ import annotations
import xarray as xr
from typing import Literal
from ...core.utility import transfer_data_multiple_units
from .vapor_pressure import calc_saturation_vapor_pressure

__all__ = ["calc_mixing_ratio", "calc_saturation_mixing_ratio"]


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
    The mixing ratio ( :math:`\\mathrm{g/g}` ).
        :py:class:`xarray.DataArray<xarray.DataArray>`.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.mixing_ratio.html
    """
    return_data = (
        molecular_weight_ratio
        * partial_pressure_data
        / (total_pressure_data - partial_pressure_data)
    )  # dimensionless, g/g

    # clean other attrs
    return_data.attrs = dict()
    return_data.attrs["units"] = "g/g"
    return_data.attrs["standard_name"] = "humidity_mixing_ratio"
    return_data.name = "mixing_ratio"
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
    The saturation mixing ratio ( :math:`\\mathrm{g/g}` ).
        :py:class:`xarray.DataArray<xarray.DataArray>`.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.saturation_mixing_ratio.html
    """
    total_pressure_data = transfer_data_multiple_units(
        total_pressure_data, total_pressure_data_units, "hPa"
    )  # hPa

    partial_pressure_data = calc_saturation_vapor_pressure(
        temperature_data=temperature_data, temperature_data_units=temperature_data_units
    )  # hPa
    return_data = calc_mixing_ratio(
        partial_pressure_data=partial_pressure_data,
        total_pressure_data=total_pressure_data,
    )  # dimensionless, g/g

    # clean other attrs
    return_data.attrs = dict()
    return_data.attrs["units"] = "g/g"
    return_data.name = "saturation_mixing_ratio"
    return return_data
