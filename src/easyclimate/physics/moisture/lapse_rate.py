"""
Moist Adiabatic Lapse Rate
"""

from __future__ import annotations
import xarray as xr
from typing import Literal
from ...core.utility import (
    transfer_data_multiple_units,
    transfer_data_temperature_units,
)
from ..transfer import transfer_dewpoint_2_mixing_ratio

__all__ = ["calc_moist_adiabatic_lapse_rate"]


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
    dtdp : :py:class:`xarray.DataArray<xarray.DataArray>` ( :math:`\\mathrm{K/hPa}` ).
        Moist adiabatic lapse rate.
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
        dewpoint_data=temperature_data,
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
