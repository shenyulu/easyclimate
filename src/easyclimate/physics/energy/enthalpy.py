"""
Atmospheric Enthalpy
"""

from __future__ import annotations
import xarray as xr
from typing import Literal
from ...core.utility import (
    transfer_data_temperature_units,
    transfer_data_multiple_units,
)
from ...core.utility import compare_multi_dataarray_coordinate

__all__ = ["calc_enthalpy"]


def calc_enthalpy(
    temperature_data: xr.DataArray,
    mixing_ratio_data: xr.DataArray,
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    mixing_ratio_data_units: Literal["kg/kg", "g/g", "g/kg"],
) -> xr.DataArray:
    """
    Calculate atmospheric enthalpy from temperature and humidity mixing ratio.

    Enthalpy is a thermodynamic quantity equivalent to the internal energy plus
    the energy the system exerts on its surroundings. The enthalpy is a constant
    pressure function. As such, it includes the work term for expansion
    against the atmosphere.

    .. math::

        T \\cdot (1.01 + 0.00189 \\cdot W) + 2.5 \\cdot W

    where the unit of :math:`T` (atmospheric temperature) is ``degC``, and :math:`W` (mixing ratio) is ``g/kg``.

    Parameters
    -----------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    mixing_ratio_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The mixing ratio of a gas.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    mixing_ratio_data_units: :py:class:`str <str>`.
        The unit corresponding to ``mixing_ratio_data`` value. Optional values are :math:`\\mathrm{kg/kg}`, :math:`\\mathrm{g/g}`, :math:`\\mathrm{g/kg}` and so on.

    Returns
    --------
    Atmospheric enthalpy ( :math:`\\mathrm{kJ/kg}` ).
        :py:class:`xarray.DataArray<xarray.DataArray>`.

    .. seealso::

        - https://www.ncl.ucar.edu/Document/Functions/Contributed/enthalpy.shtml
    """
    compare_multi_dataarray_coordinate([temperature_data, mixing_ratio_data])

    mixing_ratio_data = transfer_data_multiple_units(
        mixing_ratio_data, mixing_ratio_data_units, "g/kg"
    )  # g/kg
    temperature_data = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "degC"
    )  # degC

    # Calculate enthalpy
    e = (
        temperature_data * (1.01 + 0.00189 * mixing_ratio_data)
        + 2.5 * mixing_ratio_data
    )  # kJ/kg

    # clean other attrs
    e.attrs = dict()
    e.name = "enthalpy"
    e.attrs["long_name"] = "Enthalpy"
    e.attrs["units"] = "kJ/kg"
    e.attrs["equation"] = (
        "Enthalpy (h) kJ/kg = T*(1.01 + 0.00189*W) + 2.5*W; T=>degC, W=>g/kg"
    )

    return e
