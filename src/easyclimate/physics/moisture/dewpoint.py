"""
Dewpoint
"""

from __future__ import annotations
import numpy as np
import xarray as xr
from typing import Literal
from ...core.utility import transfer_data_multiple_units

__all__ = ["calc_dewpoint"]


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
