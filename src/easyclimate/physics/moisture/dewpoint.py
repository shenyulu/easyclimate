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
    Calculate the ambient dew point temperature given the vapor pressure.

    This function inverts the Bolton (1980) formula for saturation vapor
    pressure to instead calculate the temperature. This yields the following formula for
    dewpoint in degrees Celsius, where :math:`e` is the ambient vapor pressure in millibars:

    .. math::

        T = \frac{243.5 \log(e / 6.112)}{17.67 - \log(e / 6.112)}

    Parameters
    ----------
    vapor_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Water vapor partial pressure.
    vapor_pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `vapor_pressure_data` value. Optional values are `hPa`, `Pa`.

    Returns
    -------
    The dew point ( :math:`\\mathrm{degC}` ).
        :py:class:`xarray.DataArray<xarray.DataArray>`

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.dewpoint.html
        - Bolton, D. (1980). The Computation of Equivalent Potential Temperature. Monthly Weather Review, 108(7), 1046-1053. https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_1046_tcoept_2_0_co_2.xml
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
