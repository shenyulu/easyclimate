"""
Potential Temperature
"""

from __future__ import annotations
import numpy as np
import xarray as xr
import warnings
from typing import Literal
from ...core.utility import transfer_data_multiple_units, transfer_units_coeff


__all__ = ["calc_potential_temperature", "calc_potential_temperature_vertical"]


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
    Potential temperature, units according to ``temper_data``.
        :py:class:`xarray.DataArray<xarray.DataArray>`.

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

    # clean other attrs
    pt.attrs = dict()
    try:
        pt.attrs["units"] = temper_data.attrs["units"]
    except:
        pt.attrs["units"] = "[T]"
        warnings.warn(
            "The variable of `temper_data` do not have `units` attribution, so the attribution of units in it is assigned to the symbol for dimension!"
        )
    return pt


def calc_potential_temperature_vertical(
    temper_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
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
    Potential temperature, units according to ``temper_data``.
        :py:class:`xarray.DataArray<xarray.DataArray>`.

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
