"""
Virtual Temperature
"""

from __future__ import annotations
import numpy as np
import xarray as xr
import warnings
from typing import Literal
from ...core.utility import transfer_data_multiple_units

__all__ = ["calc_virtual_temperature", "calc_virtual_temperature_Hobbs2006"]


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

    Returns
    -------
    The virtual temperature, units according to ``temper_data``.
        :py:class:`xarray.DataArray<xarray.DataArray>`.

    Reference
    --------------
    - Doswell, C. A., and E. N. Rasmussen, 1994: The Effect of Neglecting the Virtual Temperature Correction on CAPE Calculations. Wea. Forecasting, 9, 625–629, https://journals.ametsoc.org/view/journals/wefo/9/4/1520-0434_1994_009_0625_teontv_2_0_co_2.xml
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

    where :math:`\\epsilon \\approx 0.622` when the mixing ratio (specific humidity) :math:`q` is expressed in :math:`\\mathrm{g \\ g^{-1}}`.

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

    Returns
    -------
    The virtual temperature, units according to ``temper_data``.
        :py:class:`xarray.DataArray<xarray.DataArray>`.

    Reference
    --------------
    - Hobbs, P. V., and J. M. Wallace, 2006: Atmospheric Science: An Introductory Survey. 2nd ed. Academic Press, 504 pp. https://www.sciencedirect.com/book/9780127329512/atmospheric-science
    - Doswell, C. A., and E. N. Rasmussen, 1994: The Effect of Neglecting the Virtual Temperature Correction on CAPE Calculations. Wea. Forecasting, 9, 625–629, https://journals.ametsoc.org/view/journals/wefo/9/4/1520-0434_1994_009_0625_teontv_2_0_co_2.xml
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
