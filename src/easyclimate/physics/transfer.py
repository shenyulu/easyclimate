"""
Unit conversion for Meteorological variables
"""

from __future__ import annotations
import numpy as np
import xarray as xr
import warnings
from typing import Literal
from ..core.utility import transfer_data_multiple_units, transfer_data_temperature_units
from .moisture.mix import calc_saturation_mixing_ratio
from .moisture.dewpoint import calc_dewpoint
from .moisture.vapor_pressure import calc_saturation_vapor_pressure

__all__ = [
    "transfer_mixing_ratio_2_specific_humidity",
    "transfer_specific_humidity_2_mixing_ratio",
    "transfer_dewpoint_2_specific_humidity",
    "transfer_dewpoint_2_mixing_ratio",
    "transfer_specific_humidity_2_dewpoint",
    "transfer_dewpoint_2_relative_humidity",
    "transfer_mixing_ratio_2_relative_humidity",
    "transfer_specific_humidity_2_relative_humidity",
    "transfer_relative_humidity_2_dewpoint",
]


def transfer_mixing_ratio_2_specific_humidity(
    mixing_ratio_data: xr.DataArray,
    mixing_ratio_data_units: Literal["kg/kg", "g/g", "g/kg"],
) -> xr.DataArray:
    """
    Calculate the specific humidity from mixing ratio.

    Parameters
    ----------
    mixing_ratio_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The mixing ratio of a gas.
    mixing_ratio_data_units: :py:class:`str <str>`.
        The unit corresponding to ``mixing_ratio_data`` value. Optional values are :math:`\\mathrm{kg/kg}`, :math:`\\mathrm{g/g}`, :math:`\\mathrm{g/kg}` and so on.


    Returns
    -------
    The specific humidity, dimensionless (e.g. :math:`\\mathrm{kg/kg}`, :math:`\\mathrm{g/g}`).
        :py:class:`xarray.DataArray<xarray.DataArray>`.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.specific_humidity_from_mixing_ratio.html
    """
    mixing_ratio_data = transfer_data_multiple_units(
        mixing_ratio_data, mixing_ratio_data_units, "g/g"
    )  # g/g
    return_data = mixing_ratio_data / (1 + mixing_ratio_data)

    # clean other attrs
    return_data.attrs = dict()
    return_data.attrs["units"] = "g/g"
    return_data.name = "specific_humidity"
    return return_data


def transfer_specific_humidity_2_mixing_ratio(
    specific_humidity_data: xr.DataArray,
    specific_humidity_data_units: Literal["kg/kg", "g/g", "g/kg"],
) -> xr.DataArray:
    """
    Calculate the mixing ratio from specific humidity.

    Parameters
    ----------
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The Specific humidity of air.
    specific_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are :math:`\\mathrm{kg/kg}`, :math:`\\mathrm{g/g}`, :math:`\\mathrm{g/kg}` and so on.

    Returns
    -------
    The mixing ratio, dimensionless (e.g. :math:`\\mathrm{kg/kg}`, :math:`\\mathrm{g/g}`).
        :py:class:`xarray.DataArray<xarray.DataArray>`.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.mixing_ratio_from_specific_humidity.html
    """
    specific_humidity_data = transfer_data_multiple_units(
        specific_humidity_data, specific_humidity_data_units, "g/g"
    )  # g/g

    return_data = specific_humidity_data / (1 - specific_humidity_data)  # g/g

    # clean other attrs
    return_data.attrs = dict()
    return_data.attrs["units"] = "g/g"
    return_data.name = "mixing_ratio"
    return return_data


def transfer_dewpoint_2_specific_humidity(
    dewpoint_data: xr.DataArray,
    pressure_data: xr.DataArray,
    dewpoint_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
) -> xr.DataArray:
    """
    Calculate the specific humidity from the dew point temperature and pressure.

    Parameters
    ----------
    dewpoint_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The dew point temperature.
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    dewpoint_data_units: :py:class:`str <str>`.
        The unit corresponding to `dewpoint_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.

    Returns
    -------
    The specific humidity, dimensionless (e.g. :math:`\\mathrm{kg/kg}`, :math:`\\mathrm{g/g}`).
        :py:class:`xarray.DataArray<xarray.DataArray>`.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.specific_humidity_from_dewpoint.html
    """
    mixing_ratio = calc_saturation_mixing_ratio(
        total_pressure_data=pressure_data,
        temperature_data=dewpoint_data,
        temperature_data_units=dewpoint_data_units,
        total_pressure_data_units=pressure_data_units,
    )  # g/g
    return_data = transfer_mixing_ratio_2_specific_humidity(mixing_ratio, "g/g")  # g/g
    # clean other attrs
    return_data.attrs = dict()
    return_data.attrs["units"] = "g/g"
    return_data.name = "specific_humidity"
    return return_data


def transfer_dewpoint_2_mixing_ratio(
    dewpoint_data: xr.DataArray,
    pressure_data: xr.DataArray,
    dewpoint_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
):
    """
    Calculate the mixing ratio from the dew point temperature and pressure.

    Parameters
    ----------
    dewpoint_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The dew point temperature.
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    dewpoint_data_units: :py:class:`str <str>`.
        The unit corresponding to `dewpoint_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.

    Returns
    -------
    The mixing ratio, dimensionless (e.g. :math:`\\mathrm{kg/kg}`, :math:`\\mathrm{g/g}`).
        :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    q = transfer_dewpoint_2_specific_humidity(
        dewpoint_data=dewpoint_data,
        pressure_data=pressure_data,
        dewpoint_data_units=dewpoint_data_units,
        pressure_data_units=pressure_data_units,
    )  # g/g

    mixing_ratio = transfer_specific_humidity_2_mixing_ratio(
        specific_humidity_data=q, specific_humidity_data_units="g/g"
    )  # g/g

    return mixing_ratio


def transfer_specific_humidity_2_dewpoint(
    specific_humidity_data: xr.DataArray,
    pressure_data: xr.DataArray,
    specific_humidity_data_units: Literal["kg/kg", "g/g", "g/kg"],
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
    epsilon: float = 0.6219569100577033,
) -> xr.DataArray:
    """
    Calculate the dew point temperature from specific humidity and pressure.

    Parameters
    ----------
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    specific_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are :math:`\\mathrm{kg/kg}`, :math:`\\mathrm{g/g}`, :math:`\\mathrm{g/kg}` and so on.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
    epsilon: :py:class:`float <float>`.
        The molecular weight ratio, which is molecular weight of the constituent gas to that assumed for air.
        Defaults to the ratio for water vapor to dry air. (:math:`\\epsilon \\approx 0.622`)

    Returns
    -------
    The dew point temperature ( :math:`\\mathrm{degC}` ).
        :py:class:`xarray.DataArray<xarray.DataArray>`.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.dewpoint_from_specific_humidity.html
    """
    specific_humidity_data = transfer_data_multiple_units(
        specific_humidity_data, specific_humidity_data_units, "g/g"
    )  # g/g
    w = transfer_specific_humidity_2_mixing_ratio(
        specific_humidity_data=specific_humidity_data,
        specific_humidity_data_units="g/g",
    )  # g/g
    e = pressure_data * w / (epsilon + w)
    return_data = calc_dewpoint(
        vapor_pressure_data=e, vapor_pressure_data_units=pressure_data_units
    )  # Celsius

    # clean other attrs
    return_data.attrs = dict()
    return_data.attrs["units"] = "Celsius"
    return_data.name = "dewpoint"
    return return_data


def transfer_dewpoint_2_relative_humidity(
    temperature_data: xr.DataArray,
    dewpoint_data: xr.DataArray,
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    dewpoint_data_units: Literal["celsius", "kelvin", "fahrenheit"],
) -> xr.DataArray:
    """
    Calculate the relative humidity from dew point temperature.

    Uses temperature and dew point temperature to calculate relative humidity as the ratio of vapor pressure to saturation vapor pressures.

    Parameters
    ----------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    dewpoint_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The dew point temperature.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    dewpoint_data_units: :py:class:`str <str>`.
        The unit corresponding to `dewpoint_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.

    Returns
    -------
    The relative humidity, dimensionless.
        :py:class:`xarray.DataArray<xarray.DataArray>`.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.relative_humidity_from_dewpoint.html
    """
    temperature_data = transfer_data_temperature_units(
        input_data=temperature_data,
        input_units=temperature_data_units,
        output_units="celsius",
    )  # Celsius
    dewpoint_data = transfer_data_temperature_units(
        input_data=dewpoint_data,
        input_units=dewpoint_data_units,
        output_units="celsius",
    )  # Celsius

    e = calc_saturation_vapor_pressure(
        temperature_data=dewpoint_data, temperature_data_units="celsius"
    )  # hPa
    e_s = calc_saturation_vapor_pressure(
        temperature_data=temperature_data, temperature_data_units="celsius"
    )  # hPa
    return_data = e / e_s

    # clean other attrs
    return_data.attrs = dict()
    return_data.attrs["units"] = "dimensionless"
    return_data.name = "relative_humidity"
    return return_data


def transfer_mixing_ratio_2_relative_humidity(
    pressure_data: xr.DataArray,
    temperature_data: xr.DataArray,
    mixing_ratio_data: xr.DataArray,
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    mixing_ratio_data_units: Literal["kg/kg", "g/g", "g/kg"],
    epsilon: float = 0.6219569100577033,
) -> xr.DataArray:
    """
    Calculate the relative humidity from mixing ratio, temperature, and pressure.

    Parameters
    ----------
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    mixing_ratio_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The mixing ratio of a gas.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    mixing_ratio_data_units: :py:class:`str <str>`.
        The unit corresponding to ``mixing_ratio_data`` value. Optional values are :math:`\\mathrm{kg/kg}`, :math:`\\mathrm{g/g}`, :math:`\\mathrm{g/kg}` and so on.
    epsilon: :py:class:`float <float>`.
        The molecular weight ratio, which is molecular weight of the constituent gas to that assumed for air.
        Defaults to the ratio for water vapor to dry air. (:math:`\\epsilon \\approx 0.622`)

    Returns
    -------
    The relative humidity, dimensionless.
        :py:class:`xarray.DataArray<xarray.DataArray>`.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.relative_humidity_from_mixing_ratio.html
    """
    w_s = calc_saturation_mixing_ratio(
        total_pressure_data=pressure_data,
        temperature_data=temperature_data,
        total_pressure_data_units=pressure_data_units,
        temperature_data_units=temperature_data_units,
    )  # g/g
    mixing_ratio_data = transfer_data_multiple_units(
        mixing_ratio_data, mixing_ratio_data_units, "g/g"
    )  # g/g
    return_data = (
        mixing_ratio_data / (epsilon + mixing_ratio_data) * (epsilon + w_s) / w_s
    )  # dimensionless

    # clean other attrs
    return_data.attrs = dict()
    return_data.attrs["units"] = "dimensionless"
    return_data.name = "relative_humidity"
    return return_data


def transfer_specific_humidity_2_relative_humidity(
    pressure_data: xr.DataArray,
    temperature_data: xr.DataArray,
    specific_humidity_data: xr.DataArray,
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    specific_humidity_data_units: Literal["kg/kg", "g/g", "g/kg"],
) -> xr.DataArray:
    """
    Calculate the relative humidity from specific humidity, temperature, and pressure.

    Parameters
    ----------
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    specific_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are :math:`\\mathrm{kg/kg}`, :math:`\\mathrm{g/g}`, :math:`\\mathrm{g/kg}` and so on.

    Returns
    -------
    The relative humidity, dimensionless.
        :py:class:`xarray.DataArray<xarray.DataArray>`.

    .. seealso::
        - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.relative_humidity_from_specific_humidity.html
    """
    mixing_ratio_data = transfer_specific_humidity_2_mixing_ratio(
        specific_humidity_data=specific_humidity_data,
        specific_humidity_data_units=specific_humidity_data_units,
    )  # g/g
    return_data = transfer_mixing_ratio_2_relative_humidity(
        pressure_data=pressure_data,
        temperature_data=temperature_data,
        mixing_ratio_data=mixing_ratio_data,
        pressure_data_units=pressure_data_units,
        temperature_data_units=temperature_data_units,
        mixing_ratio_data_units="g/g",
    )  # dimensionless

    # clean other attrs
    return_data.attrs = dict()
    return_data.attrs["units"] = "dimensionless"
    return_data.name = "relative_humidity"
    return return_data


def transfer_relative_humidity_2_dewpoint(
    relative_humidity_data: xr.DataArray,
    temperature_data: xr.DataArray,
    relative_humidity_data_units: Literal["%", "dimensionless"],
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
) -> xr.DataArray:
    """
    Calculate dew point temperature from temperature and relative humidity.

    The dew point temperature given temperature and relative humidity using the equations from John Dutton's "Ceaseless Wind" (pp 273-274). Missing values are ignored.

    The dew point temperature :math:`T_d` (in Kelvin) is calculated from temperature :math:`T` and relative humidity :math:`RH` using the formula:

    .. math::

        T_d = \\frac{T \\cdot L}{L - T \\cdot \\ln(RH/100)}, \\quad \\text{where} \\quad L = \\frac{597.3 - 0.57(T - 273.0)}{GCX} \\quad \\text{and} \\quad GCX = \\frac{461.5}{4186}.


    Parameters
    ----------
    relative_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The relative humidity.
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    relative_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `vapor_pressure_data` value. Optional values are ``%``, ``dimensionless``.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to ``temperature_data`` value. Optional values are ``celsius``, ``kelvin``, ``fahrenheit``.

    Returns
    -------
    dewpoint : :py:class:`xarray.DataArray<xarray.DataArray>` ( :math:`\\mathrm{K}` )
        Dew point temperature.

    Reference
    --------------
    - Dutton, J. A. (1976). The Ceaseless Wind: An introduction to the theory of atmospheric motion. McGraw-Hill, Inc. https://libgen.li/file.php?md5=7154c67714c7aee56152cfda528e2080

    .. seealso::

        https://www.ncl.ucar.edu/Document/Functions/Built-in/dewtemp_trh.shtml
    """
    # Constants
    GC = 461.5  # Gas constant for water vapor [J/(kg·K)]
    GCX = GC / (1000.0 * 4.186)  # Convert to [cal/(g·K)]

    tk = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "K"
    )  # K
    rh = transfer_data_multiple_units(
        relative_humidity_data, relative_humidity_data_units, "%"
    )  # %

    # Calculate latent heat of vaporization
    lhv = (597.3 - 0.57 * (tk - 273.0)) / GCX

    # Calculate dew point temperature
    tdk = tk * lhv / (lhv - tk * np.log(rh * 0.01))  # K

    # clean other attrs
    tdk.attrs = dict()
    tdk.name = "dewpoint"
    tdk.attrs["standard_name"] = "dew_point_temperature"
    tdk.attrs["units"] = "K"
    return tdk
