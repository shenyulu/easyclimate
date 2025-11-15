"""
Functions for handling calculations of units
"""

from __future__ import annotations

import xarray as xr
import numpy as np

__all__ = [
    "transfer_units_coeff",
    "transfer_data_multiple_units",
    "transfer_data_difference_units",
    "transfer_data_temperature_units",
    "transfer_data_units",
]


def transfer_units_coeff(input_units: str, output_units: str) -> float:
    """
    Compute the unit conversion factor from input units to output units.

    Parameters
    ----------
    input_units : :py:class:`str <str>`
        The input unit string (e.g., 'm').
    output_units : :py:class:`str <str>`
        The output unit string (e.g., 'km').

    Returns
    -------
    :py:class:`float <float>`
        The conversion factor to multiply the input values by to obtain output values.

    Example
    -------
    >>> import easyclimate as ecl
    >>> result1 = ecl.transfer_units_coeff("m/s", "km/h")
    >>> result2 = ecl.transfer_units_coeff("hPa", "mbar")
    >>> result3 = ecl.transfer_units_coeff("mm/day", "m/month")
    >>> print(result1, result2, result3)
    3.6 1.0 0.0304375
    """
    from pint import UnitRegistry

    ureg = UnitRegistry()

    base_unitmul = ureg(input_units).to(output_units).to_tuple()
    base = base_unitmul[0]
    return base


def transfer_data_multiple_units(
    input_data: xr.DataArray | xr.Dataset, input_units: str, output_units: str
) -> xr.DataArray | xr.Dataset:
    """
    Convert data units for multiplicative transitions (e.g., m to km).

    Parameters
    ----------
    input_data : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`
        The input data with units to convert.
    input_units : :py:class:`str <str>`
        The input unit string. Below is a table of supported temperature unit standards and aliases:
    output_units : :py:class:`str <str>`
        The output unit string.

    Returns
    -------
    :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`
        The converted data with updated 'units' attribute.

    Example
    -------
    >>> import xarray as xr
    >>> import numpy as np
    >>> import easyclimate as ecl
    >>> ds = xr.DataArray(
    ...     np.array([[3e-5, 5.4e-5], [5.2, -75.5]]),
    ...     dims=("lon", "lat"),
    ...     coords={"lon": np.array([160, 70]), "lat": np.array([87.5, -87.5])}
    ... )
    >>> result = ecl.transfer_data_multiple_units(ds, "mm/day", "m/day")
    >>> print(result)
    <xarray.DataArray (lon: 2, lat: 2)> Size: 32B
    array([[ 3.00e-08,  5.40e-08],
        [ 5.20e-03, -7.55e-02]])
    Coordinates:
    * lon      (lon) int64 16B 160 70
    * lat      (lat) float64 16B 87.5 -87.5
    Attributes:
        units:    m/day
    """
    base = transfer_units_coeff(input_units, output_units)
    output_data = input_data * base
    output_data.attrs["units"] = output_units
    return output_data


def transfer_data_difference_units(
    input_data: xr.DataArray | xr.Dataset, input_units: str, output_units: str
) -> xr.DataArray | xr.Dataset:
    """
    Convert data units for difference-based transitions (e.g., adjustments requiring offset).

    Parameters
    ----------
    input_data : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`
        The input data with units to convert.
    input_units : :py:class:`str <str>`
        The input unit string.
    output_units : :py:class:`str <str>`
        The output unit string.

    Returns
    -------
    :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`
        The converted data with updated 'units' attribute.

    Example
    -------
    >>> import xarray as xr
    >>> import easyclimate as ecl
    >>> result = ecl.transfer_data_difference_units(xr.DataArray(15), "celsius", "kelvin")
    >>> print(result)
    <xarray.DataArray ()> Size: 8B
    array(288.15)
    Attributes:
        units:    kelvin
    """
    coeff = transfer_units_coeff(input_units, output_units)
    output_data = input_data + coeff - 1
    output_data.attrs["units"] = output_units
    return output_data


def transfer_data_temperature_units(
    input_data: xr.DataArray | xr.Dataset, input_units: str, output_units: str
) -> xr.DataArray | xr.Dataset:
    """
    Convert temperature data from one unit to another, supporting aliases.

    Parameters
    ----------
    input_data : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`
        The input temperature data to convert.
    input_units : :py:class:`str <str>`
        The input temperature unit (e.g., 'degC', 'K', '摄氏度'). The available values are as follows:

    .. tab-set::

        .. tab-item:: Kelvin scale

            K, kelvin, degK, 开氏度, ケルビン

        .. tab-item:: Celsius scale

            degC, celsius, degC, 摄氏度, セルシウス度

        .. tab-item:: Fahrenheit scale

            degF, fahrenheit, 华氏度, カ氏度, 華氏度, ファーレンハイト度

        .. tab-item:: Rankine scale

            degR, rankine, 兰氏度, ランキン度

        .. tab-item:: Réaumur scale

            degRe, reaumur, 列氏度, レオミュール度

    output_units : :py:class:`str <str>`
        The output temperature unit (e.g., 'degF', 'kelvin'). The available values are as follows:

    .. tab-set::

        .. tab-item:: Kelvin scale

            K, kelvin, degK, 开氏度, ケルビン

        .. tab-item:: Celsius scale

            degC, celsius, degC, 摄氏度, セルシウス度

        .. tab-item:: Fahrenheit scale

            degF, fahrenheit, 华氏度, カ氏度, 華氏度, ファーレンハイト度

        .. tab-item:: Rankine scale

            degR, rankine, 兰氏度, ランキン度

        .. tab-item:: Réaumur scale

            degRe, reaumur, 列氏度, レオミュール度

    Returns
    -------
    :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`
        The converted temperature data with updated 'units' attribute.

    Example
    -------
    >>> import xarray as xr
    >>> import easyclimate as ecl
    >>> result = ecl.transfer_data_temperature_units(
    ...     xr.DataArray([104, 100, 92, 92, 86, 80, 80, 60, 30]), "摄氏度", "カ氏度"
    ... )
    >>> print(result)
    <xarray.DataArray (dim_0: 9)> Size: 72B
    array([219.2, 212. , 197.6, 197.6, 186.8, 176. , 176. , 140. ,  86. ])
    Dimensions without coordinates: dim_0
    Attributes:
        units:    カ氏度

    .. seealso::
        - https://en.wikipedia.org/wiki/Kelvin
        - https://en.wikipedia.org/wiki/Celsius
        - https://en.wikipedia.org/wiki/Fahrenheit
        - https://en.wikipedia.org/wiki/Rankine_scale
        - https://en.wikipedia.org/wiki/R%C3%A9aumur_scale
        - :py:func:`transfer_data_units <transfer_data_units>`
    """
    # 定义标准单位及其别名
    unit_aliases = {
        "K": ["kelvin", "K", "degK", "开氏度", "ケルビン"],
        "degC": ["celsius", "degC", "摄氏度", "セルシウス度"],
        "degF": [
            "fahrenheit",
            "degF",
            "华氏度",
            "カ氏度",
            "華氏度",
            "ファーレンハイト度",
        ],
        "degR": ["rankine", "degR", "兰氏度", "ランキン度"],
        "degRe": ["reaumur", "degRe", "列氏度", "レオミュール度"],
    }

    def normalize_unit(unit, unit_aliases):
        """将用户输入的单位名称标准化为标准单位。"""
        unit_lower = unit.lower()
        for std_unit, aliases in unit_aliases.items():
            if unit_lower in [alias.lower() for alias in aliases]:
                return std_unit
        raise ValueError(f"Unsupported units: {unit}")

    # Standardized unit name
    # 标准化单位名称
    from_unit_std = normalize_unit(input_units, unit_aliases)
    to_unit_std = normalize_unit(output_units, unit_aliases)

    # Define unit conversion table (scale factor, offset)
    # 定义单位转换表（缩放因子, 偏移量）
    unit_conversions = {
        "degC": {
            "degC": (1, 0),
            "degF": (1.8, 32),
            "K": (1, 273.15),
            "degR": (1.8, 491.67),  # °R = (°C + 273.15) * 1.8
            "degRe": (0.8, 0),  # °Re = °C * 0.8
        },
        "degF": {
            "degF": (1, 0),
            "degC": (5 / 9, -32 * 5 / 9),
            "K": (5 / 9, 255.372222),  # K = (°F + 459.67) * 5/9
            "degR": (1, 459.67),  # °R = °F + 459.67
            "degRe": (4 / 9, -32 * 4 / 9),  # °Re = (°F - 32) * 4/9
        },
        "K": {
            "K": (1, 0),
            "degC": (1, -273.15),
            "degF": (1.8, -459.67),
            "degR": (1.8, 0),  # °R = K * 1.8
            "degRe": (0.8, -273.15 * 0.8),  # °Re = (K - 273.15) * 0.8
        },
        "degR": {
            "degR": (1, 0),
            "K": (5 / 9, 0),  # K = °R * 5/9
            "degC": (5 / 9, -273.15),  # °C = (°R - 491.67) * 5/9
            "degF": (1, -459.67),  # °F = °R - 459.67
            "degRe": (4 / 9, -491.67 * 4 / 9),  # °Re = (°R - 491.67) * 4/9
        },
        "degRe": {
            "degRe": (1, 0),
            "degC": (1.25, 0),  # °C = °Re * 1.25
            "degF": (2.25, 32),  # °F = °Re * 2.25 + 32
            "K": (1.25, 273.15),  # K = °Re * 1.25 + 273.15
            "degR": (2.25, 491.67),  # °R = (°Re * 1.25 + 273.15) * 1.8
        },
    }

    # Check if the conversion is supported
    # 检查是否支持该转换
    if (
        from_unit_std not in unit_conversions
        or to_unit_std not in unit_conversions[from_unit_std]
    ):
        raise ValueError(
            f"Conversion from {input_units} to {output_units} is not supported"
        )

    # Get scale factor and offset
    # 获取缩放因子和偏移量
    scaling_factor, offset = unit_conversions[from_unit_std][to_unit_std]

    # Perform conversion
    # 执行转换
    result = np.multiply(input_data, scaling_factor) + offset
    if isinstance(result, (xr.DataArray, xr.Dataset)):
        result.attrs["units"] = output_units
    return result


def transfer_data_units(
    input_data: xr.DataArray | xr.Dataset, input_units: str, output_units: str
) -> xr.DataArray | xr.Dataset:
    """
    Convert data units for any type of transition using Pint.

    .. warning::

        Does NOT support ``dask``.

    Parameters
    ----------
    input_data : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`
        The input data with units to convert.
    input_units : :py:class:`str <str>`
        The input unit string.
    output_units : :py:class:`str <str>`
        The output unit string.

    Returns
    -------
    :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`xarray.Dataset <xarray.Dataset>`
        The converted data with updated 'units' attribute.

    Example
    -------
    >>> import xarray as xr
    >>> import easyclimate as ecl
    >>> result = ecl.transfer_data_units(
    ...     xr.DataArray([104, 100, 92, 92, 86, 80, 80, 60, 30]), "degC", "degF"
    ... )
    >>> print(result)
    <xarray.DataArray (dim_0: 9)> Size: 72B
    array([219.2, 212. , 197.6, 197.6, 186.8, 176. , 176. , 140. ,  86. ])
    Dimensions without coordinates: dim_0
    Attributes:
        units:    degF

    .. seealso::
        :py:func:`transfer_data_temperature_units <transfer_data_temperature_units>`
    """
    # Create a global unit registry
    import pint

    ureg = pint.UnitRegistry()
    Q_ = ureg.Quantity

    def transfer_unit(value, from_unit, to_unit):
        """
        Converts values from one unit to another.

        Parameters
        ----------
        value (float):
            The value to be converted.
        from_unit (str):
            The original unit (e.g. 'degC').
        to_unit (str):
            The target unit (e.g. 'degF').

        Returns
        -------
            float: The converted value.
        """
        quantity = Q_(value, ureg(from_unit))
        converted_quantity = quantity.to(to_unit)
        return converted_quantity.magnitude

    result = xr.apply_ufunc(
        transfer_unit,
        input_data,
        input_units,
        output_units,
        vectorize=True,
        dask="parallelized",
    )

    result.attrs = input_data.attrs
    result.attrs["units"] = output_units
    return result
