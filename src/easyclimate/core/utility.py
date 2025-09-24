"""
Functions for package utility.
"""

from __future__ import annotations

import numpy as np
import xarray as xr
import warnings

from packaging import version
from functools import wraps
from typing import Union, List, Tuple, Optional
from .datanode import DataNode
from ..version import __version__

__all__ = [
    "assert_compared_version",
    "find_dims_axis",
    "transfer_int2datetime",
    "split_datetime2yearday",
    "transfer_deg2rad",
    "transfer_inf2nan",
    "transfer_nan2value",
    "get_weighted_spatial_data",
    "get_compress_xarraydata",
    "transfer_dFdp2dFdz",
    "sort_ascending_latlon_coordinates",
    "transfer_units_coeff",
    "transfer_data_multiple_units",
    "transfer_data_difference_units",
    "transfer_data_temperature_units",
    "transfer_data_units",
    "generate_dataset_dispatcher",
    "generate_datatree_dispatcher",
    "transfer_xarray_lon_from180TO360",
    "transfer_xarray_lon_from360TO180",
    "module_available",
    "dequantify_metpy_xarraydata",
    "reverse_bool_xarraydata",
    "compare_two_dataarray_coordinate",
    "compare_multi_dataarray_coordinate",
    "validate_dataarrays",
    "deprecated",
    "datetime_to_numeric",
    "numeric_to_datetime",
    "calculate_time_steps",
    "clean_extra_coords",
]


def assert_compared_version(ver1: float, ver2: float) -> int:
    """
    Compare python library versions.

    .. attention::
        - Only for incoming version numbers without alphabetic characters.
        - Based on this method, the version number comparison should result in the following `"10.12.2.6.5">"10.12.2.6"`.

    Parameters
    ----------
    - ver1: :py:class:`float <float>`, version number 1
    - ver2: :py:class:`float <float>`, version number 2

    Returns
    -------
    :py:class:`int<int>`.

    .. note::
        If `ver1<ver2`, return `-1`; If `ver1=ver2`, return `0`; If `ver1>ver2`, return `1`.

    Examples
    --------

    .. code:: python

        >>> import easyclimate as ecl
        >>> result = ecl.assert_compared_version("10.12.2.6.5", "10.12.2.6")
        >>> print(result)
        1
    """
    list1 = str(ver1).split(".")
    list2 = str(ver2).split(".")
    # print(list1)
    # print(list2)
    # `len` of a list with a short loop count.
    for i in range(len(list1)) if len(list1) < len(list2) else range(len(list2)):
        if int(list1[i]) == int(list2[i]):
            pass
        elif int(list1[i]) < int(list2[i]):
            return -1
        else:
            return 1
    # End of loop, which list is long which version number is high.
    if len(list1) == len(list2):
        return 0
    elif len(list1) < len(list2):
        return -1
    else:
        return 1


def find_dims_axis(data: xr.DataArray, dim: str) -> int:
    """
    Find the index of `dim` in the xarray DataArray.

    Parameters
    ----------
    - data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
    - dim : :py:class:`str <str>`
        Dimension(s) over which to find axis.

    Returns
    -------
    :py:class:`int <int>`.
    """
    return data.dims.index(dim)


def transfer_int2datetime(data: np.array) -> np.datetime64:
    """
    Convert a numpy array of years of type integer to `np.datetime64` type.

    Parameters
    ----------
    - data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

    Examples
    --------

    .. code:: python

        >>> import easyclimate as ecl
        >>> import numpy as np
        >>> intyear = np.array([2054, 2061, 2062, 2067, 2071, 2075, 2076, 2078, 2085, 2089, 2096])
        >>> ecl.transfer_int2datetime(intyear)
        array(['2054-01-01T00:00:00.000000000', '2061-01-01T00:00:00.000000000',
               '2062-01-01T00:00:00.000000000', '2067-01-01T00:00:00.000000000',
               '2071-01-01T00:00:00.000000000', '2075-01-01T00:00:00.000000000',
               '2076-01-01T00:00:00.000000000', '2078-01-01T00:00:00.000000000',
               '2085-01-01T00:00:00.000000000', '2089-01-01T00:00:00.000000000',
               '2096-01-01T00:00:00.000000000'], dtype='datetime64[ns]')

    .. seealso::
        `Python(pandas)整数类型数据转换为时间类型 <https://www.jianshu.com/p/d12d95fbc90c>`__.
    """
    import pandas as pd

    # xarray coordinate axis does not accept DatetimeIndex, so use `.to_numpy()` to convert it to numpy array.
    return pd.to_datetime(data, format="%Y").to_numpy()


def split_datetime2yearday(ds: xr.DataArray, time_dim: str = "time") -> xr.DataArray:
    """
    Convert `np.datetime64` type with years and days to `year` and `day` coordinates.

    Parameters
    ----------
    data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.

    .. seealso::
        `Function in xarray to regroup monthly data into months and # of years <https://github.com/pydata/xarray/discussions/5119>`__.
    """
    year = ds[time_dim].dt.year
    month = ds[time_dim].dt.month

    # assign new coords
    ds = ds.assign_coords(year=(time_dim, year.data), month=(time_dim, month.data))

    # reshape the array to (..., "month", "year")
    return ds.set_index({time_dim: ("year", "month")}).unstack(time_dim)


def transfer_deg2rad(ds: xr.DataArray) -> xr.DataArray:
    """
    Convert Degrees to Radians.

    Parameters
    ----------
    - ds: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Degrees data.

    Returns
    -------
    - Radians data.: :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return ds * np.pi / 180


def transfer_inf2nan(ds: xr.DataArray) -> xr.DataArray:
    """
    Convert `np.inf` in `ds` to `np.nan`, respectively.

    Parameters
    ----------
    - ds: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Data include `np.inf`.

    Returns
    -------
    - Data include `np.nan`.: :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return ds.where(np.isfinite(ds), np.nan)


def transfer_nan2value(
    ds: xr.DataArray,
    value: float,
) -> xr.DataArray:
    """
    Convert `np.inf` in `ds` to `np.nan`, respectively.

    Parameters
    ----------
    - ds: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Data include `np.inf`.

    Returns
    -------
    - Data include `np.nan`.: :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return ds.fillna(value)


def get_weighted_spatial_data(
    data_input: xr.DataArray,
    lat_dim: str = "lat",
    lon_dim: str = "lon",
    method: str = "cos_lat",
) -> xr.DataArray:
    """
    Get the area-weighting data.

    Parameters
    ----------
    - data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
    - lat_dim: :py:class:`str <str>`.
        Latitude dimension over which to apply. By default is applied over the `lat` dimension.
    - lon_dim: :py:class:`str <str>`.
        Longitude dimension over which to apply. By default is applied over the `lon` dimension.
    - method: {`'cos_lat'`, `'area'`}.
        area-weighting methods.

        1. `'cos_lat'`: weighting data by the cosine of latitude.
        2. `'area'`: weighting data by area, where you weight each data point by the area of each grid cell.

    .. Caution::
        - `data_input` must be **regular lonlat grid**.
        - If you are calculating global average temperature just on land,
          then you need to mask out the ocean in your area dataset at first.

    .. seealso::
        - `The Correct Way to Average the Globe (Why area-weighting your data is important) <https://towardsdatascience.com/the-correct-way-to-average-the-globe-92ceecd172b7>`__.
        - Kevin Cowtan, Peter Jacobs, Peter Thorne, Richard Wilkinson,
          Statistical analysis of coverage error in simple global temperature estimators,
          Dynamics and Statistics of the Climate System, Volume 3, Issue 1, 2018, dzy003, https://doi.org/10.1093/climsys/dzy003.

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_basic_statistical_analysis.py
    """
    if method == "cos_lat":
        weights = np.cos(np.deg2rad(data_input[lat_dim]))
        return data_input.weighted(weights)
    elif method == "area":
        # Source: https://towardsdatascience.com/the-correct-way-to-average-the-globe-92ceecd172b7
        # Author: Luke Gloege
        def earth_radius(lat):
            """
            # Source: https://gist.githubusercontent.com/lgloege/6377b0d418982d2ec1c19d17c251f90e/raw/4742b104a250bfd9114436f6ef7135fe76a8ab3f/earth_radius.py

            calculate radius of Earth assuming oblate spheroid
            defined by WGS84

            Input
            ---------
            lat: vector or latitudes in degrees

            Output
            ----------
            r: vector of radius in meters

            Notes
            -----------
            WGS84: https://earth-info.nga.mil/GandG/publications/tr8350.2/tr8350.2-a/Chapter%203.pdf
            """
            # define oblate spheroid from WGS84
            a = 6378137
            b = 6356752.3142
            e2 = 1 - (b**2 / a**2)

            # convert from geodecic to geocentric
            # see equation 3-110 in WGS84
            lat = np.deg2rad(lat)
            lat_gc = np.arctan((1 - e2) * np.tan(lat))

            # radius equation
            # see equation 3-107 in WGS84
            r = (a * (1 - e2) ** 0.5) / (1 - (e2 * np.cos(lat_gc) ** 2)) ** 0.5
            return r

        def area_grid(lat, lon, lat_dim, lon_dim):
            """
            # Source: https://gist.githubusercontent.com/lgloege/6377b0d418982d2ec1c19d17c251f90e/raw/4742b104a250bfd9114436f6ef7135fe76a8ab3f/area_grid.py

            Calculate the area of each grid cell
            Area is in square meters

            Input
            -----------
            lat: vector of latitude in degrees
            lon: vector of longitude in degrees

            Output
            -----------
            area: grid-cell area in square-meters with dimensions, [lat,lon]

            Notes
            -----------
            Based on the function in
            https://github.com/chadagreene/CDT/blob/master/cdt/cdtarea.m
            """
            xlon, ylat = np.meshgrid(lon, lat)
            R = earth_radius(ylat)

            dlat = np.deg2rad(np.gradient(ylat, axis=0))
            dlon = np.deg2rad(np.gradient(xlon, axis=1))

            dy = dlat * R
            dx = dlon * R * np.cos(np.deg2rad(ylat))
            area = dy * dx

            xda = xr.DataArray(
                area,
                dims=[lat_dim, lon_dim],
                coords={lat_dim: lat, lon_dim: lon},
                attrs={
                    "long_name": "area_per_pixel",
                    "description": "area per pixel",
                    "units": "m^2",
                },
            )
            return xda

        # area dataArray
        lat_array = data_input[lat_dim].astype("float64").data
        lon_array = data_input[lon_dim].astype("float64").data
        da_area = area_grid(lat_array, lon_array, lat_dim, lon_dim)
        # total area
        total_area = da_area.sum(dim=(lat_dim, lon_dim))
        weights = da_area / total_area
        return data_input.weighted(weights)


def get_compress_xarraydata(
    data: xr.DataArray | xr.Dataset, complevel: int = 5
) -> xr.DataArray | xr.Dataset:
    """
    Export compressible netCDF files from xarray data (:py:class:`xarray.DataArray<xarray.DataArray>`, :py:class:`xarray.Dataset<xarray.Dataset>`)
    """
    comp = dict(zlib=True, complevel=complevel)
    if type(data) is xr.Dataset:
        for var in data:
            data[var].encoding.update(comp)
    if type(data) is xr.DataArray:
        data.encoding.update(comp)
    return data


def transfer_dFdp2dFdz(
    dFdp_data: xr.DataArray | xr.Dataset, rho_d: float = 1.2928e3, g: float = 9.8
):
    """

    The transformation relationship between the z coordinate system and the p coordinate system.

    .. math::
        \\frac{\\partial F}{\\partial z} = \\frac{\\partial F}{\\partial p} \\frac{\\partial p}{\\partial z} = - \\rho g \\frac{\\partial F}{\\partial p}
    """
    return -rho_d * g * dFdp_data


def sort_ascending_latlon_coordinates(
    data: xr.DataArray | xr.Dataset, lat_dim: str = "lat", lon_dim: str = "lon"
) -> xr.DataArray | xr.Dataset:
    """
    Sort the dimensions `lat`, `lon` in ascending order.
    """
    if (lat_dim is None) and (lon_dim is None):
        return data
    elif lat_dim is None:
        return data.sortby([lon_dim], ascending=True)
    elif lon_dim is None:
        return data.sortby([lat_dim], ascending=True)
    else:
        return data.sortby([lat_dim, lon_dim], ascending=True)


def transfer_units_coeff(input_units: str, output_units: str) -> float:
    """
    Unit conversion factor
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
    Data unit conversion for multiple transition.
    """
    base = transfer_units_coeff(input_units, output_units)
    output_data = input_data * base
    output_data.attrs["units"] = output_units
    return output_data


def transfer_data_difference_units(
    input_data: xr.DataArray | xr.Dataset, input_units: str, output_units: str
) -> xr.DataArray | xr.Dataset:
    """
    Data unit conversion for difference transition.
    """
    coeff = transfer_units_coeff(input_units, output_units)
    output_data = input_data + coeff - 1
    output_data.attrs["units"] = output_units
    return output_data


def transfer_data_temperature_units(
    input_data: xr.DataArray | xr.Dataset, input_units: str, output_units: str
) -> xr.DataArray | xr.Dataset:
    """
    Converts a value from one temperature unit to another. Aliases are supported.

    将值从一个温度单位转换为另一个温度单位，支持别名。
    """
    # 定义标准单位及其别名
    unit_aliases = {
        "K": ["kelvin", "K", "degK", "开氏度"],
        "degC": ["celsius", "degC", "摄氏度"],
        "degF": ["fahrenheit", "degF"],
        "degR": ["rankine", "degR", "兰氏度"],  # 兰氏度
        "degRe": ["reaumur", "degRe", "列氏度"],  # 列氏度
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
    Data unit conversion for ANY type transition.

    .. warning::

        NOT support ``dask``.
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


def generate_dataset_dispatcher(func):
    """
    Function Dispensers: Iterate over the variables in the `xarray.Dataset` data using a function that only supports `xarray.DataArray` data
    """

    def apply_to_dataset(data, func, *args, **kwargs):
        # Apply the function to each variable in `xarray.Dataset`.
        result = {}
        for var_name, data_array in data.data_vars.items():
            try:
                result[var_name] = func(data_array, *args, **kwargs)
            except Exception as e:
                warnings.warn(
                    f"Variable '{var_name}' cannot be processed with the provided function. Ignoring."
                )
        return xr.Dataset(result)

    def wrapper(data, *args, **kwargs):
        # Type determination
        if isinstance(data, xr.Dataset):
            return apply_to_dataset(data, func, *args, **kwargs)
        elif isinstance(data, xr.DataArray):
            return func(data, *args, **kwargs)
        else:
            raise ValueError("Unsupported input type. Expected DataArray or Dataset.")

    return wrapper


# def generate_datatree_dispatcher(func):
#     """
#     Function Dispensers: Iterate over the variables in the `xarray.Dataset` data using a function that only supports `xarray.DataArray` data
#     """

#     def apply_to_dataset(data, func, returns_type, *args, **kwargs):
#         # Apply the function to each variable in `xarray.Dataset`.
#         result = {}
#         for var_name, data_array in data.data_vars.items():
#             try:
#                 result[var_name] = func(data_array, *args, **kwargs)
#             except Exception as e:
#                 warnings.warn(
#                     f"Variable '{var_name}' cannot be processed with the provided function. Ignoring."
#                 )

#         if returns_type == "dataset_returns":
#             return sort_datatree_by_dataset_returns(result)
#         elif returns_type == "dataset_vars":
#             return sort_datatree_by_dataset_vars(result)
#         raise ValueError(
#             "Unsupported input type. Expected 'dataset_returns' or 'dataset_vars'."
#         )

#     def sort_datatree_by_dataset_returns(dataset_list):
#         """
#         The dictionary containing `xarray.Dataset` data is returned as a Datatree according to the returns.
#         """

#         def are_sets_equal(lists):
#             """
#             Check that multiple lists that do not depend on order are identical
#             """
#             if len(lists) < 2:
#                 return True

#             first_list = lists[0]
#             for i in range(1, len(lists)):
#                 if len(first_list) != len(lists[i]):
#                     return False

#                 if (set(first_list) == set(lists[i])) == False:
#                     return False
#                 else:
#                     return True

#         # dataset var list
#         key_var_list = list(dataset_list.keys())

#         # dataset returns list
#         key_type_list = []
#         for type_out in list(dataset_list.keys()):
#             tmp = list(dataset_list[type_out].keys())
#             key_type_list.append(tmp)

#         # Checking dimensional consistency
#         if are_sets_equal(key_type_list) == False:
#             raise ValueError(
#                 "Make sure that the Dataset for each key contains the exact same variables!"
#             )

#         dt = DataTree(name="root")
#         for key_type, _ in dataset_list[key_var_list[0]].items():
#             # Get data for each returns about all `xarray.Dataset` variables
#             tmp = xr.Dataset()
#             for key_var, _ in dataset_list.items():
#                 tmp[key_var] = dataset_list[key_var][key_type]

#             # Add each returns to Datatree
#             dt[key_type] = tmp

#         return dt

#     def sort_datatree_by_dataset_vars(dataset_list):
#         """
#         The dictionary containing `xarray.Dataset` data is returned as a Datatree according to the variable name.
#         """
#         # Create DataTree object
#         dt = DataTree(name="root")
#         # Add items to the object
#         for key, value in dataset_list.items():
#             dt[key] = DataTree(name=key, dataset=value)

#         return dt

#     def wrapper(data, returns_type="dataset_returns", *args, **kwargs):
#         """
#         Closure core functions
#         """
#         # Type determination
#         if isinstance(data, xr.Dataset):
#             return apply_to_dataset(data, func, returns_type, *args, **kwargs)
#         elif isinstance(data, xr.DataArray):
#             return func(data, *args, **kwargs)
#         else:
#             raise ValueError("Unsupported input type. Expected DataArray or Dataset.")

#     return wrapper


def generate_datanode_dispatcher(func):
    """
    Function Dispensers: Iterate over the variables in the `xarray.Dataset` data using a function that only supports `xarray.DataArray` data
    """

    def apply_to_dataset(data, func, returns_type, *args, **kwargs):
        # Apply the function to each variable in `xarray.Dataset`.
        result = {}
        for var_name, data_array in data.data_vars.items():
            try:
                result[var_name] = func(data_array, *args, **kwargs)
            except Exception as e:
                warnings.warn(
                    f"Variable '{var_name}' cannot be processed with the provided function. Ignoring."
                )

        if returns_type == "dataset_returns":
            return sort_datatree_by_dataset_returns(result)
        elif returns_type == "dataset_vars":
            return sort_datatree_by_dataset_vars(result)
        raise ValueError(
            "Unsupported input type. Expected 'dataset_returns' or 'dataset_vars'."
        )

    def sort_datatree_by_dataset_returns(dataset_list):
        """
        The dictionary containing `xarray.Dataset` data is returned as a Datatree according to the returns.
        """

        def are_sets_equal(lists):
            """
            Check that multiple lists that do not depend on order are identical
            """
            if len(lists) < 2:
                return True

            first_list = lists[0]
            for i in range(1, len(lists)):
                if len(first_list) != len(lists[i]):
                    return False

                if (set(first_list) == set(lists[i])) == False:
                    return False
                else:
                    return True

        # dataset var list
        key_var_list = list(dataset_list.keys())

        # dataset returns list
        key_type_list = []
        for type_out in list(dataset_list.keys()):
            tmp = list(dataset_list[type_out].keys())
            key_type_list.append(tmp)

        # Checking dimensional consistency
        if are_sets_equal(key_type_list) == False:
            raise ValueError(
                "Make sure that the Dataset for each key contains the exact same variables!"
            )

        dt = DataNode(name="root")
        for key_type, _ in dataset_list[key_var_list[0]].items():
            # Get data for each returns about all `xarray.Dataset` variables
            tmp = xr.Dataset()
            for key_var, _ in dataset_list.items():
                tmp[key_var] = dataset_list[key_var][key_type]

            # Add each returns to Datatree
            dt[key_type] = tmp

        return dt

    def sort_datatree_by_dataset_vars(dataset_list):
        """
        The dictionary containing `xarray.Dataset` data is returned as a Datatree according to the variable name.
        """
        # Create DataTree object
        dt = DataNode(name="root")
        # Add items to the object
        for key, value in dataset_list.items():
            dt[key] = value

        return dt

    def wrapper(data, returns_type="dataset_returns", *args, **kwargs):
        """
        Closure core functions
        """
        # Type determination
        if isinstance(data, xr.Dataset):
            return apply_to_dataset(data, func, returns_type, *args, **kwargs)
        elif isinstance(data, xr.DataArray):
            return func(data, *args, **kwargs)
        else:
            raise ValueError("Unsupported input type. Expected DataArray or Dataset.")

    return wrapper


@generate_dataset_dispatcher
def transfer_xarray_lon_from180TO360(
    data_input: xr.DataArray | xr.Dataset, lon_dim: str = "lon"
) -> xr.DataArray | xr.Dataset:
    """
    Longitude conversion -180-180 to 0-360.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The spatio-temporal data to be calculated.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.

    .. seealso::
        :py:func:`transfer_xarray_lon_from360TO180 <transfer_xarray_lon_from360TO180>`
    """

    def transfer_lon_180TO360(lon_array):
        return (lon_array + 360) % 360

    lon_array = data_input[lon_dim].data

    if (lon_array > 180).any():
        warnings.warn(
            "It seems that the input data longitude range is not from -180° to 180°. Please carefully check your data."
        )

    lon_array = transfer_lon_180TO360(lon_array)
    tmp = data_input.assign_coords({lon_dim: lon_array}).sortby(lon_dim)
    return tmp


@generate_dataset_dispatcher
def transfer_xarray_lon_from360TO180(
    data_input: xr.DataArray | xr.Dataset, lon_dim: str = "lon"
) -> xr.DataArray | xr.Dataset:
    """
    Longitude conversion 0-360 to -180-180.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
         The spatio-temporal data to be calculated.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.

    .. seealso::
        :py:func:`transfer_xarray_lon_from180TO360 <transfer_xarray_lon_from180TO360>`
    """

    def transfer_lon_360TO180(lon_array):
        return (lon_array + 180) % 360 - 180

    lon_array = data_input[lon_dim].data

    if (lon_array < 0).any():
        warnings.warn(
            "It seems that the input data longitude range is not from 0° to 360°. Please carefully check your data."
        )

    lon_array = transfer_lon_360TO180(lon_array)
    tmp = data_input.assign_coords({lon_dim: lon_array}).sortby(lon_dim)
    return tmp


def module_available(module: str) -> bool:
    """Checks whether a module is installed without importing it.

    Use this for a lightweight check and lazy imports.

    Parameters
    ----------
    module : dim: :py:class:`str <str>`
        Name of the module.

    Returns
    -------
    available : :py:class:`bool <bool>`
        Whether the module is installed.
    """
    import importlib

    return importlib.util.find_spec(module) is not None


def dequantify_metpy_xarraydata(data: xr.DataArray) -> xr.DataArray:
    """
    Return a new DataArray with the data as magnitude and the units as an attribute (Metpy).

    .. note::

        https://unidata.github.io/MetPy/latest/api/generated/metpy.xarray.html#metpy.xarray.MetPyDataArrayAccessor.dequantify
    """
    return data.metpy.dequantify()


def reverse_bool_xarraydata(
    data_input: xr.DataArray | xr.Dataset,
) -> xr.DataArray | xr.Dataset:
    """
    Reverse the bool type in the `data_input`. i.e., `True` -> `False`, and `False` -> `True`.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        Input dataset.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
    """
    return xr.where(data_input, False, True)


def compare_two_dataarray_coordinate(
    data_input1: xr.DataArray,
    data_input2: xr.DataArray,
    time_dim: str = "time",
    exclude_dims: list[str] = [],
):
    """
    Compare two DataArray data whether they have the same dimensions (without comparing internal data)

    Parameters
    ----------
    data_input1 : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        Input dataset 1.
    data_input2 : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        Input dataset 2.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.
    exclude_dims: :py:class:`list <list>`, default: `[]`.
        The exclude comparison dimensions.
    """
    a_dims = data_input1.dims
    b_dims = data_input2.dims

    # remove dimensions that do not participate in comparison
    a_dims = tuple(x for x in a_dims if not x in exclude_dims)
    b_dims = tuple(x for x in b_dims if not x in exclude_dims)

    for item_dims in a_dims:
        # dims
        if item_dims in b_dims:
            pass
        else:
            raise AssertionError(
                f"The coodinnate `{item_dims}` does not exist in two datasets simultaneously."
            )

        # dims shape
        a_item_dims_shape = data_input1[item_dims].shape[0]
        b_item_dims_shape = data_input2[item_dims].shape[0]
        if a_item_dims_shape == b_item_dims_shape:
            pass
        else:
            raise AssertionError(
                f"The shape of coodinnate `{item_dims}` (i.e., {a_item_dims_shape} vs {b_item_dims_shape}) is not same in two datasets."
            )

        # dims data
        a_item_dims_data = data_input1[item_dims].data
        b_item_dims_data = data_input2[item_dims].data
        if item_dims != time_dim:
            if np.allclose(a_item_dims_data, b_item_dims_data, equal_nan=True):
                pass
            else:
                raise AssertionError(
                    f"The data in the coodinnate `{item_dims}` is not same in two datasets. If it is the precision problem, you can specify the value of one data dimension as the value of another corresponding data, i.e., `a[{item_dims}] = b[{item_dims}]`."
                )
        else:
            if (data_input1[item_dims].data == data_input2[item_dims].data).all():
                pass
            else:
                warnings.warn(
                    f"The data in the coodinnate `{item_dims}` is not same in two datasets. If it is the precision problem, you can specify the value of one data dimension as the value of another corresponding data, i.e., `a[{item_dims}] = b[{item_dims}]`.",
                    RuntimeWarning,
                )


def compare_multi_dataarray_coordinate(
    data_input_list: list[xr.DataArray],
    time_dim: str = "time",
    exclude_dims: list[str] = [],
):
    """
    Compare multi-DataArray data whether they have the same dimensions (without comparing internal data)

    Parameters
    ----------
    data_input_list : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        Input dataset list.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.
    exclude_dims: :py:class:`list <list>`, default: `[]`.
        The exclude comparison dimensions.
    """
    item0 = data_input_list[0]
    for dataarray_item in data_input_list[1:]:
        compare_two_dataarray_coordinate(
            item0, dataarray_item, time_dim=time_dim, exclude_dims=exclude_dims
        )


def validate_dataarrays(
    dataarrays: Union[xr.DataArray, List[xr.DataArray], Tuple[xr.DataArray, ...]],
    dims: Optional[List[str]] = None,
    time_dims: Union[str, List[str], None] = "time",
) -> bool:
    """
    Validate consistency of multiple DataArrays across specified dimensions.

    Parameters
    -------------
    dataarrays : xarray.DataArray or list/tuple of xarray.DataArray
        DataArray(s) to be validated
    dims : list of str, optional
        List of dimension names to validate. If None, validates all dimensions.
    time_dims : str or list of str, optional
        Dimension names where only size is validated (coordinate values are not checked).
        Default is "time".

    Returns
    ----------
    bool
        Returns True if all validations pass, otherwise raises an exception.

    Raises
    ---------
    TypeError
        If input is not a DataArray or list/tuple of DataArrays
    ValueError
        If fewer than 2 DataArrays are provided

    Examples
    -----------
    >>> import xarray as xr
    >>> da1 = xr.DataArray(np.random.rand(10, 5), dims=['time', 'x'])
    >>> da2 = xr.DataArray(np.random.rand(10, 5), dims=['time', 'x'])
    >>> validate_dataarrays([da1, da2])  # Validates all dimensions
    True

    >>> da3 = xr.DataArray(np.random.rand(10, 6), dims=['time', 'x'])
    >>> validate_dataarrays([da1, da3])  # Raises ValueError for dimension mismatch
    """
    # Handle single DataArray input
    if isinstance(dataarrays, xr.DataArray):
        raise TypeError(
            "At least two DataArrays required for comparison, but got a single DataArray"
        )

    # Ensure dataarrays is a list or tuple
    if not isinstance(dataarrays, (list, tuple)):
        raise TypeError("dataarrays parameter must be a list or tuple of DataArrays")

    # Check number of DataArrays
    if len(dataarrays) < 2:
        raise ValueError("At least two DataArrays required for comparison")

    # Convert time_dims to list format
    if time_dims is None:
        time_dims = []
    elif isinstance(time_dims, str):
        time_dims = [time_dims]

    # Get first DataArray as reference
    ref_da = dataarrays[0]

    # If dims is None, validate all dimensions
    if dims is None:
        dims = list(ref_da.dims)

    # Check if each DataArray has the specified dimensions
    for da in dataarrays:
        missing_dims = set(dims) - set(da.dims)
        if missing_dims:
            raise ValueError(f"DataArray missing required dimensions: {missing_dims}")

    # Validate each dimension
    for dim in dims:
        # Check dimension sizes match
        ref_size = ref_da.sizes[dim]
        for i, da in enumerate(dataarrays[1:], 1):
            if da.sizes[dim] != ref_size:
                raise ValueError(
                    f"Dimension size mismatch: on dimension '{dim}', "
                    f"DataArray 0 has size {ref_size}, "
                    f"DataArray {i} has size {da.sizes[dim]}"
                )

        # Check coordinate values (skip for time_dims dimensions)
        if dim in time_dims:
            # Special handling for time dimensions
            ref_coords = ref_da[dim].values
            for i, da in enumerate(dataarrays[1:], 1):
                try:
                    if not np.array_equal(da[dim].values, ref_coords):
                        warnings.warn(
                            f"Time coordinate values not identical: on dimension '{dim}', "
                            f"DataArray 0 and DataArray {i} have different time coordinates"
                        )
                except AttributeError:
                    # Skip if DataArray doesn't have this dimension
                    pass
        else:
            # Strict check for non-time dimensions
            ref_coords = ref_da[dim].values
            for i, da in enumerate(dataarrays[1:], 1):
                if not np.array_equal(da[dim].values, ref_coords):
                    raise ValueError(
                        f"Coordinate value mismatch: on dimension '{dim}', "
                        f"DataArray 0 and DataArray {i} have different coordinates"
                    )

    return True


def check_deprecation_status(current_version, removal_version):
    """Check whether a deprecation warning should be turned into an error"""
    return version.parse(current_version) >= version.parse(removal_version)


def deprecated(version, removal_version, replacement=None):
    """
    Deprecate decorator

    :param version: Current version (deprecated version)
    :param removal_version: indicates the version to be removed
    :param replacement: Name of the replacement function (optional)
    """

    def decorator(old_func):
        @wraps(old_func)
        def wrapped(*args, **kwargs):
            if check_deprecation_status(__version__, removal_version):
                raise RuntimeError(
                    f"{old_func.__name__} was removed in version {removal_version}"
                )

            message = f"{old_func.__name__} is deprecated since {version} and will be removed in {removal_version}."
            if replacement:
                message += f" Use `{replacement}` instead."

            warnings.warn(message, DeprecationWarning, stacklevel=2)
            return old_func(*args, **kwargs)

        return wrapped

    return decorator


def datetime_to_numeric(datetime_array, unit="ns"):
    """
    Convert datetime64[ns] array to numeric values

    Args:
        datetime_array: numpy datetime64[ns] array
        unit: conversion unit, options: 'ns'(nanoseconds), 'us'(microseconds),
              'ms'(milliseconds), 's'(seconds), 'm'(minutes), 'h'(hours), 'D'(days)

    Returns:
        Numeric array (typically int64)
    """
    if unit not in ["ns", "us", "ms", "s", "m", "h", "D", "W", "M", "Y"]:
        raise ValueError(
            "Unit must be one of: 'ns', 'us', 'ms', 's', 'm', 'h', 'D', 'W', 'M', 'Y'"
        )

    # First convert to specified unit, then to numeric
    return datetime_array.astype(f"datetime64[{unit}]").astype("int64")


def numeric_to_datetime(numeric_array, unit="ns"):
    """
    Convert numeric array back to datetime64[ns]

    Args:
        numeric_array: numeric array
        unit: unit of the numeric values (must match conversion unit)

    Returns:
        datetime64[ns] array
    """
    if unit not in ["ns", "us", "ms", "s", "m", "h", "D", "W", "M", "Y"]:
        raise ValueError(
            "Unit must be one of: 'ns', 'us', 'ms', 's', 'm', 'h', 'D', 'W', 'M', 'Y'"
        )

    # Convert to timedelta, then add epoch time (1970-01-01)
    return (
        numeric_array.astype(f"timedelta64[{unit}]") + np.datetime64("1970-01-01")
    ).astype("datetime64[ns]")


def calculate_time_steps(datetime_array, unit="D"):
    """
    Calculate time steps between consecutive datetime values

    Args:
        datetime_array: datetime64[ns] array
        unit: unit for the returned time steps

    Returns:
        Array of time differences between consecutive points (numeric)
    """
    numeric_values = datetime_to_numeric(datetime_array, unit)
    return np.diff(numeric_values)


def clean_extra_coords(
    data_input: xr.DataArray,
) -> xr.DataArray:
    """
    Remove coordinate dimensions that exist in coordinates but not in dimensions.

    This function cleans up a DataArray by removing any coordinate variables that
    are not listed in the array's dimensions. This is useful for ensuring
    consistent data structures and avoiding potential issues with operations
    that expect dimension coordinates to match actual array dimensions.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        The input DataArray to be cleaned.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`
        The cleaned DataArray with only coordinates that match its dimensions.

    Examples
    --------
    >>> import xarray as xr
    >>> import numpy as np
    >>> da = xr.DataArray(np.random.rand(3, 4),
    ...                   dims=['x', 'y'],
    ...                   coords={'x': [1,2,3],
    ...                           'y': [1,2,3,4],
    ...                           'time': [1,2,3],
    ...                           'extra': ('z', [1.1, 2.2])})
    >>> cleaned = clean_extra_coords(da)
    >>> 'time' in cleaned.coords  # Returns False
    False
    >>> 'extra' in cleaned.coords  # Returns False
    False
    """
    # Find coordinates that aren't in dimensions
    extra_dims = [dim for dim in data_input.coords if dim not in data_input.dims]

    # Remove extra coordinates if any exist
    if extra_dims:
        return data_input.reset_coords(extra_dims, drop=True)
    return data_input
