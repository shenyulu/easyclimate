easyclimate.core.utility
========================

.. py:module:: easyclimate.core.utility

.. autoapi-nested-parse::

   Functions for package utility.



Functions
---------

.. autoapisummary::

   easyclimate.core.utility.assert_compared_version
   easyclimate.core.utility.find_dims_axis
   easyclimate.core.utility.transfer_int2datetime
   easyclimate.core.utility.split_datetime2yearday
   easyclimate.core.utility.transfer_deg2rad
   easyclimate.core.utility.transfer_inf2nan
   easyclimate.core.utility.transfer_nan2value
   easyclimate.core.utility.get_weighted_spatial_data
   easyclimate.core.utility.get_compress_xarraydata
   easyclimate.core.utility.transfer_dFdp2dFdz
   easyclimate.core.utility.sort_ascending_latlon_coordinates
   easyclimate.core.utility.transfer_units_coeff
   easyclimate.core.utility.transfer_data_multiple_units
   easyclimate.core.utility.transfer_data_difference_units
   easyclimate.core.utility.generate_dataset_dispatcher
   easyclimate.core.utility.generate_datatree_dispatcher
   easyclimate.core.utility.transfer_xarray_lon_from180TO360
   easyclimate.core.utility.transfer_xarray_lon_from360TO180
   easyclimate.core.utility.module_available
   easyclimate.core.utility.dequantify_metpy_xarraydata
   easyclimate.core.utility.reverse_bool_xarraydata
   easyclimate.core.utility.compare_two_dataarray_coordinate
   easyclimate.core.utility.compare_multi_dataarray_coordinate


Module Contents
---------------

.. py:function:: assert_compared_version(ver1: float, ver2: float) -> int

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


.. py:function:: find_dims_axis(data: xarray.DataArray, dim: str) -> int

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


.. py:function:: transfer_int2datetime(data: numpy.array) -> numpy.datetime64

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


.. py:function:: split_datetime2yearday(ds: xarray.DataArray) -> xarray.DataArray

   Convert `np.datetime64` type with years and days to `year` and `day` coordinates.

   Parameters
   ----------
   - data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

   .. seealso::
       `Function in xarray to regroup monthly data into months and # of years <https://github.com/pydata/xarray/discussions/5119>`__.


.. py:function:: transfer_deg2rad(ds: xarray.DataArray) -> xarray.DataArray

   Convert Degrees to Radians.

   Parameters
   ----------
   - ds: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Degrees data.

   Returns
   -------
   - Radians data.: :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: transfer_inf2nan(ds: xarray.DataArray) -> xarray.DataArray

   Convert `np.inf` in `ds` to `np.nan`, respectively.

   Parameters
   ----------
   - ds: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Data include `np.inf`.

   Returns
   -------
   - Data include `np.nan`.: :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: transfer_nan2value(ds: xarray.DataArray, value: float) -> xarray.DataArray

   Convert `np.inf` in `ds` to `np.nan`, respectively.

   Parameters
   ----------
   - ds: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Data include `np.inf`.

   Returns
   -------
   - Data include `np.nan`.: :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: get_weighted_spatial_data(data_input: xarray.DataArray, lat_dim: str = 'lat', lon_dim: str = 'lon', method: str = 'cos_lat') -> xarray.DataArray

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


.. py:function:: get_compress_xarraydata(data: xarray.DataArray | xarray.Dataset, complevel: int) -> xarray.DataArray | xarray.Dataset

   Export compressible netCDF files from xarray data (:py:class:`xarray.DataArray<xarray.DataArray>`, :py:class:`xarray.Dataset<xarray.Dataset>`)


.. py:function:: transfer_dFdp2dFdz(dFdp_data: xarray.DataArray | xarray.Dataset, rho_d: float = 1292.8, g: float = 9.8)

   The transformation relationship between the z coordinate system and the p coordinate system.

   .. math::
       \frac{\partial F}{\partial z} = \frac{\partial F}{\partial p} \frac{\partial p}{\partial z} = - \rho g \frac{\partial F}{\partial p}


.. py:function:: sort_ascending_latlon_coordinates(data: xarray.DataArray | xarray.Dataset, lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray | xarray.Dataset

   Sort the dimensions `lat`, `lon` in ascending order.


.. py:function:: transfer_units_coeff(input_units: str, output_units: str) -> float

   Unit conversion factor


.. py:function:: transfer_data_multiple_units(input_data: xarray.DataArray | xarray.Dataset, input_units: str, output_units: str) -> xarray.DataArray | xarray.Dataset

   Data unit conversion for multiple transition.


.. py:function:: transfer_data_difference_units(input_data: xarray.DataArray | xarray.Dataset, input_units: str, output_units: str) -> xarray.DataArray | xarray.Dataset

   Data unit conversion for difference transition.


.. py:function:: generate_dataset_dispatcher(func)

   Function Dispensers: Iterate over the variables in the `xarray.Dataset` data using a function that only supports `xarray.DataArray` data


.. py:function:: generate_datatree_dispatcher(func)

   Function Dispensers: Iterate over the variables in the `xarray.Dataset` data using a function that only supports `xarray.DataArray` data


.. py:function:: transfer_xarray_lon_from180TO360(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon') -> xarray.DataArray | xarray.Dataset

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


.. py:function:: transfer_xarray_lon_from360TO180(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon') -> xarray.DataArray | xarray.Dataset

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


.. py:function:: module_available(module: str) -> bool

   Checks whether a module is installed without importing it.

   Use this for a lightweight check and lazy imports.

   Parameters
   ----------
   module : dim: :py:class:`str <str>`
       Name of the module.

   Returns
   -------
   available : :py:class:`bool <bool>`
       Whether the module is installed.


.. py:function:: dequantify_metpy_xarraydata(data: xarray.DataArray) -> xarray.DataArray

   Return a new DataArray with the data as magnitude and the units as an attribute (Metpy).

   .. note::

       https://unidata.github.io/MetPy/latest/api/generated/metpy.xarray.html#metpy.xarray.MetPyDataArrayAccessor.dequantify


.. py:function:: reverse_bool_xarraydata(data_input: xarray.DataArray | xarray.Dataset) -> xarray.DataArray | xarray.Dataset

   Reverse the bool type in the `data_input`. i.e., `True` -> `False`, and `False` -> `True`.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       Input dataset.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.


.. py:function:: compare_two_dataarray_coordinate(data_input1: xarray.DataArray, data_input2: xarray.DataArray, time_dim: str = 'time', exclude_dims: list[str] = [])

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


.. py:function:: compare_multi_dataarray_coordinate(data_input_list: list[xarray.DataArray], time_dim: str = 'time', exclude_dims: list[str] = [])

   Compare multi-DataArray data whether they have the same dimensions (without comparing internal data)

   Parameters
   ----------
   data_input_list : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       Input dataset list.
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name.
   exclude_dims: :py:class:`list <list>`, default: `[]`.
       The exclude comparison dimensions.


