easyclimate.core.units
======================

.. py:module:: easyclimate.core.units

.. autoapi-nested-parse::

   Functions for handling calculations of units



Functions
---------

.. autoapisummary::

   easyclimate.core.units.transfer_units_coeff
   easyclimate.core.units.transfer_data_multiple_units
   easyclimate.core.units.transfer_data_difference_units
   easyclimate.core.units.transfer_data_temperature_units
   easyclimate.core.units.transfer_data_units


Module Contents
---------------

.. py:function:: transfer_units_coeff(input_units: str, output_units: str) -> float

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


.. py:function:: transfer_data_multiple_units(input_data: xarray.DataArray | xarray.Dataset, input_units: str, output_units: str) -> xarray.DataArray | xarray.Dataset

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


.. py:function:: transfer_data_difference_units(input_data: xarray.DataArray | xarray.Dataset, input_units: str, output_units: str) -> xarray.DataArray | xarray.Dataset

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


.. py:function:: transfer_data_temperature_units(input_data: xarray.DataArray | xarray.Dataset, input_units: str, output_units: str) -> xarray.DataArray | xarray.Dataset

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


.. py:function:: transfer_data_units(input_data: xarray.DataArray | xarray.Dataset, input_units: str, output_units: str) -> xarray.DataArray | xarray.Dataset

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


