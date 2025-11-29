easyclimate.field.teleconnection.index_cgt
==========================================

.. py:module:: easyclimate.field.teleconnection.index_cgt

.. autoapi-nested-parse::

   Circumglobal teleconnection pattern (CGT) Index



Functions
---------

.. autoapisummary::

   easyclimate.field.teleconnection.index_cgt.calc_index_CGT_1point_Ding_Wang_2005
   easyclimate.field.teleconnection.index_cgt.calc_index_CGT_NH_Ding_Wang_2005


Module Contents
---------------

.. py:function:: calc_index_CGT_1point_Ding_Wang_2005(z200_monthly_data: xarray.DataArray, time_range: slice = slice(None, None), lon_dim: str = 'lon', lat_dim: str = 'lat', time_dim: str = 'time', normalized: bool = True) -> xarray.DataArray

   The calculation of monthly mean circumglobal teleconnection pattern (CGT) index is constructed by following method:

       Z200 anomalies averaged over the area (:math:`\mathrm{ 35 ^{\circ}N - 40 ^{\circ}N; 60^{\circ}-70^{\circ}E }`).

   Parameters
   ----------
   z200_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The monthly 200-hPa geopotential height.
   time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
       The time range of seasonal cycle means to be calculated. The default value is the entire time range.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name.
   normalized: :py:class:`bool <bool>`, default `True`, optional.
       Whether to standardize the index based on standard deviation over `time_range`.

   Returns
   -------
   The monthly mean CGT index (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - `Ding, Q., & Wang, B. (2005). Circumglobal Teleconnection in the Northern Hemisphere Summer. Journal of Climate, 18(17), 3483-3505. <https://doi.org/10.1175/JCLI3473.1>`__
   - `Liu Y. Relationship between the Silk Road and Circumglobal Teleconnection Patterns on the Interannual and Interdecadal Timescales. Atmosphere. 2023; 14(11):1626. <https://doi.org/10.3390/atmos14111626>`__



.. py:function:: calc_index_CGT_NH_Ding_Wang_2005(z200_monthly_data: xarray.DataArray, output_freq: Literal['monthly', 'seasonally'], time_range: slice = slice(None, None), lon_dim: str = 'lon', lat_dim: str = 'lat', lat_range: slice = slice(20, 85), time_dim: str = 'time') -> xarray.DataArray

   The calculation of monthly mean circumglobal teleconnection pattern (CGT) index using Ding, Q., & Wang, B. (2005) method

   Parameters
   ----------
   z200_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The monthly 200-hPa geopotential height.
   output_freq: :py:class:`str <str>`.
       The output frequency. Optional values are `monthly`, `seasonally`.
   time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
       The time range of seasonal cycle means to be calculated. The default value is the entire time range.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lat_range: :py:class:`slice <slice>`, default: `slice(20, 85)`.
       The latitude range of computation using EOFs over the Northern Hemisphere. The default value is from :math:`\mathrm{20^{\circ}N}` to :math:`\mathrm{85^{\circ}N}`.
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name.

   Returns
   -------
   The monthly/seasonally mean CGT index (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - `Ding, Q., & Wang, B. (2005). Circumglobal Teleconnection in the Northern Hemisphere Summer. Journal of Climate, 18(17), 3483-3505. <https://doi.org/10.1175/JCLI3473.1>`__
   - `Liu Y. Relationship between the Silk Road and Circumglobal Teleconnection Patterns on the Interannual and Interdecadal Timescales. Atmosphere. 2023; 14(11):1626. <https://doi.org/10.3390/atmos14111626>`__

   .. seealso::
       :py:func:`get_EOF_model <easyclimate.core.eof.get_EOF_model>`


