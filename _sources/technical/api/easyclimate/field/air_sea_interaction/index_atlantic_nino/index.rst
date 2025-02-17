easyclimate.field.air_sea_interaction.index_atlantic_nino
=========================================================

.. py:module:: easyclimate.field.air_sea_interaction.index_atlantic_nino

.. autoapi-nested-parse::

   Atlantic Niños Index



Functions
---------

.. autoapisummary::

   easyclimate.field.air_sea_interaction.index_atlantic_nino.calc_index_ATL3


Module Contents
---------------

.. py:function:: calc_index_ATL3(sst_monthly_data: xarray.DataArray | xarray.Dataset, time_range: slice = slice(None, None), lat_dim: str = 'lat', lon_dim: str = 'lon', time_dim: str = 'time', normalized: bool = False) -> xarray.DataArray | xarray.Dataset

   Calculate ATL3 index.

   In some years, the cold tongue formation in summer is weak, leading to warm SST anomalies in the cold tongue region.
   This is what is called an Atlantic Niño. One way to gauge the strength of these events is to calculate the
   area average of SST in the cold tongue region, defined as 20°W to 0° and 3°S to 3°N.
   This is called the ATL3 index.

   Parameters
   ----------
   sst_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The monthly sea surface temperature (SST) dataset.
   time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
       The time range of seasonal cycle means to be calculated. The default value is the entire time range.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name.
   normalized: :py:class:`bool <bool>`, default `False`, optional.
       Whether to standardize the index based on standard deviation over `time_range`.

   Returns
   -------
   ATL3 index.

   Reference
   --------------
   - `Lee, S.-K., Lopez, H., Tuchen, F. P., Kim, D., Foltz, G. R., & Wittenberg, A. T. (2023). On the genesis of the 2021 Atlantic Niño. Geophysical Research Letters, 50, e2023GL104452. <https://doi.org/10.1029/2023GL104452>`__
   - Atlantic Niños. Website: https://www.jamstec.go.jp/aplinfo/climate/?page_id=1566


