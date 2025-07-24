easyclimate.field.air_sea_interaction.index_iod
===============================================

.. py:module:: easyclimate.field.air_sea_interaction.index_iod

.. autoapi-nested-parse::

   Indian Ocean Dipole (IOD) Index



Functions
---------

.. autoapisummary::

   easyclimate.field.air_sea_interaction.index_iod.calc_index_IOD_Saji_1999


Module Contents
---------------

.. py:function:: calc_index_IOD_Saji_1999(sst_monthly_data: xarray.DataArray, time_range: slice = slice(None, None), lon_dim: str = 'lon', lat_dim: str = 'lat', time_dim: str = 'time', normalized: bool = False) -> xarray.DataArray

   The calculation of monthly mean Indian Ocean Dipole (IOD) index (i.e., Dipole Mode Index; DMI) is constructed by following method:

       The difference in SST anomaly between the tropical western Indian Ocean (50°E-70°E, 10°S-10°N)
       and the tropical south-eastern Indian Ocean (90°E-110°E, 10°S-Equator).

   Parameters
   ----------
   sst_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The monthly sea surface temperature (SST) dataset.
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
   The monthly mean IOD index (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - `Saji, N., Goswami, B., Vinayachandran, P. et al. A dipole mode in the tropical Indian Ocean. Nature 401, 360–363 (1999). <https://doi.org/10.1038/43854>`__


