easyclimate.field.air_sea_interaction.index_amm
===============================================

.. py:module:: easyclimate.field.air_sea_interaction.index_amm

.. autoapi-nested-parse::

   Atlantic Meridional Mode (AMM) Index



Functions
---------

.. autoapisummary::

   easyclimate.field.air_sea_interaction.index_amm.calc_index_AMM_Doi_2009


Module Contents
---------------

.. py:function:: calc_index_AMM_Doi_2009(sst_monthly_data: xarray.DataArray, time_range: slice = slice(None, None), lon_dim: str = 'lon', lat_dim: str = 'lat', time_dim: str = 'time', normalized: bool = False) -> xarray.DataArray

   The calculation of monthly mean Atlantic Meridional Mode (AMM) index is constructed by following method:

       The difference difference between the northern index (SSTA in 5–15°N, 50–20°W) and the southern index (SSTA in 5–15°S, 20°W–10°E).
       See Fig. 7 in Doi et al. 2009.

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
   The monthly mean AMM index (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - `Doi, T., Tozuka, T. & Yamagata, T. Interannual variability of the Guinea Dome and its possible link with the Atlantic Meridional Mode. Clim Dyn 33, 985–998 (2009). <https://doi.org/10.1007/s00382-009-0574-z>`__
   - `Doi, T., Tozuka, T., & Yamagata, T. (2010). The Atlantic Meridional Mode and Its Coupled Variability with the Guinea Dome. Journal of Climate, 23(2), 455-475. <https://doi.org/10.1175/2009JCLI3198.1>`__
   - `Lee, S.-K., Lopez, H., Tuchen, F. P., Kim, D., Foltz, G. R., & Wittenberg, A. T. (2023). On the genesis of the 2021 Atlantic Niño. Geophysical Research Letters, 50, e2023GL104452. <https://doi.org/10.1029/2023GL104452>`__


