easyclimate.field.air_sea_interaction.index_enso
================================================

.. py:module:: easyclimate.field.air_sea_interaction.index_enso

.. autoapi-nested-parse::

   ENSO Index



Functions
---------

.. autoapisummary::

   easyclimate.field.air_sea_interaction.index_enso.calc_index_nino1and2
   easyclimate.field.air_sea_interaction.index_enso.calc_index_nino3
   easyclimate.field.air_sea_interaction.index_enso.calc_index_nino34
   easyclimate.field.air_sea_interaction.index_enso.calc_index_OMI
   easyclimate.field.air_sea_interaction.index_enso.calc_index_nino4


Module Contents
---------------

.. py:function:: calc_index_nino1and2(sst_monthly_data: xarray.DataArray | xarray.Dataset, time_range: slice = slice(None, None), lat_dim: str = 'lat', lon_dim: str = 'lon', time_dim: str = 'time', normalized: bool = False) -> xarray.DataArray | xarray.Dataset

   Calculate the Niño 1+2 index.

   The Niño 1+2 region (0°-10°S, 90°W-80°W) is the smallest and eastern-most of the Niño SST regions,
   and corresponds with the region of coastal South America where El Niño was
   first recognized by the local populations. This index tends to have the largest variance of the Niño SST indices.

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
   Niño 1+2 index.

   Reference
   --------------
   - Trenberth, Kevin & National Center for Atmospheric Research Staff (Eds). Last modified 2023-07-25 "The Climate Data Guide: Nino SST Indices (Nino 1+2, 3, 3.4, 4; ONI and TNI)." Retrieved from https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni on 2023-11-12.
   - El Niño Index Dashboard. Website: https://psl.noaa.gov/enso/dashboard.html


.. py:function:: calc_index_nino3(sst_monthly_data: xarray.DataArray | xarray.Dataset, time_range: slice = slice(None, None), lat_dim: str = 'lat', lon_dim: str = 'lon', time_dim: str = 'time', normalized: bool = False) -> xarray.DataArray | xarray.Dataset

   Calculate the Niño 3 index.

   This region (5°N-5°S, 150°W-90°W) was once the primary focus for monitoring and predicting El Niño,
   but researchers later learned that the key region for coupled ocean-atmosphere interactions for ENSO lies further west (Trenberth, 1997).
   Hence, the Niño 3.4 and ONI became favored for defining El Niño and La Niña events.

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
   Niño 3 index.

   Reference
   --------------
   - Trenberth, Kevin & National Center for Atmospheric Research Staff (Eds). Last modified 2023-07-25 "The Climate Data Guide: Nino SST Indices (Nino 1+2, 3, 3.4, 4; ONI and TNI)." Retrieved from https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni on 2023-11-12.
   - El Niño Index Dashboard. Website: https://psl.noaa.gov/enso/dashboard.html
   - Trenberth, K. E., 1997: The Definition of El Niño. Bull. Amer. Meteor. Soc., 78, 2771–2778, https://doi.org/10.1175/1520-0477(1997)078<2771:TDOENO>2.0.CO;2.


.. py:function:: calc_index_nino34(sst_monthly_data: xarray.DataArray | xarray.Dataset, time_range: slice = slice(None, None), lat_dim: str = 'lat', lon_dim: str = 'lon', running_mean=5, time_dim: str = 'time', normalized: bool = False) -> xarray.DataArray | xarray.Dataset

   Calculate the Niño 3.4 index.

   The  Niño 3.4 (5°N-5°S, 170°W-120°W) anomalies may be thought of as representing the average equatorial SSTs
   across the Pacific from about the dateline to the South American coast.
   The Niño 3.4 index typically uses a 5-month running mean, and El Niño or
   La  Niña events are defined when the  Niño 3.4 SSTs exceed +/- 0.4℃ for a period of six months or more.

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
   running_mean: :py:class:`int <int>`, default: `5`.
       Running mean value. If `running_mean` is `None` or `0`, it will not perform running average operation.
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name.
   normalized: :py:class:`bool <bool>`, default `False`, optional.
       Whether to standardize the index based on standard deviation over `time_range`.

   Returns
   -------
   Niño 3.4 index.

   Reference
   --------------
   - Trenberth, Kevin & National Center for Atmospheric Research Staff (Eds). Last modified 2023-07-25 "The Climate Data Guide: Nino SST Indices (Nino 1+2, 3, 3.4, 4; ONI and TNI)." Retrieved from https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni on 2023-11-12.
   - El Niño Index Dashboard. Website: https://psl.noaa.gov/enso/dashboard.html


.. py:function:: calc_index_OMI(sst_monthly_data: xarray.DataArray | xarray.Dataset, time_range: slice = slice(None, None), lat_dim: str = 'lat', lon_dim: str = 'lon', running_mean=3, time_dim: str = 'time', normalized: bool = False) -> xarray.DataArray | xarray.Dataset

   Calculate the ONI (Oceanic Niño Index) index.

   The ONI (5°N-5°S, 170°W-120°W) uses the same region as the Niño 3.4 index.
   The ONI uses a 3-month running mean, and to be classified as a full-fledged El Niño or La Niña,
   the anomalies must exceed +0.5℃ or -0.5℃ for at least five consecutive months.
   This is the operational definition used by NOAA.

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
   running_mean: :py:class:`int <int>`, default: `3`.
       Running mean value. If `running_mean` is `None` or `0`, it will not perform running average operation.
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name.
   normalized: :py:class:`bool <bool>`, default `False`, optional.
       Whether to standardize the index based on standard deviation over `time_range`.

   Returns
   -------
   ONI index.

   Reference
   --------------
   - Trenberth, Kevin & National Center for Atmospheric Research Staff (Eds). Last modified 2023-07-25 "The Climate Data Guide: Nino SST Indices (Nino 1+2, 3, 3.4, 4; ONI and TNI)." Retrieved from https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni on 2023-11-12.
   - El Niño Index Dashboard. Website: https://psl.noaa.gov/enso/dashboard.html


.. py:function:: calc_index_nino4(sst_monthly_data: xarray.DataArray | xarray.Dataset, time_range: slice = slice(None, None), lat_dim: str = 'lat', lon_dim: str = 'lon', time_dim: str = 'time', normalized: bool = False) -> xarray.DataArray | xarray.Dataset

   Calculate the Niño 4 index.

   The Niño 4 index (5°N-5°S, 160°E-150°W) captures SST anomalies in the central equatorial Pacific.
   This region tends to have less variance than the other Niño regions.

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
   Niño 4 index.

   Reference
   --------------
   - Trenberth, Kevin & National Center for Atmospheric Research Staff (Eds). Last modified 2023-07-25 "The Climate Data Guide: Nino SST Indices (Nino 1+2, 3, 3.4, 4; ONI and TNI)." Retrieved from https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni on 2023-11-12.
   - El Niño Index Dashboard. Website: https://psl.noaa.gov/enso/dashboard.html


