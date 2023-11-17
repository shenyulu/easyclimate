:py:mod:`easyclimate.index`
===========================

.. py:module:: easyclimate.index


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   enso/index.rst
   npwi/index.rst
   oceanic_front/index.rst
   pna/index.rst


Package Contents
----------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.index.calc_index_NPWI
   easyclimate.index.find_PW_monsoon_region
   easyclimate.index.cal_NPWI_monsoon_onset
   easyclimate.index.cal_NPWI_monsoon_detreat
   easyclimate.index.get_weighted_spatial_data
   easyclimate.index.sort_ascending_latlon_coordinates
   easyclimate.index.calc_intensity_STFZ
   easyclimate.index.calc_intensity_SAFZ
   easyclimate.index.calc_location_STFZ
   easyclimate.index.calc_location_SAFZ
   easyclimate.index.calc_location_line_STFZ
   easyclimate.index.calc_location_line_SAFZ



.. py:function:: calc_index_NPWI(precipitable_water_daily: xarray.DataArray, time_dim='time') -> xarray.DataArray

   Calculate the normalized precipitable water index (NPWI).

   .. math::
       \mathrm{NPWI} = \frac{\mathrm{PW} - \mathrm{PW_{min}}}{\mathrm{PW_{max}} - \mathrm{PW_{min}}}

   Parameters
   ----------
   precipitable_water_daily: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Daily precipitable water data. 
   time_dim: :py:class:`str<python.str>`.
       The time coordinate dimension name.

   Returns
   -------
   Normalized precipitable water index (NPWI).

   Reference
   --------------
   Zeng, X., and E. Lu, 2004: Globally Unified Monsoon Onset and Retreat Indexes. J. Climate, 17, 2241–2248, https://doi.org/10.1175/1520-0442(2004)017<2241:GUMOAR>2.0.CO;2.


.. py:function:: find_PW_monsoon_region(precipitable_water_daily: xarray.DataArray, time_dim='time') -> xarray.DataArray

   The refined monsoon regions.

   .. note::
       To refine the definition of monsoon regions on a grid- cell-by-cell basis, 
       we first compute the 10-yr-averaged monthly PW over each cell. 
       Then we obtain the maximum monthly PW during the three summer months [e.g., June–August in the Northern Hemisphere (NH), denoted as PWw], 
       and the maximum monthly PW during the three winter months (e.g., December–February for the NH, denoted as PWc). 
       The refined monsoon regions are simply defined as grid cells that are within the monsoon 
       regions given in the above studies and have a difference between PWw and PWc greater than 12 mm. 
       Initially we have also tried to use the annual maximum and minimum monthly PW values.

   Parameters
   ----------
   precipitable_water_daily: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Daily precipitable water data.
   time_dim: :py:class:`str<python.str>`.
       The time coordinate dimension name.

   Reference
   --------------
   Zeng, X., and E. Lu, 2004: Globally Unified Monsoon Onset and Retreat Indexes. J. Climate, 17, 2241–2248, https://doi.org/10.1175/1520-0442(2004)017<2241:GUMOAR>2.0.CO;2.


.. py:function:: cal_NPWI_monsoon_onset(NPWI, thresh=0.618, consecutive_days=3, n=7, lon_dim='lon', lat_dim='lat', time_dim='time') -> xarray.DataArray

   Calculate the summer monsoon onset date.

   The summer monsoon onset date for grid cell G is defined as the first day (:math:`d`) 
   when NWPI is greater than the Golden Ratio (0.618) for three consecutive days
   in seven of the nine cells centered at cell G in day :math:`d` or (:math:`d \pm 1`).

   .. note::
       If one or more of the nine grids are undefined, for example, at the edge of monsoon regions, 
       the required number of seven is correspondingly reduced. 
       For instance, if only seven grid cells are defined, the required number is five.

   Parameters
   ----------
   NPWI: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Normalized precipitable water index (NPWI). 

       .. attention::
           It must include three dimensions: `time`, `longitude`, and `latitude`.

   thresh: :py:class:`float<python.float>`, default: `0.618`.
       Golden Ratio value for the threshold value.
   consecutive_days: :py:class:`int<python.int>`, default: `3`.
       Consecutive days values.
   n: :py:class:`int<python.int>`, default: `7`.
       :math:`n` of the nine cells centered at cell G in day :math:`d` or (:math:`d \pm 1`).
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   time_dim: :py:class:`str<python.str>`.
       The time coordinate dimension name.

   Returns
   -------
   Summer monsoon onset date.

   Reference
   --------------
   Zeng, X., and E. Lu, 2004: Globally Unified Monsoon Onset and Retreat Indexes. J. Climate, 17, 2241–2248, https://doi.org/10.1175/1520-0442(2004)017<2241:GUMOAR>2.0.CO;2.


.. py:function:: cal_NPWI_monsoon_detreat(NPWI, monsoon_onset_date, thresh=0.618, consecutive_days=3, n=7, lon_dim='lon', lat_dim='lat', time_dim='time') -> xarray.DataArray

   Calculate the summer monsoon retreat date.

   The summer monsoon retreat date for grid cell G is defined as the first day (:math:`d`) 
   when NWPI is less than the Golden Ratio (0.618) for three consecutive days
   in seven of the nine cells centered at cell G in day :math:`d` or (:math:`d \pm 1`).

   .. note::
       If one or more of the nine grids are undefined, for example, at the edge of monsoon regions, 
       the required number of seven is correspondingly reduced.
       For instance, if only seven grid cells are defined, the required number is five.

   Parameters
   ----------
   NPWI: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Normalized precipitable water index (NPWI). 

       .. attention::
           It must include three dimensions: `time`, `longitude`, and `latitude`.

   monsoon_onset_date: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Summer monsoon onset date. The results is generated by :py:func:`easyclimate.index.cal_NPWI_monsoon_onset <easyclimate.index.cal_NPWI_monsoon_onset>`.
   thresh: :py:class:`float<python.float>`, default: `0.618`.
       Golden Ratio value for the threshold value.
   consecutive_days: :py:class:`int<python.int>`, default: `3`.
       Consecutive days values.
   n: :py:class:`int<python.int>`, default: `7`.
       :math:`n` of the nine cells centered at cell G in day :math:`d` or (:math:`d \pm 1`).
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   time_dim: :py:class:`str<python.str>`.
       The time coordinate dimension name.

   Returns
   -------
   Summer monsoon retreat date.

   Reference
   --------------
   Zeng, X., and E. Lu, 2004: Globally Unified Monsoon Onset and Retreat Indexes. J. Climate, 17, 2241–2248, https://doi.org/10.1175/1520-0442(2004)017<2241:GUMOAR>2.0.CO;2.


.. py:function:: get_weighted_spatial_data(data_input: xarray.DataArray, lat_dim: str = 'lat', lon_dim: str = 'lon', method: str = 'cos_lat') -> xarray.DataArray

   Get the area-weighting data.

   Parameters
   ----------
   - data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   - lat_dim: :py:class:`str<python.str>`.
       Latitude dimension over which to apply. By default is applied over the `lat` dimension.
   - lon_dim: :py:class:`str<python.str>`.
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


.. py:function:: sort_ascending_latlon_coordinates(data: xr.DataArray | xr.Dataset, lat_dim: str = 'lat', lon_dim: str = 'lon') -> xr.DataArray | xr.Dataset

   Sort the dimensions `lat`, `lon` in ascending order.


.. py:function:: calc_intensity_STFZ(data_sst_DtDy: xarray.DataArray, criteria=0.45 * 1e-05, lat_range=[24, 32], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray

   Calculate the intensity of the subtropical frontal zone (STFZ). 
   The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

   .. math::
       \mathrm{ITS} = \sum_{i=1}^{N} \frac{G_i}{N}

   where :math:`G_i` is the value of zonally-averaged SST meridional gradient that is no less than 
   an empirically-given critical value (here, :math:`0.45 \times 10^{-5} \mathrm{km^{-1}}` for STFZ) at the :math:`i`-th latitudinal grid point within the zone, 
   and :math:`N` is the number of total grid points that satisfy the criteria above.

   Parameters
   ----------
   data_sst_DtDy : :py:class:`xarray.DataArray<xarray.DataArray>` 
       The SST meridional gradient data.
   criteria: :py:class:`float<python.float>`, default: `0.45 *1e-5`.
       Empirically-given critical value.
   lat_range: :py:class:`list<python.list>`, default: `[24, 32]`.
       The latitude range of the oceanic frontal zone.
   lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
       The longitude range of the oceanic frontal zone.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   The intensity of the STFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017), 
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts 
   and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766


.. py:function:: calc_intensity_SAFZ(data_sst_DtDy: xarray.DataArray, criteria=0.8 * 1e-05, lat_range=[36, 44], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray

   Calculate the intensity of the subarctic frontal zone (SAFZ). 
   The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

   .. math::
       \mathrm{ITS} = \sum_{i=1}^{N} \frac{G_i}{N}

   where :math:`G_i` is the value of zonally-averaged SST meridional gradient that is no less than 
   an empirically-given critical value (here, :math:`0.80 \times 10^{-5} \mathrm{km^{-1}}` for SAFZ) at the :math:`i`-th latitudinal grid point within the zone, 
   and :math:`N` is the number of total grid points that satisfy the criteria above.

   Parameters
   ----------
   data_sst_DtDy : :py:class:`xarray.DataArray<xarray.DataArray>` 
       The SST meridional gradient data.
   criteria: :py:class:`float<python.float>`, default: `0.80 *1e-5`.
       Empirically-given critical value.
   lat_range: :py:class:`list<python.list>`, default: `[36, 44]`.
       The latitude range of the oceanic frontal zone.
   lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
       The longitude range of the oceanic frontal zone.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   The intensity of the SAFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017), 
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts 
   and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766


.. py:function:: calc_location_STFZ(data_sst_DtDy: xarray.DataArray, criteria=0.45 * 1e-05, lat_range=[24, 32], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray

   Calculate the location index of the subtropical frontal zone (STFZ). 
   The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

   .. math::
       \mathrm{LCT} = \sum_{i=1}^{N} (G_i \times \mathrm{LAT}_i) / \sum_{i=1}^{N} G_i

   where :math:`\mathrm{LAT}_i` is the latitude at the :math:`i`-th grid point within the front zone. 
   Obviously, this definition reflects a weighted-average of :math:`\mathrm{LAT}_i` with respect to :math:`G_i`, 
   indicating that the location of a front is mainly determined by larger SST meridional gradients within the frontal zone.

   Parameters
   ----------
   data_sst_DtDy : :py:class:`xarray.DataArray<xarray.DataArray>` 
       The SST meridional gradient data.
   criteria: :py:class:`float<python.float>`, default: `0.45 *1e-5`.
       Empirically-given critical value.
   lat_range: :py:class:`list<python.list>`, default: `[24, 32]`.
       The latitude range of the oceanic frontal zone.
   lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
       The longitude range of the oceanic frontal zone.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   The location index of the STFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017), 
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts 
   and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766


.. py:function:: calc_location_SAFZ(data_sst_DtDy: xarray.DataArray, criteria=0.8 * 1e-05, lat_range=[36, 44], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray

   Calculate the location index of the subarctic frontal zone (SAFZ).
   The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

   .. math::
       \mathrm{LCT} = \sum_{i=1}^{N} (G_i \times \mathrm{LAT}_i) / \sum_{i=1}^{N} G_i

   where :math:`\mathrm{LAT}_i` is the latitude at the :math:`i`-th grid point within the front zone. 
   Obviously, this definition reflects a weighted-average of :math:`\mathrm{LAT}_i` with respect to :math:`G_i`, 
   indicating that the location of a front is mainly determined by larger SST meridional gradients within the frontal zone.

   Parameters
   ----------
   data_sst_DtDy : :py:class:`xarray.DataArray<xarray.DataArray>` 
       The SST meridional gradient data.
   criteria: :py:class:`float<python.float>`, default: `0.80 *1e-5`.
       Empirically-given critical value.
   lat_range: :py:class:`list<python.list>`, default: `[36, 44]`.
       The latitude range of the oceanic frontal zone.
   lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
       The longitude range of the oceanic frontal zone.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   The location index of the SAFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017), 
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts 
   and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766


.. py:function:: calc_location_line_STFZ(data_sst_DtDy: xarray.DataArray, criteria=0.45 * 1e-05, lat_range=[24, 32], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray

   Calculate the location of the subtropical frontal zone (STFZ). 
   The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

   .. math::
       \mathrm{LCT} = (\sum_{i=1}^{N} (G_i \times \mathrm{LAT}_i)) / G_i

   where :math:`\mathrm{LAT}_i` is the latitude at the :math:`i`-th grid point within the front zone. 
   Obviously, this definition reflects a weighted-average of :math:`\mathrm{LAT}_i` with respect to :math:`G_i`, 
   indicating that the location of a front is mainly determined by larger SST meridional gradients within the frontal zone.

   Parameters
   ----------
   data_sst_DtDy : :py:class:`xarray.DataArray<xarray.DataArray>` 
       The SST meridional gradient data.
   criteria: :py:class:`float<python.float>`, default: `0.45 *1e-5`.
       Empirically-given critical value.
   lat_range: :py:class:`list<python.list>`, default: `[24, 32]`.
       The latitude range of the oceanic frontal zone.
   lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
       The longitude range of the oceanic frontal zone.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   The location of the STFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017), 
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts 
   and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766    


.. py:function:: calc_location_line_SAFZ(data_sst_DtDy: xarray.DataArray, criteria=0.8 * 1e-05, lat_range=[36, 44], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon') -> xarray.DataArray

   Calculate the location of the subarctic frontal zone (SAFZ).
   The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

   .. math::
       \mathrm{LCT} = (\sum_{i=1}^{N} (G_i \times \mathrm{LAT}_i)) / G_i

   where :math:`\mathrm{LAT}_i` is the latitude at the :math:`i`-th grid point within the front zone. 
   Obviously, this definition reflects a weighted-average of :math:`\mathrm{LAT}_i` with respect to :math:`G_i`, 
   indicating that the location of a front is mainly determined by larger SST meridional gradients within the frontal zone.

   Parameters
   ----------
   data_sst_DtDy : :py:class:`xarray.DataArray<xarray.DataArray>` 
       The SST meridional gradient data.
   criteria: :py:class:`float<python.float>`, default: `0.80 *1e-5`.
       Empirically-given critical value.
   lat_range: :py:class:`list<python.list>`, default: `[36, 44]`.
       The latitude range of the oceanic frontal zone.
   lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
       The longitude range of the oceanic frontal zone.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   The location of the SAFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017), 
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts 
   and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766    


