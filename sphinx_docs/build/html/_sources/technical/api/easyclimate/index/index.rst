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
   easyclimate.index.monsoon_onsetdate_cal
   easyclimate.index.monsoon_detreatdate_cal
   easyclimate.index.monsoonregion
   easyclimate.index.get_weighted_spatial_data
   easyclimate.index.sort_ascending_latlon_coordinates
   easyclimate.index.calc_intensity_STFZ
   easyclimate.index.calc_intensity_SAFZ
   easyclimate.index.calc_location_STFZ
   easyclimate.index.calc_location_SAFZ
   easyclimate.index.calc_location_line_STFZ
   easyclimate.index.calc_location_line_SAFZ
   easyclimate.index.remove_seasonal_cycle_mean
   easyclimate.index.calc_index_PNA



.. py:function:: calc_index_NPWI(precipitable_water_data)

   https://journals.ametsoc.org/view/journals/clim/17/11/1520-0442_2004_017_2241_gumoar_2.0.co_2.xml


.. py:function:: monsoon_onsetdate_cal(NPWI, thresh=0.618, rollingday=3, n=7)

   NPWI: 三维数组


.. py:function:: monsoon_detreatdate_cal(data, monsoon_onset, thresh=0.618, rollingday=3, n=7)

   NPWI: 三维数组


.. py:function:: monsoonregion(PW)


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


.. py:function:: sort_ascending_latlon_coordinates(data: xr.DataArray | xr.Dataset, lat_dim: str = 'lat', lon_dim: str = 'lon')

   Sort the dimensions `lat`, `lon` in ascending order.


.. py:function:: calc_intensity_STFZ(data_sst_DtDy: xarray.DataArray, criteria=0.45 * 1e-05, lat_range=[24, 32], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon')


.. py:function:: calc_intensity_SAFZ(data_sst_DtDy: xarray.DataArray, criteria=0.8 * 1e-05, lat_range=[36, 44], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon')


.. py:function:: calc_location_STFZ(data_sst_DtDy: xarray.DataArray, criteria=0.45 * 1e-05, lat_range=[24, 32], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon')


.. py:function:: calc_location_SAFZ(data_sst_DtDy: xarray.DataArray, criteria=0.8 * 1e-05, lat_range=[36, 44], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon')


.. py:function:: calc_location_line_STFZ(data_sst_DtDy: xarray.DataArray, criteria=0.45 * 1e-05, lat_range=[24, 32], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon')


.. py:function:: calc_location_line_SAFZ(data_sst_DtDy: xarray.DataArray, criteria=0.8 * 1e-05, lat_range=[36, 44], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon')


.. py:function:: remove_seasonal_cycle_mean(data_input: xarray.DataArray, dim='time', **kwargs) -> xarray.DataArray

   Remove of the seasonal cycle means over the entire time range.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        The data of :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.

   .. caution:: `data_input` must be **monthly** data.

   dim : :py:class:`str<python.str>`
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data. 
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`.


.. py:function:: calc_index_PNA(z500_monthly_data, time_dim, lat_dim='lat', remove_seasonal_cycle=True, save_analysis_path=None, load_analysis_path=None)

       
       


