:py:mod:`easyclimate.index.npwi`
================================

.. py:module:: easyclimate.index.npwi

.. autoapi-nested-parse::

   Objective globally unified summer monsoon onset or retreat dates.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.index.npwi.calc_index_NPWI
   easyclimate.index.npwi.find_PW_monsoon_region
   easyclimate.index.npwi.cal_NPWI_monsoon_onset
   easyclimate.index.npwi.cal_NPWI_monsoon_detreat



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


