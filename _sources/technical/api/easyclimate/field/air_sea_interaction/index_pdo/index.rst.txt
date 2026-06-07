easyclimate.field.air_sea_interaction.index_pdo
===============================================

.. py:module:: easyclimate.field.air_sea_interaction.index_pdo

.. autoapi-nested-parse::

   Pacific Decadal Oscillation (PDO) Index



Functions
---------

.. autoapisummary::

   easyclimate.field.air_sea_interaction.index_pdo.calc_index_PDO_EOF1


Module Contents
---------------

.. py:function:: calc_index_PDO_EOF1(sst_monthly_data: xarray.DataArray, time_range: slice = slice(None, None), lon_dim: str = 'lon', lat_dim: str = 'lat', time_dim: str = 'time', random_state: int | None = None, solver: Literal['auto', 'full', 'randomized'] = 'auto', solver_kwargs: dict = {}, normalized: bool = True, detrend_spatial: bool = True) -> xarray.DataArray

   The calculation of monthly mean Pacific Decadal Oscillation (PDO) index using empirical orthogonal functions (EOFs) method over the North Pacific basin.

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
   random_state: :py:class:`int<int>`, default `None`.
       Seed for the random number generator.
   solver: {"auto", "full", "randomized"}, default: "auto".
       Solver to use for the EOFs computation.
   solver_kwargs: :py:class:`dict<dict>`, default `{}`.
       Additional keyword arguments to be passed to the EOFs solver.
   normalized: :py:class:`bool <bool>`, default `True`, optional.
       Whether to standardize the index based on standard deviation over `time_range`.
   normalized: :py:class:`bool <bool>`, default `True`, optional.
       Remove linear trend along time coordinate dimension from data.

   Returns
   -------
   The monthly mean PDO index (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - Zhang, Y., Wallace, J. M., & Battisti, D. S. (1997). ENSO-like Interdecadal Variability: 1900–93. Journal of Climate, 10(5), 1004-1020. https://journals.ametsoc.org/view/journals/clim/10/5/1520-0442_1997_010_1004_eliv_2.0.co_2.xml
   - Mantua, N. J., Hare, S. R., Zhang, Y., Wallace, J. M., & Francis, R. C. (1997). A Pacific Interdecadal Climate Oscillation with Impacts on Salmon Production*. Bulletin of the American Meteorological Society, 78(6), 1069-1080. https://journals.ametsoc.org/view/journals/bams/78/6/1520-0477_1997_078_1069_apicow_2_0_co_2.xml
   - Mantua, N.J., Hare, S.R. The Pacific Decadal Oscillation. Journal of Oceanography 58, 35–44 (2002). https://doi.org/10.1023/A:1015820616384
   - Trenberth, K.E. and Fasullo, J.T. (2013), An apparent hiatus in global warming?. Earth's Future, 1: 19-32. https://doi.org/10.1002/2013EF000165
   - Deser, Clara &, Trenberth, Kevin & National Center for Atmospheric Research Staff (Eds). Last modified 2025-04-29 "The Climate Data Guide: Pacific Decadal Oscillation (PDO): Definition and Indices." Retrieved from https://climatedataguide.ucar.edu/climate-data/pacific-decadal-oscillation-pdo-definition-and-indices
   - Schneider, D. P., C. Deser, J. Fasullo, and K. E. Trenberth (2013), Climate Data Guide Spurs Discovery and Understanding, Eos Trans. AGU, 94(13), 121. https://doi.org/10.1002/2013EO130001

   .. seealso::
       :py:func:`get_EOF_model <easyclimate.core.eof.get_EOF_model>`


