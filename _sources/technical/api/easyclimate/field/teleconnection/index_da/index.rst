easyclimate.field.teleconnection.index_da
=========================================

.. py:module:: easyclimate.field.teleconnection.index_da

.. autoapi-nested-parse::

   Arctic Dipole Anomaly (DA/AD)



Functions
---------

.. autoapisummary::

   easyclimate.field.teleconnection.index_da.calc_index_DA_EOF2_Wu_2006


Module Contents
---------------

.. py:function:: calc_index_DA_EOF2_Wu_2006(slp_monthly_data: xarray.DataArray, time_range: slice = slice(None, None), lon_dim: str = 'lon', lat_dim: str = 'lat', lat_range: slice = slice(70, 90), time_dim: str = 'time', random_state: int | None = None, solver: Literal['auto', 'full', 'randomized'] = 'auto', solver_kwargs: dict = {}, normalized: bool = True) -> xarray.DataArray

   The calculation of monthly mean Arctic Dipole Anomaly (DA/AD) index using empirical orthogonal functions (EOFs) method

   .. tip::

       The second EOF mode of SLP anomaly north of 70°N is Arctic Dipole Anomaly (DA/AD) pattern.

   Parameters
   ----------
   slp_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The monthly data of sea level pressure (SLP).
   time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
       The time range of seasonal cycle means to be calculated. The default value is the entire time range.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lat_range: :py:class:`slice <slice>`, default: `slice(20, 90)`.
       The latitude range of computation using EOFs over the Northern Hemisphere. The default value is from :math:`\mathrm{20^{\circ}N}` to :math:`\mathrm{90^{\circ}N}`.
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

   Returns
   -------
   The monthly mean DA/AD index (:py:class:`xarray.DataArray<xarray.DataArray>`).

   Reference
   --------------
   - Wu, B., Wang, J., & Walsh, J. E. (2006). Dipole Anomaly in the Winter Arctic Atmosphere and Its Association with Sea Ice Motion. Journal of Climate, 19(2), 210-225. https://doi.org/10.1175/JCLI3619.1
   - Wu, B., and M. A. Johnson (2007), A seesaw structure in SLP anomalies between the Beaufort Sea and the Barents Sea, Geophys. Res. Lett., 34, L05811, doi: https://doi.org/10.1029/2006GL028333.
   - Wang, J., J. Zhang, E. Watanabe, M. Ikeda, K. Mizobata, J. E. Walsh, X. Bai, and B. Wu (2009), Is the Dipole Anomaly a major driver to record lows in Arctic summer sea ice extent? Geophys. Res. Lett., 36, L05706, doi: https://doi.org/10.1029/2008GL036706.
   - R. Zhang, R. Zhang, Mechanisms for low-frequency variability of summer Arctic sea ice extent, Proc. Natl. Acad. Sci. U.S.A. 112 (15) 4570-4575, https://doi.org/10.1073/pnas.1422296112 (2015).
   - Kapsch, ML., Skific, N., Graversen, R.G. et al. Summers with low Arctic sea ice linked to persistence of spring atmospheric circulation patterns. Clim Dyn 52, 2497–2512 (2019). https://doi.org/10.1007/s00382-018-4279-z
   - Bi, H., Wang, Y., Liang, Y., Sun, W., Liang, X., Yu, Q., Zhang, Z., & Xu, X. (2021). Influences of Summertime Arctic Dipole Atmospheric Circulation on Sea Ice Concentration Variations in the Pacific Sector of the Arctic during Different Pacific Decadal Oscillation Phases. Journal of Climate, 34(8), 3003-3019. https://doi.org/10.1175/JCLI-D-19-0843.1
   - Bi, H., Liang, Y., & Chen, X. (2023). Distinct role of a spring atmospheric circulation mode in the Arctic sea ice decline in summer. Journal of Geophysical Research: Atmospheres, 128, e2022JD037477. https://doi.org/10.1029/2022JD037477

   .. seealso::

       :py:func:`get_EOF_model <easyclimate.core.eof.get_EOF_model>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_da_bbo.py


