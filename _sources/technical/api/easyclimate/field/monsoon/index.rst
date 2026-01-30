easyclimate.field.monsoon
=========================

.. py:module:: easyclimate.field.monsoon


Submodules
----------

.. toctree::
   :maxdepth: 1

   /technical/api/easyclimate/field/monsoon/bsiso/index
   /technical/api/easyclimate/field/monsoon/index_npwi/index


Functions
---------

.. autoapisummary::

   easyclimate.field.monsoon.calc_index_NPWI
   easyclimate.field.monsoon.find_PW_monsoon_region
   easyclimate.field.monsoon.calc_NPWI_monsoon_onset
   easyclimate.field.monsoon.calc_NPWI_monsoon_retreat
   easyclimate.field.monsoon.draw_bsiso1_phase_space_basemap
   easyclimate.field.monsoon.draw_bsiso2_phase_space_basemap
   easyclimate.field.monsoon.draw_bsiso_phase_space
   easyclimate.field.monsoon.draw_bsiso2_phase_space
   easyclimate.field.monsoon.calc_bsiso1_phase
   easyclimate.field.monsoon.calc_bsiso2_phase
   easyclimate.field.monsoon.get_times_for_phase
   easyclimate.field.monsoon.calc_bsiso_analysis


Package Contents
----------------

.. py:function:: calc_index_NPWI(precipitable_water_daily_data: xarray.DataArray, time_dim: str = 'time') -> xarray.DataArray

   Calculate the normalized precipitable water index (NPWI).

   .. math::
       \mathrm{NPWI} = \frac{\mathrm{PW} - \mathrm{PW_{min}}}{\mathrm{PW_{max}} - \mathrm{PW_{min}}}

   Parameters
   ----------
   precipitable_water_daily_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Daily precipitable water data.
   time_dim: :py:class:`str <str>`.
       The time coordinate dimension name.

   Returns
   -------
   Normalized precipitable water index (NPWI).

   Reference
   --------------
   - Zeng, X., and E. Lu, 2004: Globally Unified Monsoon Onset and Retreat Indexes. J. Climate, 17, 2241–2248, https://journals.ametsoc.org/view/journals/clim/17/11/1520-0442_2004_017_2241_gumoar_2.0.co_2.xml.
   - Tang Xu, Chen Baode, Liang Ping, Qian Weihong. Definition and features of the north edge of Asian summer monsoon. Acta Meteorologica Sinica (Chinese), 2009, (1): 83-89. doi: http://dx.doi.org/10.11676/qxxb2009.009

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_monsoon_npwi.py


.. py:function:: find_PW_monsoon_region(precipitable_water_daily_data: xarray.DataArray, time_dim: str = 'time') -> xarray.DataArray

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
   precipitable_water_daily_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Daily precipitable water data.
   time_dim: :py:class:`str <str>`.
       The time coordinate dimension name.

   Reference
   --------------
   - Zeng, X., and E. Lu, 2004: Globally Unified Monsoon Onset and Retreat Indexes. J. Climate, 17, 2241–2248, https://journals.ametsoc.org/view/journals/clim/17/11/1520-0442_2004_017_2241_gumoar_2.0.co_2.xml.
   - Tang Xu, Chen Baode, Liang Ping, Qian Weihong. Definition and features of the north edge of Asian summer monsoon. Acta Meteorologica Sinica (Chinese), 2009, (1): 83-89. doi: http://dx.doi.org/10.11676/qxxb2009.009

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_monsoon_npwi.py


.. py:function:: calc_NPWI_monsoon_onset(NPWI: xarray.DataArray, thresh: float = 0.618, consecutive_days: int = 3, n: int = 7, lon_dim: str = 'lon', lat_dim: str = 'lat', time_dim: str = 'time') -> xarray.DataArray

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

   thresh: :py:class:`float <float>`, default: `0.618`.
       Golden Ratio value for the threshold value.
   consecutive_days: :py:class:`int<int>`, default: `3`.
       Consecutive days values.
   n: :py:class:`int<int>`, default: `7`.
       :math:`n` of the nine cells centered at cell G in day :math:`d` or (:math:`d \pm 1`).
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name.

   Returns
   -------
   Summer monsoon onset date.

   Reference
   --------------
   - Zeng, X., and E. Lu, 2004: Globally Unified Monsoon Onset and Retreat Indexes. J. Climate, 17, 2241–2248, https://journals.ametsoc.org/view/journals/clim/17/11/1520-0442_2004_017_2241_gumoar_2.0.co_2.xml.
   - Tang Xu, Chen Baode, Liang Ping, Qian Weihong. Definition and features of the north edge of Asian summer monsoon. Acta Meteorologica Sinica (Chinese), 2009, (1): 83-89. doi: http://dx.doi.org/10.11676/qxxb2009.009

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_monsoon_npwi.py


.. py:function:: calc_NPWI_monsoon_retreat(NPWI: xarray.DataArray, monsoon_onset_date: xarray.DataArray, thresh: float = 0.618, consecutive_days: int = 3, n: int = 7, lon_dim: str = 'lon', lat_dim: str = 'lat', time_dim: str = 'time') -> xarray.DataArray

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
       Summer monsoon onset date. The results is generated by :py:func:`easyclimate.index.calc_NPWI_monsoon_onset <easyclimate.index.calc_NPWI_monsoon_onset>`.
   thresh: :py:class:`float <float>`, default: `0.618`.
       Golden Ratio value for the threshold value.
   consecutive_days: :py:class:`int<int>`, default: `3`.
       Consecutive days values.
   n: :py:class:`int<int>`, default: `7`.
       :math:`n` of the nine cells centered at cell G in day :math:`d` or (:math:`d \pm 1`).
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name.

   Returns
   -------
   Summer monsoon retreat date.

   Reference
   --------------
   - Zeng, X., and E. Lu, 2004: Globally Unified Monsoon Onset and Retreat Indexes. J. Climate, 17, 2241–2248, https://journals.ametsoc.org/view/journals/clim/17/11/1520-0442_2004_017_2241_gumoar_2.0.co_2.xml.
   - Tang Xu, Chen Baode, Liang Ping, Qian Weihong. Definition and features of the north edge of Asian summer monsoon. Acta Meteorologica Sinica (Chinese), 2009, (1): 83-89. doi: http://dx.doi.org/10.11676/qxxb2009.009

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_monsoon_npwi.py


.. py:function:: draw_bsiso1_phase_space_basemap(ax=None, add_phase_text: bool = True, add_location_text: bool = True)

   Create a phase space basemap for BSISO1 visualization.

   This function generates a standardized background diagram for BSISO1 analysis,
   plotting phase division lines (1-8), amplitude circles, and geographical annotations.
   The diagram uses normalized principal components (PC1 and PC2) as axes, representing
   the BSISO1 mode in a Cartesian coordinate system (-4 to 4 range). The BSISO1 mode
   captures the northward and eastward propagation of convective anomalies over the
   Asian monsoon region.

   Parameters
   ----------
   ax : :py:class:`matplotlib.axes.Axes`, optional
       Target axes for plotting. If None, uses the current axes (plt.gca()).
   add_phase_text : :py:class:`bool`, default True
       Whether to annotate BSISO phase numbers (1-8) on the diagram.
   add_location_text : :py:class:`bool`, default True
       Whether to annotate geographical regions corresponding to BSISO phases, such as:
       - Bay of Bengal & South China Sea
       - Indian Ocean & East Asia
       - Equatorial Indian Ocean
       - Western North Pacific
       - India & Maritime Continent

   Returns
   -------
   :py:class:`matplotlib.axes.Axes`
       The configured axes object with the BSISO1 phase space basemap.

   Notes
   -----
   - The diagram includes 8 radial lines at 45° intervals to divide phases, a unit circle
     for amplitude threshold (amplitude=1), and labeled axes (PC2, PC1).

   Examples
   --------
   >>> import matplotlib.pyplot as plt
   >>> fig, ax = plt.subplots(figsize=(8, 8))
   >>> draw_bsiso1_phase_space_basemap(ax=ax)
   >>> plt.show()


.. py:function:: draw_bsiso2_phase_space_basemap(ax=None, add_phase_text: bool = True, add_location_text: bool = True)

   Create a phase space basemap for BSISO2 visualization.

   This function generates a standardized background diagram for BSISO2 analysis,
   plotting phase division lines (1-8), amplitude circles, and geographical annotations.
   The diagram uses normalized principal components (PC3 and PC4) as axes, representing
   the BSISO2 mode, which captures higher-frequency intraseasonal variability compared
   to BSISO1.

   Parameters
   ----------
   ax : :py:class:`matplotlib.axes.Axes`, optional
       Target axes for plotting. If None, uses the current axes (plt.gca()).
   add_phase_text : :py:class:`bool`, default ``True``
       Whether to annotate BSISO phase numbers (1-8) on the diagram.
   add_location_text : :py:class:`bool`, default ``True``
       Whether to annotate geographical regions corresponding to BSISO phases, such as:
       - North East Asia
       - South East Asia
       - Philippine Sea
       - India & South China Sea
       - Indian Ocean
       - Western North Pacific
       - Bay of Bengal

   Returns
   -------
   :py:class:`matplotlib.axes.Axes`
       The configured axes object with the BSISO2 phase space basemap.

   Notes
   -----
   - The diagram includes 8 radial lines at 45° intervals, a unit circle for amplitude
     threshold (amplitude=1), and labeled axes (PC4, PC3).
   - BSISO2 typically represents shorter-period oscillations (10-20 days) compared to
     BSISO1 (30-60 days).

   Examples
   --------
   >>> import matplotlib.pyplot as plt
   >>> fig, ax = plt.subplots(figsize=(8, 8))
   >>> draw_bsiso2_phase_space_basemap(ax=ax)
   >>> plt.show()


.. py:function:: draw_bsiso_phase_space(bsiso_data: xarray.Dataset, y_dim: str = 'PC1', x_dim: str = 'PC2', time_dim: str = 'time', ax=None, color='blue', start_text='START', add_start_text: bool = True)

   Visualize BSISO1 phase space trajectory using PC1 and PC2.

   Plots the temporal evolution of BSISO1 phases as a parametric curve in PC1-PC2 space,
   with an optional marker for the starting point. Combines scatter points and connecting
   lines to show the progression of BSISO phases, typically representing 30-60 day
   oscillations.

   Parameters
   ----------
   bsiso_data : :py:class:`xarray.Dataset`
       Dataset containing normalized principal components (PC1, PC2) and a time coordinate.
   y_dim : :py:class:`str`, default "PC1"
       Variable name for the y-axis component (PC1).
   x_dim : :py:class:`str`, default "PC2"
       Variable name for the x-axis component (PC2).
   time_dim : :py:class:`str`, default "time"
       Name of the time coordinate dimension.
   ax : :py:class:`matplotlib.axes.Axes`, optional
       Target axes for plotting. If None, uses the current axes (plt.gca()).
   color : :py:class:`str` or :py:class:`tuple`, default "blue"
       Color for the trajectory and points.
   start_text : :py:class:`str`, default "START"
       Text annotation for the trajectory starting point.
   add_start_text : :py:class:`bool`, default True
       Whether to display the start point annotation.

   Returns
   -------
   :py:class:`matplotlib.axes.Axes`
       Configured axes object with the BSISO1 phase space trajectory.

   Notes
   -----
   - Recommended to use with `draw_bsiso1_phase_space_basemap` for the background diagram.
   - Daily data is recommended for smooth trajectory visualization.

   Examples
   --------
   >>> import xarray as xr
   >>> import matplotlib.pyplot as plt
   >>> fig, ax = plt.subplots(figsize=(8, 8))
   >>> draw_bsiso1_phase_space_basemap(ax=ax)
   >>> draw_bsiso_phase_space(bsiso_data=ds, ax=ax, color="red")
   >>> plt.title("BSISO1 Phase Space Trajectory")
   >>> plt.show()


.. py:function:: draw_bsiso2_phase_space(bsiso_data: xarray.Dataset, pc3_dim: str = 'PC3', pc4_dim: str = 'PC4', time_dim: str = 'time', ax=None, color='blue', start_text='START', add_start_text: bool = True)

   Visualize BSISO2 phase space trajectory using PC3 and PC4.

   Plots the temporal evolution of BSISO2 phases as a parametric curve in PC3-PC4 space,
   with an optional marker for the starting point. Combines scatter points and connecting
   lines to show the progression of BSISO phases, typically representing 10-20 day
   oscillations.

   Parameters
   ----------
   bsiso_data : :py:class:`xarray.Dataset`
       Dataset containing normalized principal components (PC3, PC4) and a time coordinate.
   pc3_dim : :py:class:`str`, default "PC3"
       Variable name for the y-axis component (PC3).
   pc4_dim : :py:class:`str`, default "PC4"
       Variable name for the x-axis component (PC4).
   time_dim : :py:class:`str`, default "time"
       Name of the time coordinate dimension.
   ax : :py:class:`matplotlib.axes.Axes`, optional
       Target axes for plotting. If None, uses the current axes (plt.gca()).
   color : :py:class:`str` or tuple, default "blue"
       Color for the trajectory and points.
   start_text : :py:class:`str`, default "START"
       Text annotation for the trajectory starting point.
   add_start_text : :py:class:`bool`, default True
       Whether to display the start point annotation.

   Returns
   -------
   :py:class:`matplotlib.axes.Axes`
       Configured axes object with the BSISO2 phase space trajectory.

   Notes
   -----
   - Recommended to use with `draw_bsiso2_phase_space_basemap` for the background diagram.
   - Daily data is recommended for smooth trajectory visualization.

   Examples
   --------
   >>> import xarray as xr
   >>> import matplotlib.pyplot as plt
   >>> fig, ax = plt.subplots(figsize=(8, 8))
   >>> draw_bsiso2_phase_space_basemap(ax=ax)
   >>> draw_bsiso2_phase_space(bsiso_data=ds, ax=ax, color="red")
   >>> plt.title("BSISO2 Phase Space Trajectory")
   >>> plt.show()


.. py:function:: calc_bsiso1_phase(ds: xarray.DataArray, amplitude_threshold: float = 1.5, pc1_dim: str = 'PC1', pc2_dim: str = 'PC2')

   Calculate BSISO1 phase based on PC1 and PC2 values.

   Computes the BSISO1 phase (0-8) based on the principal components PC1 and PC2,
   following the classification rules defined by Lee et al. (2013). Adds a boolean
   variable 'is_event' to indicate if the amplitude (:math:`\sqrt{\mathrm{PC1}^2 + \mathrm{PC2}^2}`) exceeds
   the threshold, and a ``'phase'`` variable for the BSISO phase.

   Parameters
   ----------
   ds : :py:class:`xarray.Dataset`
       Dataset containing PC1 and PC2 variables with a time dimension.
   amplitude_threshold : :py:class:`float`, default 1.5
       Threshold for event detection based on amplitude.
   pc1_dim : :py:class:`str`, default "PC1"
       Variable name for PC1 in the dataset.
   pc2_dim : :py:class:`str`, default "PC2"
       Variable name for PC2 in the dataset.

   Returns
   -------
   :py:class:`xarray.Dataset`
       Dataset with added 'phase' (int, 0-8) and 'is_event' (boolean) variables.

   Notes
   -----
   Phase classification:

   - Phase 0: Non-event (amplitude <= threshold)
   - Phase 1: :math:`\mathrm{PC1} < 0, \mathrm{PC2} < 0, \mathrm{PC1} > \mathrm{PC2}`
   - Phase 2: :math:`\mathrm{PC1} < 0, \mathrm{PC2} < 0, \mathrm{PC1} < \mathrm{PC2}`
   - Phase 3: :math:`\mathrm{PC1} < 0, \mathrm{PC2} > 0, |\mathrm{PC1}| > \mathrm{PC2}`
   - Phase 4: :math:`\mathrm{PC1} < 0, \mathrm{PC2} > 0, |\mathrm{PC1}| < \mathrm{PC2}`
   - Phase 5: :math:`\mathrm{PC1} > 0, \mathrm{PC2} > 0, \mathrm{PC1} < \mathrm{PC2}`
   - Phase 6: :math:`\mathrm{PC1} > 0, \mathrm{PC2} > 0, \mathrm{PC1} > \mathrm{PC2}`
   - Phase 7: :math:`\mathrm{PC1} > 0, \mathrm{PC2} < 0, \mathrm{PC1} > |\mathrm{PC2}|`
   - Phase 8: :math:`\mathrm{PC1} > 0, \mathrm{PC2} < 0, \mathrm{PC1} < |\mathrm{PC2}|`

   Examples
   --------
   >>> ds = calc_bsiso1_phase(ds, amplitude_threshold=1.5)
   >>> print(ds['phase'])


.. py:function:: calc_bsiso2_phase(ds: xarray.DataArray, amplitude_threshold: float = 1.5, pc3_dim: str = 'PC3', pc4_dim: str = 'PC4')

   Calculate BSISO2 phase based on PC3 and PC4 values.

   Wrapper function for `calc_bsiso1_phase` to compute BSISO2 phases (0-8) using
   principal components PC3 and PC4, which capture higher-frequency (10-20 day)
   intraseasonal variability.

   Parameters
   ----------
   ds : :py:class:`xarray.Dataset`
       Dataset containing PC3 and PC4 variables with a time dimension.
   amplitude_threshold : :py:class:`float`, default 1.5
       Threshold for event detection based on amplitude.
   pc3_dim : :py:class:`str`, default "PC3"
       Variable name for PC3 in the dataset.
   pc4_dim : :py:class:`str`, default "PC4"
       Variable name for PC4 in the dataset.

   Returns
   -------
   :py:class:`xarray.Dataset`
       Dataset with added 'phase' (int, 0-8) and 'is_event' (boolean) variables.

   Notes
   -----
   - Delegates to `calc_bsiso1_phase` with PC3 and PC4 as inputs.
   - BSISO2 phases follow the same classification rules as BSISO1 but use different PCs.

   Examples
   --------
   >>> ds = calc_bsiso2_phase(ds, amplitude_threshold=1.5)
   >>> print(ds['phase'])


.. py:function:: get_times_for_phase(ds: xarray.Dataset, phase_value: int, time_dim: str = 'time', phase_dim: str = 'phase')

   Retrieve time points for a specified BSISO phase.

   Extracts time coordinates from the dataset where the 'phase' variable matches the
   specified phase value (0-8), useful for composite analysis of BSISO events.

   Parameters
   ----------
   ds : :py:class:`xarray.Dataset`
       Dataset containing 'phase' variable with a time dimension.
   phase_value : :py:class:`int`
       Target phase value to query (0-8).
   time_dim : :py:class:`str`, default "time"
       Name of the time coordinate dimension.
   phase_dim : :py:class:`str`, default "phase"
       Name of the phase variable.

   Returns
   -------
   :py:class:`xarray.DataArray`
       Array of time coordinates where the phase equals ``phase_value``.

   Notes
   -----
   - Returns an empty DataArray if no time points match the specified phase.
   - Assumes the 'phase' variable exists in the dataset.

   Examples
   --------
   >>> times = get_times_for_phase(ds, phase_value=1)
   >>> print(times)


.. py:function:: calc_bsiso_analysis(olr_data: xarray.DataArray, u850_data: xarray.DataArray, v850_data: xarray.DataArray, daily_cycle_mean_time_range: slice = slice(None, None), extract_time_range: slice = slice(None, None), harmonics_num: int = 3, threshold: float = 0.05, time_dim: str = 'time', lat_dim: str = 'lat', lon_dim: str = 'lon') -> easyclimate.core.datanode.DataNode

   Perform BSISO analysis using OLR and 850-hPa wind data.

   Computes multivariate BSISO indices (BSISO1 and BSISO2) using outgoing longwave
   radiation (OLR) and 850-hPa zonal (u850) and meridional (v850) winds, following
   the methodology of Lee et al. (2013). The function processes data to remove
   seasonal cycles and low-frequency signals, performs multivariate EOF analysis,
   and derives phase composites for BSISO1 and BSISO2.

   Parameters
   ----------
   olr_data : :py:class:`xarray.DataArray`
       Outgoing longwave radiation data with time, latitude, and longitude dimensions.
   u850_data : :py:class:`xarray.DataArray`
       850-hPa zonal wind data with the same dimensions as ``olr_data``.
   v850_data : :py:class:`xarray.DataArray`
       850-hPa meridional wind data with the same dimensions as ``olr_data``.
   daily_cycle_mean_time_range : :py:class:`slice`, default ``slice(None, None)``
       Time range for computing the daily annual cycle mean.
   extract_time_range : :py:class:`slice`, default ``slice(None, None)``
       Time range for extracting data for analysis.
   harmonics_num : :py:class:`int`, default 3
       Number of harmonics for smoothing the daily annual cycle.
   threshold : :py:class:`float`, default 0.05
       P-value threshold for statistical significance in composite analysis.
   time_dim : :py:class:`str`, default "time"
       Name of the time dimension.
   lat_dim : :py:class:`str`, default "lat"
       Name of the latitude dimension.
   lon_dim : :py:class:`str`, default "lon"
       Name of the longitude dimension.

   Returns
   -------
   :py:class:`easyclimate.DataNode`
       Hierarchical data structure containing:

       - Explained variance ratios
       - Principal components (PC1-PC4)
       - Normalized PCs
       - EOF patterns for OLR and winds
       - Lead-lag correlation coefficients
       - Phase information
       - Composite analysis results for BSISO1 and BSISO2 phases

   Notes
   -----
   - The analysis focuses on May-October data to capture the boreal summer monsoon season.
   - BSISO1 represents 30-60 day oscillations, while BSISO2 captures 10-20 day oscillations.
   - Progress indicators are logged to track major computational steps.

   References
   ----------
   - Lee, JY. (이준이), Wang, B., Wheeler, M.C. et al. Real-time multivariate indices for the boreal summer intraseasonal oscillation over the Asian summer monsoon region. Clim Dyn 40, 493–509 (2013). https://doi.org/10.1007/s00382-012-1544-4
   - Wheeler, M. C., & Hendon, H. H. (2004). An All-Season Real-Time Multivariate MJO Index: Development of an Index for Monitoring and Prediction. Monthly Weather Review, 132(8), 1917-1932. https://doi.org/10.1175/1520-0493(2004)132<1917:AARMMI>2.0.CO;2

   Examples
   --------
   >>> result = calc_bsiso_analysis(olr_data, u850_data, v850_data)
   >>> print(result["Phase/PC1_2"])

   .. seealso::

       - https://iprc.soest.hawaii.edu/users/jylee/bsiso/


