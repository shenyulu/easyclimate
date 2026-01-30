easyclimate.field.equatorial_wave.wk_spectra
============================================

.. py:module:: easyclimate.field.equatorial_wave.wk_spectra

.. autoapi-nested-parse::

   Wheeler-Kiladis Space-Time Spectra

   This module provides functions for analyzing and visualizing Wheeler-Kiladis space-time spectra,
   including signal processing, symmetric/asymmetric decomposition, spectral analysis, and plotting
   of equatorial wave dispersion relationships.

   .. seealso::

       - https://github.com/mmaiergerber/wk_spectra
       - Wheeler, M., & Kiladis, G. N. (1999). Convectively Coupled Equatorial Waves: Analysis of Clouds and Temperature in the Wavenumber–Frequency Domain. Journal of the Atmospheric Sciences, 56(3), 374-399. https://journals.ametsoc.org/view/journals/atsc/56/3/1520-0469_1999_056_0374_ccewao_2.0.co_2.xml
       - Kiladis, G. N., M. C. Wheeler, P. T. Haertel, K. H. Straub, and P. E. Roundy (2009), Convectively coupled equatorial waves, Rev. Geophys., 47, RG2003, doi:https://doi.org/10.1029/2008RG000266
       - Wheeler, M. C., & Nguyen, H. (2015). TROPICAL METEOROLOGY AND CLIMATE | Equatorial Waves. In Encyclopedia of Atmospheric Sciences (pp. 102–112). Elsevier. https://doi.org/10.1016/B978-0-12-382225-3.00414-X
       - Yoshikazu Hayashi, A Generalized Method of Resolving Disturbances into Progressive and Retrogressive Waves by Space Fourier and Time Cross-Spectral Analyses, Journal of the Meteorological Society of Japan. Ser. II, 1971, Volume 49, Issue 2, Pages 125-128, Released on J-STAGE May 27, 2008, Online ISSN 2186-9057, Print ISSN 0026-1165, https://doi.org/10.2151/jmsj1965.49.2_125, https://www.jstage.jst.go.jp/article/jmsj1965/49/2/49_2_125/_article/-char/en



Functions
---------

.. autoapisummary::

   easyclimate.field.equatorial_wave.wk_spectra.remove_dominant_signals
   easyclimate.field.equatorial_wave.wk_spectra.decompose_symasym
   easyclimate.field.equatorial_wave.wk_spectra.calc_spectral_coefficients
   easyclimate.field.equatorial_wave.wk_spectra.draw_wk_anti_analysis
   easyclimate.field.equatorial_wave.wk_spectra.draw_wk_sym_analysis


Module Contents
---------------

.. py:function:: remove_dominant_signals(data: xarray.DataArray, spd: float, nDayWin: float, nDaySkip: float, time_dim: str = 'time', lon_dim: str = 'lon', lat_dim: str = 'lat') -> xarray.DataArray

   Removes the dominant signals by removing the long term linear trend (conserving the mean) and
   by eliminating the annual cycle by removing all time periods less than a corresponding critical frequency.

   Parameters
   ----------
   data : :py:class:`xarray.DataArray<xarray.DataArray>`
       Input data array with time, lat, lon dimensions

   .. caution:: `data` must be **daily** data, and the length of time should be longer than one year (i.e., >365 day).

   spd : :py:class:`float <float>`
       Samples per day (frequency of observations)
   nDayWin : :py:class:`float <float>`
       Number of days in each analysis window
   nDaySkip : :py:class:`float <float>`
       Number of days to skip between windows
   time_dim : :py:class:`str <str>`, optional
       Name of time dimension, default 'time'
   lon_dim : :py:class:`str <str>`, optional
       Name of longitude dimension, default 'lon'
   lat_dim : :py:class:`str <str>`, optional
       Name of latitude dimension, default 'lat'

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`
       Data with dominant signals removed

   Example
   -------
   >>> data = xr.DataArray(np.random.rand(730, 10, 20),
   ...                    dims=['time', 'lat', 'lon'])
   >>> cleaned = remove_dominant_signals(data, spd=1, nDayWin=90, nDaySkip=30)

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_wk_spectra.py


.. py:function:: decompose_symasym(da, lat_dim='lat')

   Decompose data into symmetric and asymmetric parts about the equator.

   The symmetric part is stored in the Southern Hemisphere and the asymmetric
   part in the Northern Hemisphere.

   Parameters
   ----------
   da : :py:class:`xarray.DataArray<xarray.DataArray>`
       Input array with latitude dimension
   lat_dim : :py:class:`str <str>`, optional
       Name of latitude dimension, default 'lat'

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`
       Array with decomposed components (symmetric in SH, asymmetric in NH)

   Example
   -------
   >>> data = xr.DataArray(np.random.rand(10, 20), dims=['lat', 'lon'])
   >>> decomposed = decompose_symasym(data)

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_wk_spectra.py


.. py:function:: calc_spectral_coefficients(data: xarray.DataArray, spd: float, nDayWin: float, nDaySkip: float, time_dim: str = 'time', lon_dim: str = 'lon', lat_dim: str = 'lat', max_freq: float = 0.5, max_wn: float = 15)

   Calculate Wheeler-Kiladis spectral coefficients.

   Parameters
   ----------
   data : :py:class:`xarray.DataArray<xarray.DataArray>`
       Input data array.

   .. caution:: `data` must be **daily** data, and the length of time should be longer than one year (i.e., >365 day).

   spd : :py:class:`float <float>`
       Samples per day
   nDayWin : :py:class:`float <float>`
       Number of days in analysis window
   nDaySkip : :py:class:`float <float>`
       Number of days to skip between windows
   time_dim : :py:class:`str <str>`, optional
       Time dimension name, default 'time'
   lon_dim : :py:class:`str <str>`, optional
       Longitude dimension name, default 'lon'
   lat_dim : :py:class:`str <str>`, optional
       Latitude dimension name, default 'lat'
   max_freq : :py:class:`float <float>`, optional
       Maximum frequency to return (CPD), default 0.5
   max_wn : :py:class:`float <float>`, optional
       Maximum wavenumber to return, default 15

   Returns
   -------
   :py:class:`xarray.Dataset<xarray.Dataset>`
       Dataset containing:

       - ``psumanti_r``: Antisymmetric power spectrum (background removed)
       - ``psumsym_r``: Symmetric power spectrum (background removed)

   Example
   -------
   >>> data = xr.DataArray(np.random.rand(730, 10, 20),
   ...                    dims=['time', 'lat', 'lon'])
   >>> spectra = calc_spectral_coefficients(data, spd=1, nDayWin=90, nDaySkip=30)

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_wk_spectra.py


.. py:function:: draw_wk_anti_analysis(max_freq: float = 0.5, max_wn: float = 15, ax=None, add_xylabel: bool = True, add_central_line: bool = True, add_westward_and_eastward: bool = True, auto_determine_xyrange: bool = True, freq_lines: bool = True, matsuno_modes_labels: bool = True, cpd_lines_levels: list = [3, 6, 30], matsuno_lines: bool = True, he: list = [12, 25, 50], meridional_modes: list = [1])

   Plot antisymmetric Wheeler-Kiladis analysis with Matsuno dispersion curves.

   Parameters
   ----------
   max_freq : :py:class:`float <float>`, optional
       Maximum frequency to plot (CPD), default 0.5
   max_wn : :py:class:`float <float>`, optional
       Maximum wavenumber to plot, default 15
   ax : matplotlib.axes.Axes, optional
       Axes to plot on, creates new if None
   add_xylabel : :py:class:`bool <bool>`, optional
       Add x/y labels, default True
   add_central_line : :py:class:`bool <bool>`, optional
       Add central vertical line, default True
   add_westward_and_eastward : :py:class:`bool <bool>`, optional
       Add eastward/westward labels, default True
   auto_determine_xyrange : :py:class:`bool <bool>`, optional
       Auto-set axis ranges, default True
   freq_lines : :py:class:`bool <bool>`, optional
       Add frequency lines, default True
   matsuno_modes_labels : :py:class:`bool <bool>`, optional
       Add Matsuno mode labels, default True
   cpd_lines_levels : list, optional
       Periods (days) for frequency lines, default [3, 6, 30]
   matsuno_lines : :py:class:`bool <bool>`, optional
       Plot Matsuno dispersion curves, default True
   he : list, optional
       Equivalent depths for Matsuno curves, default [12, 25, 50]
   meridional_modes : list, optional
       Meridional mode numbers, default [1]

   Example
   -------
   >>> fig, ax = plt.subplots()
   >>> draw_wk_anti_analysis(ax=ax)

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_wk_spectra.py


.. py:function:: draw_wk_sym_analysis(max_freq: float = 0.5, max_wn: float = 15, ax=None, add_xylabel: bool = True, add_central_line: bool = True, add_westward_and_eastward: bool = True, auto_determine_xyrange: bool = True, freq_lines: bool = True, matsuno_modes_labels: bool = True, cpd_lines_levels: list = [3, 6, 30], matsuno_lines: bool = True, he: list = [12, 25, 50], meridional_modes: list = [1])

   Plot symmetric Wheeler-Kiladis analysis with Matsuno dispersion curves.

   Parameters
   ----------
   max_freq : :py:class:`float <float>`, optional
       Maximum frequency to plot (CPD), default 0.5
   max_wn : :py:class:`float <float>`, optional
       Maximum wavenumber to plot, default 15
   ax : matplotlib.axes.Axes, optional
       Axes to plot on, creates new if None
   add_xylabel : :py:class:`bool <bool>`, optional
       Add x/y labels, default True
   add_central_line : :py:class:`bool <bool>`, optional
       Add central vertical line, default True
   add_westward_and_eastward : :py:class:`bool <bool>`, optional
       Add eastward/westward labels, default True
   auto_determine_xyrange : :py:class:`bool <bool>`, optional
       Auto-set axis ranges, default True
   freq_lines : :py:class:`bool <bool>`, optional
       Add frequency lines, default True
   matsuno_modes_labels : :py:class:`bool <bool>`, optional
       Add Matsuno mode labels, default True
   cpd_lines_levels : list, optional
       Periods (days) for frequency lines, default [3, 6, 30]
   matsuno_lines : :py:class:`bool <bool>`, optional
       Plot Matsuno dispersion curves, default True
   he : list, optional
       Equivalent depths for Matsuno curves, default [12, 25, 50]
   meridional_modes : list, optional
       Meridional mode numbers, default [1]

   Example
   -------
   >>> fig, ax = plt.subplots()
   >>> draw_wk_sym_analysis(ax=ax)

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_wk_spectra.py


