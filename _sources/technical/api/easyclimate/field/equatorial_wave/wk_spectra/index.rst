easyclimate.field.equatorial_wave.wk_spectra
============================================

.. py:module:: easyclimate.field.equatorial_wave.wk_spectra

.. autoapi-nested-parse::

   Wheeler-Kiladis Space-Time Spectra



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

.. py:function:: decompose_symasym(da, lat_dim='lat')

   Decompose an xarray DataArray into symmetric and asymmetric parts about the equator
   along the 'lat' dimension. Supports arrays with 'lat' and optional dimensions like
   'lon', 'time', 'level', etc.

   Parameters
   ----------
   da : xarray.DataArray
       Input array with 'lat' dimension

   Returns
   -------
   xarray.DataArray
       Decomposed array with symmetric part in Southern Hemisphere and asymmetric
       part in Northern Hemisphere


.. py:function:: calc_spectral_coefficients(data: xarray.DataArray, spd: float, nDayWin: float, nDaySkip: float, time_dim: str = 'time', lon_dim: str = 'lon', lat_dim: str = 'lat', max_freq: float = 0.5, max_wn: float = 15)

.. py:function:: draw_wk_anti_analysis(max_freq: float = 0.5, max_wn: float = 15, ax=None, add_xylabel: bool = True, add_central_line: bool = True, add_westward_and_eastward: bool = True, auto_determine_xyrange: bool = True, freq_lines: bool = True, matsuno_modes_labels: bool = True, cpd_lines_levels=[3, 6, 30], matsuno_lines: bool = True, he=[12, 25, 50], meridional_modes=[1])

.. py:function:: draw_wk_sym_analysis(max_freq: float = 0.5, max_wn: float = 15, ax=None, add_xylabel: bool = True, add_central_line: bool = True, add_westward_and_eastward: bool = True, auto_determine_xyrange: bool = True, freq_lines: bool = True, matsuno_modes_labels: bool = True, cpd_lines_levels=[3, 6, 30], matsuno_lines: bool = True, he=[12, 25, 50], meridional_modes=[1])

