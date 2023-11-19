:py:mod:`easyclimate.core.tutorial`
===================================

.. py:module:: easyclimate.core.tutorial

.. autoapi-nested-parse::

   Useful for:

   * users learning easyclimate
   * building tutorials in the documentation.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.core.tutorial._construct_cache_dir
   easyclimate.core.tutorial._check_netcdf_engine_installed
   easyclimate.core.tutorial.open_tutorial_dataset



Attributes
~~~~~~~~~~

.. autoapisummary::

   easyclimate.core.tutorial._default_cache_dir_name
   easyclimate.core.tutorial.base_url
   easyclimate.core.tutorial.version
   easyclimate.core.tutorial.external_urls
   easyclimate.core.tutorial.file_formats


.. py:data:: _default_cache_dir_name
   :value: 'easylimate_tutorial_data'

   

.. py:data:: base_url
   :value: 'https://github.com/shenyulu/easyclimate-data'

   

.. py:data:: version
   :value: 'main'

   

.. py:function:: _construct_cache_dir(path)


.. py:data:: external_urls
   :type: dict

   

.. py:data:: file_formats

   

.. py:function:: _check_netcdf_engine_installed(name)


.. py:function:: open_tutorial_dataset(name: str, cache: bool = True, cache_dir: None | str | os.PathLike = None, *, engine: xarray.backends.api.T_Engine = None, **kws) -> xarray.Dataset

   Open a dataset from the online repository (requires internet).

   If a local copy is found then always use that to avoid network traffic.

   Available datasets:

   * ``"air_202201_mon_mean"``: 2m air temperature of the NCEP reanalysis subset
   * ``"hgt_202201_mon_mean"``: Geopotential height of the NCEP reanalysis subset
   * ``"precip_202201_mon_mean"``: Precipitation of the NCEP reanalysis subset
   * ``"pressfc_202201_mon_mean"``: Mean sea surface pressure of the NCEP reanalysis subset
   * ``"shum_202201_mon_mean"``: Absolute humidity of the NCEP reanalysis subset
   * ``"uwnd_202201_mon_mean"``: Zonal wind of the NCEP reanalysis subset
   * ``"vwnd_202201_mon_mean"``: Meridional wind of the NCEP reanalysis subset
   * ``"mini_HadISST_ice"``: Hadley Centre Sea Ice and Sea Surface Temperature data set (HadISST) subset
   * ``"PressQFF_202007271200_872"``: Observational data from European stations (from https://github.com/EXCITED-CO2/xarray-regrid)


   Parameters
   ----------
   name : str
       Name of the file containing the dataset.
       e.g. 'air_202201_mon_mean'
   cache_dir : path-like, optional
       The directory in which to search for and write cached data.
   cache : bool, optional
       If True, then cache data locally for use on subsequent calls
   **kws : dict, optional
       Passed to xarray.open_dataset

   Returns
   -------
   :py:class:`xarray.Dataset<xarray.Dataset>`

   Reference
   --------------
   - Kalnay et al.,The NCEP/NCAR 40-year reanalysis project, Bull. Amer. Meteor. Soc., 77, 437-470, 1996
   - Rayner, N. A.; Parker, D. E.; Horton, E. B.; Folland, C. K.; Alexander, L. V.; Rowell, D. P.; Kent, E. C.; Kaplan, A. (2003) Global analyses of sea surface temperature, sea ice, and night marine air temperature since the late nineteenth century J. Geophys. Res.Vol. 108, No. D14, 4407 10.1029/2002JD002670  (pdf ~9Mb)

   .. seealso::
       - :py:func:`xarray.tutorial.load_dataset<xarray.tutorial.load_dataset>`
       - :py:func:`xarray.open_dataset<xarray.open_dataset>`
       - :py:func:`xarray.load_dataset<xarray.load_dataset>`


