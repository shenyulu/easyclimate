:py:mod:`easyclimate.index.oceanic_front`
=========================================

.. py:module:: easyclimate.index.oceanic_front

.. autoapi-nested-parse::

   Definition of basin-scale oceanic front indexes

   Ref: Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017), 
   Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic 
   fronts and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.index.oceanic_front.calc_intensity_STFZ
   easyclimate.index.oceanic_front.calc_intensity_SAFZ
   easyclimate.index.oceanic_front.calc_location_STFZ
   easyclimate.index.oceanic_front.calc_location_SAFZ
   easyclimate.index.oceanic_front.calc_location_line_STFZ
   easyclimate.index.oceanic_front.calc_location_line_SAFZ



.. py:function:: calc_intensity_STFZ(data_sst_DtDy: xarray.DataArray, criteria=0.45 * 1e-05, lat_range=[24, 32], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon')


.. py:function:: calc_intensity_SAFZ(data_sst_DtDy: xarray.DataArray, criteria=0.8 * 1e-05, lat_range=[36, 44], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon')


.. py:function:: calc_location_STFZ(data_sst_DtDy: xarray.DataArray, criteria=0.45 * 1e-05, lat_range=[24, 32], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon')


.. py:function:: calc_location_SAFZ(data_sst_DtDy: xarray.DataArray, criteria=0.8 * 1e-05, lat_range=[36, 44], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon')


.. py:function:: calc_location_line_STFZ(data_sst_DtDy: xarray.DataArray, criteria=0.45 * 1e-05, lat_range=[24, 32], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon')


.. py:function:: calc_location_line_SAFZ(data_sst_DtDy: xarray.DataArray, criteria=0.8 * 1e-05, lat_range=[36, 44], lon_range=[145, 215], lat_dim: str = 'lat', lon_dim: str = 'lon')


