:py:mod:`easyclimate.index.enso`
================================

.. py:module:: easyclimate.index.enso

.. autoapi-nested-parse::

   ENSO Index

   https://psl.noaa.gov/enso/dashboard.html
   https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.index.enso.calc_index_nino1plus2
   easyclimate.index.enso.calc_index_nino3
   easyclimate.index.enso.calc_index_nino34
   easyclimate.index.enso.calc_index_OMI
   easyclimate.index.enso.calc_index_nino4



.. py:function:: calc_index_nino1plus2(data_sst_anomaly: xarray.DataArray, lat_dim: str = 'lat', lon_dim: str = 'lon')


.. py:function:: calc_index_nino3(data_sst_anomaly: xarray.DataArray, lat_dim: str = 'lat', lon_dim: str = 'lon')


.. py:function:: calc_index_nino34(data_sst_anomaly: xarray.DataArray, lat_dim: str = 'lat', lon_dim: str = 'lon')


.. py:function:: calc_index_OMI(data_sst_anomaly: xarray.DataArray, lat_dim: str = 'lat', lon_dim: str = 'lon')


.. py:function:: calc_index_nino4(data_sst_anomaly: xarray.DataArray, lat_dim: str = 'lat', lon_dim: str = 'lon')


