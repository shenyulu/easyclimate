"""
ENSO Index

https://psl.noaa.gov/enso/dashboard.html
https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni
"""
from ..core.utility import sort_ascending_latlon_coordinates
import xarray as xr

def calc_index_nino1plus2(data_sst_anomaly: xr.DataArray, lat_dim: str = 'lat', lon_dim: str = 'lon'):
    data_sst_anomaly = sort_ascending_latlon_coordinates(data_sst_anomaly, lat_dim = lat_dim, lon_dim = lon_dim)
    nino1plus2_index = data_sst_anomaly.sel({lon_dim: slice(-10, 0), lat_dim: slice(270, 280)})
    nino1plus2_index.name = 'Nino1+2_index'
    return nino1plus2_index

def calc_index_nino3(data_sst_anomaly: xr.DataArray, lat_dim: str = 'lat', lon_dim: str = 'lon'):
    data_sst_anomaly = sort_ascending_latlon_coordinates(data_sst_anomaly, lat_dim = lat_dim, lon_dim = lon_dim)
    nino3_index = data_sst_anomaly.sel({lon_dim: slice(-5, 5), lat_dim: slice(210, 270)})
    nino3_index.name = 'Nino3_index'
    return nino3_index

def calc_index_nino34(data_sst_anomaly: xr.DataArray, lat_dim: str = 'lat', lon_dim: str = 'lon'):
    data_sst_anomaly = sort_ascending_latlon_coordinates(data_sst_anomaly, lat_dim = lat_dim, lon_dim = lon_dim)
    nino34_index = data_sst_anomaly.sel({lon_dim: slice(-5, 5), lat_dim: slice(190, 240)})
    nino34_index = nino34_index.rolling(time = 5).mean()
    nino34_index.name = 'Nino34_index'
    return nino34_index

def calc_index_OMI(data_sst_anomaly: xr.DataArray, lat_dim: str = 'lat', lon_dim: str = 'lon'):
    data_sst_anomaly = sort_ascending_latlon_coordinates(data_sst_anomaly, lat_dim = lat_dim, lon_dim = lon_dim)
    omi_index = data_sst_anomaly.sel({lon_dim: slice(-5, 5), lat_dim: slice(190, 240)})
    omi_index = omi_index.rolling(time = 3).mean()
    omi_index.name = 'Oceanic_Nino_Index_index'
    return omi_index

def calc_index_nino4(data_sst_anomaly: xr.DataArray, lat_dim: str = 'lat', lon_dim: str = 'lon'):
    data_sst_anomaly = sort_ascending_latlon_coordinates(data_sst_anomaly, lat_dim = lat_dim, lon_dim = lon_dim)
    nino4_index = data_sst_anomaly.sel({lon_dim: slice(-5, 5), lat_dim: slice(160, 210)})
    nino4_index.name = 'Nino4_index'
    return nino4_index