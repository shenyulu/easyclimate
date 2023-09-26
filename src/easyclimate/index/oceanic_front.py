"""
Definition of basin-scale oceanic front indexes

Ref: Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017), 
Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic 
fronts and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766
"""

from ..core.utility import get_weighted_spatial_data
from ..core.utility import sort_ascending_latlon_coordinates
import xarray as xr

def calc_intensity_STFZ(data_sst_DtDy: xr.DataArray,
                        criteria = 0.45*1e-5, lat_range = [24, 32], lon_range = [145, 215],
                        lat_dim: str = 'lat', lon_dim: str = 'lon'):
    data_sst_DtDy = sort_ascending_latlon_coordinates(data_sst_DtDy, lat_dim = lat_dim, lon_dim = lon_dim)
    data = data_sst_DtDy.sel(lon = slice(lon_range[0], lon_range[1]), lat = slice(lat_range[0], lat_range[1])).mean(dim = lon_dim)
    return get_weighted_spatial_data(data.where(data > criteria)).mean(dim = lat_dim)

def calc_intensity_SAFZ(data_sst_DtDy: xr.DataArray, 
                        criteria = 0.8*1e-5, lat_range = [36, 44], lon_range = [145, 215],
                        lat_dim: str = 'lat', lon_dim: str = 'lon'):
    data_sst_DtDy = sort_ascending_latlon_coordinates(data_sst_DtDy, lat_dim = lat_dim, lon_dim = lon_dim)
    data = data_sst_DtDy.sel(lon = slice(lon_range[0], lon_range[1]), lat = slice(lat_range[0], lat_range[1])).mean(dim = lon_dim)
    return get_weighted_spatial_data(data.where(data > criteria)).mean(dim = lat_dim)

def calc_location_STFZ(data_sst_DtDy: xr.DataArray, 
                       criteria = 0.45*1e-5, lat_range = [24, 32], lon_range = [145, 215],
                       lat_dim: str = 'lat', lon_dim: str = 'lon'):
    data_sst_DtDy = sort_ascending_latlon_coordinates(data_sst_DtDy, lat_dim = lat_dim, lon_dim = lon_dim)
    data_tmp = data_sst_DtDy.sel(lon = slice(lon_range[0], lon_range[1]), lat = slice(lat_range[0], lat_range[1])).mean(dim = lon_dim)
    data_tmp = data_tmp.where(data_tmp > criteria)
    g_i = get_weighted_spatial_data(data_tmp).sum(dim = lat_dim)
    g_i_times_lat = get_weighted_spatial_data((data_tmp *data_tmp[lat_dim])).sum(dim = lat_dim)
    lct = g_i_times_lat /g_i
    return lct

def calc_location_SAFZ(data_sst_DtDy: xr.DataArray, 
                       criteria = 0.8*1e-5, lat_range = [36, 44], lon_range = [145, 215],
                       lat_dim: str = 'lat', lon_dim: str = 'lon'):
    data_sst_DtDy = sort_ascending_latlon_coordinates(data_sst_DtDy, lat_dim = lat_dim, lon_dim = lon_dim)
    data_tmp = data_sst_DtDy.sel(lon = slice(lon_range[0], lon_range[1]), lat = slice(lat_range[0], lat_range[1])).mean(dim = lon_dim)
    data_tmp = data_tmp.where(data_tmp > criteria)
    g_i = get_weighted_spatial_data(data_tmp).sum(dim = lat_dim)
    g_i_times_lat = get_weighted_spatial_data((data_tmp *data_tmp[lat_dim])).sum(dim = lat_dim)
    lct = g_i_times_lat /g_i
    return lct

def calc_location_line_STFZ(data_sst_DtDy: xr.DataArray, 
                            criteria = 0.45*1e-5, lat_range = [24, 32], lon_range = [145, 215],
                            lat_dim: str = 'lat', lon_dim: str = 'lon'):
    data_sst_DtDy = sort_ascending_latlon_coordinates(data_sst_DtDy, lat_dim = lat_dim, lon_dim = lon_dim)
    data_tmp = data_sst_DtDy.sel(lon = slice(lon_range[0], lon_range[1]), lat = slice(lat_range[0], lat_range[1]))
    data_tmp = data_tmp.where(data_tmp > criteria)
    g_i = get_weighted_spatial_data(data_tmp).sum(dim = lat_dim)
    g_i_times_lat = get_weighted_spatial_data((data_tmp *data_tmp[lat_dim])).sum(dim = lat_dim)
    lct = g_i_times_lat /g_i
    return lct

def calc_location_line_SAFZ(data_sst_DtDy: xr.DataArray, 
                            criteria = 0.8*1e-5, lat_range = [36, 44], lon_range = [145, 215],
                            lat_dim: str = 'lat', lon_dim: str = 'lon'):
    data_sst_DtDy = sort_ascending_latlon_coordinates(data_sst_DtDy, lat_dim = lat_dim, lon_dim = lon_dim)
    data_tmp = data_sst_DtDy.sel(lon = slice(lon_range[0], lon_range[1]), lat = slice(lat_range[0], lat_range[1]))
    data_tmp = data_tmp.where(data_tmp > criteria)
    g_i = get_weighted_spatial_data(data_tmp).sum(dim = lat_dim)
    g_i_times_lat = get_weighted_spatial_data((data_tmp *data_tmp[lat_dim])).sum(dim = lat_dim)
    lct = g_i_times_lat /g_i
    return lct