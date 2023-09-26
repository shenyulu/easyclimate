"""
Functions for ocean instability.
"""
from __future__ import annotations
import xarray as xr
import numpy as np
import gsw_xarray

def cal_N2_from_t_salt(data_input_t, data_input_salt, time_dim = 'time', depth_dim = 'depth', lat_dim = 'lat', lon_dim = 'lon', p_mid_output = False):
    '''
    http://www.teos-10.org/pubs/gsw/html/gsw_Nsquared.html
    http://www.teos-10.org/pubs/gsw/html/gsw_contents.html
    '''
    data_input_t = data_input_t.transpose(time_dim, depth_dim, lat_dim, lon_dim)
    data_input_salt = data_input_salt.transpose(time_dim, depth_dim, lat_dim, lon_dim)
    time_length, depth_length, lat_length, lon_length = data_input_t.shape

    ds = xr.Dataset()
    ds['z'] = data_input_t[depth_dim]
    ds['lat'] = data_input_t[lat_dim]
    ds['SP'] = data_input_salt
    ds['t'] = data_input_t    # ITS-90 温度（摄氏度）

    # 高度 -> 海水压力
    ds['p'] = gsw_xarray.p_from_z(z = ds['z'] *(-1), lat = ds['lat'])

    # 实用盐度 -> 绝对盐度
    ds['SA'] = gsw_xarray.SA_from_SP(SP = ds['SP'], p = ds['p'], lon = ds['lon'], lat = ds['lat'])

    # 保守温度
    ds['CT'] = gsw_xarray.CT_from_t(SA = ds['SA'], t = ds['t'], p = ds['p'])

    p_tmp = ds['p'].depth.data
    p_tmp_new = p_tmp[np.newaxis, :, np.newaxis, np.newaxis]
    p_tmp_new1 = np.broadcast_to(p_tmp_new, shape = (time_length, depth_length, lat_length, lon_length))
    p_needed = ds['CT'].copy(data = p_tmp_new1, deep = True).where(~np.isnan(ds['t']))

    lat_tmp = ds['lat'].data
    lat_tmp_new = lat_tmp[np.newaxis, np.newaxis, :, np.newaxis]
    lat_tmp_new1 = np.broadcast_to(lat_tmp_new, shape = (time_length, depth_length, lat_length, lon_length))
    lat_needed = ds['CT'].copy(data = lat_tmp_new1, deep = True).where(~np.isnan(ds['t']))

    [N2, p_mid] = gsw_xarray.Nsquared(SA = ds['SA'], CT = ds['CT'], p = p_needed, lat = lat_needed, axis=1)

    N2_dataarray = xr.DataArray(
        N2,
        dims = [time_dim, depth_dim, lat_dim, lon_dim],
        coords = {
            'time': ds['time'].data,
            'depth': ds['depth'].data[:-1],
            'lat': ds['lat'].data,
            'lon': ds['lon'].data,
        }
    )

    p_mid_dataarray = xr.DataArray(
        p_mid,
        dims = [time_dim, depth_dim, lat_dim, lon_dim],
        coords = {
            'time': ds['time'].data,
            'depth': ds['depth'].data[:-1],
            'lat': ds['lat'].data,
            'lon': ds['lon'].data,
        }
    )

    Nsquared = xr.Dataset()
    Nsquared['N2'] = N2_dataarray
    Nsquared['N2'].attrs['name'] = 'Brunt-Vaisala Frequency squared (M-1xN)'
    Nsquared['N2'].attrs['units'] = 'rad^2 s^-2'

    if p_mid_output == True:
        Nsquared['p_mid'] = p_mid_dataarray
        Nsquared['p_mid'].attrs['name'] = 'mid pressure between p grid (M-1xN)'
        Nsquared['p_mid'].attrs['units'] = 'dbar'

    Nsquared.attrs['Attention'] = 'The units of N2 are radians^2 s^-2 however in may textbooks this is abreviated to s^-2 as radians does not have a unit. To convert the frequency to hertz, cycles sec^-1, divide the frequency by 2π, ie N/(2π).'
    return Nsquared

def cal_N2_from_t_salt(data_input_t, data_input_salt):
    '''
    http://www.teos-10.org/pubs/gsw/html/gsw_Nsquared.html
    http://www.teos-10.org/pubs/gsw/html/gsw_contents.html
    '''
    data_input_t = data_input_t.transpose('time', 'depth', 'lat', 'lon')
    data_input_salt = data_input_salt.transpose('time', 'depth', 'lat', 'lon')
    time_length, depth_length, lat_length, lon_length = data_input_t.shape

    ds = xr.Dataset()
    ds['z'] = data_input_t['depth']
    ds['lat'] = data_input_t['lat']
    ds['SP'] = data_input_salt
    ds['t'] = data_input_t    # ITS-90 温度（摄氏度）

    # 高度 -> 海水压力
    ds['p'] = gsw_xarray.p_from_z(
        z = ds['z']*(-1), lat = ds['lat']
    )

    # 实用盐度 -> 绝对盐度
    ds['SA'] = gsw_xarray.SA_from_SP(SP = ds['SP'], p = ds['p'], lon = ds['lon'], lat = ds['lat'])

    # 保守温度
    ds['CT'] = gsw_xarray.CT_from_t(SA = ds['SA'], t = ds['t'], p = ds['p'])

    p_tmp = ds['p'].depth.data
    p_tmp_new = p_tmp[np.newaxis, :, np.newaxis, np.newaxis]
    p_tmp_new1 = np.broadcast_to(p_tmp_new, shape = (time_length, depth_length, lat_length, lon_length))
    p_needed = ds['CT'].copy(data = p_tmp_new1, deep = True).where(~np.isnan(ds['t']))

    prho = gsw_xarray.pot_rho_t_exact(SA = ds['SA'], t = ds['t'], 
                                p = p_needed, p_ref = 0)

    potential_density = xr.Dataset()
    potential_density['prho'] = prho

    potential_density['prho'].attrs['name'] = 'Potential density (not potential density anomaly)'
    potential_density['prho'].attrs['units'] = 'kg/m^3'

    potential_density.attrs['References'] = 'IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.'
    return potential_density