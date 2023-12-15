"""
Functions for calculation of ocean mixed layer variables.
"""
from __future__ import annotations
import xarray as xr
import numpy as np
from ...core.diff import calc_gradient

def calc_MLD():
    pass
    # https://pyoceans.github.io/python-oceans/_modules/oceans/ocfis.html: https://pyoceans.github.io/python-oceans/_modules/oceans/ocfis.html
    # https://github.com/pyoceans/python-oceans/blob/main/oceans/sw_extras/sw_extras.py: def computemld(fieldso, fieldthetao):


def calc_MLD_temper_tendency(temper_anomaly_vertical, mld, depth_weight, depth_dim = 'depth', time_dim = 'month'):
    """
    
    """
    # Use `xarray.broadcast_to` to extend the dimensionality of the mixed layer depth data by one dimension to correspond to the depth dimension of `temper_anomaly_vertical`
    mld_expanded = xr.broadcast(mld, temper_anomaly_vertical)[0]
    
    # Slice the `temper_anomaly_vertical` data using the `xarray.where` function, keeping only the parts of the data whose depth is less than or equal to the depth of the mixed layer
    temper_mld_vertical = temper_anomaly_vertical.where(temper_anomaly_vertical[depth_dim] <= mld_expanded).compute()
    
    # Temperature mean in the mixed layer
    mixed_layer_temperature = (temper_mld_vertical *depth_weight).sum(dim = depth_dim)

    # Calculate dT'/dt
    mixed_layer_temperature_tendency = calc_gradient(mixed_layer_temperature, dim = time_dim)

    return mixed_layer_temperature_tendency.compute()

def calc_MLD_depth_weighted(temper_vertical, mld, depth_dim = 'depth'):
    """
    
    """
    # Use `xarray.broadcast_to` to extend the dimensionality of the mixed layer depth data by one dimension to correspond to the depth dimension of `temper_vertical`
    mld_expanded = xr.broadcast(mld, temper_vertical)[0]
    depth_expanded = xr.broadcast(temper_vertical[depth_dim], temper_vertical)[0]
    
    # Slice the `temper_vertical` data using the `xarray.where` function, keeping only the parts of the data whose depth is less than or equal to the depth of the mixed layer
    depth_mld = depth_expanded.where(depth_expanded[depth_dim] <= mld_expanded)
    depth_mld_weighted = depth_mld /depth_mld.sum(dim = depth_dim)

    return depth_mld_weighted.compute()

def calc_temper_MLD(temper_vertical, mld, depth_dim = 'depth'):
    """
    
    """
    # Use `xarray.broadcast_to` to extend the dimensionality of the mixed layer depth data by one dimension to correspond to the depth dimension of `temper_vertical`
    mld_expanded = xr.broadcast(mld, temper_vertical)[0]
    
    # Slice the `temper_vertical` data using the `xarray.where` function, keeping only the parts of the data whose depth is less than or equal to the depth of the hybrid layer
    temper_mld = temper_vertical.where(temper_vertical[depth_dim] <= mld_expanded)

    return temper_mld.compute()

def calc_Horizontal_advection(u, v, temper_mld, depth_weight, lat_dim = 'lat', lon_dim = 'lon', depth_dim = 'depth', R = 6370000):
    """
    
    """
    lat = u[lat_dim]
    lon = v[lon_dim]
    dlon = - (np.gradient(lon) *np.pi /180).reshape((1,lon.shape[0]))
    dlat = - (np.gradient(lat) *np.pi /180).reshape((lat.shape[0]))
    coslat = (np.cos(np.array(lat) *np.pi /180)).reshape((lat.shape[0],1))
    dx = R *coslat *dlon
    dy = R *dlat

    dx_array = xr.DataArray(dx, dims = (lat_dim, lon_dim), coords = {lon_dim: lon, lat_dim: lat})
    dy_array = xr.DataArray(dy, dims = lat_dim, coords = {lat_dim: lat})

    # Calculate $u \frac{\partial T}{\partial x}$
    u_advection = u *calc_gradient(temper_mld, dim = lon_dim) *2626560 /dx_array *depth_weight

    # Calculate $v \frac{\partial T}{\partial y}$
    v_advection = v *calc_gradient(temper_mld, dim = lat_dim) *2626560 /dy_array *depth_weight

    return u_advection.sum(dim = depth_dim).compute(), v_advection.sum(dim = depth_dim).compute()

def calc_Vertical_advection(w, temper_vertical, mld, depth_weight, depth_dim = 'depth'):
    """
    
    """ 
    # Calculate $w \frac{\partial T}{\partial z}$
    w_advection = w *calc_gradient(temper_vertical, dim = depth_dim) *2626560 *depth_weight

    # Extracting data within the mixed layer
    w_advection_mld = calc_temper_MLD(w_advection, mld)

    return w_advection_mld.sum(dim = depth_dim).compute()

def calc_Heat_flux(qnet_anomaly, mld, rho_0 = 1027, c_p = 4007):
    """
    
    """
    # Conversion unit: W/m^2 -> degree/month
    qnet_anomaly_degreepermon = qnet_anomaly *2592000

    # Calculate $\frac{q_\mathrm{net}}{\rho_0 c_p MLD}$
    heat_flux_anomaly = qnet_anomaly_degreepermon /(rho_0 *c_p *mld)

    return heat_flux_anomaly.compute()