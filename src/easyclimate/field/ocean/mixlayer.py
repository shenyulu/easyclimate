"""
The calculation of ocean mixed layer variables.
"""
from __future__ import annotations
import xarray as xr
import numpy as np
import gsw_xarray
from ...core.diff import calc_gradient
from oceans import ocfis

def calc_mixed_layer_depth(
    seawater_temperature_data: xr.DataArray,
    seawater_practical_salinity_data: xr.DataArray,
    criterion: ['temperature', 'density', 'pdvar'] = 'pdvar', 
    depth_dim: str = 'depth',
    lat_dim: str = 'lat',
    lon_dim: str = 'lon',
) -> xr.DataArray:
    """
    Calculate the mixed layer depth.
    
    Parameters
    ----------
    seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{^\circ C}`)
        Vertical seawater temperature data.
    seawater_practical_salinity_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{PSU}`)
        Vertical seawater salinity data (practical salinity).
    criterion: ['temperature', 'density', 'pdvar'], default `'pdvar'`.
        Mixed layer depth criteria

        - **temperature** : Computed based on constant temperature difference criterion, i.e., :math:`CT(0) - T[mld] = 0.5 \\mathrm{^\circ C}`.
        - **density** : Computed based on the constant potential density difference criterion, i.e., :math:`pd[0] - pd[mld] = 0.125` in sigma units.
        - **temperature** : Computed based on variable potential density criterion :math:`pd[0] - pd[mld] = var(T[0], S[0])`, where var is a variable potential density difference which corresponds to constant temperature difference of :math:`0.5 \\mathrm{^\circ C}`.

    depth_dim: :py:class:`str<python.str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.
    lon_dim: :py:class:`str<python.str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str<python.str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

    .. seealso::

        https://pyoceans.github.io/python-oceans/ocfis.html#oceans.ocfis.mld
    """
    ds = xr.Dataset()
    ds['z'] = seawater_temperature_data[depth_dim]
    ds['lat'] = seawater_temperature_data[lat_dim]
    ds['lon'] = seawater_temperature_data[lon_dim]
    ds['SP'] = seawater_practical_salinity_data
    ds['t'] = seawater_temperature_data    # ITS-90 Temperature (Celsius)

    # Height -> seawater pressure
    ds['p'] = gsw_xarray.p_from_z(z = ds['z'] *(-1), lat = ds['lat'])

    # Practical salinity -> Absolute salinity
    ds['SA'] = gsw_xarray.SA_from_SP(SP = ds['SP'], p = ds['p'], lon = ds['lon'], lat = ds['lat'])

    # Conservative temperature
    ds['CT'] = gsw_xarray.CT_from_t(SA = ds['SA'], t = ds['t'], p = ds['p'])

    sea_pressure = ds['p'].broadcast_like(ds['CT'])
    absolute_salinity = ds['SA']
    conservative_temperature = ds['CT']

    def _calc_mld(sa, ct, p):
        if (np.isnan(sa).all() == True) or (np.isnan(ct).all() == True):
            return np.array([np.nan])
        else:
            tmp = ocfis.mld(sa, ct, p, criterion = criterion)[0]
            return np.array([tmp])   
        
    result = xr.apply_ufunc(
        _calc_mld,
        absolute_salinity, conservative_temperature, sea_pressure,
        input_core_dims=[[depth_dim], [depth_dim], [depth_dim]],
        output_core_dims = [["parameter"]],
        output_dtypes=["float64"],
        dask = "parallelized",
        vectorize=True,
        dask_gufunc_kwargs = {"output_sizes": {"parameter": 1}},
    )
    result = result[...,0]
    return result

def calc_MLD_depth_weighted(
    temper_data: xr.DataArray,
    mixed_layer_depth: xr.DataArray,
    depth_dim: str = 'depth'
) -> xr.DataArray:
    """
    Calculate the weights of mixed layer depth. The weights required by the ocean model under non-uniform distribution grids in the depth direction.
    
    Parameters
    ----------
    depth_dim: :py:class:`str<python.str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

    """
    # Use `xarray.broadcast_to` to extend the dimensionality of the mixed layer depth data by one dimension to correspond to the depth dimension of `temper_data`
    mld_expanded = xr.broadcast(mixed_layer_depth, temper_data)[0]
    depth_expanded = xr.broadcast(temper_data[depth_dim], temper_data)[0]
    
    # Slice the `temper_data` data using the `xarray.where` function, keeping only the parts of the data whose depth is less than or equal to the depth of the mixed layer
    depth_mld = depth_expanded.where(depth_expanded[depth_dim] <= mld_expanded)
    depth_mld_weighted = depth_mld /depth_mld.sum(dim = depth_dim)

    return depth_mld_weighted

def calc_MLD_temper_tendency(
    temper_anomaly_data: xr.DataArray, 
    mixed_layer_depth: xr.DataArray, 
    depth_weight: xr.DataArray, 
    depth_dim = 'depth', 
    time_dim = 'time'
) -> xr.DataArray:
    """
    Calculate the tendency of the mixing layer temperature.

    Parameters
    ----------


    depth_dim: :py:class:`str<python.str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.
    time_dim: :py:class:`str<python.str>`, default: `time`.
        The time coordinate dimension name.
    """
    # Use `xarray.broadcast_to` to extend the dimensionality of the mixed layer depth data by one dimension to correspond to the depth dimension of `temper_anomaly_data`
    mld_expanded = xr.broadcast(mixed_layer_depth, temper_anomaly_data)[0]

    # Slice the `temper_anomaly_data` data using the `xarray.where` function, keeping only the parts of the data whose depth is less than or equal to the depth of the mixed layer
    temper_mld_vertical = temper_anomaly_data.where(temper_anomaly_data[depth_dim] <= mld_expanded).compute()
    
    # Temperature mean in the mixed layer
    mixed_layer_temperature = (temper_mld_vertical *depth_weight).sum(dim = depth_dim)
    
    # Calculate dT'/dt
    mixed_layer_temperature_tendency = calc_gradient(mixed_layer_temperature, dim = time_dim)

    return mixed_layer_temperature_tendency

def get_data_within_MLD(
    data_input: xr.DataArray,
    mixed_layer_depth: xr.DataArray,
    depth_dim: str = 'depth'
) -> xr.DataArray:
    """
    Obtain data within the mixed layer


    """
    # Use `xarray.broadcast_to` to extend the dimensionality of the mixed layer depth data by one dimension to correspond to the depth dimension of `data_input`
    mld_expanded = xr.broadcast(mixed_layer_depth, data_input)[0]
    
    # Slice the `data_input` data using the `xarray.where` function, keeping only the parts of the data whose depth is less than or equal to the depth of the hybrid layer
    data_winin_mld = data_input.where(data_input[depth_dim] <= mld_expanded)

    return data_winin_mld

def get_temper_within_MLD(
    temper_data: xr.DataArray,
    mixed_layer_depth: xr.DataArray,
    depth_dim: str = 'depth'    
) -> xr.DataArray:
    """
    Obtain seawater temperature data within the mixing layer.


    """
    return get_data_within_MLD(data_input = temper_data, mixed_layer_depth = mixed_layer_depth, depth_dim = depth_dim)


def calc_MLD_average_horizontal_advection(
    u_monthly_data: xr.DataArray,
    v_monthly_data: xr.DataArray,
    temper_data: xr.DataArray,
    mixed_layer_depth: xr.DataArray,
    depth_weight: xr.DataArray, 
    lat_dim = 'lat', 
    lon_dim = 'lon', 
    depth_dim = 'depth', 
    R = 6370000
) -> xr.DataArray:
    """
    
    """
    lat = u_monthly_data[lat_dim]
    lon = u_monthly_data[lon_dim]
    dlon = - (np.gradient(lon) *np.pi /180).reshape((1,lon.shape[0]))
    dlat = - (np.gradient(lat) *np.pi /180).reshape((lat.shape[0]))
    coslat = (np.cos(np.array(lat) *np.pi /180)).reshape((lat.shape[0],1))
    dx = R *coslat *dlon
    dy = R *dlat

    dx_array = xr.DataArray(dx, dims = (lat_dim, lon_dim), coords = {lon_dim: lon, lat_dim: lat})
    dy_array = xr.DataArray(dy, dims = lat_dim, coords = {lat_dim: lat})

    # Calculate $u \frac{\partial T}{\partial x}$
    u_advection = u_monthly_data *calc_gradient(temper_data, dim = lon_dim) *2626560 /dx_array *depth_weight
    u_advection = get_data_within_MLD(data_input = u_advection, mixed_layer_depth = mixed_layer_depth, depth_dim = depth_dim)

    # Calculate $v \frac{\partial T}{\partial y}$
    v_advection = v_monthly_data *calc_gradient(temper_data, dim = lat_dim) *2626560 /dy_array *depth_weight
    v_advection = get_data_within_MLD(data_input = v_advection, mixed_layer_depth = mixed_layer_depth, depth_dim = depth_dim)

    dataset = xr.Dataset()
    dataset['u_advection'] = u_advection.sum(dim = depth_dim)
    dataset['v_advection'] = v_advection.sum(dim = depth_dim)
    return dataset

def calc_MLD_average_vertical_advection(
    w_monthly_data,
    temper_vertical, 
    mixed_layer_depth, 
    depth_weight, 
    depth_dim = 'depth'
) -> xr.DataArray:
    """
    
    """ 
    # Calculate $w \frac{\partial T}{\partial z}$
    w_advection = w_monthly_data *calc_gradient(temper_vertical, dim = depth_dim) *2626560 *depth_weight

    # Extracting data within the mixed layer
    w_advection_mld = get_data_within_MLD(data_input = w_advection, mixed_layer_depth = mixed_layer_depth, depth_dim = depth_dim)

    return w_advection_mld.sum(dim = depth_dim)

def calc_ocean_surface_heat_flux(
    qnet_anomaly_data: xr.DataArray,
    mixed_layer_depth: xr.DataArray,
    rho_0 = 1027, 
    c_p = 4007,
) -> xr.DataArray:
    """
    
    """
    # Conversion unit: W/m^2 -> degree/month
    qnet_anomaly_degreepermon = qnet_anomaly_data *2592000

    # Calculate $\frac{q_\mathrm{net}}{\rho_0 c_p MLD}$
    heat_flux_anomaly = qnet_anomaly_degreepermon /(rho_0 *c_p *mixed_layer_depth)

    return heat_flux_anomaly