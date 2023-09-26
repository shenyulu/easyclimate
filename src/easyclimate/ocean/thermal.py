"""
Functions for calculation of ocean thermocline variables.
"""
from __future__ import annotations
import xarray as xr
import numpy as np
from ..core.diff import calc_gradient

def calc_temp_thermocline_depth(temp_input, depth_dim = 'depth'):
    """
    Caculate thermocline depth of ocean temperature.

    Parameters
    ----------
    - temp_input: :py:class:`xarray.DataArray<xarray.DataArray>`
        ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
    - depth_dim: str
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    # Caculate gradient
    dTdz = calc_gradient(temp_input, dim = depth_dim)
    depth_weight = calc_gradient(temp_input[depth_dim], dim = depth_dim)
    dTdz_dataarray = dTdz /depth_weight

    # Find the minimum gradient
    dTdz_dataarray_index = dTdz_dataarray.where(dTdz_dataarray == dTdz_dataarray.min(dim = depth_dim))

    # Mismatch
    dTdz_dataarray_bool = dTdz_dataarray.isin(dTdz_dataarray_index)
    dTdz_dataarray_bool_plus1 = dTdz_dataarray_bool.isel({depth_dim: slice(None, -1)}).assign_coords({depth_dim: dTdz_dataarray_bool[depth_dim].data[1:]})
    dTdz_dataarray_bool_minus1 = dTdz_dataarray_bool.isel({depth_dim: slice(1, None)}).assign_coords({depth_dim: dTdz_dataarray_bool[depth_dim].data[:-1]})

    # Find the nearest point to linear interpolation
    depth_3d_array = xr.broadcast(temp_input[depth_dim], temp_input)[0]
    thermal_depth_1 = depth_3d_array.where(dTdz_dataarray_bool).min(dim = depth_dim)
    thermal_depth_2 = depth_3d_array.where(dTdz_dataarray_bool_plus1).min(dim = depth_dim)
    thermal_depth_3 = depth_3d_array.where(dTdz_dataarray_bool_minus1).min(dim = depth_dim)
    thermal_depth = (thermal_depth_1 + thermal_depth_2 + thermal_depth_3) /3

    return thermal_depth

def calc_Dx_depth(data_input, var, depth_dim = 'depth', mask = True):
    """
    
    """
    def _calc_linear_interpolate(x1, x2, y1, y2):
        return (x2 *y1 - x1 *y2)/ (x2 - x1)
    
    data_anormaly = data_input - var
    abs_value = np.abs(data_anormaly)

    # Find minimum index
    min_index = abs_value.fillna(9999).argmin(dim = depth_dim, skipna=True)

    # The closest value to `var` in `data_input` and the depth value
    data_middle = data_anormaly.isel({depth_dim: min_index})
    data_middle_depth = data_middle[depth_dim]

    # The 2nd, 3rd closest value to `var` in `data_input` and the depth value
    data_up = data_input.isel({depth_dim: min_index - 1})
    data_down = data_input.isel({depth_dim: min_index + 1})

    # Find 2nd closest value to `var` in `data_input` and the depth value
    data_nearest = xr.concat([data_up.assign_coords({'type': 0}).expand_dims({'type': 1}),
                          data_down.assign_coords({'type': 1}).expand_dims({'type': 1})],
                          dim = 'type').min(dim = 'type')
    tmp1_anormaly = xr.broadcast(data_input, data_nearest)[1]
    tmp2_deptharray = xr.broadcast(data_input, data_input[depth_dim])[1]
    data_nearest_depth = tmp2_deptharray.where(tmp1_anormaly == data_input).min(dim = depth_dim)

    # Linear interpolation
    isotherm_depth = _calc_linear_interpolate(data_nearest, data_middle, data_middle_depth, data_nearest_depth)

    if mask == True:
        # Create mask
        mask = (data_anormaly.fillna(-9999) < 0).all(dim = depth_dim)
        mask = 1 - mask
        return isotherm_depth.where(mask)
    elif mask == False:
        return isotherm_depth

def calc_D14_depth(data_input, var = 14, depth_dim = 'depth'):
    """
    
    """
    return calc_Dx_depth(data_input, var, depth_dim = depth_dim)

def calc_D17_depth(data_input, var = 17, depth_dim = 'depth'):
    """
    
    """
    return calc_Dx_depth(data_input, var, depth_dim = depth_dim)

def calc_D20_depth(data_input, var = 20, depth_dim = 'depth'):
    """
    
    """
    return calc_Dx_depth(data_input, var, depth_dim = depth_dim)

def calc_D26_depth(data_input, var = 26, depth_dim = 'depth'):
    """
    
    """
    return calc_Dx_depth(data_input, var, depth_dim = depth_dim)

def calc_D28_depth(data_input, var = 28, depth_dim = 'depth'):
    """
    
    """
    return calc_Dx_depth(data_input, var, depth_dim = depth_dim)