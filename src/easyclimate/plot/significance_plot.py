"""
Functions for draw significant area.
"""
from __future__ import annotations

import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt

def draw_significant_area(
    p_value_data,
    threshold = 0.05,
    lon_dim = 'lon',
    lat_dim = 'lat',
    ax = None,
    hatch_colors = 'k',
    point_density = '...',
    reverse_level_plot = False,
):
    """
    绘制显著性区域（contourf hatch 方法）
    data: 包含 p 值的 DataArray
    ax: 绘制的 axes
    hatch_colors: 更改 hatch 图案颜色，类似于 hatch_colors = ['maroon', 'red', 'darkorange']
    point_density: 点密度或者点的绘制类型，可选值有 '.', '..', '...'
    """

    # Check the projection parameters
    if (type(ax).__name__ == 'GeoAxes'):
        pass
    else:
        raise TypeError("The projection type of the Axes should be `cartopy.crs`, consider to specify the parameter `projection` as `cartopy.crs.xxxxxx`. E.g. `fig, ax = plt.subplots(subplot_kw = {'projection': ccrs.NorthPolarStereo()})`.")

    if ax == None:
        ax = plt.gca()

    if 0 < threshold < 1:
        pass
    else:
        raise ValueError('The parameter `threshold` should be between 0 and 1.')

    proj_trans = ccrs.PlateCarree()
    if(reverse_level_plot == False):
        cs = p_value_data.plot.contourf(
                x = lon_dim, 
                y = lat_dim,
                ax = ax, 
                levels = [0, threshold],
                hatches = [point_density, None], 
                colors = 'none',
                add_colorbar = False, 
                transform = proj_trans,
            )
    elif(reverse_level_plot == True):
        cs = p_value_data.plot.contourf(
                x = lon_dim, 
                y = lat_dim,
                ax = ax, 
                levels = [0, threshold],
                hatches = [None, point_density], 
                colors = 'none',
                add_colorbar = False, 
                transform = proj_trans,
            )

    # 更改 hatch 图案颜色
    # https://github.com/matplotlib/matplotlib/issues/2789/
    # For each level, we set the color of its hatch 
    for i, collection in enumerate(cs.collections):
        collection.set_edgecolor(hatch_colors[i % len(hatch_colors)])
    # Doing this also colors in the box around each level
    # We can remove the colored line around the levels by setting the linewidth to 0
    for collection in cs.collections:
        collection.set_linewidth(0.)

def get_significance_point(
    p_value_data, 
    threshold = 0.05,
    lon_dim = 'lon',
    lat_dim = 'lat',
):  
    """
    
    """
    index = np.where(p_value_data.pvalue <  threshold)
    point_lat = p_value_data[lat_dim][index[0]].data
    point_lon = p_value_data[lon_dim][index[1]].data
    return point_lon, point_lat

# point_lon_op2np, point_lat_op2np = get_significance_point(ttest_result.thin(lat = 2, lon = 2))
# point_lon_op2np_6, point_lat_op2np_6 = get_significance_point(ttest_result.thin(lat = 4, lon = 4))
# ax.scatter(point_lon_op2np_6, point_lat_op2np_6, c = "r", alpha=1, s=1,transform = ccrs.PlateCarree())