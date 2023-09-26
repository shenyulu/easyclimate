"""
Functions for cartopy projection.
"""
from __future__ import annotations

import cartopy
import matplotlib.pyplot as plt
from ..core.utility import *

# cartopy version check
check_return = assert_compared_version(cartopy.__version__, '0.20')
if check_return == 1:
    pass
else:
    print("Cartopy version is not greater than 0.20, please update cartopy package. You can use Conda to update: `conda install -c conda-forge cartopy`")

import cartopy.crs as ccrs
import matplotlib.ticker as ticker
import numpy as np
from geocat.viz import util as gvutil


def draw_Circlemap_PolarStereo(*, 
                          ax = None, 
                          lonstep = 30, 
                          latstep = 20, 
                          lat_range = [0, 40], 
                          gridcolor = 'black',
                          linestyle = "--",
                          x_inline=False, y_inline=True,
                          xlabel_style = {},
                          ylabel_style = {},
                          draw_labels=True,
                          correct_pad = {},
                        ):
    # Get Axes
    if(ax == None):
        ax = plt.gca()
    else:
        pass    
  
    # Check the projection parameters
    if (type(ax.projection) == 'cartopy.crs.NorthPolarStereo') or (type(ax.projection) == 'cartopy.crs.SouthPolarStereo'):
        pass
    else:
        raise TypeError("The projection type of the Axes should be `cartopy.crs.NorthPolarStereo` or `cartopy.crs.SouthPolarStereo`, consider to specify the parameter `projection` as `cartopy.crs.NorthPolarStereo` or `cartopy.crs.SouthPolarStereo`. E.g. `fig, ax = plt.subplots(subplot_kw = {'projection': ccrs.NorthPolarStereo()})`.")


    gvutil.set_map_boundary(ax, [-180, 180], lat_range, **correct_pad) #west_pad=2, , east_pad=2, ,north_pad=2,south_pad=2
    
    gl = ax.gridlines(ccrs.PlateCarree(),
                    draw_labels=draw_labels,
                    linestyle=linestyle,
                    color=gridcolor,
                    x_inline=x_inline, y_inline=y_inline,
                    xlocs = ticker.FixedLocator(np.arange(-180, 180, lonstep)),
                    ylocs = ticker.FixedLocator(np.arange(0, 90, latstep))
                    )
    gl.xlabel_style = xlabel_style
    gl.ylabel_style = ylabel_style

def add_cyclic(data2d, inter = 1.125):
    # SINTEX Model: 1.125
    # GPCP: 2.5
    # ERA5: 0.25
    data_temp = data2d.pad(lon=(0, 1), mode="wrap")
    draw_data = data_temp.assign_coords(lon = np.arange(0, 360 + inter, inter))
    return draw_data