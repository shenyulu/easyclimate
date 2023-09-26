# -*- coding: utf-8 -*-
"""
Quick start
===================================
"""

import xarray as xr
import numpy as np
import easyclimate as ecl

#%%
# Open dataset
uwnd = xr.open_dataset('uwnd_mean.nc').uwnd.rename({'latitude': 'lat', 'longitude': 'lon'}).isel(time = 0)
vwnd = xr.open_dataset('vwnd_mean.nc').vwnd.rename({'latitude': 'lat', 'longitude': 'lon'}).isel(time = 0)

#%%
# Open dataset
uvwnd = xr.Dataset()
uvwnd['u'] = uwnd
uvwnd['v'] = vwnd
uvwnd

#%%
# Open dataset
uwnd_dx = ecl.calc_gradient(uwnd, dim = 'lon')
uwnd_dx

#%%
# Open dataset
uwnd_dx.plot()