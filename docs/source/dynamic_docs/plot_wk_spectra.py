# -*- coding: utf-8 -*-
"""
.. _wk_spectra_example:

Wheeler-Kiladis Space-Time Spectra
============================================

Extract equatorial waves by filtering in the Wheeler-Kiladis wavenumber-frequency domain.

Before proceeding with all the steps, first import some necessary libraries and packages
"""
import xarray as xr
import matplotlib.pyplot as plt
import easyclimate as ecl


# %%
#
data = xr.open_dataset('olr_smooth_data.nc')['olr'].sel(lat = slice(-15, 15))
data


# %%
#
spd=1
nDayWin=96
nDaySkip=-71

data_dt = ecl.field.equatorial_wave.remove_dominant_signals(data, spd,nDayWin,nDaySkip)
data_dt.isel(time =0).plot.contourf(levels = 21)

# %%
#
data_as = ecl.field.equatorial_wave.decompose_symasym(data_dt)
data_as.isel(time = 0).plot.contourf(levels = 21)

# %%
#
psum = ecl.field.equatorial_wave.calc_spectral_coefficients(data_as,spd,nDayWin,nDaySkip)
psum

# %%
#
fig, ax = plt.subplots()

psum.psumanti_r.plot.contourf(ax = ax, levels=21, cmap = 'YlGnBu')
ecl.field.equatorial_wave.draw_wk_anti_analysis()

# %%
#
fig, ax = plt.subplots()

psum.psumsym_r.plot.contourf(ax = ax, levels=21, cmap = 'YlGnBu')
ecl.field.equatorial_wave.draw_wk_sym_analysis()
