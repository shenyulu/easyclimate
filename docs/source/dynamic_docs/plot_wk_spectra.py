# -*- coding: utf-8 -*-
"""
.. _wk_spectra_example:

Wheeler-Kiladis Space-Time Spectra
============================================

The **Wheeler-Kiladis (WK) Space-Time Spectra** is a diagnostic tool used in meteorology to analyze tropical atmospheric waves by decomposing their variability in both **wavenumber (space) and frequency (time)** domains. It extends traditional Fourier analysis by applying a **two-dimensional (2D) spectral decomposition** to isolate and characterize different wave modes (e.g., Kelvin waves, Rossby waves, and mixed Rossby-gravity waves) based on their dispersion relations.

Key Features
--------------------------------------------

- Uses **symmetrical (eastward/westward) and antisymmetrical (north-south) Fourier transforms** to separate wave types.
- Compares observed spectra with **theoretical dispersion curves** of shallow-water waves to identify dominant modes.
- Helps distinguish **convectively coupled waves** (linked to tropical rainfall) from uncoupled waves.

Applications in Meteorology
--------------------------------------------

1. Tropical Wave Analysis: Identifies and tracks **Kelvin waves, equatorial Rossby waves, and Madden-Julian Oscillation (MJO)** signals.
2. Convective Coupling Studies: Examines how waves interact with tropical convection (e.g., in monsoon systems).
3. Model Validation: Evaluates whether climate models correctly simulate tropical wave dynamics.
4. Extreme Weather Prediction: Helps understand precursors to tropical cyclogenesis and organized convection.

The WK method is particularly useful for studying **large-scale tropical variability** and remains a fundamental tool in **tropical meteorology and climate research**.

Before proceeding with all the steps, first import some necessary libraries and packages
"""
import xarray as xr
import matplotlib.pyplot as plt
import easyclimate as ecl


# %%
# The example here is to avoid longer calculations, thus we open the pre-processed result data directly.
#
# .. tip::
#
#   You can download following datasets here: :download:`Download olr_smooth_data.nc <https://raw.githubusercontent.com/shenyulu/easyclimate/refs/heads/main/docs/source/dynamic_docs/olr_smooth_data.nc>`
#

data = xr.open_dataset('olr_smooth_data.nc')['olr'].sel(lat = slice(-15, 15))
data


# %%
# Setting the basic parameters and removing the dominant signal
spd=1
nDayWin=96
nDaySkip=-71

data_dt = ecl.field.equatorial_wave.remove_dominant_signals(data, spd,nDayWin,nDaySkip)
data_dt.isel(time =0).plot.contourf(levels = 21)

# %%
# Separation of symmetric and asymmetric parts
data_as = ecl.field.equatorial_wave.decompose_symasym(data_dt)
data_as.isel(time = 0).plot.contourf(levels = 21)

# %%
# Calculation of spectral coefficients
psum = ecl.field.equatorial_wave.calc_spectral_coefficients(data_as,spd,nDayWin,nDaySkip)
psum

# %%
# Drawing asymmetric parts
fig, ax = plt.subplots()

psum.psumanti_r.plot.contourf(ax = ax, levels=21, cmap = 'YlGnBu')
ecl.field.equatorial_wave.draw_wk_anti_analysis()

# %%
# And drawing symmetric parts

# sphinx_gallery_thumbnail_number = -1
fig, ax = plt.subplots()

psum.psumsym_r.plot.contourf(ax = ax, levels=21, cmap = 'YlGnBu')
ecl.field.equatorial_wave.draw_wk_sym_analysis()
