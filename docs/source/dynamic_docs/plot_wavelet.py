# -*- coding: utf-8 -*-
"""
.. _wavelet_example:

Wavelets
======================

.. seealso::

    - Torrence, C., & Compo, G. P. (1998). A Practical Guide to Wavelet Analysis. Bulletin of the American Meteorological Society, 79(1), 61-78. https://doi.org/10.1175/1520-0477(1998)079<0061:APGTWA>2.0.CO;2
    - Torrence, C., & Webster, P. J. (1999). Interdecadal Changes in the ENSO–Monsoon System. Journal of Climate, 12(8), 2679-2690. https://doi.org/10.1175/1520-0442(1999)012<2679:ICITEM>2.0.CO;2
    - Grinsted, A., Moore, J. C., and Jevrejeva, S.: Application of the cross wavelet transform and wavelet coherence to geophysical time series, Nonlin. Processes Geophys., 11, 561–566, https://doi.org/10.5194/npg-11-561-2004, 2004.

Before proceeding with all the steps, first import some necessary libraries and packages
"""
import xarray as xr
import matplotlib.pyplot as plt
import easyclimate as ecl

# %%
# The

data_nino3 = xr.open_dataset('test_input_nino3_wavelet.nc')['nino3']
result_data = ecl.filter.calc_timeseries_wavelet_transform(data_nino3, dt = 0.25)
result_data

# %%
# The

fig, ax = plt.subplots()
ecl.filter.draw_global_wavelet_spectrum(result_data, ax = ax)

# %%
# The

fig, ax = plt.subplots()
ecl.filter.draw_wavelet_transform(result_data, ax = ax)
