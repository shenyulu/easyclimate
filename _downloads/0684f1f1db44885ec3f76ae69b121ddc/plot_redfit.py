# -*- coding: utf-8 -*-
"""
.. _redfit_example:

Estimate Red-noise Spectrum
============================================

Before proceeding with all the steps, first import some necessary libraries and packages
"""
import pooch
import xarray as xr
import numpy as np
import pandas as pd
import easyclimate as ecl
import matplotlib.pyplot as plt
from matplotlib import ticker

# %%
# You can download following datasets here:
#
# .. code-block:: python
#
#       pooch.retrieve(
#           "https://psl.noaa.gov/data/correlation/nina34.anom.data",
#           known_hash=None,
#           fname = "nina34.anom.data",
#           path = ".",
#       )
#
# Now we begin to read the txt file

# Read data from txt file
file_path = "nina34.anom.data"
with open(file_path, 'r') as file:
    lines = file.readlines()

# Parse data
years = []
values = []
for line in lines:
    # Skip empty lines or comment lines
    if not line.strip() or line.startswith('#'):
        continue
    parts = line.split()

    # Check if the first column is a year (integer)
    try:
        year = int(parts[0])  # Attempt to convert the first column to an integer
    except ValueError:
        print(f"Skipping invalid line (first column is not a year): {line.strip()}")
        continue

    # Check if there are 12 monthly values
    if len(parts[1:]) != 12:
        print(f"Skipping invalid line (missing 12 monthly values): {line.strip()}")
        continue

    # Convert the 12 monthly values to floats
    try:
        monthly_values = list(map(float, parts[1:]))
    except ValueError:
        print(f"Skipping invalid line (contains non-numeric data): {line.strip()}")
        continue

    years.append(year)
    values.append(monthly_values)

# Create time index
time = pd.date_range(start=f'{years[0]}-01', periods=len(years) * 12, freq='ME')

# Flatten the data into a 1D array
flat_values = np.array(values).flatten()

# Create xarray.DataArray
nino34 = xr.DataArray(
    flat_values,
    dims=['time'],
    coords={'time': time},
    attrs={'description': 'Nino 3.4 Index', 'units': 'Celsius'}
)

# Replace -99.99 with NaN
nino34 = nino34.where(nino34 != -99.99)

# %%
# Filter needed time range
nino34 = nino34.isel(time = slice(24,-12))
nino34

# %%
# Calculate red noise
result_redfit = ecl.filter.calc_redfit(nino34)
result_redfit = result_redfit.assign_coords({"period_month": (result_redfit.period)/12})
result_redfit

# %%
# Draw the red noise graph
fig, ax = plt.subplots()

result_redfit.gxx.plot(ax = ax, x = 'period_month', color = 'black')
result_redfit.chi2_95.plot(ax = ax, ls = '--', x = 'period_month', label = '95% CI', color = 'red')
result_redfit.chi2_90.plot(ax = ax, ls = '--', x = 'period_month', label = '90% CI', color = 'blue')
ax.legend()

ax.set_xscale('log', base = 2, subs = None)
ax.xaxis.set_major_formatter(ticker.ScalarFormatter())

ax.set_xlim(2, 16)

ax.set_xlabel('Period (years)')
ax.set_ylabel('Spectral Amplitude')
