# -*- coding: utf-8 -*-
"""
Visualization for Spherical Harmonics
==============================================


This example transforms a gridded field to packed triangular spectral
coefficients, visualizes the spectral coefficients, and transforms the field
back to latitude-longitude space. The standard backend and Rust-backed
functions are compared with the same input data.

"""
import numpy as np
import matplotlib.pyplot as plt
import easyclimate as ecl

# %%
# Open a tutorial specific-humidity field. The spectral transform expects latitude and longitude dimensions on a regular or Gaussian global grid.
q_grid_data = ecl.open_tutorial_dataset("shum_2022_day5").shum
q_grid_data

# %%
# Transform grid-space data to packed triangular spherical-harmonic coefficients. The output has a ``spec_dim`` dimension whose length is determined by the triangular truncation.
q_spec_data = ecl.spec.transfer_grid2spectral_transform(
    q_grid_data,
    grid_data_type = "regular"
)
q_spec_rs_data = ecl.spec.transfer_grid2spectral_transform_rs(
    q_grid_data,
    grid_data_type = "regular"
)
q_spec_data

# %%
# The plotting helpers expect a one-dimensional spectral coefficient array, so select one time and level before visualizing the spectrum.
q_spec_data_sample = q_spec_data.isel(time = 0).sel(level = 850)
q_spec_data_sample

# %%
# Plot the total-wavenumber power spectrum, computed as the sum of squared coefficient magnitude over zonal wavenumber.
ecl.spec.plot_spectral_power(q_spec_data_sample)

# %%
# Plot the packed coefficients as an ``m`` by ``n`` triangular matrix. By default this displays the base-10 logarithm of coefficient magnitude.
ecl.spec.plot_spectral_matrix(q_spec_data_sample)

# %%
# Plot the real and imaginary parts of the complex spectral coefficients in two side-by-side panels.
ecl.spec.plot_spectral_real_imag(q_spec_data_sample)

# %%
# Transform the spectral coefficients back to grid space. The output grid size and grid type are supplied explicitly, matching the original regular grid.
q_grid_spectrans = ecl.spec.transfer_spectral_transform2grid(
    q_spec_data,
    nlon = 144, nlat = 73,
    grid_data_type = "regular"
)
q_grid_rs_spectrans = ecl.spec.transfer_spectral_transform2grid_rs(
    q_spec_rs_data,
    nlon = 144, nlat = 73,
    grid_data_type = "regular"
)
q_grid_spectrans

# %%
# Compare the original field, the grid field reconstructed after one spectral transform cycle, and the backend difference between Fortran and Rust transforms.
draw_q_grid = q_grid_data.isel(time = 2).sel(level = 850)
draw_q_grid_spectrans = q_grid_spectrans.isel(time = 2).sel(level = 850)
draw_diff_grid_minus_grid_spectrans = draw_q_grid - draw_q_grid_spectrans
draw_diff_rust = (q_grid_spectrans - q_grid_rs_spectrans).isel(time = 2).sel(level = 850)

# ---------------------------
fig, ax = plt.subplots(2, 2, figsize = (10, 10))

axi = ax[0, 0]
draw_q_grid.plot.contourf(
    ax = axi,
    levels= np.linspace(0, 0.015, 11),
    cmap = "plasma",
    cbar_kwargs = {'location': 'bottom', 'aspect': 60},
)
axi.set_title("Orignal [O]")

axi = ax[0, 1]
draw_q_grid_spectrans.plot.contourf(
    ax = axi,
    levels= np.linspace(0, 0.015, 11),
    cmap = "plasma",
    cbar_kwargs = {'location': 'bottom', 'aspect': 60},
)
axi.set_title("Single spectral transformation [ST]")

axi = ax[1, 0]
draw_diff_grid_minus_grid_spectrans.plot(
    ax = axi,
    vmax = 5e-4,
    cmap = "RdBu_r",
    cbar_kwargs = {'location': 'bottom', 'aspect': 60},
)
axi.set_title("Diff: [O] - [ST]")

axi = ax[1, 1]
draw_diff_rust.plot(
    ax = axi,
    vmax = 1e-8,
    cmap = "RdBu_r",
    cbar_kwargs = {'location': 'bottom', 'aspect': 60},
)
axi.set_title("Diff: (Fortran - Rust) [ST]")
