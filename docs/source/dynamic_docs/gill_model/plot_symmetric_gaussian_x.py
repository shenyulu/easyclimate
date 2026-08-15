# -*- coding: utf-8 -*-
"""Gill Model: Symmetric Heating with Gaussian Zonal Structure
==============================================================

This experiment applies equatorially symmetric heating with a Gaussian zonal
envelope. The plots reveal how the localized forcing shapes the horizontal
response and its vertical diagnostics.

"""
import numpy as np
from easyclimate.field.simulation.gill import (
    simulate,
    plot_snapshot,
    plot_cross_sections,
)

# %%
# Run the symmetric-heating experiment with a Gaussian zonal envelope.
result_symmetric_gaussian_x = simulate(case = "gill_sym_gaussian_x")
result_symmetric_gaussian_x

# %%
# Plot eta and the horizontal wind vector field.
plot_snapshot(
    result_symmetric_gaussian_x,
    quiver_kwargs={"scale": 10},
    contourf_kwargs = {"levels": np.linspace(-2, 2, 21)}
)

# %%
# Plot the derived vertical cross-section diagnostics.
plot_cross_sections(
    result_symmetric_gaussian_x,
    contourf_uw_kwargs = {"levels": np.linspace(-0.2, 0.2, 21)},
    contourf_vw_kwargs = {"levels": np.linspace(-0.02, 0.02, 21)},
    contourf_u_kwargs = {"levels": np.linspace(-0.2, 0.2, 21)}
)
