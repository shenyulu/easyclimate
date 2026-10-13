# -*- coding: utf-8 -*-
"""Gill Model: Gaussian Initial Perturbation
==============================================

This example evolves a Gaussian-shaped initial perturbation without damping.
The horizontal and vertical plots show the transient adjustment of the
idealized circulation.

"""
import numpy as np
from easyclimate.field.simulation.gill import (
    simulate,
    plot_snapshot,
    plot_cross_sections,
)

# %%
# Run the Gaussian initial-perturbation experiment.
result_gaussian_initial = simulate(case = "gauss_pert")
result_gaussian_initial

# %%
# Plot eta and the horizontal wind vector field.
plot_snapshot(
    result_gaussian_initial,
    quiver_kwargs={"scale": 2.5},
    contourf_kwargs = {"levels": np.linspace(-0.2, 0.2, 21)}
)

# %%
# Plot the derived vertical cross-section diagnostics.
plot_cross_sections(
    result_gaussian_initial,
    contourf_uw_kwargs = {"levels": np.linspace(-0.05, 0.05, 21)},
    contourf_vw_kwargs = {"levels": np.linspace(-0.05, 0.05, 21)},
    contourf_u_kwargs = {"levels": np.linspace(-0.05, 0.05, 21)}
)
