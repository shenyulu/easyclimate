# -*- coding: utf-8 -*-
"""Gill Model: Forced Stationary Motion
==============================================

This example reproduces the Matsuno forced stationary-motion experiment. The
plots show the horizontal response to the forcing and its diagnosed vertical
cross sections.

"""
import numpy as np
from easyclimate.field.simulation.gill import (
    simulate,
    plot_snapshot,
    plot_cross_sections,
)

# %%
# Run the forced stationary-motion experiment.
result_forced_stationary = simulate(case = "matsuno_heating")
result_forced_stationary

# %%
# Plot eta and the horizontal wind vector field.
plot_snapshot(
    result_forced_stationary,
    quiver_kwargs={"scale": 10},
    contourf_kwargs = {"levels": np.linspace(-2, 2, 21)}
)

# %%
# Plot the derived vertical cross-section diagnostics.
plot_cross_sections(
    result_forced_stationary,
    contourf_uw_kwargs = {"levels": np.linspace(-0.5, 0.5, 21)},
    contourf_vw_kwargs = {"levels": np.linspace(-0.05, 0.05, 21)},
    contourf_u_kwargs = {"levels": np.linspace(-0.05, 0.05, 21)}
)
