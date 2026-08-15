# -*- coding: utf-8 -*-
"""Gill Model: Delta Heating at y = 1
==============================================

This example applies an idealized delta-function heating source at ``y = 1``.
The horizontal and vertical diagnostics highlight the circulation forced by
this narrowly localized meridional heating.

"""
import numpy as np
from easyclimate.field.simulation.gill import (
    simulate,
    plot_snapshot,
    plot_cross_sections,
)

# %%
# Run the delta-heating experiment at y = 1.
result_delta_y1 = simulate(case = "gill_delta_y1")
result_delta_y1

# %%
# Plot eta and the horizontal wind vector field.
plot_snapshot(
    result_delta_y1,
    quiver_kwargs={"scale": 25},
    contourf_kwargs = {"levels": np.linspace(-5, 5, 21)}
)

# %%
# Plot the derived vertical cross-section diagnostics.
plot_cross_sections(
    result_delta_y1,
    contourf_uw_kwargs = {"levels": np.linspace(-0.5, 0.5, 21)},
    contourf_vw_kwargs = {"levels": np.linspace(-0.5, 0.5, 21)},
    contourf_u_kwargs = {"levels": np.linspace(-2, 2, 21)}
)
