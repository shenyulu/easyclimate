# -*- coding: utf-8 -*-
"""Gill Model: Symmetric and Asymmetric Heating
================================================

This example combines symmetric and asymmetric heating components. Comparing
the horizontal circulation and vertical diagnostics shows the response to the
combined forcing.

"""
import numpy as np
from easyclimate.field.simulation.gill import (
    simulate,
    plot_snapshot,
    plot_cross_sections,
)

# %%
# Run the combined symmetric- and asymmetric-heating experiment.
result_symmetric_asymmetric = simulate(case = "gill_sym_asym")
result_symmetric_asymmetric

# %%
# Plot eta and the horizontal wind vector field.
plot_snapshot(
    result_symmetric_asymmetric,
    quiver_kwargs={"scale": 25},
    contourf_kwargs = {"levels": np.linspace(-3, 3, 21)}
)

# %%
# Plot the derived vertical cross-section diagnostics.
plot_cross_sections(
    result_symmetric_asymmetric,
    contourf_uw_kwargs = {"levels": np.linspace(-0.2, 0.2, 21)},
    contourf_vw_kwargs = {"levels": np.linspace(-0.1, 0.1, 21)},
    contourf_u_kwargs = {"levels": np.linspace(-0.2, 0.2, 21)}
)
