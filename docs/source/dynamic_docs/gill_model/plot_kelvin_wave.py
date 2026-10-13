# -*- coding: utf-8 -*-
"""Gill Model: Kelvin Wave
==============================================

This example simulates the Kelvin-wave experiment. The horizontal snapshot
shows its equatorially trapped circulation, while the cross sections show the
associated vertical structure.

"""
import numpy as np
from easyclimate.field.simulation.gill import (
    simulate,
    plot_snapshot,
    plot_cross_sections,
)

# %%
# Run the Kelvin-wave experiment.
result_kelvin_wave = simulate(case = "kelvin_wave")
result_kelvin_wave

# %%
# Plot eta and the horizontal wind vector field.
plot_snapshot(
    result_kelvin_wave,
    quiver_kwargs={"scale": 8},
    contourf_kwargs = {"levels": np.linspace(-1, 1, 21)}
)

# %%
# Plot the derived vertical cross-section diagnostics.
plot_cross_sections(
    result_kelvin_wave,
    contourf_uw_kwargs = {"levels": np.linspace(-0.2, 0.2, 21)},
    contourf_vw_kwargs = {"levels": np.linspace(-0.1, 0.1, 21)},
    contourf_u_kwargs = {"levels": np.linspace(-0.1, 0.1, 21)}
)
