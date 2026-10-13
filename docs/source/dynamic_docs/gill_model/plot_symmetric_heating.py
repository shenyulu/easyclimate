# -*- coding: utf-8 -*-
"""Gill Model: Symmetric Heating
==============================================

This example simulates the steady linear response to heating that is symmetric
about the equator. The horizontal snapshot shows the geopotential-height
perturbation and horizontal circulation; the cross sections diagnose the
corresponding vertical structure.

"""
import numpy as np
from easyclimate.field.simulation.gill import (
    simulate,
    plot_snapshot,
    plot_cross_sections,
)

# %%
# Run the symmetric Gill-heating experiment. The returned xarray Dataset holds
# the final horizontal fields and the coordinates needed by the plotting helpers.
result_symmetric_heating = simulate(case = "gill_sym")
result_symmetric_heating

# %%
# Plot eta as filled contours with the horizontal wind vector field.
plot_snapshot(
    result_symmetric_heating,
    quiver_kwargs={"scale": 30},
)

# %%
# Plot the derived u-w and v-w streamfunctions and the zonal-mean zonal wind.
plot_cross_sections(
    result_symmetric_heating,
    contourf_uw_kwargs = {"levels": np.linspace(-0.2, 0.2, 21)},
    contourf_vw_kwargs = {"levels": np.linspace(-0.02, 0.02, 21)},
    contourf_u_kwargs = {"levels": np.linspace(-0.2, 0.2, 21)}
)
