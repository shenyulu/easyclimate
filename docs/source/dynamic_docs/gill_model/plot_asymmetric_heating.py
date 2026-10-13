# -*- coding: utf-8 -*-
"""Gill Model: Asymmetric Heating
==============================================

This example simulates the response to heating that is asymmetric about the
equator. It uses the same horizontal and vertical diagnostics as the symmetric
case, making the meridional asymmetry of the circulation easy to inspect.

"""
import numpy as np
from easyclimate.field.simulation.gill import (
    simulate,
    plot_snapshot,
    plot_cross_sections,
)

# %%
# Run the asymmetric Gill-heating experiment.
result_asymmetric_heating = simulate(case = "gill_asym")
result_asymmetric_heating

# %%
# Plot eta and the horizontal wind vector field.
plot_snapshot(
    result_asymmetric_heating,
    quiver_kwargs={"scale": 25},
    contourf_kwargs = {"levels": np.linspace(-2, 2, 21)}
)

# %%
# Plot the derived vertical cross-section diagnostics.
plot_cross_sections(
    result_asymmetric_heating,
    contourf_uw_kwargs = {"levels": np.linspace(-0.01, 0.01, 21)},
    contourf_vw_kwargs = {"levels": np.linspace(-0.1, 0.1, 21)},
    contourf_u_kwargs = {"levels": np.linspace(-0.4, 0.4, 21)}
)
