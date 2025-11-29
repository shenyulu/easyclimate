# -*- coding: utf-8 -*-
"""
.. _corr_reg_example:

Regression and Correlation Analyses
============================================

Want to dive into the mysteries of climate science using Python?
Today, we’ll walk through a piece of code to explore the impact of **El Niño-Southern Oscillation (ENSO)** on global sea surface temperature (SST)!
This code is simple, beginner-friendly, and perfect for anyone curious about climate science! 😄 Whether you’re learning data analysis or eager to uncover the secrets of climate change, this tutorial will leave you inspired and empowered! 🎉

We’ll use Python’s `easyclimate` library, along with `xarray`, `matplotlib`, and `cartopy`, to step through the code, create stunning visualizations, and reveal how ENSO influences global sea temperatures! In particular, we’ll dive deep into **regression analysis**, using the :py:func:`easyclimate.calc_corr_spatial <easyclimate.calc_corr_spatial>` function, and explain the math behind it! 📐 Ready to get started? Let’s dive in! 🚀


Background: What is Regression Analysis?
++++++++++++++++++++++++++++++++++++++++++++++++

Before we jump into the code, let’s talk about **regression analysis**, the core tool of this tutorial! Regression analysis is a statistical method used to study relationships between variables, particularly how a **dependent variable** (e.g., sea temperature anomalies) changes with an **independent variable** (e.g., ENSO index). 🌡️ In climate science, it helps us explore how ENSO “drives” global sea temperatures!

> Ensure data is preprocessed, e.g., by removing the seasonal cycle to focus on non-seasonal variability.

Steps of Regression Analysis
***************************************

In regression analysis, we assume a **linear relationship** between the dependent variable :math:`y(t)` (e.g., sea temperature anomaly at a grid point) and the independent variable :math:`x(t)` (e.g., Niño 3.4 index), expressed as:

.. math::

    \\hat{y}(t) = a_1 x(t) + a_0

- :math:`\\hat{y}(t)`: Predicted sea temperature anomaly.
- :math:`a_1`: **Regression coefficient** (slope), indicating how much :math:`y(t)` changes per unit change in :math:`x(t)`.
- :math:`a_0`: Intercept, the value of :math:`y(t)` when :math:`x(t) = 0`.

To find the optimal :math:`a_1` and :math:`a_0`, we use the **least squares method** to minimize the sum of squared errors between predicted and actual values:

.. math::

    Q = \\sum_{i=1}^N \\left( \\hat{y}(i) - y(i) \\right)^2 = \\sum_{i=1}^N \\left( (a_1 x(i) + a_0) - y(i) \\right)^2

By taking partial derivatives with respect to :math:`a_0` and :math:`a_1` and setting them to zero, we get:

.. math::

    \\frac { \\mathrm{d} Q } { \\mathrm{d} a _ { 0 } } \\, = \\, 0 \\, = \\, 2 \\sum _ { i \\, = \\, 1 } ^ { N } \\bigl ( a _ { 1 } x ( i ) + a _ { 0 } - y ( i ) \\bigr ) \\, = \\, a _ { 1 } \\sum _ { i \\, = \\, 1 } ^ { N } x ( i ) + a _ { 0 } N - \\sum _ { i \\, = \\, 1 } ^ { N } y ( i )


.. math::

    \\frac { \\mathrm{d} Q } { \\mathrm{d} a _ { 1 } } = \\; 0 \\; = \\; 2 \\sum _ { i \\; = \\; 1 } ^ { N } \\left( a _ { 1 } x ( i ) + a _ { 0 } - y ( i ) \\right) \\cdot x ( i ) \\; = \\; a _ { 1 } \\sum _ { i \\; = \\; 1 } ^ { N } x ( i ) ^ { 2 } + a _ { 0 } \\sum _ { i \\; = \\; 1 } ^ { N } x ( i ) - \\sum _ { i \\; = \\; 1 } ^ { N } x ( i ) y ( i )


Given :math:`\\bar{x} = \\frac{1}{N} \\sum_{i=1}^N x(i)`, we have:

.. math::

    \\bar{y} = a_1 \\bar{x} + a_0

.. math::

    \\overline{xy} = a_1 \\overline{x^2} + a_0 \\bar{x}

Solving these equations yields the slope :math:`a_1` and intercept :math:`a_0`:

.. math::

    a_1 = \\frac{\\overline{xy} - \\bar{x}\\bar{y}}{\\overline{x^2} - \\bar{x}^2},

.. math::

    a_0 = \\bar{y} - a_1 \\bar{x}

This gives us the best-fitting line using the least squares method.
Typically, :math:`a_1` is expressed in terms of deviations from the mean.
To do this, we decompose :math:`x(t)` and :math:`y(t)` into their time mean and deviations:

.. math::

    x(t) = \\bar{x} + x'(t)

.. math::

    y(t) = \\bar{y} + y'(t)

From this, we derive:

.. math::

    \\overline{xy} = \\overline{(\\bar{x} + x')(\\bar{y} + y')} = \\bar{x}\\bar{y} + \\overline{x'y'}

.. math::

    \\overline{x^2} = \\bar{x}^2 + \\overline{x'^2}

Using these and the formula for :math:`a_1`, we get:

.. math::

    a_1 = \\frac{\\overline{x'y'}}{\\overline{x'^2}} = \\frac{\\text{Cov}(x, y)}{\\text{Var}(x)}

- :math:`\\text{Cov}(x, y)`: Covariance of :math:`x(t)` and :math:`y(t)`, measuring how they vary together.
- :math:`\\text{Var}(x)`: Variance of :math:`x(t)`, measuring its own variability.
- :math:`\\bar{x}`, :math:`\\bar{y}`: Means of :math:`x(t)` and :math:`y(t)`.

The slope :math:`a_1` indicates the expected change in :math:`y(t)` per unit change in :math:`x(t)`.

Relationship Between Correlation and Regression Coefficients
******************************************************************************

In regression analysis, we also calculate the **correlation coefficient** :math:`r` to measure the strength of the linear relationship:

.. math::

    r = \\frac{\\overline{x'y'}}{\\sqrt{\\overline{x'^2}} \\sqrt{\\overline{y'^2}}}

.. math::

    r = \\frac{\\text{Cov}(x, y)}{\\sqrt{\\text{Var}(x)} \\sqrt{\\text{Var}(y)}}

- :math:`r` ranges from :math:`[-1, 1]`:
    - :math:`r = 1` indicates perfect positive correlation.
    - :math:`r = -1` indicates perfect negative correlation.
    - :math:`r = 0` indicates no linear relationship.
- The relationship between :math:`a_1` and :math:`r` is:

.. math::

    a_1 = r \\cdot \\frac{\\sqrt{\\text{Var}(y)}}{\\sqrt{\\text{Var}(x)}}

If :math:`x(t)` is standardized (mean = 0, standard deviation = 1), then :math:`\\text{Var}(x) = 1`, and :math:`a_1 = r \\cdot \\sqrt{\\text{Var}(y)}`.

Significance Testing
***************************************

To assess the reliability of the regression results, we use the :math:`t`-statistic to test the significance of the correlation coefficient:

.. math::

    t = \\frac{r \\sqrt{N-2}}{\\sqrt{1-r^2}}

- :math:`N`: Sample size (length of the time series).
- The **p-value** is calculated from the :math:`t`-statistic. A :math:`p`-value < 0.05 (or 0.01) indicates significance at the 95% (or 99%) confidence level.

.. note::

    - Regression analysis is like finding the “best-fit line” to show how sea temperatures “follow” the ENSO index!
    - The correlation coefficient :math:`r` tells you how “close” the relationship is.
    - The :math:`p`-value tells you if it’s “real”! 😉

Preparation: Importing Required Libraries
++++++++++++++++++++++++++++++++++++++++++++++++
"""
import easyclimate as ecl
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# %%
#
# - **easyclimate**: A powerful climate analysis tool with functions like :py:func:`easyclimate.calc_corr_spatial <easyclimate.calc_corr_spatial>` to streamline data processing! 🌍
# - **xarray**: Ideal for handling multi-dimensional data (time, longitude, latitude).
# - **numpy**: The backbone for numerical computations.
# - **matplotlib.pyplot**: For plotting line graphs and maps.
# - **cartopy**: A professional library for geospatial visualizations, perfect for global maps! 🗺️
#
# These libraries are your “climate lab” toolkit—let’s start experimenting! 💪


# %%
# Loading Sea Surface Temperature Data
# ++++++++++++++++++++++++++++++++++++++++++++++++
#
# - ``ecl.open_tutorial_dataset("sst_mnmean_oisst")``: Loads a global sea surface temperature (SST) dataset with time, longitude, and latitude dimensions.
# - ``.sst``: Extracts the sea temperature variable (units: °C).
# - ``sst_data``: An :py:class:`xarray.DataArray<xarray.DataArray>`, like a 3D “data cube”! 🧊
#
# .. note::
#
#     SST stands for sea surface temperature. ENSO can make some regions warmer or cooler, affecting global weather! 🌦️

sst_data = ecl.open_tutorial_dataset("sst_mnmean_oisst").sst
sst_data

# %%
# Removing the Seasonal Cycle to Extract Anomalies
# ++++++++++++++++++++++++++++++++++++++++++++++++
#
# - :py:func:`easyclimate.remove_seasonal_cycle_mean <easyclimate.remove_seasonal_cycle_mean>`: Removes seasonal variations (e.g., warmer summers, cooler winters) to get **sea temperature anomalies**, focusing on non-seasonal signals like ENSO.
# - ``sst_data_anomaly``: A new dataset containing only “unusual” temperature changes.
#
# .. note::
#
#     - Removing the seasonal cycle is like filtering out “background noise” to hear ENSO’s “solo performance”! 🎶
#     - Heads up! ⚠️ If you don’t remove the seasonal cycle before correlation or regression analysis, most climate variables will show high correlations with distant regions—simply because they’re all influenced by the annual solar radiation cycle! 🌍☀️ Pretty cool, right?
#
sst_data_anormaly = ecl.remove_seasonal_cycle_mean(sst_data)
sst_data_anormaly

# %%
# Calculating the Niño 3.4 Index
# ++++++++++++++++++++++++++++++++++++++++++++++++
#
# 1. Calculating the Niño 3.4 Index:
#
#    - :py:func:`easyclimate.field.air_sea_interaction.calc_index_nino34 <easyclimate.field.air_sea_interaction.calc_index_nino34>`: Computes the Niño 3.4 index from sea temperature anomalies in the eastern tropical Pacific (5°S–5°N, 170°W–120°W), measuring ENSO strength.
#    - ``running_mean = 0``: Uses raw data without smoothing.
#
# 2. Plotting the Line Graph:
#
#    - ``plt.figure(figsize=(10, 3))``: Creates a canvas 10 units wide and 3 units high.
#    - :py:func:`easyclimate.plot.line_plot_with_threshold <easyclimate.plot.line_plot_with_threshold>`: Plots the Niño 3.4 index time series.
#    - :py:func:`plt.title <matplotlib.pyplot.title>`: Sets the title to “Niño 3.4 index”.
#
# .. note::
#
#     The Niño 3.4 index is like a “weather gauge” for ENSO—positive values indicate El Niño (warm), and negative values indicate La Niña (cool). 📈

nino34_index = ecl.field.air_sea_interaction.calc_index_nino34(
    sst_data,
    running_mean = 0
)

plt.figure(figsize=(10, 3))
ecl.plot.line_plot_with_threshold(nino34_index)
plt.title("Niño 3.4 index")

# %%
# Standardizing the Niño 3.4 Index
# ++++++++++++++++++++++++++++++++++++++++++++++++
#
# For the next steps, if we perform linear regression between the SST anomaly time series :math:`y(t)` at each global grid point and the standardized ENSO index :math:`x(t)`, we calculate the regression coefficient :math:`a_1`:
#

# %%
# .. math::
#
#     a_1 = \frac{\text{Cov} (x,y)}{\text{Var}(x)}

# %%
# Since :math:`x(t)` is standardized, :math:`\mathrm{Var}(x) = 1`, so :math:`a_1` directly reflects the covariance, with units of :math:`\frac{\text{ K (local SST change) }}{ \sigma \text{(ENSO index)} } = \mathrm{K}`.
#
# - :py:func:`easyclimate.normalized.timeseries_normalize_zscore <easyclimate.normalized.timeseries_normalize_zscore>`: Standardizes the Niño 3.4 index (mean = 0, standard deviation = 1) using:
#
# .. math::
#
#     z = \frac{x - \bar{x}}{\sqrt{\mathrm{Var}(x)}}
#
# - Plots the standardized time series with the title “Normalized Niño 3.4 index”.
#
# .. note::
#
#     Standardization puts the index in “standard units,” making it easier to compare and use in regression analysis! 📏
#

nino34_index_normalized = ecl.normalized.timeseries_normalize_zscore(nino34_index)

plt.figure(figsize=(10, 3))
ecl.plot.line_plot_with_threshold(nino34_index_normalized)
plt.title("Normalized Niño 3.4 index")

# %%
# Calculating Correlation Between SST Anomalies and Niño 3.4 Index
# ++++++++++++++++++++++++++++++++++++++++++++++++
#
# - ``ecl.calc_corr_spatial``: Performs regression analysis between global SST anomalies (`sst_data_anomaly`) and the standardized Niño 3.4 index (`nino34_index_normalized`).
# - Output: ``sst_reg_nino34_result``, containing:
#
#   - ``reg_coeff``: Regression coefficient, in units of °C per standard deviation of the index.
#   - ``corr``: Correlation coefficient, ranging from :math:`[-1, 1]`.
#   - ``pvalue``: Significance p-value.
#
# .. note::
#
#     This step is like giving global sea temperatures a “check-up” to see which regions are closely tied to ENSO! 🔍
sst_reg_nino34_result = ecl.calc_corr_spatial(
    sst_data_anormaly,
    x = nino34_index_normalized
)
sst_reg_nino34_result

# %%
# Regression Coefficient Map
# ++++++++++++++++++++++++++++++++++++++++++++++++
#
# - **Base Map**: ``quick_draw_spatial_basemap`` creates a global map centered at 205° longitude.
# - **Regression Coefficient Plot**: ``reg_coeff.plot.contourf`` visualizes the amplitude of SST changes (°C) per standard deviation change in the Niño 3.4 index.
# - **Significance Markers**: ``draw_significant_area_contourf`` highlights areas with :math:`p < 0.01`.
# - **Land**: Adds gray landmasses, with the title indicating regression results.
#
# The map shows the global SST anomaly regressed onto the standardized ENSO index, displayed as a spatial distribution of SST anomalies in units of K.
fig, ax = ecl.plot.quick_draw_spatial_basemap(central_longitude=205)

sst_reg_nino34_result.reg_coeff.plot.contourf(
    ax = ax,
    levels = np.linspace(-1.5, 1.5, 21),
    transform = ccrs.PlateCarree(),
    cbar_kwargs={"location": "bottom", "aspect": 30, "label": "℃"},
)

ecl.plot.draw_significant_area_contourf(
    sst_reg_nino34_result.pvalue,
    thresh = 0.01,
    ax = ax,
    transform = ccrs.PlateCarree(),
)

ax.add_feature(cfeature.LAND, facecolor = '#DDDDDD', zorder = 1)
ax.set_title("Niño 3.4 index Reg SST anormaly")

# %%
#
# - In the eastern tropical Pacific (20°N–20°S, 160°E to South American coast), a one-standard-deviation increase in the ENSO warm phase (or cold tongue) index corresponds to significant positive SST anomalies, exceeding 1 K in the eastern tropical Pacific.
# - In the western tropical Pacific, SST anomalies are near 0 K, indicating minimal ENSO influence.
# - In the North Pacific, the ENSO warm phase corresponds to negative SST anomalies of about -0.2 K, indicating cooler temperatures.
# - The regression coefficients largely reflect the spatial pattern of global SST changes during the ENSO warm phase.
# - The regression coefficients show the magnitude of change but do not indicate statistical significance. Significance can be assessed by plotting confidence levels (e.g., 95%).
#
# .. note::
#
#     - Linear Assumption
#
#       Regression analysis is a **linear** method 📈. In the example above, the cold phase of the ENSO cycle (La Niña) is assumed to exhibit a pattern exactly opposite to the anomalies shown in the regression results. It’s like playing a “mirror symmetry” game!
#
#       To capture more complex **non-linear** relationships 🤹‍♂️ between ENSO and global SST anomalies, scientists often perform **composite analysis** (i.e., averaging separately for warm (El Niño) and cold (La Niña) events) to reveal their distinct characteristics.
#
#     - Causality
#
#       Regression analysis only reveals statistical associations, not causation 🔗 ≠ 🧪.
#
#       In the example, North Atlantic SST anomalies may be a remote “response” or “result” of tropical climate anomalies, as determined by numerical simulations. However, regression results alone cannot determine cause and effect—interpret with caution! ⚠️
#
#     - Effective Sample Size
#
#       The significance of correlations depends on the “effective sample size” 🧠, which is almost always smaller than the number of data points used in the analysis.
#
#       For example 🌰: A 50-year record of monthly mean Atlantic basin SST anomalies might seem to have :math:`50 \\cdot 12 = 600` time points, but due to the strong “memory” (persistence) of regionally averaged SSTs, the number of truly independent samples may be much smaller!
#
#       So, when interpreting correlations based on limited sample sizes, stay cautious and critical—numbers can be deceptive! 🕵️‍♀️
#
# .. note::
#
#     This map shows where sea temperatures rise (red) or fall (blue) when ENSO “acts up”! 🌡️

# %%
# Correlation Coefficient Map
# ++++++++++++++++++++++++++++++++++++++++++++++++
#
# - Similar to above, but plots the **correlation coefficient** map, showing the strength of the relationship between SST and the ENSO index (-1 to 1).
# - Title changed to “corr Niño 3.4 index & SST anomaly”.
#
# .. note::
#
#     This map highlights the “connection strength” between ENSO and sea temperatures—red means they “move together,” blue means they “move opposite”! 😎
fig, ax = ecl.plot.quick_draw_spatial_basemap(central_longitude=205)

sst_reg_nino34_result.corr.plot.contourf(
    ax = ax,
    levels = np.linspace(-1, 1, 11),
    transform = ccrs.PlateCarree(),
    cbar_kwargs={"location": "bottom", "aspect": 30, "label": ""},
)

ecl.plot.draw_significant_area_contourf(
    sst_reg_nino34_result.pvalue,
    thresh = 0.01,
    ax = ax,
    transform = ccrs.PlateCarree(),
)

ax.add_feature(cfeature.LAND, facecolor = '#DDDDDD', zorder = 1)
ax.set_title("corr Niño 3.4 index & SST anormaly")

# %%
# Summary: What Did We Learn?
# ++++++++++++++++++++++++++++++++++++++++++++++++
#
# Through this code, we embarked on a climate science adventure! 🌍 We:
#
# 1. Loaded and preprocessed SST data to extract anomalies.
# 2. Calculated and standardized the Niño 3.4 index, plotting its time series.
# 3. Used :py:func:`easyclimate.calc_corr_spatial <easyclimate.calc_corr_spatial>` to analyze the relationship between global SSTs and ENSO, obtaining regression coefficients, correlation coefficients, and p-values.
# 4. Created two maps showing the amplitude and strength of ENSO’s impact on SSTs.
#
# Pro Tips
# ***************************************
#
# - The regression coefficient :math:`a_1` shows the “actual impact” of ENSO changes on SST anomalies (°C), while the correlation coefficient :math:`r` shows the “strength of the connection.”
# - :py:func:`easyclimate.calc_corr_spatial <easyclimate.calc_corr_spatial>` efficiently handles multi-dimensional data, automatically skipping NaN values, making it ideal for large-scale climate analysis.
# - Areas with :math:`p < 0.01/0.05` are “rock-solid” results worth focusing on! 🔥
#
# Beginner Tips
# ***************************************
#
# - ENSO is like the ocean’s “mood,” affecting global climate, and this code helps you “hear” its story! 😉
# - Want to go deeper? Check the `easyclimate` documentation or try other climate indices!
#
# Run the code and explore the wonders of climate science! 🚀
