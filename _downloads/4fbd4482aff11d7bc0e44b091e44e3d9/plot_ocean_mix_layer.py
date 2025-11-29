# -*- coding: utf-8 -*-
"""
Ocean Surface Mixed Layers Analysis
=========================================================================================================

Before proceeding with all the steps, first import some necessary libraries and packages
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import easyclimate as ecl

# %%
# Preprocessed data
#
#
# .. code-block:: python
#
#     temper_data = xr.open_dataset('temp_soda3.4.2_mn_ocean_reg_2020_EN4.nc', chunks="auto").temp.rename({'st_ocean': 'depth'})
#     slt_data = xr.open_dataset('salt_soda3.4.2_mn_ocean_reg_2020_EN4.nc', chunks="auto").salt.rename({'st_ocean': 'depth'})
#
#     u_data = xr.open_dataset('u_soda3.4.2_mn_ocean_reg_2020_EN4.nc', chunks="auto").u.rename({'st_ocean': 'depth'})
#     v_data = xr.open_dataset('v_soda3.4.2_mn_ocean_reg_2020_EN4.nc', chunks="auto").v.rename({'st_ocean': 'depth'})
#     net_heating_data = xr.open_dataset('net_heating_soda3.4.2_mn_ocean_reg_2020_EN4.nc', chunks="auto").net_heating
#
#     wt_data = xr.open_dataset('wt_soda3.4.2_mn_ocean_reg_2020_EN4.nc', chunks="auto").wt.rename({'sw_ocean': 'depth'})

# %%
# The following data is the depth data of the mixed layer output by the oceanic models (we will compare it with it later).
mld_data = ecl.open_tutorial_dataset('mlp_soda3_4_2_mn_ocean_reg_2020_EN4').mlp

# %%
#
# .. tip::
#
#   You can download following datasets here:
#
#   - :download:`Download temp_soda3.4.2_mn_ocean_reg_2020_EN4.nc <https://huggingface.co/datasets/shenyulu/easyclimate/resolve/main/tutorial_data/temp_soda3.4.2_mn_ocean_reg_2020_EN4.nc>`
#   - :download:`Download salt_soda3.4.2_mn_ocean_reg_2020_EN4.nc <https://huggingface.co/datasets/shenyulu/easyclimate/resolve/main/tutorial_data/salt_soda3.4.2_mn_ocean_reg_2020_EN4.nc>`
#   - :download:`Download u_soda3.4.2_mn_ocean_reg_2020_EN4.nc <https://huggingface.co/datasets/shenyulu/easyclimate/resolve/main/tutorial_data/u_soda3.4.2_mn_ocean_reg_2020_EN4.nc>`
#   - :download:`Download v_soda3.4.2_mn_ocean_reg_2020_EN4.nc <https://huggingface.co/datasets/shenyulu/easyclimate/resolve/main/tutorial_data/v_soda3.4.2_mn_ocean_reg_2020_EN4.nc>`
#   - :download:`Download net_heating_soda3.4.2_mn_ocean_reg_2020_EN4.nc <https://huggingface.co/datasets/shenyulu/easyclimate/resolve/main/tutorial_data/net_heating_soda3.4.2_mn_ocean_reg_2020_EN4.nc>`
#   - :download:`Download wt_soda3.4.2_mn_ocean_reg_2020_EN4.nc <https://huggingface.co/datasets/shenyulu/easyclimate/resolve/main/tutorial_data/wt_soda3.4.2_mn_ocean_reg_2020_EN4.nc>`
#   - :download:`Download mlp_soda3.4.2_mn_ocean_reg_2020_EN4.nc <https://huggingface.co/datasets/shenyulu/easyclimate/resolve/main/tutorial_data/mlp_soda3.4.2_mn_ocean_reg_2020_EN4.nc>`
#
# .. warning::
#
#     - Here we are using only the `SODA <https://www2.atmos.umd.edu/~ocean/>`__ 3.4.2 reanalysis data during 2024; the actual analysis will need to be analyzed using multiple years of data and removing the seasonal cycle by :py:func:`easyclimate.remove_seasonal_cycle_mean <easyclimate.remove_seasonal_cycle_mean>`.
#     - **Citation**: Carton, J. A., Chepurin, G. A., & Chen, L. (2018). SODA3: A New Ocean Climate Reanalysis. Journal of Climate, 31(17), 6967-6983. https://doi.org/10.1175/JCLI-D-18-0149.1
#
#
# Mix-Layer Depth
# -------------------------------------
#
# Here, we use :py:func:`easyclimate.field.ocean.calc_mixed_layer_depth <easyclimate.field.ocean.calc_mixed_layer_depth>` to calculate the depth of mix-layer (MLD).
#
# .. code-block:: python
#
#     mixed_layer_depth = ecl.field.ocean.calc_mixed_layer_depth(
#         seawater_temperature_data = temper_data,
#         seawater_practical_salinity_data = slt_data,
#     ).to_netcdf("sample_mixed_layer_depth.nc")
#
#
# .. seealso::
#
#     - https://github.com/pyoceans/oceans
#     - https://pyoceans.github.io/python-oceans/ocfis.html#oceans.ocfis.mld
#
# Next, we open the dataset containing the results
mixed_layer_depth = ecl.open_tutorial_dataset("sample_mixed_layer_depth")["mixed_layer_depth"]
mixed_layer_depth

# %%
# Draw the figure of mix-layer depth
#
fig, ax = ecl.plot.quick_draw_spatial_basemap()

mixed_layer_depth.isel(time = 0).plot(
    vmax = 300,
    cmap = "viridis_r",
    cbar_kwargs = {'location': 'bottom'},
)
ax.set_title("Mixed Layer Depth (Jan., 2021)")

# %%
# Compare with the depth data of the mixed layer output by the oceanic models
#
fig, ax = ecl.plot.quick_draw_spatial_basemap()
diff = mld_data.isel(time = 0) - mixed_layer_depth.isel(time = 0)
diff.plot(
    vmax = 200,
    cbar_kwargs = {'location': 'bottom', 'label': 'units: m'},
)
ax.set_title("SODA minus Easyclimate (pyoceans)")

# %%
#  MLD Internal Temperature
# -------------------------------------
#
# We use :py:func:`easyclimate.field.ocean.get_temper_within_MLD <easyclimate.field.ocean.get_temper_within_MLD>` to receive MLD internal temperature.
#
# .. code-block:: python
#
#     mld_t = ecl.field.ocean.get_temper_within_MLD(
#         seawater_temperature_data = temper_data,
#         mixed_layer_depth = mld_data,
#     ).to_netcdf("sample_mld_t.nc")
#
# Next, we open the dataset containing the results
mld_t = ecl.open_tutorial_dataset("sample_mld_t").temp
mld_t

# %%
# Now, Plotting the temperature of a marine model layer within a mixed layer
ax = plt.figure().add_subplot(projection='3d')

for depth_value in np.arange(50):
    mld_t.isel(time = 5).sel(lon = slice(160, 180), lat = slice(-10, 10)).isel(depth = depth_value).plot.surface(
        alpha=0.3,
    )

# %%
#  MLD Internal Average Temperature
# -------------------------------------
#
# We use :py:func:`easyclimate.field.ocean.calc_MLD_depth_weighted <easyclimate.field.ocean.calc_MLD_depth_weighted>` to calculate MLD internal average temperature.
#
# .. code-block:: python
#
#     weight = ecl.field.ocean.calc_MLD_depth_weighted(
#         seawater_temperature_data = temper_data,
#         mixed_layer_depth = mld_data
#     )
#
#     mld_t_ave = ecl.field.ocean.get_data_average_within_MLD(
#         data_input = temper_data,
#         mixed_layer_depth = mld_data,
#         depth_weight = weight
#     ).to_netcdf("sample_mld_t_ave.nc")
#
#
# .. warning::
#
#     You can **NOT** use above result (i.e., ``mld_t``) to directly calculate MLD internal average temperature by ``mld_t.mean(dim = "depth")``, because the vertical layers in ocean models are usually not uniformly distributed.
#
# Next, we open the dataset containing the results
mld_t_ave = ecl.open_tutorial_dataset("sample_mld_t_ave")["__xarray_dataarray_variable__"]
mld_t_ave

# %%
# Now, Plotting the  MLD internal average temperature
fig, ax = ecl.plot.quick_draw_spatial_basemap()

mld_t_ave.isel(time = 0).plot.contourf(
    levels = 21,
    cbar_kwargs = {'location': 'bottom', 'label': 'degC'},
)
ax.set_title("Mixed Layer Temperature (Jan., 2021)")

# %%
# MLD Average Temperature Tendency
# -------------------------------------
#
# We use :py:func:`easyclimate.field.ocean.calc_MLD_temper_tendency <easyclimate.field.ocean.calc_MLD_temper_tendency>` to calculate MLD internal average temperature tendency.
#
# .. code-block:: python
#
#     weight = ecl.field.ocean.calc_MLD_depth_weighted(
#         seawater_temperature_data = temper_data,
#         mixed_layer_depth = mld_data
#     )
#
#     mld_t_tendency = ecl.field.ocean.calc_MLD_temper_tendency(
#         seawater_temperature_anomaly_data = temper_data,
#         mixed_layer_depth = mld_data,
#         depth_weight = weight
#     ).to_netcdf("sample_mld_t_tendency.nc")
#
mld_t_tendency = ecl.open_tutorial_dataset("sample_mld_t_tendency")["__xarray_dataarray_variable__"]
mld_t_tendency

# %%
# Now, Plotting the result.
fig, ax = ecl.plot.quick_draw_spatial_basemap()

mld_t_tendency.isel(time = 5).plot.contourf(
    vmax = 8,
    levels = 21,
    cbar_kwargs = {'location': 'bottom'},
)
ax.set_title("Mixed Layer Temperature Tendency (Jun., 2021)")

# %%
# MLD Average Horizontal Advection
# -------------------------------------
#
# We use :py:func:`easyclimate.field.ocean.calc_MLD_average_horizontal_advection <easyclimate.field.ocean.calc_MLD_average_horizontal_advection>` to calculate MLD average horizontal advection.
#
# .. code-block:: python
#
#     weight = ecl.field.ocean.calc_MLD_depth_weighted(
#         seawater_temperature_data = temper_data,
#         mixed_layer_depth = mld_data
#     )
#
#     mld_horizontal_advection = ecl.field.ocean.calc_MLD_average_horizontal_advection(
#         u_monthly_data = u_data,
#         v_monthly_data = v_data,
#         seawater_temperature_data = temper_data,
#         mixed_layer_depth = mld_data,
#         depth_weight = weight
#     ).to_netcdf("sample_mld_horizontal_advection.nc")
#
mld_horizontal_advection = ecl.open_tutorial_dataset("sample_mld_horizontal_advection")

# %%
# Now, Plotting the result for zonal advection.
fig, ax = ecl.plot.quick_draw_spatial_basemap()

mld_horizontal_advection["u_advection"].isel(time = 5).plot(
    vmax = 2,
    levels = 21,
    cbar_kwargs = {'location': 'bottom'},
)
ax.set_title("Mixed Layer U advection (Jun., 2021)")

# %%
# And meridional advection.
fig, ax = ecl.plot.quick_draw_spatial_basemap()

mld_horizontal_advection["v_advection"].isel(time = 5).plot(
    vmax = 2,
    levels = 21,
    cbar_kwargs = {'location': 'bottom'},
)
ax.set_title("Mixed Layer V advection (Jun., 2021)")

# %%
# MLD Average Vertical Advection
# -------------------------------------
#
# We use :py:func:`easyclimate.field.ocean.calc_MLD_average_vertical_advection <easyclimate.field.ocean.calc_MLD_average_vertical_advection>` to calculate MLD average vertical advection.
#
# .. code-block:: python
#
#     weight = ecl.field.ocean.calc_MLD_depth_weighted(
#         seawater_temperature_data = temper_data,
#         mixed_layer_depth = mld_data
#     )
#
#     mld_vertical_advection = ecl.field.ocean.calc_MLD_average_vertical_advection(
#         w_monthly_data = wt_data,
#         seawater_temperature_data = temper_data,
#         mixed_layer_depth = mld_data,
#         depth_weight = weight
#     ).to_netcdf("sample_mld_vertical_advection.nc")
#
mld_vertical_advection = ecl.open_tutorial_dataset("sample_mld_vertical_advection")["__xarray_dataarray_variable__"]
mld_vertical_advection

# %%
# Now, Plotting the result.
fig, ax = ecl.plot.quick_draw_spatial_basemap()

mld_vertical_advection.isel(time = 5).plot(
    vmax = 2,
    levels = 21,
    cbar_kwargs = {'location': 'bottom'},
)
ax.set_title("Mixed Layer vertical advection (Jun., 2021)")

# %%
# Surface Heat Flux
# -------------------------------------
#
# We use :py:func:`easyclimate.field.ocean.calc_ocean_surface_heat_flux <easyclimate.field.ocean.calc_ocean_surface_heat_flux>` to calculate ocean surface heat flux.
#
# .. code-block:: python
#
#     surface_heat_flux = ecl.field.ocean.calc_ocean_surface_heat_flux(
#         qnet_monthly_anomaly_data = net_heating_data,
#         mixed_layer_depth = mld_data,
#     ).to_netcdf("sample_surface_heat_flux.nc")
#
surface_heat_flux = ecl.open_tutorial_dataset("sample_surface_heat_flux")["__xarray_dataarray_variable__"]
surface_heat_flux

# %%
# Now, Plotting the result.
fig, ax = ecl.plot.quick_draw_spatial_basemap()

surface_heat_flux.isel(time = 5).plot(
    vmax = 10,
    levels = 21,
    cbar_kwargs = {'location': 'bottom'},
)
ax.set_title("Surface Heat Flux (Jun., 2021)")
