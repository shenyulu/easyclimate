# -*- coding: utf-8 -*-
"""
Geographic Finite Difference
===================================

Before proceeding with all the steps, first import some necessary libraries and packages
"""
import easyclimate as ecl
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# %%
# Then consider obtaining meridional and zonal wind variables in tutorial data

u_data = ecl.tutorial.open_tutorial_dataset("uwnd_202201_mon_mean").sortby("lat").uwnd
v_data = ecl.tutorial.open_tutorial_dataset("vwnd_202201_mon_mean").sortby("lat").vwnd
z_data = ecl.tutorial.open_tutorial_dataset("hgt_202201_mon_mean").sortby("lat").hgt
temp_data = ecl.tutorial.open_tutorial_dataset("air_202201_mon_mean").sortby("lat").air
q_data = ecl.tutorial.open_tutorial_dataset("shum_202201_mon_mean").sortby("lat").shum
msl_data = (
    ecl.tutorial.open_tutorial_dataset("pressfc_202201_mon_mean").sortby("lat").pres
)
pr_data = (
    ecl.tutorial.open_tutorial_dataset("precip_202201_mon_mean").sortby("lat").precip
)

uvdata = xr.Dataset()
uvdata["uwnd"] = u_data
uvdata["vwnd"] = v_data

# %%
# Obtain data slices on 500hPa isobars for January 2022

uvdata_500_202201 = uvdata.sel(level=500, time="2022-01-01")
z_data_500_202201 = z_data.sel(level=500, time="2022-01-01")
temp_data_500_202201 = temp_data.sel(level=500, time="2022-01-01")

# %%
# Plotting a sample `quiver` plot of this data slice

fig, ax = plt.subplots(
    subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)}
)

ax.stock_img()
ax.gridlines(draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--")
ax.coastlines(edgecolor="black", linewidths=0.5)

uvdata_500_202201.thin(lon=3, lat=3).plot.quiver(
    ax=ax,
    u="uwnd",
    v="vwnd",
    x="lon",
    y="lat",
    # projection on data
    transform=ccrs.PlateCarree(),
)

# %%
# First-order Partial Derivative
# -------------------------------------
#
# Consider the function :py:func:`easyclimate.calc_gradient <easyclimate.calc_gradient>` to compute the gradient of the zonal wind with respect to longitude.
#
# .. math::
#
#   \frac{\partial u}{\partial \lambda}
#
# The argument `dim` to the function :py:func:`easyclimate.calc_gradient <easyclimate.calc_gradient>` specifies that the direction of the solution is `longitude`.

uwnd_dx = ecl.calc_gradient(uvdata_500_202201.uwnd, dim="lon")

uwnd_dx

# %%

fig, ax = plt.subplots(
    subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)}
)

ax.stock_img()
ax.gridlines(draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--")
ax.coastlines(edgecolor="black", linewidths=0.5)

uwnd_dx.plot.contourf(
    ax=ax,
    # projection on data
    transform=ccrs.PlateCarree(),
    # Colorbar is placed at the bottom
    cbar_kwargs={"location": "bottom"},
    levels=21,
)

# %%
# Of course, it is also possible to pass in :py:class:`xarray.Dataset<xarray.Dataset>` directly into the function :py:func:`easyclimate.calc_gradient <easyclimate.calc_gradient>` to iterate through all the variables, so that you can get the gradient of both the zonal and meridional winds with respect to longitude at the same time.

uvwnd_dx = ecl.calc_gradient(uvdata_500_202201, dim="lon")

uvwnd_dx

# %%
# However, if one is required to solve for the gradient of the zonal wind with respect to the corresponding distance at each longitude, the function `calc_lon_gradient` should be used to calculate.
#
# .. math::
#   \frac{\partial F}{\partial x} = \frac{1}{R \cos\varphi} \cdot \frac{\partial F}{\partial \lambda}

uwnd_dlon = ecl.calc_lon_gradient(uvdata_500_202201.uwnd, lon_dim="lon", lat_dim="lat")

fig, ax = plt.subplots(
    subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)}
)

ax.gridlines(draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--")
ax.coastlines(edgecolor="black", linewidths=0.5)

uwnd_dlon.plot.contourf(
    ax=ax,
    # projection on data
    transform=ccrs.PlateCarree(),
    # Colorbar is placed at the bottom
    cbar_kwargs={"location": "bottom"},
    levels=21,
)

# %%
# Similarly, use :py:func:`easyclimate.calc_lat_gradient <easyclimate.calc_lat_gradient>` to solve for the gradient of the meridional wind with respect to the corresponding distance at each latitude.

# %%
# Second-order Partial Derivative
# ------------------------------------
#
# The solution of the second-order partial derivative relies on three functional calculations
#
# - :py:func:`easyclimate.calc_lon_laplacian <easyclimate.calc_lon_laplacian>`: calculation of the second-order partial derivative term (Laplace term) along longitude.
#
# .. math::
#   \frac{\partial^2 F}{\partial x^2} = \frac{1}{(R \cos\varphi)^2} \cdot \frac{\partial^2 F}{\partial \lambda^2}

uwnd_dlon2 = ecl.calc_lon_laplacian(
    uvdata_500_202201.uwnd, lon_dim="lon", lat_dim="lat"
)

# %%
#
# - :py:func:`easyclimate.calc_lat_laplacian <easyclimate.calc_lat_laplacian>`: calculation of the second-order partial derivative term (Laplace term) along latitude.
#
# .. math::
#   \frac{\partial^2 F}{\partial y^2} = \frac{1}{R^2} \cdot \frac{\partial^2 F}{\partial \varphi^2}

uwnd_dlat2 = ecl.calc_lat_laplacian(uvdata_500_202201.uwnd, lat_dim="lat")

# %%
#
# - :py:func:`easyclimate.calc_lon_lat_mixed_derivatives <easyclimate.calc_lon_lat_mixed_derivatives>`: second-order mixed partial derivative terms along longitude and latitude.
#
# .. math::
#   \frac{\partial^2 F}{\partial x \partial y} = \frac{1}{R^2 \cos\varphi} \cdot \frac{\partial^2 F}{\partial \lambda \partial \varphi}
#

uwnd_dlonlat = ecl.calc_lon_lat_mixed_derivatives(
    uvdata_500_202201.uwnd, lon_dim="lon", lat_dim="lat"
)

# %%
# Second-order partial derivative term along longitude.

fig, ax = plt.subplots(
    subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)}
)

ax.gridlines(draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--")
ax.coastlines(edgecolor="black", linewidths=0.5)

uwnd_dlon2.plot.contourf(
    ax=ax, transform=ccrs.PlateCarree(), cbar_kwargs={"location": "bottom"}, levels=21
)
ax.set_title("$\\frac{\\partial^2 F}{\\partial x^2}$", fontsize=20)

# %%
# Second-order partial derivative term along latitude.

fig, ax = plt.subplots(
    subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)}
)

ax.gridlines(draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--")
ax.coastlines(edgecolor="black", linewidths=0.5)

uwnd_dlat2.plot.contourf(
    ax=ax, transform=ccrs.PlateCarree(), cbar_kwargs={"location": "bottom"}, levels=21
)
ax.set_title("$\\frac{\\partial^2 F}{\\partial y^2}$", fontsize=20)

# %%
# Second-order mixed partial derivative terms along longitude and latitude.

fig, ax = plt.subplots(
    subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)}
)

ax.gridlines(draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--")
ax.coastlines(edgecolor="black", linewidths=0.5)

uwnd_dlonlat.plot.contourf(
    ax=ax, transform=ccrs.PlateCarree(), cbar_kwargs={"location": "bottom"}, levels=21
)
ax.set_title("$\\frac{\\partial^2 F}{\\partial x \\partial y}$", fontsize=20)

# %%
# Vorticity and Divergence
# ------------------------------------
#
# Vorticity and divergence are measures of the degree of atmospheric rotation and volumetric flux per unit volume respectively. For vorticity and divergence in the quasi-geostrophic case, the potential height is used as input data for the calculations. In general, we first calculate the quasi-geostrophic wind.
#
# - :py:func:`easyclimate.calc_geostrophic_wind <easyclimate.calc_geostrophic_wind>`: calculate the geostrophic wind.
#
# .. math::
#   u_g = - \frac{g}{f} \frac{\partial H}{\partial y}, \ v_g = \frac{g}{f} \frac{\partial H}{\partial x}
#

geostrophic_wind_data_500_202201 = ecl.calc_geostrophic_wind(
    z_data_500_202201, lon_dim="lon", lat_dim="lat"
)

# %%
# The function :py:func:`easyclimate.calc_vorticity <easyclimate.calc_vorticity>` is then used to compute the quasi-geostrophic vorticity.
#
# - :py:func:`easyclimate.calc_vorticity <easyclimate.calc_vorticity>`: calculate the horizontal relative vorticity term.
#
# .. math::
#   \zeta = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y} + \frac{u}{R} \tan \varphi
#

qg_vor_data_500_202201 = ecl.calc_vorticity(
    u_data=geostrophic_wind_data_500_202201.ug,
    v_data=geostrophic_wind_data_500_202201.vg,
    lon_dim="lon",
    lat_dim="lat",
)

qg_vor_data_500_202201.sel(lat=slice(20, 80)).plot.contourf(levels=21)

# %%
# Similar vorticity for actual winds, but for actual winds rather than quasi-geostrophic winds.

vor_data_500_202201 = ecl.calc_vorticity(
    u_data=uvdata_500_202201["uwnd"],
    v_data=uvdata_500_202201["vwnd"],
    lon_dim="lon",
    lat_dim="lat",
)

vor_data_500_202201.sel(lat=slice(20, 80)).plot.contourf(levels=21)

# %%
# In addition, the function :py:func:`easyclimate.calc_divergence <easyclimate.calc_divergence>` calculate the quasi-geostrophic divergence.
#
# .. math::
#   \mathrm{D} = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} - \frac{v}{R} \tan \varphi
#
# - :py:func:`easyclimate.calc_divergence <easyclimate.calc_divergence>`: calculate the horizontal divergence term.
#
# Quasi-geostrophic divergence

qg_div_data_500_202201 = ecl.calc_divergence(
    u_data=geostrophic_wind_data_500_202201.ug,
    v_data=geostrophic_wind_data_500_202201.vg,
    lon_dim="lon",
    lat_dim="lat",
)

qg_div_data_500_202201.sel(lat=slice(20, 80)).plot.contourf(levels=21)

# %%
# Actual divergence

div_data_500_202201 = ecl.calc_divergence(
    u_data=uvdata_500_202201["uwnd"],
    v_data=uvdata_500_202201["vwnd"],
    lon_dim="lon",
    lat_dim="lat",
)

div_data_500_202201.sel(lat=slice(20, 80)).plot.contourf(levels=21)

# %%
# Of course, in addition to the built-in finite difference method, the spherical harmonic function mothod can be solved, but you must ensure that it is **Global** and **Regular or Gaussian grid** type data.
#
# - :py:func:`easyclimate.windspharm.calc_relative_vorticity <easyclimate.windspharm.top.calc_relative_vorticity>`: calculate the relative vorticity term with the spherical harmonic function mothod.
# - :py:func:`easyclimate.windspharm.calc_divergence <easyclimate.windspharm.top.calc_divergence>`: calculate the horizontal divergence term with the spherical harmonic function mothod.

vor_data_500_202201_windspharm = ecl.windspharm.calc_relative_vorticity(
    u_data=uvdata_500_202201["uwnd"],
    v_data=uvdata_500_202201["vwnd"],
)

vor_data_500_202201_windspharm.sortby("lat").sel(lat=slice(20, 80)).plot.contourf(
    levels=21
)

# %%

div_data_500_202201_windspharm = ecl.windspharm.calc_divergence(
    u_data=uvdata_500_202201["uwnd"],
    v_data=uvdata_500_202201["vwnd"],
)

div_data_500_202201_windspharm.sortby("lat").sel(lat=slice(20, 80)).plot.contourf(
    levels=21
)

# %%
# Generally speaking, the calculation results of the finite difference method and the spherical harmonic function method are similar. The former does not require global regional data, but the calculation results of the latter are more accurate for high latitude regions.

# %%
# Advection
# -------------------
# `Advection <https://glossary.ametsoc.org/wiki/Advection>`__ is the process of transport of an atmospheric property solely by the mass motion (velocity field) of the atmosphere; also, the rate of change of the value of the advected property at a given point.
#
# For zonal advection, we can calculate as follows.
#
# .. math::
#   -u \frac{\partial T}{\partial x}
#

u_advection_500_202201 = ecl.calc_u_advection(
    u_data=uvdata_500_202201["uwnd"], temper_data=temp_data_500_202201
)

u_advection_500_202201.sortby("lat").sel(lat=slice(20, 80)).plot.contourf(levels=21)

# %%
# Similarly, the meridional advection can acquire as follows.
#
# .. math::
#   -v \frac{\partial T}{\partial y}
#

v_advection_500_202201 = ecl.calc_v_advection(
    v_data=uvdata_500_202201["vwnd"], temper_data=temp_data_500_202201
)

v_advection_500_202201.sortby("lat").sel(lat=slice(20, 80)).plot.contourf(levels=21)

# %%
# Water Flux
# -------------------------
#
# - :py:func:`easyclimate.calc_horizontal_water_flux <easyclimate.calc_horizontal_water_flux>`: calculate horizontal water vapor flux at each vertical level.
#
# .. math::
#   \frac{1}{g} q \mathbf{V} = \frac{1}{g} (u q\ \mathbf{i} + vq\ \mathbf{j})
#
# - :py:func:`easyclimate.calc_vertical_water_flux <easyclimate.calc_vertical_water_flux>`: calculate vertical water vapor flux.
#
# .. math::
#   -\omega \frac{q}{g}
#
# - :py:func:`easyclimate.calc_water_flux_top2surface_integral <easyclimate.calc_water_flux_top2surface_integral>`: calculate the water vapor flux across the vertical level.
#
# :py:func:`easyclimate.calc_horizontal_water_flux <easyclimate.calc_horizontal_water_flux>` can calculate the horizontal water flux of single layers.

ecl.calc_horizontal_water_flux(
    specific_humidity_data=q_data,
    u_data=uvdata.uwnd,
    v_data=uvdata.vwnd,
)

# %%
# The whole layer integral needs to consider the function :py:func:`easyclimate.calc_water_flux_top2surface_integral <easyclimate.calc_water_flux_top2surface_integral>` to calculate.

water_flux_top2surface_integral = ecl.calc_water_flux_top2surface_integral(
    specific_humidity_data=q_data,
    u_data=u_data,
    v_data=v_data,
    surface_pressure_data=msl_data,
    surface_pressure_data_units="millibars",
    vertical_dim="level",
    vertical_dim_units="hPa",
)

water_flux_top2surface_integral

# %%
# Extracting the entire layer water vapor flux at mid and low latitudes at the 0th time level.

draw_water_flux = (
    water_flux_top2surface_integral.isel(time=0)
    .thin(lon=3, lat=3)
    .sel(lat=slice(-60, 60))
)
draw_pr = pr_data.isel(time=0).sel(lat=slice(-60, 60))

# %%
fig, ax = plt.subplots(
    subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)}
)

ax.gridlines(draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--")
ax.coastlines(edgecolor="black", linewidths=0.5)

draw_water_flux.plot.quiver(
    ax=ax,
    u="qu",
    v="qv",
    x="lon",
    y="lat",
    transform=ccrs.PlateCarree(),
    zorder=2,
)

draw_pr.plot.contourf(
    ax=ax,
    transform=ccrs.PlateCarree(),
    levels=21,
    cmap="Greens",
    zorder=1,
    cbar_kwargs={"location": "bottom"},
    vmax=20,
)

# %%
# Water Vapor Flux Divergence
# --------------------------------------
#
# Water vapor flux divergence represents the convergence and divergence of water vapor. There are also two built-in functions to calculate the results of single-layers and whole-layer integration respectively.
#
# - :py:func:`easyclimate.calc_divergence_watervaporflux <easyclimate.calc_divergence_watervaporflux>`: calculate water vapor flux divergence at each vertical level.
#
# .. math::
#   \nabla \left( \frac{1}{g} q \mathbf{V} \right) = \frac{1}{g} \nabla \cdot \left( q \mathbf{V} \right)
#
# - :py:func:`easyclimate.calc_divergence_watervaporflux_top2surface_integral <easyclimate.calc_divergence_watervaporflux_top2surface_integral>`: calculate water vapor flux divergence across the vertical level.

divergence_watervaporflux_top2surface_integral = (
    ecl.calc_divergence_watervaporflux_top2surface_integral(
        specific_humidity_data=q_data,
        u_data=u_data,
        v_data=v_data,
        surface_pressure_data=msl_data,
        surface_pressure_data_units="millibars",
        specific_humidity_units="grams/kg",
        vertical_dim="level",
        vertical_dim_units="hPa",
    )
)

divergence_watervaporflux_top2surface_integral

# %%
# Extracting the entire layer water vapor flux at mid and low latitudes at the 0th time level.

draw_data = divergence_watervaporflux_top2surface_integral.isel(time=0).sel(
    lat=slice(-60, 60)
)

# %%
#

fig, ax = plt.subplots(
    subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)}
)

ax.gridlines(draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--")
ax.coastlines(edgecolor="black", linewidths=0.5)

draw_data.plot.contourf(
    ax=ax, transform=ccrs.PlateCarree(), cbar_kwargs={"location": "bottom"}, levels=21
)
