# -*- coding: utf-8 -*-
"""
.. _curved_quiver_plot:

Curved Quiver Plot
============================================
Here is an example of how to draw a curved quiver
"""
import easyclimate as ecl
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

# %%
# Import the required sample wind farm data
udata = ecl.open_tutorial_dataset("uwnd_2022_day5")["uwnd"].sortby("lat")
vdata = ecl.open_tutorial_dataset("vwnd_2022_day5")["vwnd"].sortby("lat")
uvdata = xr.Dataset(data_vars = {"u": udata,"v": vdata})
uvdata



# %%
# East Asia Sample
# -----------------------------------------------
# Consideration of surface wind field data for the East Asian region
lonlat_range = [90, 180, 0, 70]
gap_value = 5

drawdata_gap1 = uvdata.isel(time = 2).sel(level = 850).sel(
    lon = slice(lonlat_range[0] - gap_value, lonlat_range[1] + gap_value),
    lat = slice(lonlat_range[2] - gap_value, lonlat_range[3] + gap_value)
)

# %%
# Plotting the results of a traditional quiver
fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree(central_longitude=120)})

ax.coastlines()
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=["bottom", "left"], alpha = 0.4)
ax.set_extent(lonlat_range, crs = ccrs.PlateCarree())

q = drawdata_gap1.plot.quiver(
    x = 'lon', y = 'lat', u = 'u', v = 'v',
    ax = ax,
    color = '#f6631c',  zorder = 0,
    width = 0.007,
    regrid_shape = 15,
    scale = 200,
    add_guide = False,
    transform = ccrs.PlateCarree(),
)

qk = ax.quiverkey(
    q, 0.9, 1.02,
    10, "10", labelpos = "N",
    color = "k",
    coordinates = "axes",
    fontproperties = {"size": 12},
    labelsep = 0.05,
    transform = ccrs.PlateCarree(),
)

ax.set_title("UV850 standard quiver (2022-01-03)")

# %%
# Now use a curved quiver to plot the wind field
fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree(central_longitude=120)})

ax.coastlines()
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=["bottom", "left"], alpha = 0.4)
ax.set_extent(lonlat_range, crs = ccrs.PlateCarree())

cq = ecl.plot.curved_quiver(
    drawdata_gap1,
    x = 'lon', y = 'lat', u = 'u', v = 'v',
    ax = ax,
    density = 5,
    color = '#f6631c',
    transform = ccrs.PlateCarree(),
)

ecl.plot.add_curved_quiverkey(
    cq, ax = ax,
    pos = (0.9, 1.05),
    U = 10,
    label = "10",
    color = "black",
    labelpos = "N",
    fontproperties = {"size": 12},
    ref_point = (120, 30)
)

ax.set_title("UV850 curved quiver (2022-01-03)")

# %%
# Maritime Continent Sample
# -----------------------------------------------
# Here, we will use an oceanic-continental region as an example of a wind farm in a tropical area.
gap_value = 5
lonlat_range = [80, 160, -30, 30]
drawdata_gap6 = uvdata.isel(time = 2).sel(level = 850).sel(
    lon = slice(lonlat_range[0] - gap_value, lonlat_range[1] + gap_value),
    lat = slice(lonlat_range[2] - gap_value, lonlat_range[3] + gap_value)
)

# %%
# Similarly, first plot the results of the traditional Quiver wind field visualization.
fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree(central_longitude=120)})

ax.coastlines()
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=["bottom", "left"], alpha = 0.4)
ax.set_extent(lonlat_range, crs = ccrs.PlateCarree())

q = drawdata_gap6.plot.quiver(
    x = 'lon', y = 'lat', u = 'u', v = 'v',
    ax = ax,
    color = '#f6631c',  zorder = 0,
    width = 0.007,
    regrid_shape = 13,
    scale = 150,
    add_guide = False,
    transform = ccrs.PlateCarree(),
)

qk = ax.quiverkey(
    q, 0.9, 1.02,
    10, "10", labelpos = "N",
    color = "k",
    coordinates = "axes",
    fontproperties = {"size": 12},
    labelsep = 0.05,
    transform = ccrs.PlateCarree(),
)
qk.set_zorder(7)

ax.set_title("UV850 quiver (2022-01-03)")

# %%
# However, this is the result of drawing a curved quiver
fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree(central_longitude=120)})

ax.coastlines()
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=["bottom", "left"], alpha = 0.4)
ax.set_extent(lonlat_range, crs = ccrs.PlateCarree())

cq = ecl.plot.curved_quiver(
    drawdata_gap6,
    x = 'lon', y = 'lat', u = 'u', v = 'v',
    ax = ax,
    density = 5,
    color = '#f6631c',
    transform = ccrs.PlateCarree(),
)

ecl.plot.add_curved_quiverkey(
    cq, ax = ax,
    pos = (0.9, 1.05),
    U = 10,
    label = "10",
    color = "black",
    labelpos = "N",
    fontproperties = {"size": 12},
    ref_point = (120, 30),
)

ax.set_title("UV850 curved quiver (2022-01-03)")

# %%
# Pole Sample
# -----------------------------------------------
# Mapping the Arctic region requires additional parameter settings; here, we will first extract the data for the Arctic region.
gap_value = 5
drawdata_gap3 = uvdata.isel(time = 2).sel(level = 500).sel(
    lat = slice(0 - gap_value, None)
)
drawdata_gap3 = ecl.plot.add_lon_cyclic(drawdata_gap3, inter = 2.5)

# %%
# Traditional quiver data requires an additional setting for the ``regrid_shape`` value.
fig, ax = plt.subplots(
    figsize = (6, 6),
    subplot_kw={"projection": ccrs.NorthPolarStereo(central_longitude=150)}
)
ax.coastlines()

gl, meta = ecl.plot.draw_polar_basemap(
    ax=ax,
    lon_step=30,
    lat_step=30,
    lat_range=[0, 90],
    draw_labels=True,
    lat_label_lon=-30,
    set_map_boundary_kwargs = {"south_pad": 0.9}
)

q = drawdata_gap3.plot.quiver(
    x = 'lon', y = 'lat', u = 'u', v = 'v',
    ax = ax,
    color = '#f6631c',  zorder = 0,
    width = 0.007,
    regrid_shape = 18,
    scale = 300,
    add_guide = False,
    transform = ccrs.PlateCarree(),
)

qk = ax.quiverkey(
    q, 0.05, 0.05,
    20, "20", labelpos = "N",
    color = "k",
    coordinates = "axes",
    fontproperties = {"size": 12},
    labelsep = 0.05,
    transform = ccrs.PlateCarree(),
)

ax.set_title("")
ecl.plot.set_polar_title(ax, "UV850 quiver(2022-01-03)", meta = meta, size = 15)

# %%
# However, in a curved quiver, you need to set not only the `regrid_shape`
# but also the magnitude (``ref_magnitude``) and length (``ref_length``) of the reference
# vector to better plot the wind field.

# sphinx_gallery_thumbnail_number = 2
fig, ax = plt.subplots(
    figsize = (6, 6),
    subplot_kw={"projection": ccrs.NorthPolarStereo(central_longitude=150)}
)
ax.coastlines()

gl, meta = ecl.plot.draw_polar_basemap(
    ax=ax,
    lon_step=30,
    lat_step=30,
    lat_range=[0, 90],
    draw_labels=True,
    lat_label_lon=-30,
    set_map_boundary_kwargs = {"south_pad": 0.9}
)

cq = ecl.plot.curved_quiver(
    drawdata_gap3,
    x = 'lon', y = 'lat', u = 'u', v = 'v',
    ax = ax,
    density = 20,
    color = '#f6631c',
    regrid_shape = 11,
    ref_magnitude = 20,
    ref_length = 0.3,
    transform = ccrs.PlateCarree(),
)

ecl.plot.add_curved_quiverkey(
    cq, ax = ax,
    pos = (0.05, 0.05),
    U = 10,
    label = "10",
    color = "black",
    labelpos = "N",
    fontproperties = {"size": 12},
    ref_point = (120, 30),
)

ax.set_title("")
ecl.plot.set_polar_title(ax, "UV850 curved quiver (2022-01-03)", meta = meta, size = 15)

# %%
# Tropics Sample
# -----------------------------------------------
# The following example covers the entire tropical region.
# For cases where the area is too large,
# additional parameters must be configured;
# one cannot simply rely on the default results.
gap_value = 5

lonlat_range = [0, 360, -40, 40]
drawdata_gap4 = uvdata.isel(time = 2).sel(level = 850).sel(
    lon = slice(lonlat_range[0] - gap_value, lonlat_range[1] + gap_value),
    lat = slice(lonlat_range[2] - gap_value, lonlat_range[3] + gap_value)
)
drawdata_gap4 = ecl.plot.add_lon_cyclic(drawdata_gap4, inter = 2.5)

# %%
# Let’s start by looking at the rendering results for quiver.
fig, ax = plt.subplots(figsize = (15, 3) ,subplot_kw={"projection": ccrs.PlateCarree(central_longitude=200)})

ax.coastlines()
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=["bottom", "left"], alpha = 0.4)
ax.set_extent(lonlat_range, crs = ccrs.PlateCarree())

q = drawdata_gap4.plot.quiver(
    x = 'lon', y = 'lat', u = 'u', v = 'v',
    ax = ax,
    color = '#f6631c', zorder = 0,
    regrid_shape = 10,
    scale = 400,
    add_guide = False,
    transform = ccrs.PlateCarree(),
)

qk = ax.quiverkey(
    q, 0.9, 1.02,
    10, "10", labelpos = "N",
    color = "k",
    coordinates = "axes",
    fontproperties = {"size": 12},
    labelsep = 0.05,
    transform = ccrs.PlateCarree(),
)

ax.set_title("UV850 quiver(2022-01-03)")

# %%
# But here is the rendering result for the curved quiver.
fig, ax = plt.subplots(figsize = (15, 3),subplot_kw={"projection": ccrs.PlateCarree(central_longitude=200)})

ax.coastlines()
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=["bottom", "left"], alpha = 0.4)
ax.set_extent(lonlat_range, crs = ccrs.PlateCarree())

cq = ecl.plot.curved_quiver(
    drawdata_gap4,
    x = 'lon', y = 'lat', u = 'u', v = 'v',
    ax = ax,
    ref_magnitude = 15,
    density = 10,
    color = '#f6631c',
    regrid_shape = 15,
    transform = ccrs.PlateCarree(),
)

ecl.plot.add_curved_quiverkey(
    cq, ax = ax,
    pos = (0.9, 1.05),
    U = 2,
    label = "2",
    color = "black",
    labelpos = "N",
    fontproperties = {"size": 12},
    ref_point = (120, 30),
)

ax.set_title("UV850 curved quiver (2022-01-03)")
