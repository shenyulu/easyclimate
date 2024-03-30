import xarray as xr
import numpy as np
import dask
import pandas as pd
import os
import sys
import numpy as np
import sys
import importlib

sys.path.append("../src")
import easyclimate as ecl

# u_data_sample = ecl.open_tutorial_dataset('uwnd_202201_mon_mean').uwnd.isel(time = 0).sel(level = 500)
# v_data_sample = ecl.open_tutorial_dataset('vwnd_202201_mon_mean').vwnd.isel(time = 0).sel(level = 500)
# u_data_sample.to_netcdf('test_input_uwnd_202201_mon_mean_500hPa_sampledata.nc', format = "NETCDF3_64BIT")
# v_data_sample.to_netcdf('test_input_vwnd_202201_mon_mean_500hPa_sampledata.nc', format = "NETCDF3_64BIT")
u_data_sample = xr.open_dataset("test_input_uwnd_202201_mon_mean_500hPa_sampledata.nc")[
    "uwnd"
]
v_data_sample = xr.open_dataset("test_input_vwnd_202201_mon_mean_500hPa_sampledata.nc")[
    "vwnd"
]

lon_start, lon_end, lat_start, lat_end = 20, 30, -10, 10

# calc_wind_speed
result_data = ecl.windspharm.calc_wind_speed(
    u_data=u_data_sample,
    v_data=v_data_sample,
)
result_data.sel(lon=slice(20, 30), lat=slice(10, -10)).to_dataset(
    name="result"
).to_netcdf("../data/test_output_calc_wind_speed.nc", format="NETCDF3_64BIT")

result_data.sel(
    lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start)
).data.flatten()
result_data.sel(lon=slice(20, 30), lat=slice(10, -10)).to_dataset(
    name="result"
).result.data

# calc_relative_vorticity_and_horizontal_divergence
result_data = ecl.windspharm.calc_relative_vorticity_and_horizontal_divergence(
    u_data=u_data_sample,
    v_data=v_data_sample,
)
result_data.sel(lon=slice(20, 30), lat=slice(10, -10)).to_netcdf(
    "../data/test_output_calc_relative_vorticity_and_horizontal_divergence.nc",
    format="NETCDF3_64BIT",
)

# calc_relative_vorticity
result_data = ecl.windspharm.calc_relative_vorticity(
    u_data=u_data_sample,
    v_data=v_data_sample,
)

# calc_divergence
result_data = ecl.windspharm.calc_divergence(
    u_data=u_data_sample,
    v_data=v_data_sample,
)

# calc_planetary_vorticity
result_data = ecl.windspharm.calc_planetary_vorticity(
    u_data=u_data_sample,
    v_data=v_data_sample,
)
result_data.sel(lon=slice(20, 30), lat=slice(10, -10)).to_dataset(
    name="result"
).to_netcdf("../data/test_output_calc_planetary_vorticity.nc", format="NETCDF3_64BIT")

# calc_absolute_vorticity
result_data = ecl.windspharm.calc_absolute_vorticity(
    u_data=u_data_sample,
    v_data=v_data_sample,
)
result_data.sel(lon=slice(20, 30), lat=slice(10, -10)).to_dataset(
    name="result"
).to_netcdf("../data/test_output_calc_absolute_vorticity.nc", format="NETCDF3_64BIT")

# calc_streamfunction_and_velocity_potential
result_data = ecl.windspharm.calc_streamfunction_and_velocity_potential(
    u_data=u_data_sample,
    v_data=v_data_sample,
)
result_data.sel(lon=slice(20, 30), lat=slice(10, -10)).to_netcdf(
    "../data/test_output_calc_streamfunction_and_velocity_potential.nc",
    format="NETCDF3_64BIT",
)

# calc_helmholtz
result_data = ecl.windspharm.calc_helmholtz(
    u_data=u_data_sample,
    v_data=v_data_sample,
)
result_data.sel(lon=slice(20, 30), lat=slice(10, -10)).to_netcdf(
    "../data/test_output_calc_helmholtz.nc", format="NETCDF3_64BIT"
)

# calc_irrotational_component
result_data = ecl.windspharm.calc_irrotational_component(
    u_data=u_data_sample,
    v_data=v_data_sample,
)
result_data.sel(lon=slice(20, 30), lat=slice(10, -10)).to_netcdf(
    "../data/test_output_calc_irrotational_component.nc", format="NETCDF3_64BIT"
)

# calc_nondivergent_component
result_data = ecl.windspharm.calc_nondivergent_component(
    u_data=u_data_sample,
    v_data=v_data_sample,
)
result_data.sel(lon=slice(20, 30), lat=slice(10, -10)).to_netcdf(
    "../data/test_output_calc_nondivergent_component.nc", format="NETCDF3_64BIT"
)

# calc_rossby_wave_source
result_data = ecl.windspharm.calc_rossby_wave_source(
    u_data=u_data_sample,
    v_data=v_data_sample,
)
result_data.sel(lon=slice(20, 30), lat=slice(10, -10)).to_dataset(
    name="result"
).to_netcdf("../data/test_output_calc_rossby_wave_source.nc", format="NETCDF3_64BIT")

# calc_gradient
result_data = ecl.windspharm.calc_gradient(
    data_input=u_data_sample,
)
result_data.sel(lon=slice(20, 30), lat=slice(10, -10)).to_netcdf(
    "../data/test_output_calc_gradient.nc", format="NETCDF3_64BIT"
)
