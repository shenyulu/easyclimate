"""
pytest for windspharm.top.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import TEST_DATA_PATH
from .util import (
    round_sf_np_new,
)  # Intel fortran outputs for Windows and linux are quit different

u_data_sample = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_uwnd_202201_mon_mean_500hPa_sampledata.nc"))
)["uwnd"]
v_data_sample = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_vwnd_202201_mon_mean_500hPa_sampledata.nc"))
)["vwnd"]

lon_start, lon_end, lat_start, lat_end = 20, 30, -10, 10


def test_calc_wind_speed():
    result_data = ecl.windspharm.calc_wind_speed(
        u_data=u_data_sample,
        v_data=v_data_sample,
    )
    result_data = result_data.sel(
        lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start)
    ).data
    refer_data = (
        xr.open_dataset(str(Path(TEST_DATA_PATH, "test_output_calc_wind_speed.nc")))[
            "result"
        ]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    assert np.isclose(
        round_sf_np_new(result_data.flatten()), round_sf_np_new(refer_data.flatten())
    ).all()


def test_calc_relative_vorticity_and_horizontal_divergence():
    result_data = ecl.windspharm.calc_relative_vorticity_and_horizontal_divergence(
        u_data=u_data_sample,
        v_data=v_data_sample,
    )
    result_data1 = (
        result_data["vrt"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    result_data2 = (
        result_data["div"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    refer_data = xr.open_dataset(
        str(
            Path(
                TEST_DATA_PATH,
                "test_output_calc_relative_vorticity_and_horizontal_divergence.nc",
            )
        )
    )
    refer_data1 = (
        refer_data["vrt"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    refer_data2 = (
        refer_data["div"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    assert np.isclose(
        round_sf_np_new(result_data1.flatten()), round_sf_np_new(refer_data1.flatten())
    ).all()
    assert np.isclose(
        round_sf_np_new(result_data2.flatten()), round_sf_np_new(refer_data2.flatten())
    ).all()


def test_calc_relative_vorticity():
    result_data = ecl.windspharm.calc_relative_vorticity(
        u_data=u_data_sample,
        v_data=v_data_sample,
    )
    result_data1 = result_data.sel(
        lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start)
    ).data
    refer_data = xr.open_dataset(
        str(
            Path(
                TEST_DATA_PATH,
                "test_output_calc_relative_vorticity_and_horizontal_divergence.nc",
            )
        )
    )
    refer_data1 = (
        refer_data["vrt"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    assert np.isclose(
        round_sf_np_new(result_data1.flatten()), round_sf_np_new(refer_data1.flatten())
    ).all()


def test_calc_divergence():
    result_data = ecl.windspharm.calc_divergence(
        u_data=u_data_sample,
        v_data=v_data_sample,
    )
    result_data1 = result_data.sel(
        lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start)
    ).data
    refer_data = xr.open_dataset(
        str(
            Path(
                TEST_DATA_PATH,
                "test_output_calc_relative_vorticity_and_horizontal_divergence.nc",
            )
        )
    )
    refer_data1 = (
        refer_data["div"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    assert np.isclose(
        round_sf_np_new(result_data1.flatten()), round_sf_np_new(refer_data1.flatten())
    ).all()


def test_calc_planetary_vorticity():
    result_data = ecl.windspharm.calc_planetary_vorticity(
        u_data=u_data_sample,
        v_data=v_data_sample,
    )
    result_data = result_data.sel(
        lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start)
    ).data
    refer_data = (
        xr.open_dataset(
            str(Path(TEST_DATA_PATH, "test_output_calc_planetary_vorticity.nc"))
        )["result"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    assert np.isclose(
        round_sf_np_new(result_data.flatten()), round_sf_np_new(refer_data.flatten())
    ).all()


def test_calc_absolute_vorticity():
    result_data = ecl.windspharm.calc_absolute_vorticity(
        u_data=u_data_sample,
        v_data=v_data_sample,
    )
    result_data = result_data.sel(
        lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start)
    ).data
    refer_data = (
        xr.open_dataset(
            str(Path(TEST_DATA_PATH, "test_output_calc_absolute_vorticity.nc"))
        )["result"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    assert np.isclose(
        round_sf_np_new(result_data.flatten()), round_sf_np_new(refer_data.flatten())
    ).all()


def test_calc_streamfunction_and_velocity_potential():
    result_data = ecl.windspharm.calc_streamfunction_and_velocity_potential(
        u_data=u_data_sample,
        v_data=v_data_sample,
    )
    result_data1 = (
        result_data["stream"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    result_data2 = (
        result_data["pv"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    refer_data = xr.open_dataset(
        str(
            Path(
                TEST_DATA_PATH,
                "test_output_calc_streamfunction_and_velocity_potential.nc",
            )
        )
    )
    refer_data1 = (
        refer_data["stream"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    refer_data2 = (
        refer_data["pv"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    assert np.isclose(
        round_sf_np_new(result_data1.flatten()), round_sf_np_new(refer_data1.flatten())
    ).all()
    assert np.isclose(
        round_sf_np_new(result_data2.flatten()), round_sf_np_new(refer_data2.flatten())
    ).all()


def test_calc_streamfunction():
    result_data = ecl.windspharm.calc_streamfunction(
        u_data=u_data_sample,
        v_data=v_data_sample,
    )
    result_data1 = result_data.sel(
        lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start)
    ).data
    refer_data = xr.open_dataset(
        str(
            Path(
                TEST_DATA_PATH,
                "test_output_calc_streamfunction_and_velocity_potential.nc",
            )
        )
    )
    refer_data1 = (
        refer_data["stream"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    assert np.isclose(
        round_sf_np_new(result_data1.flatten()), round_sf_np_new(refer_data1.flatten())
    ).all()


def test_calc_velocity_potential():
    result_data = ecl.windspharm.calc_velocity_potential(
        u_data=u_data_sample,
        v_data=v_data_sample,
    )
    result_data1 = result_data.sel(
        lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start)
    ).data
    refer_data = xr.open_dataset(
        str(
            Path(
                TEST_DATA_PATH,
                "test_output_calc_streamfunction_and_velocity_potential.nc",
            )
        )
    )
    refer_data1 = (
        refer_data["pv"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    assert np.isclose(
        round_sf_np_new(result_data1.flatten()), round_sf_np_new(refer_data1.flatten())
    ).all()


def test_calc_helmholtz():
    result_data = ecl.windspharm.calc_helmholtz(
        u_data=u_data_sample,
        v_data=v_data_sample,
    )
    result_data1 = (
        result_data["uchi"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data.flatten()
    )
    result_data1 = round_sf_np_new(result_data1)
    result_data2 = (
        result_data["vchi"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data.flatten()
    )
    result_data2 = round_sf_np_new(result_data2)
    result_data3 = (
        result_data["upsi"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data.flatten()
    )
    result_data3 = round_sf_np_new(result_data3)
    result_data4 = (
        result_data["vpsi"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data.flatten()
    )
    result_data4 = round_sf_np_new(result_data4)
    refer_data = xr.open_dataset(
        str(Path(TEST_DATA_PATH, "test_output_calc_helmholtz.nc"))
    )
    refer_data1 = (
        refer_data["uchi"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data.flatten()
    )
    refer_data1 = round_sf_np_new(refer_data1)
    refer_data2 = (
        refer_data["vchi"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data.flatten()
    )
    refer_data2 = round_sf_np_new(refer_data2)
    refer_data3 = (
        refer_data["upsi"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data.flatten()
    )
    refer_data3 = round_sf_np_new(refer_data3)
    refer_data4 = (
        refer_data["vpsi"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data.flatten()
    )
    refer_data4 = round_sf_np_new(refer_data4)
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()


def test_calc_irrotational_component():
    result_data = ecl.windspharm.calc_irrotational_component(
        u_data=u_data_sample,
        v_data=v_data_sample,
    )
    result_data1 = (
        result_data["uchi"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data.flatten()
    )
    result_data2 = (
        result_data["vchi"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data.flatten()
    )
    refer_data = xr.open_dataset(
        str(Path(TEST_DATA_PATH, "test_output_calc_irrotational_component.nc"))
    )
    refer_data1 = (
        refer_data["uchi"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data.flatten()
    )
    refer_data2 = (
        refer_data["vchi"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data.flatten()
    )
    assert np.isclose(
        round_sf_np_new(result_data1.flatten()), round_sf_np_new(refer_data1.flatten())
    ).all()
    assert np.isclose(
        round_sf_np_new(result_data2.flatten()), round_sf_np_new(refer_data2.flatten())
    ).all()


def test_calc_nondivergent_component():
    result_data = ecl.windspharm.calc_nondivergent_component(
        u_data=u_data_sample,
        v_data=v_data_sample,
    )
    result_data1 = (
        result_data["upsi"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data.flatten()
    )
    result_data2 = (
        result_data["vpsi"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data.flatten()
    )
    refer_data = xr.open_dataset(
        str(Path(TEST_DATA_PATH, "test_output_calc_nondivergent_component.nc"))
    )
    refer_data1 = (
        refer_data["upsi"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data.flatten()
    )
    refer_data2 = (
        refer_data["vpsi"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data.flatten()
    )
    assert np.isclose(
        round_sf_np_new(result_data1.flatten()), round_sf_np_new(refer_data1.flatten())
    ).all()
    assert np.isclose(
        round_sf_np_new(result_data2.flatten()), round_sf_np_new(refer_data2.flatten())
    ).all()


def test_calc_rossby_wave_source():
    result_data = ecl.windspharm.calc_rossby_wave_source(
        u_data=u_data_sample,
        v_data=v_data_sample,
    )
    result_data = result_data.sel(
        lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start)
    ).data
    refer_data = (
        xr.open_dataset(
            str(Path(TEST_DATA_PATH, "test_output_calc_rossby_wave_source.nc"))
        )["result"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    assert np.isclose(
        round_sf_np_new(result_data.flatten()), round_sf_np_new(refer_data.flatten())
    ).all()


def test_calc_gradient():
    result_data = ecl.windspharm.calc_gradient(
        data_input=u_data_sample,
    )
    result_data1 = (
        result_data["zonal_gradient"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    result_data2 = (
        result_data["meridional_gradient"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    refer_data = xr.open_dataset(
        str(Path(TEST_DATA_PATH, "test_output_calc_gradient.nc"))
    )
    refer_data1 = (
        refer_data["zonal_gradient"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    refer_data2 = (
        refer_data["meridional_gradient"]
        .sel(lon=slice(lon_start, lon_end), lat=slice(lat_end, lat_start))
        .data
    )
    assert np.isclose(
        round_sf_np_new(result_data1.flatten()), round_sf_np_new(refer_data1.flatten())
    ).all()
    assert np.isclose(
        round_sf_np_new(result_data2.flatten()), round_sf_np_new(refer_data2.flatten())
    ).all()
