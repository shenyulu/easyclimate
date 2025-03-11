"""
pytest field.heat_stress.humanindexmod_2020.py
"""

import pytest

import easyclimate as ecl
import xarray as xr
import numpy as np


def test_calc_apparent_temperature():
    temperature_data = xr.DataArray([[30.5, 31], [28, 29]])
    vapor_pressure_data = xr.DataArray([[20.0, 21.0], [20.5, 19.8]])
    wind_10m_data = xr.DataArray([[10.0, 11], [9, 12]])

    result_data = ecl.field.heat_stress.calc_apparent_temperature(
        temperature_data, vapor_pressure_data, wind_10m_data, "degC", "Pa", "m/s"
    ).data.flatten()
    refer_data = np.array([19.56599998, 19.36930084, 17.7676506, 16.66534042])

    assert np.isclose(result_data, refer_data, atol=0.01).all()


def test_calc_simplified_human_discomfort_index():
    temperature_data = xr.DataArray([[30.5, 31], [28, 29]])
    vapor_pressure_data = xr.DataArray([[24.25, 21.0], [20.5, 19.8]])
    result_data = ecl.field.heat_stress.calc_simplified_human_discomfort_index(
        temperature_data, vapor_pressure_data, "degC", "Pa"
    ).data.flatten()
    refer_data = np.array([27.375, 26.0, 24.25, 24.39999962])

    assert np.isclose(result_data, refer_data, atol=0.01).all()


def test_calc_simplified_human_discomfort_index_stull():
    temperature_data = xr.DataArray([[30.5, 31], [28, 29]])
    stull_wet_bulb_temperature_data = xr.DataArray([[24.25, 21.0], [20.5, 19.8]])
    rh = xr.DataArray([[75.3, 86.0], [79, 90.8]])

    result_data = ecl.field.heat_stress.calc_simplified_human_discomfort_index_stull(
        temperature_data, stull_wet_bulb_temperature_data, rh, "degC", "degC", "%"
    ).data.flatten()
    refer_data = np.array([52.90000153, 58.5, 53.5, 59.90000153])

    assert np.isclose(result_data, refer_data, atol=0.01).all()


def test_calc_swamp_cooler_temperatures():
    result_data = ecl.field.heat_stress.calc_swamp_cooler_temperatures(
        xr.DataArray([30.25]), xr.DataArray([25.5]), "degC", "degC"
    )
    result_data1 = result_data["swp80"].data
    result_data2 = result_data["swp65"].data
    refer_data1 = np.array([26.45000076])
    refer_data2 = np.array([27.16250038])

    assert np.isclose(result_data1, refer_data1, atol=0.01).all()
    assert np.isclose(result_data2, refer_data2, atol=0.01).all()


def test_calc_heat_thic_thip():
    result_data = ecl.field.heat_stress.calc_heat_thic_thip(
        xr.DataArray([30.25]), xr.DataArray([25.5]), "degC", "degC"
    )
    result_data1 = result_data["thic"].data
    result_data2 = result_data["thip"].data
    refer_data1 = np.array([80.74000549316406])
    refer_data2 = np.array([83.45750427246094])

    assert np.isclose(result_data1, refer_data1, atol=0.01).all()
    assert np.isclose(result_data2, refer_data2, atol=0.01).all()


def test_calc_simplified_wbgt_index():
    result_data = ecl.field.heat_stress.calc_simplified_wbgt_index(
        xr.DataArray([30.25]), xr.DataArray([25.5]), "degC", "Pa"
    ).data.flatten()
    refer_data = np.array([21.1919651])

    assert np.isclose(result_data, refer_data, atol=0.01).all()


def test_calc_human_feels_temperature():
    result_data = ecl.field.heat_stress.calc_human_feels_temperature(
        xr.DataArray([30.0]), xr.DataArray([2968.321]), "degC", "Pa"
    ).data.flatten()
    refer_data = np.array([40.93511962890625])

    assert np.isclose(result_data, refer_data, atol=0.01).all()
