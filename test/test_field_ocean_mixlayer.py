"""
pytest for field.ocean.mixlayer.py
"""

import pytest

import easyclimate as ecl
import xarray as xr
import numpy as np
from pathlib import Path
from .const_define import TEST_DATA_PATH

temper_data = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_field_ocean_mixlayer_temper.nc"))
).temp
slt_data = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_field_ocean_mixlayer_slt.nc"))
).salt
mld_data = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_field_ocean_mixlayer_mld.nc"))
).mlp
u_data = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_field_ocean_mixlayer_u.nc"))
).u
v_data = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_field_ocean_mixlayer_v.nc"))
).v
wt_data = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_field_ocean_mixlayer_wt.nc"))
).wt
net_heating_data = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_field_ocean_mixlayer_net_heating.nc"))
).net_heating

weight = ecl.field.ocean.calc_MLD_depth_weighted(
    seawater_temperature_data=temper_data, mixed_layer_depth=mld_data
)


def test_calc_mixed_layer_depth():
    tmp = ecl.field.ocean.calc_mixed_layer_depth(
        seawater_temperature_data=temper_data, seawater_practical_salinity_data=slt_data
    )
    result_data = tmp.isel(time=0).data.flatten()
    refer_data = np.array(
        [
            15.19069251,
            35.57104821,
            15.19069251,
            15.19069251,
            15.19069251,
            5.06369366,
            5.06369366,
            15.19145504,
            25.37163017,
            np.nan,
            5.06396135,
            5.06396135,
            np.nan,
            np.nan,
            np.nan,
            15.19310089,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            35.57875188,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)


def test_calc_MLD_depth_weighted():
    tmp = ecl.field.ocean.calc_MLD_depth_weighted(
        seawater_temperature_data=temper_data,
        mixed_layer_depth=mld_data,
    )
    result_data = tmp.isel(time=0).isel(lat=0, lon=0).data.flatten()
    refer_data = np.array(
        [
            1.0,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)


def test_calc_MLD_temper_tendency():
    # weight = ecl.field.ocean.calc_MLD_depth_weighted(seawater_temperature_data = temper_data, mixed_layer_depth = mld_data)
    tmp = ecl.field.ocean.calc_MLD_temper_tendency(
        seawater_temperature_anomaly_data=temper_data,
        mixed_layer_depth=mld_data,
        depth_weight=weight,
    )
    result_data = tmp.isel(time=1).data.flatten()
    refer_data = np.array(
        [
            -13.70223808,
            -13.78194427,
            -13.6286974,
            -13.63245869,
            -14.06196022,
            -13.44122601,
            -13.69342899,
            -13.59566402,
            -13.67652416,
            0.0,
            -13.2894516,
            -13.75876141,
            0.0,
            0.0,
            0.0,
            -13.48992443,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)


def test_get_temper_within_MLD_AND_get_data_within_MLD():
    tmp = ecl.field.ocean.get_temper_within_MLD(
        seawater_temperature_data=temper_data,
        mixed_layer_depth=mld_data,
    )
    result_data = tmp.isel(time=1, lon=0, lat=0).data
    refer_data = np.array(
        [
            26.409521,
            26.49613,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)


def test_get_data_average_within_MLD_AND_get_temper_average_within_MLD():
    # weight = ecl.field.ocean.calc_MLD_depth_weighted(seawater_temperature_data = temper_data, mixed_layer_depth = mld_data)
    tmp = ecl.field.ocean.get_data_average_within_MLD(
        data_input=temper_data, mixed_layer_depth=mld_data, depth_weight=weight
    )
    result_data = tmp.isel(time=1).data.flatten()
    refer_data = np.array(
        [
            26.47447777,
            27.09043503,
            27.32452202,
            27.59849548,
            28.62206459,
            26.08473492,
            26.77102661,
            27.6846981,
            27.36911392,
            np.nan,
            25.60866547,
            26.83507347,
            np.nan,
            np.nan,
            np.nan,
            26.65126991,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)


def test_calc_MLD_average_horizontal_advection():
    # weight = ecl.field.ocean.calc_MLD_depth_weighted(seawater_temperature_data = temper_data, mixed_layer_depth = mld_data)
    tmp = ecl.field.ocean.calc_MLD_average_horizontal_advection(
        u_monthly_data=u_data,
        v_monthly_data=v_data,
        seawater_temperature_data=temper_data,
        mixed_layer_depth=mld_data,
        depth_weight=weight,
    ).isel(time=1)
    result_data1 = tmp["u_advection"].data.flatten()
    result_data2 = tmp["v_advection"].data.flatten()

    refer_data1 = np.array(
        [
            1.30081884,
            0.4999351,
            -0.19766832,
            -0.40578124,
            -0.3136038,
            0.93783085,
            2.74617153,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    refer_data2 = np.array(
        [
            -1.02387063,
            0.06824966,
            np.nan,
            np.nan,
            np.nan,
            -1.54469189,
            -0.16930563,
            np.nan,
            np.nan,
            np.nan,
            0.56896169,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            -0.15559094,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )

    assert np.allclose(result_data1, refer_data1, equal_nan=True)
    assert np.allclose(result_data2, refer_data2, equal_nan=True)


def test_calc_MLD_average_vertical_advection():
    # weight = ecl.field.ocean.calc_MLD_depth_weighted(seawater_temperature_data = temper_data, mixed_layer_depth = mld_data)
    tmp = ecl.field.ocean.calc_MLD_average_vertical_advection(
        w_monthly_data=wt_data,
        seawater_temperature_data=temper_data,
        mixed_layer_depth=mld_data,
        depth_weight=weight,
    ).isel(time=1)
    result_data = tmp.data.flatten()
    refer_data = np.array(
        [
            0.1827295,
            -2.73277808,
            -4.22751977,
            -0.71883941,
            -1.19724071,
            1.22931125,
            83.20039965,
            -1.62991016,
            -0.50281687,
            np.nan,
            0.12738364,
            91.5155874,
            np.nan,
            np.nan,
            np.nan,
            -3.6676629,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)


def test_calc_ocean_surface_heat_flux():
    tmp = ecl.field.ocean.calc_ocean_surface_heat_flux(
        qnet_monthly_anomaly_data=net_heating_data, mixed_layer_depth=mld_data
    ).isel(time=1)
    result_data = tmp.data.flatten()
    refer_data = np.array(
        [
            -0.50035418,
            1.41873461,
            1.86946731,
            2.70730791,
            3.17968158,
            -0.8347024,
            3.91468393,
            3.22870301,
            5.18895496,
            np.nan,
            -0.39295859,
            89.63778746,
            np.nan,
            np.nan,
            np.nan,
            -0.69614864,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            6.13829344,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)
