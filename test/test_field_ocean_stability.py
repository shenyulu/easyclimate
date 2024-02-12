"""
pytest for field.ocean.stability.py
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


def test_calc_N2_from_temp_salt_1():
    temper_data_need = temper_data.isel(time=0)
    slt_data_need = slt_data.isel(time=0)

    result_data = ecl.field.ocean.calc_N2_from_temp_salt(
        seawater_temperature_data=temper_data_need,
        seawater_practical_salinity_data=slt_data_need,
        time_dim=None,
    )
    result_data1 = result_data.N2.isel(depth=2).data.flatten()
    result_data2 = result_data.p_mid.isel(depth=2).data.flatten()

    refer_data1 = np.array(
        [
            4.53344119e-04,
            1.33144867e-05,
            5.16107056e-04,
            6.45770019e-04,
            2.23375396e-04,
            4.52393316e-04,
            5.02642128e-04,
            1.62480486e-04,
            1.28272873e-04,
            np.nan,
            7.11007684e-04,
            1.59181926e-03,
            np.nan,
            np.nan,
            np.nan,
            3.96496549e-04,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            3.14042028e-07,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    refer_data2 = np.array(
        [
            30.28890133,
            30.28890133,
            30.28890133,
            30.28890133,
            30.28890133,
            30.28890133,
            30.28890133,
            30.28890133,
            30.28890133,
            np.nan,
            30.28890133,
            30.28890133,
            np.nan,
            np.nan,
            np.nan,
            30.28890133,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            30.28890133,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )

    assert np.allclose(result_data1, refer_data1, equal_nan=True)
    assert np.allclose(result_data2, refer_data2, equal_nan=True)


def test_calc_N2_from_temp_salt_2():
    temper_data_need = temper_data.isel(time=slice(0, 2))
    slt_data_need = slt_data.isel(time=slice(0, 2))

    result_data = ecl.field.ocean.calc_N2_from_temp_salt(
        seawater_temperature_data=temper_data_need,
        seawater_practical_salinity_data=slt_data_need,
        time_dim="time",
    ).isel(time=0)
    result_data1 = result_data.N2.isel(depth=2).data.flatten()
    result_data2 = result_data.p_mid.isel(depth=2).data.flatten()

    refer_data1 = np.array(
        [
            4.53344119e-04,
            1.33144867e-05,
            5.16107056e-04,
            6.45770019e-04,
            2.23375396e-04,
            4.52393316e-04,
            5.02642128e-04,
            1.62480486e-04,
            1.28272873e-04,
            np.nan,
            7.11007684e-04,
            1.59181926e-03,
            np.nan,
            np.nan,
            np.nan,
            3.96496549e-04,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            3.14042028e-07,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    refer_data2 = np.array(
        [
            30.28890133,
            30.28890133,
            30.28890133,
            30.28890133,
            30.28890133,
            30.28890133,
            30.28890133,
            30.28890133,
            30.28890133,
            np.nan,
            30.28890133,
            30.28890133,
            np.nan,
            np.nan,
            np.nan,
            30.28890133,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            30.28890133,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )

    assert np.allclose(result_data1, refer_data1, equal_nan=True)
    assert np.allclose(result_data2, refer_data2, equal_nan=True)


def test_calc_potential_density_from_temp_salt_1():
    temper_data_need = temper_data.isel(time=0)
    slt_data_need = slt_data.isel(time=0)

    result_data = (
        ecl.field.ocean.calc_potential_density_from_temp_salt(
            seawater_temperature_data=temper_data_need,
            seawater_practical_salinity_data=slt_data_need,
            time_dim=None,
        )
        .prho.isel(depth=2)
        .data.flatten()
    )

    refer_data = np.array(
        [
            1020.52295558,
            1019.49839832,
            1020.19481173,
            1020.14431607,
            1020.24774695,
            1020.36358272,
            1020.88569038,
            1020.24960776,
            1020.45356427,
            np.nan,
            1020.30535071,
            1019.63130614,
            np.nan,
            np.nan,
            np.nan,
            1019.95530528,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            1018.22028207,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)


def test_calc_potential_density_from_temp_salt_2():
    temper_data_need = temper_data.isel(time=slice(0, 2))
    slt_data_need = slt_data.isel(time=slice(0, 2))

    result_data = (
        ecl.field.ocean.calc_potential_density_from_temp_salt(
            seawater_temperature_data=temper_data_need,
            seawater_practical_salinity_data=slt_data_need,
            time_dim="time",
        )
        .prho.isel(time=0)
        .isel(depth=2)
        .data.flatten()
    )

    refer_data = np.array(
        [
            1020.52295558,
            1019.49839832,
            1020.19481173,
            1020.14431607,
            1020.24774695,
            1020.36358272,
            1020.88569038,
            1020.24960776,
            1020.45356427,
            np.nan,
            1020.30535071,
            1019.63130614,
            np.nan,
            np.nan,
            np.nan,
            1019.95530528,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            1018.22028207,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)
