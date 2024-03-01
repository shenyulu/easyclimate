"""
pytest for field.ocean.thermal.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import TEST_DATA_PATH

temper_data = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_field_ocean_mixlayer_temper.nc"))
).temp.isel(time=1)


def test_calc_seawater_thermocline_depth():
    result_data = ecl.field.ocean.calc_seawater_thermocline_depth(
        seawater_temperature_data=temper_data
    ).data.flatten()
    refer_data = np.array(
        [
            76.88051605,
            np.nan,
            66.30594889,
            55.89711634,
            45.5960172,
            66.30594889,
            25.22615083,
            np.nan,
            66.30594889,
            np.nan,
            66.30594889,
            25.22615083,
            np.nan,
            np.nan,
            np.nan,
            66.30594889,
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


def test_calc_D14_depth():
    result_data = ecl.field.ocean.calc_D14_depth(
        seawater_temperature_data=temper_data
    ).data.flatten()
    refer_data = np.array(
        [
            201.39993804,
            26.07418473,
            101.47570178,
            114.16131459,
            36.55500891,
            201.31653075,
            23.6723689,
            35.00541862,
            90.16086258,
            np.nan,
            201.3215514,
            23.9220486,
            np.nan,
            np.nan,
            np.nan,
            148.10938009,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            43.42732817,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)


def test_calc_D17_depth():
    result_data = ecl.field.ocean.calc_D17_depth(
        seawater_temperature_data=temper_data
    ).data.flatten()
    refer_data = np.array(
        [
            135.04039348,
            29.48559621,
            105.03672927,
            117.77883839,
            39.81546714,
            135.22242671,
            20.56572622,
            31.50449012,
            93.62098688,
            np.nan,
            135.15840681,
            20.71109501,
            np.nan,
            np.nan,
            np.nan,
            135.14180481,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            40.21432621,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)


def test_calc_D20_depth():
    result_data = ecl.field.ocean.calc_D20_depth(
        seawater_temperature_data=temper_data
    ).data.flatten()
    refer_data = np.array(
        [
            98.91157353,
            31.88121449,
            107.54939955,
            120.32214872,
            42.13344297,
            109.73471159,
            18.34673806,
            29.05087076,
            96.05762344,
            np.nan,
            98.73247183,
            18.4332164,
            np.nan,
            np.nan,
            np.nan,
            98.58349809,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            37.96518136,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)


def test_calc_D26_depth():
    result_data = ecl.field.ocean.calc_D26_depth(
        seawater_temperature_data=temper_data
    ).data.flatten()
    refer_data = np.array(
        [
            55.95681065,
            35.02347415,
            87.52184342,
            76.59972076,
            45.21002112,
            15.11376721,
            15.38814205,
            25.83859492,
            87.2669918,
            np.nan,
            55.90108701,
            15.41612968,
            np.nan,
            np.nan,
            np.nan,
            55.70906107,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            14.78040122,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)


def test_calc_D28_depth():
    result_data = ecl.field.ocean.calc_D28_depth(
        seawater_temperature_data=temper_data
    ).data.flatten()
    refer_data = np.array(
        [
            45.11465309,
            14.77262182,
            45.49125558,
            25.14717116,
            55.74399892,
            25.7507407,
            5.11929397,
            25.32994262,
            66.17781028,
            np.nan,
            25.90189441,
            5.22396715,
            np.nan,
            np.nan,
            np.nan,
            45.22050821,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            14.08347404,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)
