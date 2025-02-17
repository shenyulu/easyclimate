"""
pytest for filter/lanczos_filter.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import TEST_DATA_PATH

data_lanczos = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_field_lanczos_bandpass.nc"))
).vo


def test_calc_lanczos_lowpass():
    result_data1 = (
        ecl.filter.calc_lanczos_lowpass(data_lanczos, 20, 3, method="rolling")
        .isel(time=slice(50, 55))
        .data
    )
    result_data2 = (
        ecl.filter.calc_lanczos_lowpass(data_lanczos, 20, 3, method="convolve")
        .isel(time=slice(50, 55))
        .data
    )
    refer_data = np.array(
        [
            -1.43267264e-05,
            -2.05010938e-05,
            -7.22385184e-06,
            -5.29967042e-06,
            -3.43212725e-06,
        ]
    )

    assert np.isclose(result_data1, refer_data, atol=0.01).all()
    assert np.isclose(result_data2, refer_data, atol=0.01).all()


def test_calc_lanczos_bandpass():
    result_data1 = (
        ecl.filter.calc_lanczos_bandpass(data_lanczos, 20, [3, 10], method="rolling")
        .isel(time=slice(50, 55))
        .data
    )
    result_data2 = (
        ecl.filter.calc_lanczos_bandpass(data_lanczos, 20, [3, 10], method="convolve")
        .isel(time=slice(50, 55))
        .data
    )
    refer_data = np.array(
        [
            8.21783158e-06,
            1.33842534e-05,
            3.07898814e-07,
            -5.61406041e-07,
            -1.11504102e-06,
        ]
    )

    assert np.isclose(result_data1, refer_data, atol=0.01).all()
    assert np.isclose(result_data2, refer_data, atol=0.01).all()


def test_calc_lanczos_highpass():
    result_data1 = (
        ecl.filter.calc_lanczos_highpass(data_lanczos, 20, 10, method="rolling")
        .isel(time=slice(50, 55))
        .data
    )
    result_data2 = (
        ecl.filter.calc_lanczos_highpass(data_lanczos, 20, 10, method="convolve")
        .isel(time=slice(50, 55))
        .data
    )
    refer_data = np.array(
        [
            -7.70544478e-06,
            -1.61722139e-05,
            7.45221104e-07,
            4.93167542e-06,
            -8.88079159e-06,
        ]
    )

    assert np.isclose(result_data1, refer_data, atol=0.01).all()
    assert np.isclose(result_data2, refer_data, atol=0.01).all()
