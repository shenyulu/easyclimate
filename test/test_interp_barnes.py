"""
pytest for interp.barnes.py
"""

import pytest

import easyclimate as ecl
import numpy as np

data = ecl.open_tutorial_dataset("PressQFF_202007271200_872.csv")


def test_interp_spatial_barnes():
    result1 = ecl.interp.interp_spatial_barnes(
        data,
        var_name="qff",
        grid_x=12,
        grid_y=12,
        point=[-9, 47],
        resolution=32,
        sigma=1.0,
    )
    result_data = result1.sel(lon=slice(-4, -3.9), lat=slice(50, 50.1)).data.flatten()
    refer_data = np.array(
        [
            1005.87036,
            1005.9002,
            1005.93,
            1005.9598,
            1005.7945,
            1005.82465,
            1005.8548,
            1005.88495,
            1005.7168,
            1005.7473,
            1005.7778,
            1005.8082,
            1005.6374,
            1005.6681,
            1005.69885,
            1005.7296,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_interp_spatial_barnesS2():
    result2 = ecl.interp.interp_spatial_barnesS2(
        data,
        var_name="qff",
        grid_x=12,
        grid_y=12,
        point=[-9, 47],
        resolution=32,
        sigma=1.0,
    )
    result_data = result2.sel(lon=slice(-4, -3.9), lat=slice(50, 50.1)).data.flatten()
    refer_data = np.array(
        [
            1005.8769,
            1005.8949,
            1005.9129,
            1005.9309,
            1005.8007,
            1005.8189,
            1005.837,
            1005.8551,
            1005.72345,
            1005.7417,
            1005.75995,
            1005.7782,
            1005.645,
            1005.6633,
            1005.6817,
            1005.7001,
        ]
    )
    assert np.isclose(result_data, refer_data).all()
