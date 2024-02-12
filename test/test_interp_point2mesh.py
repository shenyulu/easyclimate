"""
pytest for interp.interp_point2mesh.py
"""

import pytest

import easyclimate as ecl
import numpy as np

data = ecl.open_tutorial_dataset("PressQFF_202007271200_872.csv")


def test_interp_point2mesh():
    result1 = ecl.interp.interp_point2mesh(
        data,
        var_name="qff",
        grid_x=37.5,
        grid_y=75.0,
        point=[-26.0, 34.5],
        resolution=32,
        sigma=1,
    )
    result_data = result1.sel(lon=slice(20, 20.1), lat=slice(50, 50.1)).data.flatten()
    refer_data = np.array(
        [
            1017.0498,
            1017.0547,
            1017.0596,
            1017.06445,
            1017.05566,
            1017.0605,
            1017.0653,
            1017.0701,
            1017.0612,
            1017.066,
            1017.07074,
            1017.07544,
            1017.06647,
            1017.07117,
            1017.07587,
            1017.0805,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_interp_point2mesh_S2():
    result2 = ecl.interp.interp_point2mesh_S2(
        data,
        var_name="qff",
        grid_x=37.5,
        grid_y=75.0,
        point=[-26.0, 34.5],
        resolution=32,
        sigma=1,
    )
    result_data = result2.sel(lon=slice(20, 20.1), lat=slice(50, 50.1)).data.flatten()
    refer_data = np.array(
        [
            1017.0279,
            1017.0305,
            1017.03296,
            1017.0354,
            1017.0329,
            1017.0354,
            1017.03784,
            1017.04016,
            1017.03766,
            1017.0401,
            1017.0424,
            1017.0447,
            1017.0421,
            1017.04443,
            1017.04675,
            1017.04895,
        ]
    )
    assert np.isclose(result_data, refer_data).all()
