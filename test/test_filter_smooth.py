"""
pytest for interp.interp_mesh2mesh.py
"""

import pytest
import xarray as xr
import numpy as np
import easyclimate as ecl
from pathlib import Path
from .const_define import TEST_DATA_PATH

data_u = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_interp_mesh2mesh.nc"))
).uwnd.sortby("lat")


def test_calc_spatial_smooth_gaussian():
    result_data = (
        ecl.filter.calc_spatial_smooth_gaussian(data_u)
        .sel(lon=slice(110, 120), lat=slice(30, 40))
        .data.flatten()
    )
    refer_data = np.array(
        [
            6.2560167,
            6.397078,
            6.2896986,
            5.9525213,
            5.596702,
            8.067476,
            8.712493,
            8.981655,
            8.932116,
            8.73082,
            9.382056,
            10.385473,
            11.075889,
            11.325369,
            11.306392,
            10.271242,
            11.173315,
            11.870452,
            12.238711,
            12.421031,
            11.301823,
            11.445985,
            11.715482,
            11.9432,
            12.154338,
        ]
    )

    assert np.isclose(result_data, refer_data).all()


def test_calc_spatial_smooth_rectangular():
    result_data = (
        ecl.filter.calc_spatial_smooth_rectangular(data_u)
        .sel(lon=slice(110, 120), lat=slice(30, 40))
        .data.flatten()
    )
    refer_data = np.array(
        [
            6.2500005,
            6.354167,
            6.229167,
            5.9375,
            5.6458335,
            8.083334,
            8.666667,
            8.958334,
            8.9375,
            8.854167,
            9.395834,
            10.395834,
            11.083334,
            11.395834,
            11.479167,
            10.270834,
            11.166667,
            11.875,
            12.3125,
            12.541667,
            11.270834,
            11.4375,
            11.666667,
            11.9375,
            12.145834,
        ]
    )

    assert np.isclose(result_data, refer_data).all()


def test_calc_spatial_smooth_5or9_point1():
    result_data = (
        ecl.filter.calc_spatial_smooth_5or9_point(data_u, n=5)
        .sel(lon=slice(110, 120), lat=slice(30, 40))
        .data.flatten()
    )
    refer_data = np.array(
        [
            6.234375,
            6.3828125,
            6.28125,
            5.9609375,
            5.609375,
            8.046875,
            8.671875,
            8.953125,
            8.90625,
            8.7265625,
            9.359375,
            10.34375,
            11.0078125,
            11.2578125,
            11.2578125,
            10.28125,
            11.140625,
            11.8203125,
            12.1875,
            12.3671875,
            11.3203125,
            11.4765625,
            11.7265625,
            11.9453125,
            12.140625,
        ]
    )

    assert np.isclose(result_data, refer_data).all()


def test_calc_spatial_smooth_5or9_point2():
    result_data = (
        ecl.filter.calc_spatial_smooth_5or9_point(data_u, n=9)
        .sel(lon=slice(110, 120), lat=slice(30, 40))
        .data.flatten()
    )
    refer_data = np.array(
        [
            6.1601562,
            6.3242188,
            6.25,
            5.96875,
            5.65625,
            7.9648438,
            8.550781,
            8.8359375,
            8.816406,
            8.6875,
            9.292969,
            10.183594,
            10.785156,
            11.03125,
            11.078125,
            10.300781,
            11.0546875,
            11.652344,
            12.003906,
            12.1796875,
            11.394531,
            11.5546875,
            11.761719,
            11.953125,
            12.1015625,
        ]
    )

    assert np.isclose(result_data, refer_data).all()


def test_calc_forward_smooth1():
    result_data = (
        ecl.filter.calc_forward_smooth(data_u, n=5)
        .sel(lon=slice(110, 120), lat=slice(30, 40))
        .data.flatten()
    )
    refer_data = np.array(
        [
            6.234375,
            6.3828125,
            6.28125,
            5.9609375,
            5.609375,
            8.046875,
            8.671875,
            8.953125,
            8.90625,
            8.7265625,
            9.359375,
            10.34375,
            11.0078125,
            11.2578125,
            11.2578125,
            10.28125,
            11.140625,
            11.8203125,
            12.1875,
            12.3671875,
            11.3203125,
            11.4765625,
            11.7265625,
            11.9453125,
            12.140625,
        ]
    )

    assert np.isclose(result_data, refer_data).all()


def test_calc_forward_smooth2():
    with pytest.raises(ValueError):
        ecl.filter.calc_forward_smooth(data_u, n=11111).sel(
            lon=slice(110, 120), lat=slice(30, 40)
        ).data.flatten()
        assert 1 == 1


def test_calc_reverse_smooth():
    result_data = (
        ecl.filter.calc_reverse_smooth(data_u, n=9)
        .sel(lon=slice(110, 120), lat=slice(30, 40))
        .data.flatten()
    )
    refer_data = np.array(
        [
            8.050781,
            8.152344,
            7.953125,
            7.359375,
            6.859375,
            10.308594,
            11.316406,
            11.5390625,
            11.441406,
            10.96875,
            11.964844,
            13.433594,
            14.566406,
            14.875,
            14.65625,
            12.738281,
            14.3046875,
            15.371094,
            15.847656,
            16.101562,
            13.925781,
            13.9921875,
            14.527344,
            14.90625,
            15.3359375,
        ]
    )

    assert np.isclose(result_data, refer_data).all()


def test_calc_spatial_smooth_circular():
    result_data = (
        ecl.filter.calc_spatial_smooth_circular(data_u, radius=6)
        .sel(lon=slice(110, 120), lat=slice(30, 40))
        .data.flatten()
    )
    refer_data = np.array(
        [
            4.2582965,
            4.3849564,
            4.507744,
            4.602876,
            4.655974,
            6.1161513,
            6.2096243,
            6.286505,
            6.3296466,
            6.3263283,
            7.896571,
            7.956858,
            7.982301,
            7.96405,
            7.8993373,
            9.286504,
            9.310839,
            9.292035,
            9.232301,
            9.134955,
            10.105089,
            10.079093,
            10.019911,
            9.936946,
            9.841814,
        ]
    )

    assert np.isclose(result_data, refer_data).all()
