"""
pytest for field.ocean.ocean_oceanic_front.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import TEST_DATA_PATH

# sst_data = ecl.open_tutorial_dataset("sst_mnmean_oisst").sst.isel(time=slice(100, 150))
# sst_data.sortby('lat').sel(lon = slice(130, 230), lat = slice(20, 50)).to_netcdf('data/test_input_sst_mnmean_oisst_oceanic_front.nc', format = "NETCDF3_64BIT")
sst_data = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_sst_mnmean_oisst_oceanic_front.nc"))
)["sst"]


def test_calc_intensity_STFZ():
    result_data = ecl.field.ocean.calc_intensity_STFZ(sst_data).data
    refer_data = np.array(
        [
            20.795668,
            22.233015,
            24.089832,
            25.465496,
            26.642035,
            26.999508,
            26.132774,
            24.369974,
            22.543621,
            21.15074,
            19.907673,
            20.054033,
            21.654413,
            22.981853,
            24.9258,
            26.581474,
            27.076017,
            27.268692,
            26.124674,
            24.691538,
            22.67498,
            21.185406,
            20.125143,
            20.019165,
            20.654062,
            21.38998,
            23.841387,
            25.69947,
            26.388285,
            26.385769,
            25.196552,
            23.659357,
            22.106482,
            20.848034,
            19.924639,
            19.958597,
            20.829454,
            22.322617,
            24.170927,
            26.082241,
            26.432878,
            26.576866,
            25.362524,
            23.732878,
            22.373325,
            20.877903,
            20.258863,
            20.166018,
            20.73346,
            22.282312,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)


def test_calc_intensity_SAFZ():
    result_data = ecl.field.ocean.calc_intensity_SAFZ(sst_data).data
    refer_data = np.array(
        [
            10.878135,
            12.396331,
            14.345946,
            18.063534,
            20.892447,
            20.355253,
            18.23466,
            16.277767,
            14.401862,
            12.68262,
            11.44366,
            11.110788,
            11.548044,
            12.521644,
            15.298601,
            18.099924,
            20.345972,
            19.448929,
            17.38222,
            14.900941,
            12.7527075,
            11.037891,
            10.375448,
            10.252701,
            10.606403,
            11.825271,
            13.806818,
            16.959625,
            19.893276,
            19.823277,
            17.295025,
            14.521939,
            12.598631,
            11.3634405,
            10.340692,
            10.063781,
            10.134882,
            10.855593,
            12.909844,
            16.241829,
            18.804844,
            18.888727,
            16.922417,
            14.731883,
            12.527191,
            10.871724,
            10.281294,
            10.199626,
            10.376377,
            11.3605585,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)


def test_calc_location_STFZ():
    result_data = ecl.field.ocean.calc_location_STFZ(sst_data).data
    refer_data = np.array(
        [
            27.747995,
            27.800121,
            27.847696,
            27.892632,
            27.925243,
            27.917028,
            27.889063,
            27.852787,
            27.81767,
            27.786947,
            27.757584,
            27.746325,
            27.771917,
            27.795351,
            27.844662,
            27.897873,
            27.932999,
            27.92573,
            27.89346,
            27.8442,
            27.804647,
            27.7701,
            27.741747,
            27.7287,
            27.726788,
            27.747356,
            27.795578,
            27.856598,
            27.913599,
            27.912598,
            27.86658,
            27.834763,
            27.804028,
            27.785181,
            27.75268,
            27.728058,
            27.722708,
            27.743116,
            27.780666,
            27.869278,
            27.897781,
            27.891571,
            27.856947,
            27.83039,
            27.802177,
            27.764042,
            27.742432,
            27.737612,
            27.742489,
            27.766163,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)


def test_calc_location_SAFZ():
    result_data = ecl.field.ocean.calc_location_SAFZ(sst_data).data
    refer_data = np.array(
        [
            39.36473,
            39.415245,
            39.444813,
            39.5287,
            39.624695,
            39.62453,
            39.560814,
            39.539936,
            39.527435,
            39.486744,
            39.46102,
            39.422844,
            39.39016,
            39.389503,
            39.463783,
            39.532993,
            39.58482,
            39.575512,
            39.54517,
            39.4996,
            39.426487,
            39.368237,
            39.36987,
            39.392414,
            39.397408,
            39.46782,
            39.532,
            39.58699,
            39.604168,
            39.63651,
            39.588104,
            39.494125,
            39.430283,
            39.390804,
            39.372566,
            39.37742,
            39.36397,
            39.370895,
            39.472893,
            39.560665,
            39.606697,
            39.61709,
            39.58999,
            39.508316,
            39.427097,
            39.392982,
            39.373184,
            39.372894,
            39.38131,
            39.38961,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)


def test_calc_location_line_STFZ():
    result_data = (
        ecl.field.ocean.calc_location_line_STFZ(sst_data)
        .isel(lon=slice(20, 25), time=slice(0, 5))
        .data.flatten()
    )
    refer_data = np.array(
        [
            27.755348,
            27.752598,
            27.746893,
            27.74205,
            27.736536,
            27.795942,
            27.796278,
            27.794508,
            27.787716,
            27.785145,
            27.828197,
            27.829203,
            27.82746,
            27.820818,
            27.82114,
            27.881891,
            27.885736,
            27.888783,
            27.897547,
            27.89636,
            27.918003,
            27.922112,
            27.92792,
            27.938538,
            27.94055,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)


def test_calc_location_line_SAFZ():
    result_data = (
        ecl.field.ocean.calc_location_line_SAFZ(sst_data)
        .isel(lon=slice(20, 25), time=slice(0, 5))
        .data.flatten()
    )
    refer_data = np.array(
        [
            39.12821,
            39.16662,
            39.221855,
            39.2772,
            39.31831,
            39.264515,
            39.302525,
            39.3418,
            39.347855,
            39.370026,
            39.36368,
            39.38455,
            39.400692,
            39.398354,
            39.41293,
            39.517666,
            39.52378,
            39.52336,
            39.491528,
            39.49229,
            39.612602,
            39.620857,
            39.624573,
            39.628662,
            39.62745,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)
