"""
pytest for eof.py
"""

import pytest
import warnings

import easyclimate as ecl
import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
from .const_define import TEST_DATA_PATH

data_time_series = xr.DataArray(
    np.array(
        [
            [57413.867, 57411.96],
            [57404.35, 57404.35],
            [57478.6, 57482.406],
            [57531.906, 57531.906],
            [57564.273, 57569.984],
            [57474.79, 57472.887],
            [57309.15, 57307.246],
            [57343.42, 57343.42],
            [57408.152, 57402.44],
            [57442.426, 57442.426],
            [57478.6, 57480.504],
            [57585.215, 57583.312],
            [57568.082, 57569.984],
            [57625.195, 57627.1],
            [57653.758, 57653.758],
            [57697.547, 57697.547],
            [57613.773, 57619.49],
            [57501.445, 57501.445],
            [57446.23, 57446.23],
            [57440.523, 57438.617],
            [57438.617, 57440.523],
            [57387.21, 57381.5],
            [57507.156, 57509.062],
            [57469.08, 57469.08],
            [57349.133, 57349.133],
            [57371.984, 57373.883],
            [57499.54, 57503.35],
            [57512.867, 57516.676],
            [57421.48, 57427.195],
            [57259.65, 57261.555],
        ]
    ),
    dims=("time", "lon"),
    coords={
        "time": pd.date_range("1982-01-01", periods=30, freq="ME"),
        "lon": np.array([100.125, 101.25]),
    },
)

ds1 = xr.DataArray(
    np.array([3e-5, 5.6e-5, 7e-5]), dims="lat", coords={"lat": np.array([10, -20, 30])}
)
ds2 = xr.DataArray(
    np.array([3e-5, 5.4e-5, 7.3e-5]), dims="lon", coords={"lon": np.array([0, 320, 78])}
)
ds3 = xr.DataArray(
    np.array([[3e-5, 5.4e-5], [5.2, -75.5]]),
    dims=("lon", "lat"),
    coords={"lon": np.array([160, 70]), "lat": np.array([87.5, -87.5])},
)
ds4 = xr.DataArray(
    np.array([3e-5, 5.4e-5, 7.3e-5]),
    dims="lon",
    coords={"lon": np.array([-178, -50, 120])},
)
ds5 = xr.DataArray(
    np.array([3e-5, 5.4e-5, 7.3e-5]),
    dims="lon",
    coords={"lon": np.array([120, 180, 310])},
)

comparedata1 = xr.DataArray(
    dims=("time", "level", "lat", "lon"),
    coords={
        "time": np.array(
            [
                "1948-01-01T00:00:00.000000000",
                "1948-02-01T00:00:00.000000000",
                "1948-03-01T00:00:00.000000000",
            ],
            dtype="datetime64[ns]",
        ),
        "level": np.array([1000.0, 925.0, 850.0], dtype=np.float32),
        "lat": np.array([10.0, 7.5, 5.0, 2.5, 0.0], dtype=np.float32),
        "lon": np.array([100.0, 102.5, 105.0, 107.5, 110.0], dtype=np.float32),
    },
)

comparedata2 = xr.DataArray(
    dims=("level", "lat", "lon"),
    coords={
        "level": np.array([1000.0, 925.0, 850.0], dtype=np.float32),
        "lat": np.array([10.0, 7.5, 5.0, 2.5, 0.0], dtype=np.float32),
        "lon": np.array([100.0, 102.5, 105.0, 107.5, 110.0], dtype=np.float32),
    },
)

comparedata3 = xr.DataArray(
    dims=("time", "level", "lat", "lon"),
    coords={
        "time": np.array(
            [
                "1948-01-01T00:00:00.000000000",
                "1948-02-01T00:00:00.000000000",
                "1948-03-01T00:00:00.000000000",
            ],
            dtype="datetime64[ns]",
        ),
        "level": np.array([1000.0, 925.0, 850.0], dtype=np.float32),
        "lat": np.array([10.0, 7.5, 5.0, 2.5, 0.0, -2.5], dtype=np.float32),
        "lon": np.array([100.0, 102.5, 105.0, 107.5, 110.0], dtype=np.float32),
    },
)

comparedata4 = xr.DataArray(
    dims=("time", "level", "lat", "lon"),
    coords={
        "time": np.array(
            [
                "1948-01-01T00:00:00.000000000",
                "1948-02-01T00:00:00.000000000",
                "1948-03-01T00:00:00.000000000",
            ],
            dtype="datetime64[ns]",
        ),
        "level": np.array([1000.0, 925.0, 500.0], dtype=np.float32),
        "lat": np.array([10.0, 7.5, 5.0, 2.5, 0.0], dtype=np.float32),
        "lon": np.array([100.0, 102.5, 105.0, 107.5, 110.0], dtype=np.float32),
    },
)

comparedata5 = xr.DataArray(
    dims=("time", "level", "lat", "lon"),
    coords={
        "time": np.array([1, 2, 3]),
        "level": np.array([1000.0, 925.0, 850.0], dtype=np.float32),
        "lat": np.array([10.0, 7.5, 5.0, 2.5, 0.0], dtype=np.float32),
        "lon": np.array([100.0, 102.5, 105.0, 107.5, 110.0], dtype=np.float32),
    },
)


def test_assert_compared_version():
    result_data1 = ecl.utility.assert_compared_version("0.0.1", "0.0.2")
    result_data2 = ecl.utility.assert_compared_version("0.0.2", "0.0.1")
    result_data3 = ecl.utility.assert_compared_version("0.0.2", "0.0.2")

    refer_data1 = -1
    refer_data2 = 1
    refer_data3 = 0
    assert result_data1 == refer_data1
    assert result_data2 == refer_data2
    assert result_data3 == refer_data3


def test_find_dims_axis():
    result_data = ecl.utility.find_dims_axis(data_time_series, dim="time")
    refer_data = 0
    assert result_data == refer_data


def test_transfer_int2datetime():
    intyear = np.array(
        [2054, 2061, 2062, 2067, 2071, 2075, 2076, 2078, 2085, 2089, 2096]
    )
    result_data = ecl.utility.transfer_int2datetime(intyear)
    result_data = np.datetime_as_string(result_data)

    refer_data = np.array(
        [
            "2054-01-01T00:00:00.000000000",
            "2061-01-01T00:00:00.000000000",
            "2062-01-01T00:00:00.000000000",
            "2067-01-01T00:00:00.000000000",
            "2071-01-01T00:00:00.000000000",
            "2075-01-01T00:00:00.000000000",
            "2076-01-01T00:00:00.000000000",
            "2078-01-01T00:00:00.000000000",
            "2085-01-01T00:00:00.000000000",
            "2089-01-01T00:00:00.000000000",
            "2096-01-01T00:00:00.000000000",
        ]
    )

    assert (result_data == refer_data).all()


def test_split_datetime2yearday():
    result_data = (
        ecl.utility.split_datetime2yearday(data_time_series.isel(lon=0))
        .sel(year=1983)
        .data
    )
    refer_data = np.array(
        [
            57568.082,
            57625.195,
            57653.758,
            57697.547,
            57613.773,
            57501.445,
            57446.23,
            57440.523,
            57438.617,
            57387.21,
            57507.156,
            57469.08,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_transfer_deg2rad():
    ds = xr.DataArray(
        np.array([30, 45, 60]), dims="lon", coords={"lon": np.array([10, 20, 30])}
    )
    result_data = ecl.utility.transfer_deg2rad(ds).data
    refer_data = np.array([0.52359878, 0.78539816, 1.04719755])
    assert np.isclose(result_data, refer_data).all()


def test_transfer_inf2nan():
    ds = xr.DataArray(
        np.array([np.inf, 45, 60]), dims="lon", coords={"lon": np.array([10, 20, 30])}
    )
    result_data = ecl.utility.transfer_inf2nan(ds).data
    refer_data = np.array([np.nan, 45.0, 60.0])
    # Compare arrays include np.nan
    assert np.allclose(result_data, refer_data, equal_nan=True)


def test_transfer_nan2value():
    ds = xr.DataArray(
        np.array([np.nan, 45, 60]), dims="lon", coords={"lon": np.array([10, 20, 30])}
    )
    result_data = ecl.utility.transfer_nan2value(ds, 2).data
    refer_data = np.array([2.0, 45.0, 60.0])
    assert np.isclose(result_data, refer_data).all()


def test_get_weighted_spatial_data():
    inputdata = xr.open_dataset(
        str(Path(TEST_DATA_PATH, "test_input_core_eof.nc"))
    ).z.isel(time=0)
    result_data1 = (
        ecl.utility.get_weighted_spatial_data(inputdata, method="cos_lat").mean().data
    )
    result_data2 = (
        ecl.utility.get_weighted_spatial_data(inputdata, method="area").mean().data
    )
    refer_data1 = 57410.45244069
    refer_data2 = 57410.45261421
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()


def test_get_compress_xarraydata():
    pass


def test_transfer_dFdp2dFdz():
    ds = xr.DataArray(
        np.array([3e-5, 5.6e-5, 7e-5]),
        dims="lon",
        coords={"lon": np.array([10, 20, 30])},
    )
    result_data = ecl.utility.transfer_dFdp2dFdz(ds).data
    refer_data = np.array([-0.3800832, -0.70948864, -0.8868608])
    assert np.isclose(result_data, refer_data).all()


def test_sort_ascending_latlon_coordinates():
    result_data1 = ecl.utility.sort_ascending_latlon_coordinates(
        ds1, lat_dim="lat", lon_dim=None
    ).data
    result_data2 = ecl.utility.sort_ascending_latlon_coordinates(
        ds2, lon_dim="lon", lat_dim=None
    ).data
    result_data3 = ecl.utility.sort_ascending_latlon_coordinates(
        ds3, lon_dim="lon", lat_dim="lat"
    ).data.flatten()
    refer_data1 = np.array([5.6e-05, 3.0e-05, 7.0e-05])
    refer_data2 = np.array([3.0e-05, 7.3e-05, 5.4e-05])
    refer_data3 = np.array([-7.55e01, 5.20e00, 5.40e-05, 3.00e-05])
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()


def test_transfer_units_coeff():
    result_data1 = ecl.utility.transfer_units_coeff("m/s", "km/h")
    result_data2 = ecl.utility.transfer_units_coeff("hPa", "mbar")
    result_data3 = ecl.utility.transfer_units_coeff("mm/day", "m/month")
    refer_data1 = 3.5999999999999996
    refer_data2 = 1.0000000000000002
    refer_data3 = 0.0304375
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()


def test_transfer_data_units():
    result_data = ecl.utility.transfer_data_units(ds3, "mm/day", "m/day").data.flatten()
    refer_data = np.array([3.00e-08, 5.40e-08, 5.20e-03, -7.55e-02])
    assert np.isclose(result_data, refer_data).all()


def test_generate_dataset_dispatcher():
    @ecl.utility.generate_dataset_dispatcher
    def func(testdata: xr.DataArray | xr.Dataset) -> xr.DataArray | xr.Dataset:
        return xr.DataArray

    assert isinstance(func(ds3.to_dataset(name="test")), xr.Dataset)


def test_generate_datatree_dispatcher():
    # see `test_calc_linregress_spatial_datatree` in test_core_stat.py
    pass


def test_transfer_xarray_lon_from180TO360():
    result_data = ecl.utility.transfer_xarray_lon_from180TO360(
        ds4, lon_dim="lon"
    ).lon.data
    refer_data = np.array([120, 182, 310])
    assert np.isclose(result_data, refer_data).all()


def test_transfer_xarray_lon_from360TO180():
    result_data = ecl.utility.transfer_xarray_lon_from360TO180(
        ds5, lon_dim="lon"
    ).lon.data
    refer_data = np.array([-180, -50, 120])
    assert np.isclose(result_data, refer_data).all()


def test_module_available():
    result_data1 = ecl.utility.module_available("xarray")
    result_data2 = ecl.utility.module_available("This is sample text!&&*+")
    refer_data1 = True
    refer_data2 = False
    assert result_data1 is refer_data1
    assert result_data2 is refer_data2


def test_reverse_bool_xarraydata():
    data = xr.DataArray(np.array([True, False]), dims="lon")
    result_data = ecl.utility.reverse_bool_xarraydata(data)
    refer_data = xr.DataArray(np.array([False, True]), dims="lon")
    assert (result_data == refer_data).all().data


def test_compare_two_dataarray_coordinate1():
    ecl.utility.compare_two_dataarray_coordinate(
        data_input1=comparedata1,
        data_input2=comparedata1,
    )
    assert 1


@pytest.mark.xfail(raises=AssertionError)
def test_compare_two_dataarray_coordinate2():
    ecl.utility.compare_two_dataarray_coordinate(
        data_input1=comparedata1,
        data_input2=comparedata2,
    )


@pytest.mark.xfail(raises=AssertionError)
def test_compare_two_dataarray_coordinate3():
    ecl.utility.compare_two_dataarray_coordinate(
        data_input1=comparedata1,
        data_input2=comparedata3,
    )


@pytest.mark.xfail(raises=AssertionError)
def test_compare_two_dataarray_coordinate4():
    ecl.utility.compare_two_dataarray_coordinate(
        data_input1=comparedata1,
        data_input2=comparedata4,
    )


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_compare_two_dataarray_coordinate5():
    ecl.utility.compare_two_dataarray_coordinate(
        data_input1=comparedata1,
        data_input2=comparedata5,
    )


def test_compare_multi_dataarray_coordinate():
    ecl.utility.compare_multi_dataarray_coordinate(
        data_input_list=[comparedata1, comparedata1, comparedata1]
    )
