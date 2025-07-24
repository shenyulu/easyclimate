"""
pytest for core/utility.py

Part2
"""

import pytest
import xarray as xr
import numpy as np
import pandas as pd
from easyclimate.core.datanode import DataNode
from easyclimate.core.mk_test import sens_slope
from easyclimate.core.utility import (
    get_compress_xarraydata,
    assert_compared_version,
    generate_dataset_dispatcher,
    datetime_to_numeric,
    numeric_to_datetime,
    calculate_time_steps,
    clean_extra_coords,
)


@pytest.fixture
def sample_data_array():
    return xr.DataArray(
        np.random.rand(10, 10),
        dims=["x", "y"],
        coords={"x": np.arange(10), "y": np.arange(10)},
    )


@pytest.fixture
def sample_dataset():
    return xr.Dataset(
        {
            "var1": (["x", "y"], np.random.rand(10, 10)),
            "var2": (["x", "y"], np.random.rand(10, 10)),
        }
    )


@pytest.fixture
def sample_time_series():
    return xr.DataArray(
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


@pytest.fixture
def datetime_array():
    return np.array(["2023-01-01", "2023-01-02", "2023-01-03"], dtype="datetime64[ns]")


class TestGetCompressXarrayData:
    def test_dataarray_compression(self, sample_data_array):
        result = get_compress_xarraydata(sample_data_array, complevel=5)
        assert result.encoding.get("zlib") == True
        assert result.encoding.get("complevel") == 5
        assert isinstance(result, xr.DataArray)

    def test_dataset_compression(self, sample_dataset):
        result = get_compress_xarraydata(sample_dataset, complevel=3)
        assert isinstance(result, xr.Dataset)
        for var in result:
            assert result[var].encoding.get("zlib") == True
            assert result[var].encoding.get("complevel") == 3


class TestAssertComparedVersion:
    @pytest.mark.parametrize(
        "ver1, ver2, expected",
        [
            ("1.0", "2.0", -1),
            ("2.0", "1.0", 1),
            ("1.0", "1.0", 0),
            ("1.2.3", "1.2", 1),
            ("1.2", "1.2.3", -1),
        ],
    )
    def test_version_comparison(self, ver1, ver2, expected):
        assert assert_compared_version(ver1, ver2) == expected


class TestGenerateDatasetDispatcher:
    def test_dataarray_processing(self, sample_data_array):
        @generate_dataset_dispatcher
        def test_func(data):
            return data * 2

        result = test_func(sample_data_array)
        assert isinstance(result, xr.DataArray)
        assert np.allclose(result, sample_data_array * 2)

    def test_dataset_processing(self, sample_dataset):
        @generate_dataset_dispatcher
        def test_func(data):
            return data * 2

        result = test_func(sample_dataset)
        assert isinstance(result, xr.Dataset)
        assert set(result.data_vars) == set(sample_dataset.data_vars)
        for var in result:
            assert np.allclose(result[var], sample_dataset[var] * 2)

    def test_unsupported_type(self):
        @generate_dataset_dispatcher
        def test_func(data):
            return data

        with pytest.raises(ValueError):
            test_func(np.array([1, 2, 3]))


class TestGenerateDatanodeDispatcher:
    def test_dataset_returns_type(self, sample_time_series):
        result1 = sens_slope(
            sample_time_series.isel(lon=0).to_dataset(name="data1"),
            dim="time",
            returns_type="dataset_returns",
        )
        assert isinstance(result1, DataNode)

    def test_dataset_vars_type(self, sample_time_series):
        result2 = sens_slope(
            sample_time_series.isel(lon=0).to_dataset(name="data1"),
            dim="time",
            returns_type="dataset_vars",
        )
        assert isinstance(result2, DataNode)

    def test_invalid_returns_type(self, sample_time_series):
        with pytest.raises(ValueError):
            sens_slope(
                sample_time_series.isel(lon=0).to_dataset(name="data1"),
                dim="time",
                returns_type="invalid_type",
            )


class TestDatetimeConversion:
    def test_datetime_to_numeric(self, datetime_array):
        result = datetime_to_numeric(datetime_array, unit="D")
        assert isinstance(result, np.ndarray)
        assert result.dtype == np.int64
        assert len(result) == 3

    def test_numeric_to_datetime(self, datetime_array):
        numeric = datetime_to_numeric(datetime_array, unit="D")
        result = numeric_to_datetime(numeric, unit="D")
        assert np.all(result == datetime_array)

    def test_invalid_unit(self):
        with pytest.raises(ValueError):
            datetime_to_numeric(
                np.array(["2023-01-01"], dtype="datetime64[ns]"), unit="invalid"
            )


class TestCalculateTimeSteps:
    def test_time_steps(self, datetime_array):
        result = calculate_time_steps(datetime_array, unit="D")
        assert isinstance(result, np.ndarray)
        assert len(result) == 2
        assert np.all(result == 1)  # One day difference between consecutive dates

    def test_single_element_array(self):
        single_date = np.array(["2023-01-01"], dtype="datetime64[ns]")
        result = calculate_time_steps(single_date, unit="D")
        assert len(result) == 0


class TestCleanExtraCoords:
    @pytest.fixture
    def sample_dataarray(self):
        """Create a test DataArray with extra coordinates"""
        data = np.random.rand(3, 4)
        da = xr.DataArray(
            data,
            dims=["x", "y"],
            coords={
                "x": [1, 2, 3],
                "y": [10, 20, 30, 40],
                "extra": 1.1,
            },
        )
        return da

    def test_removes_extra_coords(self, sample_dataarray):
        """Test that extra coordinates are removed"""
        cleaned = clean_extra_coords(sample_dataarray)

        # Verify extra coordinates are gone
        assert "time" not in cleaned.coords
        assert "extra" not in cleaned.coords

        # Verify original coordinates remain
        assert "x" in cleaned.dims and "x" in cleaned.coords
        assert "y" in cleaned.dims and "y" in cleaned.coords

        # Verify data is unchanged
        assert np.isclose(cleaned.data, sample_dataarray.data).all()

    def test_no_extra_coords(self):
        """Test with DataArray that has no extra coordinates"""
        da = xr.DataArray([1, 2, 3], dims=["x"], coords={"x": [1, 2, 3]})
        cleaned = clean_extra_coords(da)
        xr.testing.assert_equal(cleaned, da)

    def test_empty_dataarray(self):
        """Test with empty DataArray"""
        da = xr.DataArray([])
        cleaned = clean_extra_coords(da)
        xr.testing.assert_equal(cleaned, da)
