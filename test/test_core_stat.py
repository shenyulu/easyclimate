"""
pytest for stat.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
import pandas as pd
import dask.array as da
from .util import round_sf_np_new
from xarray import DataTree

sst_data = ecl.tutorial.open_tutorial_dataset("mini_HadISST_sst").sst
sic_data_Barents_Sea = ecl.tutorial.open_tutorial_dataset("mini_HadISST_ice").sic
sic_data_Barents_Sea_12 = ecl.get_specific_months_data(sic_data_Barents_Sea, 12)


def test_calc_linregress_spatial():
    result_data = ecl.calc_linregress_spatial(
        sic_data_Barents_Sea_12.sel(lon=slice(34.5, 36.5), lat=slice(78.5, 80.5)),
        dim="time",
    ).compute()
    result_data1 = result_data["slope"].data
    result_data2 = result_data["intercept"].data
    result_data3 = result_data["rvalue"].data
    result_data4 = result_data["pvalue"].data
    result_data5 = result_data["stderr"].data
    result_data6 = result_data["intercept_stderr"].data

    refer_data1 = np.array(
        [
            [-0.01689814, -0.01618345, -0.01640629],
            [-0.01011993, -0.00922373, -0.0091192],
            [-0.00641115, -0.0054169, -0.00600519],
        ]
    )
    refer_data2 = np.array(
        [
            [1.05593576, 1.04461794, 1.04132889],
            [1.01817285, 1.01218153, 1.01741975],
            [0.89857147, 0.91509416, 0.94763013],
        ]
    )
    refer_data3 = np.array(
        [
            [-0.58339207, -0.57217978, -0.57376992],
            [-0.47457306, -0.4485609, -0.45343254],
            [-0.32601495, -0.28470031, -0.33127693],
        ]
    )
    refer_data4 = np.array(
        [
            [5.01586941e-05, 7.52887062e-05, 7.11385575e-05],
            [1.49647207e-03, 2.88846392e-03, 2.56361219e-03],
            [3.51172733e-02, 6.76380713e-02, 3.21088821e-02],
        ]
    )
    refer_data5 = np.array(
        [
            [0.00371969, 0.00366767, 0.00370284],
            [0.00296779, 0.00290584, 0.00283422],
            [0.00293946, 0.00288389, 0.00270435],
        ]
    )
    refer_data6 = np.array(
        [
            [0.08858534, 0.08734656, 0.08818416],
            [0.07067876, 0.06920341, 0.06749764],
            [0.07000405, 0.06868049, 0.06440477],
        ]
    )

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()
    assert np.isclose(result_data6, refer_data6).all()


def test_calc_linregress_spatial_datatree():
    data = sic_data_Barents_Sea_12.sel(
        lon=slice(34.5, 36.5), lat=slice(78.5, 80.5)
    ).to_dataset(name="sic")
    result_data = ecl.calc_linregress_spatial(data, dim="time")
    assert isinstance(result_data, DataTree)


def test_chunk_handling_with_dask():
    """Test that function properly handles chunked dask arrays"""
    # Create a dask-backed DataArray
    data = da.random.random((10, 5, 5), chunks=(2, 5, 5))
    data_input = xr.DataArray(data, dims=["time", "lat", "lon"])

    # Call function - should rechunk time dimension to -1
    result = ecl.calc_linregress_spatial(data_input, dim="time")

    # Verify the result is correct type and contains expected variables
    assert isinstance(result, xr.Dataset)
    assert all(
        var in result.data_vars
        for var in [
            "slope",
            "intercept",
            "rvalue",
            "pvalue",
            "stderr",
            "intercept_stderr",
        ]
    )


def test_warning_for_numpy_array_x():
    """Test that warning is raised when passing numpy array for x"""
    data = np.random.random((10, 5, 5))
    data_input = xr.DataArray(data, dims=["time", "lat", "lon"])
    x = np.arange(10)

    with pytest.warns(
        UserWarning, match="Assuming that the coordinate value of 'time' in"
    ):
        ecl.calc_linregress_spatial(data_input, dim="time", x=x)


def test_error_for_mismatched_dim_sizes():
    """Test ValueError when x size doesn't match data_input dim size"""
    data = np.random.random((10, 5, 5))
    data_input = xr.DataArray(data, dims=["time", "lat", "lon"])
    x = np.arange(8)  # Wrong size

    with pytest.raises(ValueError, match="`data_input` array size along dimension"):
        ecl.calc_linregress_spatial(data_input, dim="time", x=x)


def test_error_for_mismatched_xarray_dims():
    """Test ValueError when x DataArray has wrong dimension name"""
    data = np.random.random((10, 5, 5))
    data_input = xr.DataArray(data, dims=["time", "lat", "lon"])
    x = xr.DataArray(np.arange(10), dims=["wrong_dim"])

    with pytest.raises(ValueError, match="The coordinate name of `data_input` array"):
        ecl.calc_linregress_spatial(data_input, dim="time", x=x)


def test_error_for_mismatched_xarray_coords():
    """Test ValueError when x DataArray has mismatched coordinates"""
    data = np.random.random((10, 5, 5))
    time_coords = np.arange(10)
    data_input = xr.DataArray(
        data, dims=["time", "lat", "lon"], coords={"time": time_coords}
    )
    x = xr.DataArray(
        np.arange(10), dims=["time"], coords={"time": time_coords + 1}
    )  # Offset coords

    with pytest.raises(ValueError, match="Coordinate value of 'time' in"):
        ecl.calc_linregress_spatial(data_input, dim="time", x=x)


def test_calc_detrend_spatial():
    result_data = (
        ecl.calc_detrend_spatial(
            sic_data_Barents_Sea_12.sel(lon=slice(34.5, 36.5), lat=slice(78.5, 80.5)),
            time_dim="time",
        )
        .mean(dim=("lat", "lon"))
        .data
    )
    result_data = round_sf_np_new(result_data)
    refer_data = np.array(
        [
            -5.9e-02,
            -4.9e-02,
            -4.5e-02,
            -4.9e-01,
            -2.2e-02,
            -4.7e-03,
            2.0e-02,
            3.3e-02,
            1.3e-02,
            2.9e-02,
            -3.5e-04,
            2.5e-02,
            5.3e-02,
            3.0e-02,
            5.6e-02,
            1.5e-01,
            5.8e-02,
            1.5e-01,
            1.5e-01,
            1.5e-01,
            2.2e-02,
            1.2e-01,
            1.8e-01,
            1.5e-01,
            1.2e-01,
            -3.3e-01,
            -9.8e-02,
            1.6e-01,
            -1.2e-01,
            2.9e-01,
            1.3e-01,
            -4.2e-01,
            1.8e-01,
            2.7e-01,
            -2.0e-01,
            -6.2e-01,
            -2.3e-01,
            -6.0e-01,
            3.9e-01,
            5.8e-02,
            3.9e-01,
            -7.8e-02,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calculate_corr_pvalue():
    rng = np.random.default_rng(seed=42)

    times = pd.date_range("2000-01-01", periods=5)
    lats = [10, 20, 30]
    lons = [100, 110, 120, 130]
    data = rng.random((5, 3, 4))
    da = xr.DataArray(
        data,
        dims=["time", "lat", "lon"],
        coords={"time": times, "lat": lats, "lon": lons},
    )
    x = xr.DataArray(
        np.array([0.0276367, 0.67968458, 0.61436011, 0.587718, 0.88766532]),
        dims=["time"],
        coords={"time": times},
    )
    result = ecl.calc_corr_spatial(da, x)
    result_data1 = result["corr"].data.flatten()
    result_data2 = result["pvalue"].data.flatten()
    refer_data1 = np.array(
        [
            -0.18430726,
            -0.18548118,
            -0.9047632,
            -0.74794148,
            0.80681898,
            -0.52102108,
            -0.19712648,
            -0.07580053,
            0.57095664,
            0.1390009,
            0.00566705,
            -0.65137758,
        ]
    )
    refer_data2 = np.array(
        [
            0.76666814,
            0.76519923,
            0.03477264,
            0.14603029,
            0.09891851,
            0.36798864,
            0.75064591,
            0.90358027,
            0.31473309,
            0.82359014,
            0.99278452,
            0.23374293,
        ]
    )
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()


def test_time_dimension_mismatch():
    """Test that ValueError is raised when time dimensions don't match"""
    # Create test data with different time lengths
    data_input = xr.DataArray(
        np.random.rand(10, 2, 2),
        dims=["time", "lat", "lon"],
        coords={"time": range(10)},
    )
    x = xr.DataArray(np.random.rand(8), dims=["time"])

    with pytest.raises(ValueError) as excinfo:
        ecl.calc_corr_spatial(data_input, x)
    assert "time dimension is not consistent" in str(excinfo.value)


def test_fewer_than_2_valid_points():
    """Test that nan is returned when fewer than 2 valid points exist"""
    # Create test data with only 1 valid point
    data_input = xr.DataArray(
        np.array(
            [
                [[1.0, np.nan], [np.nan, np.nan]],  # Only one valid point at [0,0,0]
                [[np.nan, np.nan], [np.nan, np.nan]],
            ]
        ),
        dims=["time", "lat", "lon"],
    )
    x = xr.DataArray([1.0, np.nan], dims=["time"])

    result = ecl.calc_corr_spatial(data_input, x)

    # Check all outputs are nan since we don't have enough valid points
    assert np.isnan(result["corr"].values).all()
    assert np.isnan(result["pvalue"].values).all()


def test_exactly_2_valid_points():
    """Test that correlation is calculated when exactly 2 valid points exist"""
    # Create test data with exactly 2 valid points
    data_input = xr.DataArray(
        np.array([[[1.0, 2.0], [np.nan, np.nan]], [[3.0, 4.0], [np.nan, np.nan]]]),
        dims=["time", "lat", "lon"],
    )
    x = xr.DataArray([1.0, 2.0], dims=["time"])

    result = ecl.calc_corr_spatial(data_input, x)

    # Check correlation is calculated for first two points (perfect correlation)
    assert result["corr"].values[0, 0] == pytest.approx(1.0)
    assert result["corr"].values[0, 1] == pytest.approx(1.0)
    # Check pvalue is 0 for perfect correlation
    assert result["pvalue"].values[0, 0] == pytest.approx(1.0)
    assert result["pvalue"].values[0, 1] == pytest.approx(1.0)
    # Other points should be nan
    assert np.isnan(result["corr"].values[1, 0])
    assert np.isnan(result["corr"].values[1, 1])
    assert np.isnan(result["pvalue"].values[1, 0])
    assert np.isnan(result["pvalue"].values[1, 1])


def test_mixed_valid_and_invalid_points():
    """Test behavior with some valid and some invalid grid points"""
    data_input = xr.DataArray(
        np.array(
            [
                [[1.0, 2.0, np.nan], [3.0, np.nan, np.nan]],
                [[4.0, 5.0, np.nan], [6.0, np.nan, np.nan]],
            ]
        ),
        dims=["time", "lat", "lon"],
    )
    x = xr.DataArray([1.0, 2.0], dims=["time"])

    result = ecl.calc_corr_spatial(data_input, x)

    # Check correlations are calculated where there are enough points
    assert not np.isnan(result["corr"].values[0, 0])
    assert not np.isnan(result["corr"].values[0, 1])
    assert not np.isnan(result["pvalue"].values[0, 0])
    assert not np.isnan(result["pvalue"].values[0, 1])

    # Check nan where not enough points
    assert np.isnan(result["corr"].values[0, 2])
    assert np.isnan(result["corr"].values[1, 1])
    assert np.isnan(result["corr"].values[1, 2])
    assert np.isnan(result["pvalue"].values[0, 2])
    assert np.isnan(result["pvalue"].values[1, 1])
    assert np.isnan(result["pvalue"].values[1, 2])


def test_calc_ttestSpatialPattern_spatial():
    sic_data_Barents_Sea_12_spatial_mean = sic_data_Barents_Sea_12.sel(
        lon=slice(34.5, 36.5), lat=slice(78.5, 80.5)
    ).mean(dim=("lat", "lon"))
    test1 = sic_data_Barents_Sea_12_spatial_mean.isel(time=slice(0, 20))
    test2 = sic_data_Barents_Sea_12_spatial_mean.isel(time=slice(21, None))
    result_data = ecl.calc_ttestSpatialPattern_spatial(test1, test2)
    result_data1 = result_data["statistic"].data
    result_data2 = result_data["pvalue"].data
    refer_data1 = 3.43382768
    refer_data2 = 0.00142438
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()


def test_calc_levenetestSpatialPattern_spatial():
    sic_data_Barents_Sea_12_spatial_mean = sic_data_Barents_Sea_12.sel(
        lon=slice(34.5, 36.5), lat=slice(78.5, 80.5)
    ).mean(dim=("lat", "lon"))
    test1 = sic_data_Barents_Sea_12_spatial_mean.isel(time=slice(0, 20))
    test2 = sic_data_Barents_Sea_12_spatial_mean.isel(time=slice(21, None))
    result_data = ecl.calc_levenetestSpatialPattern_spatial(test1, test2)
    result_data1 = result_data["statistic"].data
    result_data2 = result_data["pvalue"].data
    refer_data1 = 13.79602081
    refer_data2 = 0.00063671
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()


def test_calc_skewness_spatial():
    sic_data_Barents_Sea_12_detrend = ecl.calc_detrend_spatial(
        sic_data_Barents_Sea_12.sel(lon=slice(34.5, 36.5), lat=slice(78.5, 80.5)),
        time_dim="time",
    )
    result_data = ecl.calc_skewness_spatial(sic_data_Barents_Sea_12_detrend, dim="time")
    result_data1 = result_data["skewness"].data
    result_data2 = result_data["pvalue"].data
    refer_data1 = np.array(
        [
            [-0.24333526, -0.32189173, -0.27430525],
            [-1.3397094, -1.5588326, -1.6165946],
            [-1.8677251, -2.209491, -2.330299],
        ]
    )
    refer_data2 = np.array(
        [
            [7.70400089e-01, 6.26686378e-01, 7.18524740e-01],
            [4.16622922e-04, 3.66448314e-05, 1.56112432e-05],
            [1.88511919e-06, 6.86471937e-08, 1.30767304e-08],
        ]
    )
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()


def test_calc_kurtosis_spatial():
    sic_data_Barents_Sea_12_detrend = ecl.calc_detrend_spatial(
        sic_data_Barents_Sea_12.sel(lon=slice(34.5, 36.5), lat=slice(78.5, 80.5)),
        time_dim="time",
    )
    result_data = ecl.calc_kurtosis_spatial(sic_data_Barents_Sea_12_detrend, dim="time")
    refer_data = np.array(
        [
            [2.7300231, 2.8368442, 2.7555804],
            [4.667509, 5.464381, 5.85138],
            [6.32586, 7.463421, 8.428185],
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_theilslopes_spatial():
    result_data = ecl.calc_theilslopes_spatial(
        sic_data_Barents_Sea_12.sel(lon=slice(34.5, 36.5), lat=slice(78.5, 80.5)),
        dim="time",
    ).compute()
    result_data1 = result_data["slope"].data.flatten()
    result_data2 = result_data["intercept"].data.flatten()
    result_data3 = result_data["low_slope"].data.flatten()
    result_data4 = result_data["high_slope"].data.flatten()

    refer_data1 = np.array(
        [
            -0.00999999,
            -0.00999999,
            -0.01,
            -0.0035,
            -0.00296296,
            -0.00258064,
            -0.0025,
            -0.00179487,
            -0.001875,
        ]
    )
    refer_data2 = np.array(
        [
            1.10499978,
            1.09499979,
            1.08499995,
            1.00675,
            0.99574073,
            0.99290321,
            0.91125007,
            0.92179486,
            0.93343744,
        ]
    )
    refer_data3 = np.array(
        [
            -0.01833333,
            -0.018,
            -0.019,
            -0.00853659,
            -0.0075,
            -0.00736842,
            -0.005625,
            -0.00416667,
            -0.00333334,
        ]
    )
    refer_data4 = np.array(
        [
            -0.003,
            -0.00285714,
            -0.00333333,
            -0.00090909,
            -0.00058823,
            -0.00058823,
            0.0,
            0.0,
            -0.000625,
        ]
    )

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()


def test_chunk_handling_with_dask_theilslopes_spatial():
    """Test that function properly handles chunked dask arrays"""
    # Create a dask-backed DataArray
    data = da.random.random((10, 5, 5), chunks=(2, 5, 5))
    data_input = xr.DataArray(data, dims=["time", "lat", "lon"])

    # Call function - should rechunk time dimension to -1
    result = ecl.calc_theilslopes_spatial(data_input, dim="time")

    # Verify the result is correct type and contains expected variables
    assert isinstance(result, xr.Dataset)
    assert all(
        var in result.data_vars
        for var in ["slope", "intercept", "low_slope", "high_slope"]
    )


def test_warning_for_numpy_array_x_theilslopes_spatial():
    """Test that warning is raised when passing numpy array for x"""
    data = np.random.random((10, 5, 5))
    data_input = xr.DataArray(data, dims=["time", "lat", "lon"])
    x = np.arange(10)

    with pytest.warns(
        UserWarning, match="Assuming that the coordinate value of 'time' in"
    ):
        ecl.calc_theilslopes_spatial(data_input, dim="time", x=x)


def test_error_for_mismatched_dim_sizes_theilslopes_spatial():
    """Test ValueError when x size doesn't match data_input dim size"""
    data = np.random.random((10, 5, 5))
    data_input = xr.DataArray(data, dims=["time", "lat", "lon"])
    x = np.arange(8)  # Wrong size

    with pytest.raises(ValueError, match="`data_input` array size along dimension"):
        ecl.calc_theilslopes_spatial(data_input, dim="time", x=x)


def test_error_for_mismatched_xarray_dims_theilslopes_spatial():
    """Test ValueError when x DataArray has wrong dimension name"""
    data = np.random.random((10, 5, 5))
    data_input = xr.DataArray(data, dims=["time", "lat", "lon"])
    x = xr.DataArray(np.arange(10), dims=["wrong_dim"])

    with pytest.raises(ValueError, match="The coordinate name of `data_input` array"):
        ecl.calc_theilslopes_spatial(data_input, dim="time", x=x)


def test_error_for_mismatched_xarray_coords_theilslopes_spatial():
    """Test ValueError when x DataArray has mismatched coordinates"""
    data = np.random.random((10, 5, 5))
    time_coords = np.arange(10)
    data_input = xr.DataArray(
        data, dims=["time", "lat", "lon"], coords={"time": time_coords}
    )
    x = xr.DataArray(
        np.arange(10), dims=["time"], coords={"time": time_coords + 1}
    )  # Offset coords

    with pytest.raises(ValueError, match="Coordinate value of 'time' in"):
        ecl.calc_theilslopes_spatial(data_input, dim="time", x=x)
