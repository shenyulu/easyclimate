"""
pytest for interp.interp_mesh2point.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
from .const_define import DOCS_DATA_PATH
from .util import round_sf_np_new
from easyclimate.interp import interp_mesh2point, interp_mesh2point_withtime

t2m_data = ecl.open_tutorial_dataset("js_t2m_ERA5_2025052000").t2m

data = {
    "Site": [
        "Nanjing",
        "Suzhou",
        "Shanghai",
        "Chuzhou",
        "Changzhou",
        "Xuzhou",
        "Yancheng",
    ],
    "lon": [
        118.7788631,
        120.6212881,
        121.4700152,
        118.3139455,
        119.9691539,
        117.1810431,
        120.1577019,
    ],
    "lat": [
        32.0438284,
        31.311123,
        31.2312707,
        32.3027377,
        31.8122623,
        34.2665258,
        33.349559,
    ],
}
df = pd.DataFrame(data)


def test_interp_mesh2point():
    df_interp = interp_mesh2point(
        t2m_data,
        df,
        lon_dim_mesh="lon",
        lat_dim_mesh="lat",
        lon_dim_df="lon",
        lat_dim_df="lat",
    )
    result_data = df_interp["interpolated_value"].values
    refer_data = np.array(
        [
            27.43410485,
            27.05690739,
            25.86055861,
            27.13935815,
            27.48956649,
            28.62585815,
            26.3315014,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_interp_mesh2point_no_valid_stations():
    """
    Test that interp_mesh2point raises a ValueError when no stations
    fall within the grid's spatial extent.
    """
    # Create a small test grid (10x10 degrees)
    lats = np.linspace(0, 10, 11)
    lons = np.linspace(0, 10, 11)
    data = np.random.rand(11, 11)  # Random values for testing
    grid_data = xr.DataArray(
        data, dims=["lat", "lon"], coords={"lat": lats, "lon": lons}
    )

    # Create test points that are all outside the grid extent
    test_points = pd.DataFrame(
        {
            "lat": [15, 20, -5],  # All latitudes outside 0-10 range
            "lon": [5, 5, 5],  # Longitudes within range but latitudes invalid
        }
    )

    # Verify that the function raises ValueError with the expected message
    with pytest.raises(ValueError) as excinfo:
        interp_mesh2point(
            grid_data,
            test_points,
            lon_dim_mesh="lon",
            lat_dim_mesh="lat",
            lon_dim_df="lon",
            lat_dim_df="lat",
        )

    # Check the error message
    assert "No stations fall within the grid data extent" in str(excinfo.value)


def test_interp_mesh2point_edge_cases():
    """
    Test edge cases where points are just outside the grid bounds.
    """
    # Create a small test grid
    lats = np.linspace(0, 10, 11)
    lons = np.linspace(0, 10, 11)
    data = np.random.rand(11, 11)
    grid_data = xr.DataArray(
        data, dims=["lat", "lon"], coords={"lat": lats, "lon": lons}
    )

    # Test case 1: Points just outside grid bounds (should raise ValueError)
    edge_points = pd.DataFrame(
        {"lat": [10.0001, -0.0001], "lon": [5, 5]}  # Just outside bounds
    )

    with pytest.raises(ValueError):
        interp_mesh2point(
            grid_data,
            edge_points,
            lon_dim_mesh="lon",
            lat_dim_mesh="lat",
            lon_dim_df="lon",
            lat_dim_df="lat",
        )

    # Test case 2: Points exactly on grid bounds (should work)
    exact_edge_points = pd.DataFrame(
        {"lat": [10.0, 0.0], "lon": [5, 5]}  # Exactly on bounds
    )

    # This should not raise an exception
    result = interp_mesh2point(
        grid_data,
        exact_edge_points,
        lon_dim_mesh="lon",
        lat_dim_mesh="lat",
        lon_dim_df="lon",
        lat_dim_df="lat",
    )
    assert len(result) == 2
    assert not result["interpolated_value"].isna().any()


def test_interp_mesh2point_withtime():
    # Set the random number seed
    np.random.seed(42)
    # Create sample grid
    times = pd.date_range("2020-01-01", periods=5)
    lats = np.linspace(-45, -10, 100)
    lons = np.linspace(110, 156, 120)
    data = xr.DataArray(
        np.random.rand(5, 100, 120),
        dims=["time", "lat", "lon"],
        coords={"time": times, "lat": lats, "lon": lons},
    )
    # Create stations
    stations = pd.DataFrame(
        {
            "station_id_col": [1001, 1002, 1003],
            "lat_col": [-15.5, -20.3, 0.0],
            "lon_col": [125.5, 130.2, 112.0],
        }
    )
    # Interpolate
    result = interp_mesh2point_withtime(
        data,
        stations,
        lon_dim_df="lon_col",
        lat_dim_df="lat_col",
        station_dim_df="station_id_col",
    )
    result_data = result.interpolated_value.data.flatten()
    refer_data = np.array(
        [
            0.22494999,
            0.71982919,
            0.38652211,
            0.53725334,
            0.39874745,
            0.640323,
            0.65750595,
            0.51837553,
            0.51428686,
            0.78778669,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert np.isclose(result_data, refer_data, equal_nan=True).all()
