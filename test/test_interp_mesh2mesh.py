"""
pytest for interp.interp_mesh2mesh.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import TEST_DATA_PATH

data_u = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_interp_mesh2mesh.nc"))
).uwnd
data_output = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_output_interp_mesh2mesh.nc"))
)
data_output_linear = data_output["linear"]
data_output_nearest = data_output["nearest"]
data_output_cubic = data_output["cubic"]
data_output_conservative = data_output["conservative"]


def test_interp_mesh2mesh1():
    target_grid = xr.DataArray(
        dims=("lat", "lon"),
        coords={"lon": np.linspace(0, 357.5, 72), "lat": np.linspace(-90, 90, 36)},
    )
    result_data1 = ecl.interp.interp_mesh2mesh(
        data_u, target_grid, lon_dim="lon", lat_dim="lat", method="linear"
    ).data.flatten()
    result_data2 = ecl.interp.interp_mesh2mesh(
        data_u, target_grid, lon_dim="lon", lat_dim="lat", method="nearest"
    ).data.flatten()
    result_data3 = ecl.interp.interp_mesh2mesh(
        data_u, target_grid, lon_dim="lon", lat_dim="lat", method="cubic"
    ).data.flatten()
    result_data4 = ecl.interp.interp_mesh2mesh(
        data_u, target_grid, lon_dim="lon", lat_dim="lat", method="conservative"
    ).data.flatten()

    refer_data1 = data_output_linear.data.flatten()
    refer_data2 = data_output_nearest.data.flatten()
    refer_data3 = data_output_cubic.data.flatten()
    refer_data4 = data_output_conservative.data.flatten()

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()


def test_interp_mesh2mesh2():
    target_grid = xr.DataArray(
        dims=("lat", "lon", "time"),
        coords={
            "lon": np.linspace(0, 357.5, 72),
            "lat": np.linspace(-90, 90, 36),
            "time": np.array([1, 2]),
        },
    )
    with pytest.raises(ValueError):
        ecl.interp.interp_mesh2mesh(
            data_u, target_grid, lon_dim="lon", lat_dim="lat", method="linear"
        )
        assert 1 == 1


def test_interp_mesh2mesh3():
    target_grid = xr.DataArray(
        dims=("lat1", "lon1"),
        coords={"lon1": np.linspace(0, 357.5, 72), "lat1": np.linspace(-90, 90, 36)},
    )
    with pytest.raises(ValueError):
        ecl.interp.interp_mesh2mesh(
            data_u, target_grid, lon_dim="lon", lat_dim="lat", method="linear"
        )
        assert 1 == 1


def test_interp_mesh2mesh4():
    target_grid = xr.DataArray(
        dims=("lat", "lon"),
        coords={"lon": np.linspace(0, 357.5, 72), "lat": np.linspace(-90, 90, 36)},
    )
    target_grid = ecl.utility.transfer_xarray_lon_from360TO180(target_grid)

    result_data = ecl.interp.interp_mesh2mesh(
        data_u, target_grid, lon_dim="lon", lat_dim="lat", method="linear"
    ).data.flatten()
    refer_data = data_output_linear.data.flatten()
    assert np.isclose(result_data, refer_data).all()


def test_interp_mesh2mesh4():
    target_grid = xr.DataArray(
        dims=("lat", "lon"),
        coords={"lon": np.linspace(0, 357.5, 72), "lat": np.linspace(-90, 90, 36)},
    )
    data_u1 = ecl.utility.transfer_xarray_lon_from360TO180(data_u)

    result_data = ecl.interp.interp_mesh2mesh(
        data_u1, target_grid, lon_dim="lon", lat_dim="lat", method="linear"
    ).data.flatten()
    refer_data = data_output_linear.data.flatten()
    assert np.isclose(result_data, refer_data).all()
