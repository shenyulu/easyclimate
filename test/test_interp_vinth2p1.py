"""
pytest for interp/vinth2p.py

Part2
"""

import pytest

import pytest
import xarray as xr
import numpy as np

from easyclimate.interp.vinth2p import interp_vintp2p_ecmwf


class TestInterpVintp2pECMWF:

    @pytest.fixture
    def example_data(self):
        time = 1
        lev = 3
        lat = 4
        lon = 5
        shape_4d = (time, lev, lat, lon)
        shape_3d = (time, lat, lon)

        data_input = xr.DataArray(
            np.linspace(250, 300, np.prod(shape_4d)).reshape(shape_4d),
            dims=["time", "lev", "lat", "lon"],
            coords={
                "time": [0],
                "lev": [1000, 850, 700],
                "lat": np.linspace(-90, 90, lat),
                "lon": np.linspace(0, 360, lon, endpoint=False),
            },
        )

        pressure_data = xr.DataArray(
            np.tile(
                np.array([1000, 850, 700]).reshape((1, 3, 1, 1)), (time, 1, lat, lon)
            ),
            dims=["time", "lev", "lat", "lon"],
            coords=data_input.coords,
        )

        surface_pressure_data = xr.DataArray(
            np.full(shape_3d, 101325),
            dims=["time", "lat", "lon"],
            coords={k: v for k, v in data_input.coords.items() if k != "lev"},
        )

        temperature_bottom_data = xr.DataArray(
            np.full(shape_3d, 290),
            dims=["time", "lat", "lon"],
            coords=surface_pressure_data.coords,
        )

        surface_geopotential_data = xr.DataArray(
            np.full(shape_3d, 100),
            dims=["time", "lat", "lon"],
            coords=surface_pressure_data.coords,
        )

        return {
            "data_input": data_input,
            "pressure_data": pressure_data,
            "pressure_data_units": "hPa",
            "surface_pressure_data": surface_pressure_data,
            "surface_pressure_data_units": "Pa",
            "temperature_bottom_data": temperature_bottom_data,
            "surface_geopotential_data": surface_geopotential_data,
            "vertical_output_level": [925, 850, 700],
            "vertical_input_dim": "lev",
            "vertical_output_dim": "plev",
            "vertical_output_dim_units": "hPa",
            "lat_dim": "lat",
            "lon_dim": "lon",
        }

    def test_interp_linear_no_extrapolation(self, example_data):
        output = interp_vintp2p_ecmwf(
            **example_data,
            variable_flag="other",
            interp_method="linear",
            extrapolation=False,
        )
        assert isinstance(output, xr.DataArray)
        assert "plev" in output.dims
        assert not output.isnull().all(), "Output should contain valid data"

    def test_interp_invalid_method_raises(self, example_data):
        with pytest.raises(
            ValueError, match="interp_method must be 'linear', 'log', or 'loglog'"
        ):
            interp_vintp2p_ecmwf(
                **example_data,
                variable_flag="T",
                interp_method="invalid",
                extrapolation=False,
            )

    def test_missing_extrapolation_inputs(self, example_data):
        example_data.pop("temperature_bottom_data")
        example_data.pop("surface_geopotential_data")
        with pytest.raises(
            ValueError, match="must be provided when extrapolation=True"
        ):
            interp_vintp2p_ecmwf(
                **example_data,
                variable_flag="Z",
                extrapolation=True,
            )
