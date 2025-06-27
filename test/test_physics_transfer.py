"""
pytest for physics.transfer
"""

import pytest

import numpy as np
import xarray as xr
import pandas as pd

from easyclimate.physics import (
    transfer_mixing_ratio_2_specific_humidity,
    transfer_specific_humidity_2_mixing_ratio,
    transfer_dewpoint_2_specific_humidity,
    transfer_specific_humidity_2_dewpoint,
    transfer_dewpoint_2_relative_humidity,
    transfer_mixing_ratio_2_relative_humidity,
    transfer_specific_humidity_2_relative_humidity,
    transfer_relative_humidity_2_dewpoint,
)
from .util import round_sf_np


def test_transfer_mixing_ratio_2_specific_humidity():
    result_data = transfer_mixing_ratio_2_specific_humidity(
        mixing_ratio_data=xr.DataArray([0.01594761, 0.01923578])
    ).data
    refer_data = np.array([0.01569728, 0.01887275])
    assert np.isclose(result_data, refer_data).all()


def test_transfer_specific_humidity_2_mixing_ratio():
    result_data = transfer_specific_humidity_2_mixing_ratio(
        specific_humidity_data=xr.DataArray([0.01569728, 0.01887275]),
        specific_humidity_data_units="g/g",
    ).data
    refer_data = np.array([0.01594761, 0.01923578])
    assert np.isclose(result_data, refer_data).all()


def test_transfer_dewpoint_2_specific_humidity():
    result_data = transfer_dewpoint_2_specific_humidity(
        dewpoint_data=xr.DataArray([15, 18]),
        pressure_data=xr.DataArray([988, 980]),
        dewpoint_data_units="celsius",
        pressure_data_units="hPa",
    ).data
    refer_data = np.array([0.01079758, 0.01319518])
    assert np.isclose(result_data, refer_data).all()


def test_transfer_specific_humidity_2_dewpoint():
    result_data = transfer_specific_humidity_2_dewpoint(
        specific_humidity_data=xr.DataArray([0.01079, 0.01319]),
        pressure_data=xr.DataArray([988, 980]),
        specific_humidity_data_units="g/g",
        pressure_data_units="hPa",
    ).data
    refer_data = np.array([14.98916116, 17.9938138])
    assert np.isclose(result_data, refer_data).all()


def test_transfer_dewpoint_2_relative_humidity():
    result_data = transfer_dewpoint_2_relative_humidity(
        temperature_data=xr.DataArray([25, 28]),
        dewpoint_data=xr.DataArray([12, 13]),
        temperature_data_units="celsius",
        dewpoint_data_units="celsius",
    ).data
    refer_data = np.array([0.44248477, 0.3958322])
    assert np.isclose(result_data, refer_data).all()


def test_transfer_mixing_ratio_2_relative_humidity():
    result_data = transfer_mixing_ratio_2_relative_humidity(
        pressure_data=xr.DataArray([1013.25, 1014]),
        temperature_data=xr.DataArray([30, 28]),
        mixing_ratio_data=xr.DataArray([18 / 1000, 16 / 1000]),
        pressure_data_units="hPa",
        temperature_data_units="celsius",
    ).data
    refer_data = np.array([0.67127708, 0.67260414])
    assert np.isclose(result_data, refer_data).all()


def test_transfer_specific_humidity_2_relative_humidity():
    result_data = transfer_specific_humidity_2_relative_humidity(
        pressure_data=xr.DataArray([1013.25, 1014]),
        temperature_data=xr.DataArray([30, 28]),
        specific_humidity_data=xr.DataArray([18 / 1000, 16 / 1000]),
        pressure_data_units="hPa",
        temperature_data_units="degC",
        specific_humidity_data_units="g/g",
    ).data
    refer_data = np.array([0.6832293, 0.68326215])
    assert np.isclose(result_data, refer_data).all()


def test_transfer_relative_humidity_2_dewpoint():
    # Create sample xarray DataArrays
    temperature = xr.DataArray(
        np.linspace(280, 300, 5),
        dims=["time"],
        coords={"time": pd.date_range("2023-01-01", periods=5)},
        attrs={"units": "K"},
    )
    relative_humidity = xr.DataArray(
        np.linspace(30, 70, 5),
        dims=["time"],
        coords={"time": pd.date_range("2023-01-01", periods=5)},
        attrs={"units": "%"},
    )
    result_data = transfer_relative_humidity_2_dewpoint(
        temperature,
        relative_humidity,
        relative_humidity_data_units="%",
        temperature_data_units="K",
    )
    refer_data = np.array(
        [30.13948923, 40.25464736, 50.40840651, 60.60341513, 70.84238873]
    )
    assert np.isclose(result_data, refer_data).all()
