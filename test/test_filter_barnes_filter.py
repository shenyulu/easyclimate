"""
pytest for filter\barnes_filter.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import TEST_DATA_PATH

# sst_data = ecl.open_tutorial_dataset("sst_mnmean_oisst").sst.sortby('lat').sel(lon = slice(0, 20), lat = slice(-10, 10))
# sst_data.to_netcdf('data/test_input_sst_mnmean_oisst_barnes_filter.nc', format = "NETCDF3_64BIT")
sst_data = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_sst_mnmean_oisst_barnes_filter.nc"))
)["sst"]


def test_calc_barnes_lowpass():
    tmp = ecl.filter.barnes_filter.calc_barnes_lowpass(sst_data)
    result_data = tmp.isel(time=5).isel(lon=slice(5, 8), lat=slice(5, 8)).data.flatten()
    refer_data = np.array(
        [
            25.837467,
            25.59172,
            25.33062,
            25.847834,
            25.617031,
            25.377975,
            25.921867,
            25.711851,
            25.500973,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)


def test_calc_barnes_bandpass():
    tmp = ecl.filter.barnes_filter.calc_barnes_bandpass(sst_data)
    result_data = tmp.isel(time=5).isel(lon=slice(5, 8), lat=slice(5, 8)).data.flatten()
    refer_data = np.array(
        [
            0.36594316,
            0.26737976,
            0.10052491,
            -0.03947296,
            -0.15506744,
            -0.33001328,
            -0.4155304,
            -0.54795,
            -0.7124222,
        ]
    )
    assert np.allclose(result_data, refer_data, equal_nan=True)
