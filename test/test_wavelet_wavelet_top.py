"""
pytest for wavelet_top.py
"""
import pytest

import easyclimate as ecl
import xarray as xr
import numpy as np
import os
from .const_define import TEST_DATA_PATH

def test_timeseries_wavelet_transform():
    inputdata = xr.open_dataset(os.path.join(TEST_DATA_PATH, 'test_input_waveletpytest.nc')).sst
    refer_data = xr.open_dataset(os.path.join(TEST_DATA_PATH, 'test_output_waveletpytest.nc'))

    result_data = ecl.wavelet.timeseries_wavelet_transform(inputdata, dt = 0.25)
    result_data1 = result_data.power.data.flatten()
    result_data2 = result_data.sig.data.flatten()
    result_data3 = result_data.global_ws.data.flatten()
    result_data4 = result_data.global_signif.data.flatten()
    result_data5 = result_data.coi.data.flatten()
    result_data6 = result_data.coi_bottom.data.flatten()

    refer_data1 = refer_data.power.data.flatten()
    refer_data2 = refer_data.sig.data.flatten()
    refer_data3 = refer_data.global_ws.data.flatten()
    refer_data4 = refer_data.global_signif.data.flatten()
    refer_data5 = refer_data.coi.data.flatten()
    refer_data6 = refer_data.coi_bottom.data.flatten()

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()
    assert np.isclose(result_data6, refer_data6).all()