"""
pytest for field/teleconnection/index_ao_nam.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import DOCS_DATA_PATH

slp_data = ecl.open_tutorial_dataset("slp_monmean_NH")["slp"]


def test_calc_index_AO_EOF_Thompson_Wallace_1998_AND_calc_index_NAH_zonal_lat_Li_Wang_2003():
    slp_data_DJF_mean = ecl.calc_seasonal_mean(slp_data, extract_season="DJF")
    index_ao = ecl.field.teleconnection.calc_index_AO_EOF_Thompson_Wallace_1998(
        slp_data_DJF_mean, random_state=0
    )
    index_ao_point = ecl.field.teleconnection.calc_index_NAH_zonal_lat_Li_Wang_2003(
        slp_data_DJF_mean
    )
    result_data = np.corrcoef(index_ao_point[:-1], index_ao[:-1])[0, 1]
    refer_data = np.array([0.9146496629261469])
    assert np.isclose(np.abs(result_data), refer_data, atol=0.1, equal_nan=True).all()
