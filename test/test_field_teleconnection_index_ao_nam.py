"""
pytest for field/teleconnection/index_ao_nam.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import DOCS_DATA_PATH

slp_data = xr.open_dataset(str(Path(DOCS_DATA_PATH, "slp_monmean_NH.nc")))["slp"]


def test_calc_index_AO_EOF_Thompson_Wallace_1998_AND_calc_index_NAH_zonal_lat_Li_Wang_2003():
    slp_data_DJF_mean = ecl.calc_seasonal_mean(slp_data, extract_season="DJF")
    index_ao = ecl.field.teleconnection.calc_index_AO_EOF_Thompson_Wallace_1998(
        slp_data_DJF_mean
    )
    index_ao_point = ecl.field.teleconnection.calc_index_NAH_zonal_lat_Li_Wang_2003(
        slp_data_DJF_mean
    )
    result = np.corrcoef(index_ao_point, index_ao)[0, 1]
    assert np.abs(result) > 0.9
