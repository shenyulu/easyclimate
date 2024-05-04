"""
pytest for field.air_sea_interaction.index_amm.py
"""

import pytest

import easyclimate as ecl
import xarray as xr
import numpy as np
from pathlib import Path
from .const_define import TEST_DATA_PATH

data_sst = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_oisst_data.nc")))["sst"]


def test_calc_index_AMM_Doi_2009_1():
    result_data = ecl.field.air_sea_interaction.calc_index_AMM_Doi_2009(data_sst).data[
        :20
    ]
    refer_data = np.array(
        [
            -1.5725382e-04,
            1.6361520e-01,
            1.8975957e-01,
            3.5046536e-01,
            3.8568866e-01,
            7.5838298e-01,
            4.3417996e-01,
            2.0991376e-01,
            4.6612772e-01,
            1.5847254e-01,
            -2.6524976e-01,
            -2.1791753e-01,
            -4.4400340e-01,
            -6.3130760e-01,
            -5.8240391e-02,
            9.3769622e-01,
            1.3219962e00,
            1.2628446e00,
            1.1743338e00,
            6.2622142e-01,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_index_AMM_Doi_2009_2():
    result_data = ecl.field.air_sea_interaction.calc_index_AMM_Doi_2009(
        data_sst, normalized=True
    ).data[:20]
    refer_data = np.array(
        [
            -3.2382132e-04,
            3.3692086e-01,
            3.9075804e-01,
            7.2168773e-01,
            7.9422057e-01,
            1.5616828e00,
            8.9407516e-01,
            4.3226010e-01,
            9.5986283e-01,
            3.2633093e-01,
            -5.4620951e-01,
            -4.4874167e-01,
            -9.1430384e-01,
            -1.3000057e00,
            -1.1993019e-01,
            1.9309294e00,
            2.7222905e00,
            2.6004837e00,
            2.4182200e00,
            1.2895321e00,
        ],
        dtype=np.float32,
    )
    assert np.isclose(result_data, refer_data).all()
