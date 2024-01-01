"""
pytest for field/monsoon/index_npwi.py
"""
import pytest

import easyclimate as ecl
import numpy as np

lon_start, lon_end, lat_start, lat_end = 120, 125, 10, 20
pw_data = ecl.open_tutorial_dataset('pr_wtr_eatm_2022').pr_wtr
NPWI_index = ecl.field.monsoon.calc_index_NPWI(pw_data)
monsoon_onset_date = ecl.field.monsoon.calc_NPWI_monsoon_onset(NPWI_index)

def test_find_PW_monsoon_region():
    result_data = ecl.field.monsoon.find_PW_monsoon_region(pw_data).sel(lon = slice(lon_start, lon_end), lat = slice(lat_end, lat_start)).data.flatten()
    refer_data = np.array([ True, False, False, False, False, False, False, False, False,
       False, False, False, False, False, False])
    assert np.isclose(result_data, refer_data).all()

def test_calc_index_NPWI():
    result_data = NPWI_index.isel(time = 2).sel(lon = slice(lon_start, lon_end), lat = slice(lat_end, lat_start)).data.flatten()
    refer_data = np.array([0.30897176, 0.34349486, 0.312401  , 0.22050983, 0.14875537,
       0.11090223, 0.21444692, 0.14592266, 0.17253728, 0.4068323 ,
       0.35567716, 0.31826746, 0.65482634, 0.45058355, 0.30517498])
    assert np.isclose(result_data, refer_data).all()

def test_calc_NPWI_monsoon_onset():
    result_data = ecl.field.monsoon.calc_NPWI_monsoon_onset(NPWI_index).sel(lon = slice(lon_start, lon_end), lat = slice(lat_end, lat_start)).data.flatten()
    refer_data = np.array([137., 136., 134., 137., 141., 142., 137., 142., 179., 141., 101.,68.,  69.,  66.,  25.])
    assert np.isclose(result_data, refer_data).all()

def test_calc_NPWI_monsoon_detreat():
    result_data = ecl.field.monsoon.calc_NPWI_monsoon_detreat(NPWI_index, monsoon_onset_date).sel(lon = slice(lon_start, lon_end), lat = slice(lat_end, lat_start)).data.flatten()
    refer_data = np.array([192., 173., 164., 157., 158., 151., 145., 145., 202., 145., 106.,77.,  83.,  77.,  28.])
    assert np.isclose(result_data, refer_data).all()