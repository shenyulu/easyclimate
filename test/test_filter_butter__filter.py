"""
pytest for yearstat.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import TEST_DATA_PATH

data_nino34_area = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_nino34.nc"))
).sst


def test_calc_butter_bandpass():
    result_data = ecl.filter.calc_butter_bandpass(
        data_nino34_area, sampling_frequency=1 / 12, period=[3, 10]
    )
    result_data = ecl.calc_yearly_climatological_mean(result_data).data
    refer_data = np.array(
        [
            0.88364814,
            1.12401429,
            -0.58469931,
            -0.7815492,
            0.67832723,
            1.02367853,
            -0.58776442,
            -1.10863814,
            -0.01244383,
            0.58176857,
            0.21445867,
            0.00334761,
            0.09636685,
            -0.4594995,
            -0.28476831,
            1.0236293,
            0.50226878,
            -0.96091926,
            -0.81687546,
            0.07020657,
            0.51648347,
            0.33800556,
            0.05252858,
            0.06410267,
            0.06374163,
            -0.44170631,
            -0.30253664,
            0.60143959,
            0.16784339,
            -0.48528246,
            -0.15203788,
            -0.19446214,
            0.07205399,
            0.89042826,
            0.22708011,
            -0.82872356,
            -0.26053022,
            0.41860927,
            0.02154499,
            -0.15823524,
            0.04580367,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_butter_lowpass():
    result_data = ecl.filter.calc_butter_lowpass(
        data_nino34_area, sampling_frequency=1 / 12, period=10
    )
    result_data = ecl.calc_yearly_climatological_mean(result_data).data
    refer_data = np.array(
        [
            26.26998031,
            26.70600132,
            26.95900997,
            27.06218374,
            27.07031735,
            27.03318847,
            27.00213554,
            27.02631035,
            27.11671195,
            27.23591063,
            27.32910936,
            27.35890634,
            27.31697658,
            27.21615982,
            27.07999974,
            26.93201064,
            26.80018071,
            26.73114669,
            26.76403943,
            26.8900334,
            27.0519318,
            27.17777965,
            27.21860984,
            27.16367514,
            27.03659814,
            26.88178891,
            26.74788864,
            26.67033885,
            26.66806135,
            26.75397889,
            26.92392561,
            27.14188418,
            27.34836116,
            27.48001626,
            27.49866775,
            27.41349219,
            27.2642859,
            27.09224082,
            26.93409448,
            26.82198982,
            26.76820941,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_calc_butter_highpass():
    result_data = ecl.filter.calc_butter_highpass(
        data_nino34_area, sampling_frequency=1 / 12, period=1.5
    )
    result_data = ecl.calc_yearly_climatological_mean(result_data).data
    refer_data = np.array(
        [
            -0.06685445,
            -0.05564417,
            0.09081622,
            0.05948681,
            -0.03326621,
            -0.00228173,
            -0.08090729,
            0.04276665,
            0.03057548,
            -0.07353261,
            -0.02076468,
            0.11494913,
            -0.01651844,
            -0.07578275,
            0.12437629,
            -0.03503192,
            -0.12090074,
            0.10903428,
            -0.02588454,
            0.01310894,
            0.06333923,
            -0.04060118,
            -0.01043272,
            0.02507494,
            0.0071918,
            -0.06310646,
            0.09922859,
            -0.09663018,
            -0.07625251,
            0.09578758,
            0.0459007,
            -0.02939997,
            0.04412075,
            -0.0421145,
            -0.1399516,
            0.1290685,
            0.02595735,
            -0.09239391,
            -0.01291601,
            0.08025345,
            0.01768816,
        ]
    )
    assert np.isclose(result_data, refer_data).all()
