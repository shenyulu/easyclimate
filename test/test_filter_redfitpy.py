"""
pytest for filter/redfitpy.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import TEST_DATA_PATH

sst_data = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_nino34.nc"))).sst


def test_calc_redfit1():
    result_data = ecl.filter.calc_redfit(sst_data)

    result_data1 = result_data.gxx[20:30].data
    result_data3 = result_data.gred_th[20:30].data
    result_data6 = result_data.chi2_80[20:30].data
    result_data7 = result_data.chi2_90[20:30].data
    result_data8 = result_data.chi2_95[20:30].data
    result_data9 = result_data.chi2_99[20:30].data

    refer_data1 = np.array(
        [
            2.6100595,
            1.3072749,
            7.6913514,
            1.6619556,
            5.278525,
            3.2057655,
            0.66721755,
            17.080835,
            15.754893,
            4.182542,
        ]
    )
    refer_data3 = np.array(
        [
            4.374933,
            3.9981651,
            3.6673093,
            3.3753517,
            3.1165347,
            2.8861244,
            2.6801665,
            2.4953728,
            2.3289793,
            2.1786544,
        ]
    )
    refer_data6 = np.array(
        [
            7.0427766,
            6.436255,
            5.9036427,
            5.4336486,
            5.0170045,
            4.6460896,
            4.314538,
            4.0170565,
            3.749196,
            3.5072029,
        ]
    )
    refer_data7 = np.array(
        [
            10.071838,
            9.204455,
            8.442769,
            7.7706327,
            7.1747923,
            6.644348,
            6.1701984,
            5.744772,
            5.361706,
            5.0156326,
        ]
    )
    refer_data8 = np.array(
        [
            13.109244,
            11.980281,
            10.98889,
            10.114055,
            9.338524,
            8.648112,
            8.030971,
            7.4772468,
            6.9786577,
            6.5282173,
        ]
    )
    refer_data9 = np.array(
        [
            20.143677,
            18.40891,
            16.885538,
            15.5412655,
            14.349585,
            13.288696,
            12.340397,
            11.489544,
            10.723412,
            10.031265,
        ]
    )

    assert np.isclose(result_data1, refer_data1, atol=0.01).all()
    assert np.isclose(result_data3, refer_data3, atol=0.01).all()
    assert np.isclose(result_data6, refer_data6, atol=0.01).all()
    assert np.isclose(result_data7, refer_data7, atol=0.01).all()
    assert np.isclose(result_data8, refer_data8, atol=0.01).all()
    assert np.isclose(result_data9, refer_data9, atol=0.01).all()


def test_calc_redfit2():
    result_data = ecl.filter.calc_redfit(sst_data, iwin="welch")
    result_data1 = result_data.gxx[20:30].data
    refer_data1 = np.array(
        [
            1.1525935,
            1.5818456,
            6.1552815,
            5.997016,
            4.5128694,
            3.7235422,
            0.19082911,
            21.362434,
            19.166662,
            2.2907004,
        ]
    )
    assert np.isclose(result_data1, refer_data1, atol=0.01).all()


def test_calc_redfit3():
    result_data = ecl.filter.calc_redfit(sst_data, iwin="hanning")
    result_data1 = result_data.gxx[20:30].data
    refer_data1 = np.array(
        [
            0.48377505,
            1.1684577,
            5.568493,
            7.5657167,
            6.037566,
            3.4565885,
            1.446536,
            21.131693,
            16.22255,
            0.6097581,
        ]
    )
    assert np.isclose(result_data1, refer_data1, atol=0.01).all()


def test_calc_redfit4():
    result_data = ecl.filter.calc_redfit(sst_data, iwin="triangular")
    result_data1 = result_data.gxx[20:30].data
    refer_data1 = np.array(
        [
            0.53083134,
            1.0754433,
            5.9494467,
            6.823731,
            5.0147853,
            3.5979056,
            0.74510956,
            19.944061,
            16.900396,
            1.1186556,
        ]
    )
    assert np.isclose(result_data1, refer_data1, atol=0.01).all()


def test_calc_redfit5():
    result_data = ecl.filter.calc_redfit(sst_data, iwin="blackmanharris")
    result_data1 = result_data.gxx[20:30].data
    refer_data1 = np.array(
        [
            0.25909367,
            0.9457828,
            5.4999194,
            8.3161955,
            6.943804,
            3.0793598,
            2.3809066,
            19.174746,
            13.924949,
            0.18201135,
        ]
    )
    assert np.isclose(result_data1, refer_data1, atol=0.01).all()


def test_calc_redfit6():
    result_data = ecl.filter.calc_redfit(sst_data, iwin="rectangular", mctest=True)
    result_data1 = result_data.gxx[20:30].data
    refer_data1 = np.array(
        [
            2.6100595,
            1.3072749,
            7.6913514,
            1.6619556,
            5.278525,
            3.2057655,
            0.66721755,
            17.080835,
            15.754893,
            4.182542,
        ]
    )
    assert np.isclose(result_data1, refer_data1, atol=0.01).all()


def test_calc_redfit_cross1():
    gxx, gyy, gxy, cxy, phxy = ecl.filter.calc_redfit_cross(
        sst_data.data[:20], sst_data.data[:20]
    )
    result_gxx_data1 = gxx.gxx[5:].data
    refer_gxx_data1 = np.array(
        [0.49393213, 0.04329514, 0.10684587, 0.03644858, 0.0762934, 0.00138481]
    )

    assert np.isclose(result_gxx_data1, refer_gxx_data1, atol=0.01).all()


def test_calc_redfit_cross2():
    gxx, gyy, gxy, cxy, phxy = ecl.filter.calc_redfit_cross(
        sst_data.data[:20], sst_data.data[:20], iwin="welch"
    )
    result_gxx_data1 = gxx.gxx[5:].data
    refer_gxx_data1 = np.array(
        [0.15475138, 0.06211971, 0.02606629, 0.02504181, 0.0246939, 0.02436013]
    )

    assert np.isclose(result_gxx_data1, refer_gxx_data1, atol=0.01).all()


def test_calc_redfit_cross3():
    gxx, gyy, gxy, cxy, phxy = ecl.filter.calc_redfit_cross(
        sst_data.data[:20], sst_data.data[:20], iwin="hanning"
    )
    result_gxx_data1 = gxx.gxx[5:].data
    refer_gxx_data1 = np.array(
        [0.10175404, 0.06451696, 0.01783044, 0.00387535, 0.04815307, 0.00995838]
    )

    assert np.isclose(result_gxx_data1, refer_gxx_data1, atol=0.01).all()


def test_calc_redfit_cross4():
    gxx, gyy, gxy, cxy, phxy = ecl.filter.calc_redfit_cross(
        sst_data.data[:20], sst_data.data[:20], iwin="triangular"
    )
    result_gxx_data1 = gxx.gxx[5:].data
    refer_gxx_data1 = np.array(
        [0.08821693, 0.0617842, 0.01330115, 0.02466889, 0.02898028, 0.02473436]
    )

    assert np.isclose(result_gxx_data1, refer_gxx_data1, atol=0.01).all()


def test_calc_redfit_cross5():
    gxx, gyy, gxy, cxy, phxy = ecl.filter.calc_redfit_cross(
        sst_data.data[:20], sst_data.data[:20], iwin="blackmanharris"
    )
    result_gxx_data1 = gxx.gxx[5:].data
    refer_gxx_data1 = np.array(
        [0.06423467, 0.06125515, 0.02167186, 0.00830468, 0.03833103, 0.01371629]
    )

    assert np.isclose(result_gxx_data1, refer_gxx_data1, atol=0.01).all()


def test_calc_redfit_cross6():
    gxx, gyy, gxy, cxy, phxy = ecl.filter.calc_redfit_cross(
        sst_data.data[:20], sst_data.data[:20], mctest=False, mctest_phi=False
    )
    result_gxx_data1 = gxx.gxx[5:].data
    refer_gxx_data1 = np.array(
        [0.49393213, 0.04329514, 0.10684587, 0.03644858, 0.0762934, 0.00138481]
    )

    assert np.isclose(result_gxx_data1, refer_gxx_data1, atol=0.01).all()
