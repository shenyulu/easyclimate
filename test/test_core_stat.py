"""
pytest for stat.py
"""
import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
import pandas as pd

sst_data = ecl.tutorial.open_tutorial_dataset('mini_HadISST_sst').sst
sic_data_Barents_Sea = ecl.tutorial.open_tutorial_dataset('mini_HadISST_ice').sic
sic_data_Barents_Sea_12 = ecl.get_specific_months_data(sic_data_Barents_Sea, 12)

def test_calc_linregress_spatial():
    result_data = ecl.calc_linregress_spatial(sic_data_Barents_Sea_12.sel(lon = slice(34.5, 36.5), lat = slice(78.5, 80.5)), dim = 'time').compute()
    result_data1 = result_data['slope'].data
    result_data2 = result_data['intercept'].data
    result_data3 = result_data['rvalue'].data
    result_data4 = result_data['pvalue'].data
    result_data5 = result_data['stderr'].data
    result_data6 = result_data['intercept_stderr'].data

    refer_data1 = np.array([[-0.01689814, -0.01618345, -0.01640629],
       [-0.01011993, -0.00922373, -0.0091192 ],
       [-0.00641115, -0.0054169 , -0.00600519]])
    refer_data2 = np.array([[1.05593576, 1.04461794, 1.04132889],
       [1.01817285, 1.01218153, 1.01741975],
       [0.89857147, 0.91509416, 0.94763013]])
    refer_data3 = np.array([[-0.58339207, -0.57217978, -0.57376992],
       [-0.47457306, -0.4485609 , -0.45343254],
       [-0.32601495, -0.28470031, -0.33127693]])
    refer_data4 = np.array([[5.01586941e-05, 7.52887062e-05, 7.11385575e-05],
       [1.49647207e-03, 2.88846392e-03, 2.56361219e-03],
       [3.51172733e-02, 6.76380713e-02, 3.21088821e-02]])
    refer_data5 = np.array([[0.00371969, 0.00366767, 0.00370284],
       [0.00296779, 0.00290584, 0.00283422],
       [0.00293946, 0.00288389, 0.00270435]])
    refer_data6 = np.array([[0.08858534, 0.08734656, 0.08818416],
       [0.07067876, 0.06920341, 0.06749764],
       [0.07000405, 0.06868049, 0.06440477]])

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()
    assert np.isclose(result_data6, refer_data6).all()

def test_calc_linregress_spatial():
    result_data = ecl.calc_linregress_spatial(sic_data_Barents_Sea_12.sel(lon = slice(34.5, 36.5), lat = slice(78.5, 80.5)), dim = 'time').compute()
    result_data1 = result_data['slope'].data
    result_data2 = result_data['intercept'].data
    result_data3 = result_data['rvalue'].data
    result_data4 = result_data['pvalue'].data
    result_data5 = result_data['stderr'].data
    result_data6 = result_data['intercept_stderr'].data

    refer_data1 = np.array([[-0.01689814, -0.01618345, -0.01640629],
       [-0.01011993, -0.00922373, -0.0091192 ],
       [-0.00641115, -0.0054169 , -0.00600519]])
    refer_data2 = np.array([[1.05593576, 1.04461794, 1.04132889],
       [1.01817285, 1.01218153, 1.01741975],
       [0.89857147, 0.91509416, 0.94763013]])
    refer_data3 = np.array([[-0.58339207, -0.57217978, -0.57376992],
       [-0.47457306, -0.4485609 , -0.45343254],
       [-0.32601495, -0.28470031, -0.33127693]])
    refer_data4 = np.array([[5.01586941e-05, 7.52887062e-05, 7.11385575e-05],
       [1.49647207e-03, 2.88846392e-03, 2.56361219e-03],
       [3.51172733e-02, 6.76380713e-02, 3.21088821e-02]])
    refer_data5 = np.array([[0.00371969, 0.00366767, 0.00370284],
       [0.00296779, 0.00290584, 0.00283422],
       [0.00293946, 0.00288389, 0.00270435]])
    refer_data6 = np.array([[0.08858534, 0.08734656, 0.08818416],
       [0.07067876, 0.06920341, 0.06749764],
       [0.07000405, 0.06868049, 0.06440477]])

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()
    assert np.isclose(result_data6, refer_data6).all()

def test_calc_linregress_spatial2():
    result_data = ecl.calc_linregress_spatial(sic_data_Barents_Sea_12.sel(lon = slice(34.5, 36.5), lat = slice(78.5, 80.5)), dim = 'time', engine = 'xarray').compute()
    result_data1 = result_data['slope'].data
    result_data2 = result_data['intercept'].data
    result_data3 = result_data['rvalue'].data
    result_data4 = result_data['pvalue'].data
    result_data5 = result_data['stderr'].data
    result_data6 = result_data['intercept_stderr'].data

    refer_data1 = np.array([[-0.01689814, -0.01618345, -0.01640629],
       [-0.01011993, -0.00922373, -0.0091192 ],
       [-0.00641115, -0.0054169 , -0.00600519]])
    refer_data2 = np.array([[1.0559357 , 1.04461782, 1.04132889],
       [1.01817267, 1.01218159, 1.01741975],
       [0.89857135, 0.91509404, 0.94763007]])
    refer_data3 = np.array([[-0.5833921 , -0.57217977, -0.57376993],
       [-0.47457306, -0.4485609 , -0.45343254],
       [-0.32601492, -0.2847003 , -0.33127695]])
    refer_data4 = np.array([[0.54670677, 0.55950058, 0.55297443],
       [0.70697605, 0.73026712, 0.73443278],
       [0.78751658, 0.82292746, 0.81043421]])
    refer_data5 = np.array([[0.02779874, 0.02749916, 0.02741885],
       [0.02672881, 0.02656663, 0.02669478],
       [0.02362677, 0.02404782, 0.02487059]])
    refer_data6 = np.array([[0.66203424, 0.6548997 , 0.65298707],
       [0.63655362, 0.63269127, 0.63574306],
       [0.5626778 , 0.57270522, 0.59229961]])

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()
    assert np.isclose(result_data6, refer_data6).all()

def test_calc_detrend_data():
    result_data = ecl.calc_detrend_data(sic_data_Barents_Sea_12.sel(lon = slice(34.5, 36.5), lat = slice(78.5, 80.5)), time_dim = 'time').mean(dim = ('lat', 'lon')).data
    refer_data = np.array([-5.89946248e-02, -4.94630672e-02, -4.54870537e-02, -4.94844377e-01,
       -2.19795033e-02, -4.67017619e-03,  2.04169489e-02,  3.32818367e-02,
        1.28133763e-02,  2.90115941e-02, -3.45740060e-04,  2.47413721e-02,
        5.31618260e-02,  3.04711461e-02,  5.55582643e-02,  1.49534270e-01,
        5.79547025e-02,  1.48597360e-01,  1.53684452e-01,  1.53216019e-01,
        2.16364730e-02,  1.21168017e-01,  1.79588467e-01,  1.50231123e-01,
        1.15318246e-01, -3.25150162e-01, -9.78408679e-02,  1.57246247e-01,
       -1.24333329e-01,  2.94087142e-01,  1.32507533e-01, -4.22405332e-01,
        1.76015079e-01,  2.68880010e-01, -2.03810692e-01, -6.22056901e-01,
       -2.34747574e-01, -6.00771606e-01,  3.94315481e-01,  5.82915097e-02,
        3.93378615e-01, -7.82009140e-02])
    assert np.isclose(result_data, refer_data).all()

def test_calc_ttestSpatialPattern_spatial():
    sic_data_Barents_Sea_12_spatial_mean = sic_data_Barents_Sea_12.sel(lon = slice(34.5, 36.5), lat = slice(78.5, 80.5)).mean(dim = ('lat', 'lon'))
    test1 = sic_data_Barents_Sea_12_spatial_mean.isel(time = slice(0, 20))
    test2 = sic_data_Barents_Sea_12_spatial_mean.isel(time = slice(21, None))
    result_data = ecl.calc_ttestSpatialPattern_spatial(test1, test2)
    result_data1 = result_data['statistic'].data
    result_data2 = result_data['pvalue'].data
    refer_data1 = 3.43382768
    refer_data2 = 0.00142438
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()

def test_calc_skewness_spatial():
    sic_data_Barents_Sea_12_detrend = ecl.calc_detrend_data(sic_data_Barents_Sea_12.sel(lon = slice(34.5, 36.5), lat = slice(78.5, 80.5)), time_dim = 'time')
    result_data = ecl.calc_skewness_spatial(sic_data_Barents_Sea_12_detrend, dim = 'time')
    result_data1 = result_data['skewness'].data
    result_data2 = result_data['pvalue'].data
    refer_data1 = np.array([[-0.24333526, -0.32189173, -0.27430525],
       [-1.3397094 , -1.5588326 , -1.6165946 ],
       [-1.8677251 , -2.209491  , -2.330299  ]])
    refer_data2 = np.array([[7.70400089e-01, 6.26686378e-01, 7.18524740e-01],
       [4.16622922e-04, 3.66448314e-05, 1.56112432e-05],
       [1.88511919e-06, 6.86471937e-08, 1.30767304e-08]])
    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    
def test_calc_kurtosis_spatial():
    sic_data_Barents_Sea_12_detrend = ecl.calc_detrend_data(sic_data_Barents_Sea_12.sel(lon = slice(34.5, 36.5), lat = slice(78.5, 80.5)), time_dim = 'time')
    result_data = ecl.calc_kurtosis_spatial(sic_data_Barents_Sea_12_detrend, dim = 'time')
    refer_data = np.array([[2.7300231, 2.8368442, 2.7555804],
       [4.667509 , 5.464381 , 5.85138  ],
       [6.32586  , 7.463421 , 8.428185 ]])
    assert np.isclose(result_data, refer_data).all()
