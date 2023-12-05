"""
pytest for diff.py
"""
import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from .util import round_sf_np

# Sample Data Declaration
t_data = xr.DataArray(
    np.array([[268.03, 268.13, 268.19],
       [268.28, 268.33, 268.28],
       [268.46, 268.42, 268.25]]),
    dims = ('lat', 'lon'),
    coords = {'lat': np.array([5.0, 2.5, 0.0]), 'lon': np.array([0.0, 2.5, 5.0])}
)

t_level_data = xr.DataArray(
    np.array([299.80002, 295.4, 291.2, 282.7, 275.80002, 268., 258.1, 243.40001, 232.7, 220.3, 205.1, 192.40001, 196. , 203. ,212.3, 220.5, 232.1]),
    dims = 'level',
    coords = {'level': np.array([1000., 925., 850., 700., 600., 500., 400., 300., 250., 200., 150., 100., 70., 50., 30., 20., 10.])}
)

def test_calc_gradient():
    result_data = ecl.calc_gradient(t_data, dim = 'lon').data
    refer_data = np.array([[ 0.12 ,  0.08 ,  0.04 ], [ 0.1  ,  0.   , -0.1  ], [ 0.025, -0.105, -0.235]])
    assert np.isclose(result_data.data, refer_data).all()

def test_calc_lon_gradient():
    result_data = ecl.calc_lon_gradient(t_data).data
    refer_data = np.array([[ 4.33391322e-07,  2.88927548e-07,  1.44463774e-07],
       [ 3.60127877e-07,  0.00000000e+00, -3.60127877e-07],
       [ 8.99462787e-08, -3.77774370e-07, -8.45495020e-07]])
    assert np.isclose(result_data, refer_data).all()

def test_calc_lat_gradient():
    result_data = ecl.calc_lat_gradient(t_data).data
    refer_data = np.array([[-1.02538758e-06, -9.17452042e-07, -5.39677672e-07],
       [-7.73537997e-07, -5.21688416e-07, -1.07935534e-07],
       [-5.21688416e-07, -1.25924790e-07,  3.23806603e-07]])
    assert np.isclose(result_data, refer_data).all()

def test_calc_lon_laplacian():
    result_data = ecl.calc_lon_laplacian(t_data).data
    refer_data = np.array([[-5.21744551e-13, -5.21744551e-13, -5.21744551e-13],
       [-1.29692088e-12, -1.29692088e-12, -1.29692088e-12],
       [-1.68278927e-12, -1.68278927e-12, -1.68278927e-12]])
    assert np.isclose(result_data, refer_data).all()

def test_calc_lat_laplacian():
    result_data = ecl.calc_lat_laplacian(t_data).data
    refer_data = np.array([[-9.06117301e-13, -1.42389862e-12, -1.55334394e-12],
       [-9.06117301e-13, -1.42389862e-12, -1.55334394e-12],
       [-9.06117301e-13, -1.42389862e-12, -1.55334394e-12]])
    assert np.isclose(result_data, refer_data).all()

def test_calc_lon_lat_mixed_derivatives():
    result_data = ecl.calc_lon_lat_mixed_derivatives(t_data).data
    refer_data = np.array([[-9.74548417e-14,  8.77093575e-13,  1.85164199e-12],
       [ 6.15451085e-13,  1.19851001e-12,  1.78156893e-12],
       [ 1.32681462e-12,  1.52098261e-12,  1.71515061e-12]])
    assert np.isclose(result_data, refer_data).all()

def test_calc_p_gradient():
    result_data = ecl.calc_p_gradient(t_level_data, vertical_dim = 'level', vertical_dim_units = 'hPa').data
    refer_data = np.array([ 0.0006,  0.00057333,  0.00056444,  0.000616,  0.000735,
        0.000885,  0.00123,  0.00169333,  0.00231,  0.00276,
        0.00279,  0.0011375, -0.00212, -0.004075, -0.00583333,
       -0.0099, -0.0133])
    assert np.isclose(result_data, refer_data).all()