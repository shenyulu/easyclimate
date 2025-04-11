"""
pytest for field/typhoon/potential_intensity
"""

import pytest

import easyclimate as ecl
import xarray as xr
import numpy as np

testdata_sst = xr.DataArray(
    np.array([[28.861475, 28.552972], [29.157007, 28.649042]]),
    dims=("lat", "lon"),
    coords={"lat": np.array([2.5, 0.0]), "lon": np.array([-177.5, -175.0])},
    attrs={"standard_name": "Sea Surface Temperature", "units": "degrees C"},
)
testdata_msl = xr.DataArray(
    np.array([[1008.075708, 1008.171802], [1008.005922, 1008.088748]]),
    dims=("lat", "lon"),
    coords={"lat": np.array([2.5, 0.0]), "lon": np.array([-177.5, -175.0])},
    attrs={"standard_name": "Mean Sea Level Pressure", "units": "hPa"},
)
testdata_t = xr.DataArray(
    np.array(
        [
            [[27.353636, 27.222672], [27.21099, 26.938538]],
            [[25.190057, 25.069938], [25.007941, 24.750216]],
            [[23.092178, 22.964785], [22.902211, 22.671951]],
            [[21.464064, 21.194161], [21.422638, 21.28721]],
            [[20.135445, 19.904204], [20.2931, 20.135855]],
            [[19.135903, 19.107665], [19.392363, 19.440976]],
            [[18.297735, 18.37927], [18.58249, 18.641281]],
            [[17.206641, 17.328192], [17.450027, 17.518263]],
            [[15.9993, 16.088557], [16.159143, 16.248171]],
            [[14.678321, 14.700328], [14.817254, 14.830638]],
            [[13.276337, 13.264193], [13.373714, 13.333164]],
            [[11.805724, 11.77128], [11.80218, 11.800547]],
            [[10.245422, 10.189297], [10.17729, 10.205829]],
            [[6.958592, 6.955302], [6.958394, 6.957223]],
            [[3.667473, 3.660736], [3.809784, 3.702164]],
            [[0.086626, 0.099989], [0.178932, 0.153677]],
            [[-4.108487, -3.972281], [-3.976002, -3.855817]],
            [[-8.5458, -8.51175], [-8.630228, -8.564163]],
            [[-13.561634, -13.553428], [-13.84696, -13.798673]],
            [[-20.077641, -19.99213], [-20.159432, -20.108628]],
            [[-28.503128, -28.444674], [-28.420612, -28.401702]],
            [[-38.741623, -38.853656], [-38.611279, -38.734565]],
            [[-51.723332, -51.758043], [-51.649847, -51.732749]],
            [[-67.683181, -67.657047], [-67.849234, -67.727188]],
            [[-84.802469, -84.811532], [-84.599578, -84.656952]],
            [[-82.980426, -83.06524], [-83.076922, -83.306107]],
            [[-72.109836, -71.961065], [-72.277637, -72.140224]],
            [[-65.229228, -65.242002], [-65.227432, -65.291973]],
            [[-56.258264, -56.437082], [-55.930041, -56.037194]],
            [[-51.276108, -51.220998], [-51.082832, -51.050359]],
            [[-44.492267, -44.431645], [-44.665341, -44.653366]],
        ]
    ),
    dims=("level", "lat", "lon"),
    coords={
        "level": np.array(
            [
                1000.0,
                975.0,
                950.0,
                925.0,
                900.0,
                875.0,
                850.0,
                825.0,
                800.0,
                775.0,
                750.0,
                725.0,
                700.0,
                650.0,
                600.0,
                550.0,
                500.0,
                450.0,
                400.0,
                350.0,
                300.0,
                250.0,
                200.0,
                150.0,
                100.0,
                70.0,
                50.0,
                40.0,
                30.0,
                20.0,
                10.0,
            ]
        ),
        "lat": np.array([2.5, 0.0]),
        "lon": np.array([-177.5, -175.0]),
    },
    attrs={"standard_name": "Atmospheric Temperature", "units": "degrees C"},
)
testdata_q = xr.DataArray(
    np.array(
        [
            [[1.752169e01, 1.739373e01], [1.755779e01, 1.747079e01]],
            [[1.711666e01, 1.697109e01], [1.715569e01, 1.708690e01]],
            [[1.677095e01, 1.662145e01], [1.673994e01, 1.663940e01]],
            [[1.593267e01, 1.590978e01], [1.558633e01, 1.544861e01]],
            [[1.473882e01, 1.455964e01], [1.387122e01, 1.373549e01]],
            [[1.304239e01, 1.267559e01], [1.199529e01, 1.146268e01]],
            [[1.127705e01, 1.066488e01], [1.025696e01, 9.864766e00]],
            [[9.847610e00, 9.396217e00], [9.291184e00, 8.872944e00]],
            [[8.918181e00, 8.675946e00], [8.745527e00, 8.273892e00]],
            [[8.272269e00, 8.040189e00], [8.070926e00, 7.926245e00]],
            [[7.787410e00, 7.493204e00], [7.475906e00, 7.490718e00]],
            [[7.253088e00, 6.989394e00], [6.969766e00, 6.952205e00]],
            [[6.702615e00, 6.516205e00], [6.319682e00, 6.375766e00]],
            [[5.709916e00, 5.521406e00], [5.232806e00, 5.305737e00]],
            [[4.398887e00, 4.284769e00], [4.266110e00, 4.302661e00]],
            [[3.194687e00, 2.911649e00], [3.216351e00, 3.126061e00]],
            [[2.166236e00, 1.930091e00], [2.258709e00, 2.104025e00]],
            [[1.585033e00, 1.382272e00], [1.726041e00, 1.579070e00]],
            [[1.132454e00, 9.975798e-01], [1.327680e00, 1.184414e00]],
            [[7.866390e-01, 7.202522e-01], [8.988403e-01, 8.505270e-01]],
            [[4.513145e-01, 4.453679e-01], [5.425988e-01, 5.050453e-01]],
            [[2.165959e-01, 1.997725e-01], [2.575166e-01, 2.287391e-01]],
            [[7.002809e-02, 6.306999e-02], [8.013048e-02, 6.932322e-02]],
            [[1.401548e-02, 1.282958e-02], [1.460552e-02, 1.323063e-02]],
            [[1.314765e-03, 1.300477e-03], [1.292135e-03, 1.284880e-03]],
            [[1.717125e-03, 1.726488e-03], [1.672255e-03, 1.674200e-03]],
            [[2.587957e-03, 2.586587e-03], [2.588366e-03, 2.586200e-03]],
            [[2.572437e-03, 2.572403e-03], [2.572174e-03, 2.572327e-03]],
            [[2.530370e-03, 2.530549e-03], [2.530676e-03, 2.530960e-03]],
            [[2.556868e-03, 2.556884e-03], [2.555793e-03, 2.555478e-03]],
            [[2.825350e-03, 2.825275e-03], [2.810023e-03, 2.809112e-03]],
        ]
    ),
    dims=("level", "lat", "lon"),
    coords={
        "level": np.array(
            [
                1000.0,
                975.0,
                950.0,
                925.0,
                900.0,
                875.0,
                850.0,
                825.0,
                800.0,
                775.0,
                750.0,
                725.0,
                700.0,
                650.0,
                600.0,
                550.0,
                500.0,
                450.0,
                400.0,
                350.0,
                300.0,
                250.0,
                200.0,
                150.0,
                100.0,
                70.0,
                50.0,
                40.0,
                30.0,
                20.0,
                10.0,
            ]
        ),
        "lat": np.array([2.5, 0.0]),
        "lon": np.array([-177.5, -175.0]),
    },
    attrs={"standard_name": "Specific Humidity", "units": "g/kg"},
)


def test_calc_potential_intensity_Bister_Emanuel_2002():
    result_data = ecl.field.typhoon.calc_potential_intensity_Bister_Emanuel_2002(
        sst_data=testdata_sst,
        sst_data_units="degC",
        surface_pressure_data=testdata_msl,
        surface_pressure_data_units="hPa",
        temperature_data=testdata_t,
        temperature_data_units="degC",
        specific_humidity_data=testdata_q,
        specific_humidity_data_units="g/kg",
        vertical_dim="level",
        vertical_dim_units="hPa",
    )
    result_data1 = result_data.vmax.data.flatten()[1:3]
    result_data2 = result_data.pmin.data.flatten()[1:3]
    result_data3 = result_data.ifl.data.flatten()[1:3]
    result_data4 = result_data.t0.data.flatten()[1:3]
    result_data5 = result_data.otl.data.flatten()[1:3]
    result_data6 = result_data.eff.data.flatten()[1:3]
    result_data7 = result_data.diseq.data.flatten()[1:3]
    result_data8 = result_data.lnpi.data.flatten()[1:3]
    result_data9 = result_data.lneff.data.flatten()[1:3]
    result_data10 = result_data.lndiseq.data.flatten()[1:3]
    result_data11 = result_data.lnCKCD.data.flatten()[1:3]

    refer_data1 = np.array([90.93510812, 98.74747325])
    refer_data2 = np.array([900.4934007, 882.44831663])
    refer_data3 = np.array([1, 1])
    refer_data4 = np.array([189.45530961, 189.75865991])
    refer_data5 = np.array([80.81349033, 76.1947957])
    refer_data6 = np.array([0.59247567, 0.5931131])
    refer_data7 = np.array([15507.79821296, 18267.19889297])
    refer_data8 = np.array([9.02029231, 9.18513163])
    refer_data9 = np.array([-0.52344546, -0.52237017])
    refer_data10 = np.array([9.64909829, 9.81286232])
    refer_data11 = np.array([-0.10536052, -0.10536052])

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()
    assert np.isclose(result_data6, refer_data6).all()
    assert np.isclose(result_data7, refer_data7).all()
    assert np.isclose(result_data8, refer_data8).all()
    assert np.isclose(result_data9, refer_data9).all()
    assert np.isclose(result_data10, refer_data10).all()
    assert np.isclose(result_data11, refer_data11).all()
