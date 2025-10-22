"""
pytest for eof.py
"""

import pytest

import easyclimate as ecl
import xarray as xr
import numpy as np
import xeofs
from pathlib import Path
from .const_define import TEST_DATA_PATH
from .const_define import TEST_TMP_PATH
from .util import assert_path_dir_exist
from easyclimate.core.eof import calc_eof_projection_coefficient


def test_get_EOF_model_and_calc_EOF_analysis():
    inputdata = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_core_eof.nc"))).z
    model = ecl.eof.get_EOF_model(
        inputdata,
        lat_dim="lat",
        lon_dim="lon",
        remove_seasonal_cycle_mean=True,
        use_coslat=True,
    )
    result_data = ecl.eof.calc_EOF_analysis(model)
    result_data1 = result_data.EOF.isel(mode=0).data.flatten()
    result_data2 = result_data.PC.isel(mode=0, time=slice(None, 20)).data.flatten()
    result_data3 = result_data.explained_variance.data.flatten()
    result_data4 = result_data.explained_variance_ratio.data.flatten()
    result_data5 = result_data.singular_values.data.flatten()

    refer_data1 = np.array(
        [
            0.21979125,
            0.22049542,
            0.22112887,
            0.22203068,
            0.22288253,
            0.22176682,
            0.22243614,
            0.22330636,
            0.2243775,
            0.22592835,
            0.22290957,
            0.22398918,
            0.22526811,
            0.22591317,
            0.22645566,
            0.22210256,
            0.2235009,
            0.2250623,
            0.22565883,
            0.22695101,
        ]
    )
    refer_data2 = np.array(
        [
            -0.04923435,
            -0.053984,
            -0.0365616,
            -0.0285534,
            0.00723713,
            0.00680275,
            -0.03870684,
            -0.02295436,
            -0.01665111,
            -0.02548954,
            -0.02032592,
            0.02496027,
            0.03082447,
            0.05606671,
            0.04795404,
            0.05304259,
            0.03197721,
            0.02105718,
            0.03165233,
            0.02392981,
        ]
    )
    refer_data3 = np.array(
        [
            1.64010367e05,
            1.63808525e03,
            2.99641421e02,
            1.87287196e01,
            1.67479146e01,
            1.03994458e01,
            6.61782085e00,
            6.29723182e00,
            6.11247419e00,
            5.65868326e00,
        ]
    )
    refer_data4 = np.array(
        [
            9.87641956e-01,
            9.86426494e-03,
            1.80438861e-03,
            1.12781098e-04,
            1.00853034e-04,
            6.26236568e-05,
            3.98513682e-05,
            3.79208367e-05,
            3.68082583e-05,
            3.40756081e-05,
        ]
    )
    refer_data5 = np.array(
        [
            8973.80021646,
            896.82766221,
            383.56738361,
            95.89474092,
            90.68200514,
            71.45717509,
            57.00307044,
            55.605223,
            54.78343572,
            52.71065815,
        ]
    )

    assert np.isclose(result_data1, refer_data1, atol=0.01).all()
    assert np.isclose(result_data2, refer_data2, atol=0.01).all()
    assert np.isclose(result_data3, refer_data3, atol=0.01).all()
    assert np.isclose(result_data4, refer_data4, atol=0.01).all()
    assert np.isclose(result_data5, refer_data5, atol=0.01).all()


# def test_save_EOF_model_and_load_EOF_model():
#     assert_path_dir_exist(TEST_TMP_PATH)
#     inputdata = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_core_eof.nc"))).z
#     model = ecl.eof.get_EOF_model(
#         inputdata,
#         lat_dim="lat",
#         lon_dim="lon",
#         remove_seasonal_cycle_mean=True,
#         use_coslat=True,
#     )
#     result_data = ecl.eof.calc_EOF_analysis(model)
#     outputfile_path = str(Path(TEST_TMP_PATH, "test_output_core_eof_model.zarr"))
#     ecl.eof.save_EOF_model(
#         model, path=outputfile_path, overwrite=True, save_data=False, engine="zarr"
#     )
#     mymodel = ecl.eof.load_EOF_model(path=outputfile_path, engine="zarr")
#     assert isinstance(mymodel, xeofs.single.eof.EOF)


def test_get_REOF_model_and_calc_REOF_analysis():
    inputdata = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_core_eof.nc"))).z
    rof_model = ecl.eof.get_REOF_model(
        inputdata,
        lat_dim="lat",
        lon_dim="lon",
        remove_seasonal_cycle_mean=True,
        use_coslat=True,
    )
    result_data = ecl.eof.calc_REOF_analysis(rof_model)
    result_data1 = result_data.EOF.isel(mode=0).data.flatten()
    result_data2 = result_data.PC.isel(mode=0, time=slice(None, 20)).data.flatten()
    result_data3 = result_data.explained_variance.data.flatten()
    result_data4 = result_data.explained_variance_ratio.data.flatten()
    result_data5 = result_data.singular_values.data.flatten()

    refer_data1 = np.array(
        [
            0.24682069,
            0.24770336,
            0.24866589,
            0.24919784,
            0.25025373,
            0.23188296,
            0.23375969,
            0.23396476,
            0.23576607,
            0.23693779,
            0.21245328,
            0.2131522,
            0.21435789,
            0.21489508,
            0.21542574,
            0.19088355,
            0.19268268,
            0.19392525,
            0.19414491,
            0.19573786,
        ]
    )
    refer_data2 = np.array(
        [
            -0.01686419,
            -0.02530423,
            -0.04233954,
            0.03120813,
            -0.01915748,
            0.00317191,
            -0.0159518,
            -0.00773108,
            0.03129432,
            -0.02648422,
            0.03334937,
            0.03512614,
            0.12212657,
            0.13570949,
            0.14529659,
            0.05060206,
            0.02823904,
            0.01518247,
            -0.01690386,
            0.01378143,
        ]
    )
    refer_data3 = np.array([83455.58916613, 82192.86377077])
    refer_data4 = np.array([0.50255507, 0.49495115])
    refer_data5 = np.array([6401.3041078, 6352.69203657])

    assert np.isclose(result_data1, refer_data1, atol=0.01).all()
    assert np.isclose(result_data2, refer_data2, atol=0.01).all()
    assert np.isclose(result_data3, refer_data3, atol=0.01).all()
    assert np.isclose(result_data4, refer_data4, atol=0.01).all()
    assert np.isclose(result_data5, refer_data5, atol=0.01).all()


# def test_save_REOF_model_and_load_REOF_model():
#     assert_path_dir_exist(TEST_TMP_PATH)
#     inputdata = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_core_eof.nc"))).z
#     model = ecl.eof.get_REOF_model(
#         inputdata,
#         lat_dim="lat",
#         lon_dim="lon",
#         remove_seasonal_cycle_mean=True,
#         use_coslat=True,
#     )
#     result_data = ecl.eof.calc_REOF_analysis(model)
#     outputfile_path = str(Path(TEST_TMP_PATH, "test_output_core_reof_model.zarr"))
#     ecl.eof.save_REOF_model(
#         model, path=outputfile_path, overwrite=True, save_data=False, engine="zarr"
#     )
#     mymodel = ecl.eof.load_REOF_model(path=outputfile_path, engine="zarr")
#     assert isinstance(mymodel, xeofs.single.EOFRotator)


def test_get_MCA_model_and_calc_MCA_analysis():
    t2m = xr.tutorial.load_dataset("air_temperature")["air"]
    da1 = t2m.isel(lon=slice(0, 26))
    da2 = t2m.isel(lon=slice(27, None))

    mca_model = ecl.eof.get_MCA_model(
        da1,
        da2,
        lat_dim="lat",
        lon_dim="lon",
        n_modes=10,
        use_coslat=True,
        solver="full",
        random_state=0,
    )
    result_datatree = ecl.eof.calc_MCA_analysis(mca_model)

    result_data1 = (
        result_datatree["EOF"]["left_EOF"]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )
    result_data2 = (
        result_datatree["EOF"]["right_EOF"]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )
    result_data3 = (
        result_datatree["PC"]["left_PC"]
        .isel(mode=0, time=slice(None, 20))
        .data.flatten()
    )
    result_data4 = (
        result_datatree["PC"]["right_PC"]
        .isel(mode=0, time=slice(None, 20))
        .data.flatten()
    )
    result_data5 = result_datatree["correlation_coefficients_X"].data.flatten()
    result_data6 = result_datatree["correlation_coefficients_Y"].data.flatten()
    result_data7 = result_datatree["covariance_fraction"].data.flatten()
    result_data8 = result_datatree["cross_correlation_coefficients"].data.flatten()
    result_data9 = result_datatree["fraction_variance_X_explained_by_X"].data.flatten()
    result_data10 = result_datatree["fraction_variance_Y_explained_by_X"].data.flatten()
    result_data11 = result_datatree["fraction_variance_Y_explained_by_Y"].data.flatten()
    result_data12 = result_datatree["squared_covariance_fraction"].data.flatten()
    result_data13 = (
        result_datatree["heterogeneous_patterns/left_heterogeneous_patterns"]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )
    result_data14 = (
        result_datatree["heterogeneous_patterns/right_heterogeneous_patterns"]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )
    result_data15 = (
        result_datatree["heterogeneous_patterns/pvalues_of_left_heterogeneous_patterns"]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )
    result_data16 = (
        result_datatree[
            "heterogeneous_patterns/pvalues_of_right_heterogeneous_patterns"
        ]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )
    result_data17 = (
        result_datatree["homogeneous_patterns/left_homogeneous_patterns"]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )
    result_data18 = (
        result_datatree["homogeneous_patterns/right_homogeneous_patterns"]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )
    result_data19 = (
        result_datatree["homogeneous_patterns/pvalues_of_left_homogeneous_patterns"]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )
    result_data20 = (
        result_datatree["homogeneous_patterns/pvalues_of_right_homogeneous_patterns"]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )

    refer_data1 = np.array(
        [0.01558326, 0.01405792, 0.01504114, 0.01474562, 0.01370998, 0.01430686]
    )
    refer_data2 = np.array(
        [0.05806844, 0.05429851, 0.0604402, 0.06080861, 0.06308315, 0.0634863]
    )
    refer_data3 = np.array(
        [
            -0.01965964,
            -0.02188687,
            -0.02300512,
            -0.02149792,
            -0.01936105,
            -0.02230554,
            -0.02288072,
            -0.021101,
            -0.01745216,
            -0.01957284,
            -0.01852802,
            -0.01836741,
            -0.01614946,
            -0.0205901,
            -0.02224295,
            -0.02124052,
            -0.01838362,
            -0.02135036,
            -0.02138998,
            -0.0209474,
        ]
    )
    refer_data4 = np.array(
        [
            -0.02089254,
            -0.02286041,
            -0.02326339,
            -0.02265824,
            -0.02363329,
            -0.02712326,
            -0.02675353,
            -0.0261337,
            -0.0254113,
            -0.02813793,
            -0.02868098,
            -0.02645875,
            -0.02586578,
            -0.02594402,
            -0.02520787,
            -0.02279183,
            -0.02296499,
            -0.02568254,
            -0.02591385,
            -0.02441013,
        ]
    )
    refer_data5 = np.array(
        [
            1.00034258,
            0.05091213,
            0.0477693,
            -0.04618039,
            0.08341818,
            0.01019111,
            0.03106142,
            0.08745134,
            -0.06169519,
            0.0323779,
            0.05091213,
            1.00034258,
            -0.2698223,
            -0.04627493,
            0.11332768,
            0.12340422,
            0.28155392,
            0.11583862,
            -0.0025707,
            -0.10371667,
            0.0477693,
            -0.2698223,
            1.00034258,
            0.11594853,
            -0.25484756,
            -0.24710089,
            -0.43387763,
            0.0381828,
            0.24572371,
            0.3747766,
            -0.04618039,
            -0.04627493,
            0.11594853,
            1.00034258,
            -0.17266141,
            0.38627309,
            -0.11352632,
            -0.00329691,
            0.27015383,
            0.09440547,
            0.08341818,
            0.11332768,
            -0.25484756,
            -0.17266141,
            1.00034258,
            0.03187743,
            0.1435686,
            0.33747335,
            -0.27991527,
            -0.17946749,
            0.01019111,
            0.12340422,
            -0.24710089,
            0.38627309,
            0.03187743,
            1.00034258,
            0.14073894,
            -0.04791318,
            -0.06978607,
            -0.26116732,
            0.03106142,
            0.28155392,
            -0.43387763,
            -0.11352632,
            0.1435686,
            0.14073894,
            1.00034258,
            0.00306663,
            -0.25236892,
            -0.28796784,
            0.08745134,
            0.11583862,
            0.0381828,
            -0.00329691,
            0.33747335,
            -0.04791318,
            0.00306663,
            1.00034258,
            -0.06696252,
            -0.03856645,
            -0.06169519,
            -0.0025707,
            0.24572371,
            0.27015383,
            -0.27991527,
            -0.06978607,
            -0.25236892,
            -0.06696252,
            1.00034258,
            0.14860438,
            0.0323779,
            -0.10371667,
            0.3747766,
            0.09440547,
            -0.17946749,
            -0.26116732,
            -0.28796784,
            -0.03856645,
            0.14860438,
            1.00034258,
        ]
    )
    refer_data6 = np.array(
        [
            1.00034258,
            -0.03323755,
            -0.06621311,
            0.02149172,
            -0.04363363,
            -0.05428197,
            -0.00640223,
            0.06088193,
            -0.06661587,
            0.11255033,
            -0.03323755,
            1.00034258,
            0.02362015,
            0.05605379,
            -0.03993013,
            -0.18056184,
            -0.33470027,
            -0.02887211,
            -0.1287425,
            -0.05885655,
            -0.06621311,
            0.02362015,
            1.00034258,
            0.13241564,
            0.27994298,
            0.37009574,
            -0.02332147,
            -0.10510192,
            0.1352765,
            0.03899048,
            0.02149172,
            0.05605379,
            0.13241564,
            1.00034258,
            0.16297803,
            -0.05573436,
            -0.17095611,
            -0.23545458,
            0.01153948,
            -0.1293934,
            -0.04363363,
            -0.03993013,
            0.27994298,
            0.16297803,
            1.00034258,
            0.01784411,
            -0.11589252,
            -0.01809376,
            0.10523593,
            -0.06055696,
            -0.05428197,
            -0.18056184,
            0.37009574,
            -0.05573436,
            0.01784411,
            1.00034258,
            -0.00906442,
            -0.31905356,
            0.43199383,
            0.15733733,
            -0.00640223,
            -0.33470027,
            -0.02332147,
            -0.17095611,
            -0.11589252,
            -0.00906442,
            1.00034258,
            0.07692467,
            0.02321068,
            0.15483428,
            0.06088193,
            -0.02887211,
            -0.10510192,
            -0.23545458,
            -0.01809376,
            -0.31905356,
            0.07692467,
            1.00034258,
            -0.25836711,
            0.15994814,
            -0.06661587,
            -0.1287425,
            0.1352765,
            0.01153948,
            0.10523593,
            0.43199383,
            0.02321068,
            -0.25836711,
            1.00034258,
            0.03261397,
            0.11255033,
            -0.05885655,
            0.03899048,
            -0.1293934,
            -0.06055696,
            0.15733733,
            0.15483428,
            0.15994814,
            0.03261397,
            1.00034258,
        ]
    )
    refer_data7 = np.array(
        [
            9.28776726e-01,
            2.81577930e-02,
            1.50080474e-02,
            9.27542244e-03,
            7.68112471e-03,
            4.28324636e-03,
            2.74762904e-03,
            2.14268177e-03,
            1.20557269e-03,
            7.21756237e-04,
        ]
    )
    refer_data8 = np.array(
        [
            0.9641005,
            0.7657287,
            0.49969426,
            0.45858186,
            0.41303811,
            0.23081255,
            0.16580367,
            0.18137918,
            0.09319696,
            0.06168472,
        ]
    )
    refer_data9 = np.array(
        [
            0.82936789,
            0.02959719,
            0.03261848,
            0.02330888,
            0.01837483,
            0.01807617,
            0.01839542,
            0.01028853,
            0.01073144,
            0.00924117,
        ]
    )
    refer_data10 = np.array(
        [
            9.51870941e-01,
            2.78700018e-02,
            9.39324447e-03,
            4.63173060e-03,
            3.76155974e-03,
            1.29245449e-03,
            4.94903686e-04,
            4.69950767e-04,
            1.51283996e-04,
            6.39289999e-05,
        ]
    )
    refer_data11 = np.array(
        [
            0.85637781,
            0.03496484,
            0.02116462,
            0.01343219,
            0.01440391,
            0.01457993,
            0.01142491,
            0.01038058,
            0.01193328,
            0.01133793,
        ]
    )
    refer_data12 = np.array(
        [
            9.98615922e-01,
            9.17852850e-04,
            2.60749943e-04,
            9.95963067e-05,
            6.83007496e-05,
            2.12384074e-05,
            8.73960962e-06,
            5.31485215e-06,
            1.68252932e-06,
            6.03055008e-07,
        ]
    )
    refer_data13 = np.array(
        [0.81719904, 0.78632968, 0.79406758, 0.78797284, 0.76287803, 0.77685482]
    )
    refer_data14 = np.array(
        [0.88580087, 0.90172098, 0.88207484, 0.89022451, 0.88250947, 0.88197606]
    )
    refer_data15 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    refer_data16 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    refer_data17 = np.array(
        [0.7763152, 0.74608813, 0.75346575, 0.75260881, 0.72411332, 0.74683124]
    )
    refer_data18 = np.array(
        [0.9298704, 0.94731727, 0.93471545, 0.94020059, 0.94094888, 0.93843361]
    )
    refer_data19 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    refer_data20 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    assert np.isclose(result_data1, refer_data1, atol=0.01).all()
    assert np.isclose(result_data2, refer_data2, atol=0.01).all()
    assert np.isclose(result_data3, refer_data3, atol=0.01).all()
    assert np.isclose(result_data4, refer_data4, atol=0.01).all()
    assert np.isclose(result_data5, refer_data5, atol=0.01).all()
    assert np.isclose(result_data6, refer_data6, atol=0.01).all()
    assert np.isclose(result_data7, refer_data7, atol=0.01).all()
    assert np.isclose(result_data8, refer_data8, atol=0.01).all()
    assert np.isclose(result_data9, refer_data9, atol=0.01).all()
    assert np.isclose(result_data10, refer_data10, atol=0.01).all()
    assert np.isclose(result_data11, refer_data11, atol=0.01).all()
    assert np.isclose(result_data12, refer_data12, atol=0.01).all()
    assert np.isclose(result_data13, refer_data13, atol=0.01).all()
    assert np.isclose(result_data14, refer_data14, atol=0.01).all()
    assert np.isclose(result_data15, refer_data15, atol=0.01).all()
    assert np.isclose(result_data16, refer_data16, atol=0.01).all()
    assert np.isclose(result_data17, refer_data17, atol=0.01).all()
    assert np.isclose(result_data18, refer_data18, atol=0.01).all()
    assert np.isclose(result_data19, refer_data19, atol=0.01).all()
    assert np.isclose(result_data20, refer_data20, atol=0.01).all()


# def test_save_MCA_model_and_load_MCA_model():
#     assert_path_dir_exist(TEST_TMP_PATH)

#     t2m = xr.tutorial.load_dataset("air_temperature")["air"]
#     da1 = t2m.isel(lon=slice(0, 26))
#     da2 = t2m.isel(lon=slice(27, None))
#     mca_model = ecl.eof.get_MCA_model(
#         da1, da2, lat_dim="lat", lon_dim="lon", n_modes=2, use_coslat=True
#     )
#     result_datatree = ecl.eof.calc_MCA_analysis(mca_model)

#     outputfile_path = str(Path(TEST_TMP_PATH, "test_output_core_mca_model.zarr"))
#     ecl.eof.save_MCA_model(
#         mca_model, path=outputfile_path, overwrite=True, save_data=False, engine="zarr"
#     )
#     mymodel = ecl.eof.load_MCA_model(path=outputfile_path, engine="zarr")
#     assert isinstance(mymodel, xeofs.cross.MCA)


def test_calc_eof_projection_coefficient():
    rng = np.random.default_rng(42)
    field = xr.DataArray(rng.random((2, 3)), dims=["lat", "lon"])
    eof_v = xr.DataArray(rng.random((2, 3)), dims=["lat", "lon"])
    result_data1 = calc_eof_projection_coefficient(field, eof_v)

    time = xr.DataArray(np.arange(4), dims=["time"])
    timed_field = xr.DataArray(rng.random((4, 2, 3)), dims=["time", "lat", "lon"])
    result_data2 = calc_eof_projection_coefficient(timed_field, eof_v)

    refer_data1 = np.array([0.95208032])
    refer_data2 = np.array([0.64684219, 1.06549741, 0.62797062, 0.78151191])
    assert np.isclose(result_data1, refer_data1, atol=0.01).all()
    assert np.isclose(result_data2, refer_data2, atol=0.01).all()


def test_clean_tmp_file():
    # Delete temp file
    import shutil

    assert_path_dir_exist(TEST_TMP_PATH)
    shutil.rmtree(TEST_TMP_PATH)
