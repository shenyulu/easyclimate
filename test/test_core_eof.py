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

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()


def test_save_EOF_model_and_load_EOF_model():
    assert_path_dir_exist(TEST_TMP_PATH)
    inputdata = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_core_eof.nc"))).z
    model = ecl.eof.get_EOF_model(
        inputdata,
        lat_dim="lat",
        lon_dim="lon",
        remove_seasonal_cycle_mean=True,
        use_coslat=True,
    )
    result_data = ecl.eof.calc_EOF_analysis(model)
    outputfile_path = str(Path(TEST_TMP_PATH, "test_output_core_eof_model.zarr"))
    ecl.eof.save_EOF_model(
        model, path=outputfile_path, overwrite=True, save_data=False, engine="zarr"
    )
    mymodel = ecl.eof.load_EOF_model(path=outputfile_path, engine="zarr")
    assert isinstance(mymodel, xeofs.models.eof.EOF)


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

    assert np.isclose(result_data1, refer_data1).all()
    assert np.isclose(result_data2, refer_data2).all()
    assert np.isclose(result_data3, refer_data3).all()
    assert np.isclose(result_data4, refer_data4).all()
    assert np.isclose(result_data5, refer_data5).all()


def test_save_EOF_model_and_load_EOF_model():
    assert_path_dir_exist(TEST_TMP_PATH)
    inputdata = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_core_eof.nc"))).z
    model = ecl.eof.get_REOF_model(
        inputdata,
        lat_dim="lat",
        lon_dim="lon",
        remove_seasonal_cycle_mean=True,
        use_coslat=True,
    )
    result_data = ecl.eof.calc_REOF_analysis(model)
    outputfile_path = str(Path(TEST_TMP_PATH, "test_output_core_reof_model.zarr"))
    ecl.eof.save_REOF_model(
        model, path=outputfile_path, overwrite=True, save_data=False, engine="zarr"
    )
    mymodel = ecl.eof.load_REOF_model(path=outputfile_path, engine="zarr")
    assert isinstance(mymodel, xeofs.models.EOFRotator)


def test_get_MCA_model_and_calc_MCA_analysis():
    t2m = xr.tutorial.load_dataset("air_temperature")["air"]
    da1 = t2m.isel(lon=slice(0, 26))
    da2 = t2m.isel(lon=slice(27, None))

    mca_model = ecl.eof.get_MCA_model(
        da1, da2, lat_dim="lat", lon_dim="lon", n_modes=2, use_coslat=True
    )
    result_datatree = ecl.eof.calc_MCA_analysis(mca_model)

    result_data1 = (
        result_datatree["EOF"]["left_EOF"]["components1"]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )
    result_data2 = (
        result_datatree["EOF"]["right_EOF"]["components2"]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )
    result_data3 = (
        result_datatree["PC"]["left_PC"]["scores1"]
        .isel(mode=0, time=slice(None, 20))
        .data.flatten()
    )
    result_data4 = (
        result_datatree["PC"]["right_PC"]["scores2"]
        .isel(mode=0, time=slice(None, 20))
        .data.flatten()
    )
    result_data5 = result_datatree["covariance_fraction"][
        "covariance_fraction"
    ].data.flatten()
    result_data6 = result_datatree["singular_values"]["singular_values"].data.flatten()
    result_data7 = result_datatree["squared_covariance"][
        "squared_covariance"
    ].data.flatten()
    result_data8 = result_datatree["squared_covariance_fraction"][
        "squared_covariance_fraction"
    ].data.flatten()

    result_data9 = (
        result_datatree["heterogeneous_patterns/left_heterogeneous_patterns"][
            "left_heterogeneous_patterns"
        ]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )
    result_data10 = (
        result_datatree["heterogeneous_patterns/right_heterogeneous_patterns"][
            "right_heterogeneous_patterns"
        ]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )
    result_data11 = (
        result_datatree[
            "heterogeneous_patterns/pvalues_of_left_heterogeneous_patterns"
        ]["pvalues_of_left_heterogeneous_patterns"]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )
    result_data12 = (
        result_datatree[
            "heterogeneous_patterns/pvalues_of_right_heterogeneous_patterns"
        ]["pvalues_of_right_heterogeneous_patterns"]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )

    result_data13 = (
        result_datatree["homogeneous_patterns/left_homogeneous_patterns"][
            "left_homogeneous_patterns"
        ]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )
    result_data14 = (
        result_datatree["homogeneous_patterns/right_homogeneous_patterns"][
            "right_homogeneous_patterns"
        ]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )
    result_data15 = (
        result_datatree["homogeneous_patterns/pvalues_of_left_homogeneous_patterns"][
            "pvalues_of_left_homogeneous_patterns"
        ]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )
    result_data16 = (
        result_datatree["homogeneous_patterns/pvalues_of_right_homogeneous_patterns"][
            "pvalues_of_right_homogeneous_patterns"
        ]
        .isel(mode=0)
        .isel(lat=slice(10, 13), lon=slice(6, 8))
        .data.flatten()
    )

    refer_data1 = np.array(
        [0.01424639, 0.01283669, 0.01373628, 0.01355581, 0.01252516, 0.01323858]
    )
    refer_data2 = np.array(
        [0.05879741, 0.05500811, 0.06178727, 0.0619469, 0.06488786, 0.06516425]
    )
    refer_data3 = np.array(
        [
            0.00361445,
            0.00359823,
            0.00359024,
            0.00360107,
            0.00361677,
            0.00359558,
            0.00359158,
            0.00360423,
            0.00363048,
            0.00361504,
            0.0036224,
            0.00362335,
            0.00363947,
            0.00360717,
            0.00359529,
            0.00360251,
            0.00362341,
            0.00360205,
            0.00360161,
            0.00360466,
        ]
    )
    refer_data4 = np.array(
        [
            0.00387221,
            0.00385854,
            0.00385573,
            0.00385987,
            0.00385299,
            0.00382875,
            0.00383138,
            0.00383567,
            0.00384058,
            0.00382164,
            0.00381791,
            0.00383335,
            0.00383727,
            0.00383647,
            0.00384156,
            0.00385833,
            0.00385721,
            0.00383832,
            0.0038367,
            0.0038472,
        ]
    )
    refer_data5 = np.array([9.99321686e-01, 6.78314305e-04])
    refer_data6 = np.array([1.2383819e12, 8.4058234e08])
    refer_data7 = np.array([1.53358974e24, 7.06578670e17])
    refer_data8 = np.array([9.99999539e-01, 4.60734919e-07])
    refer_data9 = np.array(
        [
            2643.94120679,
            2814.18914907,
            2606.29513332,
            2660.039444,
            2692.63997955,
            2649.38745258,
        ]
    )
    refer_data10 = np.array(
        [
            766.11631185,
            839.15093393,
            715.23096131,
            717.78585558,
            667.5871293,
            661.89358203,
        ]
    )
    refer_data11 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    refer_data12 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    refer_data13 = np.array(
        [
            2383.9696907,
            2537.47378184,
            2350.02760208,
            2398.49181239,
            2427.88325588,
            2388.89479008,
        ]
    )
    refer_data14 = np.array(
        [
            849.62210588,
            930.62236392,
            793.19624281,
            796.02573919,
            740.36144092,
            734.04497586,
        ]
    )
    refer_data15 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    refer_data16 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

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
    assert np.isclose(result_data12, refer_data12).all()
    assert np.isclose(result_data13, refer_data13).all()
    assert np.isclose(result_data14, refer_data14).all()
    assert np.isclose(result_data15, refer_data15).all()
    assert np.isclose(result_data16, refer_data16).all()


def test_save_MCA_model_and_load_MCA_model():
    assert_path_dir_exist(TEST_TMP_PATH)

    t2m = xr.tutorial.load_dataset("air_temperature")["air"]
    da1 = t2m.isel(lon=slice(0, 26))
    da2 = t2m.isel(lon=slice(27, None))
    mca_model = ecl.eof.get_MCA_model(
        da1, da2, lat_dim="lat", lon_dim="lon", n_modes=2, use_coslat=True
    )
    result_datatree = ecl.eof.calc_MCA_analysis(mca_model)

    outputfile_path = str(Path(TEST_TMP_PATH, "test_output_core_mca_model.zarr"))
    ecl.eof.save_MCA_model(
        mca_model, path=outputfile_path, overwrite=True, save_data=False, engine="zarr"
    )
    mymodel = ecl.eof.load_MCA_model(path=outputfile_path, engine="zarr")
    assert isinstance(mymodel, xeofs.models.MCA)


def test_clean_tmp_file():
    # Delete temp file
    import shutil

    assert_path_dir_exist(TEST_TMP_PATH)
    shutil.rmtree(TEST_TMP_PATH)
