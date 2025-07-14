"""
pytest for interp/vinth2p.py

Part1
"""

import pytest

from easyclimate.interp.vinth2p import interp_vinth2p_dp, interp_vinth2p_ecmwf
import numpy as np
import xarray as xr


@pytest.fixture
def sample_data1():
    t_cesm_ndata = xr.DataArray(
        np.array(
            [
                [[240.38052]],
                [[229.63916]],
                [[223.04153]],
                [[218.12868]],
                [[209.68248]],
                [[205.21547]],
                [[195.03557]],
                [[183.25945]],
                [[187.07957]],
                [[193.25035]],
                [[200.69344]],
                [[208.84563]],
                [[217.64043]],
                [[226.82707]],
                [[236.31105]],
                [[245.6203]],
                [[254.20686]],
                [[261.92285]],
                [[269.39594]],
                [[274.62244]],
                [[280.0596]],
                [[284.87875]],
                [[288.72394]],
                [[291.2716]],
                [[292.85422]],
                [[294.20374]],
                [[295.60782]],
                [[297.20935]],
                [[298.7263]],
                [[299.97076]],
            ],
            dtype=np.float32,
        ),
        dims=("lev", "lat", "lon"),
        coords={
            "lev": np.array(
                [
                    3.64346569,
                    7.59481965,
                    14.35663225,
                    24.61222,
                    38.26829977,
                    54.59547974,
                    72.01245055,
                    87.82123029,
                    103.31712663,
                    121.54724076,
                    142.99403876,
                    168.22507977,
                    197.9080867,
                    232.82861896,
                    273.91081676,
                    322.24190235,
                    379.10090387,
                    445.9925741,
                    524.68717471,
                    609.77869481,
                    691.38943031,
                    763.40448111,
                    820.85836865,
                    859.53476653,
                    887.02024892,
                    912.64454694,
                    936.19839847,
                    957.48547954,
                    976.32540739,
                    992.55609512,
                ]
            ),
            "lat": np.array([0.47120419]),
            "lon": np.array([150.0]),
        },
    )

    z_cesm_ndata = xr.DataArray(
        np.array(
            [
                [[37682.703]],
                [[32597.73]],
                [[28347.734]],
                [[24839.229]],
                [[22050.93]],
                [[19877.744]],
                [[18240.92]],
                [[17129.004]],
                [[16247.974]],
                [[15343.029]],
                [[14405.623]],
                [[13431.077]],
                [[12412.256]],
                [[11347.556]],
                [[10239.184]],
                [[9086.737]],
                [[7892.1274]],
                [[6658.683]],
                [[5388.1597]],
                [[4183.013]],
                [[3155.8298]],
                [[2330.3606]],
                [[1716.3195]],
                [[1322.1869]],
                [[1050.8429]],
                [[803.8812]],
                [[581.6391]],
                [[384.38214]],
                [[212.45535]],
                [[66.37027]],
            ],
            dtype=np.float32,
        ),
        dims=("lev", "lat", "lon"),
        coords={
            "lev": np.array(
                [
                    3.64346569,
                    7.59481965,
                    14.35663225,
                    24.61222,
                    38.26829977,
                    54.59547974,
                    72.01245055,
                    87.82123029,
                    103.31712663,
                    121.54724076,
                    142.99403876,
                    168.22507977,
                    197.9080867,
                    232.82861896,
                    273.91081676,
                    322.24190235,
                    379.10090387,
                    445.9925741,
                    524.68717471,
                    609.77869481,
                    691.38943031,
                    763.40448111,
                    820.85836865,
                    859.53476653,
                    887.02024892,
                    912.64454694,
                    936.19839847,
                    957.48547954,
                    976.32540739,
                    992.55609512,
                ]
            ),
            "lat": np.array([0.47120419]),
            "lon": np.array([150.0]),
        },
    )

    msl_cesm_ndata = xr.DataArray(
        np.array([[100635.1]], dtype=np.float32),
        dims=("lat", "lon"),
        coords={"lat": np.array([0.47120419]), "lon": np.array([150.0])},
    )

    phis_cesm_ndata = xr.DataArray(
        np.array([[2.4943347e-14]], dtype=np.float32),
        dims=("lat", "lon"),
        coords={"lat": np.array([0.47120419]), "lon": np.array([150.0])},
    )

    my_hybrid_An_coefficients = xr.DataArray(
        np.array(
            [
                0.00364347,
                0.00759482,
                0.01435663,
                0.02461222,
                0.0382683,
                0.05459548,
                0.07201245,
                0.08782123,
                0.10331713,
                0.12154724,
                0.14299404,
                0.16822508,
                0.17823067,
                0.17032433,
                0.16102291,
                0.15008029,
                0.13720686,
                0.12206194,
                0.10424471,
                0.08497915,
                0.0665017,
                0.05019679,
                0.03718866,
                0.02843195,
                0.02220898,
                0.01640738,
                0.01107456,
                0.00625495,
                0.00198941,
                0.0,
            ]
        ),
        dims=("lev"),
        coords={
            "lev": np.array(
                [
                    3.64346569,
                    7.59481965,
                    14.35663225,
                    24.61222,
                    38.26829977,
                    54.59547974,
                    72.01245055,
                    87.82123029,
                    103.31712663,
                    121.54724076,
                    142.99403876,
                    168.22507977,
                    197.9080867,
                    232.82861896,
                    273.91081676,
                    322.24190235,
                    379.10090387,
                    445.9925741,
                    524.68717471,
                    609.77869481,
                    691.38943031,
                    763.40448111,
                    820.85836865,
                    859.53476653,
                    887.02024892,
                    912.64454694,
                    936.19839847,
                    957.48547954,
                    976.32540739,
                    992.55609512,
                ]
            )
        },
    )

    my_hybrid_Bn_coefficients = xr.DataArray(
        np.array(
            [
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.01967741,
                0.06250429,
                0.11288791,
                0.17216162,
                0.24189404,
                0.32393064,
                0.42044246,
                0.52479954,
                0.62488773,
                0.71320769,
                0.78366971,
                0.83110282,
                0.86481127,
                0.89623716,
                0.92512384,
                0.95123053,
                0.974336,
                0.9925561,
            ]
        ),
        dims=("lev"),
        coords={
            "lev": np.array(
                [
                    3.64346569,
                    7.59481965,
                    14.35663225,
                    24.61222,
                    38.26829977,
                    54.59547974,
                    72.01245055,
                    87.82123029,
                    103.31712663,
                    121.54724076,
                    142.99403876,
                    168.22507977,
                    197.9080867,
                    232.82861896,
                    273.91081676,
                    322.24190235,
                    379.10090387,
                    445.9925741,
                    524.68717471,
                    609.77869481,
                    691.38943031,
                    763.40448111,
                    820.85836865,
                    859.53476653,
                    887.02024892,
                    912.64454694,
                    936.19839847,
                    957.48547954,
                    976.32540739,
                    992.55609512,
                ]
            )
        },
    )
    return (
        t_cesm_ndata,
        z_cesm_ndata,
        msl_cesm_ndata,
        phis_cesm_ndata,
        my_hybrid_An_coefficients,
        my_hybrid_Bn_coefficients,
    )


def test_interp_vinth2p_dp1(sample_data1):
    (
        t_cesm_ndata,
        z_cesm_ndata,
        msl_cesm_ndata,
        phis_cesm_ndata,
        my_hybrid_An_coefficients,
        my_hybrid_Bn_coefficients,
    ) = sample_data1
    result_data = interp_vinth2p_dp(
        data_input=t_cesm_ndata,
        surface_pressure_data=msl_cesm_ndata,
        surface_pressure_data_units="Pa",
        hybrid_A_coefficients=my_hybrid_An_coefficients,
        hybrid_B_coefficients=my_hybrid_Bn_coefficients,
        vertical_input_dim="lev",
        vertical_output_dim="plev",
        vertical_output_dim_units="hPa",
        vertical_output_level=np.array(
            [1000.0, 900.0, 800.0, 700.0, 600.0, 500.0, 400.0, 100.0]
        ),
    ).data.flatten()
    refer_data = np.array(
        [
            np.nan,
            293.24549387,
            287.00822564,
            280.3678298,
            273.82331246,
            266.81806934,
            256.42311278,
            186.26181855,
        ]
    )
    assert np.isclose(result_data, refer_data, equal_nan=True).all()


def test_interp_vinth2p_ecmwf1(sample_data1):
    (
        t_cesm_ndata,
        z_cesm_ndata,
        msl_cesm_ndata,
        phis_cesm_ndata,
        my_hybrid_An_coefficients,
        my_hybrid_Bn_coefficients,
    ) = sample_data1
    result_data = interp_vinth2p_ecmwf(
        data_input=t_cesm_ndata,
        surface_pressure_data=msl_cesm_ndata,
        surface_pressure_data_units="Pa",
        hybrid_A_coefficients=my_hybrid_An_coefficients,
        hybrid_B_coefficients=my_hybrid_Bn_coefficients,
        vertical_input_dim="lev",
        vertical_output_dim="plev",
        vertical_output_dim_units="hPa",
        vertical_output_level=np.array(
            [1000.0, 900.0, 800.0, 700.0, 600.0, 500.0, 400.0]
        ),
        variable_flag="T",
        temperature_bottom_data=t_cesm_ndata.isel(lev=-1),
        surface_geopotential_data=phis_cesm_ndata,
    ).data.flatten()
    refer_data = np.array(
        [
            300.03717378,
            293.24549387,
            287.00822564,
            280.3678298,
            273.82331246,
            266.81806934,
            256.42311278,
        ]
    )
    assert np.isclose(result_data, refer_data, equal_nan=True).all()
