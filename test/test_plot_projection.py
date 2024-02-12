"""
pytest for plot.projection.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from .const_define import TEST_DATA_PATH
from pathlib import Path

data_u = xr.open_dataset(
    str(Path(TEST_DATA_PATH, "test_input_interp_mesh2mesh.nc"))
).uwnd


@pytest.mark.mpl_image_compare
def test_draw_Circlemap_PolarStereo1():
    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.NorthPolarStereo()})
    ax.coastlines(edgecolor="black", linewidths=0.5)
    ecl.plot.draw_Circlemap_PolarStereo(
        ax=ax,
        lon_step=30,
        lat_step=10,
        lat_range=[50, 90],
        draw_labels=True,
        set_map_boundary_kwargs={"north_pad": 0.3, "south_pad": 0.4},
        gridlines_kwargs={"color": "grey", "alpha": 0.5, "linestyle": "--"},
    )
    return fig


@pytest.mark.mpl_image_compare
def test_draw_Circlemap_PolarStereo2():
    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.NorthPolarStereo()})

    ax.coastlines(edgecolor="black", linewidths=0.5)

    ecl.plot.draw_Circlemap_PolarStereo(
        ax=ax,
        lon_step=30,
        lat_step=10,
        lat_range=[50, 90],
        add_gridlines=False,
        draw_labels=True,
        set_map_boundary_kwargs={"north_pad": 0.3, "south_pad": 0.4},
        gridlines_kwargs={"color": "grey", "alpha": 0.5, "linestyle": "--"},
    )
    return fig


@pytest.mark.mpl_image_compare
def test_draw_Circlemap_PolarStereo3():
    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.NorthPolarStereo()})
    ax.coastlines(edgecolor="black", linewidths=0.5)
    ecl.plot.draw_Circlemap_PolarStereo(
        ax=ax,
        lon_step=30,
        lat_step=10,
        lat_range=[50, 90],
        add_gridlines=True,
        draw_labels=False,
        set_map_boundary_kwargs={"north_pad": 0.3, "south_pad": 0.4},
        gridlines_kwargs={"color": "grey", "alpha": 0.5, "linestyle": "--"},
    )
    return fig


def test_draw_Circlemap_PolarStereo4():
    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.NorthPolarStereo()})
    ax.coastlines(edgecolor="black", linewidths=0.5)
    with pytest.raises(ValueError):
        ecl.plot.draw_Circlemap_PolarStereo(
            ax=ax,
            lon_step=None,  # Breakpoint
            lat_step=10,
            lat_range=[50, 90],
            add_gridlines=True,
            draw_labels=True,
            set_map_boundary_kwargs={"north_pad": 0.3, "south_pad": 0.4},
            gridlines_kwargs={"color": "grey", "alpha": 0.5, "linestyle": "--"},
        )
        assert 1 == 1


def test_draw_Circlemap_PolarStereo5():
    fig, ax = plt.subplots()  # Breakpoint
    with pytest.raises(AttributeError):
        ecl.plot.draw_Circlemap_PolarStereo(
            ax=ax,
            lon_step=30,
            lat_step=10,
            lat_range=[50, 90],
            add_gridlines=False,
            draw_labels=True,
            set_map_boundary_kwargs={"north_pad": 0.3, "south_pad": 0.4},
            gridlines_kwargs={"color": "grey", "alpha": 0.5, "linestyle": "--"},
        )
        assert 1 == 1


# def test_draw_Circlemap_PolarStereo6():
#     fig, ax = plt.subplots(subplot_kw = {'projection': ccrs.PlateCarree()}) # Breakpoint
#     with pytest.raises(TypeError):
#         ecl.plot.draw_Circlemap_PolarStereo(
#             ax = ax,
#             lon_step = 30,
#             lat_step = 10,
#             lat_range = [50, 90],
#             add_gridlines = False,
#             draw_labels = True,
#             set_map_boundary_kwargs = {'north_pad': 0.3, 'south_pad': 0.4},
#             gridlines_kwargs = {'color': 'grey', 'alpha': 0.5, 'linestyle' : '--'}
#         )
#         assert 1 == 1


def test_add_lon_cyclic1():
    result_data = ecl.plot.add_lon_cyclic(data_u, inter=2.5).lon.data
    refer_data = np.array(
        [
            0.0,
            2.5,
            5.0,
            7.5,
            10.0,
            12.5,
            15.0,
            17.5,
            20.0,
            22.5,
            25.0,
            27.5,
            30.0,
            32.5,
            35.0,
            37.5,
            40.0,
            42.5,
            45.0,
            47.5,
            50.0,
            52.5,
            55.0,
            57.5,
            60.0,
            62.5,
            65.0,
            67.5,
            70.0,
            72.5,
            75.0,
            77.5,
            80.0,
            82.5,
            85.0,
            87.5,
            90.0,
            92.5,
            95.0,
            97.5,
            100.0,
            102.5,
            105.0,
            107.5,
            110.0,
            112.5,
            115.0,
            117.5,
            120.0,
            122.5,
            125.0,
            127.5,
            130.0,
            132.5,
            135.0,
            137.5,
            140.0,
            142.5,
            145.0,
            147.5,
            150.0,
            152.5,
            155.0,
            157.5,
            160.0,
            162.5,
            165.0,
            167.5,
            170.0,
            172.5,
            175.0,
            177.5,
            180.0,
            182.5,
            185.0,
            187.5,
            190.0,
            192.5,
            195.0,
            197.5,
            200.0,
            202.5,
            205.0,
            207.5,
            210.0,
            212.5,
            215.0,
            217.5,
            220.0,
            222.5,
            225.0,
            227.5,
            230.0,
            232.5,
            235.0,
            237.5,
            240.0,
            242.5,
            245.0,
            247.5,
            250.0,
            252.5,
            255.0,
            257.5,
            260.0,
            262.5,
            265.0,
            267.5,
            270.0,
            272.5,
            275.0,
            277.5,
            280.0,
            282.5,
            285.0,
            287.5,
            290.0,
            292.5,
            295.0,
            297.5,
            300.0,
            302.5,
            305.0,
            307.5,
            310.0,
            312.5,
            315.0,
            317.5,
            320.0,
            322.5,
            325.0,
            327.5,
            330.0,
            332.5,
            335.0,
            337.5,
            340.0,
            342.5,
            345.0,
            347.5,
            350.0,
            352.5,
            355.0,
            357.5,
            360.0,
        ]
    )
    assert np.isclose(result_data, refer_data).all()


def test_add_lon_cyclic2():
    data_u_180 = ecl.utility.transfer_xarray_lon_from360TO180(data_u)
    result_data = ecl.plot.add_lon_cyclic(data_u_180, inter=2.5).lon.data
    refer_data = np.array(
        [
            0.0,
            2.5,
            5.0,
            7.5,
            10.0,
            12.5,
            15.0,
            17.5,
            20.0,
            22.5,
            25.0,
            27.5,
            30.0,
            32.5,
            35.0,
            37.5,
            40.0,
            42.5,
            45.0,
            47.5,
            50.0,
            52.5,
            55.0,
            57.5,
            60.0,
            62.5,
            65.0,
            67.5,
            70.0,
            72.5,
            75.0,
            77.5,
            80.0,
            82.5,
            85.0,
            87.5,
            90.0,
            92.5,
            95.0,
            97.5,
            100.0,
            102.5,
            105.0,
            107.5,
            110.0,
            112.5,
            115.0,
            117.5,
            120.0,
            122.5,
            125.0,
            127.5,
            130.0,
            132.5,
            135.0,
            137.5,
            140.0,
            142.5,
            145.0,
            147.5,
            150.0,
            152.5,
            155.0,
            157.5,
            160.0,
            162.5,
            165.0,
            167.5,
            170.0,
            172.5,
            175.0,
            177.5,
            180.0,
            182.5,
            185.0,
            187.5,
            190.0,
            192.5,
            195.0,
            197.5,
            200.0,
            202.5,
            205.0,
            207.5,
            210.0,
            212.5,
            215.0,
            217.5,
            220.0,
            222.5,
            225.0,
            227.5,
            230.0,
            232.5,
            235.0,
            237.5,
            240.0,
            242.5,
            245.0,
            247.5,
            250.0,
            252.5,
            255.0,
            257.5,
            260.0,
            262.5,
            265.0,
            267.5,
            270.0,
            272.5,
            275.0,
            277.5,
            280.0,
            282.5,
            285.0,
            287.5,
            290.0,
            292.5,
            295.0,
            297.5,
            300.0,
            302.5,
            305.0,
            307.5,
            310.0,
            312.5,
            315.0,
            317.5,
            320.0,
            322.5,
            325.0,
            327.5,
            330.0,
            332.5,
            335.0,
            337.5,
            340.0,
            342.5,
            345.0,
            347.5,
            350.0,
            352.5,
            355.0,
            357.5,
            360.0,
        ]
    )
    assert np.isclose(result_data, refer_data).all()
