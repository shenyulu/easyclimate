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

data_u = xr.open_dataset(str(Path(TEST_DATA_PATH, 'test_input_interp_mesh2mesh.nc'))).uwnd

@pytest.mark.mpl_image_compare
def test_draw_Circlemap_PolarStereo1():
    fig, ax = plt.subplots(subplot_kw = {'projection': ccrs.NorthPolarStereo()})
    ax.coastlines(edgecolor = 'black', linewidths = 0.5)
    ecl.plot.draw_Circlemap_PolarStereo(
        ax = ax,
        lon_step = 30,
        lat_step = 10,
        lat_range = [50, 90],
        draw_labels = True,
        set_map_boundary_kwargs = {'north_pad': 0.3, 'south_pad': 0.4},
        gridlines_kwargs = {'color': 'grey', 'alpha': 0.5, 'linestyle' : '--'}
    )
    return fig

@pytest.mark.mpl_image_compare
def test_draw_Circlemap_PolarStereo2():
    fig, ax = plt.subplots(subplot_kw = {'projection': ccrs.NorthPolarStereo()})

    ax.coastlines(edgecolor = 'black', linewidths = 0.5)

    ecl.plot.draw_Circlemap_PolarStereo(
        ax = ax,
        lon_step = 30,
        lat_step = 10,
        lat_range = [50, 90],
        add_gridlines = False,
        draw_labels = True,
        set_map_boundary_kwargs = {'north_pad': 0.3, 'south_pad': 0.4},
        gridlines_kwargs = {'color': 'grey', 'alpha': 0.5, 'linestyle' : '--'}
    )
    return fig

@pytest.mark.mpl_image_compare
def test_draw_Circlemap_PolarStereo3():
    fig, ax = plt.subplots(subplot_kw = {'projection': ccrs.NorthPolarStereo()})
    ax.coastlines(edgecolor = 'black', linewidths = 0.5)
    ecl.plot.draw_Circlemap_PolarStereo(
        ax = ax,
        lon_step = 30,
        lat_step = 10,
        lat_range = [50, 90],
        add_gridlines = True,
        draw_labels = False,
        set_map_boundary_kwargs = {'north_pad': 0.3, 'south_pad': 0.4},
        gridlines_kwargs = {'color': 'grey', 'alpha': 0.5, 'linestyle' : '--'}
    )
    return fig

def test_draw_Circlemap_PolarStereo4():
    fig, ax = plt.subplots(subplot_kw = {'projection': ccrs.NorthPolarStereo()})
    ax.coastlines(edgecolor = 'black', linewidths = 0.5)
    with pytest.raises(ValueError):
        ecl.plot.draw_Circlemap_PolarStereo(
            ax = ax,
            lon_step = None, # Breakpoint
            lat_step = 10,
            lat_range = [50, 90],
            add_gridlines = True,
            draw_labels = True,
            set_map_boundary_kwargs = {'north_pad': 0.3, 'south_pad': 0.4},
            gridlines_kwargs = {'color': 'grey', 'alpha': 0.5, 'linestyle' : '--'}
        )
        assert 1 == 1

def test_draw_Circlemap_PolarStereo5():
    fig, ax = plt.subplots() # Breakpoint
    with pytest.raises(AttributeError):
        ecl.plot.draw_Circlemap_PolarStereo(
            ax = ax,
            lon_step = 30,
            lat_step = 10,
            lat_range = [50, 90],
            add_gridlines = False,
            draw_labels = True,
            set_map_boundary_kwargs = {'north_pad': 0.3, 'south_pad': 0.4},
            gridlines_kwargs = {'color': 'grey', 'alpha': 0.5, 'linestyle' : '--'}
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
    result_data = ecl.plot.add_lon_cyclic(data_u, inter = 2.5).lon.data
    refer_data = np.array([  0. ,   2.5,   5. ,   7.5,  10. ,  12.5,  15. ,  17.5,  20. ,
        22.5,  25. ,  27.5,  30. ,  32.5,  35. ,  37.5,  40. ,  42.5,
        45. ,  47.5,  50. ,  52.5,  55. ,  57.5,  60. ,  62.5,  65. ,
        67.5,  70. ,  72.5,  75. ,  77.5,  80. ,  82.5,  85. ,  87.5,
        90. ,  92.5,  95. ,  97.5, 100. , 102.5, 105. , 107.5, 110. ,
       112.5, 115. , 117.5, 120. , 122.5, 125. , 127.5, 130. , 132.5,
       135. , 137.5, 140. , 142.5, 145. , 147.5, 150. , 152.5, 155. ,
       157.5, 160. , 162.5, 165. , 167.5, 170. , 172.5, 175. , 177.5,
       180. , 182.5, 185. , 187.5, 190. , 192.5, 195. , 197.5, 200. ,
       202.5, 205. , 207.5, 210. , 212.5, 215. , 217.5, 220. , 222.5,
       225. , 227.5, 230. , 232.5, 235. , 237.5, 240. , 242.5, 245. ,
       247.5, 250. , 252.5, 255. , 257.5, 260. , 262.5, 265. , 267.5,
       270. , 272.5, 275. , 277.5, 280. , 282.5, 285. , 287.5, 290. ,
       292.5, 295. , 297.5, 300. , 302.5, 305. , 307.5, 310. , 312.5,
       315. , 317.5, 320. , 322.5, 325. , 327.5, 330. , 332.5, 335. ,
       337.5, 340. , 342.5, 345. , 347.5, 350. , 352.5, 355. , 357.5,
       360. ])
    assert np.isclose(result_data, refer_data).all()

def test_add_lon_cyclic2():
    data_u_180 = ecl.utility.transfer_xarray_lon_from360TO180(data_u)
    result_data = ecl.plot.add_lon_cyclic(data_u_180, inter = 2.5).lon.data
    refer_data = np.array([  0. ,   2.5,   5. ,   7.5,  10. ,  12.5,  15. ,  17.5,  20. ,
        22.5,  25. ,  27.5,  30. ,  32.5,  35. ,  37.5,  40. ,  42.5,
        45. ,  47.5,  50. ,  52.5,  55. ,  57.5,  60. ,  62.5,  65. ,
        67.5,  70. ,  72.5,  75. ,  77.5,  80. ,  82.5,  85. ,  87.5,
        90. ,  92.5,  95. ,  97.5, 100. , 102.5, 105. , 107.5, 110. ,
       112.5, 115. , 117.5, 120. , 122.5, 125. , 127.5, 130. , 132.5,
       135. , 137.5, 140. , 142.5, 145. , 147.5, 150. , 152.5, 155. ,
       157.5, 160. , 162.5, 165. , 167.5, 170. , 172.5, 175. , 177.5,
       180. , 182.5, 185. , 187.5, 190. , 192.5, 195. , 197.5, 200. ,
       202.5, 205. , 207.5, 210. , 212.5, 215. , 217.5, 220. , 222.5,
       225. , 227.5, 230. , 232.5, 235. , 237.5, 240. , 242.5, 245. ,
       247.5, 250. , 252.5, 255. , 257.5, 260. , 262.5, 265. , 267.5,
       270. , 272.5, 275. , 277.5, 280. , 282.5, 285. , 287.5, 290. ,
       292.5, 295. , 297.5, 300. , 302.5, 305. , 307.5, 310. , 312.5,
       315. , 317.5, 320. , 322.5, 325. , 327.5, 330. , 332.5, 335. ,
       337.5, 340. , 342.5, 345. , 347.5, 350. , 352.5, 355. , 357.5,
       360. ])
    assert np.isclose(result_data, refer_data).all()