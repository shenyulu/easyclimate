"""
pytest for plot.velovect.py
"""

import pytest

import easyclimate as ecl
import xarray as xr
import cartopy.crs as ccrs

udata = ecl.open_tutorial_dataset("uwnd_202201_mon_mean").uwnd
vdata = ecl.open_tutorial_dataset("vwnd_202201_mon_mean").vwnd
uvdata = xr.Dataset(data_vars={"u": udata, "v": vdata}).sortby("lat")
uvdata850 = uvdata.isel(time=0).sel(lon=slice(90, 180), lat=slice(0, 70), level=850)


@pytest.mark.mpl_image_compare(remove_text=True)
def test_curved_quiver():
    fig, ax = ecl.plot.quick_draw_spatial_basemap(central_longitude=180)
    cq = uvdata850.easyclimate.plot.curved_quiver(
        x="lon",
        y="lat",
        u="u",
        v="v",
        ax=ax,
        density=2,
        color="mediumpurple",
        transform=ccrs.PlateCarree(),
    )
    ecl.plot.add_curved_quiverkey(
        cq, 0.95, 1.02, 3, label="", color="mediumpurple", labelpos="N", labelsep=0.04
    )
    return fig
