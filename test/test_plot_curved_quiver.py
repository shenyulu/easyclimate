"""
pytest for plot.velovect.py
"""

import pytest

import easyclimate as ecl
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


@pytest.fixture(scope="module")
def uvdata():
    udata = ecl.open_tutorial_dataset("uwnd_2022_day5")["uwnd"].sortby("lat")
    vdata = ecl.open_tutorial_dataset("vwnd_2022_day5")["vwnd"].sortby("lat")
    uvdata = xr.Dataset(data_vars={"u": udata, "v": vdata})
    return uvdata


class TestCurvedQuiver:
    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_curved_quiver1(self, uvdata):
        lonlat_range = [90, 180, 0, 70]
        gap_value = 5

        drawdata_gap1 = (
            uvdata.isel(time=2)
            .sel(level=850)
            .sel(
                lon=slice(lonlat_range[0] - gap_value, lonlat_range[1] + gap_value),
                lat=slice(lonlat_range[2] - gap_value, lonlat_range[3] + gap_value),
            )
        )

        fig, ax = plt.subplots(
            subplot_kw={"projection": ccrs.PlateCarree(central_longitude=120)}
        )

        ax.coastlines()
        ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=["bottom", "left"], alpha=0.4)
        ax.set_extent(lonlat_range, crs=ccrs.PlateCarree())

        cq = ecl.plot.curved_quiver(
            drawdata_gap1,
            x="lon",
            y="lat",
            u="u",
            v="v",
            ax=ax,
            density=5,
            color="#f6631c",
            transform=ccrs.PlateCarree(),
        )

        ecl.plot.add_curved_quiverkey(
            cq,
            ax=ax,
            pos=(0.9, 1.05),
            U=10,
            label="10",
            color="black",
            labelpos="N",
            fontproperties={"size": 12},
            ref_point=(120, 30),
        )

        ax.set_title("")
        return fig

    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_curved_quiver2(self, uvdata):
        lonlat_range = [90, 180, 0, 70]
        gap_value = 5

        drawdata_gap1 = (
            uvdata.isel(time=2)
            .sel(level=850)
            .sel(
                lon=slice(lonlat_range[0] - gap_value, lonlat_range[1] + gap_value),
                lat=slice(lonlat_range[2] - gap_value, lonlat_range[3] + gap_value),
            )
        )

        fig, ax = plt.subplots()

        cq = ecl.plot.curved_quiver(
            drawdata_gap1,
            x="lon",
            y="lat",
            u="u",
            v="v",
            ax=ax,
            density=5,
            color="#f6631c",
        )

        ecl.plot.add_curved_quiverkey(
            cq,
            ax=ax,
            pos=(0.9, 1.05),
            U=10,
            label="10",
            color="black",
            labelpos="N",
            fontproperties={"size": 12},
            ref_point=(120, 30),
        )

        ax.set_title("")
        return fig
