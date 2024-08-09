"""
pytest for plot\quick_draw.py
"""

import pytest

import easyclimate as ecl
import cartopy.crs as ccrs

t_data_origin = (
    ecl.tutorial.open_tutorial_dataset("air_202201_mon_mean")
    .air.isel(time=0)
    .sel(level=1000)
)


@pytest.mark.mpl_image_compare
def test_quick_draw_spatial_basemap1():
    fig, ax = ecl.plot.quick_draw_spatial_basemap(central_longitude=180)

    t_data_origin.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cbar_kwargs={"location": "bottom", "label": ""},
        levels=21,
    )
    ax.set_title("")
    return fig


@pytest.mark.mpl_image_compare
def test_quick_draw_spatial_basemap2():
    fig, ax = ecl.plot.quick_draw_spatial_basemap(nrows=2, central_longitude=180)
    fig.subplots_adjust(hspace=0.3)

    # Figure 1
    t_data_origin.plot.contourf(
        ax=ax[0],
        transform=ccrs.PlateCarree(),
        cbar_kwargs={"location": "right", "label": ""},
        levels=21,
    )
    ax[0].set_title("")

    # Figure 2
    t_data_origin.plot.contourf(
        ax=ax[1],
        transform=ccrs.PlateCarree(),
        cbar_kwargs={"location": "right", "label": ""},
        levels=11,
    )
    ax[1].set_title("")
    return fig


@pytest.mark.mpl_image_compare
def test_quick_draw_rectangular_box1():
    fig, ax = ecl.plot.quick_draw_spatial_basemap(central_longitude=180)
    ax.set_extent([120, 290, -30, 30], crs=ccrs.PlateCarree())
    ecl.plot.quick_draw_rectangular_box(
        -120, 170, 5, -5, ax, ec="r", fc="none", lw=2, transform=ccrs.PlateCarree()
    )
    return fig
