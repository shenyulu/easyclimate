"""
pytest for filter/mask.py
"""

import matplotlib.pyplot as plt
import pytest
import xarray as xr

import cartopy.crs as ccrs
import easyclimate as ecl


def _build_test_dataarray():
    lon = xr.DataArray(range(260, 351), dims="lon", name="lon")
    lat = xr.DataArray(range(20, 71), dims="lat", name="lat")
    data = xr.DataArray(
        [[0.0] * lon.size for _ in range(lat.size)],
        coords={"lat": lat, "lon": lon},
        dims=("lat", "lon"),
        name="sample",
    )
    return data


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_mask_custom_rectangular_box1():
    da = _build_test_dataarray()
    mask = ecl.filter.mask_custom_rectangular_box(
        da, 280, 340, 38, 47, angle=20, center="center"
    )

    assert mask.dtype == bool
    assert mask.dims == ("lat", "lon")

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()})
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=["bottom", "left"], alpha=0.3)
    ax.set_extent([260, 350, 20, 70], crs=ccrs.PlateCarree())

    mask.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap="Greys",
        add_colorbar=False,
    )
    ecl.plot.quick_draw_custom_rectangular_box(
        280, 340, 38, 47, ax=ax, angle=20, ec="r", transform=ccrs.PlateCarree()
    )
    ax.set_title("")
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_mask_custom_rectangular_box2():
    da = _build_test_dataarray()
    mask = ecl.filter.mask_custom_rectangular_box(
        da,
        280,
        340,
        38,
        47,
        angle=20,
        center="lowerleft",
        use_coslat=False,
    )

    assert mask.sel(lat=38, lon=280).item() is True

    fig, ax = plt.subplots()
    mask.plot(ax=ax, cmap="Greys", add_colorbar=False)
    ecl.plot.quick_draw_custom_rectangular_box(
        280, 340, 38, 47, ax=ax, angle=20, center="lowerleft", ec="r"
    )
    ax.set_xlim(260, 350)
    ax.set_ylim(20, 70)
    ax.set_title("")
    return fig
