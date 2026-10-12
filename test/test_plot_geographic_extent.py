"""Regression checks for the geographic windows used by the gallery."""

import matplotlib
import numpy as np
import pytest
import importlib
import sys
import types
from pathlib import Path

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr


# Keep these plotting-only checks runnable without the optional compiled
# analysis backends, including in environments with an older Cartopy.
_root = Path(__file__).parents[1] / "src" / "easyclimate" / "plot"
for _suffix in ["", ".icon", ".mpas"]:
    _name = "_test_extent_plot" + _suffix
    _package = types.ModuleType(_name)
    _package.__path__ = [str(_root / _suffix.lstrip("."))]
    sys.modules[_name] = _package

_set_cartopy_extent = importlib.import_module(
    "_test_extent_plot.icon.common"
)._set_cartopy_extent
_decorate_vertex_contour_axis = importlib.import_module(
    "_test_extent_plot.mpas.common"
)._decorate_vertex_contour_axis
draw_polar_basemap = importlib.import_module(
    "_test_extent_plot.polar"
).draw_polar_basemap
plot_triangular_grid = importlib.import_module(
    "_test_extent_plot.icon.triangular_grid"
).plot_triangular_grid
set_map_extent = importlib.import_module("_test_extent_plot._geo")._set_cartopy_extent
set_geoextent = importlib.import_module("_test_extent_plot.geoextent").set_geoextent


@pytest.mark.parametrize("projection_type", [ccrs.PlateCarree, ccrs.Mercator])
@pytest.mark.parametrize("central_longitude", [0, 140])
def test_full_longitude_window_does_not_collapse(projection_type, central_longitude):
    projection = projection_type(central_longitude=central_longitude)
    fig, ax = plt.subplots(subplot_kw={"projection": projection})
    try:
        assert set_geoextent([0, 360, 10, 85], ax=ax) is ax
        np.testing.assert_allclose(ax.get_xlim(), projection.x_limits)
        expected_y = [
            projection.transform_point(central_longitude, lat, ccrs.PlateCarree())[1]
            for lat in [10, 85]
        ]
        np.testing.assert_allclose(ax.get_ylim(), expected_y)
    finally:
        plt.close(fig)


@pytest.mark.parametrize("extents", [[170, -170, 10, 60], [170, 190, 10, 60]])
def test_geoextent_dateline_window(extents):
    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree(180)})
    try:
        set_geoextent(extents, ax=ax)
        np.testing.assert_allclose(ax.get_xlim(), [-10, 10], atol=1e-6)
        np.testing.assert_allclose(ax.get_ylim(), [10, 60], atol=1e-6)
    finally:
        plt.close(fig)


def test_geoextent_accepts_shifted_source_crs_and_current_axes():
    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree(-120)})
    try:
        assert set_geoextent([-130, -20, 25, 70], crs=ccrs.PlateCarree(-120)) is ax
        np.testing.assert_allclose(ax.get_xlim(), [-130, -20], atol=1e-6)
    finally:
        plt.close(fig)


def test_geoextent_window_crossing_target_seam_includes_both_edges():
    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()})
    try:
        set_geoextent([170, -170, 10, 60], ax=ax)
        np.testing.assert_allclose(ax.get_xlim(), ax.projection.x_limits)
        np.testing.assert_allclose(ax.get_ylim(), [10, 60])
    finally:
        plt.close(fig)


@pytest.mark.parametrize(
    "projection",
    [
        ccrs.NorthPolarStereo(140),
        ccrs.SouthPolarStereo(140),
        ccrs.Robinson(140),
    ],
)
def test_geoextent_full_longitude_non_cylindrical(projection):
    fig, ax = plt.subplots(subplot_kw={"projection": projection})
    try:
        south = isinstance(projection, ccrs.SouthPolarStereo)
        set_geoextent([0, 360, -90, -10] if south else [0, 360, 10, 90], ax=ax)
        assert ax.get_xlim()[0] < 0 < ax.get_xlim()[1]
        if isinstance(projection, (ccrs.NorthPolarStereo, ccrs.SouthPolarStereo)):
            np.testing.assert_allclose(ax.get_xlim(), ax.get_ylim())
            np.testing.assert_allclose(ax.get_xlim()[0], -ax.get_xlim()[1])
    finally:
        plt.close(fig)


@pytest.mark.parametrize("projection", [ccrs.PlateCarree(140), ccrs.Mercator(140)])
def test_geoextent_global_including_poles(projection):
    fig, ax = plt.subplots(subplot_kw={"projection": projection})
    try:
        set_geoextent([-180, 180, -90, 90], ax=ax)
        np.testing.assert_allclose(ax.get_xlim(), projection.x_limits)
        np.testing.assert_allclose(ax.get_ylim(), projection.y_limits)
    finally:
        plt.close(fig)


@pytest.mark.parametrize(
    "extents",
    [
        [0, 0, 10, 85],
        [0, 720, 10, 85],
        [0, 360, 85, 10],
        [0, 360, 10, 91],
        [0, np.nan, 10, 85],
        [0, 360, 10],
    ],
)
def test_geoextent_rejects_invalid_geographic_window(extents):
    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()})
    try:
        with pytest.raises(ValueError):
            set_geoextent(extents, ax=ax)
    finally:
        plt.close(fig)


def test_geoextent_requires_map_axes():
    fig, ax = plt.subplots()
    try:
        with pytest.raises(TypeError, match="GeoAxes"):
            set_geoextent([0, 360, 10, 85], ax=ax)
    finally:
        plt.close(fig)


def test_geoextent_delegates_native_projection_coordinates():
    projection = ccrs.Mercator(140)
    fig, ax = plt.subplots(subplot_kw={"projection": projection})
    try:
        set_geoextent([-1e6, 1e6, 1e6, 2e6], ax=ax, crs=projection)
        np.testing.assert_allclose(ax.get_xlim(), [-1e6, 1e6])
        np.testing.assert_allclose(ax.get_ylim(), [1e6, 2e6])
    finally:
        plt.close(fig)


@pytest.mark.parametrize("backend", ["icon", "mpas"])
def test_shifted_plate_carree_keeps_regional_window(backend):
    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree(-120)})
    try:
        if backend == "icon":
            _set_cartopy_extent(ax, 110, 220, 25, 70)
        else:
            _decorate_vertex_contour_axis(
                ax,
                cartopy_axis=True,
                lon_min_plot=110,
                lon_max_plot=220,
                lat_min_plot=25,
                lat_max_plot=70,
            )
        np.testing.assert_allclose(ax.get_xlim(), [-130, -20], atol=1e-5)
        np.testing.assert_allclose(ax.get_ylim(), [25, 70], atol=1e-5)
    finally:
        plt.close(fig)


def test_icon_projected_grid_splits_at_target_longitude_seam():
    # This triangle crosses the seam of a map centered on 180 degrees.
    ds = xr.Dataset(
        {
            "clon": ("ncells", [0.0]),
            "clat": ("ncells", [30.0]),
            "clon_bnds": (("ncells", "vertices"), [[-1.0, 1.0, 0.0]]),
            "clat_bnds": (("ncells", "vertices"), [[29.0, 29.0, 31.0]]),
        }
    )
    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree(180)})
    try:
        lines = plot_triangular_grid(
            ds,
            ax=ax,
            transform=ccrs.PlateCarree(),
            input_radians=False,
            global_plot=True,
        )
        segments = np.asarray(lines.get_segments())
        assert np.max(np.abs(np.diff(segments[:, :, 0], axis=1))) < 5
    finally:
        plt.close(fig)


@pytest.mark.parametrize("south", [False, True])
@pytest.mark.parametrize("central_longitude", [0, 140])
def test_polar_basemap_contains_complete_latitude_circle(south, central_longitude):
    projection = (ccrs.SouthPolarStereo if south else ccrs.NorthPolarStereo)(
        central_longitude=central_longitude
    )
    fig, ax = plt.subplots(subplot_kw={"projection": projection})
    try:
        draw_polar_basemap(
            ax=ax,
            lat_range=[-90, -10] if south else [10, 90],
            lon_step=30,
            lat_step=20,
        )
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        np.testing.assert_allclose([x0, y0], [-x1, -y1], atol=1e-5)
        np.testing.assert_allclose(x1, y1, rtol=1e-5)
        lons = np.linspace(-180, 180, 361)
        points = projection.transform_points(
            ccrs.PlateCarree(), lons, np.full_like(lons, -10 if south else 10)
        )
        assert np.max(np.abs(points[:, :2])) <= x1
    finally:
        plt.close(fig)
