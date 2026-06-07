"""
Tests for ICON plotting helpers.
"""

import numpy as np

import importlib.util
import matplotlib
import pytest
import sys
import types
from pathlib import Path

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt


def _load_icon_module(module_name):
    root = Path(__file__).parents[1] / "src" / "easyclimate" / "plot"
    plot_package_name = "_test_easyclimate_plot"
    icon_package_name = f"{plot_package_name}.icon"
    icon_package_path = root / "icon"

    if plot_package_name not in sys.modules:
        plot_package = types.ModuleType(plot_package_name)
        plot_package.__path__ = [str(root)]
        sys.modules[plot_package_name] = plot_package

    if icon_package_name not in sys.modules:
        icon_package = types.ModuleType(icon_package_name)
        icon_package.__path__ = [str(icon_package_path)]
        sys.modules[icon_package_name] = icon_package

    full_name = f"{icon_package_name}.{module_name}"
    if full_name in sys.modules:
        return sys.modules[full_name]

    spec = importlib.util.spec_from_file_location(
        full_name,
        icon_package_path / f"{module_name}.py",
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[full_name] = module
    spec.loader.exec_module(module)
    return module


cell_streamplot = _load_icon_module("cell_streamplot")
cell_curved_quiver = _load_icon_module("cell_curved_quiver")
cell_quiver = _load_icon_module("cell_quiver")
cell_barbs = _load_icon_module("cell_barbs")
cell_contour = _load_icon_module("cell_contour")
cell_triangular = _load_icon_module("cell_triangular")
triangular_grid = _load_icon_module("triangular_grid")


class _Values:
    def __init__(self, values, *, name=None, attrs=None):
        self.values = np.asarray(values)
        self.name = name
        self.attrs = attrs or {}


def _vector_dataset():
    return {
        "clon": _Values([100.0, 110.0, 120.0, 130.0, 140.0]),
        "clat": _Values([10.0, 20.0, 30.0, 40.0, 50.0]),
    }


def _triangular_dataset():
    return {
        "clon": _Values([1.0 / 3.0, 2.0 / 3.0, 1.5]),
        "clat": _Values([1.0 / 3.0, 2.0 / 3.0, 0.4]),
        "clon_bnds": _Values(
            [
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
                [1.0, 2.0, 1.5],
            ]
        ),
        "clat_bnds": _Values(
            [
                [0.0, 0.0, 1.0],
                [0.0, 1.0, 1.0],
                [0.0, 0.0, 1.0],
            ]
        ),
    }


def _scalar_values():
    return _Values(
        [2.0, 4.0, 1.0],
        name="temperature",
        attrs={"units": "K", "long_name": "air temperature"},
    )


def test_icon_cell_quiver_thins_vectors_and_adds_key(monkeypatch):
    recorded = {}

    def fake_thin_vectors_by_bins(x, y, u, v, **kwargs):
        recorded.update(kwargs)
        return (
            np.asarray([105.0, 135.0]),
            np.asarray([15.0, 45.0]),
            np.asarray([1.0, 3.0]),
            np.asarray([0.0, 4.0]),
            np.asarray([1.0, 5.0]),
        )

    monkeypatch.setattr(cell_quiver, "thin_vectors_by_bins", fake_thin_vectors_by_bins)

    fig, ax = plt.subplots()
    result = cell_quiver.plot_cell_quiver(
        _vector_dataset(),
        _Values([1.0, 2.0, 3.0, 4.0, 5.0]),
        _Values([0.0, 1.0, 0.0, 1.0, 0.0]),
        ax=ax,
        lon_min=100.0,
        lon_max=140.0,
        lat_min=10.0,
        lat_max=50.0,
        input_radians=False,
        nx_bins=8,
        ny_bins=4,
        thin_method="max_speed",
        quiverkey_value=5.0,
        quiverkey_label="5 m/s",
    )

    assert result in ax.collections
    assert recorded["nx_bins"] == 8
    assert recorded["ny_bins"] == 4
    assert recorded["method"] == "max_speed"
    assert recorded["use_cache"] is True
    assert np.allclose(ax.get_xlim(), (100.0, 140.0))
    assert np.allclose(ax.get_ylim(), (10.0, 50.0))
    assert ax.get_title() == "Cell-centered wind"
    plt.close(fig)


def test_icon_cell_barbs_thins_vectors_without_cache(monkeypatch):
    recorded = {}

    def fake_thin_vectors_by_bins(x, y, u, v, **kwargs):
        recorded.update(kwargs)
        return (
            np.asarray([105.0, 135.0]),
            np.asarray([15.0, 45.0]),
            np.asarray([1.0, 3.0]),
            np.asarray([0.0, 4.0]),
            np.asarray([1.0, 5.0]),
        )

    monkeypatch.setattr(cell_barbs, "thin_vectors_by_bins", fake_thin_vectors_by_bins)

    fig, ax = plt.subplots()
    result = cell_barbs.plot_cell_barbs(
        _vector_dataset(),
        _Values([1.0, 2.0, 3.0, 4.0, 5.0]),
        _Values([0.0, 1.0, 0.0, 1.0, 0.0]),
        ax=ax,
        lon_min=100.0,
        lon_max=140.0,
        lat_min=10.0,
        lat_max=50.0,
        input_radians=False,
        nx_bins=6,
        ny_bins=3,
        thin_method="mean",
        use_thinning_cache=False,
    )

    assert result in ax.collections
    assert recorded["nx_bins"] == 6
    assert recorded["ny_bins"] == 3
    assert recorded["method"] == "mean"
    assert recorded["cache_key"] is None
    assert recorded["use_cache"] is False
    assert ax.get_title() == "Cell-centered wind barbs"
    plt.close(fig)


def test_icon_cell_contours_return_contour_sets():
    ds = _triangular_dataset()
    da = _scalar_values()

    fig, (ax1, ax2) = plt.subplots(ncols=2)
    contourf = cell_contour.plot_cell_contourf(
        ds,
        da,
        ax=ax1,
        lon_min=0.0,
        lon_max=2.0,
        lat_min=0.0,
        lat_max=1.0,
        input_radians=False,
        levels=[0.0, 2.0, 4.0, 6.0],
        add_colorbar=False,
        use_triangulation_cache=False,
    )
    contour = cell_contour.plot_cell_contour(
        ds,
        da,
        ax=ax2,
        lon_min=0.0,
        lon_max=2.0,
        lat_min=0.0,
        lat_max=1.0,
        input_radians=False,
        levels=[0.0, 2.0, 4.0, 6.0],
        use_triangulation_cache=False,
    )

    assert contourf.axes is ax1
    assert contour.axes is ax2
    assert np.allclose(ax1.get_xlim(), (0.0, 2.0))
    assert np.allclose(ax1.get_ylim(), (0.0, 1.0))
    assert ax1.get_title() == "air temperature"
    assert ax2.get_title() == "air temperature"
    plt.close(fig)


def test_icon_cell_triangular_returns_poly_collection():
    fig, ax = plt.subplots()
    collection = cell_triangular.plot_cell_triangular(
        _triangular_dataset(),
        _scalar_values(),
        ax=ax,
        lon_min=0.0,
        lon_max=2.0,
        lat_min=0.0,
        lat_max=1.0,
        input_radians=False,
        add_colorbar=False,
        edgecolor="black",
        linewidth=0.2,
    )

    assert collection in ax.collections
    assert collection.get_array().size == 3
    assert np.allclose(ax.get_xlim(), (0.0, 2.0))
    assert np.allclose(ax.get_ylim(), (0.0, 1.0))
    assert ax.get_title() == "air temperature"
    plt.close(fig)


def test_icon_triangular_grid_deduplicates_shared_edges():
    ds = _triangular_dataset()

    fig, (ax1, ax2) = plt.subplots(ncols=2)
    deduped = triangular_grid.plot_triangular_grid(
        ds,
        ax=ax1,
        lon_min=0.0,
        lon_max=1.0,
        lat_min=0.0,
        lat_max=1.0,
        input_radians=False,
        deduplicate_edges=True,
    )
    raw = triangular_grid.plot_triangular_grid(
        ds,
        ax=ax2,
        lon_min=0.0,
        lon_max=1.0,
        lat_min=0.0,
        lat_max=1.0,
        input_radians=False,
        deduplicate_edges=False,
    )

    assert deduped in ax1.collections
    assert raw in ax2.collections
    assert len(deduped.get_segments()) == 5
    assert len(raw.get_segments()) == 6
    assert np.allclose(ax1.get_xlim(), (0.0, 1.0))
    assert np.allclose(ax1.get_ylim(), (0.0, 1.0))
    plt.close(fig)


def test_icon_cell_streamplot_pads_internal_grid_not_visible_axes(monkeypatch):
    recorded = {}

    def fake_interpolate_vectors_to_grid(
        x,
        y,
        u,
        v,
        *,
        lon_min,
        lon_max,
        lat_min,
        lat_max,
        nx,
        ny,
        **kwargs,
    ):
        recorded.update(
            lon_min=lon_min,
            lon_max=lon_max,
            lat_min=lat_min,
            lat_max=lat_max,
        )
        lon_grid = np.linspace(lon_min, lon_max, int(nx))
        lat_grid = np.linspace(lat_min, lat_max, int(ny))
        xx, yy = np.meshgrid(lon_grid, lat_grid)
        u_grid = np.ma.masked_invalid(np.ones_like(xx))
        v_grid = np.ma.masked_invalid(np.zeros_like(xx))
        return lon_grid, lat_grid, xx, yy, u_grid, v_grid

    monkeypatch.setattr(
        cell_streamplot,
        "_interpolate_vectors_to_grid",
        fake_interpolate_vectors_to_grid,
    )

    ds = {
        "clon": _Values([180.0, 200.0, 240.0, 290.0, 310.0]),
        "clat": _Values([0.0, 20.0, 40.0, 60.0, 80.0]),
    }
    u = _Values([1.0, 1.0, 1.0, 1.0, 1.0])
    v = _Values([0.0, 0.0, 0.0, 0.0, 0.0])

    fig, ax = plt.subplots()
    cell_streamplot.plot_cell_streamplot(
        ds,
        u,
        v,
        ax=ax,
        lon_min=190.0,
        lon_max=-60.0,
        lat_min=10.0,
        lat_max=70.0,
        input_radians=False,
        nx=12,
        ny=8,
    )

    assert recorded["lon_min"] < 190.0
    assert recorded["lon_max"] > 300.0
    assert recorded["lat_min"] < 10.0
    assert recorded["lat_max"] > 70.0
    assert np.allclose(ax.get_xlim(), (190.0, 300.0))
    assert np.allclose(ax.get_ylim(), (10.0, 70.0))
    plt.close(fig)


def test_icon_cell_curved_quiver_uses_velovect_backend(monkeypatch):
    calls = {}

    def fake_interpolate_vectors_to_grid(
        x,
        y,
        u,
        v,
        *,
        lon_min,
        lon_max,
        lat_min,
        lat_max,
        nx,
        ny,
        **kwargs,
    ):
        lon_grid = np.linspace(lon_min, lon_max, int(nx))
        lat_grid = np.linspace(lat_min, lat_max, int(ny))
        xx, yy = np.meshgrid(lon_grid, lat_grid)
        u_grid = np.ma.masked_invalid(np.ones_like(xx))
        v_grid = np.ma.masked_invalid(np.zeros_like(xx))
        return lon_grid, lat_grid, xx, yy, u_grid, v_grid

    def fake_velovect(ax, x, y, u, v, **kwargs):
        calls.update(kwargs)
        calls["x"] = x
        calls["y"] = y
        return "curved-quiver-result"

    monkeypatch.setattr(
        cell_curved_quiver,
        "_interpolate_vectors_to_grid",
        fake_interpolate_vectors_to_grid,
    )
    monkeypatch.setattr(cell_curved_quiver, "velovect", fake_velovect)

    ds = {
        "clon": _Values([140.0, 160.0, 180.0, 200.0, 220.0]),
        "clat": _Values([25.0, 35.0, 45.0, 55.0, 70.0]),
    }
    u = _Values([1.0, 1.0, 1.0, 1.0, 1.0])
    v = _Values([0.0, 0.0, 0.0, 0.0, 0.0])

    fig, ax = plt.subplots()
    result = cell_curved_quiver.plot_cell_curved_quiver(
        ds,
        u,
        v,
        ax=ax,
        lon_min=140.0,
        lon_max=220.0,
        lat_min=25.0,
        lat_max=70.0,
        input_radians=False,
        nx=12,
        ny=8,
        grains=9,
        mask_density=6,
        arrow_position=0.7,
    )

    assert result == "curved-quiver-result"
    assert calls["grains"] == 9
    assert calls["mask_density"] == 6
    assert calls["arrow_position"] == 0.7
    assert calls["glyph_mode"] is False
    assert np.allclose(ax.get_xlim(), (140.0, 220.0))
    assert np.allclose(ax.get_ylim(), (25.0, 70.0))
    plt.close(fig)


def test_icon_vector_plots_require_explicit_extent():
    ds = _vector_dataset()
    u = _Values([1.0, 1.0, 1.0, 1.0, 1.0])
    v = _Values([0.0, 0.0, 0.0, 0.0, 0.0])

    with pytest.raises(ValueError, match="lon_min"):
        cell_quiver.plot_cell_quiver(ds, u, v, input_radians=False)

    with pytest.raises(ValueError, match="Streamplot needs an explicit"):
        cell_streamplot.plot_cell_streamplot(ds, u, v, input_radians=False)

    with pytest.raises(ValueError, match="Curved quiver needs an explicit"):
        cell_curved_quiver.plot_cell_curved_quiver(ds, u, v, input_radians=False)
