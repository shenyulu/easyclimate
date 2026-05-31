"""
pytest for plot/mpas
"""

import pytest

import matplotlib
import numpy as np

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import easyclimate as ecl

from easyclimate.plot.mpas import (
    plot_cell_barbs,
    plot_cell_contour,
    plot_cell_contourf,
    plot_cell_curved_quiver,
    plot_cell_quiver,
    plot_cell_streamplot,
    plot_cell_voronoi,
    plot_vertex_contour,
    plot_vertex_contourf,
    plot_vertex_voronoi,
    plot_voronoi_grid,
)
from easyclimate.plot.mpas.voronoi_extract import extract_cell_latlon


def _open_tutorial_dataset_or_skip(name):
    try:
        return ecl.open_tutorial_dataset(name)
    except Exception as exc:
        pytest.skip(f"tutorial dataset {name!r} is unavailable: {exc}")


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_plot_mpas_voronoi_grid():
    data = _open_tutorial_dataset_or_skip("x1_2562_grid")

    fig = plt.figure(figsize=(8, 4))
    ax = plt.axes(projection=ccrs.PlateCarree(180))
    ax.coastlines(resolution="110m", linewidth=0.6)

    plot_voronoi_grid(data, ax=ax, transform=ccrs.PlateCarree())
    ax.set_global()

    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_plot_mpas_cell_voronoi():
    data = _open_tutorial_dataset_or_skip("mpas_JWwave_T10_nVertLevels10")

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree(160)})

    plot_cell_voronoi(
        data,
        data.divergence,
        ax=ax,
        transform=ccrs.PlateCarree(),
        lon_min=120,
        lon_max=200,
        lat_min=20,
        lat_max=70,
        vmax=1e-6,
        linewidth=0.1,
        edgecolor="grey",
        cbar_kwargs={"location": "bottom", "aspect": 60},
    )

    ax.coastlines(resolution="110m", linewidth=0.6)
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_plot_mpas_cell_contourf_and_contour():
    data = _open_tutorial_dataset_or_skip("mpas_JWwave_T10_nVertLevels10")

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree(160)})

    plot_cell_contourf(
        data,
        data.divergence,
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap="RdBu_r",
        levels=np.linspace(-1e-6, 1e-6, 11),
        lon_min=120,
        lon_max=200,
        lat_min=20,
        lat_max=70,
        cbar_kwargs={"location": "bottom", "aspect": 60},
    )

    plot_cell_contour(
        data,
        data.divergence,
        ax=ax,
        transform=ccrs.PlateCarree(),
        levels=np.linspace(-1e-6, 1e-6, 5),
        colors="black",
        lon_min=120,
        lon_max=200,
        lat_min=20,
        lat_max=70,
    )

    ax.coastlines(resolution="110m", linewidth=0.6)
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_plot_mpas_vertex_contourf_and_contour():
    data = _open_tutorial_dataset_or_skip("mpas_JWwave_T10_nVertLevels10")

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree(-120)})

    plot_vertex_contourf(
        data,
        data.vorticity,
        ax=ax,
        transform=ccrs.PlateCarree(),
        levels=np.linspace(-4e-5, 4e-5, 11),
        lon_min=130,
        lon_max=320,
        lat_min=10,
        lat_max=80,
        cbar_kwargs={"location": "bottom", "aspect": 60},
    )

    plot_vertex_contour(
        data,
        data.vorticity,
        ax=ax,
        transform=ccrs.PlateCarree(),
        levels=np.linspace(-4e-5, 4e-5, 5),
        colors="black",
        lon_min=130,
        lon_max=320,
        lat_min=10,
        lat_max=80,
    )

    ax.coastlines(resolution="110m", linewidth=0.6)
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_plot_mpas_cell_quiver_and_barbs():
    data = _open_tutorial_dataset_or_skip("mpas_JWwave_T10_nVertLevels10")

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree(-120)})

    plot_cell_quiver(
        data,
        data.uReconstructZonal,
        data.uReconstructMeridional,
        ax=ax,
        lon_min=190,
        lon_max=-60,
        lat_min=10,
        lat_max=80,
        nx_bins=18,
        ny_bins=12,
        thin_method="nearest_center",
        scale=900,
        width=0.004,
    )

    plot_cell_barbs(
        data,
        data.uReconstructZonal,
        data.uReconstructMeridional,
        ax=ax,
        lon_min=190,
        lon_max=-60,
        lat_min=10,
        lat_max=80,
        nx_bins=12,
        ny_bins=8,
        thin_method="mean",
        length=4,
        linewidth=0.5,
        alpha=0.8,
    )

    ax.coastlines(resolution="110m", linewidth=0.6)
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_plot_mpas_vertex_voronoi_from_vertex_values():
    data = _open_tutorial_dataset_or_skip("mpas_JWwave_T10_nVertLevels10")

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree(160)})

    plot_vertex_voronoi(
        data,
        data.vorticity,
        ax=ax,
        transform=ccrs.PlateCarree(),
        lon_min=200,
        lon_max=270,
        lat_min=30,
        lat_max=70,
        linewidth=0.1,
        edgecolor="grey",
        cbar_kwargs={"location": "bottom", "aspect": 60},
    )

    ax.coastlines(resolution="110m", linewidth=0.6)
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_plot_mpas_vertex_voronoi_from_cell_values_plain_axes():
    data = _open_tutorial_dataset_or_skip("mpas_JWwave_T10_nVertLevels10")

    fig, ax = plt.subplots(figsize=(6, 4))

    plot_vertex_voronoi(
        data,
        data.divergence,
        ax=ax,
        lon_min=120,
        lon_max=200,
        lat_min=20,
        lat_max=70,
        symmetric=False,
        add_colorbar=False,
        linewidth=0.1,
        edgecolor="grey",
    )

    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_plot_mpas_cell_streamplot_plain_axes():
    data = _open_tutorial_dataset_or_skip("mpas_JWwave_T10_nVertLevels10")

    fig, ax = plt.subplots(figsize=(6, 4))

    plot_cell_streamplot(
        data,
        data.uReconstructZonal,
        data.uReconstructMeridional,
        ax=ax,
        lon_min=190,
        lon_max=300,
        lat_min=10,
        lat_max=80,
        nx=24,
        ny=16,
        interpolation_padding=8,
        density=0.6,
        linewidth=0.8,
    )

    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_plot_mpas_cell_curved_quiver_plain_axes():
    data = _open_tutorial_dataset_or_skip("mpas_JWwave_T10_nVertLevels10")

    fig, ax = plt.subplots(figsize=(6, 4))

    plot_cell_curved_quiver(
        data,
        data.uReconstructZonal,
        data.uReconstructMeridional,
        ax=ax,
        lon_min=190,
        lon_max=300,
        lat_min=10,
        lat_max=80,
        nx=18,
        ny=12,
        interpolation_padding=(8, 8),
        density=0.7,
        arrowsize=0.8,
        grains=8,
        mask_density=6,
    )

    return fig


def test_plot_mpas_extract_cell_latlon_global_and_no_auto_extent():
    data = _open_tutorial_dataset_or_skip("mpas_JWwave_T10_nVertLevels10")

    global_geo = extract_cell_latlon(
        data,
        da=data.divergence,
        global_plot=True,
        center_lon=180,
    )
    assert global_geo["lon_min_plot"] == 0
    assert global_geo["lon_max_plot"] == 360
    assert global_geo["lat_min_plot"] == -90
    assert global_geo["lat_max_plot"] == 90

    no_extent_geo = extract_cell_latlon(data, auto_extent=False, center_lon=0)
    assert no_extent_geo["lon_min_plot"] is None
    assert no_extent_geo["lon_max_plot"] is None
    assert no_extent_geo["lat_min_plot"] is None
    assert no_extent_geo["lat_max_plot"] is None


def test_plot_mpas_invalid_scalar_shapes_raise():
    data = _open_tutorial_dataset_or_skip("mpas_JWwave_T10_nVertLevels10")
    bad_values = np.zeros((2, 2))

    with pytest.raises(ValueError, match="must be 1D"):
        plot_cell_voronoi(data, bad_values)

    with pytest.raises(ValueError, match="must be 1D"):
        plot_cell_contourf(data, bad_values)

    with pytest.raises(ValueError, match="must be 1D"):
        plot_vertex_voronoi(data, bad_values)

    with pytest.raises(ValueError, match="does not match nVertices or nCells"):
        plot_vertex_voronoi(data, np.zeros(3))


def test_plot_mpas_vector_extent_and_shape_errors():
    data = _open_tutorial_dataset_or_skip("mpas_JWwave_T10_nVertLevels10")
    u = data.uReconstructZonal
    v = data.uReconstructMeridional

    with pytest.raises(ValueError, match="Please provide lon_min"):
        plot_cell_quiver(data, u, v)

    with pytest.raises(ValueError, match="Please provide lon_min"):
        plot_cell_barbs(data, u, v)

    with pytest.raises(ValueError, match="Streamplot needs"):
        plot_cell_streamplot(data, u, v)

    with pytest.raises(ValueError, match="Curved quiver needs"):
        plot_cell_curved_quiver(data, u, v)

    with pytest.raises(ValueError, match="shape mismatch"):
        plot_cell_quiver(
            data,
            np.asarray(u.values)[:-1],
            v,
            lon_min=190,
            lon_max=300,
            lat_min=10,
            lat_max=80,
        )


def test_plot_mpas_regrid_shape_requires_cartopy_axes():
    data = _open_tutorial_dataset_or_skip("mpas_JWwave_T10_nVertLevels10")
    u = data.uReconstructZonal
    v = data.uReconstructMeridional
    fig, ax = plt.subplots()

    common_kwargs = dict(
        ax=ax,
        lon_min=190,
        lon_max=300,
        lat_min=10,
        lat_max=80,
        regrid_shape=8,
    )

    with pytest.raises(ValueError, match="cartopy GeoAxes"):
        plot_cell_quiver(data, u, v, **common_kwargs)

    with pytest.raises(ValueError, match="cartopy GeoAxes"):
        plot_cell_barbs(data, u, v, **common_kwargs)

    with pytest.raises(ValueError, match="cartopy GeoAxes"):
        plot_cell_streamplot(data, u, v, **common_kwargs)

    with pytest.raises(ValueError, match="cartopy GeoAxes"):
        plot_cell_curved_quiver(data, u, v, **common_kwargs)

    plt.close(fig)


def test_plot_mpas_cartopy_regrid_shape_branches():
    data = _open_tutorial_dataset_or_skip("mpas_JWwave_T10_nVertLevels10")
    u = data.uReconstructZonal
    v = data.uReconstructMeridional
    kwargs = dict(
        lon_min=190,
        lon_max=300,
        lat_min=10,
        lat_max=80,
        transform=ccrs.PlateCarree(),
        regrid_shape=8,
    )

    fig, axs = plt.subplots(
        2,
        2,
        subplot_kw={"projection": ccrs.PlateCarree(-120)},
        figsize=(8, 6),
    )

    q = plot_cell_quiver(data, u, v, ax=axs[0, 0], add_quiverkey=False, **kwargs)
    b = plot_cell_barbs(data, u, v, ax=axs[0, 1], **kwargs)
    sp = plot_cell_streamplot(
        data,
        u,
        v,
        ax=axs[1, 0],
        density=0.5,
        **kwargs,
    )
    cq = plot_cell_curved_quiver(
        data,
        u,
        v,
        ax=axs[1, 1],
        density=0.5,
        grains=6,
        mask_density=5,
        **kwargs,
    )

    assert q is not None
    assert b is not None
    assert sp is not None
    assert cq is not None
    plt.close(fig)


def test_plot_mpas_cartopy_project_to_map_branches():
    data = _open_tutorial_dataset_or_skip("mpas_JWwave_T10_nVertLevels10")
    u = data.uReconstructZonal
    v = data.uReconstructMeridional
    kwargs = dict(
        lon_min=190,
        lon_max=300,
        lat_min=10,
        lat_max=80,
        transform=ccrs.PlateCarree(),
        nx=18,
        ny=12,
        interpolation_padding=8,
    )

    fig, axs = plt.subplots(
        2,
        2,
        subplot_kw={"projection": ccrs.PlateCarree(-120)},
        figsize=(8, 6),
    )

    sp_projected = plot_cell_streamplot(
        data,
        u,
        v,
        ax=axs[0, 0],
        density=0.5,
        project_to_map=True,
        **kwargs,
    )
    sp_transform = plot_cell_streamplot(
        data,
        u,
        v,
        ax=axs[0, 1],
        density=0.5,
        project_to_map=False,
        **kwargs,
    )
    cq_projected = plot_cell_curved_quiver(
        data,
        u,
        v,
        ax=axs[1, 0],
        density=0.5,
        grains=6,
        mask_density=5,
        project_to_map=True,
        **kwargs,
    )
    cq_transform = plot_cell_curved_quiver(
        data,
        u,
        v,
        ax=axs[1, 1],
        density=0.5,
        grains=6,
        mask_density=5,
        project_to_map=False,
        **kwargs,
    )

    assert sp_projected is not None
    assert sp_transform is not None
    assert cq_projected is not None
    assert cq_transform is not None
    plt.close(fig)


def test_plot_mpas_global_and_plain_axes_branches():
    data = _open_tutorial_dataset_or_skip("mpas_JWwave_T10_nVertLevels10")

    fig, axs = plt.subplots(2, 2, figsize=(8, 6))

    ax_cell, pc_cell = plot_cell_voronoi(
        data,
        data.divergence,
        ax=axs[0, 0],
        global_plot=True,
        center_lon=180,
        add_colorbar=False,
        title="",
    )
    cs_cell = plot_cell_contourf(
        data,
        data.divergence,
        ax=axs[0, 1],
        global_plot=True,
        center_lon=180,
        add_colorbar=False,
        levels=5,
        title="",
    )
    cs_vertex = plot_vertex_contourf(
        data,
        data.vorticity,
        ax=axs[1, 0],
        global_plot=True,
        center_lon=180,
        add_colorbar=False,
        levels=5,
        title="",
    )
    pc_vertex = plot_vertex_voronoi(
        data,
        data.vorticity,
        ax=axs[1, 1],
        global_plot=True,
        center_lon=180,
        add_colorbar=False,
        title="",
    )

    assert ax_cell is axs[0, 0]
    assert pc_cell.get_array().size > 0
    assert cs_cell is not None
    assert cs_vertex is not None
    assert pc_vertex.get_array().size > 0
    plt.close(fig)


def test_plot_mpas_voronoi_grid_plain_subset_and_default_axes():
    data = _open_tutorial_dataset_or_skip("x1_2562_grid")

    fig, ax = plt.subplots()
    lc_subset = plot_voronoi_grid(
        data,
        ax=ax,
        lon_min=100,
        lon_max=180,
        lat_min=0,
        lat_max=80,
        colors="black",
    )
    assert len(lc_subset.get_segments()) > 0
    plt.close(fig)

    fig = plt.figure()
    lc_default_ax = plot_voronoi_grid(data)
    assert len(lc_default_ax.get_segments()) > 0
    plt.close(fig)


def test_plot_mpas_no_selected_data_errors():
    data = _open_tutorial_dataset_or_skip("mpas_JWwave_T10_nVertLevels10")
    u = data.uReconstructZonal
    v = data.uReconstructMeridional

    with pytest.raises(RuntimeError, match="No cells selected"):
        plot_cell_voronoi(
            data,
            data.divergence,
            lon_min=0,
            lon_max=1,
            lat_min=-90,
            lat_max=-89,
        )

    with pytest.raises(RuntimeError, match="No valid cell vectors selected"):
        plot_cell_quiver(
            data,
            u,
            v,
            lon_min=0,
            lon_max=1,
            lat_min=-90,
            lat_max=-89,
        )

    with pytest.raises(RuntimeError, match="At least 3 valid cell vectors"):
        plot_cell_streamplot(
            data,
            u,
            v,
            lon_min=0,
            lon_max=1,
            lat_min=-90,
            lat_max=-89,
        )
