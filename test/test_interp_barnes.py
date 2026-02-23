"""
pytest for interp.barnes.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


@pytest.fixture
def sample_station_data():
    return ecl.open_tutorial_dataset("PressQFF_202007271200_872.csv")


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_interp_spatial_barnes1(sample_station_data):
    result = ecl.interp.interp_spatial_barnes(
        sample_station_data,
        "qff",
        lon_dim="lon",
        lat_dim="lat",
        influence_radius_deg=None,  # Default = 4*sigma
        mask_radius_deg=None,  # Default = influence radius
        method="optimized_convolution",
    )

    fig, ax = plt.subplots(
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=0)}
    )

    ax.gridlines(
        draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--"
    )
    ax.coastlines(edgecolor="black", linewidths=0.5)
    ax.set_extent([-30, 50, 32, 75])

    # Draw interpolation results
    result.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cbar_kwargs={"location": "bottom", "aspect": 50, "shrink": 0.9},
        cmap="RdBu_r",
        levels=21,
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_interp_spatial_barnes2(sample_station_data):
    result = ecl.interp.interp_spatial_barnes(
        sample_station_data,
        "qff",
        lon_dim="lon",
        lat_dim="lat",
        grid_res_deg=0.25,
        sigma_deg=0.5,
        influence_radius_deg=None,  # Default = 4*sigma
        mask_radius_deg=None,  # Default = influence radius
        method="optimized_convolution",
        buffer_deg=5.0,
    )

    fig, ax = plt.subplots(
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=0)}
    )

    ax.gridlines(
        draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--"
    )
    ax.coastlines(edgecolor="black", linewidths=0.5)
    ax.set_extent([-30, 50, 32, 75])

    # Draw interpolation results
    result.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cbar_kwargs={"location": "bottom", "aspect": 50, "shrink": 0.9},
        cmap="RdBu_r",
        levels=21,
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_interp_spatial_barnes3(sample_station_data):
    result = ecl.interp.interp_spatial_barnes(
        sample_station_data,
        "qff",
        lon_dim="lon",
        lat_dim="lat",
        grid_res_deg=0.25,
        sigma_deg=0.5,
        influence_radius_deg=None,  # Default = 4*sigma
        mask_radius_deg=None,  # Default = influence radius
        method="convolution",
        buffer_deg=5.0,
    )

    fig, ax = plt.subplots(
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=0)}
    )

    ax.gridlines(
        draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--"
    )
    ax.coastlines(edgecolor="black", linewidths=0.5)
    ax.set_extent([-30, 50, 32, 75])

    # Draw interpolation results
    result.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cbar_kwargs={"location": "bottom", "aspect": 50, "shrink": 0.9},
        cmap="RdBu_r",
        levels=21,
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_interp_spatial_barnes4(sample_station_data):
    result = ecl.interp.interp_spatial_barnes(
        sample_station_data,
        "qff",
        lon_dim="lon",
        lat_dim="lat",
        grid_res_deg=0.25,
        sigma_deg=0.5,
        influence_radius_deg=None,  # Default = 4*sigma
        mask_radius_deg=None,  # Default = influence radius
        method="radius",
        buffer_deg=5.0,
    )

    fig, ax = plt.subplots(
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=0)}
    )

    ax.gridlines(
        draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--"
    )
    ax.coastlines(edgecolor="black", linewidths=0.5)
    ax.set_extent([-30, 50, 32, 75])

    # Draw interpolation results
    result.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cbar_kwargs={"location": "bottom", "aspect": 50, "shrink": 0.9},
        cmap="RdBu_r",
        levels=21,
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_interp_spatial_barnes_rs1(sample_station_data):
    result = ecl.interp.interp_spatial_barnes_rs(
        sample_station_data,
        "qff",
        lon_dim="lon",
        lat_dim="lat",
        influence_radius_deg=None,  # Default = 4*sigma
        mask_radius_deg=None,  # Default = influence radius
        method="optimized_convolution",
    )

    fig, ax = plt.subplots(
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=0)}
    )

    ax.gridlines(
        draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--"
    )
    ax.coastlines(edgecolor="black", linewidths=0.5)
    ax.set_extent([-30, 50, 32, 75])

    # Draw interpolation results
    result.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cbar_kwargs={"location": "bottom", "aspect": 50, "shrink": 0.9},
        cmap="RdBu_r",
        levels=21,
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_interp_spatial_barnes_rs2(sample_station_data):
    result = ecl.interp.interp_spatial_barnes_rs(
        sample_station_data,
        "qff",
        lon_dim="lon",
        lat_dim="lat",
        grid_res_deg=0.25,
        sigma_deg=0.5,
        influence_radius_deg=None,  # Default = 4*sigma
        mask_radius_deg=None,  # Default = influence radius
        method="optimized_convolution",
        buffer_deg=5.0,
    )

    fig, ax = plt.subplots(
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=0)}
    )

    ax.gridlines(
        draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--"
    )
    ax.coastlines(edgecolor="black", linewidths=0.5)
    ax.set_extent([-30, 50, 32, 75])

    # Draw interpolation results
    result.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cbar_kwargs={"location": "bottom", "aspect": 50, "shrink": 0.9},
        cmap="RdBu_r",
        levels=21,
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_interp_spatial_barnes_rs3(sample_station_data):
    result = ecl.interp.interp_spatial_barnes_rs(
        sample_station_data,
        "qff",
        lon_dim="lon",
        lat_dim="lat",
        grid_res_deg=0.25,
        sigma_deg=0.5,
        influence_radius_deg=None,  # Default = 4*sigma
        mask_radius_deg=None,  # Default = influence radius
        method="convolution",
        buffer_deg=5.0,
    )

    fig, ax = plt.subplots(
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=0)}
    )

    ax.gridlines(
        draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--"
    )
    ax.coastlines(edgecolor="black", linewidths=0.5)
    ax.set_extent([-30, 50, 32, 75])

    # Draw interpolation results
    result.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cbar_kwargs={"location": "bottom", "aspect": 50, "shrink": 0.9},
        cmap="RdBu_r",
        levels=21,
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_interp_spatial_barnes_rs4(sample_station_data):
    result = ecl.interp.interp_spatial_barnes_rs(
        sample_station_data,
        "qff",
        lon_dim="lon",
        lat_dim="lat",
        grid_res_deg=0.25,
        sigma_deg=0.5,
        influence_radius_deg=None,  # Default = 4*sigma
        mask_radius_deg=None,  # Default = influence radius
        method="radius",
        buffer_deg=5.0,
    )

    fig, ax = plt.subplots(
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=0)}
    )

    ax.gridlines(
        draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--"
    )
    ax.coastlines(edgecolor="black", linewidths=0.5)
    ax.set_extent([-30, 50, 32, 75])

    # Draw interpolation results
    result.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cbar_kwargs={"location": "bottom", "aspect": 50, "shrink": 0.9},
        cmap="RdBu_r",
        levels=21,
    )
    return fig


# --------- S2 -------------------------------------------------


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_interp_spatial_barnesS21(sample_station_data):
    result = ecl.interp.interp_spatial_barnesS2(
        sample_station_data,
        "qff",
        lon_dim="lon",
        lat_dim="lat",
        influence_radius_deg=None,  # Default = 4*sigma
        mask_radius_deg=None,  # Default = influence radius
        method="optimized_convolution_S2",
    )

    fig, ax = plt.subplots(
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=0)}
    )

    ax.gridlines(
        draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--"
    )
    ax.coastlines(edgecolor="black", linewidths=0.5)
    ax.set_extent([-30, 50, 32, 75])

    # Draw interpolation results
    result.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cbar_kwargs={"location": "bottom", "aspect": 50, "shrink": 0.9},
        cmap="RdBu_r",
        levels=21,
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_interp_spatial_barnesS22(sample_station_data):
    result = ecl.interp.interp_spatial_barnesS2(
        sample_station_data,
        "qff",
        lon_dim="lon",
        lat_dim="lat",
        grid_res_deg=0.25,
        sigma_deg=0.5,
        influence_radius_deg=None,  # Default = 4*sigma
        mask_radius_deg=None,  # Default = influence radius
        method="optimized_convolution_S2",
        buffer_deg=5.0,
    )

    fig, ax = plt.subplots(
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=0)}
    )

    ax.gridlines(
        draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--"
    )
    ax.coastlines(edgecolor="black", linewidths=0.5)
    ax.set_extent([-30, 50, 32, 75])

    # Draw interpolation results
    result.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cbar_kwargs={"location": "bottom", "aspect": 50, "shrink": 0.9},
        cmap="RdBu_r",
        levels=21,
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_interp_spatial_barnesS23(sample_station_data):
    result = ecl.interp.interp_spatial_barnesS2(
        sample_station_data,
        "qff",
        lon_dim="lon",
        lat_dim="lat",
        grid_res_deg=0.25,
        sigma_deg=0.5,
        influence_radius_deg=None,  # Default = 4*sigma
        mask_radius_deg=None,  # Default = influence radius
        method="naive_S2",
        buffer_deg=5.0,
    )

    fig, ax = plt.subplots(
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=0)}
    )

    ax.gridlines(
        draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--"
    )
    ax.coastlines(edgecolor="black", linewidths=0.5)
    ax.set_extent([-30, 50, 32, 75])

    # Draw interpolation results
    result.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cbar_kwargs={"location": "bottom", "aspect": 50, "shrink": 0.9},
        cmap="RdBu_r",
        levels=21,
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_interp_spatial_barnesS2_rs1(sample_station_data):
    result = ecl.interp.interp_spatial_barnesS2_rs(
        sample_station_data,
        "qff",
        lon_dim="lon",
        lat_dim="lat",
        influence_radius_deg=None,  # Default = 4*sigma
        mask_radius_deg=None,  # Default = influence radius
        method="optimized_convolution_S2",
    )

    fig, ax = plt.subplots(
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=0)}
    )

    ax.gridlines(
        draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--"
    )
    ax.coastlines(edgecolor="black", linewidths=0.5)
    ax.set_extent([-30, 50, 32, 75])

    # Draw interpolation results
    result.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cbar_kwargs={"location": "bottom", "aspect": 50, "shrink": 0.9},
        cmap="RdBu_r",
        levels=21,
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_interp_spatial_barnesS2_rs2(sample_station_data):
    result = ecl.interp.interp_spatial_barnesS2_rs(
        sample_station_data,
        "qff",
        lon_dim="lon",
        lat_dim="lat",
        grid_res_deg=0.25,
        sigma_deg=0.5,
        influence_radius_deg=None,  # Default = 4*sigma
        mask_radius_deg=None,  # Default = influence radius
        method="optimized_convolution_S2",
        buffer_deg=5.0,
    )

    fig, ax = plt.subplots(
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=0)}
    )

    ax.gridlines(
        draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--"
    )
    ax.coastlines(edgecolor="black", linewidths=0.5)
    ax.set_extent([-30, 50, 32, 75])

    # Draw interpolation results
    result.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cbar_kwargs={"location": "bottom", "aspect": 50, "shrink": 0.9},
        cmap="RdBu_r",
        levels=21,
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_interp_spatial_barnesS2_rs3(sample_station_data):
    result = ecl.interp.interp_spatial_barnesS2_rs(
        sample_station_data,
        "qff",
        lon_dim="lon",
        lat_dim="lat",
        grid_res_deg=0.25,
        sigma_deg=0.5,
        influence_radius_deg=None,  # Default = 4*sigma
        mask_radius_deg=None,  # Default = influence radius
        method="naive_S2",
        buffer_deg=5.0,
    )

    fig, ax = plt.subplots(
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=0)}
    )

    ax.gridlines(
        draw_labels=["bottom", "left"], color="grey", alpha=0.5, linestyle="--"
    )
    ax.coastlines(edgecolor="black", linewidths=0.5)
    ax.set_extent([-30, 50, 32, 75])

    # Draw interpolation results
    result.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cbar_kwargs={"location": "bottom", "aspect": 50, "shrink": 0.9},
        cmap="RdBu_r",
        levels=21,
    )
    return fig
