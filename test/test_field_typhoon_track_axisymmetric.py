"""
pytest for field/typhoon/track.py and field/typhoon/axisymmetric.py
"""

import pytest

import easyclimate as ecl
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pathlib import Path
from .const_define import DOCS_DATA_PATH

# Data
t_data = xr.open_dataset(str(Path(DOCS_DATA_PATH, "test_t_typhoon_201919.nc"))).t
slp_data = xr.open_dataset(
    str(Path(DOCS_DATA_PATH, "test_slp_typhoon_201919.nc"))
).msl.isel(time=0)
tp_data = (
    xr.open_dataset(str(Path(DOCS_DATA_PATH, "test_pr_typhoon_201919.nc"))).tp * 1000
)
tp_data.attrs["units"] = "mm/h"


class TestTrackCycloneCenterMslOnly:
    def test_track_cyclone_center_msl_only(self):
        df_track = ecl.field.typhoon.track_cyclone_center_msl_only(slp_data, (140, 20))
        result_data = df_track.values.flatten()
        refer_data = np.array([1.40197533e02, 1.97669414e01, 9.62637600e04])
        assert np.isclose(result_data, refer_data, atol=0.01).all()


class TestCycloneAxisymmetricAnalysis:
    @pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
    def test_temp_rotated_symmetric(self):
        temp_tc_polar_result = ecl.field.typhoon.cyclone_axisymmetric_analysis(
            t_data, (140.19753277422885, 19.76694143269665)
        )

        ds_sym = temp_tc_polar_result["rotated_symmetric"]
        r1 = ds_sym - ds_sym.isel(y=0)

        # Draw begin
        fig, ax = plt.subplots()
        ecl.plot.set_p_format_axis()

        r1.plot.contourf(ax=ax, levels=np.linspace(-8, 8, 21), yincrease=False)
        return fig

    @pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
    def test_temp_rotated_asymmetric(self):
        temp_tc_polar_result = ecl.field.typhoon.cyclone_axisymmetric_analysis(
            t_data, (140.19753277422885, 19.76694143269665)
        )

        ds_asym = temp_tc_polar_result["rotated_asymmetric"]
        ds_asym850 = ds_asym.sel(level=850)

        ds_asym850_latlon = ds_asym850.swap_dims({"y": "lat", "polar_lon": "lon"})
        ds_asym850_latlon_0360 = ecl.plot.add_lon_cyclic(ds_asym850_latlon, inter=2)

        # Draw begin
        fig, ax = plt.subplots(subplot_kw={"projection": ccrs.NorthPolarStereo()})

        ecl.plot.draw_Circlemap_PolarStereo(
            ax=ax,
            lon_step=30,
            lat_step=2.5,
            lat_range=[80, 90],
            draw_labels=False,
            set_map_boundary_kwargs={"north_pad": 0.3, "south_pad": 0.4},
            gridlines_kwargs={"color": "grey", "alpha": 0.8, "linestyle": "--"},
        )

        (ds_asym850_latlon_0360 + 1e-13).plot.contourf(
            x="lon", y="lat", transform=ccrs.PlateCarree(), cmap="RdBu_r", levels=21
        )
        return fig

    @pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
    def test_pr_rotated_symmetric(self):
        pr_tc_polar_result = ecl.field.typhoon.cyclone_axisymmetric_analysis(
            tp_data, (140.19753277422885, 19.76694143269665)
        )

        ds_sym = pr_tc_polar_result["rotated_symmetric"]
        r1 = ds_sym - ds_sym.isel(y=0)

        # Draw begin
        fig, ax = plt.subplots()
        r1.plot(ax=ax)
        return fig

    @pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
    def test_pr_rotated_asymmetric(self):
        pr_tc_polar_result = ecl.field.typhoon.cyclone_axisymmetric_analysis(
            tp_data, (140.19753277422885, 19.76694143269665)
        )
        ds_asym = pr_tc_polar_result["rotated_asymmetric"]

        ds_asym_latlon = ds_asym.swap_dims({"y": "lat", "polar_lon": "lon"})
        ds_asym_latlon_0360 = ecl.plot.add_lon_cyclic(ds_asym_latlon, inter=2)

        # Draw begin
        fig, ax = plt.subplots(subplot_kw={"projection": ccrs.NorthPolarStereo()})

        ecl.plot.draw_Circlemap_PolarStereo(
            ax=ax,
            lon_step=30,
            lat_step=2.5,
            lat_range=[80, 90],
            draw_labels=False,
            set_map_boundary_kwargs={"north_pad": 0.3, "south_pad": 0.4},
            gridlines_kwargs={"color": "grey", "alpha": 0.8, "linestyle": "--"},
        )

        (ds_asym_latlon_0360 + 1e-13).plot.contourf(
            x="lon", y="lat", transform=ccrs.PlateCarree(), cmap="RdBu_r", levels=21
        )
        return fig
