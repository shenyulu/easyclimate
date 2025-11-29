"""
pytest for field/equatorial_wave/mjo.py
"""

import pytest
import xarray as xr
import matplotlib.pyplot as plt
import easyclimate as ecl
from pathlib import Path
from .const_define import DOCS_DATA_PATH

mjo_data = ecl.open_tutorial_dataset("mjo_data").rename({"T": "time"})


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_draw_mjo_phase_space_basemap_AND_draw_mjo_phase_space():
    fig, ax = plt.subplots(figsize=(7.5, 7.5))

    ecl.field.equatorial_wave.draw_mjo_phase_space_basemap()
    ecl.field.equatorial_wave.draw_mjo_phase_space(
        mjo_data=mjo_data.sel(time=slice("2024-12-01", "2024-12-31")),
        rmm1_dim="RMM1",
        rmm2_dim="RMM2",
        time_dim="time",
    )
    return fig
