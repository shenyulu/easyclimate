"""
pytest for field/monsoon/bsiso.py
"""

import pytest

import easyclimate as ecl
from easyclimate.core.datanode import DataNode
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path
from .const_define import TEST_DATA_PATH

data_input = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_asm_area_olr_uv.nc")))
test_olr_data = data_input["olr"]
test_u850_data = data_input["u"]
test_v850_data = data_input["v"]
result_node = ecl.field.monsoon.calc_bsiso_analysis(
    test_olr_data, test_u850_data, test_v850_data
)


def test_calc_bsiso_analysis():
    assert isinstance(result_node, DataNode)


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_draw_bsiso1_phase_space_basemap():
    fig, ax = plt.subplots()

    ecl.field.monsoon.draw_bsiso1_phase_space_basemap()

    ecl.field.monsoon.draw_bsiso_phase_space(
        result_node["PC_normalized/pc12"].isel(time=slice(0, 10)),
        color="r",
        x_dim="PC2",
        y_dim="PC1",
    )

    ecl.field.monsoon.draw_bsiso_phase_space(
        result_node["PC_normalized/pc12"].isel(time=slice(100, 130)),
        color="b",
        x_dim="PC2",
        y_dim="PC1",
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_draw_bsiso2_phase_space_basemap():
    fig, ax = plt.subplots()

    ecl.field.monsoon.draw_bsiso2_phase_space_basemap()

    ecl.field.monsoon.draw_bsiso_phase_space(
        result_node["PC_normalized/pc34"].isel(time=slice(0, 10)),
        color="r",
        x_dim="PC4",
        y_dim="PC3",
    )

    ecl.field.monsoon.draw_bsiso_phase_space(
        result_node["PC_normalized/pc34"].isel(time=slice(100, 130)),
        color="b",
        x_dim="PC4",
        y_dim="PC3",
    )
    return fig
