"""
pytest for field.air_sea_interaction.index_pdo.py
"""

import pytest

import easyclimate as ecl
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from .const_define import DOCS_DATA_PATH

data_sst = xr.open_dataset(str(Path(DOCS_DATA_PATH, "test_input_oisst_data.nc")))["sst"]


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_calc_index_PDO_EOF1():
    pdo_index = ecl.field.air_sea_interaction.calc_index_PDO_EOF1(
        data_sst, normalized=True, detrend_spatial=True
    )
    pdo_index_normalized = ecl.normalized.normalize_zscore(pdo_index, dim="time")

    fig, ax = plt.subplots()
    ecl.plot.bar_plot_with_threshold(pdo_index_normalized)
    return fig
