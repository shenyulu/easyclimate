"""
pytest for plot.taylor_diagrams.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from .const_define import DOCS_DATA_PATH

js_data = xr.open_dataset(
    str(Path(DOCS_DATA_PATH, "js_H09_20250617_0500.nc")), decode_timedelta=False
)


@pytest.mark.mpl_image_compare(remove_text=True)
def test_get_stretched_rgb_data_base():
    rgb_result = ecl.satellite.get_stretched_rgb_data(
        js_data, r_band="albedo_03", g_band="albedo_02", b_band="albedo_01"
    )

    fig, ax = plt.subplots()
    rgb_result.plot.imshow(ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("")
    return fig
