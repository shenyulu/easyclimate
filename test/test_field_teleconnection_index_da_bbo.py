"""
pytest for field/teleconnection/index_da.py and field/teleconnection/index_bbo.py
"""

import pytest

import easyclimate as ecl
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path
from .const_define import DOCS_DATA_PATH
from easyclimate.field.teleconnection import (
    calc_index_DA_EOF2_Wu_2006,
    calc_index_BBO_EOF3_Wu_2007,
)

slp_data = ecl.open_tutorial_dataset("slp_monmean_NH")["slp"]


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_calc_index_DA_EOF2_Wu_2006():
    slp_data_anormaly = ecl.remove_seasonal_cycle_mean(slp_data)

    index_ad = calc_index_DA_EOF2_Wu_2006(slp_data) * (-1)

    ad_lower_time = ecl.get_time_exceed_index_lower_bound(index_ad, -index_ad.std())
    ad_upper_time = ecl.get_time_exceed_index_upper_bound(index_ad, index_ad.std())

    ad_minus = slp_data_anormaly.sel(time=ad_lower_time).mean(dim="time")
    ad_plus = slp_data_anormaly.sel(time=ad_upper_time).mean(dim="time")

    fig, ax = plt.subplots(2)
    ad_minus.plot.contourf(ax=ax[0], levels=21, cbar_kwargs={"label": ""})
    ad_plus.plot.contourf(ax=ax[1], levels=21, cbar_kwargs={"label": ""})

    for axi in ax.flat:
        axi.set_xlabel("")
        axi.set_ylabel("")
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_calc_index_BBO_EOF3_Wu_2007():
    slp_data_anormaly = ecl.remove_seasonal_cycle_mean(slp_data)

    index_bbo = calc_index_BBO_EOF3_Wu_2007(slp_data) * (-1)

    bbo_lower_time = ecl.get_time_exceed_index_lower_bound(index_bbo, -index_bbo.std())
    bbo_upper_time = ecl.get_time_exceed_index_upper_bound(index_bbo, index_bbo.std())

    bbo_minus = slp_data_anormaly.sel(time=bbo_lower_time).mean(dim="time")
    bbo_plus = slp_data_anormaly.sel(time=bbo_upper_time).mean(dim="time")

    fig, ax = plt.subplots(2)
    bbo_minus.plot.contourf(ax=ax[0], levels=21, cbar_kwargs={"label": ""})
    bbo_plus.plot.contourf(ax=ax[1], levels=21, cbar_kwargs={"label": ""})

    for axi in ax.flat:
        axi.set_xlabel("")
        axi.set_ylabel("")
    return fig
