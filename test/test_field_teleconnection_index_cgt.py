"""
pytest for field/teleconnection/index_CGT.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path
from .const_define import TEST_DATA_PATH

z200_data = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_z200_mon.nc")))["hgt"]


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_calc_index_CGT_1point_Ding_Wang_2005():
    cgt_index1 = ecl.field.teleconnection.calc_index_CGT_1point_Ding_Wang_2005(
        z200_data
    )
    cgt_index1 = ecl.get_specific_months_data(cgt_index1, [6, 7, 8, 9])
    cgt_index1_normalized = ecl.normalized.timeseries_normalize_zscore(
        cgt_index1, dim="time"
    )

    z200_anomaly_data = ecl.remove_seasonal_cycle_mean(z200_data)
    z200_anormaly_JJAS = ecl.get_specific_months_data(
        z200_anomaly_data, month_array=[6, 7, 8, 9]
    )

    z200_reg_cgt_result1 = ecl.calc_corr_spatial(
        z200_anormaly_JJAS, cgt_index1_normalized
    )

    fig, ax = plt.subplots()
    z200_reg_cgt_result1.reg_coeff.plot.contourf(levels=21)
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_calc_index_CGT_NH_Ding_Wang_2005():
    cgt_monthly_index = ecl.field.teleconnection.calc_index_CGT_NH_Ding_Wang_2005(
        z200_data, output_freq="monthly"
    )
    cgt_seasonally_index = ecl.field.teleconnection.calc_index_CGT_NH_Ding_Wang_2005(
        z200_data, output_freq="seasonally"
    )

    def get_anormaly_dataset(ds, freq):
        ds_anormaly = ecl.remove_seasonal_cycle_mean(ds)
        ds_JJAS_anormaly = ecl.get_specific_months_data(
            ds_anormaly, month_array=[6, 7, 8, 9]
        )
        if freq == "monthly":
            return ecl.calc_detrend_spatial(ds_JJAS_anormaly)
        elif freq == "seasonally":
            yearly_mean = ecl.calc_yearly_climatological_mean(ds_JJAS_anormaly)
            return ecl.calc_detrend_spatial(yearly_mean)
        else:
            return 0

    z200_anormaly_JJAS_monthly_data = get_anormaly_dataset(z200_data, freq="monthly")
    z200_anormaly_JJAS_seasonally_data = get_anormaly_dataset(
        z200_data, freq="seasonally"
    )

    z200_reg_cgt_result1 = ecl.calc_corr_spatial(
        z200_anormaly_JJAS_monthly_data, cgt_monthly_index
    )
    z200_reg_cgt_result1 = ecl.plot.add_lon_cyclic(z200_reg_cgt_result1, 2.5)

    z200_reg_cgt_result2 = ecl.calc_corr_spatial(
        z200_anormaly_JJAS_seasonally_data, cgt_seasonally_index
    )
    z200_reg_cgt_result2 = ecl.plot.add_lon_cyclic(z200_reg_cgt_result2, 2.5)
    z200_reg_cgt_result2

    fig, ax = plt.subplots(2)
    z200_reg_cgt_result1.reg_coeff.plot.contourf(ax=ax[0], levels=21)
    z200_reg_cgt_result2.reg_coeff.plot.contourf(ax=ax[1], levels=21)
    return fig
