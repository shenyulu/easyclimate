"""
pytest for stat.py
"""

import pytest
import matplotlib.pyplot as plt
import numpy as np
import easyclimate as ecl

sst_data = ecl.open_tutorial_dataset("sst_mnmean_oisst").sst
sst_data_anormaly = ecl.remove_seasonal_cycle_mean(sst_data)

nino34_index = ecl.field.air_sea_interaction.calc_index_nino34(sst_data, running_mean=0)
nino34_index_normalized = ecl.normalized.normalize_zscore(nino34_index)


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_calc_corr_spatial():
    sst_reg_nino34_result = ecl.calc_corr_spatial(
        sst_data_anormaly, x=nino34_index_normalized
    )

    fig, ax = plt.subplots()

    sst_reg_nino34_result.reg_coeff.plot.contourf(
        ax=ax,
        levels=np.linspace(-1.5, 1.5, 21),
    )

    ecl.plot.draw_significant_area_contourf(
        sst_reg_nino34_result.pvalue,
        thresh=0.01,
        ax=ax,
    )
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_calc_leadlag_corr_spatial():
    leadlag_result = ecl.calc_leadlag_corr_spatial(
        sst_data_anormaly, x=nino34_index_normalized, leadlag_array=[-6, 6]
    )

    fig, ax = plt.subplots(2)

    axi = ax[0]
    leadlag_result.isel(leadlag=0).reg_coeff.plot.contourf(
        ax=axi,
        levels=np.linspace(-1.5, 1.5, 21),
    )

    ecl.plot.draw_significant_area_contourf(
        leadlag_result.isel(leadlag=0).pvalue,
        thresh=0.01,
        ax=axi,
    )

    axi = ax[1]
    leadlag_result.isel(leadlag=1).reg_coeff.plot.contourf(
        ax=axi,
        levels=np.linspace(-1.5, 1.5, 21),
    )

    ecl.plot.draw_significant_area_contourf(
        leadlag_result.isel(leadlag=1).pvalue,
        thresh=0.01,
        ax=axi,
    )
    return fig
