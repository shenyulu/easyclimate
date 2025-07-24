"""
pytest for stat.py calc_lead_lag_correlation_coefficients
"""

import pytest

import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import easyclimate as ecl

pc1 = xr.DataArray(
    np.array(
        [
            -11.02909005,
            -9.24363696,
            -8.9240402,
            -4.05222314,
            0.36699088,
            -1.16396849,
            -3.48129717,
            -3.89368637,
            3.90733285,
            10.50904935,
            15.79242665,
            17.42215435,
            21.28036026,
            23.90369037,
            25.09963444,
            26.04435696,
            24.89421697,
            17.65926959,
            13.36064014,
            14.7511998,
        ]
    ),
    # data_in.pc1.values[:20],
    dims="time",
    coords={
        "time": pd.to_datetime(
            [
                "1986-05-01",
                "1986-05-02",
                "1986-05-03",
                "1986-05-04",
                "1986-05-05",
                "1986-05-06",
                "1986-05-07",
                "1986-05-08",
                "1986-05-09",
                "1986-05-10",
                "1986-05-11",
                "1986-05-12",
                "1986-05-13",
                "1986-05-14",
                "1986-05-15",
                "1986-05-16",
                "1986-05-17",
                "1986-05-18",
                "1986-05-19",
                "1986-05-20",
            ]
        )
    },
)

pc2 = xr.DataArray(
    np.array(
        [
            9.06924315,
            5.79843584,
            7.45809984,
            10.56795327,
            11.34306118,
            6.00196093,
            -0.69995022,
            -0.97253475,
            -5.78444372,
            -8.65616443,
            -2.98160504,
            1.54400324,
            5.7888796,
            2.64337392,
            -6.38860926,
            -9.70990684,
            -10.2995918,
            -14.12106694,
            -12.89854077,
            -6.8267084,
        ]
    ),
    # data_in.pc2.values[:20],
    dims="time",
    coords={
        "time": pd.to_datetime(
            [
                "1986-05-01",
                "1986-05-02",
                "1986-05-03",
                "1986-05-04",
                "1986-05-05",
                "1986-05-06",
                "1986-05-07",
                "1986-05-08",
                "1986-05-09",
                "1986-05-10",
                "1986-05-11",
                "1986-05-12",
                "1986-05-13",
                "1986-05-14",
                "1986-05-15",
                "1986-05-16",
                "1986-05-17",
                "1986-05-18",
                "1986-05-19",
                "1986-05-20",
            ]
        )
    },
)


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_calc_lead_lag_correlation_coefficients():
    fig, ax = plt.subplots()

    pcs = {"PC1": pc1, "PC2": pc2}
    pairs = [("PC1_vs_PC2", "PC1", "PC2")]

    corr_da = ecl.calc_lead_lag_correlation_coefficients(
        pcs=pcs, pairs=pairs, max_lag=30
    )

    for pair_name in corr_da.data_vars:
        corr_da[pair_name].plot(label=pair_name)
    return fig
