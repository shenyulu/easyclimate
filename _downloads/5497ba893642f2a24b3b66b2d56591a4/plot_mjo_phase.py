# -*- coding: utf-8 -*-
"""
.. _mjo_phase_diagram_example:

MJO Phase Space Trajectory Diagram
===================================

Before proceeding with all the steps, first import some necessary libraries and packages
"""
import xarray as xr
import pandas as pd
import easyclimate as ecl
import matplotlib.pyplot as plt

# %%
# Load MJO phase space data
#
# .. code-block:: python
#
#       mjo_ds = xr.open_dataset('http://iridl.ldeo.columbia.edu/SOURCES/.BoM/.MJO/.RMM/dods',
#                                decode_times=False)
#       T = mjo_ds.T.values
#       mjo_ds['T'] = pd.date_range("1974-06-01", periods=len(T))
#       mjo_ds = ecl.utility.get_compress_xarraydata(mjo_ds)
#       mjo_ds.to_netcdf("mjo_data.nc")
#
# Modify the name of the time parameter
#
# .. tip::
#
#   You can download following datasets here: :download:`Download mjo_data.nc <https://raw.githubusercontent.com/shenyulu/easyclimate/refs/heads/main/docs/source/dynamic_docs/mjo_data.nc>`
#

mjo_data = xr.open_dataset("mjo_data.nc").rename({"T": "time"})
mjo_data

# %%
# Draw MJO Phase Space Trajectory Diagram

fig, ax = plt.subplots(figsize = (7.5, 7.5))

ecl.field.equatorial_wave.draw_mjo_phase_space_basemap()
ecl.field.equatorial_wave.draw_mjo_phase_space(
    mjo_data = mjo_data.sel(time = slice('2024-12-01', '2024-12-31')),
    rmm1_dim = "RMM1",
    rmm2_dim = "RMM2",
    time_dim = "time"
)
