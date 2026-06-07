# -*- coding: utf-8 -*-
"""
.. _smooth_daily_cycle_example:

Smooth Mean Daily Annual Cycle
============================================

Calculates a smooth mean daily annual cycle.

.. seealso::

    https://www.ncl.ucar.edu/Document/Functions/Contributed/smthClmDayTLL.shtml

Before proceeding with all the steps, first import some necessary libraries and packages
"""
import xarray as xr
import easyclimate as ecl
import matplotlib.pyplot as plt

# %%
# Preprocessed data
#
#
# .. tip::
#
#   You can download following datasets here:
#
#   - :download:`Download olr-daily_v01r02_19800101_20231231.nc (3.34 GB) <https://huggingface.co/datasets/shenyulu/easyclimate/resolve/main/tutorial_data/olr-daily_v01r02_19800101_20231231.nc>`
#   - :download:`Download olr_daily_annual_cycle_mean.nc <https://raw.githubusercontent.com/shenyulu/easyclimate/refs/heads/main/docs/source/dynamic_docs/olr_daily_annual_cycle_mean.nc>`
#
#
# .. code-block:: python
#
#       lats, latn = -20, 20
#
#       olr_data = xr.open_dataset('olr-daily_v01r02_19800101_20231231.nc', chunks='auto').sel(lat=slice(lats,latn)).olr
#       olr_data_daily_annual_cycle_mean = ecl.calc_daily_annual_cycle_mean(olr_data).thin(lat = 2, lon = 5).compute()
#
#       olr_data_daily_annual_cycle_mean = ecl.utility.get_compress_xarraydata(olr_data_daily_annual_cycle_mean)
#       olr_data_daily_annual_cycle_mean.to_netcdf("olr_daily_annual_cycle_mean.nc")
#
# Here, we directly load the data of mean daily annual cycle for outgoing longwave radiation (OLR).
#

olr_data_daily_annual_cycle_mean = ecl.open_tutorial_dataset("olr_daily_annual_cycle_mean").olr
olr_data_daily_annual_cycle_mean

# %%
# Simply do a regional average for the tropics

olr_data_ave = olr_data_daily_annual_cycle_mean.mean(dim = ("lon", "lat"))
olr_data_ave.plot()


# %%
# By removing excess noise in the mean daily annual cycle through the way of :py:func:`easyclimate.smooth_daily_annual_cycle <easyclimate.smooth_daily_annual_cycle>`, we can plot the following:

# sphinx_gallery_thumbnail_number = -1
olr_data_ave_smoothed = ecl.smooth_daily_annual_cycle(olr_data_daily_annual_cycle_mean).mean(dim = ("lon", "lat"))

olr_data_ave.plot(label = "Daily annual cycle mean")
olr_data_ave_smoothed.plot(label = "Smoothed daily annual cycle mean")
plt.legend()
