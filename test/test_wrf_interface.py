"""
pytest for filter/wrf/interface.py
"""

import pytest
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import easyclimate as ecl
from pathlib import Path
from .const_define import DOCS_DATA_PATH


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_wrf_interface():
    path = str(Path(DOCS_DATA_PATH, "wrf", "wrfout_d01_2022-05-01_00_00_00.nc4"))
    data = xr.open_dataset(path)
    ncfile = ecl.wrf.transfer_xarray2nctype(data)

    #
    ncfile = ecl.wrf.open_wrf_data(path)

    # Get the Sea Level Pressure
    slp = ecl.wrf.getvar(ncfile, "slp")

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()})
    slp.plot(
        ax=ax, cbar_kwargs={"location": "bottom"}, transform=ecl.wrf.get_cartopy(slp)
    )
    ax.gridlines(draw_labels=True)
    return fig
