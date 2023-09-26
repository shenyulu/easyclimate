"""
tutorial
"""
from __future__ import annotations
import pooch
from pathlib import Path
import warnings

def get_tuturial_data(cache_dir = None):
    if cache_dir == None:
        cache_dir = Path.home() / '.eclcache'
        if cache_dir.exists() == False:
            # create folder
            cache_dir.mkdir(parents=True)
            warnings.warn(f"Create folder at '{str(cache_dir)}'. Ignoring.")
        cache_dir = str(cache_dir)

    # uwnd data 
    url = "https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/Monthlies/pressure/uwnd.mon.mean.nc"
    pooch.retrieve(url, known_hash = None, fname = "uwnd.mon.mean.nc", path = cache_dir)
            