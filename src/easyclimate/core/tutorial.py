"""
Useful for:

* users learning easyclimate
* building tutorials in the documentation.

"""

# idea borrowed from xarray

from __future__ import annotations

import os
import pathlib
from typing import TYPE_CHECKING
from glob import glob
import xarray as xr
from pathlib import Path
import pandas as pd

if TYPE_CHECKING:
    from xarray.backends.api import T_Engine

__all__ = ["open_tutorial_dataset"]

_default_cache_dir_name = "easylimate_tutorial_data"
base_url = "https://github.com/shenyulu/easyclimate-data"
version = "main"


def _construct_cache_dir(path):
    import pooch

    if isinstance(path, os.PathLike):
        path = os.fspath(path)
    elif path is None:
        path = pooch.os_cache(_default_cache_dir_name)

    return path


external_urls = {}  # type: dict
file_formats = {
    "air_202201_mon_mean": 4,
    "hgt_202201_mon_mean": 4,
    "precip_202201_mon_mean": 4,
    "pressfc_202201_mon_mean": 4,
    "shum_202201_mon_mean": 4,
    "uwnd_202201_mon_mean": 4,
    "vwnd_202201_mon_mean": 4,
    "omega_202201_mon_mean": 4,
    "mini_HadISST_ice": 4,
    "PressQFF_202007271200_872": "csv",
    "pr_wtr_eatm_2022": 4,
    "sst_mnmean_oisst": 4,
    "hgt_day_ltm_1991_2020_0to6day": 4,
    "uwnd_day_ltm_1991_2020_0to6day": 4,
    "vwnd_day_ltm_1991_2020_0to6day": 4,
}


def _check_netcdf_engine_installed(name):
    version = file_formats.get(name)
    if version == 3:
        try:
            import scipy  # noqa
        except ImportError:
            try:
                import netCDF4  # noqa
            except ImportError:
                raise ImportError(
                    f"opening tutorial dataset {name} requires either scipy or "
                    "netCDF4 to be installed."
                )
    if version == 4:
        try:
            import h5netcdf  # noqa
        except ImportError:
            try:
                import netCDF4  # noqa
            except ImportError:
                raise ImportError(
                    f"opening tutorial dataset {name} requires either h5netcdf "
                    "or netCDF4 to be installed."
                )
    if version == "csv":
        try:
            import pandas  # noqa
        except ImportError:
            raise ImportError(
                f"opening tutorial dataset {name} requires pandas " " to be installed."
            )


def open_tutorial_dataset(
    name: str,
    cache: bool = True,
    cache_dir: None | str | os.PathLike = None,
    progressbar: bool = False,
    *,
    engine: T_Engine = None,
    **kws,
) -> xr.Dataset:
    """
    Open a dataset from the online repository (requires internet).

    If a local copy is found then always use that to avoid network traffic.

    Available datasets:

    * ``"air_202201_mon_mean"``: 2m air temperature of the NCEP reanalysis subset
    * ``"hgt_202201_mon_mean"``: Geopotential height of the NCEP reanalysis subset
    * ``"precip_202201_mon_mean"``: Precipitation of the NCEP reanalysis subset
    * ``"pressfc_202201_mon_mean"``: Mean sea surface pressure of the NCEP reanalysis subset
    * ``"shum_202201_mon_mean"``: Absolute humidity of the NCEP reanalysis subset
    * ``"uwnd_202201_mon_mean"``: Zonal wind of the NCEP reanalysis subset
    * ``"vwnd_202201_mon_mean"``: Meridional wind of the NCEP reanalysis subset
    * ``"omega_202201_mon_mean"``: Vertical velocity of the NCEP reanalysis subset
    * ``"mini_HadISST_ice"``: Hadley Centre Sea Ice and Sea Surface Temperature data set (HadISST) subset
    * ``"PressQFF_202007271200_872"``: Observational data from European stations (from https://github.com/EXCITED-CO2/xarray-regrid)
    * ``"pr_wtr_eatm_2022"``: Precipitable water of the NCEP reanalysis subset in the 2022
    * ``"sst_mnmean_oisst"``: NOAA Optimum Interpolation (OI) SST V2 (from https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html)

    Parameters
    ----------
    name : :py:class:`str <str>`
        Name of the file containing the dataset.
        e.g. 'air_202201_mon_mean'
    cache_dir : path-like, optional
        The directory in which to search for and write cached data.
    cache : dim: :py:class:`bool <bool>`, optional
        If True, then cache data locally for use on subsequent calls
    progressbar: :py:class:`bool <bool>`, default `False`.
        If True, will print a progress bar of the download to standard error (stderr). Requires `tqdm` to be installed.
    **kws : :py:class:`dict <dict>`, optional
        Passed to xarray.open_dataset

    Returns
    -------
    :py:class:`xarray.Dataset<xarray.Dataset>`

    Reference
    --------------
    - Kalnay et al.,The NCEP/NCAR 40-year reanalysis project, Bull. Amer. Meteor. Soc., 77, 437-470, 1996
    - Rayner, N. A.; Parker, D. E.; Horton, E. B.; Folland, C. K.; Alexander, L. V.; Rowell, D. P.; Kent, E. C.; Kaplan, A. (2003) Global analyses of sea surface temperature, sea ice, and night marine air temperature since the late nineteenth century J. Geophys. Res.Vol. 108, No. D14, 4407 10.1029/2002JD002670  (pdf ~9Mb)

    .. seealso::
        - :py:func:`xarray.tutorial.load_dataset<xarray.tutorial.load_dataset>`
        - :py:func:`xarray.open_dataset<xarray.open_dataset>`
        - :py:func:`xarray.load_dataset<xarray.load_dataset>`
    """
    try:
        import pooch
        import tqdm
    except ImportError as e:
        raise ImportError(
            "tutorial.open_dataset depends on `pooch` and `tqdm` to download and manage datasets."
            " To proceed please install `pooch` and `tqdm` by the pypi."
        ) from e

    logger = pooch.get_logger()
    logger.setLevel("WARNING")

    cache_dir = _construct_cache_dir(cache_dir)
    if name in external_urls:
        url = external_urls[name]
    else:
        path = pathlib.Path(name)
        if not path.suffix:
            # process the name
            default_extension = ".nc"
            if engine is None:
                _check_netcdf_engine_installed(name)
            path = path.with_suffix(default_extension)
        elif path.suffix == ".grib":
            if engine is None:
                engine = "cfgrib"
                try:
                    import cfgrib  # noqa
                except ImportError as e:
                    raise ImportError(
                        "Reading this tutorial dataset requires the cfgrib package."
                    ) from e

        url = f"{base_url}/raw/{version}/{path.name}"

    # retrieve the file
    filepath = pooch.retrieve(
        url=url, known_hash=None, path=cache_dir, progressbar=progressbar
    )

    if Path(filepath).suffix == ".nc" or Path(filepath).suffix == ".grib":
        ds = xr.open_dataset(filepath, engine=engine, **kws)
        if not cache:
            ds = ds.load()
            pathlib.Path(filepath).unlink()
    elif Path(filepath).suffix == ".csv" or Path(filepath).suffix == ".CSV":
        ds = pd.read_csv(filepath)

    return ds
