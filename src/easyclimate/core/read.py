"""
Functions for read data.
"""
from __future__ import annotations

from glob import glob
import xarray as xr
import pooch

def open_muliti_dataset(files: str, dim: str, **kwargs) -> xr.Dataset:
    """
    Compare python library versions.

    .. attention::
        - Only for incoming version numbers without alphabetic characters.
        - Based on this method, the version number comparison should result in the following `"10.12.2.6.5">"10.12.2.6"`.

    Parameters
    ----------
    - ver1: Version number 1
    - ver2: Version number 2

    Returns
    -------
    :py:class:`int<python.int>`.

    .. note::
        If `ver1<ver2`, return `-1`; If `ver1=ver2`, return `0`; If `ver1>ver2`, return `1`.

    Examples
    --------

    .. code:: python

        >>> import easyclimate as ecl
        >>> result = assert_compared_version("10.12.2.6.5", "10.12.2.6")
        >>> print(result)
        1

    .. note::
        - https://medium.com/pangeo/accessing-netcdf-and-grib-file-collections-as-cloud-native-virtual-datasets-using-kerchunk-625a2d0a9191
        - https://github.com/fsspec/kerchunk/issues/240
    """

    # Dask is used by default for multi-data file reads for lazy computation
    kwargs = kwargs.get('chunks', 'auto')

    # glob expands paths with * to a list of files, like the unix shell
    paths = sorted(glob(files))
    datasets = [xr.open_dataset(p, **kwargs) for p in paths]
    combined = xr.concat(datasets, dim)
    return combined