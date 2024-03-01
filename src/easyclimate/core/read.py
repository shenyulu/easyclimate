"""
Functions for read data.
"""

from __future__ import annotations

from glob import glob
import xarray as xr

__all__ = ["open_muliti_dataset"]


def open_muliti_dataset(files: str, dim: str, **kwargs) -> xr.Dataset:
    """
    Open multiple netCDF files without the need for xarray's necessary dimension checks

    Parameters
    ----------
    - ver1: Version number 1
    - ver2: Version number 2

    Returns
    -------
    :py:class:`int <int>`.

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
    if "chunks" not in kwargs:
        kwargs["chunks"] = "auto"

    # glob expands paths with * to a list of files, like the unix shell
    paths = sorted(glob(files))
    datasets = [xr.open_dataset(p, **kwargs) for p in paths]
    combined = xr.concat(datasets, dim)
    return combined
