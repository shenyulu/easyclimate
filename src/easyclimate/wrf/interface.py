"""
WRF-python interface

.. seealso::

    This is the **built-in compiled** `wrf-python <https://wrf-python.readthedocs.io/en/latest/index.html>`__ package.
    Please refer to the `User API <https://wrf-python.readthedocs.io/en/latest/user_api/index.html>`__ in the `wrf-python <https://wrf-python.readthedocs.io/en/latest/index.html>`__ package for related routines.
"""

from __future__ import annotations

import xarray as xr
import os
from netCDF4 import Dataset
from typing import Literal, Any
from typing_extensions import Buffer

__all__ = ["open_wrf_data", "transfer_xarray2nctype"]


def open_wrf_data(
    filename: str | os.PathLike,
    mode: Literal["r", "w", "r+", "a", "x", "rs", "ws", "r+s", "as"] = "r",
    clobber: bool = True,
    format: Literal[
        "NETCDF4",
        "NETCDF4_CLASSIC",
        "NETCDF3_CLASSIC",
        "NETCDF3_64BIT_OFFSET",
        "NETCDF3_64BIT_DATA",
    ] = "NETCDF4",
    diskless: bool = False,
    persist: bool = False,
    keepweakref: bool = False,
    memory: Buffer | int | None = None,
    encoding: str | None = None,
    parallel: bool = False,
    comm: Any = None,
    info: Any = None,
    auto_complex: bool = False,
    **kwargs,
) -> Dataset:
    """
    Open WRF (Weather Research and Forecasting) Model Output data.

    Parameters
    ----------
    filename: :py:class:`str<str>`.
        Name of netCDF file to hold dataset. Can also be a python 3 pathlib instance or the URL of an OpenDAP dataset. When memory is set this is just used to set the `filepath()`.
    mode: ["r", "w", "r+", "a", "x", "rs", "ws", "r+s", "as"], default: "r".
        access mode. ``r`` means read-only; no data can be modified. ``w`` means write; a new file is created,
        an existing file with the same name is deleted. ``x`` means write, but fail if an existing file with the same name already exists.
        ``a`` and ``r+`` mean append; an existing file is opened for reading and writing, if file does not exist already, one is created.
        Appending ``s`` to modes ``r``, ``w``, ``r+`` or a will enable unbuffered shared access to ``NETCDF3_CLASSIC``,
        ``NETCDF3_64BIT_OFFSET`` or ``NETCDF3_64BIT_DATA`` formatted files. Unbuffered access may be useful even if you don't need shared access,
        since it may be faster for programs that don't access data sequentially. This option is ignored for ``NETCDF4`` and ``NETCDF4_CLASSIC`` formatted files.
    clobber: :py:class:`bool<bool>`, default: ``True``.
        If True (default), opening a file with ``mode='w'`` will clobber an existing file with the same name.
        If False, an exception will be raised if a file with the same name already exists. ``mode=x`` is identical to ``mode=w`` with ``clobber=False``.
    format: ["NETCDF4", "NETCDF4_CLASSIC", "NETCDF3_CLASSIC", "NETCDF3_64BIT_OFFSET", "NETCDF3_64BIT_DATA"], default: "NETCDF4".
        underlying file format (one of ``'NETCDF4'``, ``'NETCDF4_CLASSIC'``, ``'NETCDF3_CLASSIC'``, ``'NETCDF3_64BIT_OFFSET'`` or ``'NETCDF3_64BIT_DATA'``.
        Only relevant if ``mode = 'w'`` (if ``mode = 'r','a'`` or ``'r+'`` the file format is automatically detected).
        Default ``'NETCDF4'``, which means the data is stored in an HDF5 file, using netCDF 4 API features.
        Setting ``format='NETCDF4_CLASSIC'`` will create an HDF5 file, using only netCDF 3 compatible API features.
        netCDF 3 clients must be recompiled and linked against the netCDF 4 library to read files in ``NETCDF4_CLASSIC`` format.
        ``'NETCDF3_CLASSIC'`` is the classic netCDF 3 file format that does not handle 2+ Gb files.
        ``'NETCDF3_64BIT_OFFSET'`` is the 64-bit offset version of the netCDF 3 file format,
        which fully supports 2+ GB files, but is only compatible with clients linked against netCDF version 3.6.0 or later.
        ``'NETCDF3_64BIT_DATA'`` is the 64-bit data version of the netCDF 3 file format,
        which supports 64-bit dimension sizes plus unsigned and 64 bit integer data types,
        but is only compatible with clients linked against netCDF version 4.4.0 or later.
    diskless: :py:class:`bool<bool>`, default: ``False``.
        If ``True``, create diskless (in-core) file. This is a feature added to the C library after the netcdf-4.2 release.
        If you need to access the memory buffer directly, use the in-memory feature instead (see ``memory`` kwarg).
    persist:  :py:class:`bool<bool>`, default: ``False``.
        if ``diskless=True``, persist file to disk when closed (default ``False``).
    keepweakref: :py:class:`bool<bool>`, default: ``False``.
        if ``True``, child Dimension and Variable instances will keep weak references to the parent Dataset or Group object.
        Default is ``False``, which means strong references will be kept.
        Having Dimension and Variable instances keep a strong reference to the parent Dataset instance,
        which in turn keeps a reference to child Dimension and Variable instances, creates circular references.
        Circular references complicate garbage collection,
        which may mean increased memory usage for programs that create may Dataset instances with lots of Variables.
        It also will result in the Dataset object never being deleted, which means it may keep open files alive as well.
        Setting ``keepweakref=True`` allows Dataset instances to be garbage collected as soon as they go out of scope,
        potentially reducing memory usage and open file handles. However, in many cases this is not desirable,
        since the associated Variable instances may still be needed, but are rendered unusable when the parent Dataset instance is garbage collected.
    memory: Buffer or :py:class:`int<int>`.
        if not ``None``, create or open an in-memory Dataset.
        If mode = ``r``, the memory kwarg must contain a memory buffer object (an object that supports the python buffer interface).
        The Dataset will then be created with contents taken from this block of memory.
        If mode = ``w``, the memory kwarg should contain the anticipated size of the Dataset in bytes (used only for NETCDF3 files).
        A memory buffer containing a copy of the Dataset is returned by the ``Dataset.close()`` method.
        Requires netcdf-c version 4.4.1 for mode=``r`` netcdf-c 4.6.2 for mode=``w``.
        To persist the file to disk, the raw bytes from the returned buffer can be written into a binary file.
        The Dataset can also be re-opened using this memory buffer.
    encoding: :py:class:`str<str>`.
        encoding used to encode filename string into bytes. Default is None (``sys.getdefaultfileencoding()`` is used).
    parallel: :py:class:`bool<bool>`, default: ``False``.
        open for parallel access using MPI (requires mpi4py and parallel-enabled netcdf-c and hdf5 libraries).
        Default is ``False``. If ``True``, ``comm`` and ``info`` kwargs may also be specified.
    comm: Any.
        MPI_Comm object for parallel access. Default ``None``, which means MPI_COMM_WORLD will be used. Ignored if ``parallel=False``.
    info: Any.
        MPI_Info object for parallel access. Default ``None``, which means MPI_INFO_NULL will be used. Ignored if ``parallel=False``.
    auto_complex: :py:class:`bool<bool>`, default: ``False``.
        if ``True``, then automatically convert complex number types.

    Returns
    -------
    Data (:py:class:`netCDF4.Dataset<netCDF4.Dataset>`).

    .. seealso::
        https://unidata.github.io/netcdf4-python/#netCDF4.Dataset

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_wrf_tutorial.py
    """

    ncfile = Dataset(
        filename,
        mode=mode,
        clobber=clobber,
        format=format,
        diskless=diskless,
        persist=persist,
        keepweakref=keepweakref,
        memory=memory,
        encoding=encoding,
        parallel=parallel,
        comm=comm,
        info=info,
        auto_complex=auto_complex,
        **kwargs,
    )
    return ncfile


def transfer_xarray2nctype(ds: xr.Dataset) -> Dataset:
    """
    Transfer WRF data from :py:class:`xarray.Dataset<xarray.Dataset>` to :py:class:`netCDF4.Dataset<netCDF4.Dataset>`.

    Parameters
    ----------
    ds: :py:class:`xarray.Dataset<xarray.Dataset>`.
        WRF data read by xarray engine.

    Returns
    -------
    Data (:py:class:`netCDF4.Dataset<netCDF4.Dataset>`).

    .. seealso::
        - https://wrf-python.readthedocs.io/en/latest/faq.html
        - https://github.com/pydata/xarray/issues/5175
    """
    return ds._close.__self__.ds
