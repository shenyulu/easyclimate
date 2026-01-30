easyclimate.wrf
===============

.. py:module:: easyclimate.wrf


Submodules
----------

.. toctree::
   :maxdepth: 1

   /technical/api/easyclimate/wrf/interface/index


Attributes
----------

.. autoapisummary::

   easyclimate.wrf.xarray_enabled
   easyclimate.wrf.disable_xarray
   easyclimate.wrf.enable_xarray
   easyclimate.wrf.cartopy_enabled
   easyclimate.wrf.disable_cartopy
   easyclimate.wrf.enable_cartopy
   easyclimate.wrf.basemap_enabled
   easyclimate.wrf.disable_basemap
   easyclimate.wrf.enable_basemap
   easyclimate.wrf.pyngl_enabled
   easyclimate.wrf.enable_pyngl
   easyclimate.wrf.disable_pyngl
   easyclimate.wrf.set_cache_size
   easyclimate.wrf.get_cache_size
   easyclimate.wrf.omp_enabled
   easyclimate.wrf.ALL_TIMES
   easyclimate.wrf.Constants
   easyclimate.wrf.ConversionFactors
   easyclimate.wrf.ProjectionTypes
   easyclimate.wrf.default_fill
   easyclimate.wrf.OMP_SCHED_STATIC
   easyclimate.wrf.OMP_SCHED_DYNAMIC
   easyclimate.wrf.OMP_SCHED_GUIDED
   easyclimate.wrf.OMP_SCHED_AUTO
   easyclimate.wrf.destagger
   easyclimate.wrf.getvar
   easyclimate.wrf.xy
   easyclimate.wrf.interp1d
   easyclimate.wrf.interp2dxy
   easyclimate.wrf.interpz3d
   easyclimate.wrf.slp
   easyclimate.wrf.tk
   easyclimate.wrf.td
   easyclimate.wrf.rh
   easyclimate.wrf.uvmet
   easyclimate.wrf.smooth2d
   easyclimate.wrf.cape_2d
   easyclimate.wrf.cape_3d
   easyclimate.wrf.cloudfrac
   easyclimate.wrf.ctt
   easyclimate.wrf.dbz
   easyclimate.wrf.srhel
   easyclimate.wrf.udhel
   easyclimate.wrf.avo
   easyclimate.wrf.pvo
   easyclimate.wrf.eth
   easyclimate.wrf.wetbulb
   easyclimate.wrf.tvirtual
   easyclimate.wrf.omega
   easyclimate.wrf.pw
   easyclimate.wrf.DiagnosticError
   easyclimate.wrf.omp_set_num_threads
   easyclimate.wrf.omp_get_num_threads
   easyclimate.wrf.omp_get_max_threads
   easyclimate.wrf.omp_get_thread_num
   easyclimate.wrf.omp_get_num_procs
   easyclimate.wrf.omp_in_parallel
   easyclimate.wrf.omp_set_dynamic
   easyclimate.wrf.omp_get_dynamic
   easyclimate.wrf.omp_set_nested
   easyclimate.wrf.omp_get_nested
   easyclimate.wrf.omp_set_schedule
   easyclimate.wrf.omp_get_schedule
   easyclimate.wrf.omp_get_thread_limit
   easyclimate.wrf.omp_set_max_active_levels
   easyclimate.wrf.omp_get_max_active_levels
   easyclimate.wrf.omp_get_level
   easyclimate.wrf.omp_get_ancestor_thread_num
   easyclimate.wrf.omp_get_team_size
   easyclimate.wrf.omp_get_active_level
   easyclimate.wrf.omp_in_final
   easyclimate.wrf.omp_init_lock
   easyclimate.wrf.omp_init_nest_lock
   easyclimate.wrf.omp_destroy_lock
   easyclimate.wrf.omp_destroy_nest_lock
   easyclimate.wrf.omp_set_lock
   easyclimate.wrf.omp_set_nest_lock
   easyclimate.wrf.omp_unset_lock
   easyclimate.wrf.omp_unset_nest_lock
   easyclimate.wrf.omp_test_lock
   easyclimate.wrf.omp_test_nest_lock
   easyclimate.wrf.omp_get_wtime
   easyclimate.wrf.omp_get_wtick
   easyclimate.wrf.interplevel
   easyclimate.wrf.vertcross
   easyclimate.wrf.interpline
   easyclimate.wrf.vinterp
   easyclimate.wrf.xy_to_ll
   easyclimate.wrf.ll_to_xy
   easyclimate.wrf.xy_to_ll_proj
   easyclimate.wrf.ll_to_xy_proj
   easyclimate.wrf.viewitems
   easyclimate.wrf.viewkeys
   easyclimate.wrf.viewvalues
   easyclimate.wrf.isstr
   easyclimate.wrf.py2round
   easyclimate.wrf.py3range
   easyclimate.wrf.ucode
   easyclimate.wrf.to_np
   easyclimate.wrf.extract_global_attrs
   easyclimate.wrf.is_standard_wrf_var
   easyclimate.wrf.extract_dim
   easyclimate.wrf.extract_vars
   easyclimate.wrf.extract_times
   easyclimate.wrf.combine_files
   easyclimate.wrf.npbytes_to_str
   easyclimate.wrf.is_moving_domain
   easyclimate.wrf.is_staggered
   easyclimate.wrf.get_left_indexes
   easyclimate.wrf.iter_left_indexes
   easyclimate.wrf.get_right_slices
   easyclimate.wrf.get_proj_params
   easyclimate.wrf.from_args
   easyclimate.wrf.args_to_list
   easyclimate.wrf.arg_location
   easyclimate.wrf.psafilepath
   easyclimate.wrf.get_id
   easyclimate.wrf.from_var
   easyclimate.wrf.combine_dims
   easyclimate.wrf.either
   easyclimate.wrf.get_iterable
   easyclimate.wrf.IterWrapper
   easyclimate.wrf.is_coordvar
   easyclimate.wrf.latlon_coordvars
   easyclimate.wrf.is_mapping
   easyclimate.wrf.has_time_coord
   easyclimate.wrf.is_multi_file
   easyclimate.wrf.is_multi_time_req
   easyclimate.wrf.get_coord_pairs
   easyclimate.wrf.is_time_coord_var
   easyclimate.wrf.geo_bounds
   easyclimate.wrf.get_cartopy
   easyclimate.wrf.get_basemap
   easyclimate.wrf.get_pyngl
   easyclimate.wrf.cartopy_xlim
   easyclimate.wrf.cartopy_ylim
   easyclimate.wrf.latlon_coords
   easyclimate.wrf.ll_points
   easyclimate.wrf.pairs_to_latlon
   easyclimate.wrf.GeoBounds
   easyclimate.wrf.NullGeoBounds
   easyclimate.wrf.WrfProj
   easyclimate.wrf.NullProjection
   easyclimate.wrf.LambertConformal
   easyclimate.wrf.Mercator
   easyclimate.wrf.PolarStereographic
   easyclimate.wrf.LatLon
   easyclimate.wrf.RotatedLatLon
   easyclimate.wrf.getproj
   easyclimate.wrf.CoordPair
   easyclimate.wrf.to_xy_coords
   easyclimate.wrf.cache_item
   easyclimate.wrf.get_cached_item


Functions
---------

.. autoapisummary::

   easyclimate.wrf.open_wrf_data
   easyclimate.wrf.transfer_xarray2nctype


Package Contents
----------------

.. py:data:: xarray_enabled
   :value: False


.. py:data:: disable_xarray
   :value: None


.. py:data:: enable_xarray
   :value: None


.. py:data:: cartopy_enabled
   :value: False


.. py:data:: disable_cartopy
   :value: None


.. py:data:: enable_cartopy
   :value: None


.. py:data:: basemap_enabled
   :value: False


.. py:data:: disable_basemap
   :value: None


.. py:data:: enable_basemap
   :value: None


.. py:data:: pyngl_enabled
   :value: False


.. py:data:: enable_pyngl
   :value: None


.. py:data:: disable_pyngl
   :value: None


.. py:data:: set_cache_size
   :value: None


.. py:data:: get_cache_size
   :value: None


.. py:data:: omp_enabled
   :value: False


.. py:data:: ALL_TIMES
   :value: None


.. py:data:: Constants
   :value: None


.. py:data:: ConversionFactors
   :value: None


.. py:data:: ProjectionTypes
   :value: None


.. py:data:: default_fill
   :value: None


.. py:data:: OMP_SCHED_STATIC
   :value: None


.. py:data:: OMP_SCHED_DYNAMIC
   :value: None


.. py:data:: OMP_SCHED_GUIDED
   :value: None


.. py:data:: OMP_SCHED_AUTO
   :value: None


.. py:data:: destagger
   :value: None


.. py:data:: getvar

.. py:data:: xy
   :value: None


.. py:data:: interp1d
   :value: None


.. py:data:: interp2dxy
   :value: None


.. py:data:: interpz3d
   :value: None


.. py:data:: slp
   :value: None


.. py:data:: tk
   :value: None


.. py:data:: td
   :value: None


.. py:data:: rh
   :value: None


.. py:data:: uvmet
   :value: None


.. py:data:: smooth2d
   :value: None


.. py:data:: cape_2d
   :value: None


.. py:data:: cape_3d
   :value: None


.. py:data:: cloudfrac
   :value: None


.. py:data:: ctt
   :value: None


.. py:data:: dbz
   :value: None


.. py:data:: srhel
   :value: None


.. py:data:: udhel
   :value: None


.. py:data:: avo
   :value: None


.. py:data:: pvo
   :value: None


.. py:data:: eth
   :value: None


.. py:data:: wetbulb
   :value: None


.. py:data:: tvirtual
   :value: None


.. py:data:: omega
   :value: None


.. py:data:: pw
   :value: None


.. py:data:: DiagnosticError
   :value: None


.. py:data:: omp_set_num_threads
   :value: None


.. py:data:: omp_get_num_threads
   :value: None


.. py:data:: omp_get_max_threads
   :value: None


.. py:data:: omp_get_thread_num
   :value: None


.. py:data:: omp_get_num_procs
   :value: None


.. py:data:: omp_in_parallel
   :value: None


.. py:data:: omp_set_dynamic
   :value: None


.. py:data:: omp_get_dynamic
   :value: None


.. py:data:: omp_set_nested
   :value: None


.. py:data:: omp_get_nested
   :value: None


.. py:data:: omp_set_schedule
   :value: None


.. py:data:: omp_get_schedule
   :value: None


.. py:data:: omp_get_thread_limit
   :value: None


.. py:data:: omp_set_max_active_levels
   :value: None


.. py:data:: omp_get_max_active_levels
   :value: None


.. py:data:: omp_get_level
   :value: None


.. py:data:: omp_get_ancestor_thread_num
   :value: None


.. py:data:: omp_get_team_size
   :value: None


.. py:data:: omp_get_active_level
   :value: None


.. py:data:: omp_in_final
   :value: None


.. py:data:: omp_init_lock
   :value: None


.. py:data:: omp_init_nest_lock
   :value: None


.. py:data:: omp_destroy_lock
   :value: None


.. py:data:: omp_destroy_nest_lock
   :value: None


.. py:data:: omp_set_lock
   :value: None


.. py:data:: omp_set_nest_lock
   :value: None


.. py:data:: omp_unset_lock
   :value: None


.. py:data:: omp_unset_nest_lock
   :value: None


.. py:data:: omp_test_lock
   :value: None


.. py:data:: omp_test_nest_lock
   :value: None


.. py:data:: omp_get_wtime
   :value: None


.. py:data:: omp_get_wtick
   :value: None


.. py:data:: interplevel

.. py:data:: vertcross
   :value: None


.. py:data:: interpline
   :value: None


.. py:data:: vinterp
   :value: None


.. py:data:: xy_to_ll
   :value: None


.. py:data:: ll_to_xy
   :value: None


.. py:data:: xy_to_ll_proj
   :value: None


.. py:data:: ll_to_xy_proj
   :value: None


.. py:data:: viewitems
   :value: None


.. py:data:: viewkeys
   :value: None


.. py:data:: viewvalues
   :value: None


.. py:data:: isstr
   :value: None


.. py:data:: py2round
   :value: None


.. py:data:: py3range
   :value: None


.. py:data:: ucode
   :value: None


.. py:data:: to_np
   :value: None


.. py:data:: extract_global_attrs
   :value: None


.. py:data:: is_standard_wrf_var
   :value: None


.. py:data:: extract_dim
   :value: None


.. py:data:: extract_vars
   :value: None


.. py:data:: extract_times
   :value: None


.. py:data:: combine_files
   :value: None


.. py:data:: npbytes_to_str
   :value: None


.. py:data:: is_moving_domain
   :value: None


.. py:data:: is_staggered
   :value: None


.. py:data:: get_left_indexes
   :value: None


.. py:data:: iter_left_indexes
   :value: None


.. py:data:: get_right_slices
   :value: None


.. py:data:: get_proj_params
   :value: None


.. py:data:: from_args
   :value: None


.. py:data:: args_to_list
   :value: None


.. py:data:: arg_location
   :value: None


.. py:data:: psafilepath
   :value: None


.. py:data:: get_id
   :value: None


.. py:data:: from_var
   :value: None


.. py:data:: combine_dims
   :value: None


.. py:data:: either
   :value: None


.. py:data:: get_iterable
   :value: None


.. py:data:: IterWrapper
   :value: None


.. py:data:: is_coordvar
   :value: None


.. py:data:: latlon_coordvars
   :value: None


.. py:data:: is_mapping
   :value: None


.. py:data:: has_time_coord
   :value: None


.. py:data:: is_multi_file
   :value: None


.. py:data:: is_multi_time_req
   :value: None


.. py:data:: get_coord_pairs
   :value: None


.. py:data:: is_time_coord_var
   :value: None


.. py:data:: geo_bounds
   :value: None


.. py:data:: get_cartopy
   :value: None


.. py:data:: get_basemap
   :value: None


.. py:data:: get_pyngl
   :value: None


.. py:data:: cartopy_xlim
   :value: None


.. py:data:: cartopy_ylim
   :value: None


.. py:data:: latlon_coords
   :value: None


.. py:data:: ll_points
   :value: None


.. py:data:: pairs_to_latlon
   :value: None


.. py:data:: GeoBounds
   :value: None


.. py:data:: NullGeoBounds
   :value: None


.. py:data:: WrfProj
   :value: None


.. py:data:: NullProjection
   :value: None


.. py:data:: LambertConformal
   :value: None


.. py:data:: Mercator
   :value: None


.. py:data:: PolarStereographic
   :value: None


.. py:data:: LatLon
   :value: None


.. py:data:: RotatedLatLon
   :value: None


.. py:data:: getproj
   :value: None


.. py:data:: CoordPair
   :value: None


.. py:data:: to_xy_coords
   :value: None


.. py:data:: cache_item
   :value: None


.. py:data:: get_cached_item
   :value: None


.. py:function:: open_wrf_data(filename: str | os.PathLike, mode: Literal['r', 'w', 'r+', 'a', 'x', 'rs', 'ws', 'r+s', 'as'] = 'r', clobber: bool = True, format: Literal['NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_CLASSIC', 'NETCDF3_64BIT_OFFSET', 'NETCDF3_64BIT_DATA'] = 'NETCDF4', diskless: bool = False, persist: bool = False, keepweakref: bool = False, memory: typing_extensions.Buffer | int | None = None, encoding: str | None = None, parallel: bool = False, comm: Any = None, info: Any = None, auto_complex: bool = False, **kwargs) -> netCDF4.Dataset

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


.. py:function:: transfer_xarray2nctype(ds: xarray.Dataset) -> netCDF4.Dataset

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


