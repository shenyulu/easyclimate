easyclimate.backend
===================

.. py:module:: easyclimate.backend


Attributes
----------

.. autoapisummary::

   easyclimate.backend.CURRENT_PLATFORM
   easyclimate.backend.spharm
   easyclimate.backend.VectorWind
   easyclimate.backend.aeronoskin
   easyclimate.backend.aeroskin
   easyclimate.backend.human_index_mod
   easyclimate.backend.human_index_mod_old
   easyclimate.backend._ecl_redfit
   easyclimate.backend._ecl_redfit_x
   easyclimate.backend._vinth2p_dp
   easyclimate.backend._vinth2p_ecmwf
   easyclimate.backend._vintp2p_ecmwf
   easyclimate.backend._wet_bulb_temperature
   easyclimate.backend.dvibeta
   easyclimate.backend.dvrfidf
   easyclimate.backend.ddvfidf
   easyclimate.backend.xarray_enabled
   easyclimate.backend.disable_xarray
   easyclimate.backend.enable_xarray
   easyclimate.backend.cartopy_enabled
   easyclimate.backend.disable_cartopy
   easyclimate.backend.enable_cartopy
   easyclimate.backend.basemap_enabled
   easyclimate.backend.disable_basemap
   easyclimate.backend.enable_basemap
   easyclimate.backend.pyngl_enabled
   easyclimate.backend.enable_pyngl
   easyclimate.backend.disable_pyngl
   easyclimate.backend.set_cache_size
   easyclimate.backend.get_cache_size
   easyclimate.backend.omp_enabled
   easyclimate.backend.ALL_TIMES
   easyclimate.backend.Constants
   easyclimate.backend.ConversionFactors
   easyclimate.backend.ProjectionTypes
   easyclimate.backend.default_fill
   easyclimate.backend.OMP_SCHED_STATIC
   easyclimate.backend.OMP_SCHED_DYNAMIC
   easyclimate.backend.OMP_SCHED_GUIDED
   easyclimate.backend.OMP_SCHED_AUTO
   easyclimate.backend.destagger
   easyclimate.backend.getvar
   easyclimate.backend.xy
   easyclimate.backend.interp1d
   easyclimate.backend.interp2dxy
   easyclimate.backend.interpz3d
   easyclimate.backend.slp
   easyclimate.backend.tk
   easyclimate.backend.td
   easyclimate.backend.rh
   easyclimate.backend.uvmet
   easyclimate.backend.smooth2d
   easyclimate.backend.cape_2d
   easyclimate.backend.cape_3d
   easyclimate.backend.cloudfrac
   easyclimate.backend.ctt
   easyclimate.backend.dbz
   easyclimate.backend.srhel
   easyclimate.backend.udhel
   easyclimate.backend.avo
   easyclimate.backend.pvo
   easyclimate.backend.eth
   easyclimate.backend.wetbulb
   easyclimate.backend.tvirtual
   easyclimate.backend.omega
   easyclimate.backend.pw
   easyclimate.backend.DiagnosticError
   easyclimate.backend.omp_set_num_threads
   easyclimate.backend.omp_get_num_threads
   easyclimate.backend.omp_get_max_threads
   easyclimate.backend.omp_get_thread_num
   easyclimate.backend.omp_get_num_procs
   easyclimate.backend.omp_in_parallel
   easyclimate.backend.omp_set_dynamic
   easyclimate.backend.omp_get_dynamic
   easyclimate.backend.omp_set_nested
   easyclimate.backend.omp_get_nested
   easyclimate.backend.omp_set_schedule
   easyclimate.backend.omp_get_schedule
   easyclimate.backend.omp_get_thread_limit
   easyclimate.backend.omp_set_max_active_levels
   easyclimate.backend.omp_get_max_active_levels
   easyclimate.backend.omp_get_level
   easyclimate.backend.omp_get_ancestor_thread_num
   easyclimate.backend.omp_get_team_size
   easyclimate.backend.omp_get_active_level
   easyclimate.backend.omp_in_final
   easyclimate.backend.omp_init_lock
   easyclimate.backend.omp_init_nest_lock
   easyclimate.backend.omp_destroy_lock
   easyclimate.backend.omp_destroy_nest_lock
   easyclimate.backend.omp_set_lock
   easyclimate.backend.omp_set_nest_lock
   easyclimate.backend.omp_unset_lock
   easyclimate.backend.omp_unset_nest_lock
   easyclimate.backend.omp_test_lock
   easyclimate.backend.omp_test_nest_lock
   easyclimate.backend.omp_get_wtime
   easyclimate.backend.omp_get_wtick
   easyclimate.backend.interplevel
   easyclimate.backend.vertcross
   easyclimate.backend.interpline
   easyclimate.backend.vinterp
   easyclimate.backend.xy_to_ll
   easyclimate.backend.ll_to_xy
   easyclimate.backend.xy_to_ll_proj
   easyclimate.backend.ll_to_xy_proj
   easyclimate.backend.viewitems
   easyclimate.backend.viewkeys
   easyclimate.backend.viewvalues
   easyclimate.backend.isstr
   easyclimate.backend.py2round
   easyclimate.backend.py3range
   easyclimate.backend.ucode
   easyclimate.backend.to_np
   easyclimate.backend.extract_global_attrs
   easyclimate.backend.is_standard_wrf_var
   easyclimate.backend.extract_dim
   easyclimate.backend.extract_vars
   easyclimate.backend.extract_times
   easyclimate.backend.combine_files
   easyclimate.backend.npbytes_to_str
   easyclimate.backend.is_moving_domain
   easyclimate.backend.is_staggered
   easyclimate.backend.get_left_indexes
   easyclimate.backend.iter_left_indexes
   easyclimate.backend.get_right_slices
   easyclimate.backend.get_proj_params
   easyclimate.backend.from_args
   easyclimate.backend.args_to_list
   easyclimate.backend.arg_location
   easyclimate.backend.psafilepath
   easyclimate.backend.get_id
   easyclimate.backend.from_var
   easyclimate.backend.combine_dims
   easyclimate.backend.either
   easyclimate.backend.get_iterable
   easyclimate.backend.IterWrapper
   easyclimate.backend.is_coordvar
   easyclimate.backend.latlon_coordvars
   easyclimate.backend.is_mapping
   easyclimate.backend.has_time_coord
   easyclimate.backend.is_multi_file
   easyclimate.backend.is_multi_time_req
   easyclimate.backend.get_coord_pairs
   easyclimate.backend.is_time_coord_var
   easyclimate.backend.geo_bounds
   easyclimate.backend.get_cartopy
   easyclimate.backend.get_basemap
   easyclimate.backend.get_pyngl
   easyclimate.backend.cartopy_xlim
   easyclimate.backend.cartopy_ylim
   easyclimate.backend.latlon_coords
   easyclimate.backend.ll_points
   easyclimate.backend.pairs_to_latlon
   easyclimate.backend.GeoBounds
   easyclimate.backend.NullGeoBounds
   easyclimate.backend.WrfProj
   easyclimate.backend.NullProjection
   easyclimate.backend.LambertConformal
   easyclimate.backend.Mercator
   easyclimate.backend.PolarStereographic
   easyclimate.backend.LatLon
   easyclimate.backend.RotatedLatLon
   easyclimate.backend.getproj
   easyclimate.backend.CoordPair
   easyclimate.backend.to_xy_coords
   easyclimate.backend.cache_item
   easyclimate.backend.get_cached_item
   easyclimate.backend.calc_wet_bulb_temperature_rs
   easyclimate.backend.calc_sphere_laplacian_numpy_rs
   easyclimate.backend.calc_sphere_laplacian_conservative_numpy_rs
   easyclimate.backend.calc_detrend_spatial_3d_rs
   easyclimate.backend.calc_detrend_spatial_3d_chunked_rs
   easyclimate.backend.calc_detrend_spatial_flexible_rs
   easyclimate.backend.interp1d_linear_core_rs
   easyclimate.backend.interp1d_linear_2d_rs
   easyclimate.backend.interp1d_linear_3d_rs
   easyclimate.backend.interp1d_linear_4d_rs
   easyclimate.backend.spharm
   easyclimate.backend.VectorWind
   easyclimate.backend._wet_bulb_temperature
   easyclimate.backend.getvar
   easyclimate.backend.interplevel
   easyclimate.backend.RUST_AVAILABLE


Functions
---------

.. autoapisummary::

   easyclimate.backend.wk_analysis
   easyclimate.backend._dummy_function


Module Contents
---------------

.. py:data:: CURRENT_PLATFORM

.. py:data:: spharm
   :value: None


.. py:data:: VectorWind
   :value: None


.. py:data:: aeronoskin
   :value: None


.. py:data:: aeroskin
   :value: None


.. py:data:: human_index_mod
   :value: None


.. py:data:: human_index_mod_old
   :value: None


.. py:data:: _ecl_redfit
   :value: None


.. py:data:: _ecl_redfit_x
   :value: None


.. py:data:: _vinth2p_dp
   :value: None


.. py:data:: _vinth2p_ecmwf
   :value: None


.. py:data:: _vintp2p_ecmwf
   :value: None


.. py:data:: _wet_bulb_temperature
   :value: None


.. py:data:: dvibeta
   :value: None


.. py:data:: dvrfidf
   :value: None


.. py:data:: ddvfidf
   :value: None


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
   :value: None


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
   :value: None


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


.. py:data:: calc_wet_bulb_temperature_rs
   :value: None


.. py:data:: calc_sphere_laplacian_numpy_rs
   :value: None


.. py:data:: calc_sphere_laplacian_conservative_numpy_rs
   :value: None


.. py:data:: calc_detrend_spatial_3d_rs
   :value: None


.. py:data:: calc_detrend_spatial_3d_chunked_rs
   :value: None


.. py:data:: calc_detrend_spatial_flexible_rs
   :value: None


.. py:data:: interp1d_linear_core_rs
   :value: None


.. py:data:: interp1d_linear_2d_rs
   :value: None


.. py:data:: interp1d_linear_3d_rs
   :value: None


.. py:data:: interp1d_linear_4d_rs
   :value: None


.. py:function:: wk_analysis(*args, **kwargs)

.. py:function:: _dummy_function(*args, **kwargs)

   A dummy function that raises an informative error when called.


.. py:data:: spharm

.. py:data:: VectorWind

.. py:data:: _wet_bulb_temperature

.. py:data:: getvar

.. py:data:: interplevel

.. py:data:: RUST_AVAILABLE
   :value: True


