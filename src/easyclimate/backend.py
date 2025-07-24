import platform
import warnings


# Detect the current operating system
CURRENT_PLATFORM = platform.system()


# Initialize variables for backend functions
spharm = None
VectorWind = None
aeronoskin = None
aeroskin = None
human_index_mod = None
human_index_mod_old = None
_ecl_redfit = None
_ecl_redfit_x = None
_vinth2p_dp = None
_vinth2p_ecmwf = None
_vintp2p_ecmwf = None
_wet_bulb_temperature = None


# Cross platform module
from easyclimate_backend.wk_spectra import wk_analysis, matsuno_plot
from easyclimate_backend.wavelet.waveletFunctions import wave_signif, wavelet


# Attempt to import platform-specific functions only on Windows and Linux
if CURRENT_PLATFORM in ("Windows", "Linux"):
    try:
        from easyclimate_backend.pyspharm import spharm
        from easyclimate_backend.windspharm.xarray import VectorWind
        from easyclimate_backend.aerobulk import mod_aerobulk_wrap_noskin as aeronoskin
        from easyclimate_backend.aerobulk import mod_aerobulk_wrap_skin as aeroskin
        from easyclimate_backend.heat_stress import human_index_mod, human_index_mod_old
        from easyclimate_backend.redfit import _ecl_redfit
        from easyclimate_backend.redfit import _ecl_redfit_x
        from easyclimate_backend.vinth2p._vinth2p_dp import vinth2p as _vinth2p_dp
        from easyclimate_backend.vinth2p._vinth2p_ecmwf import (
            vinth2pecmwf as _vinth2p_ecmwf,
        )
        from easyclimate_backend.vinth2p._vintp2p_ecmwf import (
            vintp2pecmwf as _vintp2p_ecmwf,
        )
        from easyclimate_backend.wet_bulb import _wet_bulb_temperature
        from easyclimate_backend.wrf import (
            xarray_enabled,
            disable_xarray,
            enable_xarray,
            cartopy_enabled,
            disable_cartopy,
            enable_cartopy,
            basemap_enabled,
            disable_basemap,
            enable_basemap,
            pyngl_enabled,
            enable_pyngl,
            disable_pyngl,
            set_cache_size,
            get_cache_size,
            omp_enabled,
            ALL_TIMES,
            Constants,
            ConversionFactors,
            ProjectionTypes,
            default_fill,
            OMP_SCHED_STATIC,
            OMP_SCHED_DYNAMIC,
            OMP_SCHED_GUIDED,
            OMP_SCHED_AUTO,
            destagger,
            getvar,
            xy,
            interp1d,
            interp2dxy,
            interpz3d,
            slp,
            tk,
            td,
            rh,
            uvmet,
            smooth2d,
            cape_2d,
            cape_3d,
            cloudfrac,
            ctt,
            dbz,
            srhel,
            udhel,
            avo,
            pvo,
            eth,
            wetbulb,
            tvirtual,
            omega,
            pw,
            DiagnosticError,
            omp_set_num_threads,
            omp_get_num_threads,
            omp_get_max_threads,
            omp_get_thread_num,
            omp_get_num_procs,
            omp_in_parallel,
            omp_set_dynamic,
            omp_get_dynamic,
            omp_set_nested,
            omp_get_nested,
            omp_set_schedule,
            omp_get_schedule,
            omp_get_thread_limit,
            omp_set_max_active_levels,
            omp_get_max_active_levels,
            omp_get_level,
            omp_get_ancestor_thread_num,
            omp_get_team_size,
            omp_get_active_level,
            omp_in_final,
            omp_init_lock,
            omp_init_nest_lock,
            omp_destroy_lock,
            omp_destroy_nest_lock,
            omp_set_lock,
            omp_set_nest_lock,
            omp_unset_lock,
            omp_unset_nest_lock,
            omp_test_lock,
            omp_test_nest_lock,
            omp_get_wtime,
            omp_get_wtick,
            interplevel,
            vertcross,
            interpline,
            vinterp,
            xy_to_ll,
            ll_to_xy,
            xy_to_ll_proj,
            ll_to_xy_proj,
            viewitems,
            viewkeys,
            viewvalues,
            isstr,
            py2round,
            py3range,
            ucode,
            to_np,
            extract_global_attrs,
            is_standard_wrf_var,
            extract_dim,
            extract_vars,
            extract_times,
            combine_files,
            npbytes_to_str,
            is_moving_domain,
            is_staggered,
            get_left_indexes,
            iter_left_indexes,
            get_right_slices,
            get_proj_params,
            from_args,
            args_to_list,
            arg_location,
            psafilepath,
            get_id,
            from_var,
            combine_dims,
            either,
            get_iterable,
            IterWrapper,
            is_coordvar,
            latlon_coordvars,
            is_mapping,
            has_time_coord,
            is_multi_file,
            is_multi_time_req,
            get_coord_pairs,
            is_time_coord_var,
            geo_bounds,
            get_cartopy,
            get_basemap,
            get_pyngl,
            cartopy_xlim,
            cartopy_ylim,
            latlon_coords,
            ll_points,
            pairs_to_latlon,
            GeoBounds,
            NullGeoBounds,
            WrfProj,
            NullProjection,
            LambertConformal,
            Mercator,
            PolarStereographic,
            LatLon,
            RotatedLatLon,
            getproj,
            CoordPair,
            to_xy_coords,
            cache_item,
            get_cached_item,
        )

    except ImportError as e:
        warnings.warn(
            f"Failed to import easyclimate-backend related modules: {e}", ImportWarning
        )
else:
    warnings.warn(
        f"easyclimate-backend modules is not supported on {CURRENT_PLATFORM}. Related functionality will be disabled."
    )
