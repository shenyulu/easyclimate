import platform
import warnings
import sys


# Detect the current operating system
CURRENT_PLATFORM = platform.system()

# --------------------------------------------
# Easyclimate-backend Initialize variables
# --------------------------------------------

# Initialize variables for backend functions with default values
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
dvibeta = None
dvrfidf = None
ddvfidf = None

# WRF-related variables with default values
xarray_enabled = False
disable_xarray = None
enable_xarray = None
cartopy_enabled = False
disable_cartopy = None
enable_cartopy = None
basemap_enabled = False
disable_basemap = None
enable_basemap = None
pyngl_enabled = False
enable_pyngl = None
disable_pyngl = None
set_cache_size = None
get_cache_size = None
omp_enabled = False
ALL_TIMES = None
Constants = None
ConversionFactors = None
ProjectionTypes = None
default_fill = None
OMP_SCHED_STATIC = None
OMP_SCHED_DYNAMIC = None
OMP_SCHED_GUIDED = None
OMP_SCHED_AUTO = None
destagger = None
getvar = None
xy = None
interp1d = None
interp2dxy = None
interpz3d = None
slp = None
tk = None
td = None
rh = None
uvmet = None
smooth2d = None
cape_2d = None
cape_3d = None
cloudfrac = None
ctt = None
dbz = None
srhel = None
udhel = None
avo = None
pvo = None
eth = None
wetbulb = None
tvirtual = None
omega = None
pw = None
DiagnosticError = None
omp_set_num_threads = None
omp_get_num_threads = None
omp_get_max_threads = None
omp_get_thread_num = None
omp_get_num_procs = None
omp_in_parallel = None
omp_set_dynamic = None
omp_get_dynamic = None
omp_set_nested = None
omp_get_nested = None
omp_set_schedule = None
omp_get_schedule = None
omp_get_thread_limit = None
omp_set_max_active_levels = None
omp_get_max_active_levels = None
omp_get_level = None
omp_get_ancestor_thread_num = None
omp_get_team_size = None
omp_get_active_level = None
omp_in_final = None
omp_init_lock = None
omp_init_nest_lock = None
omp_destroy_lock = None
omp_destroy_nest_lock = None
omp_set_lock = None
omp_set_nest_lock = None
omp_unset_lock = None
omp_unset_nest_lock = None
omp_test_lock = None
omp_test_nest_lock = None
omp_get_wtime = None
omp_get_wtick = None
interplevel = None
vertcross = None
interpline = None
vinterp = None
xy_to_ll = None
ll_to_xy = None
xy_to_ll_proj = None
ll_to_xy_proj = None
viewitems = None
viewkeys = None
viewvalues = None
isstr = None
py2round = None
py3range = None
ucode = None
to_np = None
extract_global_attrs = None
is_standard_wrf_var = None
extract_dim = None
extract_vars = None
extract_times = None
combine_files = None
npbytes_to_str = None
is_moving_domain = None
is_staggered = None
get_left_indexes = None
iter_left_indexes = None
get_right_slices = None
get_proj_params = None
from_args = None
args_to_list = None
arg_location = None
psafilepath = None
get_id = None
from_var = None
combine_dims = None
either = None
get_iterable = None
IterWrapper = None
is_coordvar = None
latlon_coordvars = None
is_mapping = None
has_time_coord = None
is_multi_file = None
is_multi_time_req = None
get_coord_pairs = None
is_time_coord_var = None
geo_bounds = None
get_cartopy = None
get_basemap = None
get_pyngl = None
cartopy_xlim = None
cartopy_ylim = None
latlon_coords = None
ll_points = None
pairs_to_latlon = None
GeoBounds = None
NullGeoBounds = None
WrfProj = None
NullProjection = None
LambertConformal = None
Mercator = None
PolarStereographic = None
LatLon = None
RotatedLatLon = None
getproj = None
CoordPair = None
to_xy_coords = None
cache_item = None
get_cached_item = None

# --------------------------------------------
# Easyclimate-rust Initialize variables
# --------------------------------------------
calc_wet_bulb_temperature_rs = None
calc_sphere_laplacian_numpy_rs = None
calc_sphere_laplacian_conservative_numpy_rs = None
calc_detrend_spatial_3d_rs = None
calc_detrend_spatial_3d_chunked_rs = None
calc_detrend_spatial_flexible_rs = None
interp1d_linear_core_rs = None
interp1d_linear_2d_rs = None
interp1d_linear_3d_rs = None
interp1d_linear_4d_rs = None

# --------------------------------------------
# Easyclimate-backend Import
# --------------------------------------------

# Cross platform module - these should work on all platforms
try:
    from easyclimate_backend.wk_spectra import wk_analysis, matsuno_plot
    from easyclimate_backend.wavelet.waveletFunctions import wave_signif, wavelet
except ImportError as e:
    warnings.warn(
        f"Failed to import cross-platform modules: {e}. Some core functionality will be disabled.",
        ImportWarning,
    )

    # Define fallbacks for cross-platform modules
    def wk_analysis(*args, **kwargs):
        raise ImportError("wk_analysis is not available due to import failure")

    def matsuno_plot(*args, **kwargs):
        raise ImportError("matsuno_plot is not available due to import failure")

    def wave_signif(*args, **kwargs):
        raise ImportError("wave_signif is not available due to import failure")

    def wavelet(*args, **kwargs):
        raise ImportError("wavelet is not available due to import failure")


# Attempt to import platform-specific functions only on Windows and Linux
if CURRENT_PLATFORM in ("Windows", "Linux"):
    try:
        # print(
        #     f"Attempting to import easyclimate-backend modules on {CURRENT_PLATFORM}...",
        #     file=sys.stderr,
        # )

        # Import basic modules
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

        from easyclimate_backend.vibeta._vibeta_dp import dvibeta
        from easyclimate_backend.rvdv._rvdv import ddvfidf, dvrfidf

        from easyclimate_backend.wet_bulb import _wet_bulb_temperature

        # print(
        #     "Successfully imported basic easyclimate-backend modules", file=sys.stderr
        # )

    except ImportError as e:
        warnings.warn(
            f"Failed to import basic easyclimate-backend modules: {e}. Some functionality will be disabled.",
            ImportWarning,
        )
        # Record detailed error information to stderr
        print(f"Detailed import error for basic modules: {e}", file=sys.stderr)

    # Attempt to import WRF-related modules - handle separately, as WRF may have additional dependencies
    try:
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

        # print("Successfully imported WRF modules", file=sys.stderr)

    except ImportError as wrf_e:
        warnings.warn(
            f"Failed to import WRF-related modules: {wrf_e}. WRF functionality will be disabled.",
            ImportWarning,
        )
        print(f"Detailed WRF import error: {wrf_e}", file=sys.stderr)

else:
    warnings.warn(
        f"easyclimate-backend modules is not supported on {CURRENT_PLATFORM}. Related functionality will be disabled."
    )


# Define some basic fallback functions
def _dummy_function(*args, **kwargs):
    """A dummy function that raises an informative error when called."""
    raise ImportError(
        "This function is not available because easyclimate-backend modules failed to import. "
        f"Check that easyclimate-backend is properly installed on {CURRENT_PLATFORM}."
    )


# If the import of some key functions fails, replace them with fallback functions.
if spharm is None:
    spharm = _dummy_function

if VectorWind is None:
    VectorWind = _dummy_function

if _wet_bulb_temperature is None:
    _wet_bulb_temperature = _dummy_function

# For the WRF function, if it is None, set it to the fallback function.
if getvar is None:
    getvar = _dummy_function

if interplevel is None:
    interplevel = _dummy_function

# Print the summary of the import status
if CURRENT_PLATFORM in ("Windows", "Linux"):
    if all(func is not None for func in [spharm, VectorWind, getvar]):
        # print("easyclimate-backend modules imported successfully", file=sys.stderr)
        pass
    else:
        print("Some easyclimate-backend modules failed to import", file=sys.stderr)

# --------------------------------------------
# Easyclimate-rust Import
# --------------------------------------------

# Attempt to import platform-specific functions only on Windows and Linux
if CURRENT_PLATFORM in ("Windows", "Linux"):
    try:
        # print(
        #     f"Attempting to import easyclimate-rust modules on {CURRENT_PLATFORM}...",
        #     file=sys.stderr,
        # )

        # Import basic modules

        # wet_bulb module
        from easyclimate_rust._easyclimate_rust import (
            calc_wet_bulb_temperature as calc_wet_bulb_temperature_rs,
        )

        # sphere_laplacian module (NumPy interface)
        from easyclimate_rust._easyclimate_rust import (
            calc_sphere_laplacian_numpy as calc_sphere_laplacian_numpy_rs,
        )
        from easyclimate_rust._easyclimate_rust import (
            calc_sphere_laplacian_conservative_numpy as calc_sphere_laplacian_conservative_numpy_rs,
        )

        # detrend_spatial module - high-performance spatial detrending
        from easyclimate_rust._easyclimate_rust import (
            calc_detrend_spatial_3d as calc_detrend_spatial_3d_rs,
        )
        from easyclimate_rust._easyclimate_rust import (
            calc_detrend_spatial_3d_chunked as calc_detrend_spatial_3d_chunked_rs,
        )
        from easyclimate_rust._easyclimate_rust import (
            calc_detrend_spatial_flexible as calc_detrend_spatial_flexible_rs,
        )

        # interp1d
        from easyclimate_rust._easyclimate_rust import (
            interp1d_linear_core as interp1d_linear_core_rs,
        )
        from easyclimate_rust._easyclimate_rust import (
            interp1d_linear_2d as interp1d_linear_2d_rs,
        )
        from easyclimate_rust._easyclimate_rust import (
            interp1d_linear_3d as interp1d_linear_3d_rs,
        )
        from easyclimate_rust._easyclimate_rust import (
            interp1d_linear_4d as interp1d_linear_4d_rs,
        )

        # print(
        #     "Successfully imported basic easyclimate-rust modules", file=sys.stderr
        # )

        RUST_AVAILABLE = True

    except ImportError as e:
        warnings.warn(
            f"Failed to import basic easyclimate-rust modules: {e}. Some functionality will be disabled.",
            ImportWarning,
        )
        # Record detailed error information to stderr
        print(f"Detailed import error for basic modules: {e}", file=sys.stderr)

        RUST_AVAILABLE = False
else:
    warnings.warn(
        f"easyclimate-rust modules is not supported on {CURRENT_PLATFORM}. Related functionality will be disabled."
    )
