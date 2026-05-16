import pytest

from easyclimate import backend


def _available(name):
    return getattr(backend, name, None) is not None


def _spharm_available():
    return hasattr(getattr(backend, "spharm", None), "Spharmt")


def _vectorwind_available():
    value = getattr(backend, "VectorWind", None)
    return value is not None and getattr(value, "__name__", "") != "_dummy_function"


def _wrf_available():
    value = getattr(backend, "getvar", None)
    return value is not None and getattr(value, "__name__", "") != "_dummy_function"


BACKEND_FILE_REQUIREMENTS = {
    "test_core_spharm.py": (
        _spharm_available,
        "easyclimate-backend pyspharm is unavailable",
    ),
    "test_core_windspharm.py": (
        _vectorwind_available,
        "easyclimate-backend windspharm is unavailable",
    ),
    "test_field_boundary_layer_aerobulk.py": (
        lambda: _available("aeronoskin") and _available("aeroskin"),
        "easyclimate-backend aerobulk is unavailable",
    ),
    "test_field_heat_stress_humanindexmod.py": (
        lambda: _available("human_index_mod"),
        "easyclimate-backend heat_stress is unavailable",
    ),
    "test_filter_redfitpy.py": (
        lambda: _available("_ecl_redfit") and _available("_ecl_redfit_x"),
        "easyclimate-backend redfit is unavailable",
    ),
    "test_interp_vinth2p.py": (
        lambda: _available("_vinth2p_dp") and _available("_vinth2p_ecmwf"),
        "easyclimate-backend vinth2p is unavailable",
    ),
    "test_interp_vinth2p1.py": (
        lambda: _available("_vintp2p_ecmwf"),
        "easyclimate-backend vintp2p is unavailable",
    ),
    "test_wrf_interface.py": (_wrf_available, "easyclimate-backend WRF is unavailable"),
}

CORE_DIFF_NCL_TESTS = {
    "test_calc_top2surface_integral1",
    "test_calc_top2surface_integral2",
    "test_calc_top2surface_integral3",
    "test_calc_top2surface_integral4",
    "test_calc_top2surface_integral5",
    "test_calc_top2surface_integral6",
    "test_calc_top2surface_average1",
    "test_calc_top2surface_integral_mass_weighted1",
    "test_calc_top2surface_integral_normalize_dispatch1",
    "test_calc_water_flux_top2surface_integral1",
    "test_calc_water_flux_top2surface_integral2",
    "test_calc_water_flux_top2surface_integral3",
}

CORE_DIFF_RVDV_NCL_TESTS = {
    "test_calc_divergence2",
    "test_calc_vorticity2",
    "test_calc_geostrophic_wind_vorticity2_ncl",
    "test_calc_geostrophic_wind_vorticity3_ncl",
    "test_calc_geostrophic_wind_vorticity4_ncl",
    "test_calc_geostrophic_wind_vorticity5_ncl",
    "test_calc_divergence_watervaporflux1",
    "test_calc_divergence_watervaporflux2_ncl",
}

CORE_DIFF_VIBETA_AND_RVDV_NCL_TESTS = {
    "test_calc_divergence_watervaporflux_top2surface_integral1",
    "test_calc_divergence_watervaporflux_top2surface_integral2",
}


def pytest_collection_modifyitems(config, items):
    vibeta_available = _available("dvibeta_ncl")
    rvdv_ncl_available = _available("ddvfidf_ncl") and _available("dvrfidf_ncl")

    for item in items:
        file_name = item.path.name

        requirement = BACKEND_FILE_REQUIREMENTS.get(file_name)
        if requirement is not None:
            predicate, reason = requirement
            if not predicate():
                item.add_marker(pytest.mark.skip(reason=reason))
            continue

        if file_name != "test_core_diff.py":
            continue

        if item.name in CORE_DIFF_NCL_TESTS and not vibeta_available:
            item.add_marker(
                pytest.mark.skip(reason="easyclimate-backend vibeta is unavailable")
            )
        elif item.name in CORE_DIFF_RVDV_NCL_TESTS and not rvdv_ncl_available:
            item.add_marker(
                pytest.mark.skip(reason="easyclimate-backend rvdv is unavailable")
            )
        elif item.name in CORE_DIFF_VIBETA_AND_RVDV_NCL_TESTS and not (
            vibeta_available and rvdv_ncl_available
        ):
            item.add_marker(
                pytest.mark.skip(
                    reason="easyclimate-backend vibeta or rvdv is unavailable"
                )
            )
