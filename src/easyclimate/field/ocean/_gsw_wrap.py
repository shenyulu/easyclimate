"""
Small xarray-friendly wrapper around GSW-Python.

This module provides the xarray-aware subset of GSW-Python used by the ocean
field calculations.
"""

from __future__ import annotations

import gsw
import xarray as xr

__all__ = [
    "CT_from_t",
    "Nsquared",
    "SA_from_SP",
    "p_from_z",
    "pot_rho_t_exact",
    "rho",
]


def _apply(func, *args):
    return xr.apply_ufunc(
        func,
        *args,
        dask="parallelized",
        keep_attrs=True,
        output_dtypes=[float],
    )


def p_from_z(z, lat):
    """Calculate sea pressure from height using GSW-Python."""
    return _apply(gsw.p_from_z, z, lat)


def SA_from_SP(SP, p, lon, lat):
    """Calculate Absolute Salinity from Practical Salinity."""
    return _apply(gsw.SA_from_SP, SP, p, lon, lat)


def CT_from_t(SA, t, p):
    """Calculate Conservative Temperature from in-situ temperature."""
    return _apply(gsw.CT_from_t, SA, t, p)


def Nsquared(SA, CT, p, lat=None, axis=0):
    """Calculate buoyancy frequency squared and midpoint pressure."""
    return gsw.Nsquared(SA, CT, p, lat=lat, axis=axis)


def pot_rho_t_exact(SA, t, p, p_ref):
    """Calculate potential density from in-situ temperature."""
    return _apply(gsw.pot_rho_t_exact, SA, t, p, p_ref)


def rho(SA, CT, p):
    """Calculate in-situ density."""
    return _apply(gsw.rho, SA, CT, p)
