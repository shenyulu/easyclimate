"""
Geostrophic Wind
"""

from __future__ import annotations

from .utility import (
    transfer_inf2nan,
)
import xarray as xr
from typing import Literal

__all__ = [
    "calc_geostrophic_wind",
    "calc_geostrophic_wind_vorticity",
]


def calc_geostrophic_wind(
    z_data: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    omega: float = 7.292e-5,
    g: float = 9.8,
    R: float = 6371200.0,
) -> xr.DataArray:
    """
    Calculate the geostrophic wind.

    .. math::
        u_g = - \\frac{g}{f} \\frac{\\partial H}{\\partial y}

    .. math::
        v_g = \\frac{g}{f} \\frac{\\partial H}{\\partial x}

    Parameters
    ----------
    z_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric geopotential height.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    omega: :py:class:`float <float>`, default: `7.292e-5`.
        The angular speed of the earth.
    g: :py:class:`float <float>`, default: `9.8`.
        The acceleration of gravity.
    R: :py:class:`float <float>`, default: `6370000`.
        Radius of the Earth.

    Returns
    -------
    The geostrophic wind term. (:py:class:`xarray.DataArray<xarray.DataArray>`).
        - ug
        - vg

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    from ..physics.geo.coriolis import get_coriolis_parameter
    from .diff import calc_dx_gradient, calc_dy_gradient

    lat_array = z_data[lat_dim]
    f = get_coriolis_parameter(lat_array, omega=omega)

    dHdy = calc_dy_gradient(z_data, lat_dim=lat_dim, R=R)
    dHdx = calc_dx_gradient(z_data, lon_dim=lon_dim, lat_dim=lat_dim, R=R)

    ug = -(g / f) * dHdy
    vg = (g / f) * dHdx

    uvg = xr.Dataset(data_vars={"ug": transfer_inf2nan(ug), "vg": transfer_inf2nan(vg)})

    return uvg


def calc_geostrophic_wind_vorticity(
    z_data: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    spherical_coord: bool = True,
    omega: float = 7.292e-5,
    g: float = 9.8,
    R: float = 6371200.0,
    cyclic_boundary_setting: Literal["nan", "cyclic", "cyclic+diff", "diff"] = "nan",
    method: Literal["raw", "ncl", "rust"] = "ncl",
) -> xr.DataArray:
    """
    Calculate the geostrophic vorticity.

    Rectangular coordinates

    .. math::
        \\zeta_g = \\frac{\\partial v_g}{\\partial x} - \\frac{\\partial u_g}{\\partial y}

    Spherical coordinates

    .. math::
        \\zeta_g = \\frac{\\partial v_g}{\\partial x} - \\frac{\\partial u_g}{\\partial y} + \\frac{u_g}{R} \\tan \\varphi

    Parameters
    ----------
    z_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric geopotential height.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    spherical_coord: :py:class:`bool<bool>`, default: `True`.
        Whether or not to compute the horizontal Laplace term in spherical coordinates.
    omega: :py:class:`float <float>`, default: `7.292e-5`.
        The angular speed of the earth.
    g: :py:class:`float <float>`, default: `9.8`.
        The acceleration of gravity.
    R: :py:class:`float <float>`, default: `6370000`.
        Radius of the Earth.
    cyclic_boundary_setting: {"nan", "cyclic", "cyclic+diff", "diff"}, default: `nan`.
        A scalar integer equal to the boundary condition option:

        - ``nan``: Boundary points are set to the missing value.
        - ``cyclic``: The u and v arrays are cyclic in longitude. (The arrays should **NOT** include the cyclic point.) The upper and lower boundaries will be set to missing.
        - ``cyclic+diff``: Boundary points are estimated using one-sided difference schemes normal to the boundary.
        - ``diff``: The u and v arrays are cyclic in longitude. (The arrays should **NOT** include the cyclic points.) The upper and lower boundaries are estimated using a one-sided difference scheme normal to the boundary.

        .. note::
            The parameter is applicable only when ``method = ncl`` or ``method = rust``.

    method: {"raw", "ncl", rust}, default: `ncl`.
        The method to calculate horizontal divergence term. Optional values are ``raw``, ``ncl`` or ``rust``.

    Returns
    -------
    The geostrophic vorticity term. (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    from .rvdv import calc_vorticity, calc_vorticity_ncl, calc_vorticity_rs

    geostrophic_wind = calc_geostrophic_wind(
        z_data, lon_dim=lon_dim, lat_dim=lat_dim, omega=omega, g=g, R=R
    )
    ug, vg = geostrophic_wind["ug"], geostrophic_wind["vg"]

    if method in ["raw"]:
        vor_g = calc_vorticity(
            ug,
            vg,
            spherical_coord=spherical_coord,
            lon_dim=lon_dim,
            lat_dim=lat_dim,
            R=R,
        )
    elif method in ["ncl"]:
        vor_g = calc_vorticity_ncl(
            ug,
            vg,
            lon_dim=lon_dim,
            lat_dim=lat_dim,
            cyclic_boundary_setting=cyclic_boundary_setting,
        )
    elif method == "rust":
        vor_g = calc_vorticity_rs(
            ug,
            vg,
            lon_dim=lon_dim,
            lat_dim=lat_dim,
            R=R,
            cyclic_boundary_setting=cyclic_boundary_setting,
        )
    else:
        raise ValueError(
            f"Unsupported method={method!r}. Expected one of 'raw', 'ncl' or 'rust'."
        )

    return vor_g
