"""
Advection
"""

from __future__ import annotations

import xarray as xr
from typing import Literal

__all__ = [
    "calc_u_advection",
    "calc_v_advection",
    "calc_p_advection",
]


def calc_u_advection(
    u_data: xr.DataArray,
    temper_data: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    min_dx: float = 1.0,
    edge_order: int = 2,
    R: float = 6371200.0,
) -> xr.DataArray:
    """
    Calculate zonal temperature advection at each vertical level.

    .. math::

        -u \\frac{\\partial T}{\\partial x}

    Parameters
    ----------
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Air temperature.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    min_dx: :py:class:`float <float>`, default: `1.0`.
        The minimum acceptable value of `dx`, below which parts will set `nan` to avoid large computational errors.
        The unit is m. You can set it to a negative value in order to remove this benefit.
    edge_order: {1, 2}, optional
        Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
    R: :py:class:`float <float>`, default: `6370000`.
        Radius of the Earth.

    Returns
    -------
    The zonal temperature advection. (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        geographic_finite_difference
    """
    from .diff import calc_dx_gradient

    dTdx = calc_dx_gradient(
        temper_data,
        lon_dim=lon_dim,
        lat_dim=lat_dim,
        min_dx=min_dx,
        edge_order=edge_order,
        R=R,
    )
    u_adv = (-1) * u_data * dTdx
    return u_adv


def calc_v_advection(
    v_data: xr.DataArray,
    temper_data: xr.DataArray,
    lat_dim: str = "lat",
    min_dy: float = 1.0,
    edge_order: int = 2,
    R: float = 6371200.0,
) -> xr.DataArray:
    """
    Calculate meridional temperature advection at each vertical level.

    .. math::

        -v \\frac{\\partial T}{\\partial y}

    Parameters
    ----------
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Air temperature.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

    Returns
    -------
    The meridional temperature advection. (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        geographic_finite_difference
    """
    from .diff import calc_dy_gradient

    dTdy = calc_dy_gradient(
        temper_data, lat_dim=lat_dim, min_dy=min_dy, edge_order=edge_order, R=R
    )
    v_adv = (-1) * v_data * dTdy
    return v_adv


def calc_p_advection(
    omega_data: xr.DataArray,
    temper_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
) -> xr.DataArray:
    """
    Calculate vertical temperature transport at each vertical level.

    .. math::

        -\\omega \\frac{\\partial T}{\\partial p}

    Parameters
    ----------
    omega: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The vertical velocity data (:math:`\\frac{\\mathrm{d} p}{\\mathrm{d} t}`).
    temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Air temperature.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.

    Returns
    -------
    The vertical temperature transport. (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    from .diff import calc_p_gradient

    dTdp = calc_p_gradient(
        temper_data, vertical_dim=vertical_dim, vertical_dim_units=vertical_dim_units
    )
    p_adv = -omega_data * dTdp
    return p_adv
