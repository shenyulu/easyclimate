"""
Advection
"""

from __future__ import annotations

import xarray as xr
from typing import Literal
from .utility import compare_multi_dataarray_coordinate

__all__ = [
    "calc_u_advection",
    "calc_v_advection",
    "calc_p_advection",
]


def calc_u_advection(
    u_data: xr.DataArray,
    data_input: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    min_dx: float = 1.0,
    edge_order: int = 2,
    R: float = 6371200.0,
) -> xr.DataArray:
    """
    Calculate zonal advection of the input variable at each vertical level.

    .. math::

        -u \\frac{\\partial F}{\\partial x}

    Parameters
    ----------
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Input variable :math:`F` to be advected.
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
    The zonal advection of :math:`F`. (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        geographic_finite_difference
    """
    from .diff import calc_dx_gradient

    dFdx = calc_dx_gradient(
        data_input,
        lon_dim=lon_dim,
        lat_dim=lat_dim,
        min_dx=min_dx,
        edge_order=edge_order,
        R=R,
    )
    u_adv = (-1) * u_data * dFdx
    input_name = data_input.name or "data"
    input_units = data_input.attrs.get("units", "dimensionless")
    u_adv.attrs = {
        "long_name": f"Zonal advection of {input_name}",
        "units": f"{input_units} per second",
        "standard_name": "zonal_advection_of_data",
    }
    u_adv.name = f"{input_name}_u_advection"
    return u_adv


def calc_v_advection(
    v_data: xr.DataArray,
    data_input: xr.DataArray,
    lat_dim: str = "lat",
    min_dy: float = 1.0,
    edge_order: int = 2,
    R: float = 6371200.0,
) -> xr.DataArray:
    """
    Calculate meridional advection of the input variable at each vertical level.

    .. math::

        -v \\frac{\\partial F}{\\partial y}

    Parameters
    ----------
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Input variable :math:`F` to be advected.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

    Returns
    -------
    The meridional advection of :math:`F`. (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        geographic_finite_difference
    """
    from .diff import calc_dy_gradient

    dFdy = calc_dy_gradient(
        data_input, lat_dim=lat_dim, min_dy=min_dy, edge_order=edge_order, R=R
    )
    v_adv = (-1) * v_data * dFdy
    input_name = data_input.name or "data"
    input_units = data_input.attrs.get("units", "dimensionless")
    v_adv.attrs = {
        "long_name": f"Meridional advection of {input_name}",
        "units": f"{input_units} per second",
        "standard_name": "meridional_advection_of_data",
    }
    v_adv.name = f"{input_name}_v_advection"
    return v_adv


def calc_p_advection(
    omega_data: xr.DataArray,
    data_input: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
) -> xr.DataArray:
    """
    Calculate vertical advection of the input variable at each vertical level.

    .. math::

        -\\omega \\frac{\\partial F}{\\partial p}

    Parameters
    ----------
    omega_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The vertical velocity data (:math:`\\frac{\\mathrm{d} p}{\\mathrm{d} t}`).
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Input variable :math:`F` to be advected.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.

    Returns
    -------
    The vertical advection of :math:`F`. (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    from .diff import calc_p_gradient

    compare_multi_dataarray_coordinate([omega_data, data_input])

    dFdp = calc_p_gradient(
        data_input, vertical_dim=vertical_dim, vertical_dim_units=vertical_dim_units
    )
    p_adv = -omega_data * dFdp
    input_name = data_input.name or "data"
    input_units = data_input.attrs.get("units", "dimensionless")
    p_adv.attrs = {
        "long_name": f"Pressure-coordinate advection of {input_name}",
        "units": f"{input_units} per second",
        "standard_name": "pressure_coordinate_advection_of_data",
    }
    p_adv.name = f"{input_name}_p_advection"
    return p_adv
