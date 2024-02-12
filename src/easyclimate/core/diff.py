"""
The calculation of geographic finite difference
"""

from __future__ import annotations

import numpy as np
from .utility import (
    find_dims_axis,
    transfer_deg2rad,
    transfer_inf2nan,
    transfer_units_coeff,
    transfer_data_units,
    generate_dataset_dispatcher,
)
import xarray as xr
import dask

__all__ = [
    "calc_gradient",
    "calc_lon_gradient",
    "calc_lat_gradient",
    "calc_lon_laplacian",
    "calc_lat_laplacian",
    "calc_lon_lat_mixed_derivatives",
    "calc_p_gradient",
    "calc_time_gradient",
    "calc_delta_pressure",
    "calc_p_integral",
    "calc_top2surface_integral",
    "calc_laplacian",
    "calc_divergence",
    "calc_vorticity",
    "calc_geostrophic_wind",
    "calc_geostrophic_wind_vorticity",
    "calc_horizontal_water_flux",
    "calc_vertical_water_flux",
    "calc_water_flux_top2surface_integral",
    "calc_divergence_watervaporflux",
    "calc_divergence_watervaporflux_top2surface_integral",
    "calc_u_advection",
    "calc_v_advection",
    "calc_p_advection",
]


@generate_dataset_dispatcher
def calc_gradient(
    data_input: xr.DataArray | xr.Dataset,
    dim: str,
    varargs: int = 1,
    edge_order: int = 2,
) -> xr.DataArray | xr.Dataset:
    """
    Compute the gradient along the coordinate `dim` direction.

    The gradient is computed using **second order accurate central differences** in the interior points
    and either first or second order accurate one-sides (forward or backwards) differences at the boundaries.
    The returned gradient hence has the same shape as the input array.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
         The spatio-temporal data to be calculated.
    dim : :py:class:`str <str>`.
        Dimension(s) over which to apply gradient. By default gradient is applied over the `time` dimension.
    varargs: :py:class:`list <list>` of scalar or array, optional
        Spacing between f values. Default unitary spacing for all dimensions. Spacing can be specified using:

        1. Single scalar to specify a sample distance for all dimensions.
        2. N scalars to specify a constant sample distance for each dimension. i.e. :math:`\\mathrm{d}x, \\mathrm{d}y, \\mathrm{d}z, ...`
        3. N arrays to specify the coordinates of the values along each dimension of F.
           The length of the array must match the size of the corresponding dimension.
        4. Any combination of N scalars/arrays with the meaning of 2. and 3.

    edge_order: {1, 2}, optional
        Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.

    Returns
    -------
    The gradient along the coordinate `dim` direction (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. seealso::
        :py:func:`numpy.gradient <numpy:numpy.gradient>`
    """

    def _calc_gradient(data_input, dim, varargs, edge_order) -> xr.DataArray:
        dim_index = find_dims_axis(data_input, dim=dim)
        # `data_input.data` make dask to process Dask.array
        data_gradient = np.gradient(
            data_input.data, varargs, axis=dim_index, edge_order=edge_order
        )
        return data_input.copy(data=data_gradient, deep=True)

    return _calc_gradient(data_input, dim, varargs, edge_order)


def calc_lon_gradient(
    data_input: xr.DataArray | xr.Dataset,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    min_dx: float = 1.0,
    edge_order: int = 2,
    R: float = 6370000,
) -> xr.DataArray | xr.Dataset:
    """
    Calculate the gradient along the longitude.

    .. math::
        \\frac{\\partial F}{\\partial x} = \\frac{1}{R \\cos\\varphi} \\cdot \\frac{\\partial F}{\\partial \\lambda}

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
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
    The gradient along the longitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. seealso::
        :py:func:`calc_gradient <calc_gradient>`
    """
    # Set to `float64` for more accurate results in trigonometric calculations.
    lon_array = data_input[lon_dim].astype("float64")
    lat_array = data_input[lat_dim].astype("float64")
    dlon = transfer_deg2rad(calc_gradient(lon_array, dim=lon_dim))
    coslat = np.cos(transfer_deg2rad(lat_array))

    dx = R * coslat * dlon
    dx = dx.where(np.abs(dx) >= min_dx)

    dFdx_raw = calc_gradient(data_input, dim=lon_dim, edge_order=edge_order)
    dFdx = dFdx_raw / dx
    return dFdx


def calc_lat_gradient(
    data_input: xr.DataArray | xr.Dataset,
    lat_dim: str = "lat",
    min_dy: float = 1.0,
    edge_order: int = 2,
    R: float = 6370000,
) -> xr.DataArray | xr.Dataset:
    """
    Calculate the gradient along the latitude.

    .. math::
        \\frac{\\partial F}{\\partial y} = \\frac{1}{R} \\cdot \\frac{\\partial F}{\\partial \\varphi}

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    min_dy: :py:class:`float <float>`, default: `1.0`.
        The minimum acceptable value of `dy`, below which parts will set `nan` to avoid large computational errors.
        The unit is m. You can set it to a negative value in order to remove this benefit.
    edge_order: {1, 2}, optional
        Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
    R: :py:class:`float <float>`, default: `6370000`.
        Radius of the Earth.

    Returns
    -------
    The gradient along the latitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. seealso::
        :py:func:`calc_gradient <calc_gradient>`
    """
    # Set to `float64` for more accurate results in trigonometric calculations.
    lat_array = data_input[lat_dim].astype("float64")
    dlat = transfer_deg2rad(calc_gradient(lat_array, dim=lat_dim))
    dy = R * dlat
    dy = dy.where(np.abs(dy) >= min_dy)

    dFdy_raw = calc_gradient(data_input, dim=lat_dim, edge_order=edge_order)
    dFdy = dFdy_raw / dy
    return dFdy


def calc_lon_laplacian(
    data_input: xr.DataArray | xr.Dataset,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    min_dx2: float = 1e9,
    edge_order: int = 2,
    R: float = 6370000,
) -> xr.DataArray | xr.Dataset:
    """
    Calculation of the second-order partial derivative term (Laplace term) along longitude.

    .. math::
        \\frac{\\partial^2 F}{\\partial x^2} = \\frac{1}{(R \\cos\\varphi)^2} \\cdot \\frac{\\partial^2 F}{\\partial \\lambda^2}

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    min_dx2: :py:class:`float <float>`, default: `1e9`.
        The minimum acceptable value of :math:`(\\mathrm{d}x)^2`, below which parts will set `nan` to avoid large computational errors.
        The unit is m. You can set it to a negative value in order to remove this benefit.
    edge_order: {1, 2}, optional
        Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
    R: :py:class:`float <float>`, default: `6370000`.
        Radius of the Earth.

    Returns
    -------
    The second-order partial derivative term (Laplace term) along longitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. seealso::
        :py:func:`calc_gradient <calc_gradient>`
    """
    # Set to `float64` for more accurate results in trigonometric calculations.
    lon_array = data_input[lon_dim].astype("float64")
    lat_array = data_input[lat_dim].astype("float64")

    dlon = transfer_deg2rad(calc_gradient(lon_array, dim=lon_dim))
    # The purpose of this calculation is to compute the trigonometric function only once and reduce the computational error.
    # cos^2(x)= ( 1+cos(2x) ) /2
    cos2lat = (np.cos(2.0 * transfer_deg2rad(lat_array)) + 1.0) / 2.0
    dx2 = (dlon * R) ** 2 * cos2lat
    dx2 = dx2.where(dx2 >= min_dx2)

    dFdx_raw = calc_gradient(data_input, dim=lon_dim, edge_order=edge_order)
    d2Fdx2_raw = calc_gradient(dFdx_raw, dim=lon_dim, edge_order=edge_order)
    d2Fdx2 = d2Fdx2_raw / dx2
    return d2Fdx2


def calc_lat_laplacian(
    data_input: xr.DataArray | xr.Dataset,
    lat_dim: str = "lat",
    min_dy2: float = 1.0,
    edge_order: int = 2,
    R: float = 6370000,
) -> xr.DataArray | xr.Dataset:
    """
    Calculation of the second-order partial derivative term (Laplace term) along latitude.

    .. math::
        \\frac{\\partial^2 F}{\\partial y^2} = \\frac{1}{R^2} \\cdot \\frac{\\partial^2 F}{\\partial \\varphi^2}

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    min_dy2: :py:class:`float <float>`, default: `1.0`.
        The minimum acceptable value of :math:`(\\mathrm{d}y)^2`, below which parts will set `nan` to avoid large computational errors.
        The unit is m. You can set it to a negative value in order to remove this benefit.
    edge_order: {1, 2}, optional
        Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
    R: :py:class:`float <float>`, default: `6370000`.
        Radius of the Earth.

    Returns
    -------
    The second-order partial derivative term (Laplace term) along latitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. seealso::
        :py:func:`calc_gradient <calc_gradient>`
    """
    lat_array = data_input[lat_dim].astype("float64")

    dlat = transfer_deg2rad(calc_gradient(lat_array, dim=lat_dim))
    dy2 = (R * dlat) ** 2
    dy2 = dy2.where(dy2 >= min_dy2)

    dFdy_raw = calc_gradient(data_input, dim=lat_dim, edge_order=edge_order)
    d2Fdy2_raw = calc_gradient(dFdy_raw, dim=lat_dim, edge_order=edge_order)
    d2Fdy2 = d2Fdy2_raw / dy2
    return d2Fdy2


def calc_lon_lat_mixed_derivatives(
    data_input: xr.DataArray | xr.Dataset,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    min_dxdy: float = 1e10,
    edge_order: int = 2,
    R: float = 6370000,
) -> xr.DataArray | xr.Dataset:
    """
    Calculation of second-order mixed partial derivative terms along longitude and latitude.

    .. math::
        \\frac{\\partial^2 F}{\\partial x \\partial y} = \\frac{1}{R^2 \\cos\\varphi} \\cdot \\frac{\\partial^2 F}{\\partial \\lambda \\partial \\varphi}

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    min_dxdy: :py:class:`float <float>`, default: `1e10`.
        The minimum acceptable value of :math:`\\mathrm{d}x\\mathrm{d}y`, below which parts will set `nan` to avoid large computational errors.
        The unit is m. You can set it to a negative value in order to remove this benefit.
    edge_order: {1, 2}, optional
        Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
    R: :py:class:`float <float>`, default: `6370000`.
        Radius of the Earth.

    Returns
    -------
    The second-order mixed partial derivative terms along longitude and latitude (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. seealso::
        :py:func:`calc_gradient <calc_gradient>`
    """
    lon_array = data_input[lon_dim].astype("float64")
    lat_array = data_input[lat_dim].astype("float64")

    dlon = transfer_deg2rad(calc_gradient(lon_array, dim=lon_dim))
    dlat = transfer_deg2rad(calc_gradient(lat_array, dim=lat_dim))
    coslat = np.cos(transfer_deg2rad(lat_array))
    dxdy = R**2 * coslat * dlon * dlat
    dxdy = dxdy.where(np.abs(dxdy) >= min_dxdy)

    dFdy_raw = calc_gradient(data_input, dim=lat_dim, edge_order=edge_order)
    d2Fdxy_raw = calc_gradient(dFdy_raw, dim=lon_dim, edge_order=edge_order)
    d2Fdxy = d2Fdxy_raw / dxdy
    return d2Fdxy


def calc_p_gradient(
    data_input: xr.DataArray, vertical_dim: str, vertical_dim_units: str
) -> xr.DataArray:
    """
    Calculate the gradient along the barometric pressure direction in the p-coordinate system.

    .. math::
        \\frac{\\partial F}{\\partial p}

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.

    Returns
    -------
    The gradient along the barometric pressure direction in the p-coordinate system (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. seealso::
        :py:func:`calc_gradient <calc_gradient>`
    """
    if isinstance(data_input.data, dask.array.core.Array):
        # The vertical coordinate dimension is set to no chunking
        data_input = data_input.chunk({vertical_dim: -1})

    # Convert the pressure unit to Pascal
    dp_base = transfer_units_coeff(vertical_dim_units, "Pa")

    dp = calc_gradient(data_input[vertical_dim], dim=vertical_dim) * dp_base
    dF_dp = calc_gradient(data_input, dim=vertical_dim) / dp
    return dF_dp


def calc_time_gradient(
    data_input: xr.DataArray,
    time_units: str,
    time_dim: str = "time",
) -> xr.DataArray:
    """
    Calculate the gradient along the time direction.

    .. math::
        \\frac{\\partial F}{\\partial t}

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
    time_units: :py:class:`str <str>`.
        The unit corresponding to the time dimension value. Optional values are `seconds`, `months`, `years` and so on.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.

    Returns
    -------
    The gradient along the time direction (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. caution:: The units for partial derivative of `time` are :math:`\mathrm{s^{-1}}`.

    .. seealso::
        :py:func:`calc_gradient <calc_gradient>`
    """
    if isinstance(data_input.data, dask.array.core.Array):
        # The vertical coordinate dimension is set to no chunking
        data_input = data_input.chunk({time_dim: -1})

    # Convert time units to seconds
    dt = transfer_units_coeff(time_units, "seconds")
    dFdt = calc_gradient(data_input, dim=time_dim) / dt
    dFdt.attrs["units"] = "s^-1"
    return dFdt


def calc_delta_pressure(
    data_input: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: str,
    surface_pressure_data_units: str,
) -> xr.DataArray:
    """
    Calculates the pressure layer thickness (delta pressure) of a constant
    pressure level coordinate system.

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The spatio-temporal data to be calculated.
    surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Mean surface sea level pressure.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    surface_pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.

    Returns
    -------
    The pressure layer thickness (delta pressure) of a constant pressure level coordinate system (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. seealso::
        - :py:func:`geocat.comp.meteorology.delta_pressure <geocat-comp:geocat.comp.meteorology.delta_pressure>`
        - `dpres_plevel - NCL <https://www.ncl.ucar.edu/Document/Functions/Built-in/dpres_plevel.shtml>`__
    """
    # `vertical_dim` data is forced into descending order
    data_input = data_input.sortby(vertical_dim, ascending=False)

    # Change the vertical coordinate unit to `Pa`
    vertical_dim_base = transfer_units_coeff(vertical_dim_units, "Pa")
    data_input = data_input.assign_coords(
        {vertical_dim: data_input[vertical_dim] * vertical_dim_base}
    )

    # Change the surface pressure data unit to `Pa`
    surface_pressure_data = transfer_data_units(
        surface_pressure_data, surface_pressure_data_units, "Pa"
    )

    pressure_lev = xr.broadcast(data_input[vertical_dim], data_input)[0]
    pressure_top = pressure_lev.min(dim=vertical_dim)
    surface_pressure = surface_pressure_data

    delta_pressure = np.abs(calc_gradient(pressure_lev, dim=vertical_dim))
    # top level
    delta_pressure[{vertical_dim: -1}] = (
        pressure_lev[{vertical_dim: -1}] + pressure_lev[{vertical_dim: -2}]
    ) / 2 - pressure_top
    # bottom level
    delta_pressure[{vertical_dim: 0}] = (
        surface_pressure
        - (pressure_lev[{vertical_dim: 0}] + pressure_lev[{vertical_dim: 1}]) / 2
    )
    return delta_pressure


def calc_p_integral(
    data_input: xr.DataArray, vertical_dim: str, normalize: bool = True
) -> xr.DataArray:
    """
    Calculate the vertical integral along the barometric pressure direction in the p-coordinate system.

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The spatio-temporal data to be calculated.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    normalize: :py:class:`bool<bool>`, default: `True`.
        Whether or not the integral results are averaged over the entire layer.

    Returns
    -------
    The vertical integral along the barometric pressure direction in the p-coordinate system (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. attention::
        This method ignores the effect of topography, so it applies to altitudes **above 900hPa** and is **NOT applicable to the Tibetan Plateau region**.
        For a fully accurate vertical integration, please use the :py:func:`calc_top2surface_integral <calc_top2surface_integral>` function to calculate,
        but the speed of the calculation is slightly slowed down.
    """
    dim_index = find_dims_axis(data_input, dim=vertical_dim)
    # Get level data
    x = data_input[vertical_dim].data
    # Integrate along the given axis using the composite trapezoidal rule
    data_trapz = np.trapz(data_input.data, x, axis=dim_index)
    data_trapz = np.abs(data_trapz)

    if normalize == True:
        # Normalize
        data_trapz_normalized = data_trapz / (np.max(x) - np.min(x))
    elif normalize == False:
        data_trapz_normalized = data_trapz
    else:
        raise ValueError("The parameter `normalize` should be `True` or `False`.")

    return (
        data_input.isel({vertical_dim: 0})
        .drop_vars(vertical_dim)
        .copy(data=data_trapz_normalized, deep=True)
    )


def calc_top2surface_integral(
    data_input: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    vertical_dim: str,
    surface_pressure_data_units: str,
    vertical_dim_units: str,
    method: str = "Trenberth-vibeta",
    normalize: bool = True,
) -> xr.DataArray:
    """
    Calculate the vertical integral in the p-coordinate system from the ground to the zenith along the barometric pressure direction.

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The spatio-temporal data to be calculated.
    surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Mean surface sea level pressure.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    surface_pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    method: :py:class:`str <str>`, default: `'Trenberth-vibeta'`.
        vertical integration method. Optional values are `Boer-vibeta`, `'Trenberth-vibeta'`.

        .. note::
            The trapezoidal rule of integration is exactly equivalent to

            .. math::
                I = \\sum_{j=1,2J-1,2} (\\beta M)_j \\Delta p_j,

            where Kevin E. Trenberth (1991) define

            .. math::
                \\beta_j = \\left\\lbrace
                \\begin{array}{ll}
                1, & \\mathrm{if} \\ p_{j-1} < p_s,\\\\
                0, & \\mathrm{if} \\ p_{j+1} > p_s ,\\\\
                \\frac{p_s - p_{j+1}}{p_{j-1} - p_{j+1}}, & \\mathrm{if}  \\ p_{j-1} > p_s > p_{j+1}.
                \\end{array}
                \\right.

            While G. J. Boer (1982) define :math:`\\beta = 0, 1` only.

    normalize: :py:class:`bool<bool>`, default: `True`.
        Whether or not the integral results are averaged over the entire layer.

    Returns
    -------
    The vertical integral in the p-coordinate system from the ground to the zenith along the barometric pressure direction. (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - `Boer, G. J., 1982: Diagnostic Equations in Isobaric Coordinates. Mon. Wea. Rev., 110, 1801–1820, <https://doi.org/10.1175/1520-0493(1982)110%3C1801:DEIIC%3E2.0.CO;2>`__
    - `Trenberth, K. E., 1991: Climate Diagnostics from Global Analyses: Conservation of Mass in ECMWF Analyses. J. Climate, 4, 707–722, <https://doi.org/10.1175/1520-0442(1991)004%3C0707:CDFGAC%3E2.0.CO;2>`__

    .. seealso::
        - `vibeta - NCL <https://www.ncl.ucar.edu/Document/Functions/Built-in/vibeta.shtml>`__
        - `dpres_plevel - NCL <https://www.ncl.ucar.edu/Document/Functions/Built-in/dpres_plevel.shtml>`__
    """
    # `vertical_dim` data is forced into descending order
    data_input = data_input.sortby(vertical_dim, ascending=False)

    if method == "Boer-vibeta":
        # https://doi.org/10.1175/1520-0493(1982)110%3C1801:DEIIC%3E2.0.CO;2
        # Change the vertical coordinate unit to `Pa`
        vertical_dim_array = transfer_data_units(
            data_input[vertical_dim], vertical_dim_units, "Pa"
        )

        # Change the surface pressure data unit to `Pa`
        surface_pressure_data = transfer_data_units(
            surface_pressure_data, surface_pressure_data_units, "Pa"
        )

        # Filtering the part of the air pressure level greater than the surface pressure
        level_expanded = xr.broadcast(vertical_dim_array, data_input)[0]
        pressure2surface_data = data_input.where(
            level_expanded <= surface_pressure_data
        )
        pressure2surface_data = pressure2surface_data.fillna(0)
        fint = calc_p_integral(
            pressure2surface_data, vertical_dim=vertical_dim, normalize=normalize
        )

        # Clear unnecessary attributes
        fint.attrs = dict()
        # Add `Vertical integral method` property
        fint.attrs["Vertical integral method"] = "Boer-vibeta"
        return fint

    elif method == "Trenberth-vibeta":
        # https://doi.org/10.1175/1520-0442(1991)004%3C0707:CDFGAC%3E2.0.CO;2
        # https://www.ncl.ucar.edu/Document/Functions/Built-in/vibeta.shtml
        dp = calc_delta_pressure(
            data_input,
            surface_pressure_data,
            vertical_dim,
            surface_pressure_data_units,
            vertical_dim_units,
        )

        # Change the vertical coordinate unit to `Pa`
        vertical_dim_base = transfer_units_coeff(vertical_dim_units, "Pa")
        data_input = data_input.assign_coords(
            {vertical_dim: data_input[vertical_dim] * vertical_dim_base}
        )

        # weighted sum/sum_of_layer_thickness
        part1 = (data_input * dp).sum(dim=vertical_dim)

        if normalize == True:
            part2 = dp.sum(dim=vertical_dim)
            fint = part1 / part2
        elif normalize == False:
            fint = part1
        else:
            raise ValueError("The parameter `normalize` should be `True` or `False`.")

        # Add `Vertical integral method` property
        fint.attrs["Vertical integral method"] = "Trenberth-vibeta"
        return fint

    else:
        raise ValueError(
            "The parameter `method` should be `'Boer-vibeta'` or `'Trenberth-vibeta'`."
        )


def calc_laplacian(
    data_input: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    R: float = 6370000,
    spherical_coord: bool = True,
) -> xr.DataArray:
    """
    Calculate the horizontal Laplace term.

    rectangular coordinates

    .. math::
        \\nabla^2 F = \\frac{\\partial^2 F}{\\partial x^2} + \\frac{\\partial^2 F}{\\partial y^2}

    Spherical coordinates

    .. math::
        \\nabla^2 F = \\frac{\\partial^2 F}{\\partial x^2} + \\frac{\\partial^2 F}{\\partial y^2} - \\frac{1}{R} \\frac{\\partial F}{\\partial y} \\tan \\varphi

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The spatio-temporal data to be calculated.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    R: :py:class:`float <float>`, default: `6370000`.
        Radius of the Earth.
    spherical_coord: :py:class:`bool <bool>`, default: `True`.
        Whether or not to compute the horizontal Laplace term in spherical coordinates.

    Returns
    -------
    The horizontal Laplace term. (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    dlon = transfer_deg2rad(calc_gradient(data_input[lon_dim], dim=lon_dim))
    dlat = transfer_deg2rad(calc_gradient(data_input[lat_dim], dim=lat_dim))
    coslat = np.cos(transfer_deg2rad(data_input[lat_dim]))
    dx = R * coslat * dlon
    dy = R * dlat

    d2Hdx2 = data_input.diff(dim=lon_dim, n=2) / (dx**2)
    d2Hdy2 = data_input.diff(dim=lat_dim, n=2) / (dy**2)

    if spherical_coord == True:
        term3 = (
            1
            / R
            * calc_gradient(data_input, dim=lat_dim)
            / dy
            * np.tan(transfer_deg2rad(data_input[lat_dim]))
        )
        laplacian = d2Hdx2 + d2Hdy2 - term3
    elif spherical_coord == False:
        laplacian = d2Hdx2 + d2Hdy2
    return laplacian


def calc_divergence(
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    R: float = 6370000,
    spherical_coord=True,
) -> xr.DataArray:
    """
    Calculate the horizontal divergence term.

    rectangular coordinates

    .. math::
        \\mathrm{D} = \\frac{\\partial u}{\\partial x} + \\frac{\\partial v}{\\partial y}

    Spherical coordinates

    .. math::
        \\mathrm{D} = \\frac{\\partial u}{\\partial x} + \\frac{\\partial v}{\\partial y} - \\frac{v}{R} \\tan \\varphi

    Parameters
    ----------
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    R: :py:class:`float <float>`, default: `6370000`.
        Radius of the Earth.
    spherical_coord: :py:class:`bool<bool>`, default: `True`.
        Whether or not to compute the horizontal Laplace term in spherical coordinates.

    Returns
    -------
    The horizontal divergence term. (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    if u_data.shape != v_data.shape:
        raise ValueError(
            "`u`, `v` shape must be same! However, the u shape is ",
            u_data.shape,
            ", and v shape is ",
            v_data.shape,
            ".",
        )

    metadata = u_data
    dlon = transfer_deg2rad(calc_gradient(metadata[lon_dim], dim=lon_dim))
    dlat = transfer_deg2rad(calc_gradient(metadata[lat_dim], dim=lat_dim))
    coslat = np.cos(transfer_deg2rad(metadata[lat_dim]))
    dx = R * coslat * dlon
    dy = R * dlat

    dudx = calc_gradient(u_data, dim=lon_dim) / dx
    dvdy = calc_gradient(v_data, dim=lat_dim) / dy

    if spherical_coord == True:
        term3 = v_data / R * np.tan(transfer_deg2rad(metadata[lat_dim]))
        div = dudx + dvdy - term3
    elif spherical_coord == False:
        div = dudx + dvdy
    return div


def calc_vorticity(
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    R: float = 6370000,
    spherical_coord: bool = True,
) -> xr.DataArray:
    """
    Calculate the horizontal relative vorticity term.

    rectangular coordinates

    .. math::
        \\zeta = \\frac{\\partial v}{\\partial x} - \\frac{\\partial u}{\\partial y}

    Spherical coordinates

    .. math::
        \\zeta = \\frac{\\partial v}{\\partial x} - \\frac{\\partial u}{\\partial y} + \\frac{u}{R} \\tan \\varphi

    Parameters
    ----------
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    R: :py:class:`float <float>`, default: `6370000`.
        Radius of the Earth.
    spherical_coord: :py:class:`bool<bool>`, default: `True`.
        Whether or not to compute the horizontal Laplace term in spherical coordinates.

    Returns
    -------
    The horizontal relative vorticity term. (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    if u_data.shape != v_data.shape:
        raise ValueError(
            "`u`, `v` shape must be same! However, the u shape is ",
            u_data.shape,
            ", and v shape is ",
            v_data.shape,
            ".",
        )

    metadata = u_data
    dlon = transfer_deg2rad(calc_gradient(metadata[lon_dim], dim=lon_dim))
    dlat = transfer_deg2rad(calc_gradient(metadata[lat_dim], dim=lat_dim))
    coslat = np.cos(transfer_deg2rad(metadata[lat_dim]))
    dx = R * coslat * dlon
    dy = R * dlat

    dvdx = calc_gradient(v_data, dim=lon_dim) / dx
    dudy = calc_gradient(u_data, dim=lat_dim) / dy

    if spherical_coord == True:
        term3 = u_data / R * np.tan(transfer_deg2rad(metadata[lat_dim]))
        vor = dvdx - dudy + term3
    elif spherical_coord == False:
        vor = dvdx - dudy
    return vor


def calc_geostrophic_wind(
    z_data: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    omega: float = 7.292e-5,
    g: float = 9.8,
    R: float = 6370000,
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
    """
    dlon = transfer_deg2rad(calc_gradient(z_data[lon_dim], dim=lon_dim))
    dlat = transfer_deg2rad(calc_gradient(z_data[lat_dim], dim=lat_dim))
    coslat = np.cos(transfer_deg2rad(z_data[lat_dim]))
    dx = R * coslat * dlon
    dy = R * dlat
    f = 2 * omega * np.sin(transfer_deg2rad(z_data[lat_dim]))

    dHdy = calc_gradient(z_data, dim=lat_dim) / dy
    dHdx = calc_gradient(z_data, dim=lon_dim) / dx

    dataset = xr.Dataset()
    ug = -(g / f) * dHdy
    dataset["ug"] = transfer_inf2nan(ug)
    vg = (g / f) * dHdx
    dataset["vg"] = transfer_inf2nan(vg)

    return dataset


def calc_geostrophic_wind_vorticity(
    z_data: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    spherical_coord: bool = True,
    omega: float = 7.292e-5,
    g: float = 9.8,
    R: float = 6370000,
) -> xr.DataArray:
    """
    Calculate the geostrophic vorticity.

    rectangular coordinates

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

    Returns
    -------
    The geostrophic vorticity term. (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    geostrophic_wind = calc_geostrophic_wind(
        z_data, lon_dim=lon_dim, lat_dim=lat_dim, omega=omega, g=g, R=R
    )
    ug, vg = geostrophic_wind["ug"], geostrophic_wind["vg"]
    vor_g = calc_vorticity(
        ug, vg, spherical_coord=spherical_coord, lon_dim=lon_dim, lat_dim=lat_dim, R=R
    )
    return vor_g


def calc_horizontal_water_flux(
    specific_humidity_data: xr.DataArray,
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    g: float = 9.8,
) -> xr.Dataset:
    """
    Calculate horizontal water vapor flux at each vertical level.

    .. math::
        \\frac{1}{g} q \\mathbf{V} = \\frac{1}{g} (u q\\ \\mathbf{i} + vq\\ \\mathbf{j})

    Parameters
    ----------
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    g: :py:class:`float <float>`, default: `9.8`.
        The acceleration of gravity.

    Returns
    -------
    The water vapor flux. (:py:class:`xarray.Dataset<xarray.Dataset>`).

    - :math:`qu`: zonal water vapor flux.
    - :math:`qv`: meridional water vapor flux.
    """
    water_flux = xr.Dataset(
        data_vars={
            "qu": specific_humidity_data * u_data / g,
            "qv": specific_humidity_data * v_data / g,
        }
    )
    return water_flux


def calc_vertical_water_flux(
    specific_humidity_data: xr.DataArray, omega_data: xr.DataArray, g: float = 9.8
) -> xr.DataArray:
    """
    Calculate vertical water vapor flux.

    .. math::
        -\\omega \\frac{q}{g}

    Parameters
    ----------
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    omega_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The vertical velocity data (:math:`\\frac{\\mathrm{d} p}{\\mathrm{d} t}`).
    g: :py:class:`float <float>`, default: `9.8`.
        The acceleration of gravity.

    Returns
    -------
    The vertical water flux. (:py:class:`xarray.DataArray <xarray.DataArray>`).
    """
    water_flux = -omega_data * specific_humidity_data / g
    return water_flux


def calc_water_flux_top2surface_integral(
    specific_humidity_data: xr.DataArray,
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    surface_pressure_data_units: str,
    vertical_dim: str,
    vertical_dim_units: str,
    method: str = "Trenberth-vibeta",
    g: float = 9.8,
) -> xr.DataArray:
    """
    Calculate the water vapor flux across the vertical level.

    Parameters
    ----------
    specific_humidity: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Mean surface sea level pressure.
    surface_pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    method: :py:class:`str <str>`, default: `'Trenberth-vibeta'`.
        Vertical integration method. Optional values are `Boer-vibeta`, `'Trenberth-vibeta'`.
    g: :py:class:`float <float>`, default: `9.8`.
        The acceleration of gravity.

    Returns
    -------
    The water vapor flux. (:py:class:`xarray.Dataset<xarray.Dataset>`).

    - :math:`qu`: zonal water vapor flux.
    - :math:`qv`: meridional water vapor flux.

    .. seealso::
        :py:func:`calc_top2surface_integral <calc_top2surface_integral>`
    """
    # Calculate the single-layer water flux
    water_flux_single_layer = calc_horizontal_water_flux(
        specific_humidity_data, u_data, v_data, g=g
    )
    qu = water_flux_single_layer["qu"]
    qv = water_flux_single_layer["qv"]

    # Calculate the integral of `qu` over the whole atmosphere
    qu_top2surface_integral = calc_top2surface_integral(
        data_input=qu,
        vertical_dim=vertical_dim,
        vertical_dim_units=vertical_dim_units,
        surface_pressure_data=surface_pressure_data,
        surface_pressure_data_units=surface_pressure_data_units,
        method=method,
        normalize=False,
    )

    # Calculate the integral of `qv` over the whole atmosphere
    qv_top2surface_integral = calc_top2surface_integral(
        data_input=qv,
        vertical_dim=vertical_dim,
        vertical_dim_units=vertical_dim_units,
        surface_pressure_data=surface_pressure_data,
        surface_pressure_data_units=surface_pressure_data_units,
        method=method,
        normalize=False,
    )
    quv_top2surface_integral = xr.Dataset(
        data_vars={"qu": qu_top2surface_integral, "qv": qv_top2surface_integral}
    )

    return quv_top2surface_integral


def calc_divergence_watervaporflux(
    specific_humidity_data: xr.DataArray,
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    specific_humidity_units: str,
    spherical_coord: bool = True,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    g: float = 9.8,
    R: float = 6370000,
) -> xr.DataArray:
    """
    Calculate water vapor flux divergence at each vertical level.

    .. math::
        \\nabla \\left( \\frac{1}{g} q \\mathbf{V} \\right) = \\frac{1}{g} \\nabla \\cdot \\left( q \\mathbf{V} \\right)


    Parameters
    ----------
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    specific_humidity_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
    spherical_coord: :py:class:`bool<bool>`, default: `True`.
        Whether or not to compute the horizontal Laplace term in spherical coordinates.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    g: :py:class:`float <float>`, default: `9.8`.
        The acceleration of gravity.
    R: :py:class:`float <float>`, default: `6370000`.
        Radius of the Earth.

    Returns
    -------
    The water vapor flux divergence. (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    specific_humidity = transfer_data_units(
        specific_humidity_data, specific_humidity_units, "kg/kg"
    )
    divergence_watervaporflux = (1 / g) * calc_divergence(
        u_data=specific_humidity * u_data,
        v_data=specific_humidity * v_data,
        spherical_coord=spherical_coord,
        lon_dim=lon_dim,
        lat_dim=lat_dim,
        R=R,
    )

    return divergence_watervaporflux


def calc_divergence_watervaporflux_top2surface_integral(
    specific_humidity_data: xr.DataArray,
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    vertical_dim: str,
    specific_humidity_units: str,
    surface_pressure_data_units: str,
    vertical_dim_units: str,
    spherical_coord: bool = True,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    method: str = "Trenberth-vibeta",
    g: float = 9.8,
    R: float = 6370000,
) -> xr.DataArray:
    """
    Calculate water vapor flux divergence across the vertical level.

    Parameters
    ----------
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Mean surface sea level pressure.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    specific_humidity_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
    surface_pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    spherical_coord: :py:class:`bool<bool>`, default: `True`.
        Whether or not to compute the horizontal Laplace term in spherical coordinates.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    g: :py:class:`float <float>`, default: `9.8`.
        The acceleration of gravity.
    R: :py:class:`float <float>`, default: `6370000`.
        Radius of the Earth.

    Returns
    -------
    The water vapor flux divergence. (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    specific_humidity = transfer_data_units(
        specific_humidity_data, specific_humidity_units, "kg/kg"
    )

    # Calculation of single-layer water vapor flux divergence
    divergence_watervaporflux_single_layer = calc_divergence_watervaporflux(
        specific_humidity,
        u_data,
        v_data,
        spherical_coord=spherical_coord,
        specific_humidity_units=specific_humidity_units,
        lon_dim=lon_dim,
        lat_dim=lat_dim,
        g=g,
        R=R,
    )

    # Calculate the whole layer water vapor flux divergence
    divergence_watervaporflux_single_layer_top2surface_integral = (
        calc_top2surface_integral(
            data_input=divergence_watervaporflux_single_layer,
            vertical_dim=vertical_dim,
            vertical_dim_units=vertical_dim_units,
            surface_pressure_data=surface_pressure_data,
            surface_pressure_data_units=surface_pressure_data_units,
            method=method,
            normalize=False,
        )
    )

    return divergence_watervaporflux_single_layer_top2surface_integral


def calc_u_advection(
    u_data: xr.DataArray,
    temper_data: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    min_dx: float = 1.0,
    edge_order: int = 2,
    R: float = 6370000,
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
    """
    return (
        (-1)
        * u_data
        * calc_lon_gradient(
            temper_data,
            lon_dim=lon_dim,
            lat_dim=lat_dim,
            min_dx=min_dx,
            edge_order=edge_order,
            R=R,
        )
    )


def calc_v_advection(
    v_data: xr.DataArray,
    temper_data: xr.DataArray,
    lat_dim: str = "lat",
    min_dy: float = 1.0,
    edge_order: int = 2,
    R: float = 6370000,
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
    """
    return (
        (-1)
        * v_data
        * calc_lat_gradient(
            temper_data, lat_dim=lat_dim, min_dy=min_dy, edge_order=edge_order, R=R
        )
    )


def calc_p_advection(
    omega_data: xr.DataArray,
    temper_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: str,
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
    return -omega_data * calc_p_gradient(temper_data, vertical_dim, vertical_dim_units)
