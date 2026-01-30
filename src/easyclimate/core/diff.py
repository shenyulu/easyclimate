"""
The calculation of geographic finite difference
"""

from __future__ import annotations

import numpy as np
from .utility import (
    find_dims_axis,
    transfer_deg2rad,
    transfer_inf2nan,
    generate_dataset_dispatcher,
    compare_multi_dataarray_coordinate,
)
from .units import (
    transfer_units_coeff,
    transfer_data_multiple_units,
)
from ..backend import dvibeta, dvrfidf, ddvfidf
import xarray as xr
import dask
from typing import Literal

__all__ = [
    "calc_gradient",
    "calc_dx_gradient",
    "calc_dlon_radian_gradient",
    "calc_dlon_degree_gradient",
    "calc_dy_gradient",
    "calc_dlat_radian_gradient",
    "calc_dlat_degree_gradient",
    "calc_dx_laplacian",
    "calc_dy_laplacian",
    "calc_dxdy_mixed_derivatives",
    "calc_p_gradient",
    "calc_time_gradient",
    "calc_delta_pressure",
    "calc_p_integral",
    "calc_top2surface_integral",
    "calc_dxdy_laplacian",
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
    "calc_shear_stretch_deform",
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
        Gradient is calculated using N-th order accurate differences at the boundaries. Default: 2.

    Returns
    -------
    The gradient along the coordinate `dim` direction (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. seealso::
        :py:func:`numpy.gradient <numpy:numpy.gradient>`

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """

    def _calc_gradient(data_input, dim, varargs, edge_order) -> xr.DataArray:
        dim_index = find_dims_axis(data_input, dim=dim)
        # `data_input.data` make dask to process Dask.array
        data_gradient = np.gradient(
            data_input.data, varargs, axis=dim_index, edge_order=edge_order
        )
        return data_input.copy(data=data_gradient, deep=True)

    result = _calc_gradient(data_input, dim, varargs, edge_order)

    # clean attrs
    result.attrs = dict()
    return result


def calc_dx_gradient(
    data_input: xr.DataArray | xr.Dataset,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    min_dx: float = 1.0,
    edge_order: int = 2,
    R: float = 6371200.0,
) -> xr.DataArray | xr.Dataset:
    """
    Calculate the zonal gradient of the input data in physical units (meters).

    This function computes the partial derivative :math:`\\partial F / \\partial x`,
    where :math:`x` is the eastward distance along a parallel (in meters). It is
    the full physical zonal gradient on a sphere, given by:

    .. math::

        \\frac{\\partial F}{\\partial x} = \\frac{1}{R \\cos \\varphi} \\cdot \\frac{\\partial F}{\\partial \\lambda}

    where :math:`R` is the Earth's radius, :math:`\\varphi` is latitude, and
    :math:`\\lambda` is longitude in radians. This is essential for dynamical
    calculations like advection or wave propagation in atmospheric/oceanic models.

    The computation uses finite differences along the longitude dimension:
    :math:`\\partial F / \\partial x = (\\partial F / \\partial i) / (\\partial x / \\partial i)`,
    where :math:`i` is the grid index. Longitude is assumed in degrees and
    converted to radians; latitude is broadcasted for the cosine factor.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated. If a Dataset, the gradient is applied to all
        data variables.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default, the gradient is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. Used to compute the cosine factor; broadcasted if necessary.
    min_dx: :py:class:`float <float>`, default: `1.0` (:math:`\\mathrm{m}`).
        The minimum acceptable value of `dx` (zonal spacing in meters), below which the output
        is set to NaN to avoid numerical instabilities from very small grid spacings.
        Set to a negative value to disable this check.
    edge_order: {1, 2}, optional
        Order of the finite difference used at the boundaries. 1 uses first-order accurate
        one-sided differences; 2 uses second-order accurate one-sided differences. Default: 2.
    R: :py:class:`float <float>`, default: `6370000` (:math:`\\mathrm{m}`).
        Radius of the Earth in meters (approximate mean radius). Can be adjusted for specific models.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The zonal gradient :math:`\\partial F / \\partial x`, with the same shape and coordinates
        as the input. Units are those of the input data divided by meters (e.g., if F is in K,
        output is K/m). Invalid regions (dx < min_dx) are NaN.

    .. seealso::
        - :py:func:`calc_gradient <calc_gradient>`
        - :py:func:`calc_dlon_radian_gradient <calc_dlon_radian_gradient>`
        - :py:func:`calc_dy_gradient <calc_dy_gradient>`

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
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
    dFdx = dFdx.astype("float32")

    # Add attributes for clarity
    if isinstance(dFdx, xr.DataArray):
        input_units = data_input.attrs.get("units", "dimensionless")
        dFdx.attrs.update(
            {
                "long_name": f'Zonal gradient of {data_input.name if data_input.name else "data"}',
                "units": f"{input_units} per meter",
                "standard_name": "partial_derivative_of_data_with_respect_to_zonal_distance",
            }
        )
        dFdx.name = f'd{dFdx.name or "data"}_dx' if dFdx.name else "dF_dx"
    elif isinstance(dFdx, xr.Dataset):
        for var_name, var in dFdx.data_vars.items():
            input_units = (
                data_input[var_name].attrs.get("units", "dimensionless")
                if var_name in data_input
                else "dimensionless"
            )
            var.attrs.update(
                {
                    "long_name": f"Zonal gradient of {var_name}",
                    "units": f"{input_units} per meter",
                    "standard_name": f"partial_derivative_of_{var_name}_with_respect_to_zonal_distance",
                }
            )
            dFdx[var_name].name = f"d{var_name}_dx"

    return dFdx


def calc_dlon_radian_gradient(
    data_input: xr.DataArray | xr.Dataset,
    lon_dim: str = "lon",
    edge_order: int = 2,
) -> xr.DataArray | xr.Dataset:
    """
    Calculate the gradient of the input data with respect to longitude in radians.

    This function computes the partial derivative :math:`\\partial F / \\partial \\lambda`,
    where :math:`\\lambda` is the longitude in radians. It is useful for spherical coordinate
    calculations, such as in wave activity flux (WAF) formulations, where angular gradients
    must be in radians for consistency with trigonometric functions and the Earth's radius.

    The computation uses finite differences:

    .. math::

        \\frac{\\partial F}{\\partial \\lambda} = \\frac{\\partial F}{\\partial i } / \\frac{\\partial \\lambda}{\\partial i},

    where :math:`i` is the grid index along the longitude dimension. Longitude values are
    assumed to be in degrees initially and converted to radians for the denominator.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated. If a Dataset, the gradient is applied to all data variables.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default, the gradient is applied over the `lon` dimension.
    edge_order: {1, 2}, optional
        Order of the finite difference used at the boundaries. 1 uses first-order accurate
        one-sided differences; 2 uses second-order accurate one-sided differences. Default: 2.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The gradient :math:`\\partial F / \\partial \\lambda` (in radians), with the same shape
        and coordinates as the input. Units are inherited from the input data divided by radians
        (e.g., if F is in K, output is K/rad).
    """
    # Set to `float64` for more accurate results in trigonometric calculations.
    lon_array = data_input[lon_dim].astype("float64")
    dlon = transfer_deg2rad(calc_gradient(lon_array, dim=lon_dim))

    dFdx_raw = calc_gradient(data_input, dim=lon_dim, edge_order=edge_order)
    dFdx = dFdx_raw / dlon
    dFdx = dFdx.astype("float32")

    # Add attributes for clarity
    if isinstance(dFdx, xr.DataArray):
        dFdx.attrs.update(
            {
                "long_name": f'Gradient of {data_input.name if data_input.name else "data"} with respect to longitude (radians)',
                "units": f'{data_input.attrs.get("units", "dimensionless")} per radian',
                "standard_name": "partial_derivative_of_data_with_respect_to_longitude_radians",
            }
        )
        dFdx.name = (
            f'd{dFdx.name or "data"}_dlambda_rad' if dFdx.name else "dF_dlambda_rad"
        )
    elif isinstance(dFdx, xr.Dataset):
        for var_name, var in dFdx.data_vars.items():
            var.attrs.update(
                {
                    "long_name": f"Gradient of {var_name} with respect to longitude (radians)",
                    "units": f'{data_input[var_name].attrs.get("units", "dimensionless")} per radian',
                    "standard_name": f"partial_derivative_of_{var_name}_with_respect_to_longitude_radians",
                }
            )
            dFdx[var_name].name = f"d{var_name}_dlambda_rad"

    return dFdx


def calc_dlon_degree_gradient(
    data_input: xr.DataArray | xr.Dataset,
    lon_dim: str = "lon",
    edge_order: int = 2,
) -> xr.DataArray | xr.Dataset:
    """
    Calculate the gradient of the input data with respect to longitude in degrees.

    This function computes the partial derivative :math:`\\partial F / \\partial \\lambda`,
    where :math:`\\lambda` is the longitude in degrees. It is suitable for general-purpose
    gradient calculations where degree units are preferred for interpretability.

    The computation uses finite differences:

    .. math::

        \\frac{\\partial F}{\\partial \\lambda} = \\frac{\\partial F}{\\partial i} / \\frac{\\partial \\lambda}{\\partial i},

    where :math:`i` is the grid index along the longitude dimension. Longitude values remain
    in degrees for the denominator.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated. If a Dataset, the gradient is applied to all
        data variables.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default, the gradient is applied over the `lon` dimension.
    edge_order: {1, 2}, optional
        Order of the finite difference used at the boundaries. 1 uses first-order accurate
        one-sided differences; 2 uses second-order accurate one-sided differences. Default: 2.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The gradient :math:`\\partial F / \\partial \\lambda` (in degrees), with the same shape
        and coordinates as the input. Units are inherited from the input data divided by degrees
        (e.g., if F is in K, output is K/deg).

    .. seealso::
        - :py:func:`calc_gradient <calc_gradient>`
        - :py:func:`calc_dlon_radian_gradient <calc_dlon_radian_gradient>`
        - :py:func:`calc_dlat_degree_gradient <calc_dlat_degree_gradient>`

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    # Set to `float64` for more accurate results in trigonometric calculations.
    lon_array = data_input[lon_dim].astype("float64")
    dlon = calc_gradient(lon_array, dim=lon_dim)

    dFdx_raw = calc_gradient(data_input, dim=lon_dim, edge_order=edge_order)
    dFdx = dFdx_raw / dlon
    dFdx = dFdx.astype("float32")

    # Add attributes for clarity
    if isinstance(dFdx, xr.DataArray):
        dFdx.attrs.update(
            {
                "long_name": f'Gradient of {data_input.name if data_input.name else "data"} with respect to longitude (degrees)',
                "units": f'{data_input.attrs.get("units", "dimensionless")} per degree',
                "standard_name": "partial_derivative_of_data_with_respect_to_longitude_degrees",
            }
        )
        dFdx.name = (
            f'd{dFdx.name or "data"}_dlambda_deg' if dFdx.name else "dF_dlambda_deg"
        )
    elif isinstance(dFdx, xr.Dataset):
        for var_name, var in dFdx.data_vars.items():
            var.attrs.update(
                {
                    "long_name": f"Gradient of {var_name} with respect to longitude (degrees)",
                    "units": f'{data_input[var_name].attrs.get("units", "dimensionless")} per degree',
                    "standard_name": f"partial_derivative_of_{var_name}_with_respect_to_longitude_degrees",
                }
            )
            dFdx[var_name].name = f"d{var_name}_dlambda_deg"

    return dFdx


def calc_dy_gradient(
    data_input: xr.DataArray | xr.Dataset,
    lat_dim: str = "lat",
    min_dy: float = 1.0,
    edge_order: int = 2,
    R: float = 6371200.0,
) -> xr.DataArray | xr.Dataset:
    """
    Calculate the meridional gradient of the input data in physical units (meters).

    This function computes the partial derivative :math:`\\partial F / \\partial y`,
    where :math:`y` is the northward distance along a meridian (in meters). It is
    the full physical meridional gradient on a sphere, given by:

    .. math::

        \\frac{\\partial F}{\\partial y} = \\frac{1}{R} \\cdot \\frac{\\partial F}{\\partial \\varphi}

    where :math:`R` is the Earth's radius and :math:`\\varphi` is latitude in radians.
    This is essential for dynamical calculations like advection or vorticity in models.

    The computation uses finite differences along the latitude dimension:
    :math:`\\partial F / \\partial y = (\\partial F / \\partial j) / (\\partial y / \\partial j)`,
    where :math:`j` is the grid index. Latitude is assumed in degrees and
    converted to radians.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated. If a Dataset, the gradient is applied to all
        data variables.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default, the gradient is applied over the `lat` dimension.
    min_dy: :py:class:`float <float>`, default: `1.0`.
        The minimum acceptable value of `dy` (meridional spacing in meters), below which the output
        is set to NaN to avoid numerical instabilities from very small grid spacings.
        Set to a negative value to disable this check. Unit: meters.
    edge_order: {1, 2}, optional
        Order of the finite difference used at the boundaries. 1 uses first-order accurate
        one-sided differences; 2 uses second-order accurate one-sided differences. Default: 2.
    R: :py:class:`float <float>`, default: `6370000`.
        Radius of the Earth in meters (approximate mean radius). Can be adjusted for specific models.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The meridional gradient :math:`\\partial F / \\partial y`, with the same shape and coordinates
        as the input. Units are those of the input data divided by meters (e.g., if F is in K,
        output is K/m). Invalid regions (dy < min_dy) are NaN.

    .. seealso::
        - :py:func:`calc_gradient <calc_gradient>`
        - :py:func:`calc_dlat_radian_gradient <calc_dlat_radian_gradient>`
        - :py:func:`calc_dx_gradient <calc_dx_gradient>`

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    # Set to `float64` for more accurate results in trigonometric calculations.
    lat_array = data_input[lat_dim].astype("float64")
    dlat = transfer_deg2rad(calc_gradient(lat_array, dim=lat_dim))
    dy = R * dlat
    dy = dy.where(np.abs(dy) >= min_dy)

    dFdy_raw = calc_gradient(data_input, dim=lat_dim, edge_order=edge_order)
    dFdy = dFdy_raw / dy
    dFdy = dFdy.astype("float32")

    # Add attributes for clarity
    if isinstance(dFdy, xr.DataArray):
        input_units = data_input.attrs.get("units", "dimensionless")
        dFdy.attrs.update(
            {
                "long_name": f'Meridional gradient of {data_input.name if data_input.name else "data"}',
                "units": f"{input_units} per meter",
                "standard_name": "partial_derivative_of_data_with_respect_to_meridional_distance",
            }
        )
        dFdy.name = f'd{dFdy.name or "data"}_dy' if dFdy.name else "dF_dy"
    elif isinstance(dFdy, xr.Dataset):
        for var_name, var in dFdy.data_vars.items():
            input_units = (
                data_input[var_name].attrs.get("units", "dimensionless")
                if var_name in data_input
                else "dimensionless"
            )
            var.attrs.update(
                {
                    "long_name": f"Meridional gradient of {var_name}",
                    "units": f"{input_units} per meter",
                    "standard_name": f"partial_derivative_of_{var_name}_with_respect_to_meridional_distance",
                }
            )
            dFdy[var_name].name = f"d{var_name}_dy"

    return dFdy


def calc_dlat_radian_gradient(
    data_input: xr.DataArray | xr.Dataset,
    lat_dim: str = "lat",
    edge_order: int = 2,
) -> xr.DataArray | xr.Dataset:
    """
    Calculate the gradient of the input data with respect to latitude in radians.

    This function computes the partial derivative :math:`\\partial F / \\partial \\phi`,
    where :math:`\\phi` is the latitude in radians. It is useful for spherical coordinate
    calculations, such as in wave activity flux (WAF) formulations.

    The computation uses finite differences:

    .. math::

        \\frac{\\partial F}{\\partial \\phi} = \\frac{\\partial F}{\\partial j} / \\frac{\\partial \\phi}{\\partial j},

    where :math:`j` is the grid index along the latitude dimension. Latitude values are
    assumed to be in degrees initially and converted to radians for the denominator.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated. If a Dataset, the gradient is applied to all
        data variables.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default, the gradient is applied over the `lat` dimension.
    edge_order: {1, 2}, optional
        Order of the finite difference used at the boundaries. 1 uses first-order accurate
        one-sided differences; 2 uses second-order accurate one-sided differences. Default: 2.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The gradient :math:`\\partial F / \\partial \\phi` (in radians), with the same shape
        and coordinates as the input. Units are inherited from the input data divided by radians
        (e.g., if F is in K, output is K/rad).

    .. seealso::
        - :py:func:`calc_gradient <calc_gradient>`
        - :py:func:`calc_dlat_degree_gradient <calc_dlat_degree_gradient>`
        - :py:func:`calc_dlon_radian_gradient <calc_dlon_radian_gradient>`

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    # Set to `float64` for more accurate results in trigonometric calculations.
    lat_array = data_input[lat_dim].astype("float64")
    dlat = transfer_deg2rad(calc_gradient(lat_array, dim=lat_dim))

    dFdy_raw = calc_gradient(data_input, dim=lat_dim, edge_order=edge_order)
    dFdy = dFdy_raw / dlat
    dFdy = dFdy.astype("float32")

    # Add attributes for clarity
    if isinstance(dFdy, xr.DataArray):
        dFdy.attrs.update(
            {
                "long_name": f'Gradient of {data_input.name if data_input.name else "data"} with respect to latitude (radians)',
                "units": f'{data_input.attrs.get("units", "dimensionless")} per radian',
                "standard_name": "partial_derivative_of_data_with_respect_to_latitude_radians",
            }
        )
        dFdy.name = f'd{dFdy.name or "data"}_dphi_rad' if dFdy.name else "dF_dphi_rad"
    elif isinstance(dFdy, xr.Dataset):
        for var_name, var in dFdy.data_vars.items():
            var.attrs.update(
                {
                    "long_name": f"Gradient of {var_name} with respect to latitude (radians)",
                    "units": f'{data_input[var_name].attrs.get("units", "dimensionless")} per radian',
                    "standard_name": f"partial_derivative_of_{var_name}_with_respect_to_latitude_radians",
                }
            )
            dFdy[var_name].name = f"d{var_name}_dphi_rad"

    return dFdy.astype("float32")


def calc_dlat_degree_gradient(
    data_input: xr.DataArray | xr.Dataset,
    lat_dim: str = "lat",
    edge_order: int = 2,
) -> xr.DataArray | xr.Dataset:
    """
    Calculate the gradient of the input data with respect to latitude in degrees.

    This function computes the partial derivative :math:`\\partial F / \\partial \\phi`,
    where :math:`\\phi` is the latitude in degrees. It is suitable for general-purpose
    gradient calculations where degree units are preferred.

    The computation uses finite differences:

    .. math::

        \\frac{\\partial F}{\\partial \\phi} = \\frac{\\partial F}{\\partial j} / \\frac{\\partial \\phi}{\\partial j},

    where :math:`j` is the grid index along the latitude dimension. Latitude values remain
    in degrees for the denominator.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated. If a Dataset, the gradient is applied to all
        data variables.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default, the gradient is applied over the `lat` dimension.
    edge_order: {1, 2}, optional
        Order of the finite difference used at the boundaries. 1 uses first-order accurate
        one-sided differences; 2 uses second-order accurate one-sided differences. Default: 2.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The gradient :math:`\\partial F / \\partial \\phi` (in degrees), with the same shape
        and coordinates as the input. Units are inherited from the input data divided by degrees
        (e.g., if F is in K, output is K/deg).

    .. seealso::
        - :py:func:`calc_gradient <calc_gradient>`
        - :py:func:`calc_dlat_radian_gradient <calc_dlat_radian_gradient>`
        - :py:func:`calc_dlon_degree_gradient <calc_dlon_degree_gradient>`

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    # Set to `float64` for more accurate results in trigonometric calculations.
    lat_array = data_input[lat_dim].astype("float64")
    dlat = calc_gradient(lat_array, dim=lat_dim)

    dFdy_raw = calc_gradient(data_input, dim=lat_dim, edge_order=edge_order)
    dFdy = dFdy_raw / dlat
    dFdy = dFdy.astype("float32")

    # Add attributes for clarity
    if isinstance(dFdy, xr.DataArray):
        dFdy.attrs.update(
            {
                "long_name": f'Gradient of {data_input.name if data_input.name else "data"} with respect to latitude (degrees)',
                "units": f'{data_input.attrs.get("units", "dimensionless")} per degree',
                "standard_name": "partial_derivative_of_data_with_respect_to_latitude_degrees",
            }
        )
        dFdy.name = f'd{dFdy.name or "data"}_dphi_deg' if dFdy.name else "dF_dphi_deg"
    elif isinstance(dFdy, xr.Dataset):
        for var_name, var in dFdy.data_vars.items():
            var.attrs.update(
                {
                    "long_name": f"Gradient of {var_name} with respect to latitude (degrees)",
                    "units": f'{data_input[var_name].attrs.get("units", "dimensionless")} per degree',
                    "standard_name": f"partial_derivative_of_{var_name}_with_respect_to_latitude_degrees",
                }
            )
            dFdy[var_name].name = f"d{var_name}_dphi_deg"

    return dFdy.astype("float32")


def calc_dx_laplacian(
    data_input: xr.DataArray | xr.Dataset,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    min_dx2: float = 1e9,
    edge_order: int = 2,
    R: float = 6371200.0,
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

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
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


def calc_dy_laplacian(
    data_input: xr.DataArray | xr.Dataset,
    lat_dim: str = "lat",
    min_dy2: float = 1.0,
    edge_order: int = 2,
    R: float = 6371200.0,
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

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    lat_array = data_input[lat_dim].astype("float64")

    dlat = transfer_deg2rad(calc_gradient(lat_array, dim=lat_dim))
    dy2 = (R * dlat) ** 2
    dy2 = dy2.where(dy2 >= min_dy2)

    dFdy_raw = calc_gradient(data_input, dim=lat_dim, edge_order=edge_order)
    d2Fdy2_raw = calc_gradient(dFdy_raw, dim=lat_dim, edge_order=edge_order)
    d2Fdy2 = d2Fdy2_raw / dy2
    return d2Fdy2


def calc_dxdy_mixed_derivatives(
    data_input: xr.DataArray | xr.Dataset,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    min_dxdy: float = 1e10,
    edge_order: int = 2,
    R: float = 6371200.0,
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

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
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
    data_input: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
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
    # top → bottom
    data_input = data_input.sortby(vertical_dim, ascending=True)

    if isinstance(data_input.data, dask.array.core.Array):
        # The vertical coordinate dimension is set to no chunking
        data_input = data_input.chunk({vertical_dim: -1})

    # Convert the pressure unit to Pascal
    dp_base = transfer_units_coeff(vertical_dim_units, "Pa")

    dp = calc_gradient(data_input[vertical_dim], dim=vertical_dim) * dp_base
    dF_dp = calc_gradient(data_input, dim=vertical_dim) / dp

    # clean other attrs
    if "units" in data_input.attrs:
        original_units = data_input.attrs["units"]
        dF_dp.attrs = dict()
        dF_dp.attrs["units"] = str(original_units) + " Pa^-1"
    else:
        dF_dp.attrs = dict()
        dF_dp.attrs["units"] = "[data_input units] Pa^-1"
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

    .. caution:: The units for partial derivative of `time` are :math:`\\mathrm{s^{-1}}`.

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
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
    surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
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

    Examples
    --------
    The results in :py:func:`geocat.comp.meteorology.delta_pressure <geocat-comp:geocat.comp.meteorology.delta_pressure>`:

    >>> from geocat.comp.meteorology import delta_pressure
    >>> dp = delta_pressure(
    ...     pressure_lev= np.array([1000.,925.,850.,700.,600.,500., 400.,300.,250.,200.,150.,100., 70.,50.,30.,20.,10.]),
    ...     surface_pressure = np.array([1013]),
    ... )
    >>> print(dp)
    [[ 50.5  75.  112.5 125.  100.  100.  100.   75.   50.   50.   50.   40.
       25.   20.   15.   10.    5. ]]

    For comparison, the results in :py:func:`easyclimate.calc_delta_pressure <calc_delta_pressure>`:

    >>> temp_sample = xr.DataArray(
    ...     np.array([[292.,285.,283.,277.,270.,260., 250.,235.,225.,215.,207.,207., 213.,220.,225.,228.,230.]]),
    ...     dims = ("lat", "plev"),
    ...     coords = {"plev": np.array([1000.,925.,850.,700.,600.,500., 400.,300.,250.,200.,150.,100., 70.,50.,30.,20.,10.]),
    ...             "lat": np.array([0])}
    ... )
    >>> dp = ecl.calc_delta_pressure(
    ...     data_input = temp_sample,
    ...     surface_pressure_data = xr.DataArray([1013], dims = "lat"),
    ...     vertical_dim = "plev",
    ...     surface_pressure_data_units = "Pa",
    ...     vertical_dim_units = "Pa",
    ... ).transpose("lat", "plev")
    >>> print(dp)
    <xarray.DataArray 'plev' (lat: 1, plev: 17)> Size: 136B
    array([[ 50.5,  75. , 112.5, 125. , 100. , 100. , 100. ,  75. ,  50. ,
            50. ,  50. ,  40. ,  25. ,  20. ,  15. ,  10. ,   5. ]])
    Coordinates:
    * lat      (lat) int64 8B 0
    * plev     (plev) float64 136B 1e+03 925.0 850.0 700.0 ... 50.0 30.0 20.0 10.0
    """
    # `vertical_dim` data is forced into descending order
    data_input = data_input.sortby(vertical_dim, ascending=False)

    # Change the vertical coordinate unit to `Pa`
    vertical_dim_base = transfer_units_coeff(vertical_dim_units, "Pa")
    data_input = data_input.assign_coords(
        {vertical_dim: data_input[vertical_dim] * vertical_dim_base}
    )

    # Change the surface pressure data unit to `Pa`
    surface_pressure_data_Pa = transfer_data_multiple_units(
        surface_pressure_data, surface_pressure_data_units, "Pa"
    )

    pressure_lev = xr.broadcast(data_input[vertical_dim], data_input)[0]
    pressure_top = pressure_lev.min(dim=vertical_dim)
    surface_pressure = surface_pressure_data_Pa

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

    # Add attributes to the output DataArray
    delta_pressure.attrs = {
        "units": "Pa",
        "long_name": "pressure layer thickness",
        "description": "Thickness of pressure layers in the vertical coordinate system",
        "standard_name": "pressure_thickness",
    }
    delta_pressure[vertical_dim].attrs = {"units": "Pa"}
    return delta_pressure


def calc_p_integral(
    data_input: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
    normalize: bool = False,
) -> xr.DataArray:
    """
    Calculate the vertical integral along the barometric pressure direction in the p-coordinate system.

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The spatio-temporal data to be calculated.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
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
    # Change the vertical coordinate unit to `Pa`
    vertical_dim_base = transfer_units_coeff(vertical_dim_units, "Pa")
    data_input = data_input.assign_coords(
        {vertical_dim: data_input[vertical_dim] * vertical_dim_base}
    )

    part1 = data_input.sortby(vertical_dim).integrate(
        coord=vertical_dim
    )  # unit: [data_input] *kg *s^-2

    if normalize == False:
        fint = part1

    elif normalize == True:
        part2 = (
            xr.ones_like(data_input).sortby(vertical_dim).integrate(coord=vertical_dim)
        )
        fint = part1 / part2  # unit: [data_input] *kg *s^-2

    else:
        raise ValueError("The parameter `normalize` should be `True` or `False`.")

    # Create output with coords from surface_pressure_data
    output = xr.DataArray(
        fint,
        dims=fint.dims,
        coords=fint.coords,
        attrs={
            "long_name": f"Vertical integral of {data_input.name} from top to surface",
            "units": f'{data_input.attrs.get("units", "")} Pa',
        },
        name=data_input.name,
    )
    return output


def calc_top2surface_integral(
    data_input: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    vertical_dim: str,
    surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
    method: Literal["Boer1982", "Trenberth1991", "vibeta-ncl"] = "vibeta-ncl",
    normalize: bool = False,
) -> xr.DataArray:
    """
    Calculate the vertical integral in the p-coordinate system from the ground to the zenith along the barometric pressure direction.

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The spatio-temporal data to be calculated.
    surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Surface level pressure.

    .. warning::

        This parameter only accepts local pressure (slp) and must **NOT** be substituted
        with mean sea level pressure (msl). There is a fundamental difference in their physical meaning and numerical
        characteristics—local pressure reflects atmospheric pressure at the actual elevation of local area,
        while mean sea level pressure is a theoretical value adjusted to sea level.

    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    surface_pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    method: {"Boer1982", "Trenberth1991", "vibeta-ncl"}, default: `vibeta-ncl`.
        vertical integration method. Optional values are ``Boer1982``, ``Trenberth1991`` and ``vibeta-ncl``.

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
    if method == "Boer1982":
        # https://doi.org/10.1175/1520-0493(1982)110%3C1801:DEIIC%3E2.0.CO;2

        # Change the vertical coordinate unit to `Pa`
        vertical_dim_array = transfer_data_multiple_units(
            data_input[vertical_dim], vertical_dim_units, "Pa"
        )  # Pa

        # Change the surface pressure data unit to `Pa`
        surface_pressure_data_Pa = transfer_data_multiple_units(
            surface_pressure_data, surface_pressure_data_units, "Pa"
        )  # Pa

        # Filtering the part of the air pressure level greater than the surface pressure
        level_expanded = xr.broadcast(vertical_dim_array, data_input)[0]  # Pa
        data_input_need = data_input.where(level_expanded <= surface_pressure_data_Pa)
        data_input_need = data_input_need.fillna(0)

        # Integral
        fint = calc_p_integral(
            data_input_need,
            vertical_dim=vertical_dim,
            vertical_dim_units=vertical_dim_units,
            normalize=normalize,
        )

        if normalize == True:
            integral_type = "integral average"
        elif normalize == False:
            integral_type = "integral summation"

        # Create output with coords from surface_pressure_data
        output = xr.DataArray(
            fint,
            dims=surface_pressure_data.dims,
            coords=surface_pressure_data.coords,
            attrs={
                "long_name": f"Vertical integral of {data_input.name} from top to surface using {method}",
                "units": f'{data_input.attrs.get("units", "")} Pa',
                "integral_type": f"{integral_type}",
            },
            name=data_input.name,
        )
        return output

    elif method == "Trenberth1991":
        # https://doi.org/10.1175/1520-0442(1991)004%3C0707:CDFGAC%3E2.0.CO;2
        # https://www.ncl.ucar.edu/Document/Functions/Built-in/vibeta.shtml

        dp = calc_delta_pressure(
            data_input=data_input,
            surface_pressure_data=surface_pressure_data,
            vertical_dim=vertical_dim,
            vertical_dim_units=vertical_dim_units,
            surface_pressure_data_units=surface_pressure_data_units,
        )  # Pa

        # Change the vertical coordinate unit to `Pa`
        vertical_dim_base = transfer_units_coeff(vertical_dim_units, "Pa")
        data_input_vertical_dim_Pa = data_input.assign_coords(
            {vertical_dim: data_input[vertical_dim] * vertical_dim_base}
        )  # vertical_dim unit: Pa

        # weighted sum/sum_of_layer_thickness
        part1 = (data_input_vertical_dim_Pa * dp).sum(
            dim=vertical_dim
        )  # unit: Pa kg s^-2

        if normalize == True:
            part2 = dp.sum(dim=vertical_dim)  # kg s^-2
            fint = part1 / part2  # unit: [data_input]
            integral_type = "integral average"
        elif normalize == False:
            fint = part1  # unit: [data_input] *kg *s^-2
            integral_type = "integral summation"
        else:
            raise ValueError("The parameter `normalize` should be `True` or `False`.")

        # Create output with coords from surface_pressure_data
        output = xr.DataArray(
            fint,
            dims=surface_pressure_data.dims,
            coords=surface_pressure_data.coords,
            attrs={
                "long_name": f"Vertical integral of {data_input.name} from top to surface using {method}",
                "units": f'{data_input.attrs.get("units", "")} Pa',
                "integral_type": f"{integral_type}",
            },
            name=data_input.name,
        )
        return output

    elif method == "vibeta-ncl":
        # Change the vertical coordinate unit to `Pa`
        vertical_dim_base = transfer_units_coeff(vertical_dim_units, "Pa")
        data_input = data_input.assign_coords(
            {vertical_dim: data_input[vertical_dim] * vertical_dim_base}
        )

        # Change the surface pressure data unit to `Pa`
        surface_pressure_data = transfer_data_multiple_units(
            surface_pressure_data, surface_pressure_data_units, "Pa"
        )

        # Assume vertical coordinate is pressure levels, decreasing (surface to top)
        level = data_input.coords[vertical_dim].values
        if not np.all(np.diff(level) < 0):
            raise ValueError("Vertical levels must be in decreasing order (hPa).")

        ptop = level[-1]  # top pressure (smallest)
        nlev = len(level)
        xmsg = 1e30
        linlog = 2  # Log interpolation in pressure

        def vibeta_core(x_profile, psfc):
            if psfc <= ptop:
                return np.nan

            # Replace NaNs with xmsg for Fortran handling
            x_profile = np.where(np.isnan(x_profile), xmsg, x_profile).astype(
                np.float64
            )

            xsfc = x_profile[0] if abs(x_profile[0] - xmsg) > 1e10 else 0.0
            pbot = psfc
            plvcrt = ptop

            vint, ier = dvibeta(
                level.astype(np.float64),
                x_profile.astype(np.float64),
                np.float64(xmsg),
                np.int32(linlog),
                np.float64(psfc),
                np.float64(xsfc),
                np.float64(pbot),
                np.float64(ptop),
                np.float64(plvcrt),
                np.int32(nlev),
            )

            if ier == 0:
                return vint
            else:
                return np.nan

        # Apply the core function using apply_ufunc for vectorization
        integral = xr.apply_ufunc(
            vibeta_core,
            data_input,
            surface_pressure_data,
            input_core_dims=[(vertical_dim,), ()],
            output_core_dims=[()],
            vectorize=True,
            dask="parallelized",
            output_dtypes=[np.float64],
            join="outer",
            dataset_fill_value=np.nan,
            keep_attrs=False,  # We'll set attrs manually
        )
        integral_type = "integral summation"

        if normalize == True:
            integral_ = xr.apply_ufunc(
                vibeta_core,
                xr.ones_like(data_input),
                surface_pressure_data,
                input_core_dims=[(vertical_dim,), ()],
                output_core_dims=[()],
                vectorize=True,
                dask="parallelized",
                output_dtypes=[np.float64],
                join="outer",
                dataset_fill_value=np.nan,
                keep_attrs=False,  # We'll set attrs manually
            )
            integral = integral / integral_
            integral_type = "integral average"

        # Create output with coords from surface_pressure_data
        output = xr.DataArray(
            integral,
            dims=surface_pressure_data.dims,
            coords=surface_pressure_data.coords,
            attrs={
                "long_name": f"Vertical integral of {data_input.name} from top to surface using {method}",
                "units": f'{data_input.attrs.get("units", "")} Pa',
                "integral_type": f"{integral_type}",
            },
            name=data_input.name,
        )
        return output

    else:
        raise ValueError(
            "The parameter `method` should be `Boer1982`, `Trenberth1991` or `vibeta-ncl`."
        )


def calc_dxdy_laplacian(
    data_input: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    R: float = 6371200.0,
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
    R: float = 6371200.0,
    spherical_coord=True,
    cyclic_boundary: bool = False,
    method: Literal["easyclimate", "uv2dv_cfd-ncl"] = "uv2dv_cfd-ncl",
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
    R: :py:class:`float <float>`, default: `6371200.0`.
        Radius of the Earth.
    spherical_coord: :py:class:`bool<bool>`, default: `True`.
        Whether or not to compute the horizontal Laplace term in spherical coordinates. The parameter is applicable only when ``method = easyclimate``.
    cyclic_boundary: :py:class:`bool <bool>`, default: `False`.
        If True, assume cyclic (periodic) boundaries in longitude. The parameter is applicable only when ``method = ddvfidf-ncl``.
    method: {"easyclimate", "ddvfidf-ncl"}, default: `ddvfidf-ncl`.
        The method to calculate horizontal divergence term. Optional values are ``easyclimate`` and ``ddvfidf-ncl``.

    Returns
    -------
    The horizontal divergence term. (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. seealso::

        - https://www.ncl.ucar.edu/Document/Functions/Built-in/uv2dv_cfd.shtml
        - Howard B. Bluestein. (1992). Synoptic-Dynamic Meteorology in Midlatitudes: Principles of Kinematics and Dynamics, Vol. 1. p113-114

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    compare_multi_dataarray_coordinate([u_data, v_data])

    match method:
        case "easyclimate":
            lat_array = u_data[lat_dim]
            dudx = calc_dx_gradient(u_data, lon_dim=lon_dim, lat_dim=lat_dim, R=R)
            dvdy = calc_dy_gradient(v_data, lat_dim=lat_dim, R=R)

            if spherical_coord == True:
                term3 = v_data / R * np.tan(transfer_deg2rad(lat_array))
                div = dudx + dvdy - term3
            elif spherical_coord == False:
                div = dudx + dvdy

            div.name = "divergence"
            div.attrs["long_name"] = "divergence"
            div.attrs["units"] = "s^-1"
            if "_FillValue" in u_data.attrs:
                div.attrs["_FillValue"] = u_data.attrs["_FillValue"]
            return div

        case "uv2dv_cfd-ncl":
            if lon_dim not in u_data.dims or lat_dim not in u_data.dims:
                raise ValueError(
                    f"u_data and v_data must have dimensions {lon_dim} and {lat_dim}."
                )

            missing = None
            xmsg = (
                missing
                if missing is not None
                else u_data.attrs.get("_FillValue", np.nan)
            )
            iopt = 1 if cyclic_boundary else 0

            # Transpose to ensure lat_dim is second-to-last, lon_dim is last
            transpose_order_u = [
                d for d in u_data.dims if d not in [lat_dim, lon_dim]
            ] + [lat_dim, lon_dim]
            u_trans = u_data.transpose(*transpose_order_u)
            v_trans = v_data.transpose(*transpose_order_u)  # Assume same dims

            def _core(u_vals, v_vals, lon_vals, lat_vals, xmsg):
                # Core slices: u_vals.shape = (..., nlat, mlon)
                orig_lat_vals = lat_vals.copy()
                flip_lat = orig_lat_vals[0] > orig_lat_vals[-1]
                if flip_lat:
                    u_vals = np.flip(u_vals, axis=-2)
                    v_vals = np.flip(v_vals, axis=-2)
                    lat_vals = np.flip(lat_vals)

                # Swap to (..., mlon, nlat) for Fortran
                u_vals_t = np.swapaxes(u_vals, -2, -1)
                v_vals_t = np.swapaxes(v_vals, -2, -1)
                glat_vals = lat_vals.astype(np.float64)
                glon_vals = lon_vals.astype(np.float64)

                # Handle missing values
                if np.isnan(xmsg):
                    sentinel = 1e20
                    u_fort_batch = np.nan_to_num(u_vals_t, nan=sentinel).astype(
                        np.float64
                    )
                    v_fort_batch = np.nan_to_num(v_vals_t, nan=sentinel).astype(
                        np.float64
                    )
                    xmsg_fort = float(sentinel)
                else:
                    u_fort_batch = u_vals_t.astype(np.float64)
                    v_fort_batch = v_vals_t.astype(np.float64)
                    xmsg_fort = float(xmsg)

                # Batch processing
                batch_shape = u_fort_batch.shape[:-2]
                nlat_out = u_fort_batch.shape[-1]
                mlon_out = u_fort_batch.shape[-2]
                dv_out = np.empty((*batch_shape, nlat_out, mlon_out), dtype=np.float64)

                for idx in np.ndindex(*batch_shape):
                    u_slice = u_fort_batch[idx]
                    v_slice = v_fort_batch[idx]

                    dv_slice_fort, ier = ddvfidf(
                        u_slice, v_slice, glat_vals, glon_vals, xmsg_fort, iopt
                    )
                    if ier != 0:
                        raise ValueError(
                            f"easyclimate-backend error in ddvfidf: ier={ier}"
                        )

                    # Restore missing values if using sentinel
                    if np.isnan(xmsg):
                        dv_slice_fort = np.where(
                            dv_slice_fort == xmsg_fort, np.nan, dv_slice_fort
                        )

                    # Swap back to (nlat, mlon)
                    dv_slice = np.swapaxes(dv_slice_fort, -2, -1)
                    dv_out[idx] = dv_slice

                # Flip back if original was decreasing
                if flip_lat:
                    dv_out = np.flip(dv_out, axis=-2)

                return dv_out

            div = xr.apply_ufunc(
                _core,
                u_trans,
                v_trans,
                u_trans[lon_dim],
                u_trans[lat_dim],
                xmsg,  # Pass as scalar, broadcasted
                input_core_dims=[
                    (lat_dim, lon_dim),
                    (lat_dim, lon_dim),
                    (lon_dim,),
                    (lat_dim,),
                    [],
                ],
                output_core_dims=[(lat_dim, lon_dim)],
                output_dtypes=[np.float64],
                keep_attrs=True,
                dask="allowed",  # Allow dask on non-core dimensions
                vectorize=False,
            )

            # Transpose back to original dimension order
            div = div.transpose(*u_data.dims)

            div.name = "divergence"
            div.attrs["long_name"] = "divergence"
            div.attrs["units"] = "s^-1"
            if "_FillValue" in u_data.attrs:
                div.attrs["_FillValue"] = u_data.attrs["_FillValue"]
            return div

        case _:
            raise ValueError("Method should be `easyclimate` or `uv2dv_cfd-ncl`.")


def calc_vorticity(
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    R: float = 6371200.0,
    spherical_coord: bool = True,
    cyclic_boundary: bool = False,
    method: Literal["easyclimate", "uv2vr_cfd-ncl"] = "uv2vr_cfd-ncl",
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
    cyclic_boundary: :py:class:`bool <bool>`, default: `False`.
        If True, assume cyclic (periodic) boundaries in longitude. The parameter is applicable only when ``method = uv2vr_cfd-ncl``.
    method: {"easyclimate", "uv2vr_cfd-ncl"}, default: `uv2vr_cfd-ncl`.
        The method to calculate horizontal divergence term. Optional values are ``easyclimate`` and ``uv2vr_cfd-ncl``.

    Returns
    -------
    The horizontal relative vorticity term. (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. seealso::

        - https://www.ncl.ucar.edu/Document/Functions/Built-in/uv2vr_cfd.shtml
        - Howard B. Bluestein. (1992). Synoptic-Dynamic Meteorology in Midlatitudes: Principles of Kinematics and Dynamics, Vol. 1. p113-114

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    compare_multi_dataarray_coordinate([u_data, v_data])

    match method:
        case "easyclimate":
            dvdx = calc_dx_gradient(v_data, lon_dim=lon_dim, lat_dim=lat_dim, R=R)
            dudy = calc_dy_gradient(u_data, lat_dim=lat_dim, R=R)

            if spherical_coord == True:
                term3 = u_data / R * np.tan(transfer_deg2rad(u_data[lat_dim]))
                vor = dvdx - dudy + term3
            elif spherical_coord == False:
                vor = dvdx - dudy

            vor.name = "relative_vorticity"
            vor.attrs["long_name"] = "relative_vorticity"
            vor.attrs["units"] = "s^-1"
            if "_FillValue" in u_data.attrs:
                vor.attrs["_FillValue"] = u_data.attrs["_FillValue"]
            return vor

        case "uv2vr_cfd-ncl":
            missing = None
            xmsg = (
                missing
                if missing is not None
                else u_data.attrs.get("_FillValue", np.nan)
            )
            iopt = 1 if cyclic_boundary else 0

            # Transpose to ensure lat_dim is second-to-last, lon_dim is last
            transpose_order_u = [
                d for d in u_data.dims if d not in [lat_dim, lon_dim]
            ] + [lat_dim, lon_dim]
            u_trans = u_data.transpose(*transpose_order_u)
            v_trans = v_data.transpose(*transpose_order_u)  # Assume same dims

            def _core(u_vals, v_vals, lon_vals, lat_vals, xmsg):
                # Core slices: u_vals.shape = (..., nlat, mlon)
                orig_lat_vals = lat_vals.copy()
                flip_lat = orig_lat_vals[0] > orig_lat_vals[-1]
                if flip_lat:
                    u_vals = np.flip(u_vals, axis=-2)
                    v_vals = np.flip(v_vals, axis=-2)
                    lat_vals = np.flip(lat_vals)

                # Swap to (..., mlon, nlat) for Fortran
                u_vals_t = np.swapaxes(u_vals, -2, -1)
                v_vals_t = np.swapaxes(v_vals, -2, -1)
                glat_vals = lat_vals.astype(np.float64)
                glon_vals = lon_vals.astype(np.float64)

                # Handle missing values
                if np.isnan(xmsg):
                    sentinel = 1e20
                    u_fort_batch = np.nan_to_num(u_vals_t, nan=sentinel).astype(
                        np.float64
                    )
                    v_fort_batch = np.nan_to_num(v_vals_t, nan=sentinel).astype(
                        np.float64
                    )
                    xmsg_fort = float(sentinel)
                else:
                    u_fort_batch = u_vals_t.astype(np.float64)
                    v_fort_batch = v_vals_t.astype(np.float64)
                    xmsg_fort = float(xmsg)

                # Batch processing
                batch_shape = u_fort_batch.shape[:-2]
                nlat_out = u_fort_batch.shape[-1]
                mlon_out = u_fort_batch.shape[-2]
                rv_out = np.empty((*batch_shape, nlat_out, mlon_out), dtype=np.float64)

                for idx in np.ndindex(*batch_shape):
                    u_slice = u_fort_batch[idx]
                    v_slice = v_fort_batch[idx]

                    rv_slice_fort, ier = dvrfidf(
                        u_slice, v_slice, glat_vals, glon_vals, xmsg_fort, iopt
                    )
                    if ier != 0:
                        raise ValueError(f"Fortran error in dvrfidf: ier={ier}")

                    # Restore missing values if using sentinel
                    if np.isnan(xmsg):
                        rv_slice_fort = np.where(
                            rv_slice_fort == xmsg_fort, np.nan, rv_slice_fort
                        )

                    # Swap back to (nlat, mlon)
                    rv_slice = np.swapaxes(rv_slice_fort, -2, -1)
                    rv_out[idx] = rv_slice

                # Flip back if original was decreasing
                if flip_lat:
                    rv_out = np.flip(rv_out, axis=-2)

                return rv_out

            rv = xr.apply_ufunc(
                _core,
                u_trans,
                v_trans,
                u_trans[lon_dim],
                u_trans[lat_dim],
                xmsg,  # Pass as scalar, broadcasted
                input_core_dims=[
                    (lat_dim, lon_dim),
                    (lat_dim, lon_dim),
                    (lon_dim,),
                    (lat_dim,),
                    [],
                ],
                output_core_dims=[(lat_dim, lon_dim)],
                output_dtypes=[np.float64],
                keep_attrs=True,
                dask="allowed",  # Allow dask on non-core dimensions
                vectorize=False,
            )

            # Transpose back to original dimension order
            rv = rv.transpose(*u_data.dims)

            rv.name = "relative_vorticity"
            rv.attrs["long_name"] = "relative_vorticity"
            rv.attrs["units"] = "s^-1"
            if "_FillValue" in u_data.attrs:
                rv.attrs["_FillValue"] = u_data.attrs["_FillValue"]
            return rv

        case _:
            raise ValueError("Method should be `easyclimate` or `uv2vr_cfd-ncl`.")


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
    cyclic_boundary: bool = False,
    method: Literal["easyclimate", "uv2vr_cfd-ncl"] = "uv2vr_cfd-ncl",
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
    cyclic_boundary: :py:class:`bool <bool>`, default: `False`.
        If True, assume cyclic (periodic) boundaries in longitude. The parameter is applicable only when ``method = uv2vr_cfd-ncl``.
    method: {"easyclimate", "uv2vr_cfd-ncl"}, default: `uv2vr_cfd-ncl`.
        The method to calculate horizontal divergence term. Optional values are ``easyclimate`` and ``uv2vr_cfd-ncl``.

    Returns
    -------
    The geostrophic vorticity term. (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    geostrophic_wind = calc_geostrophic_wind(
        z_data, lon_dim=lon_dim, lat_dim=lat_dim, omega=omega, g=g, R=R
    )
    ug, vg = geostrophic_wind["ug"], geostrophic_wind["vg"]
    vor_g = calc_vorticity(
        ug,
        vg,
        spherical_coord=spherical_coord,
        lon_dim=lon_dim,
        lat_dim=lat_dim,
        R=R,
        cyclic_boundary=cyclic_boundary,
        method=method,
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

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
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

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    water_flux = -omega_data * specific_humidity_data / g
    return water_flux


def calc_water_flux_top2surface_integral(
    specific_humidity_data: xr.DataArray,
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    specific_humidity_data_units: Literal["kg/kg", "g/kg", "g/g"],
    vertical_dim: str,
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
    method: Literal["Boer1982", "Trenberth1991", "vibeta-ncl"] = "vibeta-ncl",
    g: float = 9.8,
) -> xr.DataArray:
    """
    Calculate the water vapor flux across the vertical level.

    .. math::

        \\frac{1}{g} \\int_0^{p_s} (q\\mathbf{v}),dp

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
    specific_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    method: :py:class:`str <str>`, default: `'Trenberth1991'`.
        Vertical integration method. Optional values are `Boer1982`, `'Trenberth1991'`.
    g: :py:class:`float <float>`, default: `9.8`.
        The acceleration of gravity.

    Returns
    -------
    The water vapor flux. (:py:class:`xarray.Dataset<xarray.Dataset>`, :math:`\\mathrm{kg \\cdot m^-1 \\cdot s^-1 }`).

    - :math:`qu`: zonal water vapor flux.
    - :math:`qv`: meridional water vapor flux.

    .. seealso::
        :py:func:`calc_top2surface_integral <calc_top2surface_integral>`

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    specific_humidity_kg_kg = transfer_data_multiple_units(
        specific_humidity_data, specific_humidity_data_units, "kg/kg"
    )

    # Calculate the single-layer water flux
    water_flux_single_layer = calc_horizontal_water_flux(
        specific_humidity_kg_kg, u_data, v_data, g=g
    )  # 1/g *(q \mathbf{v})
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

    quv_top2surface_integral.attrs = dict()
    quv_top2surface_integral.attrs["long_name"] = (
        "Vertical integral of water vapour flux"
    )
    quv_top2surface_integral.attrs["units"] = "kg m**-1 s**-1"

    return quv_top2surface_integral


def calc_divergence_watervaporflux(
    specific_humidity_data: xr.DataArray,
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    specific_humidity_data_units: Literal["kg/kg", "g/kg", "g/g"],
    spherical_coord: bool = True,
    cyclic_boundary: bool = False,
    method: Literal["easyclimate", "uv2dv_cfd-ncl"] = "uv2dv_cfd-ncl",
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    g: float = 9.8,
    R: float = 6371200.0,
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
    specific_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
    spherical_coord: :py:class:`bool<bool>`, default: `True`.
        Whether or not to compute the horizontal Laplace term in spherical coordinates. The parameter is applicable only when ``method = easyclimate``.
    cyclic_boundary: :py:class:`bool <bool>`, default: `False`.
        If True, assume cyclic (periodic) boundaries in longitude. The parameter is applicable only when ``method = ddvfidf-ncl``.
    method: {"easyclimate", "ddvfidf-ncl"}, default: `ddvfidf-ncl`.
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

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    specific_humidity_data_kgkg = transfer_data_multiple_units(
        specific_humidity_data, specific_humidity_data_units, "kg/kg"
    )
    divergence_watervaporflux = (1 / g) * calc_divergence(
        u_data=specific_humidity_data_kgkg * u_data,
        v_data=specific_humidity_data_kgkg * v_data,
        spherical_coord=spherical_coord,
        lon_dim=lon_dim,
        lat_dim=lat_dim,
        R=R,
        method=method,
        cyclic_boundary=cyclic_boundary,
    )

    return divergence_watervaporflux


def calc_divergence_watervaporflux_top2surface_integral(
    specific_humidity_data: xr.DataArray,
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    vertical_dim: str,
    specific_humidity_data_units: Literal["kg/kg", "g/kg", "g/g"],
    surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
    spherical_coord: bool = True,
    cyclic_boundary: bool = False,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    integral_method: Literal["Boer1982", "Trenberth1991", "vibeta-ncl"] = "vibeta-ncl",
    div_method: Literal["easyclimate", "uv2dv_cfd-ncl"] = "uv2dv_cfd-ncl",
    g: float = 9.8,
    R: float = 6371200.0,
) -> xr.DataArray:
    """
    Calculate water vapor flux divergence across the vertical level.

    .. math::

        \\nabla \\cdot \\frac{1}{g} \\int_0^{p_s} (q\\mathbf{v}),dp

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
    specific_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
    surface_pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    spherical_coord: :py:class:`bool<bool>`, default: `True`.
        Whether or not to compute the horizontal Laplace term in spherical coordinates. The parameter is applicable only when ``method = easyclimate``.
    cyclic_boundary: :py:class:`bool <bool>`, default: `False`.
        If True, assume cyclic (periodic) boundaries in longitude. The parameter is applicable only when ``method = ddvfidf-ncl``.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    integral_method: {"Boer1982", "Trenberth1991", "vibeta-ncl"}, default: `vibeta-ncl`.
        The vertical integration method. Optional values are ``Boer1982``, ``Trenberth1991`` and ``vibeta-ncl``.

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

    div_method: {"easyclimate", "ddvfidf-ncl"}, default: `ddvfidf-ncl`.
        The method to calculate horizontal divergence term. Optional values are ``easyclimate`` and ``ddvfidf-ncl``.

    g: :py:class:`float <float>`, default: `9.8`.
        The acceleration of gravity.
    R: :py:class:`float <float>`, default: `6370000`.
        Radius of the Earth.

    Returns
    -------
    The water vapor flux divergence. (:py:class:`xarray.DataArray<xarray.DataArray>`, :math:`\\mathrm{kg \\cdot m^-2 \\cdot s^-1 }`).

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    specific_humidity_kg_kg = transfer_data_multiple_units(
        specific_humidity_data, specific_humidity_data_units, "kg/kg"
    )

    # Calculation water vapor flux of the whole layer
    quv = calc_water_flux_top2surface_integral(
        specific_humidity_data=specific_humidity_kg_kg,
        u_data=u_data,
        v_data=v_data,
        surface_pressure_data=surface_pressure_data,
        surface_pressure_data_units=surface_pressure_data_units,
        specific_humidity_data_units="kg/kg",
        vertical_dim=vertical_dim,
        vertical_dim_units=vertical_dim_units,
        method=integral_method,
        g=g,
    )

    # Calculation of water vapor flux divergence
    div_quv = calc_divergence(
        quv["qu"],
        quv["qv"],
        lon_dim=lon_dim,
        lat_dim=lat_dim,
        R=R,
        spherical_coord=spherical_coord,
        cyclic_boundary=cyclic_boundary,
        method=div_method,
    )

    div_quv.attrs = dict()
    div_quv.name = "wvdiv"
    div_quv.attrs["long_name"] = "Divergence of vertical integral of water vapour flux"
    div_quv.attrs["units"] = "kg m**-2 s**-1"

    result = xr.Dataset(data_vars={"qu": quv["qu"], "qv": quv["qv"], "wvdiv": div_quv})

    return result


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

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
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

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
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
    dTdp = calc_p_gradient(
        temper_data, vertical_dim=vertical_dim, vertical_dim_units=vertical_dim_units
    )
    p_adv = -omega_data * dTdp
    return p_adv


def calc_shear_stretch_deform(
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    edge_order: {1, 2} = 2,
    R: float = 6371200.0,
):
    """
    - https://www.ncl.ucar.edu/Document/Functions/Contributed/shear_stretch_deform_cfd.shtml
    - Spensberger, C., & Spengler, T. (2014). A New Look at Deformation as a Diagnostic for Large-Scale Flow. Journal of the Atmospheric Sciences, 71(11), 4221-4234. https://doi.org/10.1175/JAS-D-14-0108.1
    """
    compare_multi_dataarray_coordinate([u_data, v_data])

    lon_array = u_data[lon_dim].astype("float64")
    lat_array = u_data[lat_dim].astype("float64")
    dlon = transfer_deg2rad(calc_gradient(lon_array, dim=lon_dim))
    dlat = transfer_deg2rad(calc_gradient(lat_array, dim=lat_dim))
    coslat = np.cos(transfer_deg2rad(lat_array))

    dx = R * coslat * dlon
    dy = R * dlat

    dudx_raw = calc_gradient(u_data, dim=lon_dim, edge_order=edge_order)
    dudx = dudx_raw / dx

    dudy_raw = calc_gradient(u_data, dim=lat_dim, edge_order=edge_order)
    dudy = dudy_raw / dy

    dvdx_raw = calc_gradient(v_data, dim=lon_dim, edge_order=edge_order)
    dvdx = dvdx_raw / dx

    dvdy_raw = calc_gradient(v_data, dim=lat_dim, edge_order=edge_order)
    dvdy = dvdy_raw / dy

    shear = dvdx + dudy
    stretch = dudx - dvdy
    deform = np.sqrt(shear**2 + stretch**2)

    result = xr.Dataset(
        data_vars={"shear": shear, "stretch": stretch, "deform": deform}
    )
    return result
