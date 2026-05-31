"""
Geographic Finite Difference
"""

from __future__ import annotations

import numpy as np
from .utility import (
    find_dims_axis,
    transfer_deg2rad,
    generate_dataset_dispatcher,
    compare_multi_dataarray_coordinate,
)
from .units import (
    transfer_units_coeff,
)
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
    "calc_dxdy_laplacian",
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

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        geographic_finite_difference
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

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        geographic_finite_difference
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

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        geographic_finite_difference
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

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        geographic_finite_difference
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

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        geographic_finite_difference
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

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        geographic_finite_difference
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

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        geographic_finite_difference
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

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        geographic_finite_difference
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

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        geographic_finite_difference
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
