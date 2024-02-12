"""
The calculation of ocean thermocline variables.
"""

from __future__ import annotations
import xarray as xr
import numpy as np
from ...core.diff import calc_gradient


def calc_seawater_thermocline_depth(
    seawater_temperature_data: xr.DataArray, depth_dim: str = "depth"
) -> xr.DataArray:
    """
    Caculate thermocline depth of ocean temperature.

    Parameters
    ----------
    seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\\mathrm{^\circ C}`).
        ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
    depth_dim: :py:class:`str <str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    # Caculate gradient
    dTdz = calc_gradient(seawater_temperature_data, dim=depth_dim)
    depth_weight = calc_gradient(seawater_temperature_data[depth_dim], dim=depth_dim)
    dTdz_dataarray = dTdz / depth_weight

    # Find the minimum gradient
    dTdz_dataarray_index = dTdz_dataarray.where(
        dTdz_dataarray == dTdz_dataarray.min(dim=depth_dim)
    )

    # Mismatch
    dTdz_dataarray_bool = dTdz_dataarray.isin(dTdz_dataarray_index)
    dTdz_dataarray_bool_plus1 = dTdz_dataarray_bool.isel(
        {depth_dim: slice(None, -1)}
    ).assign_coords({depth_dim: dTdz_dataarray_bool[depth_dim].data[1:]})
    dTdz_dataarray_bool_minus1 = dTdz_dataarray_bool.isel(
        {depth_dim: slice(1, None)}
    ).assign_coords({depth_dim: dTdz_dataarray_bool[depth_dim].data[:-1]})

    # Find the nearest point to linear interpolation
    depth_3d_array = xr.broadcast(
        seawater_temperature_data[depth_dim], seawater_temperature_data
    )[0]
    thermal_depth_1 = depth_3d_array.where(dTdz_dataarray_bool).min(dim=depth_dim)
    thermal_depth_2 = depth_3d_array.where(dTdz_dataarray_bool_plus1).min(dim=depth_dim)
    thermal_depth_3 = depth_3d_array.where(dTdz_dataarray_bool_minus1).min(
        dim=depth_dim
    )
    thermal_depth = (thermal_depth_1 + thermal_depth_2 + thermal_depth_3) / 3

    return thermal_depth


def calc_Dx_depth(
    seawater_temperature_data: xr.DataArray,
    value: float,
    depth_dim: str = "depth",
) -> xr.DataArray:
    """
    Caculate `value` depth of ocean temperature.

    Parameters
    ----------
    seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\\mathrm{^\circ C}`).
        ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
    value: float.
        The depth of ocean temperature to be calculated.
    depth_dim: :py:class:`str <str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """

    def _calc_linear_interpolate(x1, x2, y1, y2):
        return (x2 * y1 - x1 * y2) / (x2 - x1)

    data_anormaly = seawater_temperature_data - value
    abs_value = np.abs(data_anormaly)

    # Find minimum index
    min_index = abs_value.fillna(9999).argmin(dim=depth_dim, skipna=True)

    # The closest value to `value` in `seawater_temperature_data` and the depth value
    data_middle = data_anormaly.isel({depth_dim: min_index})
    data_middle_depth = data_middle[depth_dim]

    # The 2nd, 3rd closest value to `value` in `seawater_temperature_data` and the depth value
    data_up = seawater_temperature_data.isel({depth_dim: min_index - 1})
    data_down = seawater_temperature_data.isel({depth_dim: min_index + 1})

    # Find 2nd closest value to `value` in `seawater_temperature_data` and the depth value
    data_nearest = xr.concat(
        [
            data_up.assign_coords({"type": 0}).expand_dims({"type": 1}),
            data_down.assign_coords({"type": 1}).expand_dims({"type": 1}),
        ],
        dim="type",
    ).min(dim="type")
    tmp1_anormaly = xr.broadcast(seawater_temperature_data, data_nearest)[1]
    tmp2_deptharray = xr.broadcast(
        seawater_temperature_data, seawater_temperature_data[depth_dim]
    )[1]
    data_nearest_depth = tmp2_deptharray.where(
        tmp1_anormaly == seawater_temperature_data
    ).min(dim=depth_dim)

    # Linear interpolation
    isotherm_depth = _calc_linear_interpolate(
        data_nearest, data_middle, data_middle_depth, data_nearest_depth
    )

    # Create mask
    mask = (np.isnan(data_anormaly)).all(dim=depth_dim)
    mask = 1 - mask

    return isotherm_depth.where(mask)


def calc_D14_depth(
    seawater_temperature_data: xr.DataArray,
    value: float = 14,
    depth_dim: str = "depth",
) -> xr.DataArray:
    """
    Caculate 14m depth of ocean temperature.

    Parameters
    ----------
    seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\\mathrm{^\circ C}`).
        ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
    value: float (:math:`\\mathrm{m}`).
        The depth of ocean temperature to be calculated.
    depth_dim: :py:class:`str <str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return calc_Dx_depth(
        seawater_temperature_data=seawater_temperature_data,
        value=value,
        depth_dim=depth_dim,
    )


def calc_D17_depth(
    seawater_temperature_data: xr.DataArray,
    value: float = 17,
    depth_dim: str = "depth",
) -> xr.DataArray:
    """
    Caculate 17m depth of ocean temperature.

    Parameters
    ----------
    seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\\mathrm{^\circ C}`).
        ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
    value: float (:math:`\\mathrm{m}`).
        The depth of ocean temperature to be calculated.
    depth_dim: :py:class:`str <str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return calc_Dx_depth(
        seawater_temperature_data=seawater_temperature_data,
        value=value,
        depth_dim=depth_dim,
    )


def calc_D20_depth(
    seawater_temperature_data: xr.DataArray,
    value: float = 20,
    depth_dim: str = "depth",
) -> xr.DataArray:
    """
    Caculate 20m depth of ocean temperature.

    Parameters
    ----------
    seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\\mathrm{^\circ C}`).
        ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
    value: float (:math:`\\mathrm{m}`).
        The depth of ocean temperature to be calculated.
    depth_dim: :py:class:`str <str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return calc_Dx_depth(
        seawater_temperature_data=seawater_temperature_data,
        value=value,
        depth_dim=depth_dim,
    )


def calc_D26_depth(
    seawater_temperature_data: xr.DataArray,
    value: float = 26,
    depth_dim: str = "depth",
) -> xr.DataArray:
    """
    Caculate 26m depth of ocean temperature.

    Parameters
    ----------
    seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\\mathrm{^\circ C}`).
        ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
    value: float (:math:`\\mathrm{m}`).
        The depth of ocean temperature to be calculated.
    depth_dim: :py:class:`str <str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return calc_Dx_depth(
        seawater_temperature_data=seawater_temperature_data,
        value=value,
        depth_dim=depth_dim,
    )


def calc_D28_depth(
    seawater_temperature_data: xr.DataArray,
    value: float = 28,
    depth_dim: str = "depth",
) -> xr.DataArray:
    """
    Caculate 28m depth of ocean temperature.

    Parameters
    ----------
    seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`  (:math:`\\mathrm{^\circ C}`).
        ocean temperature :py:class:`xarray.DataArray<xarray.DataArray>` data to be calculated.
    value: float (:math:`\\mathrm{m}`).
        The depth of ocean temperature to be calculated.
    depth_dim: :py:class:`str <str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    return calc_Dx_depth(
        seawater_temperature_data=seawater_temperature_data,
        value=value,
        depth_dim=depth_dim,
    )
