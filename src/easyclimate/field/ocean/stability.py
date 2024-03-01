"""
The calculation of ocean instability.
"""

from __future__ import annotations
import xarray as xr
import numpy as np
import gsw_xarray
from ...core.utility import find_dims_axis


def calc_N2_from_temp_salt(
    seawater_temperature_data: xr.DataArray,
    seawater_practical_salinity_data: xr.DataArray,
    time_dim: str | None,
    depth_dim: str = "depth",
    lat_dim: str = "lat",
    lon_dim: str = "lon",
) -> xr.Dataset:
    """
    Calculate the frequency of seawater buoyancy.

    Parameters
    ----------
    seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{^\circ C}`)
        Vertical seawater temperature data.
    seawater_practical_salinity_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{PSU}`)
        Vertical seawater salinity data (practical salinity).
    time_dim: :py:class:`str <str>` or `None`, default: `time`.
        The time coordinate dimension name.
    depth_dim: :py:class:`str <str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

    Returns
    -------
    The frequency of seawater buoyancy (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. seealso::

        - http://www.teos-10.org/pubs/gsw/html/gsw_Nsquared.html
        - http://www.teos-10.org/pubs/gsw/html/gsw_contents.html
    """
    # For the convenience of subsequent broadcasting operations, it is necessary to transfer the same dimension sorting method here
    if time_dim != None:
        seawater_temperature_data = seawater_temperature_data.transpose(
            time_dim, depth_dim, lat_dim, lon_dim
        )
        seawater_practical_salinity_data = seawater_practical_salinity_data.transpose(
            time_dim, depth_dim, lat_dim, lon_dim
        )
        time_length, depth_length, lat_length, lon_length = (
            seawater_temperature_data.shape
        )
    else:
        seawater_temperature_data = seawater_temperature_data.transpose(
            depth_dim, lat_dim, lon_dim
        )
        seawater_practical_salinity_data = seawater_practical_salinity_data.transpose(
            depth_dim, lat_dim, lon_dim
        )
        depth_length, lat_length, lon_length = seawater_temperature_data.shape

    ds = xr.Dataset()
    ds["z"] = seawater_temperature_data[depth_dim]
    ds["lat"] = seawater_temperature_data[lat_dim]
    ds["lon"] = seawater_temperature_data[lon_dim]
    ds["depth"] = seawater_temperature_data[depth_dim]
    ds["SP"] = seawater_practical_salinity_data
    ds["t"] = seawater_temperature_data  # ITS-90 Temperature (Celsius)

    # Height -> seawater pressure
    ds["p"] = gsw_xarray.p_from_z(z=ds["z"] * (-1), lat=ds["lat"])

    # Practical salinity -> Absolute salinity
    ds["SA"] = gsw_xarray.SA_from_SP(
        SP=ds["SP"], p=ds["p"], lon=ds["lon"], lat=ds["lat"]
    )

    # Conservative temperature
    ds["CT"] = gsw_xarray.CT_from_t(SA=ds["SA"], t=ds["t"], p=ds["p"])

    if time_dim != None:
        p_tmp = ds["p"].depth.data
        p_tmp_new = p_tmp[np.newaxis, :, np.newaxis, np.newaxis]
        p_tmp_new1 = np.broadcast_to(
            p_tmp_new, shape=(time_length, depth_length, lat_length, lon_length)
        )
        p_needed = ds["CT"].copy(data=p_tmp_new1, deep=True).where(~np.isnan(ds["t"]))

        lat_tmp = ds["lat"].data
        lat_tmp_new = lat_tmp[np.newaxis, np.newaxis, :, np.newaxis]
        lat_tmp_new1 = np.broadcast_to(
            lat_tmp_new, shape=(time_length, depth_length, lat_length, lon_length)
        )
        lat_needed = (
            ds["CT"].copy(data=lat_tmp_new1, deep=True).where(~np.isnan(ds["t"]))
        )
    else:
        p_tmp = ds["p"].depth.data
        p_tmp_new = p_tmp[:, np.newaxis, np.newaxis]
        p_tmp_new1 = np.broadcast_to(
            p_tmp_new, shape=(depth_length, lat_length, lon_length)
        )
        p_needed = ds["CT"].copy(data=p_tmp_new1, deep=True).where(~np.isnan(ds["t"]))

        lat_tmp = ds["lat"].data
        lat_tmp_new = lat_tmp[np.newaxis, :, np.newaxis]
        lat_tmp_new1 = np.broadcast_to(
            lat_tmp_new, shape=(depth_length, lat_length, lon_length)
        )
        lat_needed = (
            ds["CT"].copy(data=lat_tmp_new1, deep=True).where(~np.isnan(ds["t"]))
        )

    depth_axis_num = find_dims_axis(ds["SA"], depth_dim)
    [N2, p_mid] = gsw_xarray.Nsquared(
        SA=ds["SA"], CT=ds["CT"], p=p_needed, lat=lat_needed, axis=depth_axis_num
    )

    if time_dim != None:

        N2_dataarray = xr.DataArray(
            N2,
            dims=[time_dim, depth_dim, lat_dim, lon_dim],
            coords={
                time_dim: ds["time"].data,
                depth_dim: ds["depth"].data[:-1],
                lat_dim: ds["lat"].data,
                lon_dim: ds["lon"].data,
            },
        )

        p_mid_dataarray = xr.DataArray(
            p_mid,
            dims=[time_dim, depth_dim, lat_dim, lon_dim],
            coords={
                time_dim: ds["time"].data,
                depth_dim: ds["depth"].data[:-1],
                lat_dim: ds["lat"].data,
                lon_dim: ds["lon"].data,
            },
        )
    else:
        N2_dataarray = xr.DataArray(
            N2,
            dims=[depth_dim, lat_dim, lon_dim],
            coords={
                depth_dim: ds["depth"].data[:-1],
                lat_dim: ds["lat"].data,
                lon_dim: ds["lon"].data,
            },
        )

        p_mid_dataarray = xr.DataArray(
            p_mid,
            dims=[depth_dim, lat_dim, lon_dim],
            coords={
                depth_dim: ds["depth"].data[:-1],
                lat_dim: ds["lat"].data,
                lon_dim: ds["lon"].data,
            },
        )

    Nsquared = xr.Dataset()
    Nsquared["N2"] = N2_dataarray
    Nsquared["N2"].attrs["name"] = "Brunt-Vaisala Frequency squared (M-1xN)"
    Nsquared["N2"].attrs["units"] = "rad^2 s^-2"

    Nsquared["p_mid"] = p_mid_dataarray
    Nsquared["p_mid"].attrs["name"] = "mid pressure between p grid (M-1xN)"
    Nsquared["p_mid"].attrs["units"] = "dbar"

    Nsquared.attrs["Attention"] = (
        "The units of N2 are radians^2 s^-2 however in may textbooks this is abreviated to s^-2 as radians does not have a unit. To convert the frequency to hertz, cycles sec^-1, divide the frequency by 2π, ie N/(2π)."
    )
    return Nsquared


def calc_potential_density_from_temp_salt(
    seawater_temperature_data: xr.DataArray,
    seawater_practical_salinity_data: xr.DataArray,
    time_dim: str | None,
    depth_dim: str = "depth",
    lat_dim: str = "lat",
    lon_dim: str = "lon",
) -> xr.Dataset:
    """
    Calculate the potential density of seawater.

    Parameters
    ----------
    seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{^\circ C}`)
        Vertical seawater temperature data.
    seawater_practical_salinity_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{PSU}`)
        Vertical seawater salinity data (practical salinity).
    time_dim: :py:class:`str <str>` or `None`, default: `time`.
        The time coordinate dimension name.
    depth_dim: :py:class:`str <str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

    Returns
    -------
    The potential density of seawater (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. seealso::

        - http://www.teos-10.org/pubs/gsw/html/gsw_Nsquared.html
        - http://www.teos-10.org/pubs/gsw/html/gsw_contents.html
    """
    # For the convenience of subsequent broadcasting operations, it is necessary to transfer the same dimension sorting method here
    if time_dim != None:
        seawater_temperature_data = seawater_temperature_data.transpose(
            time_dim, depth_dim, lat_dim, lon_dim
        )
        seawater_practical_salinity_data = seawater_practical_salinity_data.transpose(
            time_dim, depth_dim, lat_dim, lon_dim
        )
        time_length, depth_length, lat_length, lon_length = (
            seawater_temperature_data.shape
        )
    else:
        seawater_temperature_data = seawater_temperature_data.transpose(
            depth_dim, lat_dim, lon_dim
        )
        seawater_practical_salinity_data = seawater_practical_salinity_data.transpose(
            depth_dim, lat_dim, lon_dim
        )
        depth_length, lat_length, lon_length = seawater_temperature_data.shape

    ds = xr.Dataset()
    ds["z"] = seawater_temperature_data[depth_dim]
    ds["lat"] = seawater_temperature_data[lat_dim]
    ds["SP"] = seawater_practical_salinity_data
    ds["t"] = seawater_temperature_data  # ITS-90 Temperature (Celsius)

    # Height -> seawater pressure
    ds["p"] = gsw_xarray.p_from_z(z=ds["z"] * (-1), lat=ds["lat"])

    # Practical salinity -> Absolute salinity
    ds["SA"] = gsw_xarray.SA_from_SP(
        SP=ds["SP"], p=ds["p"], lon=ds["lon"], lat=ds["lat"]
    )

    # Conservative temperature
    ds["CT"] = gsw_xarray.CT_from_t(SA=ds["SA"], t=ds["t"], p=ds["p"])

    if time_dim != None:
        p_tmp = ds["p"].depth.data
        p_tmp_new = p_tmp[np.newaxis, :, np.newaxis, np.newaxis]
        p_tmp_new1 = np.broadcast_to(
            p_tmp_new, shape=(time_length, depth_length, lat_length, lon_length)
        )
        p_needed = ds["CT"].copy(data=p_tmp_new1, deep=True).where(~np.isnan(ds["t"]))
    else:
        p_tmp = ds["p"].depth.data
        p_tmp_new = p_tmp[:, np.newaxis, np.newaxis]
        p_tmp_new1 = np.broadcast_to(
            p_tmp_new, shape=(depth_length, lat_length, lon_length)
        )
        p_needed = ds["CT"].copy(data=p_tmp_new1, deep=True).where(~np.isnan(ds["t"]))

    prho = gsw_xarray.pot_rho_t_exact(SA=ds["SA"], t=ds["t"], p=p_needed, p_ref=0)

    potential_density = xr.Dataset()
    potential_density["prho"] = prho

    potential_density["prho"].attrs[
        "name"
    ] = "Potential density (not potential density anomaly)"
    potential_density["prho"].attrs["units"] = "kg/m^3"

    potential_density.attrs["References"] = (
        "IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp."
    )
    return potential_density
