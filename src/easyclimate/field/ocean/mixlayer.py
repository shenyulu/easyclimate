"""
The calculation of ocean mixed layer variables.
"""

from __future__ import annotations
import xarray as xr
import numpy as np
import gsw_xarray
from ...core.diff import calc_gradient, calc_u_advection, calc_v_advection
from oceans import ocfis


def calc_mixed_layer_depth(
    seawater_temperature_data: xr.DataArray,
    seawater_practical_salinity_data: xr.DataArray,
    criterion: {"temperature", "density", "pdvar"} = "pdvar",
    depth_dim: str = "depth",
    lon_dim: str = "lon",
    lat_dim: str = "lat",
) -> xr.DataArray:
    """
    Calculate the mixed layer depth.

    Parameters
    ----------
    seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{^\circ C}`)
        Vertical seawater temperature data.
    seawater_practical_salinity_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{PSU}`)
        Vertical seawater salinity data (practical salinity).
    criterion: {'temperature', 'density', 'pdvar'}, default `'pdvar'`.
        Mixed layer depth criteria

        - **temperature** : Computed based on constant temperature difference criterion, i.e., :math:`CT(0) - T[mld] = 0.5 \\mathrm{^\circ C}`.
        - **density** : Computed based on the constant potential density difference criterion, i.e., :math:`pd[0] - pd[mld] = 0.125` in sigma units.
        - **pdvar** : Computed based on variable potential density criterion :math:`pd[0] - pd[mld] = var(T[0], S[0])`, where var is a variable potential density difference which corresponds to constant temperature difference of :math:`0.5 \\mathrm{^\circ C}`.

    depth_dim: :py:class:`str <str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

    Returns
    -------
    The mixed layer depth (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. seealso::

        - https://github.com/pyoceans/oceans
        - https://pyoceans.github.io/python-oceans/ocfis.html#oceans.ocfis.mld
    """
    ds = xr.Dataset()
    ds["z"] = seawater_temperature_data[depth_dim]
    ds["lat"] = seawater_temperature_data[lat_dim]
    ds["lon"] = seawater_temperature_data[lon_dim]
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

    sea_pressure = ds["p"].broadcast_like(ds["CT"])
    absolute_salinity = ds["SA"]
    conservative_temperature = ds["CT"]

    def _calc_mld(sa, ct, p):
        if (np.isnan(sa).all() == True) or (np.isnan(ct).all() == True):
            return np.array([np.nan])
        else:
            tmp = ocfis.mld(sa, ct, p, criterion=criterion)[0]
            return np.array([tmp])

    result = xr.apply_ufunc(
        _calc_mld,
        absolute_salinity,
        conservative_temperature,
        sea_pressure,
        input_core_dims=[[depth_dim], [depth_dim], [depth_dim]],
        output_core_dims=[["parameter"]],
        output_dtypes=["float64"],
        dask="parallelized",
        vectorize=True,
        dask_gufunc_kwargs={"output_sizes": {"parameter": 1}, "allow_rechunk": True},
    )
    result = result[..., 0]

    result.attrs = {}
    result.name = "mixed_layer_depth"
    return result


def calc_MLD_depth_weighted(
    seawater_temperature_data: xr.DataArray | xr.Dataset,
    mixed_layer_depth: xr.DataArray,
    depth_dim: str = "depth",
) -> xr.DataArray | xr.Dataset:
    """
    Calculate the weights of mixed layer depth. The weights required by the ocean model under non-uniform distribution grids in the depth direction.

    Parameters
    ----------
    seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{^\circ C}`).
        Vertical seawater temperature data.
    mixed_layer_depth: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{m}`).
        The mixed layer depth.
    depth_dim: :py:class:`str <str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

    Returns
    -------
    The weights of the mixed layer depth (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).
    """
    # Use `xarray.broadcast_to` to extend the dimensionality of the mixed layer depth data by one dimension to correspond to the depth dimension of `seawater_temperature_data`
    mld_expanded = xr.broadcast(mixed_layer_depth, seawater_temperature_data)[0]
    depth_expanded = xr.broadcast(
        seawater_temperature_data[depth_dim], seawater_temperature_data
    )[0]

    # Slice the `seawater_temperature_data` data using the `xarray.where` function, keeping only the parts of the data whose depth is less than or equal to the depth of the mixed layer
    depth_mld = depth_expanded.where(depth_expanded[depth_dim] <= mld_expanded)
    depth_mld_weighted = depth_mld / depth_mld.sum(dim=depth_dim)

    return depth_mld_weighted


def calc_MLD_temper_tendency(
    seawater_temperature_anomaly_data: xr.DataArray,
    mixed_layer_depth: xr.DataArray,
    depth_weight: xr.DataArray,
    depth_dim="depth",
    time_dim="time",
) -> xr.DataArray:
    """
    Calculate the tendency of the mixing layer temperature.

    Parameters
    ----------
    seawater_temperature_anomaly_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{^\circ C}`).
        The anomaly of the vertical seawater temperature data.
    mixed_layer_depth: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{m}`).
        The mixed layer depth.
    depth_weight: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The weights of the mixed layer depth. The data is generated by the :py:class:`easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted <easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted>`
    depth_dim: :py:class:`str <str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.

    Returns
    -------
    The weights of the mixed layer depth (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    # Use `xarray.broadcast_to` to extend the dimensionality of the mixed layer depth data by one dimension to correspond to the depth dimension of `seawater_temperature_anomaly_data`
    mld_expanded = xr.broadcast(mixed_layer_depth, seawater_temperature_anomaly_data)[0]

    # Slice the `seawater_temperature_anomaly_data` data using the `xarray.where` function, keeping only the parts of the data whose depth is less than or equal to the depth of the mixed layer
    temper_mld_vertical = seawater_temperature_anomaly_data.where(
        seawater_temperature_anomaly_data[depth_dim] <= mld_expanded
    ).compute()

    # Temperature mean in the mixed layer
    mixed_layer_temperature = (temper_mld_vertical * depth_weight).sum(dim=depth_dim)

    # Calculate dT'/dt
    mixed_layer_temperature_tendency = calc_gradient(
        mixed_layer_temperature, dim=time_dim
    )

    return mixed_layer_temperature_tendency


def get_data_within_MLD(
    data_input: xr.DataArray, mixed_layer_depth: xr.DataArray, depth_dim: str = "depth"
) -> xr.DataArray:
    """
    Obtain data within the mixed layer.

    .. caution::

        This function sets the data outside the mixing layer as missing values, i.e. `np.nan`,
        but it does not calculate the average value for the data inside the mixing layer.

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
         The spatio-temporal data to be calculated.
    mixed_layer_depth: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{m}`).
        The mixed layer depth.
    depth_dim: :py:class:`str <str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

    Returns
    -------
    The data within the mixed layer (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    # Use `xarray.broadcast_to` to extend the dimensionality of the mixed layer depth data by one dimension to correspond to the depth dimension of `data_input`
    mld_expanded = xr.broadcast(mixed_layer_depth, data_input)[0]

    # Slice the `data_input` data using the `xarray.where` function, keeping only the parts of the data whose depth is less than or equal to the depth of the hybrid layer
    data_winin_mld = data_input.where(data_input[depth_dim] <= mld_expanded)

    return data_winin_mld


def get_temper_within_MLD(
    seawater_temperature_data: xr.DataArray,
    mixed_layer_depth: xr.DataArray,
    depth_dim: str = "depth",
) -> xr.DataArray:
    """
    Obtain seawater temperature data within the mixing layer.

    .. caution::

        This function sets the data outside the mixing layer as missing values, i.e. `np.nan`,
        but it does not calculate the average value for the data inside the mixing layer.

    Parameters
    ----------
    seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{^\circ C}`)
        Vertical seawater temperature data.
    mixed_layer_depth: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{m}`).
        The mixed layer depth.
    depth_dim: :py:class:`str <str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

    Returns
    -------
    The seawater temperature data within the mixed layer (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    return get_data_within_MLD(
        data_input=seawater_temperature_data,
        mixed_layer_depth=mixed_layer_depth,
        depth_dim=depth_dim,
    )


def get_data_average_within_MLD(
    data_input: xr.DataArray,
    mixed_layer_depth: xr.DataArray,
    depth_weight: xr.DataArray,
    depth_dim: str = "depth",
) -> xr.DataArray:
    """
    Obtain averaged data within the mixed layer.

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`
         The spatio-temporal data to be calculated.
    mixed_layer_depth: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{m}`).
        The mixed layer depth.
    depth_weight: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The weights of the mixed layer depth. The data is generated by the :py:class:`easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted <easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted>`
    depth_dim: :py:class:`str <str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

    Returns
    -------
    The averaged data within the mixed layer (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    # Use `xarray.broadcast_to` to extend the dimensionality of the mixed layer depth data by one dimension to correspond to the depth dimension of `data_input`
    mld_expanded = xr.broadcast(mixed_layer_depth, data_input)[0]

    # Slice the `data_input` data using the `xarray.where` function, keeping only the parts of the data whose depth is less than or equal to the depth of the hybrid layer
    data_winin_mld = data_input.where(data_input[depth_dim] <= mld_expanded)

    return (data_winin_mld * depth_weight).sum(dim=depth_dim, min_count=1)


def get_temper_average_within_MLD(
    seawater_temperature_data: xr.DataArray,
    mixed_layer_depth: xr.DataArray,
    depth_weight: xr.DataArray,
    depth_dim: str = "depth",
) -> xr.DataArray:
    """
    Obtain averaged seawater temperature data within the mixing layer.

    .. caution::

        This function sets the data outside the mixing layer as missing values, i.e. `np.nan`,
        but it does not calculate the average value for the data inside the mixing layer.

    Parameters
    ----------
    seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{^\circ C}`)
        Vertical seawater temperature data.
    mixed_layer_depth: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{m}`).
        The mixed layer depth.
    depth_weight: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The weights of the mixed layer depth. The data is generated by the :py:class:`easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted <easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted>`
    depth_dim: :py:class:`str <str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

    Returns
    -------
    The averaged seawater temperature data within the mixed layer (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    return get_data_average_within_MLD(
        data_input=seawater_temperature_data,
        mixed_layer_depth=mixed_layer_depth,
        depth_weight=depth_weight,
        depth_dim=depth_dim,
    )


def calc_MLD_average_horizontal_advection(
    u_monthly_data: xr.DataArray,
    v_monthly_data: xr.DataArray,
    seawater_temperature_data: xr.DataArray,
    mixed_layer_depth: xr.DataArray,
    depth_weight: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    depth_dim: str = "depth",
    min_dx: float = 1.0,
    min_dy: float = 1.0,
    edge_order: int = 2,
    R: float = 6370000,
) -> xr.DataArray:
    """
    Obtain the average horizontal advection within the mixed layer

    Parameters
    ----------
    u_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{m \cdot s^{-1}}`).
        The monthly ocean current data.
    v_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{m \cdot s^{-1}}`).
        The monthly meridional ocean current data.
    seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{^\circ C}`).
        Vertical seawater temperature data.
    mixed_layer_depth: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{m}`).
        The mixed layer depth.
    depth_weight: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The weights of the mixed layer depth. The data is generated by the :py:class:`easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted <easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted>`
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    depth_dim: :py:class:`str <str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.
    min_dx: :py:class:`float <float>`, default: `1.0`.
        The minimum acceptable value of `dx`, below which parts will set `nan` to avoid large computational errors.
        The unit is m. You can set it to a negative value in order to remove this benefit.
    min_dy: :py:class:`float <float>`, default: `1.0`.
        The minimum acceptable value of `dy`, below which parts will set `nan` to avoid large computational errors.
        The unit is m. You can set it to a negative value in order to remove this benefit.
    edge_order: {1, 2}, optional
        Gradient is calculated using N-th order accurate differences at the boundaries. Default: 1.
    R: :py:class:`float <float>`, default: `6370000`.
        Radius of the Earth.

    Returns
    -------
    The average horizontal advection within the mixed layer (:math:`\\mathrm{^\circ C} \\cdot \mathrm{month}^{-1}`, :py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    # Calculate $u \frac{\partial T}{\partial x}$
    u_advection = (
        calc_u_advection(
            u_data=u_monthly_data,
            temper_data=seawater_temperature_data,
            lon_dim=lon_dim,
            lat_dim=lat_dim,
            min_dx=min_dx,
            edge_order=edge_order,
            R=R,
        )
        * 2626560
        * depth_weight
    )
    u_advection = get_data_within_MLD(
        data_input=u_advection, mixed_layer_depth=mixed_layer_depth, depth_dim=depth_dim
    )
    u_advection.attrs = {}

    # Calculate $v \frac{\partial T}{\partial y}$
    v_advection = (
        calc_v_advection(
            v_data=v_monthly_data,
            temper_data=seawater_temperature_data,
            lat_dim=lat_dim,
            min_dy=min_dy,
            edge_order=edge_order,
            R=R,
        )
        * 2626560
        * depth_weight
    )
    v_advection = get_data_within_MLD(
        data_input=v_advection, mixed_layer_depth=mixed_layer_depth, depth_dim=depth_dim
    )
    v_advection.attrs = {}

    dataset = xr.Dataset()
    dataset["u_advection"] = u_advection.sum(dim=depth_dim, min_count=1)
    dataset["v_advection"] = v_advection.sum(dim=depth_dim, min_count=1)
    return dataset


def calc_MLD_average_vertical_advection(
    w_monthly_data: xr.DataArray,
    seawater_temperature_data: xr.DataArray,
    mixed_layer_depth: xr.DataArray,
    depth_weight: xr.DataArray,
    depth_dim: str = "depth",
) -> xr.DataArray:
    """
    Obtain the average vertical advection within the mixed layer.

    Parameters
    ----------
    w_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{m \cdot s^{-1}}`).
        The monthly vertical ocean current data.
    seawater_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{^\circ C}`).
        Vertical seawater temperature data.
    mixed_layer_depth: :py:class:`xarray.DataArray<xarray.DataArray>`(:math:`\\mathrm{m}`).
        The mixed layer depth.
    depth_weight: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The weights of the mixed layer depth. The data is generated by the :py:class:`easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted <easyclimate.field.ocean.mixlayer.calc_MLD_depth_weighted>`
    depth_dim: :py:class:`str <str>`, default: `depth`.
        `depth` like dimension over which to apply calculate. By default extracting is applied over the `depth` dimension.

    Returns
    -------
    The average vertical advection within the mixed layer (:math:`\\mathrm{^\circ C} \\cdot \mathrm{month}^{-1}`, :py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    # Calculate $w \frac{\partial T}{\partial z}$
    w_advection = (
        w_monthly_data
        * calc_gradient(seawater_temperature_data, dim=depth_dim)
        * 2626560
        * depth_weight
    )

    # Extracting data within the mixed layer
    w_advection_mld = get_data_within_MLD(
        data_input=w_advection, mixed_layer_depth=mixed_layer_depth, depth_dim=depth_dim
    )

    return w_advection_mld.sum(dim=depth_dim, min_count=1)


def calc_ocean_surface_heat_flux(
    qnet_monthly_anomaly_data: xr.DataArray,
    mixed_layer_depth: xr.DataArray | float,
    rho_0: float = 1000,
    c_p: float = 4000,
) -> xr.DataArray:
    """
    Obtain ocean surface heat flux.

    Parameters
    ----------
    qnet_monthly_anomaly_data: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{W \cdot m^{-2}}`).
        The monthly anomaly of the downward net heat flux at the ocean surface.
    mixed_layer_depth: :py:class:`xarray.DataArray<xarray.DataArray>` (:math:`\\mathrm{m}`).
        The mixed layer depth.
    rho_0: :py:class:`float <float>`, default: `1000` (:math:`\\mathrm{kg \cdot m^{-3}}`).
        The density of water.
    c_p: :py:class:`float <float>`, default: `4000` (:math:`\\mathrm{J \cdot kg \cdot K^{-1}}`).
        The specific heat of water.

    Returns
    -------
    The ocean surface heat flux (:math:`\\mathrm{^\circ C} \\cdot \mathrm{month}^{-1}`, :py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    Nnamchi, H., Li, J., Kucharski, F. et al. Thermodynamic controls of the Atlantic Niño. Nat Commun 6, 8895 (2015). https://doi.org/10.1038/ncomms9895
    """
    # Conversion unit: W/m^2 -> degree/month
    qnet_anomaly_degreepermon = qnet_monthly_anomaly_data * 2592000

    # Calculate $\frac{q_\mathrm{net}}{\rho_0 c_p MLD}$
    heat_flux_anomaly = qnet_anomaly_degreepermon / (rho_0 * c_p * mixed_layer_depth)

    return heat_flux_anomaly
