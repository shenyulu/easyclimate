"""
Quantify the intensity and location of two basin-scale oceanic frontal zones
in the wintertime North Pacific, i.e. the subtropical and subarctic frontal zones (STFZ, SAFZ).
"""

from ...core.utility import get_weighted_spatial_data
from ...core.utility import sort_ascending_latlon_coordinates
import xarray as xr


def calc_intensity_STFZ(
    sst_DtDy_data: xr.DataArray,
    criteria=0.45 * 1e-5,
    lat_range=[24, 32],
    lon_range=[145, 215],
    lat_dim: str = "lat",
    lon_dim: str = "lon",
) -> xr.DataArray:
    """
    Calculate the intensity of the subtropical frontal zone (STFZ).
    The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

    .. math::
        \\mathrm{ITS} = \\sum_{i=1}^{N} \\frac{G_i}{N}

    where :math:`G_i` is the value of zonally-averaged SST meridional gradient that is no less than
    an empirically-given critical value (here, :math:`0.45 \\times 10^{-5} \\mathrm{km^{-1}}` for STFZ) at the :math:`i`-th latitudinal grid point within the zone,
    and :math:`N` is the number of total grid points that satisfy the criteria above.

    Parameters
    ----------
    sst_DtDy_data : :py:class:`xarray.DataArray<xarray.DataArray>`
        The SST meridional gradient data.
    criteria: :py:class:`float <float>`, default: `0.45 *1e-5`.
        Empirically-given critical value.
    lat_range: :py:class:`list<python.list>`, default: `[24, 32]`.
        The latitude range of the oceanic frontal zone.
    lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
        The longitude range of the oceanic frontal zone.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

    Returns
    -------
    The intensity of the STFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017),
    Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts
    and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766
    """
    sst_DtDy_data = sort_ascending_latlon_coordinates(
        sst_DtDy_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    data = sst_DtDy_data.sel(
        lon=slice(lon_range[0], lon_range[1]), lat=slice(lat_range[0], lat_range[1])
    ).mean(dim=lon_dim)
    return get_weighted_spatial_data(data.where(data > criteria)).mean(dim=lat_dim)


def calc_intensity_SAFZ(
    sst_DtDy_data: xr.DataArray,
    criteria=0.80 * 1e-5,
    lat_range=[36, 44],
    lon_range=[145, 215],
    lat_dim: str = "lat",
    lon_dim: str = "lon",
) -> xr.DataArray:
    """
    Calculate the intensity of the subarctic frontal zone (SAFZ).
    The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

    .. math::
        \\mathrm{ITS} = \\sum_{i=1}^{N} \\frac{G_i}{N}

    where :math:`G_i` is the value of zonally-averaged SST meridional gradient that is no less than
    an empirically-given critical value (here, :math:`0.80 \\times 10^{-5} \\mathrm{km^{-1}}` for SAFZ) at the :math:`i`-th latitudinal grid point within the zone,
    and :math:`N` is the number of total grid points that satisfy the criteria above.

    Parameters
    ----------
    sst_DtDy_data : :py:class:`xarray.DataArray<xarray.DataArray>`
        The SST meridional gradient data.
    criteria: :py:class:`float <float>`, default: `0.80 *1e-5`.
        Empirically-given critical value.
    lat_range: :py:class:`list<python.list>`, default: `[36, 44]`.
        The latitude range of the oceanic frontal zone.
    lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
        The longitude range of the oceanic frontal zone.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

    Returns
    -------
    The intensity of the SAFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017),
    Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts
    and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766
    """
    sst_DtDy_data = sort_ascending_latlon_coordinates(
        sst_DtDy_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    data = sst_DtDy_data.sel(
        lon=slice(lon_range[0], lon_range[1]), lat=slice(lat_range[0], lat_range[1])
    ).mean(dim=lon_dim)
    return get_weighted_spatial_data(data.where(data > criteria)).mean(dim=lat_dim)


def calc_location_STFZ(
    sst_DtDy_data: xr.DataArray,
    criteria=0.45 * 1e-5,
    lat_range=[24, 32],
    lon_range=[145, 215],
    lat_dim: str = "lat",
    lon_dim: str = "lon",
) -> xr.DataArray:
    """
    Calculate the location index of the subtropical frontal zone (STFZ).
    The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

    .. math::
        \\mathrm{LCT} = \\sum_{i=1}^{N} (G_i \\times \\mathrm{LAT}_i) / \\sum_{i=1}^{N} G_i

    where :math:`\\mathrm{LAT}_i` is the latitude at the :math:`i`-th grid point within the front zone.
    Obviously, this definition reflects a weighted-average of :math:`\\mathrm{LAT}_i` with respect to :math:`G_i`,
    indicating that the location of a front is mainly determined by larger SST meridional gradients within the frontal zone.

    Parameters
    ----------
    sst_DtDy_data : :py:class:`xarray.DataArray<xarray.DataArray>`
        The SST meridional gradient data.
    criteria: :py:class:`float <float>`, default: `0.45 *1e-5`.
        Empirically-given critical value.
    lat_range: :py:class:`list<python.list>`, default: `[24, 32]`.
        The latitude range of the oceanic frontal zone.
    lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
        The longitude range of the oceanic frontal zone.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

    Returns
    -------
    The location index of the STFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017),
    Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts
    and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766
    """
    sst_DtDy_data = sort_ascending_latlon_coordinates(
        sst_DtDy_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    data_tmp = sst_DtDy_data.sel(
        lon=slice(lon_range[0], lon_range[1]), lat=slice(lat_range[0], lat_range[1])
    ).mean(dim=lon_dim)
    data_tmp = data_tmp.where(data_tmp > criteria)
    g_i = get_weighted_spatial_data(data_tmp).sum(dim=lat_dim)
    g_i_times_lat = get_weighted_spatial_data((data_tmp * data_tmp[lat_dim])).sum(
        dim=lat_dim
    )
    lct = g_i_times_lat / g_i
    return lct


def calc_location_SAFZ(
    sst_DtDy_data: xr.DataArray,
    criteria=0.80 * 1e-5,
    lat_range=[36, 44],
    lon_range=[145, 215],
    lat_dim: str = "lat",
    lon_dim: str = "lon",
) -> xr.DataArray:
    """
    Calculate the location index of the subarctic frontal zone (SAFZ).
    The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

    .. math::
        \\mathrm{LCT} = \\sum_{i=1}^{N} (G_i \\times \\mathrm{LAT}_i) / \\sum_{i=1}^{N} G_i

    where :math:`\\mathrm{LAT}_i` is the latitude at the :math:`i`-th grid point within the front zone.
    Obviously, this definition reflects a weighted-average of :math:`\\mathrm{LAT}_i` with respect to :math:`G_i`,
    indicating that the location of a front is mainly determined by larger SST meridional gradients within the frontal zone.

    Parameters
    ----------
    sst_DtDy_data : :py:class:`xarray.DataArray<xarray.DataArray>`
        The SST meridional gradient data.
    criteria: :py:class:`float <float>`, default: `0.80 *1e-5`.
        Empirically-given critical value.
    lat_range: :py:class:`list<python.list>`, default: `[36, 44]`.
        The latitude range of the oceanic frontal zone.
    lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
        The longitude range of the oceanic frontal zone.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

    Returns
    -------
    The location index of the SAFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017),
    Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts
    and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766
    """
    sst_DtDy_data = sort_ascending_latlon_coordinates(
        sst_DtDy_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    data_tmp = sst_DtDy_data.sel(
        lon=slice(lon_range[0], lon_range[1]), lat=slice(lat_range[0], lat_range[1])
    ).mean(dim=lon_dim)
    data_tmp = data_tmp.where(data_tmp > criteria)
    g_i = get_weighted_spatial_data(data_tmp).sum(dim=lat_dim)
    g_i_times_lat = get_weighted_spatial_data((data_tmp * data_tmp[lat_dim])).sum(
        dim=lat_dim
    )
    lct = g_i_times_lat / g_i
    return lct


def calc_location_line_STFZ(
    sst_DtDy_data: xr.DataArray,
    criteria=0.45 * 1e-5,
    lat_range=[24, 32],
    lon_range=[145, 215],
    lat_dim: str = "lat",
    lon_dim: str = "lon",
) -> xr.DataArray:
    """
    Calculate the location of the subtropical frontal zone (STFZ).
    The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

    .. math::
        \\mathrm{LCT} = (\\sum_{i=1}^{N} (G_i \\times \\mathrm{LAT}_i)) / G_i

    where :math:`\\mathrm{LAT}_i` is the latitude at the :math:`i`-th grid point within the front zone.
    Obviously, this definition reflects a weighted-average of :math:`\\mathrm{LAT}_i` with respect to :math:`G_i`,
    indicating that the location of a front is mainly determined by larger SST meridional gradients within the frontal zone.

    Parameters
    ----------
    sst_DtDy_data : :py:class:`xarray.DataArray<xarray.DataArray>`
        The SST meridional gradient data.
    criteria: :py:class:`float <float>`, default: `0.45 *1e-5`.
        Empirically-given critical value.
    lat_range: :py:class:`list<python.list>`, default: `[24, 32]`.
        The latitude range of the oceanic frontal zone.
    lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
        The longitude range of the oceanic frontal zone.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

    Returns
    -------
    The location of the STFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017),
    Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts
    and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766
    """
    sst_DtDy_data = sort_ascending_latlon_coordinates(
        sst_DtDy_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    data_tmp = sst_DtDy_data.sel(
        lon=slice(lon_range[0], lon_range[1]), lat=slice(lat_range[0], lat_range[1])
    )
    data_tmp = data_tmp.where(data_tmp > criteria)
    g_i = get_weighted_spatial_data(data_tmp).sum(dim=lat_dim)
    g_i_times_lat = get_weighted_spatial_data((data_tmp * data_tmp[lat_dim])).sum(
        dim=lat_dim
    )
    lct = g_i_times_lat / g_i
    return lct


def calc_location_line_SAFZ(
    sst_DtDy_data: xr.DataArray,
    criteria=0.80 * 1e-5,
    lat_range=[36, 44],
    lon_range=[145, 215],
    lat_dim: str = "lat",
    lon_dim: str = "lon",
) -> xr.DataArray:
    """
    Calculate the location of the subarctic frontal zone (SAFZ).
    The intensity index defined reflects an average of the SST meridional gradient within a frontal zone.

    .. math::
        \\mathrm{LCT} = (\\sum_{i=1}^{N} (G_i \\times \\mathrm{LAT}_i)) / G_i

    where :math:`\\mathrm{LAT}_i` is the latitude at the :math:`i`-th grid point within the front zone.
    Obviously, this definition reflects a weighted-average of :math:`\\mathrm{LAT}_i` with respect to :math:`G_i`,
    indicating that the location of a front is mainly determined by larger SST meridional gradients within the frontal zone.

    Parameters
    ----------
    sst_DtDy_data : :py:class:`xarray.DataArray<xarray.DataArray>`
        The SST meridional gradient data.
    criteria: :py:class:`float <float>`, default: `0.80 *1e-5`.
        Empirically-given critical value.
    lat_range: :py:class:`list<python.list>`, default: `[36, 44]`.
        The latitude range of the oceanic frontal zone.
    lon_range: :py:class:`list<python.list>`, default: `[145, 215]`.
        The longitude range of the oceanic frontal zone.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

    Returns
    -------
    The location of the SAFZ (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    Wang, L., Yang, X.-Q., Yang, D., Xie, Q., Fang, J. and Sun, X. (2017),
    Two typical modes in the variabilities of wintertime North Pacific basin-scale oceanic fronts
    and associated atmospheric eddy-driven jet. Atmos. Sci. Lett, 18: 373-380. Website: https://doi.org/10.1002/asl.766
    """
    sst_DtDy_data = sort_ascending_latlon_coordinates(
        sst_DtDy_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    data_tmp = sst_DtDy_data.sel(
        lon=slice(lon_range[0], lon_range[1]), lat=slice(lat_range[0], lat_range[1])
    )
    data_tmp = data_tmp.where(data_tmp > criteria)
    g_i = get_weighted_spatial_data(data_tmp).sum(dim=lat_dim)
    g_i_times_lat = get_weighted_spatial_data((data_tmp * data_tmp[lat_dim])).sum(
        dim=lat_dim
    )
    lct = g_i_times_lat / g_i
    return lct
