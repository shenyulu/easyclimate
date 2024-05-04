"""
Atlantic Meridional Mode (AMM) Index
"""

import xarray as xr
from ...core.utility import (
    sort_ascending_latlon_coordinates,
    transfer_xarray_lon_from360TO180,
)
from ...core.variability import remove_seasonal_cycle_mean

__all__ = [
    "calc_index_AMM_Doi_2009",
]


def calc_index_AMM_Doi_2009(
    sst_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    time_dim: str = "time",
    normalized: bool = False,
) -> xr.DataArray:
    """
    The calculation of monthly mean Atlantic Meridional Mode (AMM) index is constructed by following method:

        The difference difference between the northern index (SSTA in 5–15°N, 50–20°W) and the southern index (SSTA in 5–15°S, 20°W–10°E).
        See Fig. 7 in Doi et al. 2009.

    Parameters
    ----------
    sst_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The monthly sea surface temperature (SST) dataset.
    time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
        The time range of seasonal cycle means to be calculated. The default value is the entire time range.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.
    normalized: :py:class:`bool <bool>`, default `True`, optional.
        Whether to standardize the index based on standard deviation over `time_range`.

    Returns
    -------
    The monthly mean AMM index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - `Doi, T., Tozuka, T. & Yamagata, T. Interannual variability of the Guinea Dome and its possible link with the Atlantic Meridional Mode. Clim Dyn 33, 985–998 (2009). <https://doi.org/10.1007/s00382-009-0574-z>`__
    - `Doi, T., Tozuka, T., & Yamagata, T. (2010). The Atlantic Meridional Mode and Its Coupled Variability with the Guinea Dome. Journal of Climate, 23(2), 455-475. <https://doi.org/10.1175/2009JCLI3198.1>`__
    - `Lee, S.-K., Lopez, H., Tuchen, F. P., Kim, D., Foltz, G. R., & Wittenberg, A. T. (2023). On the genesis of the 2021 Atlantic Niño. Geophysical Research Letters, 50, e2023GL104452. <https://doi.org/10.1029/2023GL104452>`__
    """
    sst_monthly_data = sort_ascending_latlon_coordinates(
        sst_monthly_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    # anomaly
    sst_anomaly_data = remove_seasonal_cycle_mean(
        sst_monthly_data, dim=time_dim, time_range=time_range
    )
    sst_anomaly_data = transfer_xarray_lon_from360TO180(
        sst_anomaly_data, lon_dim=lon_dim
    )

    # The difference in SST anomaly between the tropical western Indian Ocean (50°E-70°E, 10°S-10°N)
    # and the tropical south-eastern Indian Ocean (90°E-110°E, 10°S-Equator)
    northern_index = sst_anomaly_data.sel(lon=slice(-50, -20), lat=slice(5, 15)).mean(
        dim=(lat_dim, lon_dim)
    )
    southern_index = sst_anomaly_data.sel(lon=slice(-20, 10), lat=slice(-15, -5)).mean(
        dim=(lat_dim, lon_dim)
    )
    index_AMM = northern_index - southern_index

    # Normalized
    if normalized == True:
        index_normalized_std = (
            index_AMM.sel({time_dim: time_range}).std(dim=time_dim).data
        )
        result = (index_AMM / index_normalized_std).drop_vars("month")
        return result
    elif normalized == False:
        return index_AMM
