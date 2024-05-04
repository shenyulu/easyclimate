"""
Atlantic Niños Index
"""

import xarray as xr
from ...core.utility import sort_ascending_latlon_coordinates
from ...core.variability import remove_seasonal_cycle_mean


def calc_index_ATL3(
    sst_monthly_data: xr.DataArray | xr.Dataset,
    time_range: slice = slice(None, None),
    lat_dim: str = "lat",
    lon_dim: str = "lon",
    time_dim: str = "time",
    normalized: bool = False,
) -> xr.DataArray | xr.Dataset:
    """
    Calculate ATL3 index.

    In some years, the cold tongue formation in summer is weak, leading to warm SST anomalies in the cold tongue region.
    This is what is called an Atlantic Niño. One way to gauge the strength of these events is to calculate the
    area average of SST in the cold tongue region, defined as 20°W to 0° and 3°S to 3°N.
    This is called the ATL3 index.

    Parameters
    ----------
    sst_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The monthly sea surface temperature (SST) dataset.
    time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
        The time range of seasonal cycle means to be calculated. The default value is the entire time range.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.
    normalized: :py:class:`bool <bool>`, default `False`, optional.
        Whether to standardize the index based on standard deviation over `time_range`.

    Returns
    -------
    ATL3 index.

    Reference
    --------------
    - `Lee, S.-K., Lopez, H., Tuchen, F. P., Kim, D., Foltz, G. R., & Wittenberg, A. T. (2023). On the genesis of the 2021 Atlantic Niño. Geophysical Research Letters, 50, e2023GL104452. <https://doi.org/10.1029/2023GL104452>`__
    - Atlantic Niños. Website: https://www.jamstec.go.jp/aplinfo/climate/?page_id=1566
    """
    sst_monthly_data = sort_ascending_latlon_coordinates(
        sst_monthly_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    # anomaly
    sst_monthly_anomaly_data = remove_seasonal_cycle_mean(
        sst_monthly_data, dim=time_dim, time_range=time_range
    )

    ATL3_index = sst_monthly_anomaly_data.sel(
        {lat_dim: slice(-3, 3), lon_dim: slice(340, 360)}
    ).mean(dim=(lat_dim, lon_dim))

    # Normalized
    if normalized == True:
        index_normalized_std = (
            ATL3_index.sel({time_dim: time_range}).std(dim=time_dim).data
        )
        result = (ATL3_index / index_normalized_std).drop_vars("month")
        return result
    elif normalized == False:
        return ATL3_index
