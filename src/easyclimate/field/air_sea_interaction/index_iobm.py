"""
Indian Ocean Basin mode (IOBM) Index
"""

import xarray as xr
from ...core.utility import sort_ascending_latlon_coordinates
from ...core.variability import remove_seasonal_cycle_mean
from ...core.eof import get_EOF_model, calc_EOF_analysis

__all__ = ["calc_index_IOBM_1point", "calc_index_IOBM_EOF1"]


def calc_index_IOBM_1point(
    sst_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_range: slice = slice(40, 110),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    time_dim: str = "time",
    normalized: bool = False,
) -> xr.DataArray:
    """
    The calculation of monthly mean Indian Ocean Basin mode (IOBM) index is constructed by following method:

        The SSTA averaged over the tropical Indian Ocean (40°E-100°E, 20°S-20°N) or (40°E-110°E, 20°S-20°N).

    Parameters
    ----------
    sst_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The monthly sea surface temperature (SST) dataset.
    time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
        The time range of seasonal cycle means to be calculated. The default value is the entire time range.
    lon_range: :py:class:`slice <slice>`, default: `slice(40, 110)`.
        The range of longitude to calculate the IOBM index. Common choices include `slice(40, 110)` and `slice(40, 100)`.
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
    The monthly mean IOBM index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - (40°E-100°E, 20°S-20°N)
        1. `Wu, Y., Tang, Y. Seasonal predictability of the tropical Indian Ocean SST in the North American multimodel ensemble. Clim Dyn 53, 3361–3372 (2019). <https://doi.org/10.1007/s00382-019-04709-0>`__
        1. `Xie, S., Hu, K., Hafner, J., Tokinaga, H., Du, Y., Huang, G., & Sampe, T. (2009). Indian Ocean Capacitor Effect on Indo–Western Pacific Climate during the Summer following El Niño. Journal of Climate, 22(3), 730-747. <https://doi.org/10.1175/2008JCLI2544.1>`__
        1. `Wang, J., Liu, Y., Cheng, F., Song, C., Li, Q., Ding, Y., and Xu, X (2023). Potential modulation of Indian Ocean basin mode on the interdecadal variations of summer precipitation over the East Asian monsoon boundary zone, EGUsphere [preprint]. <https://doi.org/10.5194/egusphere-2023-1529>`__
    - (40°E-110°E, 20°S-20°N)
        1. `Yang, J., Q. Liu, S.-P. Xie, Z. Liu, and L. Wu (2007), Impact of the Indian Ocean SST basin mode on the Asian summer monsoon, Geophys. Res. Lett., 34, L02708 <https://doi.org/10.1029/2006GL028571>`__
        2. `Sun, Y., Zhu, Z., Yang, Y. et al. Decadal change in the connection between the Pacific–Japan pattern and the Indian Ocean SST basin mode. Clim Dyn (2024). <https://doi.org/10.1007/s00382-024-07132-2>`__
    """
    sst_monthly_data = sort_ascending_latlon_coordinates(
        sst_monthly_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    # anomaly
    sst_anomaly_data = remove_seasonal_cycle_mean(
        sst_monthly_data, dim=time_dim, time_range=time_range
    )

    index_IOBM = sst_anomaly_data.sel(lon=lon_range, lat=slice(-20, 20)).mean(
        dim=(lat_dim, lon_dim)
    )

    # Normalized
    if normalized == True:
        index_normalized_std = (
            index_IOBM.sel({time_dim: time_range}).std(dim=time_dim).data
        )
        result = (index_IOBM / index_normalized_std).drop_vars("month")
        return result
    elif normalized == False:
        return index_IOBM


def calc_index_IOBM_EOF1(
    sst_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_range: slice = slice(40, 110),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    time_dim: str = "time",
    random_state=None,
    solver="auto",
    solver_kwargs={},
    normalized: bool = True,
) -> xr.DataArray:
    """
    The calculation of monthly mean PNA index using rotated empirical orthogonal functions (REOFs) method over the entire Northern Hemisphere:

    Parameters
    ----------
    sst_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The monthly sea surface temperature (SST) dataset.
    time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
        The time range of seasonal cycle means to be calculated. The default value is the entire time range.
    lon_range: :py:class:`slice <slice>`, default: `slice(40, 110)`.
        The range of longitude to calculate the IOBM index. Common choices include `slice(40, 110)` and `slice(40, 100)`.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.
    random_state: :py:class:`int<int>`, default `None`.
        Seed for the random number generator.
    solver: {"auto", "full", "randomized"}, default: "auto".
        Solver to use for the EOFs computation.
    solver_kwargs: :py:class:`dict<dict>`, default `{}`.
        Additional keyword arguments to be passed to the EOFs solver.
    normalized: :py:class:`bool <bool>`, default `True`, optional.
        Whether to standardize the index based on standard deviation over `time_range`.

    Returns
    -------
    The monthly mean IOBM index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - (40°E-100°E, 20°S-20°N)
        1. `Wu, Y., Tang, Y. Seasonal predictability of the tropical Indian Ocean SST in the North American multimodel ensemble. Clim Dyn 53, 3361–3372 (2019). <https://doi.org/10.1007/s00382-019-04709-0>`__
        1. `Xie, S., Hu, K., Hafner, J., Tokinaga, H., Du, Y., Huang, G., & Sampe, T. (2009). Indian Ocean Capacitor Effect on Indo–Western Pacific Climate during the Summer following El Niño. Journal of Climate, 22(3), 730-747. <https://doi.org/10.1175/2008JCLI2544.1>`__
        1. `Wang, J., Liu, Y., Cheng, F., Song, C., Li, Q., Ding, Y., and Xu, X (2023). Potential modulation of Indian Ocean basin mode on the interdecadal variations of summer precipitation over the East Asian monsoon boundary zone, EGUsphere [preprint]. <https://doi.org/10.5194/egusphere-2023-1529>`__
    - (40°E-110°E, 20°S-20°N)
        1. `Yang, J., Q. Liu, S.-P. Xie, Z. Liu, and L. Wu (2007), Impact of the Indian Ocean SST basin mode on the Asian summer monsoon, Geophys. Res. Lett., 34, L02708 <https://doi.org/10.1029/2006GL028571>`__
        2. `Sun, Y., Zhu, Z., Yang, Y. et al. Decadal change in the connection between the Pacific–Japan pattern and the Indian Ocean SST basin mode. Clim Dyn (2024). <https://doi.org/10.1007/s00382-024-07132-2>`__

    .. seealso::
        :py:func:`get_EOF_model <easyclimate.core.eof.get_EOF_model>`
    """
    sst_monthly_data = sort_ascending_latlon_coordinates(
        sst_monthly_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    # anomaly
    sst_anomaly_data = remove_seasonal_cycle_mean(
        sst_monthly_data, dim=time_dim, time_range=time_range
    )

    # EOF
    sst_EOF_model = get_EOF_model(
        sst_anomaly_data.sel(lon=lon_range, lat=slice(-20, 20)),
        lat_dim=lat_dim,
        lon_dim=lon_dim,
        time_dim=time_dim,
        n_modes=1,
        random_state=random_state,
        solver=solver,
        solver_kwargs=solver_kwargs,
    )
    sst_EOF_result = calc_EOF_analysis(sst_EOF_model)
    index_IOBM = sst_EOF_result["PC"].sel(mode=1)

    # Normalized
    if normalized == True:
        index_normalized_std = (
            index_IOBM.sel({time_dim: time_range}).std(dim=time_dim).data
        )
        result = (index_IOBM / index_normalized_std).drop_vars("month")
        return result
    elif normalized == False:
        return index_IOBM
