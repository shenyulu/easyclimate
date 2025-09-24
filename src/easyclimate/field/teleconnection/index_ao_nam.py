"""
The Arctic Oscillation (AO)/ Monthly Northern Hemisphere Annular Mode (NAM) Index

The Arctic Oscillation (AO) Index (or Monthly Northern Hemisphere Annular Mode (NAM) Index) is a key metric used to describe
large-scale atmospheric variability in the Northern Hemisphere,
particularly influencing mid-to-high latitude weather patterns. It is defined by the leading mode of Empirical Orthogonal Function (EOF)
analysis of sea-level pressure (SLP) anomalies north of 20°N. The AO Index quantifies fluctuations
in atmospheric pressure between the Arctic and mid-latitudes, with positive and negative phases reflecting
distinct circulation patterns. In the positive phase, lower Arctic pressure and higher mid-latitude
pressure strengthen westerly winds, confining cold air to polar regions, often leading to milder
winters in North America and Europe. The negative phase, with higher Arctic pressure and weaker winds,
allows cold air to move southward, causing colder, stormier weather in these regions.

The AO's role in climate variability is significant, as it modulates temperature and precipitation, especially in winter.
The AO Index, typically derived from monthly or seasonal SLP data, reflects the strength of the polar vortex,
with positive values indicating a stronger vortex and negative values a weaker one.
It is closely linked to the North Atlantic Oscillation (NAO) due to shared variability patterns.

The AO's fluctuations are driven by internal atmospheric dynamics, stratospheric processes,
and external forcings like sea surface temperatures. Its teleconnections make it a critical factor in seasonal weather predictions and long-term climate modeling.
In a warming climate, Arctic amplification may alter AO dynamics, making its study essential for understanding future climate trends.

.. seealso::

    - Thompson, D. W. J., & Wallace, J. M. (1998). The Arctic oscillation signature in the wintertime geopotential height and temperature fields. Geophysical Research Letters, 25(9), 1297–1300. https://doi.org/10.1029/98gl00950
    - Fang, Z., Sun, X., Yang, X.-Q., & Zhu, Z. (2024). Interdecadal variations in the spatial pattern of the Arctic Oscillation Arctic center in wintertime. Geophysical Research Letters, 51, e2024GL111380. https://doi.org/10.1029/2024GL111380
    - Li, J., and J. X. L. Wang (2003), A modified zonal index and its physical sense, Geophys. Res. Lett., 30, 1632, doi: https://doi.org/10.1029/2003GL017441, 12.
    - Thompson, D. W. J. , & Wallace, J. M. . (1944). The Arctic oscillation signature in the wintertime geopotential height and temperature fields. Geophys. Res. Lett., doi: https://10.1029/98GL00950
"""

import xarray as xr
from ...core.utility import sort_ascending_latlon_coordinates
from ...core.variability import remove_seasonal_cycle_mean
from ...core.eof import get_EOF_model, calc_EOF_analysis
from typing import Literal

__all__ = [
    "calc_index_AO_EOF_Thompson_Wallace_1998",
    "calc_index_NAH_zonal_lat_Li_Wang_2003",
]


def calc_index_AO_EOF_Thompson_Wallace_1998(
    slp_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    lat_range: slice = slice(20, 90),
    time_dim: str = "time",
    random_state: int | None = None,
    solver: Literal["auto", "full", "randomized"] = "auto",
    solver_kwargs: dict = {},
    normalized: bool = True,
) -> xr.DataArray:
    """
    The calculation of monthly mean Arctic Oscillation (AO) index using empirical orthogonal functions (EOFs) method over the entire Northern Hemisphere:

    .. tip::

        EOF analysis on SLP anomalies poleward of 20°N to obtain the winter AO pattern and AO index, which is used in Thompson and Wallace (1998)

    Parameters
    ----------
    slp_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The monthly data of sea level pressure (SLP).
    time_range: :py:class:`slice <slice>`, default: `slice(None, None)`.
        The time range of seasonal cycle means to be calculated. The default value is the entire time range.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    lat_range: :py:class:`slice <slice>`, default: `slice(20, 90)`.
        The latitude range of computation using EOFs over the Northern Hemisphere. The default value is from :math:`\\mathrm{20^{\\circ}N}` to :math:`\\mathrm{90^{\\circ}N}`.
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
    The monthly mean AO index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - Thompson, D. W. J., & Wallace, J. M. (1998). The Arctic oscillation signature in the wintertime geopotential height and temperature fields. Geophysical Research Letters, 25(9), 1297–1300. https://doi.org/10.1029/98gl00950
    - Fang, Z., Sun, X., Yang, X.-Q., & Zhu, Z. (2024). Interdecadal variations in the spatial pattern of the Arctic Oscillation Arctic center in wintertime. Geophysical Research Letters, 51, e2024GL111380. https://doi.org/10.1029/2024GL111380
    - Li, J., and J. X. L. Wang (2003), A modified zonal index and its physical sense, Geophys. Res. Lett., 30, 1632, doi: https://doi.org/10.1029/2003GL017441, 12.
    - Thompson, D. W. J. , & Wallace, J. M. . (1944). The Arctic oscillation signature in the wintertime geopotential height and temperature fields. Geophys. Res. Lett., doi: https://doi.org/10.1029/98GL00950, 12.

    .. seealso::

        :py:func:`get_EOF_model <easyclimate.core.eof.get_EOF_model>`

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_ao_index.py
        ./dynamic_docs/plot_multi_linear_reg.py
    """
    slp_monthly_data = sort_ascending_latlon_coordinates(
        slp_monthly_data, lat_dim=lat_dim, lon_dim=lon_dim
    )
    # anomaly
    slp_monthly_data_NH = slp_monthly_data.sel({lat_dim: lat_range})
    slp_monthly_data_NH = remove_seasonal_cycle_mean(
        slp_monthly_data_NH, dim=time_dim, time_range=time_range
    )

    # EOF
    slp_EOF_model = get_EOF_model(
        slp_monthly_data_NH,
        lat_dim=lat_dim,
        lon_dim=lon_dim,
        time_dim=time_dim,
        n_modes=2,
        random_state=random_state,
        solver=solver,
        solver_kwargs=solver_kwargs,
    )
    slp_EOF_result = calc_EOF_analysis(slp_EOF_model, PC_normalized=normalized)
    index_AO = slp_EOF_result["PC"].sel(mode=1)
    return index_AO


def calc_index_NAH_zonal_lat_Li_Wang_2003(
    slp_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    time_dim: str = "time",
    normalized: bool = True,
) -> xr.DataArray:
    """
    The calculation of Monthly Northern Hemisphere Annular Mode (NAM) Index using normalized monthly zonal-mean sea level pressure (SLP) between 35°N and 65°N.

    .. tip::

        The monthly NAM index (NAMI) or AO index (AOI) is defined as the idfference in the normalized monthly zonal-mean sea level pressure (SLP) between 35°N and 65°N (Li and Wang, 2003)

    Parameters
    ----------
    slp_monthly_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The monthly data of sea level pressure (SLP).
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
    The monthly mean NAH/AO index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - Li, J., and J. X. L. Wang (2003), A modified zonal index and its physical sense, Geophys. Res. Lett., 30, 1632, doi: https://doi.org/10.1029/2003GL017441, 12.
    - Thompson, D. W. J. , & Wallace, J. M. . (1944). The Arctic oscillation signature in the wintertime geopotential height and temperature fields. Geophys. Res. Lett., doi: https://doi.org/10.1029/98GL00950, 12.
    - 李建平，海气耦合涛动与中国气候变化，中国气候与环境演变（上卷）（秦大河主编），北京：气象出版社，2005，324-333.  http://lijianping.cn/dct/attach/Y2xiOmNsYjpwZGY6MTk3
    - http://lijianping.cn/dct/page/65607

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_ao_index.py
    """
    slp_data_nocycle = remove_seasonal_cycle_mean(
        slp_monthly_data, dim=time_dim, time_range=time_range
    )

    index = slp_data_nocycle.sel(
        {lat_dim: 35}, method="nearest"
    ) - slp_data_nocycle.sel({lat_dim: 65}, method="nearest")
    index = index.mean(dim=lon_dim)

    # Normalized
    if normalized == True:
        index_normalized_std = index.sel({time_dim: time_range}).std(dim=time_dim).data
        result = index / index_normalized_std
        return result
    elif normalized == False:
        return index
