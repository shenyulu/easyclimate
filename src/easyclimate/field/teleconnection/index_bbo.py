"""
Barents-Beaufort Oscillation (BBO)
"""

import xarray as xr
from typing import Literal
from easyclimate.core.utility import sort_ascending_latlon_coordinates
from easyclimate.core.variability import remove_seasonal_cycle_mean
from easyclimate.core.eof import get_EOF_model, calc_EOF_analysis

__all__ = ["calc_index_BBO_EOF3_Wu_2007"]


def calc_index_BBO_EOF3_Wu_2007(
    slp_monthly_data: xr.DataArray,
    time_range: slice = slice(None, None),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    lat_range: slice = slice(70, 90),
    time_dim: str = "time",
    random_state: int | None = None,
    solver: Literal["auto", "full", "randomized"] = "auto",
    solver_kwargs: dict = {},
    normalized: bool = True,
) -> xr.DataArray:
    """
    The calculation of monthly mean Barents-Beaufort Oscillation (BBO) index using empirical orthogonal functions (EOFs) method

    .. tip::

        The third EOF mode of SLP anomaly north of 70°N is Barents-Beaufort Oscillation (BBO) pattern.

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
    The monthly mean BBO index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - Wu, B., and M. A. Johnson (2007), A seesaw structure in SLP anomalies between the Beaufort Sea and the Barents Sea, Geophys. Res. Lett., 34, L05811, doi: https://doi.org/10.1029/2006GL028333.
    - H. Bi, K. Sun, X. Zhou, H. Huang and X. Xu, "Arctic Sea Ice Area Export Through the Fram Strait Estimated From Satellite-Based Data:1988–2012," in IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, vol. 9, no. 7, pp. 3144-3157, July 2016, doi: https://doi.org/10.1109/JSTARS.2016.2584539.
    - Bi, H., Liang, Y., & Chen, X. (2023). Distinct role of a spring atmospheric circulation mode in the Arctic sea ice decline in summer. Journal of Geophysical Research: Atmospheres, 128, e2022JD037477. https://doi.org/10.1029/2022JD037477

    .. seealso::

        :py:func:`get_EOF_model <easyclimate.core.eof.get_EOF_model>`

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_da_bbo.py
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
        n_modes=10,
        random_state=random_state,
        solver=solver,
        solver_kwargs=solver_kwargs,
    )
    slp_EOF_result = calc_EOF_analysis(slp_EOF_model, PC_normalized=normalized)
    index_BBO = slp_EOF_result["PC"].sel(mode=3)
    return index_BBO
