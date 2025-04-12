"""
Barnes filter

Barnes filter is a commonly used spatial filtering method that mainly uses two constants **g** and **c** to calculate Gaussian weights,
and performs spatial interpolation for each grid point, thus becoming a low-pass filter that filters out high-frequency fluctuations.
When using two different schemes of constant g and c schemes, both retain low-frequency fluctuations of different scales.
The difference between the filtering results of the two methods can result in mesoscale fluctuations.

.. seealso::
    - Maddox, R. A. (1980). An Objective Technique for Separating Macroscale and Mesoscale Features in Meteorological Data. Monthly Weather Review, 108(8), 1108-1121. https://journals.ametsoc.org/view/journals/mwre/108/8/1520-0493_1980_108_1108_aotfsm_2_0_co_2.xml
    - https://github.com/LinOuyang/pybarnes
"""

from __future__ import annotations
import xarray as xr
from pybarnes import BarnesFilter

__all__ = [
    "calc_barnes_lowpass",
    "calc_barnes_bandpass",
]


def calc_barnes_lowpass(
    data: xr.DataArray,
    g: float = 0.3,
    c: int = 150000,
    lon_dim="lon",
    lat_dim="lat",
    radius_degree=8,
    print_progress=True,
) -> xr.DataArray:
    """
    Selecting different parameters **g** and **c**
    will result in different filtering characteristics.

    .. warning::

        **NOT** support ``data`` contains ``np.nan``.


    Parameters
    ----------
    g : :py:class:`float <float>`, generally between (0, 1]
        Constant parameter.
    c : :py:class:`int <int>`
        Constant parameter. When *c* takes a larger value, the filter function converges
        at a larger wavelength, and the response function slowly approaches the maximum value,
        which means that high-frequency fluctuations have been filtered out.
    print_progress : :py:class:`bool <bool>`
        Whether to print the progress bar when executing computation.

    Returns
    -------
    data_vars : :py:class:`xarray.DataArray <xarray.DataArray>`
        Data field after filtering out high-frequency fluctuations

    .. seealso::
        - Maddox, R. A. (1980). An Objective Technique for Separating Macroscale and Mesoscale Features in Meteorological Data. Monthly Weather Review, 108(8), 1108-1121. https://journals.ametsoc.org/view/journals/mwre/108/8/1520-0493_1980_108_1108_aotfsm_2_0_co_2.xml
        - https://github.com/LinOuyang/pybarnes
    """
    f = BarnesFilter(
        data,
        lon=data[lon_dim].data,
        lat=data[lat_dim].data,
        radius_degree=radius_degree,
    )
    result = f.lowpass(g=g, c=c, print_progress=print_progress)
    return result


def calc_barnes_bandpass(
    data: xr.DataArray,
    g1: float = 0.3,
    g2: float = 0.3,
    c1: int = 30000,
    c2: int = 150000,
    r=1.2,
    lon_dim="lon",
    lat_dim="lat",
    radius_degree=8,
    print_progress=True,
) -> xr.DataArray:
    """
    Select two different filtering schemes 1 and 2, and perform the filtering separately.
    And then perform the difference, that means **scheme1 - scheme2**.
    The mesoscale fluctuations are thus preserved.

    .. warning::

        **NOT** support ``data`` contains ``np.nan``.

    Parameters
    ----------
    g1 : :py:class:`float <float>`, generally between (0, 1]
        Constant parameter of scheme1.
    c1 : :py:class:`int <int>`
        Constant parameterof scheme1.
    g2 : :py:class:`float <float>`, generally between (0, 1]
        Constant parameter of scheme2.
    c2 : :py:class:`int <int>`
        Constant parameterof scheme2.
    r :  :py:class:`float <float>`
        The inverse of the maximum response differenc.
        It is prevented from being unduly large and very small difference fields are not greatly amplified.
    print_progress : :py:class:`bool <bool>`
        Whether to print the progress bar when executing computation.

    Returns
    -------
    data_vars : :py:class:`xarray.DataArray <xarray.DataArray>`
        Mesoscale wave field filtered out from raw data

    .. seealso::
        - Maddox, R. A. (1980). An Objective Technique for Separating Macroscale and Mesoscale Features in Meteorological Data. Monthly Weather Review, 108(8), 1108-1121. https://journals.ametsoc.org/view/journals/mwre/108/8/1520-0493_1980_108_1108_aotfsm_2_0_co_2.xml
        - https://github.com/LinOuyang/pybarnes
    """
    f = BarnesFilter(
        data,
        lon=data[lon_dim].data,
        lat=data[lat_dim].data,
        radius_degree=radius_degree,
    )
    result = f.bandpass(g1=g1, c1=c1, g2=g2, c2=c2, r=r, print_progress=print_progress)
    return result
