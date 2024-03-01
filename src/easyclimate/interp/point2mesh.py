"""
Interpolate from points to grid
"""

from __future__ import annotations

from fastbarnes.interpolation import barnes
from fastbarnes.interpolationS2 import barnes_S2
import numpy as np
import xarray as xr
import pandas as pd

__all__ = ["interp_point2mesh", "interp_point2mesh_S2"]


def interp_point2mesh(
    data: pd.DataFrame,
    var_name: str,
    point: list[int],
    grid_x: float,
    grid_y: float,
    resolution: float,
    sigma: float,
    lon_dim_name="lon",
    lat_dim_name="lat",
    method="optimized_convolution",
    num_iter=4,
    min_weight=0.001,
) -> xr.DataArray:
    """
    Computes the Barnes interpolation for observation values `var_name` taken at sample
    points `data` using Gaussian weights for the width parameter `sigma`.
    The underlying grid embedded in a Euclidean space is given with start point
    `point`, regular x-direction length `grid_x` (degree), regular y-direction length `grid_y` (degree),
    and resolution `resolution`.

    Barnes interpolation is a method that is widely used in geospatial sciences like meteorology
    to remodel data values recorded at irregularly distributed points into a representative
    analytical field. It is defined as

    .. math::
        f(\\boldsymbol{x})=\\frac{\\sum_{k=1}^N f_k\\cdot w_k(\\boldsymbol{x})}{\\sum_{k=1}^N w_k(\\boldsymbol{x})}

    with Gaussian weights

    .. math::
        w_k(\\boldsymbol{x})=\\text{e}^{-\\frac{1}{2\\sigma^2}\\left\|x-\\boldsymbol{x}_k\\right\\|^2}

    Naive computation of Barnes interpolation leads to an algorithmic complexity of :math:`O(N \\times W \\times H)`,
    where :math:`N` is the number of sample points and :math:`W \\times H` the size of the underlying grid.

    For sufficiently large :math:`n` (in general in the range from 3 to 6) a good approximation of
    Barnes interpolation with a reduced complexity :math:`O(N + W \\times H)` can be obtained by the convolutional expression

    .. math::
        f(\\boldsymbol{x})\\approx \\frac{ (\\sum_{k=1}^{N}f_k\\cdot\\delta_{\\boldsymbol{x}_k}) *  ( r_n^{*n[x]}(x)\\cdot r_n^{*n[y]}(y) )   }{ ( \\sum_{k=1}^{N} \\delta_{\\boldsymbol{x}_k}  ) *  (  r_{n}^{*n[x]}(x)\\cdot r_{n}^{*n[y]}(y)  )   }

    where :math:`\\delta` is the Dirac impulse function and :math:`r(.)` an elementary rectangular function of a specific length that depends on :math:`\\sigma` and :math:`n`.

    data : :py:class:`pandas.DataFrame<pandas.DataFrame>`
        Longitude and latitude grid point discrete data. There should be a similar structure as follows

        +------------+------------+-----------+
        |    lon     |    lat     |    qff    |
        +============+============+===========+
        |   -3.73    |   56.33    |   995.1   |
        +------------+------------+-----------+
        |    2.64    |   47.05    |  1012.5   |
        +------------+------------+-----------+
        |    ...     |   ...      |   ...     |
        +------------+------------+-----------+

        .. note::
            Data points should contain longitude (`lon`), latitude (`lat`) and data variables (the above data variable name is `qff`).

    var_name: :py:class:`str <str>`
        The name of the data variable. This should match the one in the parameter `data`.
    lat_dim_name: :py:class:`str <str>`.
        Latitude dimension name. This should match the one in the parameter `data`. By default is `lat`.
    lon_dim_name: :py:class:`str <str>`.
        Longitude dimension name. This should match the one in the parameter `data`. By default is `lon`.
    point : :py:class:`numpy.ndarray <numpy.ndarray>`.
        A 1-dimensional array of size 2 containing the coordinates of the
        start point of the grid to be used.
    grid_x : :py:class:`int <int>`.
        Length in degrees in the x-direction of the interpolated rectangular grid.
    grid_y : :py:class:`int <int>`.
        Length in degrees in the y-direction of the interpolated rectangular grid.
    resolution: :py:class:`float <float>`
        Grid resolution. The distance between regular grid points is the reciprocal of the value. Common values: 4.0, 8.0, 16.0, 32.0, 64.0.
    sigma : :py:class:`float <float>`
        The Gaussian width parameter to be used. Common values: 0.25, 0.5, 1.0, 2.0, 4.0.
    method : {'optimized_convolution', 'convolution', 'radius', 'naive'}, default: 'optimized_convolution'.
        Designates the Barnes interpolation method to be used. The possible
        implementations that can be chosen are 'naive' for the straightforward
        implementation (algorithm A from paper), 'radius' to consider only sample
        points within a specific radius of influence, both with an algorithmic
        complexity of :math:`O(N \\times W \\times H)`.
        The choice 'convolution' implements algorithm B specified in the paper
        and 'optimized_convolution' is its optimization by appending tail values
        to the rectangular kernel. The latter two algorithms reduce the complexity
        down to :math:`O(N + W \\times H)`.
        The default is 'optimized_convolution'.
    num_iter : :py:class:`int <int>`, optional
        The number of performed self-convolutions of the underlying rect-kernel.
        Applies only if method is 'optimized_convolution' or 'convolution'.
        The default is 4. Applies only to Convol interpolations: one of 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50.
    min_weight : :py:class:`float <float>`, optional
        Choose radius of influence such that Gaussian weight of considered sample
        points is greater than `min_weight`.
        Applies only if method is 'radius'. Recommended values are 0.001 and less.
        The default is 0.001, which corresponds to a radius of 3.717 * sigma.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.

    .. seealso::
        - https://github.com/MeteoSwiss/fast-barnes-py
        - Zürcher, B. K.: Fast approximate Barnes interpolation: illustrated by Python-Numba implementation fast-barnes-py v1.0, Geosci. Model Dev., 16, 1697–1711, https://doi.org/10.5194/gmd-16-1697-2023, 2023.
    """
    # definition of a (grid_x)° x (grid_y)° grid starting at point[0]°W / point[1]°N
    step = 1.0 / resolution
    x0 = np.asarray([point[0] + step, point[1]], dtype=np.float64)
    size = (int(grid_x / step), int(grid_y / step))

    lon_lat_data = data[[lon_dim_name, lat_dim_name]].to_numpy()
    qff_values = data[[var_name]].to_numpy().flatten()

    gridX = np.arange(x0[0], x0[0] + size[1] * step, step)
    gridY = np.arange(x0[1], x0[1] + size[0] * step, step)

    # calculate Barnes interpolation from fastbarnes import interpolation
    field = barnes(
        lon_lat_data,
        qff_values,
        sigma,
        x0,
        step,
        size,
        method=method,
        num_iter=num_iter,
        min_weight=min_weight,
    )

    field_dataarray = xr.DataArray(
        field,
        dims=(lat_dim_name, lon_dim_name),
        coords={lon_dim_name: gridX, lat_dim_name: gridY},
        name=var_name,
    )
    return field_dataarray


def interp_point2mesh_S2(
    data: pd.DataFrame,
    var_name: str,
    point: list[int],
    grid_x: float,
    grid_y: float,
    resolution: float,
    sigma: float,
    lon_dim_name="lon",
    lat_dim_name="lat",
    method="optimized_convolution_S2",
    num_iter=4,
    resample=True,
) -> xr.DataArray:
    """
    Computes the Barnes interpolation for observation values `var_name` taken at sample
    points `data` using Gaussian weights for the width parameter `sigma`.

    The underlying grid embedded on the unit sphere :math:`S^2` and thus inherits the
    spherical distance measure (taken in degrees). The grid is given by the start
    point `point`, regular x-direction length `grid_x` (degree), regular y-direction length `grid_y` (degree),
    and resolution `resolution`.

    Parameters
    ----------
    data : :py:class:`pandas.DataFrame<pandas.DataFrame>`
        Longitude and latitude grid point discrete data. There should be a similar structure as follows

        +------------+------------+-----------+
        |    lon     |    lat     |    qff    |
        +============+============+===========+
        |   -3.73    |   56.33    |   995.1   |
        +------------+------------+-----------+
        |    2.64    |   47.05    |  1012.5   |
        +------------+------------+-----------+
        |    ...     |   ...      |   ...     |
        +------------+------------+-----------+

        .. note::
            Data points should contain longitude (`lon`), latitude (`lat`) and data variables (the above data variable name is `qff`).

    var_name: :py:class:`str <str>`
        The name of the data variable. This should match the one in the parameter `data`.
    lat_dim_name: :py:class:`str <str>`.
        Latitude dimension name. This should match the one in the parameter `data`. By default is `lat`.
    lon_dim_name: :py:class:`str <str>`.
        Longitude dimension name. This should match the one in the parameter `data`. By default is `lon`.
    point : :py:class:`numpy.ndarray <numpy.ndarray>`
        A 1-dimensional array of size 2 containing the coordinates of the
        start point of the grid to be used.
    grid_x : :py:class:`int <int>`.
        Length in degrees in the x-direction of the interpolated rectangular grid.
    grid_y : :py:class:`int <int>`.
        Length in degrees in the y-direction of the interpolated rectangular grid.
    resolution: :py:class:`float <float>`
        Grid resolution. The distance between regular grid points is the reciprocal of the value. Common values: 4.0, 8.0, 16.0, 32.0, 64.0.
    sigma : :py:class:`float <float>`
        The Gaussian width parameter to be used. Common values: 0.25, 0.5, 1.0, 2.0, 4.0.
    method : {'optimized_convolution_S2', 'naive_S2'}, default: 'optimized_convolution_S2'.
        Designates the Barnes interpolation method to be used. The possible
        implementations that can be chosen are 'naive_S2' for the straightforward
        implementation (algorithm A from the paper) with an algorithmic complexity
        of :math:`O(N \\times W \\times H)`.
        The choice 'optimized_convolution_S2' implements the optimized algorithm B
        specified in the paper by appending tail values to the rectangular kernel.
        The latter algorithm has a reduced complexity of :math:`O(N + W \\times H)`.
        The default is 'optimized_convolution_S2'.
    num_iter : :py:class:`int <int>`, optional, default: 4.
        The number of performed self-convolutions of the underlying rect-kernel.
        Applies only if method is 'optimized_convolution_S2'.
        The default is 4. Applies only to Convol interpolations: one of 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50.
    resample : :py:class:`bool <bool>`, optional, default: `True`.
        Specifies whether to resample Lambert grid field to lonlat grid.
        Applies only if method is 'optimized_convolution_S2'.
        The default is True.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.

    .. seealso::
        - https://github.com/MeteoSwiss/fast-barnes-py
        - Zürcher, B. K.: Fast approximate Barnes interpolation: illustrated by Python-Numba implementation fast-barnes-py v1.0, Geosci. Model Dev., 16, 1697–1711, https://doi.org/10.5194/gmd-16-1697-2023, 2023.
    """
    # definition of a (grid_x)° x (grid_y)° grid starting at point[0]°W / point[1]°N
    step = 1.0 / resolution
    x0 = np.asarray([point[0] + step, point[1]], dtype=np.float64)
    size = (int(grid_x / step), int(grid_y / step))

    lon_lat_data = data[[lon_dim_name, lat_dim_name]].to_numpy()
    qff_values = data[[var_name]].to_numpy().flatten()

    gridX = np.arange(x0[0], x0[0] + size[1] * step, step)
    gridY = np.arange(x0[1], x0[1] + size[0] * step, step)

    # calculate Barnes interpolation from fastbarnes import interpolation
    field = barnes_S2(
        lon_lat_data,
        qff_values,
        sigma,
        x0,
        step,
        size,
        method=method,
        num_iter=num_iter,
        resample=resample,
    )

    field_dataarray = xr.DataArray(
        field,
        dims=(lat_dim_name, lon_dim_name),
        coords={lon_dim_name: gridX, lat_dim_name: gridY},
        name=var_name,
    )
    return field_dataarray
