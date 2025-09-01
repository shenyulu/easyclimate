"""
Track cyclone center
"""

import xarray as xr
import xarray as xr
import numpy as np
import pandas as pd
from scipy.ndimage import minimum_filter
from typing import List, Tuple
from ...core.diff import calc_gradient

__all__ = ["track_cyclone_center_msl_only"]


def track_cyclone_center_msl_only(
    msl_data: xr.DataArray,
    sample_point: Tuple[float, float],
    index_value: List[int] | List[float] = [0],
    lon_dim: str = "lon",
    lat_dim: str = "lat",
) -> pd.DataFrame:
    """
    Tracks the center of a cyclone using biquadratic interpolation on mean sea level pressure (MSL) data.

    This function identifies local minima in the MSL data using a minimum filter with cyclic boundary conditions
    in longitude, selects the minimum closest to the provided sample point in geodesic distance, and applies
    biquadratic interpolation to estimate the precise location of the cyclone center. If interpolation fails,
    the grid-based minimum is used.

    This is a simple approach used by the author in `pytrack <https://github.com/tenomoto/pytrack>`__ to identify local minima in sea-level pressure.

    Local minima are identified using :py:func:`scipy.ndimage.minimum_filter() <scipy.ndimage.minimum_filter>`.
    To ensure consistency with the subsequent quadratic interpolation, we search for minima within a :math:`3 \\cdot 3` grid.
    Since the region is cropped, the ``mode`` parameter is set to ``nearest``,
    which extends the boundary values for both latitude and longitude dimensions.

    Given an estimated position, e.g., :math:`\\lambda = 140^\\circ, \\phi = 20^\\circ` in ``sample_point``,
    we calculate the great-circle distance :math:`d = a\\alpha` (:math:`a` is the Earth's radius) to the identified local minima using the formula:

    .. math::

        \\cos \\alpha = \\sin\\theta_0 \\sin\\theta + \\cos\\theta_0 \\cos\\theta \\cos(\\lambda - \\lambda_0).

    Since the relative magnitude remains unchanged when comparing the central angle :math:`\\alpha` on a unit sphere, the Earth's radius is omitted.

    The location of the local minima closest to the given longitude and latitude is identified.
    This point and its eight neighboring points, totaling nine points, are stored in a one-dimensional array.
    The center is indexed as 0, and the points are stored counterclockwise starting from the bottom-left corner.

    Let :math:`f` be a quadratic function of :math:`x` and :math:`y`:

    .. math::

        f(x, y) = c_0 + c_1x + c_2y + c_3xy + c_4x^2 + c_5y^2 + c_6x^2y + c_7xy^2 + c_8x^2y^2.

    A necessary condition for the quadratic function to have an extremum is that its gradient is zero.
    At grid points, assume :math:`f(x_0, y_0)` is a local minimum.
    The extremum of :math:`f(x, y)` may not necessarily lie on a grid point.
    Suppose the extremum of :math:`f(x, y)` is at :math:`x_0 + \\Delta x, y_0 + \\Delta y`. Expanding around :math:`x_0, y_0` using a Taylor series gives:

    .. math::

        f(x, y) = f(x_0, y_0) + f_x\\Delta x + f_y\\Delta y + \\frac{1}{2}f_{xx}(\\Delta x)^2 + \\frac{1}{2}f_{yy}(\\Delta y)^2 + f_{xy}\\Delta x\\Delta y.


    Define:

    .. math::

        \\begin{align}
        \\mathbf{b} &= -\\begin{bmatrix} f_x \\ f_y \\end{bmatrix}, \\\\
        \\mathbf{A} &= \\begin{bmatrix} f_{xx} & f_{xy} \\\ f_{xy} & f_{yy} \\end{bmatrix}, \\\\
        \\mathbf{x} &= \\begin{bmatrix} \\Delta x \\\ \\Delta y \\end{bmatrix}.
        \\end{align}

    Then:

    .. math::
        \\begin{align}
        f(x, y) &= f(x_0, y_0) - \\mathbf{b}^T\\mathbf{x} + \\frac{1}{2}\\mathbf{x}^T\\mathbf{A}\\mathbf{x}, \\\\
        \\nabla f &= \\mathbf{A}\\mathbf{x} - \\mathbf{b} = 0, \\\\
        \\mathbf{x} &= \\mathbf{A}^{-1}\\mathbf{b}.
        \\end{align}

    When :math:`d \\equiv f_{xx}f_{yy} - f_{xy}^2 \\ne 0`:

    .. math::

        \\mathbf{A}^{-1} = \\frac{1}{d} \\begin{pmatrix} f_{yy} & -f_{xy} \\\ -f_{xy} & f_{xx} \\end{pmatrix},

    allowing :math:`\\Delta x` and :math:`\\Delta y` to be determined, thus locating the extremum.

    Parameters
    ----------
    msl_data : :py:class:`xarray.DataArray<xarray.DataArray>`
        Input mean sea level pressure data with latitude and longitude dimensions.
    sample_point : :py:class:`tuple[float, float] <tuple>`
        Initial guess for the cyclone center as (longitude, latitude) in degrees.
    index_value : :py:class:`List[int] | List[float] <list>`, optional
        Index value(s) for the output DataFrame, default is [0].
    lon_dim : :py:class:`str <str>`, optional
        Name of the longitude dimension, default is 'lon'.
    lat_dim : :py:class:`str <str>`, optional
        Name of the latitude dimension, default is 'lat'.

    Returns
    -------
    :py:class:`pandas.DataFrame<pandas.DataFrame>`
        DataFrame containing the longitude, latitude, and minimum MSL pressure of the cyclone center.

    .. seealso::

        - https://github.com/tenomoto/pytrack
        - 台風: https://www.dpac.dpri.kyoto-u.ac.jp/enomoto/pymetds/Typhoon.html

    Example
    -------
    >>> import xarray as xr
    >>> import numpy as np
    >>> slp = xr.DataArray(np.random.rand(20, 30), dims=['lat', 'lon'],
    ...                    coords={'lat': np.linspace(-10, 10, 20), 'lon': np.linspace(100, 130, 30)})
    >>> result = track_cyclone_center_msl_only(slp, (110, 0), index_value = [0])

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_tc_track_axis.py
    """
    # Calculate spatial resolution for longitude and latitude
    hres_lon = np.abs(calc_gradient(msl_data[lon_dim], dim=lon_dim).data[0])
    hres_lat = np.abs(calc_gradient(msl_data[lat_dim], dim=lat_dim).data[0])

    # Identify local minima using a 3x3 minimum filter with cyclic boundary conditions
    loc_min = np.where(
        minimum_filter(msl_data, size=(3, 3), mode=("nearest", "nearest")) == msl_data
    )

    # Extract longitude and latitude coordinates
    lon = msl_data[lon_dim].data
    lat = msl_data[lat_dim].data

    # Compute geodesic distance from sample point to all local minima
    lon0, lat0 = sample_point
    dx = np.deg2rad(lon[loc_min[1]] - lon0)
    y1 = np.deg2rad(lat[loc_min[0]])
    y0 = np.deg2rad(lat0)
    d = np.arccos(np.sin(y0) * np.sin(y1) + np.cos(y0) * np.cos(y1) * np.cos(dx))

    # Find the index of the closest minimum
    n = np.argmin(d)
    jmin, imin = loc_min[0][n], loc_min[1][n]

    # Set up 9-point stencil centered on the minimum
    f = np.zeros(9)
    f[0] = msl_data[jmin, imin]
    f[1] = msl_data[jmin - 1, imin - 1]
    f[2] = msl_data[jmin - 1, imin]
    f[3] = msl_data[jmin - 1, imin + 1]
    f[4] = msl_data[jmin, imin + 1]
    f[5] = msl_data[jmin + 1, imin + 1]
    f[6] = msl_data[jmin + 1, imin]
    f[7] = msl_data[jmin + 1, imin - 1]
    f[8] = msl_data[jmin, imin - 1]

    # Calculate biquadratic interpolation coefficients
    #
    # f: one dimensional array storing stencil value
    #    in the following order
    #  y^
    #  1| 7 6 5
    #  0| 8 0 4
    # -1| 1 2 3
    #  ---------->
    #    -1 0 1 x
    c = np.zeros_like(f)
    c[0] = f[0]
    c[1] = f[4] - f[8]
    c[2] = f[6] - f[2]
    c[3] = f[1] - f[3] + f[5] - f[7]
    c[4] = 2.0 * (f[4] + f[8] - 2.0 * f[0])
    c[5] = 2.0 * (f[2] + f[6] - 2.0 * f[0])
    c[6] = 2.0 * (f[5] + f[7] - f[1] - f[3] - 2.0 * (f[6] - f[2]))
    c[7] = 2.0 * (f[3] + f[5] - f[1] - f[7] - 2.0 * (f[4] - f[8]))
    c[8] = 4.0 * (
        f[1] + f[3] + f[5] + f[7] - 2.0 * (f[2] + f[4] + f[6] + f[8]) + 4.0 * f[0]
    )

    # Compute derivatives and Hessian for interpolation
    #
    # b = (/-fx, -fy/)
    # d = fxx fyy - fxy^2
    # A = |fxx fxy|
    #     |fxy fyy|
    # A^(-1) = 1 | fyy -fxy|
    #          - |-fxy  fxx|
    #          d
    # x = A^(-1)b
    fx, fy, fxy = c[1], c[2], c[3]
    fxx, fyy = 2.0 * c[4], 2.0 * c[5]
    d = fxx * fyy - fxy * fxy

    # Perform biquadratic interpolation or use grid minimum if determinant is zero
    if d != 0:
        x = (fyy * (-fx) + (-fxy) * (-fy)) / d
        y = (-fxy * (-fx) + fxx * (-fy)) / d
        slpmin = (
            c[0]
            + c[1] * x
            + c[2] * y
            + c[3] * x * y
            + c[4] * x * x
            + c[5] * y * y
            + c[6] * x * x * y
            + c[7] * x * y * y
            + c[8] * x * x * y * y
        )
    else:
        x = 0
        y = 0
        slpmin = f[0]

    # Calculate final coordinates of the cyclone center
    lonmin = lon[imin] + x * hres_lon
    latmin = lat[jmin] - y * hres_lat

    # Return results as a DataFrame
    df = pd.DataFrame(
        {lon_dim: lonmin, lat_dim: latmin, "slp_min": slpmin}, index=index_value
    )
    return df
