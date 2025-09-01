"""
Axisymmetric Analysis

The axisymmetric analysis was conducted by spherical orthodrome transformation (Ritchie 1987; Nakamura et al. 1997; Yamazaki 2011).
In the transformed coordinates, the cyclonic centre is relocated to the North Pole. The coordinate transform is easily
conducted in three-dimensional Cartesian coordinates.

Consider a longitude–latitude grid, where the cyclone relocated
at the North Pole (the NP coordinates hereafter). The transformation
between the original and NP coordinates is straightforward
with rotations in Cartesian coordinates. A positional vector in the
Cartesian coordinates :math:`\\pmb { r } = ( x , y , z ) ^ { \\mathrm { T } }` is obtained from the longitude–
latitude coordinates by

.. math::

    x=\\cos\\lambda\\sin\\theta

.. math::

    y=\\sin\\lambda\\sin\\theta

.. math::

    z=\\cos\\theta

Now, :math:`\\pmb { r }` in the NP coordinates is projected to the original coordinates
by following two rotations. Assuming that the cyclonic
centre is :math:`(\\lambda_c, \\theta_c)`, where :math:`\\lambda_c` and :math:`\\theta_c` are the central longitude and
colatitude, respectively:

1. Rotate :math:`-\\theta_c` around the :math:`y` axis (negative because of a clockwise
rotation in the x–z plane).
2. Rotate :math:`\\lambda_c` around the :math:`z` axis.

The two rotations are described by the following rotation matrix:

.. math::

    A _ { 1 } = \\left[ \\begin{array} { c c c } { { \\cos \\theta _ { \\mathrm { c } } } } & { { 0 } } & { { \\sin \\theta _ { \\mathrm { c } } } } \\\\ { { 0 } } & { { 1 } } & { { 0 } } \\\\ { { - \\sin \\theta _ { \\mathrm { c } } } } & { { 0 } } & { { \\cos \\theta _ { \\mathrm { c } } } } \\end{array} \\right]

.. math::

    A _ { 2 } = \\left[ \\begin{array} { c c c } { \\cos \\lambda _ { \\mathrm { c } } } & { - \\sin \\lambda _ { \\mathrm { c } } } & { 0 } \\\\ { \\sin \\lambda _ { \\mathrm { c } } } & { \\cos \\lambda _ { \\mathrm { c } } } & { 0 } \\\\ { 0 } & { 0 } & { 1 } \\end{array}\\right]

.. math::

    A = A_2 A_1 = \\begin{bmatrix}
    \\cos \\lambda_c \\cos \\theta_c & -\\sin \\lambda_c & \\cos \\lambda_c \\sin \\theta_c \\\\
    \\sin \\lambda_c \\cos \\theta_c & \\cos \\lambda_c & \\sin \\lambda_c \\sin \\theta_c \\\\
    -\\sin \\theta_c & 0 & \\cos \\theta_c
    \\end{bmatrix}

The points in the NP coordinates are transformed to the original
coordinates and are obtained by :math:`A\\pmb{ r }`. The longitude and latitude
corresponding to :math:`A\\pmb{ r }` are obtained by

.. math::

    \\lambda = { \\left\\{ \\begin{array} { l l } { \\text{arctan} \\! { \\frac { y } { x } } { \\bmod { 2 \\pi } } , } & { x \\! \\neq \\! 0 } \\\\ { 0 , } & { x \\! = \\! 0 } \\end{array} \\right. }

.. math::

    \\theta = a \\cos z.

The scalar field is interpolated at (:math:`\\lambda, \\theta`). The zonally symmetric and
asymmetric components in the transformed coordinates represent
the axially symmetric and asymmetric components, respectively.

The vectors such as winds need additional procedures. First,
the horizontal winds expressed in the Cartesian coordinates are

.. math::

    \\dot { x } = - u \\sin \\lambda - v \\cos \\lambda \\cos \\theta

.. math::

    \\dot { y } = u \\cos \\lambda - v \\sin \\lambda \\cos \\theta

.. math::

    \\dot { z } = v \\sin \\theta

and each component of a vector in the Cartesian coordinates
is interpolated as a scalar. Then, the coordinates are rotated with
:math:`A ^ { \\mathrm { T } } \\dot { x }`, where :math:`A^T = A_1 A_2` , to obtain winds in the NP coordinates.
The winds in the longitude–latitude coordinates are obtained from
those in the Cartesian coordinates by

.. math::

    u = \\dot { y } \\cos \\lambda - \\dot { x } \\sin \\lambda

.. math::

    v = \\operatorname{sgn} { ( \\dot { z } ) } \\sqrt { { ( \\dot { x } \\cos { \\lambda } + \\dot { y } \\sin { \\lambda } ) } ^ { 2 } + \\dot { z } ^ { 2 } }

The transformed zonal and meridional winds represent tangential
and radial components, respectively. In this study, the grid
spacing is uniform in longitude (:math:`\\lambda`) and in colatitude (:math:`\\theta=\\dfrac{\\pi}{2}-\\phi`),
and the meridional extent is :math:`10^\\circ` from the North Pole.

.. seealso::

    - Ritchie, H. (1987). Semi-Lagrangian Advection on a Gaussian Grid. Monthly Weather Review, 115(2), 608-619. https://journals.ametsoc.org/view/journals/mwre/115/2/1520-0493_1987_115_0608_slaoag_2_0_co_2.xml
    - Nakamura, H., Nakamura, M., & Anderson, J. L. (1997). The Role of High- and Low-Frequency Dynamics in Blocking Formation. Monthly Weather Review, 125(9), 2074-2093. https://journals.ametsoc.org/view/journals/mwre/125/9/1520-0493_1997_125_2074_trohal_2.0.co_2.xml.

    - Yamazaki, A. (山崎 哲), 2011: The maintenance mechanism of atmospheric blocking. D.S. thesis, Kyushu University (Available online at http://hdl.handle.net/2324/21709, https://doi.org/10.15017/21709).

    - Enomoto, T. (榎本 剛) (2019). Influence of the Track Forecast of Typhoon Prapiroon on the Heavy Rainfall in Western Japan in July 2018. SOLA, 15A, 66-71. https://doi.org/10.2151/sola.15A-012.

    - Nakashita, S. (中下 早織), & Enomoto, T. (2021). Factors for an Abrupt Increase in Track Forecast Error of Typhoon Hagibis (2019). SOLA, 17A(Special_Edition), 33-37. https://doi.org/10.2151/sola.17A-006.
"""

import xarray as xr
import xarray as xr
import numpy as np
from typing import Tuple
from ...core.datanode import DataNode

__all__ = ["cyclone_axisymmetric_analysis"]


def cyclone_axisymmetric_analysis(
    data_input: xr.DataArray,
    cyclone_center_point: Tuple[float, float],
    polar_lon: np.ndarray = np.arange(0, 360, 2),
    polar_lat: np.ndarray = np.arange(80, 90.1, 1),
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    vertical_dim: str = "level",
    R: float = 6371.0087714,
) -> DataNode:
    """
    Performs axisymmetric analysis of a cyclone by transforming data into a polar coordinate system centered on the cyclone.

    This function converts input data to a polar coordinate system based on the cyclone center, interpolates the data onto a polar grid,
    and decomposes it into symmetric and asymmetric components. The symmetric component is the azimuthal mean, and the asymmetric component
    is the deviation from this mean.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        Input data with latitude, longitude, and optionally vertical dimensions, i.e., ``(lat, lon)`` or ``(level, lat, lon)``.
    cyclone_center_point : Tuple[float, float]
        Cyclone center as ``(longitude, latitude)`` in degrees.
    polar_lon : :py:class:`numpy.ndarray <numpy.ndarray>`, optional
        Array of longitudinal angles in degrees for the polar grid, default is ``np.arange(0, 360, 2)``.
    polar_lat : :py:class:`numpy.ndarray <numpy.ndarray>`, optional
        Array of latitudinal angles in degrees for the polar grid, default is ``np.arange(80, 90.1, 1)``.
    lon_dim : :py:class:`str <str>`, optional
        Name of the longitude dimension, default is 'lon'.
    lat_dim : :py:class:`str <str>`, optional
        Name of the latitude dimension, default is 'lat'.
    vertical_dim : :py:class:`str <str>`, optional
        Name of the vertical dimension, default is 'level'.
    R : :py:class:`float <float>`, optional ( :math:`\\mathrm{km}` ).
        Earth's radius, default is 6371.0087714.

    Returns
    -------
    :py:class:`easyclimate.DataNode <easyclimate.DataNode>`
        A DataNode containing three xarray.DataArray objects:

        - rotated: Data interpolated onto the polar grid.
        - rotated_symmetric: Azimuthal mean of the rotated data.
        - rotated_asymmetric: Deviation from the azimuthal mean.

    .. seealso::

        - https://www.dpac.dpri.kyoto-u.ac.jp/enomoto/pymetds/Typhoon.html
        - Enomoto, T. (榎本 剛) (2019). Influence of the Track Forecast of Typhoon Prapiroon on the Heavy Rainfall in Western Japan in July 2018. SOLA, 15A, 66-71. https://doi.org/10.2151/sola.15A-012
        - Nakashita, S. (中下 早織), & Enomoto, T. (2021). Factors for an Abrupt Increase in Track Forecast Error of Typhoon Hagibis (2019). SOLA, 17A(Special_Edition), 33-37. https://doi.org/10.2151/sola.17A-006

    Example
    -------
    >>> import xarray as xr
    >>> import numpy as np
    >>> data = xr.DataArray(np.random.rand(37, 241, 241), dims=['level', 'lat', 'lon'],
    ...                     coords={'level': np.array([1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]),
    ...                             'lat': np.arange(60, 0-0.25, -0.25),
    ...                             'lon': np.arange(110, 170 + 0.25, 0.25)
    ... )
    >>> result = cyclone_axisymmetric_analysis(data, (140.20, 19.77))
    >>> print(result)
    <easyclimate.DataNode 'root'>
    root: /
    ├── rotated: <xarray.DataArray>
    │   Dimensions:  (level: 37, y: 11, polar_lon: 180)
    │   Coordinates:
    │     * lat        (y: 11): float64
    │     * level      (level: 37): int32
    │     * lon        (polar_lon: 180): int64
    │     * polar_lat  (y: 11): float64
    │     * polar_lon  (polar_lon: 180): int64
    │     * y          (y: 11): float64
    ├── rotated_asymmetric: <xarray.DataArray>
    │   Dimensions:  (level: 37, y: 11, polar_lon: 180)
    │   Coordinates:
    │     * lat        (y: 11): float64
    │     * level      (level: 37): int32
    │     * lon        (polar_lon: 180): int64
    │     * polar_lat  (y: 11): float64
    │     * polar_lon  (polar_lon: 180): int64
    │     * y          (y: 11): float64
    └── rotated_symmetric: <xarray.DataArray>
        Dimensions:  (level: 37, y: 11)
        Coordinates:
        * lat        (y: 11): float64
        * level      (level: 37): int32
        * polar_lat  (y: 11): float64
        * y          (y: 11): float64

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_tc_track_axis.py
    """
    # Cyclone center
    lonmin, latmin = cyclone_center_point

    # Latlon -> Cartesian coordinates
    plon = polar_lon
    plat = polar_lat
    PLON, PLAT = np.meshgrid(plon, plat)
    PLON = PLON.flatten()
    PLAT = PLAT.flatten()

    X = np.cos(np.deg2rad(PLON)) * np.cos(np.deg2rad(PLAT))
    Y = np.sin(np.deg2rad(PLON)) * np.cos(np.deg2rad(PLAT))
    Z = np.sin(np.deg2rad(PLAT))

    rlonmin, rlatmin = np.deg2rad(lonmin), np.deg2rad(latmin)
    clon, slon = np.cos(rlonmin), np.sin(rlonmin)
    clat, slat = np.cos(rlatmin), np.sin(rlatmin)
    amat = np.array(
        [
            [clon * slat, -slon, clon * clat],
            [slon * slat, clon, slon * clat],
            [-clat, 0, slat],
        ]
    )
    xyz = amat @ np.vstack([X, Y, Z])

    tlon = np.rad2deg(np.arctan2(xyz[1, :], xyz[0, :]))
    tlat = np.rad2deg(np.arcsin(xyz[2, :]))

    def linint2_points_vectorized(xi, yi, fi, xo, yo):
        # Ensure xi, yi, and fi are numpy arrays
        xi = np.asarray(xi)
        yi = np.asarray(yi)
        fi = np.asarray(fi)
        xo = np.asarray(xo)
        yo = np.asarray(yo)

        # Find indices for interpolation
        i = np.argmin(np.abs(xi - xo[:, None]), axis=1)
        i = np.where(xi[i] <= xo, i, i - 1)
        j = np.argmin(np.abs(yi - yo[:, None]), axis=1)
        j = np.where(yi[j] <= yo, j, j - 1)

        # Compute interpolation weights
        t = (xo - xi[i]) / (xi[i + 1] - xi[i])
        u = (yo - yi[j]) / (yi[j + 1] - yi[j])

        # Handle multi-dimensional input (e.g., with level)
        out_shape = fi.shape[:-2] + (xo.size,)
        fo = np.zeros(out_shape)

        # Perform bilinear interpolation
        for idx in np.ndindex(fi.shape[:-2]):
            fo[idx] = (1 - u) * (
                (1 - t) * fi[idx][j, i] + t * fi[idx][j, i + 1]
            ) + u * (t * fi[idx][j + 1, i + 1] + (1 - t) * fi[idx][j + 1, i])
        return fo

    # Prepare input coordinates and data
    xi = data_input[lon_dim].data
    yi = data_input[lat_dim].data[::-1]
    fi = (
        data_input[:, ::-1].data
        if vertical_dim in data_input.dims
        else data_input[::-1].data
    )

    # Define input and output core dimensions
    input_core_dims = (
        [lat_dim, lon_dim]
        if vertical_dim not in data_input.dims
        else [vertical_dim, lat_dim, lon_dim]
    )
    output_core_dims = ["points"]

    # Apply interpolation using apply_ufunc
    data_rotated = xr.apply_ufunc(
        linint2_points_vectorized,
        xi,
        yi,
        fi,
        tlon,
        tlat,
        input_core_dims=[[], [], input_core_dims, [], []],
        output_core_dims=[output_core_dims],
        vectorize=True,
        dask="allowed",
    )

    # Reshape the output to match desired dimensions
    if vertical_dim in data_input.dims:
        data_rotated = data_rotated.reshape(
            data_input[vertical_dim].size if vertical_dim in data_input.dims else 1,
            plat.size,
            plon.size,
        )
    else:
        data_rotated = data_rotated.reshape(plat.size, plon.size)

    # Create xarray DataArray for the rotated data
    dims = (
        [vertical_dim, "y", "polar_lon"]
        if vertical_dim in data_input.dims
        else ["y", "polar_lon"]
    )
    coords = {
        "y": R * (0.5 * np.pi - np.deg2rad(plat)),
        "polar_lon": plon,
        "polar_lat": ("y", plat),
        lat_dim: ("y", plat),
        lon_dim: ("polar_lon", plon),
    }
    if vertical_dim in data_input.dims:
        coords[vertical_dim] = data_input[vertical_dim].data

    data_rotated_xarray = xr.DataArray(data_rotated, dims=dims, coords=coords)

    node = DataNode()
    node["rotated"] = data_rotated_xarray

    data_rotated_symmetric_xarray = data_rotated_xarray.mean(dim="polar_lon")
    node["rotated_symmetric"] = data_rotated_symmetric_xarray

    data_rotated_asymmetric_xarray = data_rotated_xarray - data_rotated_symmetric_xarray
    node["rotated_asymmetric"] = data_rotated_asymmetric_xarray

    return node
