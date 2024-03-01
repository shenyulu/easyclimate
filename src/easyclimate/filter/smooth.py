"""
spatial smoothing
"""

from __future__ import annotations
import numpy as np
import xarray as xr
import metpy.calc as mpcalc
from ..core.utility import dequantify_metpy_xarraydata, generate_dataset_dispatcher
from typing import Sequence

__all__ = [
    "calc_spatial_smooth_gaussian",
    "calc_spatial_smooth_rectangular",
    "calc_spatial_smooth_5or9_point",
    "calc_forward_smooth",
    "calc_reverse_smooth",
    "calc_spatial_smooth_circular",
    "calc_spatial_smooth_window",
]


@generate_dataset_dispatcher
def calc_spatial_smooth_gaussian(
    data: xr.DataArray | xr.Dataset,
    n: int = 3,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
) -> xr.DataArray | xr.Dataset:
    """
    Filter with normal distribution of weights.

    Parameters
    ----------
    data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
        Some n-dimensional latitude-longitude dataset.
    n: :py:class:`int <int>`.
        Degree of filtering.

    Returns
    -------
    The filtered latitude-longitude dataset (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. seealso:

        :py:func:`metpy.calc.smooth_gaussian <metpy:metpy.calc.smooth_gaussian>`
    """
    data_original_dims_order = data.dims
    data = data.transpose(..., lat_dim, lon_dim)
    smoothed_data = mpcalc.smooth_gaussian(data, n)
    smoothed_data = smoothed_data.transpose(*data_original_dims_order)
    return dequantify_metpy_xarraydata(smoothed_data)


@generate_dataset_dispatcher
def calc_spatial_smooth_rectangular(
    data: xr.DataArray | xr.Dataset,
    rectangle_shapes: int | Sequence[int] = 3,
    times: int = 1,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
) -> xr.DataArray | xr.Dataset:
    """
    Filter with a rectangular window smoother.

    Parameters
    ----------
    data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
        Some n-dimensional latitude-longitude dataset.
    rectangle_shapes: :py:class:`int <int>` or :py:class:`Sequence[int]`, default `3`.
        Shape of rectangle along the trailing dimension(s) of the data.
    times: :py:class:`int <int>` or , default `1`.
        The number of times to apply the filter to the data.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

    Returns
    -------
    The filtered latitude-longitude dataset (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. seealso:

        :py:func:`metpy.calc.smooth_rectangular <metpy:metpy.calc.smooth_rectangular>`
    """
    data_original_dims_order = data.dims
    data = data.transpose(..., lat_dim, lon_dim)
    smoothed_data = mpcalc.smooth_rectangular(data, size=rectangle_shapes, passes=times)
    smoothed_data = smoothed_data.transpose(*data_original_dims_order)
    return dequantify_metpy_xarraydata(smoothed_data)


@generate_dataset_dispatcher
def calc_spatial_smooth_5or9_point(
    data: xr.DataArray | xr.Dataset,
    n: {5, 9} = 5,
    times: int = 1,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
) -> xr.DataArray | xr.Dataset:
    """
    Filter with an 5-point or 9-point smoother.

    Parameters
    ----------
    data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
        Some n-dimensional latitude-longitude dataset.
    n: {5, 9}, default `5`.
        The number of points to use in smoothing, only valid inputs are 5 and 9.
    times: :py:class:`int <int>` or , default `1`.
        The number of times to apply the filter to the data.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

    Returns
    -------
    The filtered latitude-longitude dataset (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. seealso:

        :py:func:`metpy.calc.smooth_n_point <metpy:metpy.calc.smooth_n_point>`
    """
    data_original_dims_order = data.dims
    data = data.transpose(..., lat_dim, lon_dim)
    smoothed_data = mpcalc.smooth_n_point(data, n=n, passes=times)
    smoothed_data = smoothed_data.transpose(*data_original_dims_order)
    return dequantify_metpy_xarraydata(smoothed_data)


def _spatial_smooth(data, n, S, times, normalize_weights, lon_dim, lat_dim):
    if n == 9:
        weights = np.array(
            [
                [S**2 / 4.0, S * (1 - S) / 2.0, S**2 / 4.0],
                [S * (1 - S) / 2.0, 1 - 2 * S * (1 - S), S * (1 - S) / 2.0],
                [S**2 / 4.0, S * (1 - S) / 2.0, S**2 / 4.0],
            ]
        )
    elif n == 5:
        weights = np.array(
            [[0.0, S / 4.0, 0.0], [S / 4.0, 1 - S, S / 4.0], [0.0, S / 4.0, 0.0]]
        )
    else:
        raise ValueError(f"`n` should be 5 or 9, not {n}")

    return calc_spatial_smooth_window(
        data=data,
        window=weights,
        times=times,
        normalize_weights=normalize_weights,
        lon_dim=lon_dim,
        lat_dim=lat_dim,
    )


@generate_dataset_dispatcher
def calc_forward_smooth(
    data: xr.DataArray | xr.Dataset,
    n: int = 5,
    S: float = 0.5,
    times: int = 1,
    normalize_weights=False,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
) -> xr.DataArray | xr.Dataset:
    """
    Filter with the forward smooth.

    Parameters
    ----------
    data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
        Some n-dimensional latitude-longitude dataset.
    n: {5, 9}, default `5`.
        The number of points to use in smoothing, only valid inputs are 5 and 9.
    S: :py:class:`float <float>`, default `0.5`.
        The smoothing factor.
    times: :py:class:`int <int>`, default `1`.
        The number of times to apply the filter to the data.
    normalize_weights: :py:class:`bool <bool>`, default `False`.
        If `True`, divide the values in window by the sum of all values in the window to obtain the normalized smoothing weights. If `False`, use supplied values directly as the weights.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

    Returns
    -------
    The filtered latitude-longitude dataset (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. note::

        For any discrete variable, its value on the :math:`x`-axis can be described as :math:`F_i = F(x_i)`, where :math:`x_i = i \\Delta x, i = 0, \\pm 1, \\pm 2`.

        Define a one-dimensional smoothing operator

        .. math::

            \\overset{\\sim}{F}_{i}^x = (1-S)F_i + \\frac{S}{2} (F_{i+1} + F_{i-1})

        Here, :math:`S` is the spatial smoothing coefficient. Since the above formula only involves the values of three grid data points (:math:`F_i, F_{i+1}, F_{i-1}`), it is also referred to as a three-point smoothing.

        Correspondingly, for a two-dimensional variable :math:`F_{i,j}`, smoothing is performed in the :math:`x` direction as follows

        .. math::

            \\overset{\\sim}{F}_{i,j}^x = F_{i,j} + \\frac{S}{2} (F_{i+1,j} + F_{i-1,j}-2F_{i,j})

        And, smoothing is done in the :math:`y` direction separately

        .. math::

            \\overset{\\sim}{F}_{i,j}^y = F_{i,j} + \\frac{S}{2} (F_{i,j+1} + F_{i,j-1}-2F_{i,j})

        Taking the average of both results:

        .. math::

            \\overset{\\sim}{F}_{i,j}^{x,y} = \\frac{\\overset{\\sim}{F}_{i,j}^x + \\overset{\\sim}{F}_{i,j}^y}{2} = F_{i,j} + \\frac{S}{4}(F_{i+1,j}+F_{i,j+1}+F_{i-1,j}+F_{i,j-1}-4F_{i,j})

        At this point, the following "window" matrix :math:`\\mathbf{W}` can be obtained

        .. math::

            \\mathbf{W} =
            \\left(
            \\begin{array}{ccc}
                F_{i-1,j+1} &  F_{i,j+1} & F_{i+1,j+1} \\\\
                F_{i-1,j} &  F_{i,j} & F_{i+1,j} \\\\
                F_{i-1,j-1} &  F_{i,j-1} & F_{i+1,j-1}
            \\end{array}
            \\right)
            =
            \\left(
            \\begin{array}{ccc}
                0 &  \\frac{S}{4} & 0 \\\\
                \\frac{S}{4} &  (1-S) & \\frac{S}{4} \\\\
                0 &  \\frac{S}{4} & 0
            \\end{array}
            \\right)

        When :math:`S=0.5`, we have

        .. math::

            \\mathbf{W} =
            \\left(
            \\begin{array}{ccc}
                0 &  0.125 & 0 \\\\
                0.125 &  0.5 & 0.125 \\\\
                0 &  0.125 & 0
            \\end{array}
            \\right)

        This is the "window" matrix used by the function when :math:`n=5`.

        For the variable :math:`F_{i,j}`, applying a three-point smoothing in the :math:`x`-axis direction followed by a three-point smoothing in the :math:`y`-axis direction results in the following nine-point smoothing scheme

        .. math::

            \\overset{\\sim}{F}_{i,j}^{x,y} = {F}_{i,j} + \\frac{S}{2}(1-S)(F_{i+1,j}+F_{i,j+1}+F_{i-1,j}+F_{i,j-1}-4F_{i,j})+\\frac{S^2}{4}(F_{i+1,j+1}+F_{i-1,j+1}+F_{i-1,j-1}+F_{i+1,j-1}-4F_{i,j})

        At this point, the corresponding "window" matrix :math:`\\mathbf{W}` is given by

        .. math::

            \\mathbf{W} =
            \\left(
            \\begin{array}{ccc}
                F_{i-1,j+1} &  F_{i,j+1} & F_{i+1,j+1} \\\\
                F_{i-1,j} &  F_{i,j} & F_{i+1,j} \\\\
                F_{i-1,j-1} &  F_{i,j-1} & F_{i+1,j-1}
            \\end{array}
            \\right)
            =
            \\left(
            \\begin{array}{ccc}
                \\frac{S^2}{4} &  \\frac{S}{2}(1-S) & \\frac{S^2}{4} \\\\
                \\frac{S}{2}(1-S) &  1-2S(1-S) & \\frac{S}{2}(1-S) \\\\
                \\frac{S^2}{4} &  \\frac{S}{2}(1-S) & \\frac{S^2}{4}
            \\end{array}
            \\right)

        When :math:`S=0.5`, we have

        .. math::

            \\mathbf{W} =
            \\left(
            \\begin{array}{ccc}
                0.0625 &  0.125 & 0.0625 \\\\
                0.125 &  0.5 & 0.125 \\\\
                0.0625 &  0.125 & 0.0625
            \\end{array}
            \\right)

        This is the "window" matrix used by the function when :math:`n=9`.

    .. seealso::
        :py:func:`metpy.calc.smooth_n_point <metpy:metpy.calc.smooth_n_point>`
    """
    return _spatial_smooth(
        data, n, S, times, normalize_weights, lon_dim=lon_dim, lat_dim=lat_dim
    )


def calc_reverse_smooth(
    data: xr.DataArray | xr.Dataset,
    n: int = 5,
    S: float = -0.5,
    times: int = 1,
    normalize_weights=False,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
) -> xr.DataArray | xr.Dataset:
    """
    Filter with the reverse smooth.

    Parameters
    ----------
    data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
        Some n-dimensional latitude-longitude dataset.
    n: {5, 9}, default `5`.
        The number of points to use in smoothing, only valid inputs are 5 and 9.
    S: :py:class:`float <float>`, default `-0.5`.
        The smoothing factor.
    times: :py:class:`int <int>`, default `1`.
        The number of times to apply the filter to the data.
    normalize_weights: :py:class:`bool <bool>`, default `False`.
        If `True`, divide the values in window by the sum of all values in the window to obtain the normalized smoothing weights. If `False`, use supplied values directly as the weights.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

    Returns
    -------
    The filtered latitude-longitude dataset (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. note::

        For any discrete variable, its value on the :math:`x`-axis can be described as :math:`F_i = F(x_i)`, where :math:`x_i = i \\Delta x, i = 0, \\pm 1, \\pm 2`.

        Define a one-dimensional smoothing operator

        .. math::

            \\overset{\\sim}{F}_{i}^x = (1-S)F_i + \\frac{S}{2} (F_{i+1} + F_{i-1})

        Here, :math:`S` is the spatial smoothing coefficient. Since the above formula only involves the values of three grid data points (:math:`F_i, F_{i+1}, F_{i-1}`), it is also referred to as a three-point smoothing.

        Correspondingly, for a two-dimensional variable :math:`F_{i,j}`, smoothing is performed in the :math:`x` direction as follows

        .. math::

            \\overset{\\sim}{F}_{i,j}^x = F_{i,j} + \\frac{S}{2} (F_{i+1,j} + F_{i-1,j}-2F_{i,j})

        And, smoothing is done in the :math:`y` direction separately

        .. math::

            \\overset{\\sim}{F}_{i,j}^y = F_{i,j} + \\frac{S}{2} (F_{i,j+1} + F_{i,j-1}-2F_{i,j})

        Taking the average of both results:

        .. math::

            \\overset{\\sim}{F}_{i,j}^{x,y} = \\frac{\\overset{\\sim}{F}_{i,j}^x + \\overset{\\sim}{F}_{i,j}^y}{2} = F_{i,j} + \\frac{S}{4}(F_{i+1,j}+F_{i,j+1}+F_{i-1,j}+F_{i,j-1}-4F_{i,j})

        At this point, the following "window" matrix :math:`\\mathbf{W}` can be obtained

        .. math::

            \\mathbf{W} =
            \\left(
            \\begin{array}{ccc}
                F_{i-1,j+1} &  F_{i,j+1} & F_{i+1,j+1} \\\\
                F_{i-1,j} &  F_{i,j} & F_{i+1,j} \\\\
                F_{i-1,j-1} &  F_{i,j-1} & F_{i+1,j-1}
            \\end{array}
            \\right)
            =
            \\left(
            \\begin{array}{ccc}
                0 &  \\frac{S}{4} & 0 \\\\
                \\frac{S}{4} &  (1-S) & \\frac{S}{4} \\\\
                0 &  \\frac{S}{4} & 0
            \\end{array}
            \\right)

        When :math:`S=0.5`, we have

        .. math::

            \\mathbf{W} =
            \\left(
            \\begin{array}{ccc}
                0 &  0.125 & 0 \\\\
                0.125 &  0.5 & 0.125 \\\\
                0 &  0.125 & 0
            \\end{array}
            \\right)

        This is the "window" matrix used by the function when :math:`n=5`.

        For the variable :math:`F_{i,j}`, applying a three-point smoothing in the :math:`x`-axis direction followed by a three-point smoothing in the :math:`y`-axis direction results in the following nine-point smoothing scheme

        .. math::

            \\overset{\\sim}{F}_{i,j}^{x,y} = {F}_{i,j} + \\frac{S}{2}(1-S)(F_{i+1,j}+F_{i,j+1}+F_{i-1,j}+F_{i,j-1}-4F_{i,j})+\\frac{S^2}{4}(F_{i+1,j+1}+F_{i-1,j+1}+F_{i-1,j-1}+F_{i+1,j-1}-4F_{i,j})

        At this point, the corresponding "window" matrix :math:`\\mathbf{W}` is given by

        .. math::

            \\mathbf{W} =
            \\left(
            \\begin{array}{ccc}
                F_{i-1,j+1} &  F_{i,j+1} & F_{i+1,j+1} \\\\
                F_{i-1,j} &  F_{i,j} & F_{i+1,j} \\\\
                F_{i-1,j-1} &  F_{i,j-1} & F_{i+1,j-1}
            \\end{array}
            \\right)
            =
            \\left(
            \\begin{array}{ccc}
                \\frac{S^2}{4} &  \\frac{S}{2}(1-S) & \\frac{S^2}{4} \\\\
                \\frac{S}{2}(1-S) &  1-2S(1-S) & \\frac{S}{2}(1-S) \\\\
                \\frac{S^2}{4} &  \\frac{S}{2}(1-S) & \\frac{S^2}{4}
            \\end{array}
            \\right)

        When :math:`S=0.5`, we have

        .. math::

            \\mathbf{W} =
            \\left(
            \\begin{array}{ccc}
                0.0625 &  0.125 & 0.0625 \\\\
                0.125 &  0.5 & 0.125 \\\\
                0.0625 &  0.125 & 0.0625
            \\end{array}
            \\right)

        This is the "window" matrix used by the function when :math:`n=9`.

    .. seealso::
        :py:func:`metpy.calc.smooth_n_point <metpy:metpy.calc.smooth_n_point>`
    """
    return _spatial_smooth(
        data, n, S, times, normalize_weights, lon_dim=lon_dim, lat_dim=lat_dim
    )


@generate_dataset_dispatcher
def calc_spatial_smooth_circular(
    data: xr.DataArray | xr.Dataset,
    radius: int,
    times: int = 1,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
) -> xr.DataArray | xr.Dataset:
    """
    Filter with a circular window smoother.

    Parameters
    ----------
    data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
        Some n-dimensional latitude-longitude dataset.
    radius: :py:class:`int <int>`.
        Radius of the circular smoothing window. The “diameter” of the circle (width of smoothing window) is :math:`2 \\cdot \\mathrm{radius} + 1` to provide a smoothing window with odd shape.
    times: :py:class:`int <int>`, default `1`.
        The number of times to apply the filter to the data.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

    Returns
    -------
    The filtered latitude-longitude dataset (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. seealso:

        :py:func:`metpy.calc.smooth_circular <metpy:metpy.calc.smooth_circular>`
    """
    data_original_dims_order = data.dims
    data = data.transpose(..., lat_dim, lon_dim)
    smoothed_data = mpcalc.smooth_circular(data, radius=radius, passes=times)
    smoothed_data = smoothed_data.transpose(*data_original_dims_order)
    return dequantify_metpy_xarraydata(smoothed_data)


@generate_dataset_dispatcher
def calc_spatial_smooth_window(
    data: xr.DataArray | xr.Dataset,
    window: np.ndarray | xr.DataArray,
    times: int = 1,
    normalize_weights: bool = False,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
) -> xr.DataArray | xr.Dataset:
    """
    Filter with an arbitrary window smoother.

    Parameters
    ----------
    data: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.
        Some n-dimensional latitude-longitude dataset.
    window : :py:class:`numpy.ndarray <numpy.ndarray>` or :py:class:`xarray.DataArray<xarray.DataArray>`.
        Window to use in smoothing. Can have dimension less than or equal to :math:`N`. If dimension less than :math:`N`, the scalar grid will be smoothed along its trailing dimensions. Shape along each dimension must be odd.
    times: :py:class:`int <int>`, default `1`.
        The number of times to apply the filter to the data.
    normalize_weights: :py:class:`bool <bool>`, default `False`.
        If `True`, divide the values in window by the sum of all values in the window to obtain the normalized smoothing weights. If `False`, use supplied values directly as the weights.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

    Returns
    -------
    The filtered latitude-longitude dataset (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. seealso:

        :py:func:`metpy.calc.smooth_window <metpy:metpy.calc.smooth_window>`
    """
    data_original_dims_order = data.dims
    data = data.transpose(..., lat_dim, lon_dim)
    smoothed_data = mpcalc.smooth_window(
        data, window=window, passes=times, normalize_weights=normalize_weights
    )
    smoothed_data = smoothed_data.transpose(*data_original_dims_order)
    return dequantify_metpy_xarraydata(smoothed_data)
