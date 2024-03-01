"""
Graph processing related functions
"""

from __future__ import annotations

import cartopy
import xarray as xr
import warnings
import matplotlib
import matplotlib.pyplot as plt
from ..core.utility import transfer_xarray_lon_from180TO360, assert_compared_version

# cartopy version check
check_return = assert_compared_version(cartopy.__version__, "0.20")
if check_return == 1:
    pass
else:
    print(
        "Cartopy version is not greater than 0.20, please update cartopy package. You can use Conda to update: `conda install -c conda-forge cartopy`"
    )

import cartopy.crs as ccrs
import matplotlib.ticker as ticker
import numpy as np
from geocat.viz import util as gvutil


def draw_Circlemap_PolarStereo(
    *,
    lat_range: tuple | list,
    add_gridlines: bool = True,
    lon_step: float = None,
    lat_step: float = None,
    ax: matplotlib.axes.Axes = None,
    draw_labels: bool = True,
    set_map_boundary_kwargs: dict = {},
    gridlines_kwargs: dict = {},
):
    """
    Utility function to set the boundary of ax to a path that surrounds a
    given region specified by latitude and longitude coordinates. This boundary
    is drawn in the projection coordinates and therefore follows any curves
    created by the projection. As of now, this works consistently for the
    North/South Polar Stereographic Projections.

    Parameters
    ----------
    lat_range : :py:class:`tuple`, :py:class:`list`.
        The two-tuple containing the start and end of the desired range of
        latitudes. The first entry must be smaller than the second entry.
        Both entries must be between [-90 , 90].
    add_gridlines: :py:class:`bool`.
        whether or not add gridlines and tick labels to a map.
    lon_step: :py:class:`float`.
        The step of grid lines in longitude.
    lat_step: :py:class:`float`.
        The step of grid lines in latitude.
    ax : :py:class:`matplotlib.axes.Axes`
        The axes to which the boundary will be applied.
    draw_labels: :py:class:`bool`.
        Whether to draw labels. Defaults to `True`.
    **set_map_boundary_kwargs: :py:class:`dict`.
        Additional keyword arguments to wrapped :py:func:`geocat.viz.util.set_map_boundary <geocat.viz:geocat.viz.util.set_map_boundary>`.
    **gridlines_kwargs: :py:class:`dict`.
        Additional keyword arguments to wrapped :py:class:`cartopy.mpl.gridliner.Gridliner <cartopy:cartopy.mpl.gridliner.Gridliner>`.
    .. seealso
        :py:func:`geocat.viz.util.set_map_boundary <geocat.viz:geocat.viz.util.set_map_boundary>`, :py:class:`cartopy.mpl.gridliner.Gridliner <cartopy:cartopy.mpl.gridliner.Gridliner>`.
    """
    # Get Axes
    if ax == None:
        ax = plt.gca()
    else:
        pass

    # Check the projection parameters
    if (type(ax.projection).__name__ == "NorthPolarStereo") or (
        type(ax.projection).__name__ == "SouthPolarStereo"
    ):
        pass
    else:
        raise TypeError(
            "The projection type of the Axes should be `cartopy.crs.NorthPolarStereo` or `cartopy.crs.SouthPolarStereo`, consider to specify the parameter `projection` as `cartopy.crs.NorthPolarStereo` or `cartopy.crs.SouthPolarStereo`. E.g. `fig, ax = plt.subplots(subplot_kw = {'projection': ccrs.NorthPolarStereo()})`."
        )

    gvutil.set_map_boundary(ax, [-180, 180], lat_range, **set_map_boundary_kwargs)

    if add_gridlines == True:
        if lon_step == None or lat_step == None:
            raise ValueError(
                "If `add_gridlines = True`, the parameters `lon_step`, `lat_step`, and `draw_labels` should be specified."
            )

        gl = ax.gridlines(
            crs=ccrs.PlateCarree(),
            draw_labels=draw_labels,
            xlocs=ticker.FixedLocator(np.arange(-180, 180, lon_step)),
            ylocs=ticker.FixedLocator(
                np.arange(min(lat_range), max(lat_range), lat_step)
            ),
            **gridlines_kwargs,
        )
    elif add_gridlines == False:
        pass
    else:
        raise ValueError(
            "`add_gridlines` is bool type, it should be `True` or `False`."
        )


def add_lon_cyclic(data_input: xr.DataArray, inter: float, lon_dim: str = "lon"):
    """
    Add a cyclic point to an array and optionally a corresponding coordinate.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
    inter: :py:class:`float<float>`
        Longitude interval (assuming longitude is arranged in a sequence of equal differences).
    lon_dim: :py:class:`str<str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

    .. seealso
        :py:func:`xarray.DataArray.pad <xarray:xarray.DataArray.pad>`, :py:func:`cartopy.util.add_cyclic_point <cartopy:cartopy.util.add_cyclic_point>`
    """
    lon_array_data_input = data_input[lon_dim].data

    if (lon_array_data_input < 0).any():
        warnings.warn(
            "It seems that the input data longitude range is from -180° to 180°. Currently automatically converted to it from 0° to 360°."
        )
        data_input = transfer_xarray_lon_from180TO360(data_input)

    temp = data_input.pad(pad_width={lon_dim: (0, 1)}, mode="wrap")
    result_data = temp.assign_coords({lon_dim: np.arange(0, 360 + inter, inter)})
    return result_data
