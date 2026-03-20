"""
Functions for curved quiver plots.
"""

from __future__ import annotations

import warnings
from collections.abc import Hashable
from typing import TYPE_CHECKING

import matplotlib.font_manager
from .modplot import CurvedQuiverplotSet
from typing import Literal
from matplotlib.patches import FancyArrowPatch

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

if TYPE_CHECKING:
    from matplotlib.axes import Axes

__all__ = ["curved_quiver", "add_curved_quiverkey"]


def curved_quiver(
    ds: xr.Dataset,
    x: Hashable,
    y: Hashable,
    u: Hashable,
    v: Hashable,
    ax: Axes | None = None,
    density=1,
    linewidth=None,
    color=None,
    cmap=None,
    norm=None,
    arrowsize=1,
    arrowstyle="-|>",
    transform=None,
    zorder=None,
    start_points=None,
    integration_direction="both",
    grains=15,
    broken_streamlines=True,
    regrid_shape: int | tuple[int, int] | None = None,
    ref_magnitude: float | None = None,
    ref_length: float | None = None,
    min_frac_length: float = 0.0,
    length_norm: Literal["reference", "max", "percentile"] = "reference",
    mask_density: int | tuple[int, int] = 10,
    line_start_stride: int = 1,
    arrow_stride: int = 1,
    min_distance: float = 0.0,
    arrow_head_ratio: float = 1.0,
    arrow_position: float = 1.0,
) -> CurvedQuiverplotSet:
    """
    Plot streamlines of a vector flow.

    Parameters
    ----------
    ds: :py:class:`xarray.Dataset`.
        Wind dataset.
    x: Hashable or None, optional.
        Variable name for x-axis.
    y: Hashable or None, optional.
        Variable name for y-axis.
    u: Hashable or None, optional.
        Variable name for the u velocity (in `x` direction).
    v: Hashable or None, optional.
        Variable name for the v velocity (in `y` direction).
    ax: :py:class:`matplotlib.axes.Axes`, optional.
        Axes on which to plot. By default, use the current axes. Mutually exclusive with `size` and `figsize`.
    density: :py:class:`float <float>` or (:py:class:`float <float>`, :py:class:`float <float>`)
        Controls the closeness of streamlines. When ``density = 1``, the domain
        is divided into a 30x30 grid. *density* linearly scales this grid.
        Each cell in the grid can have, at most, one traversing streamline.
        For different densities in each direction, use a tuple
        (density_x, density_y).
    linewidth: :py:class:`float <float>` or 2D array
        The width of the streamlines. With a 2D array the line width can be
        varied across the grid. The array must have the same shape as *u*
        and *v*.
    color: color or 2D array
        The streamline color. If given an array, its values are converted to
        colors using *cmap* and *norm*.  The array must have the same shape
        as *u* and *v*.
    cmap, norm
        Data normalization and colormapping parameters for *color*; only used
        if *color* is an array of floats. See `~.Axes.imshow` for a detailed
        description.
    arrowsize: :py:class:`float <float>`
        Scaling factor for the arrow size.
    arrowstyle: :py:class:`str <str>`
        Arrow style specification.
        See `~matplotlib.patches.FancyArrowPatch`.
    start_points: (N, 2) array
        Coordinates of starting points for the streamlines in data coordinates
        (the same coordinates as the *x* and *y* arrays).
    zorder: :py:class:`float <float>`
        The zorder of the streamlines and arrows.
        Artists with lower zorder values are drawn first.
    integration_direction: {'forward', 'backward', 'both'}, default: 'both'
        Integrate the streamline in forward, backward or both directions.
    broken_streamlines: :py:class:`bool <bool>`, default: True
        If False, forces streamlines to continue until they
        leave the plot domain.  If True, they may be terminated if they
        come too close to another streamline.
    regrid_shape: :py:class:`int <int>` or (:py:class:`int <int>`, :py:class:`int <int>`) or None, default: None
        Regrid the vector field onto a regular grid in the target map projection,
        similar to `cartopy`'s quiver and barbs methods. This is only supported
        when plotting on a `cartopy` GeoAxes and `transform` is provided.
    ref_magnitude: :py:class:`float <float>` or None, default: None
        Reference vector magnitude used to convert field magnitude into curved
        vector length. If None, use the rule defined by `length_norm`.
    ref_length: :py:class:`float <float>` or None, default: None
        Reference curved-vector length in the internal axes-normalized scale.
        If None, keep the previous density/grains-based behavior.
    min_frac_length: :py:class:`float <float>`, default: 0.0
        Minimum fraction of `ref_length` retained for weak vectors.
    length_norm: {"reference", "max", "percentile"}, default: ``"reference"``
        Fallback normalization rule used when `ref_magnitude` is not explicitly
        provided.
    mask_density: :py:class:`int <int>` or tuple, default: 10
        Density of the occupancy mask used to avoid overly dense curved vectors.
    line_start_stride: :py:class:`int <int>`, default: 1
        Seed-point subsampling stride.
    arrow_stride: :py:class:`int <int>`, default: 1
        Segment stride for placing multiple arrow heads along a curved trajectory.
    min_distance: :py:class:`float <float>`, default: 0.0
        Minimum spacing between generated seed points in data coordinates.
    arrow_head_ratio: :py:class:`float <float>`, default: 1.0
        Multiplier applied to arrow-head size independently of line length.
    arrow_position: :py:class:`float <float>`, default: 0.8
        Relative position of the first arrow along each trajectory. Larger values
        move the arrow closer to the front/head of the curved vector.

    Returns
    -------
    :py:class:`easyclimate.plot.modplot.CurvedQuiverplotSet <easyclimate.plot.modplot.CurvedQuiverplotSet>`
        Container object with attributes

        - ``lines``: `.LineCollection` of streamlines

        - ``arrows``: `.PatchCollection` containing `.FancyArrowPatch`
          objects representing the arrows half-way along streamlines.

            This container will probably change in the future to allow changes
            to the colormap, alpha, etc. for both lines and arrows, but these
            changes should be backward compatible.

    .. seealso::
        - https://github.com/matplotlib/matplotlib/issues/20038
        - https://github.com/kieranmrhunt/curved-quivers
        - https://github.com/Deltares/dfm_tools/issues/483
        - https://github.com/NCAR/geocat-viz/issues/4
        - https://docs.xarray.dev/en/stable/generated/xarray.Dataset.plot.streamplot.html#xarray.Dataset.plot.streamplot

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_curve_quiver.py
    """
    from .modplot import velovect

    ds = ds.sortby(y)
    x = ds[x].data
    y = ds[y].data
    u = ds[u].data
    v = ds[v].data

    # https://scitools.org.uk/cartopy/docs/latest/gallery/miscellanea/logo.html#sphx-glr-gallery-miscellanea-logo-py
    if ax == None:
        ax = plt.gca()

    if regrid_shape is not None:
        if transform is None:
            raise ValueError("`regrid_shape` requires `transform` to be specified.")
        if not hasattr(ax, "projection"):
            raise ValueError(
                "`regrid_shape` is only supported when `ax` is a cartopy GeoAxes."
            )

        from cartopy.vector_transform import vector_scalar_to_grid

        scalar_fields = []
        scalar_field_names = []
        if isinstance(color, np.ndarray):
            scalar_fields.append(color)
            scalar_field_names.append("color")
        if isinstance(linewidth, np.ndarray):
            scalar_fields.append(linewidth)
            scalar_field_names.append("linewidth")

        regridded = vector_scalar_to_grid(
            transform,
            ax.projection,
            regrid_shape,
            x,
            y,
            u,
            v,
            *scalar_fields,
            target_extent=ax.get_extent(crs=ax.projection),
        )
        x, y, u, v = regridded[:4]

        for field_name, field_value in zip(scalar_field_names, regridded[4:]):
            if field_name == "color":
                color = field_value
            elif field_name == "linewidth":
                linewidth = field_value

        transform = ax.transData

    if type(transform).__name__ == "PlateCarree":
        transform = transform._as_mpl_transform(ax)

    # https://github.com/Deltares/dfm_tools/issues/294
    # https://github.com/Deltares/dfm_tools/blob/main/dfm_tools/modplot.py
    obj = velovect(
        ax,
        x,
        y,
        u,
        v,
        density=density,
        linewidth=linewidth,
        color=color,
        cmap=cmap,
        norm=norm,
        arrowsize=arrowsize,
        arrowstyle=arrowstyle,
        transform=transform,
        zorder=zorder,
        start_points=start_points,
        integration_direction=integration_direction,
        grains=grains,
        broken_streamlines=broken_streamlines,
        ref_magnitude=ref_magnitude,
        ref_length=ref_length,
        min_frac_length=min_frac_length,
        length_norm=length_norm,
        mask_density=mask_density,
        line_start_stride=line_start_stride,
        arrow_stride=arrow_stride,
        min_distance=min_distance,
        arrow_head_ratio=arrow_head_ratio,
        glyph_mode=False,
        arrow_position=arrow_position,
    )
    return obj


def add_curved_quiverkey(
    curved_quiver: CurvedQuiverplotSet,
    pos: tuple[float, float],
    U: float,
    label: str,
    ax: Axes | None = None,
    color: str = "black",
    angle: float = 0,
    labelpos: Literal["N", "S", "E", "W"] = "N",
    labelsep: float = 0.02,
    labelcolor: str = None,
    fontproperties: dict = None,
    zorder=None,
    ref_point=None,
):
    """
    Add a curved-quiver key to a plot.

    The key is drawn in axes coordinates and represents a reference vector
    magnitude for the object returned by :py:class:`CurvedQuiverplotSet <easyclimate.plot.modplot.CurvedQuiverplotSet>`.
    Its length is scaled consistently with the internal normalization used by
    the curved-quiver algorithm. When plotting on projected axes, a local
    projected length can optionally be estimated using a reference point.

    Parameters
    ----------
    curved_quiver: :py:class:`easyclimate.plot.modplot.CurvedQuiverplotSet <easyclimate.plot.modplot.CurvedQuiverplotSet>`
        The curved-quiver object.
    ax: :py:class:`matplotlib.axes.Axes`
        Axes on which to draw the quiver key.
    pos: tuple[float, float]
        Position ``(x, y)`` of the quiver key in axes coordinates.
        For example, ``(0.5, 0.5)`` places the key at the center of the axes.
        Values outside the ``[0, 1]`` range are allowed so the key can be drawn
        outside the axes frame.
    U: :py:class:`float <float>`
        Vector magnitude represented by the quiver key.
    label: :py:class:`str <str>`
        Label displayed alongside the quiver key.
    color: matplotlib color, default: ``"black"``
        Color of both the key arrow and the label text unless
        [`labelcolor`](src/easyclimate/plot/curved_quiver_plot.py:168) is given.
    angle: :py:class:`float <float>`, default: 0
        Arrow angle in degrees, measured counterclockwise from the positive
        x-direction of the axes.
    labelpos: {'N', 'S', 'E', 'W'}, default: ``'N'``
        Position of the label relative to the arrow:

        - ``'N'``: above the arrow,
        - ``'S'``: below the arrow,
        - ``'E'``: at the arrow-head side,
        - ``'W'``: at the arrow-tail side.
    labelsep: :py:class:`float <float>`, default: 0.02
        Separation between the arrow and the label in axes coordinates.
    labelcolor: matplotlib color or None, default: None
        Color of the label text.
    fontproperties: :py:class:`dict <dict>` or :py:class:`matplotlib.font_manager.FontProperties` or None, default: None
        Font properties passed to :py:func:`matplotlib.axes.Axes.text` when drawing the label.
    zorder: :py:class:`float <float>` or None, default: None
        Drawing order of the quiver key. If None, a default value above the
        curved quiver is used.
    ref_point: tuple[float, float] or None, default: None
        Optional reference point ``(x, y)`` in data coordinates.
        When provided on projected axes, the key length is converted using the
        local projection scale near this point, so the displayed key length is
        more representative of the curved-quiver length at that location.
        For non-projected axes, this parameter is usually unnecessary.

    .. note::

        For projected axes, the displayed length of a vector can vary with
        location. Supplying ``ref_point`` makes the quiver key correspond
        to a local projected length near that position.

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_curve_quiver.py
    """
    if ax == None:
        ax = plt.gca()

    # Calculate arrow length in axes coordinates using the same reference-based
    # scaling used during streamline integration.
    resolution = curved_quiver.resolution
    ref_magnitude = getattr(curved_quiver, "ref_magnitude", None)
    if ref_magnitude is None or ref_magnitude == 0:
        ref_magnitude = curved_quiver.speed_max

    if ref_magnitude == 0:
        normalized_length = 0.0
    else:
        normalized_length = U / ref_magnitude
        if curved_quiver.integration_direction == "both":
            normalized_length /= 2.0

    min_frac_length = getattr(curved_quiver, "min_frac_length", 0.0)
    normalized_length = max(min_frac_length, normalized_length)

    arrow_length = resolution * normalized_length

    # Convert angle to radians.
    angle_rad = np.radians(angle)

    # `resolution` is defined in axes-like normalized coordinates inside
    # [`_get_integrator()`](src/easyclimate/plot/modplot.py:528), so the key
    # length should be drawn directly in axes-fraction coordinates.
    dx = arrow_length * np.cos(angle_rad)
    dy = arrow_length * np.sin(angle_rad)

    if ref_point is not None:
        ref_x, ref_y = ref_point
        ref_transform = curved_quiver.transform

        data_dx = arrow_length * curved_quiver.width * np.cos(angle_rad)
        data_dy = arrow_length * curved_quiver.height * np.sin(angle_rad)

        start_disp = ref_transform.transform((ref_x, ref_y))
        end_disp = ref_transform.transform((ref_x + data_dx, ref_y + data_dy))

        dx = (end_disp[0] - start_disp[0]) / ax.bbox.width
        dy = (end_disp[1] - start_disp[1]) / ax.bbox.height
    elif curved_quiver.transform != ax.transData:
        warnings.warn(
            "Projected curved quiver key without `ref_point` uses domain-"
            "normalized scaling; pass `ref_point=(x, y)` for locally projected "
            "key length.",
            UserWarning,
            stacklevel=2,
        )

    # Define arrow properties. Using [`FancyArrowPatch`](src/easyclimate/plot/curved_quiver_plot.py:14)
    # directly avoids annotation-specific offsets and keeps the arrow geometry
    # aligned with the requested angle.
    arrow = FancyArrowPatch(
        posA=pos,
        posB=(pos[0] + dx, pos[1] + dy),
        arrowstyle="-|>",
        mutation_scale=10,
        linewidth=1.2,
        facecolor=color,
        edgecolor=color,
        transform=ax.transAxes,
        clip_on=False,
        zorder=zorder if zorder is not None else 2,
    )
    arrow.set_clip_on(False)
    arrow.set_clip_path(None)
    ax.add_patch(arrow)

    # Tangential and normal directions for label placement.
    mid_x = pos[0] + dx / 2
    mid_y = pos[1] + dy / 2
    tan_x = np.cos(angle_rad)
    tan_y = np.sin(angle_rad)
    normal_x = -np.sin(angle_rad)
    normal_y = np.cos(angle_rad)

    # Calculate the label position based on `labelpos`
    if labelpos == "N":
        label_x = mid_x + labelsep * normal_x
        label_y = mid_y + labelsep * normal_y
        horizontalalignment = "center"
        verticalalignment = "bottom"
    elif labelpos == "S":
        label_x = mid_x - labelsep * normal_x
        label_y = mid_y - labelsep * normal_y
        horizontalalignment = "center"
        verticalalignment = "top"
    elif labelpos == "E":
        label_x = pos[0] + dx + labelsep * tan_x
        label_y = pos[1] + dy + labelsep * tan_y
        horizontalalignment = "left"
        verticalalignment = "center"
    elif labelpos == "W":
        label_x = pos[0] - labelsep * tan_x
        label_y = pos[1] - labelsep * tan_y
        horizontalalignment = "right"
        verticalalignment = "center"
    else:
        raise ValueError("labelpos must be one of 'N', 'S', 'E', or 'W'")

    # Add the label text
    ax.text(
        label_x,
        label_y,
        label,
        color=labelcolor if labelcolor is not None else color,
        va=verticalalignment,
        ha=horizontalalignment,
        transform=ax.transAxes,
        clip_on=False,
        fontproperties=fontproperties,
        zorder=zorder if zorder is not None else 2,
    )
