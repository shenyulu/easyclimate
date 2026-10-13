""" """

from __future__ import annotations

from typing import Any
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np


__all__ = ["plot_snapshot", "plot_cross_sections"]


def _metadata(result: Any, name: str, default: Any = None) -> Any:
    if isinstance(result, dict):
        return result.get(name, default)
    return result.attrs.get(name, default)


def _merge_kwargs(default: dict, user: dict | None) -> dict:
    """Merge user kwargs into defaults recursively."""
    result = deepcopy(default)

    if user is None:
        return result

    for key, value in user.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _merge_kwargs(result[key], value)
        else:
            result[key] = value

    return result


def plot_snapshot(
    result: Any,
    *,
    ax: plt.Axes | None = None,
    nq: int = 5,
    title: str | None = None,
    contourf_kwargs: dict = None,
    quiver_kwargs: dict = None,
) -> plt.Axes:
    """Draw the final horizontal height and wind response of a Gill experiment.

    ``eta`` is rendered as filled contours, so its sign and horizontal
    structure can be compared with the forcing stored in ``result.q``. The
    collocated ``u`` and interpolated ``v`` fields are thinned and drawn as a
    horizontal vector field. The function does not select a time index: it
    plots the final fields returned by :func:`simulate`.

    Parameters
    ----------
    result : xarray.Dataset
        Dataset returned by :func:`easyclimate.field.simulation.gill.simulate`.
        It must contain ``eta``, ``u``, ``v``, ``x``, and ``y``.
    ax : matplotlib.axes.Axes, optional
        Axes on which to draw. A new figure and axes are created by default.
    nq : int, default: 5
        Every ``nq``-th ``x`` and ``y`` point is retained before plotting the
        velocity vectors. Smaller values show a denser field; values greater
        than one reduce overlap on fine grids.
    title : str, optional
        Plot title. By default, the case label and final simulation time are
        taken from ``result.attrs``.
    contourf_kwargs : dict, optional
        Keyword arguments merged with the contour defaults for ``eta``. The
        defaults use 21 levels, set the upper colour limit to 3, and create a
        horizontal colourbar below the axes labelled ``eta``. Supply values
        such as ``levels``, ``vmin``, ``vmax``, ``cmap``, or ``cbar_kwargs``
        to change them.
    quiver_kwargs : dict, optional
        Keyword arguments merged with the vector defaults for ``u`` and ``v``.
        Vectors use data-coordinate angles, tail pivots, and a default scale
        of 0.5. Supply Matplotlib quiver options such as ``scale``, ``color``,
        or ``width`` to change their appearance.

    Returns
    -------
    matplotlib.axes.Axes
        The axes containing the snapshot.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        gill_asymmetric_heating
        gill_delta_y1
        gill_forced_stationary
        gill_gaussian_initial
        gill_kelvin_wave
        gill_symmetric_asymmetric
        gill_symmetric_gaussian_x
        gill_symmetric_heating
    """

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 4.5))

    default_contourf_kwargs = {
        "levels": 21,
        "vmax": 3,
        "cbar_kwargs": {
            "location": "bottom",
            "pad": 0.2,
            "aspect": 50,
            "shrink": 0.8,
            "label": "$\\eta$",
        },
    }

    default_quiver_kwargs = {
        "scale": 0.5,
        "angles": "xy",
        "pivot": "tail",
        "width": 0.004,
        "headwidth": 3.5,
        "headlength": 4.5,
        "headaxislength": 4.0,
        "minlength": 0.0,
    }

    contourf_kwargs = _merge_kwargs(
        default_contourf_kwargs,
        contourf_kwargs,
    )

    quiver_kwargs = _merge_kwargs(
        default_quiver_kwargs,
        quiver_kwargs,
    )

    result["eta"].plot.contourf(ax=ax, **contourf_kwargs)

    result.thin(x=nq, y=nq).plot.quiver(
        x="x", y="y", u="u", v="v", ax=ax, **quiver_kwargs
    )

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    if title is None:
        title = (
            f"{_metadata(result, 'case_label', '')}\nt = {_metadata(result, 't'):.2f}"
        )
    ax.set_title(title)
    return ax


def plot_cross_sections(
    result: Any,
    *,
    fig: plt.Figure | None = None,
    contourf_uw_kwargs: dict = None,
    contourf_vw_kwargs: dict = None,
    contourf_u_kwargs: dict = None,
) -> plt.Figure:
    """Draw three vertical diagnostics reconstructed from the horizontal solution.

    The shallow-water integration is two-dimensional; these panels do not
    introduce a separately integrated vertical velocity. Instead, they apply
    the first vertical sine or cosine structure over ``0 <= z <= lz`` to
    horizontal means from the final state. The three plotted quantities are:

    - ``psi_uw(z, x) = lz / pi * mean_y(u) * sin(pi * z / lz)`` for the
      zonal-vertical circulation diagnostic;
    - ``psi_vw(z, y) = -lz / pi * mean_x(v) * sin(pi * z / lz)`` for the
      meridional-vertical circulation diagnostic; and
    - ``u_mean(z, y) = mean_x(u) * cos(pi * z / lz)`` for the zonal-mean
      zonal-wind diagnostic.

    Each quantity is plotted with ``z`` increasing upward, and the vertical
    axes label the lower and upper boundaries as ``0`` and ``z``.

    Parameters
    ----------
    result : xarray.Dataset
        Dataset returned by :func:`easyclimate.field.simulation.gill.simulate`.
        It must contain ``u``, ``v``, and ``z``, and have an ``lz`` attribute.
    fig : matplotlib.figure.Figure, optional
        Empty figure on which to create the three stacked panels. If omitted,
        a new ``8 x 8`` inch figure is created.
    contourf_uw_kwargs, contourf_vw_kwargs, contourf_u_kwargs : dict, optional
        Keyword arguments merged with the default ``levels=21`` and passed to
        the filled-contour calls for the ``u-w``, ``v-w``, and ``u`` panels,
        respectively. Each mapping is independent, allowing panel-specific
        limits, levels, colormaps, or colourbar options.

    Returns
    -------
    matplotlib.figure.Figure
        The figure containing the three cross-section panels.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        gill_asymmetric_heating
        gill_delta_y1
        gill_forced_stationary
        gill_gaussian_initial
        gill_kelvin_wave
        gill_symmetric_asymmetric
        gill_symmetric_gaussian_x
        gill_symmetric_heating
    """

    lz_value = result.attrs["lz"]
    u_ymean = result["u"].mean(dim="y")
    v_xmean = result["v"].mean(dim="x")
    u_xmean = result["u"].mean(dim="x")
    z_value = result["z"]

    psiuw_value = lz_value / np.pi * u_ymean * np.sin(np.pi / lz_value * z_value)
    psivw_value = -lz_value / np.pi * v_xmean * np.sin(np.pi / lz_value * z_value)
    uxm_value = u_xmean * np.cos(np.pi / lz_value * z_value)

    psiuw_value = psiuw_value.transpose("z", ...)
    psivw_value = psivw_value.transpose("z", ...)
    uxm_value = uxm_value.transpose("z", ...)

    default_contourf_kwargs = {
        "levels": 21,
    }

    contourf_uw_kwargs = _merge_kwargs(
        default_contourf_kwargs,
        contourf_uw_kwargs,
    )
    contourf_vw_kwargs = _merge_kwargs(
        default_contourf_kwargs,
        contourf_vw_kwargs,
    )
    contourf_u_kwargs = _merge_kwargs(
        default_contourf_kwargs,
        contourf_u_kwargs,
    )

    if fig is None:
        fig, axes = plt.subplots(3, 1, figsize=(8, 8))
    else:
        axes = fig.subplots(3, 1)

    psiuw_value.plot.contourf(ax=axes[0], **contourf_uw_kwargs)
    axes[0].set_title(
        f"{_metadata(result, 'case_label', '')}\nt = {_metadata(result, 't'):.2f}\nu-w"
    )
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("z")

    psivw_value.plot.contourf(ax=axes[1], **contourf_vw_kwargs)
    axes[1].set_title("v-w")
    axes[1].set_xlabel("y")
    axes[1].set_ylabel("z")

    uxm_value.plot.contourf(ax=axes[2], **contourf_u_kwargs)
    axes[2].set_title("u")
    axes[2].set_xlabel("y")
    axes[2].set_ylabel("z")

    for ax in axes:
        ax.set_yticks([0, z_value[-1]], ["0", "z"])
    fig.tight_layout()
    return fig
