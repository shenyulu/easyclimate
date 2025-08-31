"""
Boreal Summer Intraseasonal Oscillation (BSISO)

This module provides functions to calculate and visualize the BSISO, a dominant mode of intraseasonal variability
in the Asian summer monsoon region during May-October. The BSISO is characterized by northward and eastward
propagation of convective anomalies, influencing monsoon rainfall and tropical cyclone activity. The methodology
is based on Lee et al. (2013), which defines multivariate indices using outgoing longwave radiation (OLR) and
850-hPa winds to capture BSISO phases.

.. seealso::
    - Lee, JY. (이준이), Wang, B., Wheeler, M.C. et al. Real-time multivariate indices for the boreal summer intraseasonal oscillation over the Asian summer monsoon region. Clim Dyn 40, 493–509 (2013). https://doi.org/10.1007/s00382-012-1544-4
    - https://iprc.soest.hawaii.edu/users/jylee/bsiso/
"""

__all__ = [
    "draw_bsiso1_phase_space_basemap",
    "draw_bsiso2_phase_space_basemap",
    "draw_bsiso_phase_space",
    "draw_bsiso2_phase_space",
    "calc_bsiso1_phase",
    "calc_bsiso2_phase",
    "get_times_for_phase",
    "calc_bsiso_analysis",
]

import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
import matplotlib.lines as lines
from ...core.datanode import DataNode
from ...core.variability import (
    remove_smooth_daily_annual_cycle_mean,
    remove_low_frequency_signal,
)
from ...core.extract import get_specific_months_data
from ...core.eof import get_EOF_model, calc_EOF_analysis
from ...core.stat import (
    calc_ttestSpatialPattern_spatial,
    calc_corr_spatial,
    calc_lead_lag_correlation_coefficients,
)
from ...core.utility import sort_ascending_latlon_coordinates
import logging
from rich.progress import Progress, BarColumn, TimeRemainingColumn, TextColumn

# Configure logging for progress updates
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def draw_bsiso1_phase_space_basemap(
    ax=None, add_phase_text: bool = True, add_location_text: bool = True
):
    """
    Create a phase space basemap for BSISO1 visualization.

    This function generates a standardized background diagram for BSISO1 analysis,
    plotting phase division lines (1-8), amplitude circles, and geographical annotations.
    The diagram uses normalized principal components (PC1 and PC2) as axes, representing
    the BSISO1 mode in a Cartesian coordinate system (-4 to 4 range). The BSISO1 mode
    captures the northward and eastward propagation of convective anomalies over the
    Asian monsoon region.

    Parameters
    ----------
    ax : :py:class:`matplotlib.axes.Axes`, optional
        Target axes for plotting. If None, uses the current axes (plt.gca()).
    add_phase_text : :py:class:`bool`, default True
        Whether to annotate BSISO phase numbers (1-8) on the diagram.
    add_location_text : :py:class:`bool`, default True
        Whether to annotate geographical regions corresponding to BSISO phases, such as:
        - Bay of Bengal & South China Sea
        - Indian Ocean & East Asia
        - Equatorial Indian Ocean
        - Western North Pacific
        - India & Maritime Continent

    Returns
    -------
    :py:class:`matplotlib.axes.Axes`
        The configured axes object with the BSISO1 phase space basemap.

    Notes
    -----
    - The diagram includes 8 radial lines at 45° intervals to divide phases, a unit circle
      for amplitude threshold (amplitude=1), and labeled axes (PC2, PC1).

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(figsize=(8, 8))
    >>> draw_bsiso1_phase_space_basemap(ax=ax)
    >>> plt.show()
    """
    if ax is None:
        ax = plt.gca()

    ax.set_aspect("equal")

    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)
    ax.set_xticks(range(-4, 5))
    ax.set_yticks(range(-4, 5))

    # plot mjo phase diagram lines
    line1 = lines.Line2D(
        [np.cos(np.pi / 4), 4],
        [np.sin(np.pi / 4), 4],
        color="black",
        linestyle="--",
        lw=1,
    )
    ax.add_line(line1)

    line2 = lines.Line2D(
        [np.cos(3 * np.pi / 4), -4],
        [np.sin(np.pi / 4), 4],
        color="black",
        linestyle="--",
        lw=1,
    )
    ax.add_line(line2)

    line3 = lines.Line2D(
        [np.cos(np.pi / 4), 4],
        [np.sin(7 * np.pi / 4), -4],
        color="black",
        linestyle="--",
        lw=1,
    )
    ax.add_line(line3)

    line4 = lines.Line2D(
        [np.cos(3 * np.pi / 4), -4],
        [np.sin(7 * np.pi / 4), -4],
        color="black",
        linestyle="--",
        lw=1,
    )
    ax.add_line(line4)

    line5 = lines.Line2D([-4, -1], [0, 0], color="black", linestyle="--", lw=1)
    ax.add_line(line5)

    line6 = lines.Line2D([1, 4], [0, 0], color="black", linestyle="--", lw=1)
    ax.add_line(line6)

    line7 = lines.Line2D([0, 0], [1, 4], color="black", linestyle="--", lw=1)
    ax.add_line(line7)

    line8 = lines.Line2D([0, 0], [-1, -4], color="black", linestyle="--", lw=1)
    ax.add_line(line8)

    amp1_circ = plt.Circle((0, 0), 1, color="black", fill=False)
    ax.add_patch(amp1_circ)

    if add_phase_text:
        # add phase diagram texts
        ax.text(1.5, 3, "6", size="x-large", weight="semibold")
        ax.text(-1.5, 3, "7", size="x-large", weight="semibold")
        ax.text(2.8, 1.5, "5", size="x-large", weight="semibold", ha="center")
        ax.text(-2.8, 1.5, "8", size="x-large", weight="semibold", ha="center")

        ax.text(1.5, -3, "3", size="x-large", weight="semibold")
        ax.text(-1.5, -3, "2", size="x-large", weight="semibold")
        ax.text(2.8, -1.5, "4", size="x-large", weight="semibold", ha="center")
        ax.text(-2.8, -1.5, "1", size="x-large", weight="semibold", ha="center")

    if add_location_text:
        ax.text(
            0,
            3.6,
            "Bay of Bengal & South China Sea",
            ha="center",
            bbox=dict(facecolor="white", edgecolor="white"),
        )
        ax.text(
            0,
            -3.8,
            "Indian Ocean\n& East Asia",
            ha="center",
            bbox=dict(facecolor="white", edgecolor="white"),
        )
        ax.text(
            -3.8,
            -1.5,
            "Eq Indian Ocean",
            va="center",
            rotation="vertical",
            bbox=dict(facecolor="white", edgecolor="white"),
        )
        ax.text(
            -3.8,
            1.5,
            "Western North Pacific",
            va="center",
            rotation="vertical",
            bbox=dict(facecolor="white", edgecolor="white"),
        )
        ax.text(
            3.5,
            0,
            "India & Maritime Continent",
            va="center",
            rotation=270,
            bbox=dict(facecolor="white", edgecolor="white"),
        )

    ax.set_xlabel("Normalized PC2", size="large")
    ax.set_ylabel("Normalized PC1", size="large")
    return ax


def draw_bsiso2_phase_space_basemap(
    ax=None, add_phase_text: bool = True, add_location_text: bool = True
):
    """
    Create a phase space basemap for BSISO2 visualization.

    This function generates a standardized background diagram for BSISO2 analysis,
    plotting phase division lines (1-8), amplitude circles, and geographical annotations.
    The diagram uses normalized principal components (PC3 and PC4) as axes, representing
    the BSISO2 mode, which captures higher-frequency intraseasonal variability compared
    to BSISO1.

    Parameters
    ----------
    ax : :py:class:`matplotlib.axes.Axes`, optional
        Target axes for plotting. If None, uses the current axes (plt.gca()).
    add_phase_text : :py:class:`bool`, default ``True``
        Whether to annotate BSISO phase numbers (1-8) on the diagram.
    add_location_text : :py:class:`bool`, default ``True``
        Whether to annotate geographical regions corresponding to BSISO phases, such as:
        - North East Asia
        - South East Asia
        - Philippine Sea
        - India & South China Sea
        - Indian Ocean
        - Western North Pacific
        - Bay of Bengal

    Returns
    -------
    :py:class:`matplotlib.axes.Axes`
        The configured axes object with the BSISO2 phase space basemap.

    Notes
    -----
    - The diagram includes 8 radial lines at 45° intervals, a unit circle for amplitude
      threshold (amplitude=1), and labeled axes (PC4, PC3).
    - BSISO2 typically represents shorter-period oscillations (10-20 days) compared to
      BSISO1 (30-60 days).

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(figsize=(8, 8))
    >>> draw_bsiso2_phase_space_basemap(ax=ax)
    >>> plt.show()
    """
    if ax is None:
        ax = plt.gca()

    ax.set_aspect("equal")

    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)
    ax.set_xticks(range(-4, 5))
    ax.set_yticks(range(-4, 5))

    # plot mjo phase diagram lines
    line1 = lines.Line2D(
        [np.cos(np.pi / 4), 4],
        [np.sin(np.pi / 4), 4],
        color="black",
        linestyle="--",
        lw=1,
    )
    ax.add_line(line1)

    line2 = lines.Line2D(
        [np.cos(3 * np.pi / 4), -4],
        [np.sin(np.pi / 4), 4],
        color="black",
        linestyle="--",
        lw=1,
    )
    ax.add_line(line2)

    line3 = lines.Line2D(
        [np.cos(np.pi / 4), 4],
        [np.sin(7 * np.pi / 4), -4],
        color="black",
        linestyle="--",
        lw=1,
    )
    ax.add_line(line3)

    line4 = lines.Line2D(
        [np.cos(3 * np.pi / 4), -4],
        [np.sin(7 * np.pi / 4), -4],
        color="black",
        linestyle="--",
        lw=1,
    )
    ax.add_line(line4)

    line5 = lines.Line2D([-4, -1], [0, 0], color="black", linestyle="--", lw=1)
    ax.add_line(line5)

    line6 = lines.Line2D([1, 4], [0, 0], color="black", linestyle="--", lw=1)
    ax.add_line(line6)

    line7 = lines.Line2D([0, 0], [1, 4], color="black", linestyle="--", lw=1)
    ax.add_line(line7)

    line8 = lines.Line2D([0, 0], [-1, -4], color="black", linestyle="--", lw=1)
    ax.add_line(line8)

    amp1_circ = plt.Circle((0, 0), 1, color="black", fill=False)
    ax.add_patch(amp1_circ)

    if add_phase_text:
        # add phase diagram texts
        ax.text(1.5, 3, "6", size="x-large", weight="semibold")
        ax.text(-1.5, 3, "7", size="x-large", weight="semibold")
        ax.text(2.8, 1.5, "5", size="x-large", weight="semibold", ha="center")
        ax.text(-2.8, 1.5, "8", size="x-large", weight="semibold", ha="center")

        ax.text(1.5, -3, "3", size="x-large", weight="semibold")
        ax.text(-1.5, -3, "2", size="x-large", weight="semibold")
        ax.text(2.8, -1.5, "4", size="x-large", weight="semibold", ha="center")
        ax.text(-2.8, -1.5, "1", size="x-large", weight="semibold", ha="center")

    if add_location_text:
        ax.text(
            -1.5,
            3.6,
            "N. East Asia",
            ha="center",
            bbox=dict(facecolor="white", edgecolor="white"),
        )
        ax.text(
            1.5,
            3.6,
            "S. East Asia",
            ha="center",
            bbox=dict(facecolor="white", edgecolor="white"),
        )
        ax.text(
            -1.5,
            -3.8,
            "Philippine Sea",
            ha="center",
            bbox=dict(facecolor="white", edgecolor="white"),
        )
        ax.text(
            1.7,
            -3.8,
            "India & S. China Sea",
            ha="center",
            bbox=dict(facecolor="white", edgecolor="white"),
        )
        ax.text(
            -3.8,
            -1.6,
            "India Ocean",
            va="center",
            rotation="vertical",
            bbox=dict(facecolor="white", edgecolor="white"),
        )
        ax.text(
            -3.8,
            1.6,
            "Western North Pacific",
            va="center",
            rotation="vertical",
            bbox=dict(facecolor="white", edgecolor="white"),
        )
        ax.text(
            3.5,
            0,
            "Bay of Bengal",
            va="center",
            rotation=270,
            bbox=dict(facecolor="white", edgecolor="white"),
        )

    ax.set_xlabel("Normalized PC4", size="large")
    ax.set_ylabel("Normalized PC3", size="large")
    return ax


def draw_bsiso_phase_space(
    bsiso_data: xr.Dataset,
    y_dim: str = "PC1",
    x_dim: str = "PC2",
    time_dim: str = "time",
    ax=None,
    color="blue",
    start_text="START",
    add_start_text: bool = True,
):
    """
    Visualize BSISO1 phase space trajectory using PC1 and PC2.

    Plots the temporal evolution of BSISO1 phases as a parametric curve in PC1-PC2 space,
    with an optional marker for the starting point. Combines scatter points and connecting
    lines to show the progression of BSISO phases, typically representing 30-60 day
    oscillations.

    Parameters
    ----------
    bsiso_data : :py:class:`xarray.Dataset`
        Dataset containing normalized principal components (PC1, PC2) and a time coordinate.
    y_dim : :py:class:`str`, default "PC1"
        Variable name for the y-axis component (PC1).
    x_dim : :py:class:`str`, default "PC2"
        Variable name for the x-axis component (PC2).
    time_dim : :py:class:`str`, default "time"
        Name of the time coordinate dimension.
    ax : :py:class:`matplotlib.axes.Axes`, optional
        Target axes for plotting. If None, uses the current axes (plt.gca()).
    color : :py:class:`str` or :py:class:`tuple`, default "blue"
        Color for the trajectory and points.
    start_text : :py:class:`str`, default "START"
        Text annotation for the trajectory starting point.
    add_start_text : :py:class:`bool`, default True
        Whether to display the start point annotation.

    Returns
    -------
    :py:class:`matplotlib.axes.Axes`
        Configured axes object with the BSISO1 phase space trajectory.

    Notes
    -----
    - Recommended to use with `draw_bsiso1_phase_space_basemap` for the background diagram.
    - Daily data is recommended for smooth trajectory visualization.

    Examples
    --------
    >>> import xarray as xr
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(figsize=(8, 8))
    >>> draw_bsiso1_phase_space_basemap(ax=ax)
    >>> draw_bsiso_phase_space(bsiso_data=ds, ax=ax, color="red")
    >>> plt.title("BSISO1 Phase Space Trajectory")
    >>> plt.show()
    """
    if ax is None:
        ax = plt.gca()

    # Plot trajectory lines
    ax.plot(bsiso_data[x_dim], bsiso_data[y_dim], color=color)
    # Plot data points
    ax.scatter(bsiso_data[x_dim], bsiso_data[y_dim], color=color, s=15)
    # Plot start point
    ax.scatter(
        bsiso_data[x_dim].isel({time_dim: 0}),
        bsiso_data[y_dim].isel({time_dim: 0}),
        s=15 * 5,
        c=color,
        marker="*",
    )
    # Add start point text
    if add_start_text:
        ax.text(
            bsiso_data[x_dim].isel({time_dim: 0}) + 0.1,
            bsiso_data[y_dim].isel({time_dim: 0}),
            start_text,
            ha="left",
            zorder=2,
        )
    return ax


def draw_bsiso2_phase_space(
    bsiso_data: xr.Dataset,
    pc3_dim: str = "PC3",
    pc4_dim: str = "PC4",
    time_dim: str = "time",
    ax=None,
    color="blue",
    start_text="START",
    add_start_text: bool = True,
):
    """
    Visualize BSISO2 phase space trajectory using PC3 and PC4.

    Plots the temporal evolution of BSISO2 phases as a parametric curve in PC3-PC4 space,
    with an optional marker for the starting point. Combines scatter points and connecting
    lines to show the progression of BSISO phases, typically representing 10-20 day
    oscillations.

    Parameters
    ----------
    bsiso_data : :py:class:`xarray.Dataset`
        Dataset containing normalized principal components (PC3, PC4) and a time coordinate.
    pc3_dim : :py:class:`str`, default "PC3"
        Variable name for the y-axis component (PC3).
    pc4_dim : :py:class:`str`, default "PC4"
        Variable name for the x-axis component (PC4).
    time_dim : :py:class:`str`, default "time"
        Name of the time coordinate dimension.
    ax : :py:class:`matplotlib.axes.Axes`, optional
        Target axes for plotting. If None, uses the current axes (plt.gca()).
    color : :py:class:`str` or tuple, default "blue"
        Color for the trajectory and points.
    start_text : :py:class:`str`, default "START"
        Text annotation for the trajectory starting point.
    add_start_text : :py:class:`bool`, default True
        Whether to display the start point annotation.

    Returns
    -------
    :py:class:`matplotlib.axes.Axes`
        Configured axes object with the BSISO2 phase space trajectory.

    Notes
    -----
    - Recommended to use with `draw_bsiso2_phase_space_basemap` for the background diagram.
    - Daily data is recommended for smooth trajectory visualization.

    Examples
    --------
    >>> import xarray as xr
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(figsize=(8, 8))
    >>> draw_bsiso2_phase_space_basemap(ax=ax)
    >>> draw_bsiso2_phase_space(bsiso_data=ds, ax=ax, color="red")
    >>> plt.title("BSISO2 Phase Space Trajectory")
    >>> plt.show()
    """
    if ax is None:
        ax = plt.gca()

    # Plot trajectory lines
    ax.plot(bsiso_data[pc4_dim], bsiso_data[pc3_dim], color=color)
    # Plot data points
    ax.scatter(bsiso_data[pc4_dim], bsiso_data[pc3_dim], color=color, s=15)
    # Plot start point
    ax.scatter(
        bsiso_data[pc4_dim].isel({time_dim: 0}),
        bsiso_data[pc3_dim].isel({time_dim: 0}),
        s=15 * 5,
        c=color,
        marker="*",
    )
    # Add start point text
    if add_start_text:
        ax.text(
            bsiso_data[pc4_dim].isel({time_dim: 0}) + 0.1,
            bsiso_data[pc3_dim].isel({time_dim: 0}),
            start_text,
            ha="left",
            zorder=2,
        )
    return ax


def calc_bsiso1_phase(
    ds: xr.DataArray,
    amplitude_threshold: float = 1.5,
    pc1_dim: str = "PC1",
    pc2_dim: str = "PC2",
):
    """
    Calculate BSISO1 phase based on PC1 and PC2 values.

    Computes the BSISO1 phase (0-8) based on the principal components PC1 and PC2,
    following the classification rules defined by Lee et al. (2013). Adds a boolean
    variable 'is_event' to indicate if the amplitude (:math:`\\sqrt{\\mathrm{PC1}^2 + \\mathrm{PC2}^2}`) exceeds
    the threshold, and a ``'phase'`` variable for the BSISO phase.

    Parameters
    ----------
    ds : :py:class:`xarray.Dataset`
        Dataset containing PC1 and PC2 variables with a time dimension.
    amplitude_threshold : :py:class:`float`, default 1.5
        Threshold for event detection based on amplitude.
    pc1_dim : :py:class:`str`, default "PC1"
        Variable name for PC1 in the dataset.
    pc2_dim : :py:class:`str`, default "PC2"
        Variable name for PC2 in the dataset.

    Returns
    -------
    :py:class:`xarray.Dataset`
        Dataset with added 'phase' (int, 0-8) and 'is_event' (boolean) variables.

    Notes
    -----
    Phase classification:

    - Phase 0: Non-event (amplitude <= threshold)
    - Phase 1: :math:`\\mathrm{PC1} < 0, \\mathrm{PC2} < 0, \\mathrm{PC1} > \\mathrm{PC2}`
    - Phase 2: :math:`\\mathrm{PC1} < 0, \\mathrm{PC2} < 0, \\mathrm{PC1} < \\mathrm{PC2}`
    - Phase 3: :math:`\\mathrm{PC1} < 0, \\mathrm{PC2} > 0, |\\mathrm{PC1}| > \\mathrm{PC2}`
    - Phase 4: :math:`\\mathrm{PC1} < 0, \\mathrm{PC2} > 0, |\\mathrm{PC1}| < \\mathrm{PC2}`
    - Phase 5: :math:`\\mathrm{PC1} > 0, \\mathrm{PC2} > 0, \\mathrm{PC1} < \\mathrm{PC2}`
    - Phase 6: :math:`\\mathrm{PC1} > 0, \\mathrm{PC2} > 0, \\mathrm{PC1} > \\mathrm{PC2}`
    - Phase 7: :math:`\\mathrm{PC1} > 0, \\mathrm{PC2} < 0, \\mathrm{PC1} > |\\mathrm{PC2}|`
    - Phase 8: :math:`\\mathrm{PC1} > 0, \\mathrm{PC2} < 0, \\mathrm{PC1} < |\\mathrm{PC2}|`

    Examples
    --------
    >>> ds = calc_bsiso1_phase(ds, amplitude_threshold=1.5)
    >>> print(ds['phase'])
    """
    # Calculate amplitude
    amplitude = np.sqrt(ds[pc1_dim] ** 2 + ds[pc2_dim] ** 2)

    # Add is_event variable
    ds["is_event"] = amplitude > amplitude_threshold

    # Initialize phase array
    phase = xr.DataArray(
        np.zeros_like(ds[pc1_dim], dtype=int),
        coords=ds[pc1_dim].coords,
        dims=ds[pc1_dim].dims,
        name="phase",
    )

    # Event condition mask
    event_mask = ds["is_event"]

    # Phase classification logic
    phase = phase.where(
        ~(
            event_mask
            & (ds[pc1_dim] < 0)
            & (ds[pc2_dim] < 0)
            & (ds[pc1_dim] > ds[pc2_dim])
        ),
        1,
    )
    phase = phase.where(
        ~(
            event_mask
            & (ds[pc1_dim] < 0)
            & (ds[pc2_dim] < 0)
            & (ds[pc1_dim] < ds[pc2_dim])
        ),
        2,
    )
    phase = phase.where(
        ~(
            event_mask
            & (ds[pc1_dim] < 0)
            & (ds[pc2_dim] > 0)
            & (np.abs(ds[pc1_dim]) > ds[pc2_dim])
        ),
        3,
    )
    phase = phase.where(
        ~(
            event_mask
            & (ds[pc1_dim] < 0)
            & (ds[pc2_dim] > 0)
            & (np.abs(ds[pc1_dim]) < ds[pc2_dim])
        ),
        4,
    )
    phase = phase.where(
        ~(
            event_mask
            & (ds[pc1_dim] > 0)
            & (ds[pc2_dim] > 0)
            & (ds[pc1_dim] < ds[pc2_dim])
        ),
        5,
    )
    phase = phase.where(
        ~(
            event_mask
            & (ds[pc1_dim] > 0)
            & (ds[pc2_dim] > 0)
            & (ds[pc1_dim] > ds[pc2_dim])
        ),
        6,
    )
    phase = phase.where(
        ~(
            event_mask
            & (ds[pc1_dim] > 0)
            & (ds[pc2_dim] < 0)
            & (ds[pc1_dim] > np.abs(ds[pc2_dim]))
        ),
        7,
    )
    phase = phase.where(
        ~(
            event_mask
            & (ds[pc1_dim] > 0)
            & (ds[pc2_dim] < 0)
            & (ds[pc1_dim] < np.abs(ds[pc2_dim]))
        ),
        8,
    )

    # Add phase to dataset
    ds["phase"] = phase
    return ds


def calc_bsiso2_phase(
    ds: xr.DataArray,
    amplitude_threshold: float = 1.5,
    pc3_dim: str = "PC3",
    pc4_dim: str = "PC4",
):
    """
    Calculate BSISO2 phase based on PC3 and PC4 values.

    Wrapper function for `calc_bsiso1_phase` to compute BSISO2 phases (0-8) using
    principal components PC3 and PC4, which capture higher-frequency (10-20 day)
    intraseasonal variability.

    Parameters
    ----------
    ds : :py:class:`xarray.Dataset`
        Dataset containing PC3 and PC4 variables with a time dimension.
    amplitude_threshold : :py:class:`float`, default 1.5
        Threshold for event detection based on amplitude.
    pc3_dim : :py:class:`str`, default "PC3"
        Variable name for PC3 in the dataset.
    pc4_dim : :py:class:`str`, default "PC4"
        Variable name for PC4 in the dataset.

    Returns
    -------
    :py:class:`xarray.Dataset`
        Dataset with added 'phase' (int, 0-8) and 'is_event' (boolean) variables.

    Notes
    -----
    - Delegates to `calc_bsiso1_phase` with PC3 and PC4 as inputs.
    - BSISO2 phases follow the same classification rules as BSISO1 but use different PCs.

    Examples
    --------
    >>> ds = calc_bsiso2_phase(ds, amplitude_threshold=1.5)
    >>> print(ds['phase'])
    """
    return calc_bsiso1_phase(
        ds=ds, amplitude_threshold=amplitude_threshold, pc1_dim=pc3_dim, pc2_dim=pc4_dim
    )


def get_times_for_phase(
    ds: xr.Dataset, phase_value: int, time_dim: str = "time", phase_dim: str = "phase"
):
    """
    Retrieve time points for a specified BSISO phase.

    Extracts time coordinates from the dataset where the 'phase' variable matches the
    specified phase value (0-8), useful for composite analysis of BSISO events.

    Parameters
    ----------
    ds : :py:class:`xarray.Dataset`
        Dataset containing 'phase' variable with a time dimension.
    phase_value : :py:class:`int`
        Target phase value to query (0-8).
    time_dim : :py:class:`str`, default "time"
        Name of the time coordinate dimension.
    phase_dim : :py:class:`str`, default "phase"
        Name of the phase variable.

    Returns
    -------
    :py:class:`xarray.DataArray`
        Array of time coordinates where the phase equals ``phase_value``.

    Notes
    -----
    - Returns an empty DataArray if no time points match the specified phase.
    - Assumes the 'phase' variable exists in the dataset.

    Examples
    --------
    >>> times = get_times_for_phase(ds, phase_value=1)
    >>> print(times)
    """
    return ds[time_dim].where(ds[phase_dim] == phase_value, drop=True)


def calc_bsiso_analysis(
    olr_data: xr.DataArray,
    u850_data: xr.DataArray,
    v850_data: xr.DataArray,
    daily_cycle_mean_time_range: slice = slice(None, None),
    extract_time_range: slice = slice(None, None),
    harmonics_num: int = 3,
    threshold: float = 0.05,
    time_dim: str = "time",
    lat_dim: str = "lat",
    lon_dim: str = "lon",
) -> DataNode:
    """
    Perform BSISO analysis using OLR and 850-hPa wind data.

    Computes multivariate BSISO indices (BSISO1 and BSISO2) using outgoing longwave
    radiation (OLR) and 850-hPa zonal (u850) and meridional (v850) winds, following
    the methodology of Lee et al. (2013). The function processes data to remove
    seasonal cycles and low-frequency signals, performs multivariate EOF analysis,
    and derives phase composites for BSISO1 and BSISO2.

    Parameters
    ----------
    olr_data : :py:class:`xarray.DataArray`
        Outgoing longwave radiation data with time, latitude, and longitude dimensions.
    u850_data : :py:class:`xarray.DataArray`
        850-hPa zonal wind data with the same dimensions as ``olr_data``.
    v850_data : :py:class:`xarray.DataArray`
        850-hPa meridional wind data with the same dimensions as ``olr_data``.
    daily_cycle_mean_time_range : :py:class:`slice`, default ``slice(None, None)``
        Time range for computing the daily annual cycle mean.
    extract_time_range : :py:class:`slice`, default ``slice(None, None)``
        Time range for extracting data for analysis.
    harmonics_num : :py:class:`int`, default 3
        Number of harmonics for smoothing the daily annual cycle.
    threshold : :py:class:`float`, default 0.05
        P-value threshold for statistical significance in composite analysis.
    time_dim : :py:class:`str`, default "time"
        Name of the time dimension.
    lat_dim : :py:class:`str`, default "lat"
        Name of the latitude dimension.
    lon_dim : :py:class:`str`, default "lon"
        Name of the longitude dimension.

    Returns
    -------
    :py:class:`easyclimate.DataNode`
        Hierarchical data structure containing:

        - Explained variance ratios
        - Principal components (PC1-PC4)
        - Normalized PCs
        - EOF patterns for OLR and winds
        - Lead-lag correlation coefficients
        - Phase information
        - Composite analysis results for BSISO1 and BSISO2 phases

    Notes
    -----
    - The analysis focuses on May-October data to capture the boreal summer monsoon season.
    - BSISO1 represents 30-60 day oscillations, while BSISO2 captures 10-20 day oscillations.
    - Progress indicators are logged to track major computational steps.

    References
    ----------
    - Lee, JY. (이준이), Wang, B., Wheeler, M.C. et al. Real-time multivariate indices for the boreal summer intraseasonal oscillation over the Asian summer monsoon region. Clim Dyn 40, 493–509 (2013). https://doi.org/10.1007/s00382-012-1544-4
    - Wheeler, M. C., & Hendon, H. H. (2004). An All-Season Real-Time Multivariate MJO Index: Development of an Index for Monitoring and Prediction. Monthly Weather Review, 132(8), 1917-1932. https://doi.org/10.1175/1520-0493(2004)132<1917:AARMMI>2.0.CO;2

    Examples
    --------
    >>> result = calc_bsiso_analysis(olr_data, u850_data, v850_data)
    >>> print(result["Phase/PC1_2"])

    .. seealso::

        - https://iprc.soest.hawaii.edu/users/jylee/bsiso/
    """
    logger.info("Starting BSISO analysis...")
    result_node = DataNode()

    olr_data = sort_ascending_latlon_coordinates(
        olr_data, lat_dim=lat_dim, lon_dim=lon_dim
    ).sel({lon_dim: slice(40, 160), lat_dim: slice(-10, 40)})
    u850_data = sort_ascending_latlon_coordinates(
        u850_data, lat_dim=lat_dim, lon_dim=lon_dim
    ).sel({lon_dim: slice(40, 160), lat_dim: slice(-10, 40)})
    v850_data = sort_ascending_latlon_coordinates(
        v850_data, lat_dim=lat_dim, lon_dim=lon_dim
    ).sel({lon_dim: slice(40, 160), lat_dim: slice(-10, 40)})

    # 创建进度条配置
    progress_columns = [
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeRemainingColumn(),
    ]

    with Progress(*progress_columns) as progress:
        # Step 1: Remove daily annual cycle
        logger.info("Removing smooth daily annual cycle...")
        olr_daily_smoothed = remove_smooth_daily_annual_cycle_mean(
            olr_data,
            daily_cycle_mean_time_range=daily_cycle_mean_time_range,
            extract_time_range=extract_time_range,
            harmonics_num=harmonics_num,
            dim=time_dim,
        )
        u850_daily_smoothed = remove_smooth_daily_annual_cycle_mean(
            u850_data,
            daily_cycle_mean_time_range=daily_cycle_mean_time_range,
            extract_time_range=extract_time_range,
            harmonics_num=harmonics_num,
            dim=time_dim,
        )
        v850_daily_smoothed = remove_smooth_daily_annual_cycle_mean(
            v850_data,
            daily_cycle_mean_time_range=daily_cycle_mean_time_range,
            extract_time_range=extract_time_range,
            harmonics_num=harmonics_num,
            dim=time_dim,
        )

        # Step 2: Remove low-frequency signals
        logger.info("Removing low-frequency signals...")
        olr_anomaly = remove_low_frequency_signal(
            olr_daily_smoothed, window=120, center=False, time_dim=time_dim
        )
        u850_anomaly = remove_low_frequency_signal(
            u850_daily_smoothed, window=120, center=False, time_dim=time_dim
        )
        v850_anomaly = remove_low_frequency_signal(
            v850_daily_smoothed, window=120, center=False, time_dim=time_dim
        )

        # Step 3: Normalize anomalies
        logger.info("Normalizing anomalies...")
        olr_factor = olr_anomaly.std(dim=time_dim).mean().data
        u850_factor = u850_anomaly.std(dim=time_dim).mean().data
        v850_factor = v850_anomaly.std(dim=time_dim).mean().data

        olr_anomaly_normalized = olr_anomaly / olr_factor
        u850_anomaly_normalized = u850_anomaly / u850_factor
        v850_anomaly_normalized = v850_anomaly / v850_factor

        # Step 4: Extract May-October data
        logger.info("Extracting May-October data...")
        olr_anomaly_may2oct = get_specific_months_data(
            olr_anomaly_normalized, [5, 6, 7, 8, 9, 10]
        )
        u850_anomaly_may2oct = get_specific_months_data(
            u850_anomaly_normalized, [5, 6, 7, 8, 9, 10]
        )
        v850_anomaly_may2oct = get_specific_months_data(
            v850_anomaly_normalized, [5, 6, 7, 8, 9, 10]
        )

        # Step 5: Perform multivariate EOF analysis
        logger.info("Performing multivariate EOF analysis...")
        model = get_EOF_model(
            [olr_anomaly_may2oct.fillna(0), u850_anomaly_may2oct.fillna(0)],
            lat_dim=lat_dim,
            lon_dim=lon_dim,
            use_coslat=True,
        )
        meof_analysis_result = calc_EOF_analysis(model, PC_normalized=False)
        result_node["explained_variance_ratio"] = meof_analysis_result[
            "explained_variance_ratio"
        ]

        # Step 6: Adjust EOF coefficients
        logger.info("Adjusting EOF coefficients...")
        arr1 = meof_analysis_result["EOF/var0"].sel(mode=1).mean(dim=lat_dim)
        arr2 = meof_analysis_result["EOF/var0"].sel(mode=2)
        arr3 = meof_analysis_result["EOF/var0"].sel(mode=3).mean(dim=lat_dim)
        arr4 = (
            meof_analysis_result["EOF/var0"]
            .sel(mode=4)
            .sel({lon_dim: slice(100, None)})
            .mean(dim=lon_dim)
        )

        bsiso1_coeff = (
            1
            if arr1.sel({lon_dim: 80}, method="nearest") > 0
            and arr1.sel({lon_dim: 140}, method="nearest") < 0
            else -1
        )
        bsiso2_coeff = (
            1
            if arr2.sel({lon_dim: 100, lat_dim: 10}, method="nearest") < 0
            and arr2.sel({lon_dim: 140, lat_dim: 20}, method="nearest") > 0
            else -1
        )
        bsiso3_coeff = 1 if arr3.sel({lon_dim: 100}, method="nearest") > 0 else -1
        bsiso4_coeff = (
            1
            if arr4.sel({lat_dim: 10}, method="nearest") > 0
            and arr4.sel({lat_dim: 30}, method="nearest") < 0
            else -1
        )

        # Step 7: Compute principal components
        logger.info("Computing principal components...")
        pc1 = meof_analysis_result["PC"].sel(mode=1) * bsiso1_coeff
        pc2 = meof_analysis_result["PC"].sel(mode=2) * bsiso2_coeff
        pc3 = meof_analysis_result["PC"].sel(mode=3) * bsiso3_coeff
        pc4 = meof_analysis_result["PC"].sel(mode=4) * bsiso4_coeff

        result_node["PC/pc1"] = pc1
        result_node["PC/pc2"] = pc2
        result_node["PC/pc3"] = pc3
        result_node["PC/pc4"] = pc4

        # Step 8: Normalize principal components
        logger.info("Normalizing principal components...")
        pc12_normalized = xr.Dataset()
        pc12_normalized["PC1"] = pc1 / pc1.std(dim=time_dim)
        pc12_normalized["PC2"] = pc2 / pc2.std(dim=time_dim)
        pc34_normalized = xr.Dataset()
        pc34_normalized["PC3"] = pc3 / pc3.std(dim=time_dim)
        pc34_normalized["PC4"] = pc4 / pc4.std(dim=time_dim)

        result_node["PC_normalized/pc12"] = pc12_normalized
        result_node["PC_normalized/pc34"] = pc34_normalized

        # Step 9: Store EOF patterns
        logger.info("Storing EOF patterns...")
        result_node["EOF/olr1"] = (
            meof_analysis_result["EOF/var0"].sel(mode=1) * bsiso1_coeff
        )
        result_node["EOF/olr2"] = (
            meof_analysis_result["EOF/var0"].sel(mode=2) * bsiso2_coeff
        )
        result_node["EOF/olr3"] = (
            meof_analysis_result["EOF/var0"].sel(mode=3) * bsiso3_coeff
        )
        result_node["EOF/olr4"] = (
            meof_analysis_result["EOF/var0"].sel(mode=4) * bsiso4_coeff
        )

        wind_eof_task = progress.add_task("[cyan]Computing wind EOFs...", total=4)
        for mode in range(1, 5):
            uvdata = xr.Dataset()
            uvdata["u"] = calc_corr_spatial(
                u850_anomaly_may2oct, meof_analysis_result["PC"].sel(mode=mode)
            ).reg_coeff
            uvdata["v"] = calc_corr_spatial(
                v850_anomaly_may2oct, meof_analysis_result["PC"].sel(mode=mode)
            ).reg_coeff
            result_node[f"EOF/uv850_{mode}"] = uvdata * locals()[f"bsiso{mode}_coeff"]
            progress.update(
                wind_eof_task,
                advance=1,
                description=f"[cyan]Computing wind EOF {mode}/4",
            )

        # Step 10: Compute lead-lag correlations
        logger.info("Computing lead-lag correlations...")
        pcs = {"PC1": pc1, "PC2": pc2, "PC3": pc3, "PC4": pc4}
        pairs = [("PC1_vs_PC2", "PC1", "PC2"), ("PC3_vs_PC4", "PC3", "PC4")]
        max_lag = 30
        corr_da = calc_lead_lag_correlation_coefficients(
            pcs=pcs, pairs=pairs, max_lag=max_lag
        )
        result_node["lead_lag_correlation_coefficients"] = corr_da

        # Step 11: Calculate BSISO phases
        logger.info("Calculating BSISO phases...")
        dpc1_2 = xr.Dataset()
        dpc1_2["PC1"] = pc1
        dpc1_2["PC2"] = pc2
        phase1_2 = calc_bsiso1_phase(
            dpc1_2, amplitude_threshold=1.5, pc1_dim="PC1", pc2_dim="PC2"
        )
        result_node["Phase/PC1_2"] = phase1_2

        dpc3_4 = xr.Dataset()
        dpc3_4["PC3"] = pc3
        dpc3_4["PC4"] = pc4
        phase3_4 = calc_bsiso1_phase(
            dpc3_4, amplitude_threshold=1.5, pc1_dim="PC3", pc2_dim="PC4"
        )
        result_node["Phase/PC3_4"] = phase3_4

        # Step 12: Perform composite analysis for BSISO1
        logger.info("Performing composite analysis for BSISO1...")
        bsiso1_task = progress.add_task("[green]BSISO1 composite analysis...", total=8)
        for phase in range(1, 9):
            olr_anomaly_phase = olr_daily_smoothed.where(
                get_times_for_phase(phase1_2, phase, time_dim=time_dim)
            )
            uv850_anomaly_phase = xr.Dataset(
                {"u": u850_daily_smoothed, "v": v850_daily_smoothed}
            ).where(get_times_for_phase(phase1_2, phase, time_dim=time_dim))

            olr_pvalue = calc_ttestSpatialPattern_spatial(
                olr_daily_smoothed.fillna(0), olr_anomaly_phase.fillna(0)
            ).pvalue
            u850_pvalue = calc_ttestSpatialPattern_spatial(
                u850_daily_smoothed.fillna(0), uv850_anomaly_phase["u"].fillna(0)
            ).pvalue

            uv850_mean = uv850_anomaly_phase.mean(dim=time_dim).where(
                u850_pvalue < threshold
            )
            uv850_mean["olr"] = olr_anomaly_phase.mean(dim=time_dim).where(
                olr_pvalue < threshold
            )
            result_node[f"composite_analysis/bsiso1/phase{phase}"] = uv850_mean
            progress.update(
                bsiso1_task, advance=1, description=f"[green]BSISO1 phase {phase}/8"
            )

        # Step 13: Perform composite analysis for BSISO2
        logger.info("Performing composite analysis for BSISO2...")
        bsiso2_task = progress.add_task("[blue]BSISO2 composite analysis...", total=8)
        for phase in range(1, 9):
            olr_anomaly_phase = olr_daily_smoothed.where(
                get_times_for_phase(phase3_4, phase, time_dim=time_dim)
            )
            uv850_anomaly_phase = xr.Dataset(
                {"u": u850_daily_smoothed, "v": v850_daily_smoothed}
            ).where(get_times_for_phase(phase3_4, phase, time_dim=time_dim))

            olr_pvalue = calc_ttestSpatialPattern_spatial(
                olr_daily_smoothed.fillna(0), olr_anomaly_phase.fillna(0)
            ).pvalue
            u850_pvalue = calc_ttestSpatialPattern_spatial(
                u850_daily_smoothed.fillna(0), uv850_anomaly_phase["u"].fillna(0)
            ).pvalue

            uv850_mean = uv850_anomaly_phase.mean(dim=time_dim).where(
                u850_pvalue < threshold
            )
            uv850_mean["olr"] = olr_anomaly_phase.mean(dim=time_dim).where(
                olr_pvalue < threshold
            )
            result_node[f"composite_analysis/bsiso2/phase{phase}"] = uv850_mean
            progress.update(
                bsiso2_task, advance=1, description=f"[blue]BSISO2 phase {phase}/8"
            )

        # 确保所有进度条完成
        progress.update(wind_eof_task, completed=4, refresh=True)
        progress.update(bsiso1_task, completed=8, refresh=True)
        progress.update(bsiso2_task, completed=8, refresh=True)
        logger.info("BSISO analysis completed.")
    return result_node
