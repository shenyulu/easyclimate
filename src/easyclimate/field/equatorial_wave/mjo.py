"""
Madden Julian Oscillation (MJO)
"""

import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
import matplotlib.lines as lines


__all__ = ["draw_mjo_phase_space_basemap", "draw_mjo_phase_space"]


def draw_mjo_phase_space_basemap(
    ax=None, add_phase_text: bool = True, add_location_text: bool = True
):
    """
    Create a base map for Madden-Julian Oscillation (MJO) phase space visualization.

    This function generates the standardized background diagram for MJO analysis,
    including phase division lines (1-8), amplitude circles, and optional geographical
    annotations. The diagram uses RMM1/RMM2 indices as axes in a Cartesian coordinate
    system (-4 to 4 range).

    Parameters
    ----------
    ax : :py:class:`matplotlib.axes.Axes`, optional
        Target axes object for plotting. If None, uses current axes (gca()).
    add_phase_text : :py:class:`bool`, default: True
        Whether to annotate MJO phase numbers (1-8) on the diagram.
    add_location_text : :py:class:`bool`, default: True
        Whether to annotate geographical regions corresponding to MJO phases:
        - Pacific Ocean (top)
        - Indian Ocean (bottom)
        - Western Hemisphere/Africa (left)
        - Maritime Continent (right)

    Returns
    -------
    :py:class:`matplotlib.axes.Axes`
        The configured axes object with MJO phase space basemap elements.

    Notes
    -----
    The diagram includes:
    1. 8 radial lines dividing phases (45° intervals)
    2. Unit circle indicating amplitude threshold
    3. Axes labeled as RMM1 (x-axis) and RMM2 (y-axis)

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> ecl.field.equatorial_wave.draw_mjo_phase_space_basemap(ax=ax)
    >>> plt.show()

    .. seealso::

        https://www.ncl.ucar.edu/Applications/mjoclivar.shtml

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_mjo_phase.py
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
            3.2,
            "Western\nPacific",
            ha="center",
            bbox=dict(facecolor="white", edgecolor="white"),
        )
        ax.text(
            0,
            -3.8,
            "Indian\nOcean",
            ha="center",
            bbox=dict(facecolor="white", edgecolor="white"),
        )
        ax.text(
            -3.8,
            0,
            "West. Hem.\nandAfrica",
            va="center",
            rotation="vertical",
            bbox=dict(facecolor="white", edgecolor="white"),
        )
        ax.text(
            3.2,
            0,
            "Maritime\nContinent",
            va="center",
            rotation=270,
            bbox=dict(facecolor="white", edgecolor="white"),
        )

    ax.set_xlabel("RMM1")
    ax.set_ylabel("RMM2")


def draw_mjo_phase_space(
    mjo_data: xr.Dataset,
    rmm1_dim: str = "RMM1",
    rmm2_dim: str = "RMM2",
    time_dim: str = "time",
    ax=None,
    color="blue",
    start_text="START",
    add_start_text: bool = True,
):
    """
    Visualize Madden-Julian Oscillation (MJO) phase space trajectory from RMM indices.

    Plots the temporal evolution of MJO phases as a parametric curve in RMM1-RMM2 space,
    with optional starting point marker. Combines scatter points (individual observations)
    with connecting lines to show temporal progression.

    Parameters
    ----------
    mjo_data : :py:class:`xarray.Dataset`
        Input dataset containing RMM indices. Must include:
        - RMM1 component (real-time multivariate MJO index 1)
        - RMM2 component (real-time multivariate MJO index 2)
        - Time coordinate
    rmm1_dim : :py:class:`str`, default: "RMM1"
        Variable name for RMM1 component in the dataset
    rmm2_dim : :py:class:`str`, default: "RMM2"
        Variable name for RMM2 component in the dataset
    time_dim : :py:class:`str`, default: "time"
        Time coordinate dimension name
    ax : :py:class:`matplotlib.axes.Axes`, optional
        Target axes for plotting. Creates new axes if None.
    color : :py:class:`str` or tuple, default: "blue"
        Color specification for trajectory and points
    start_text : :py:class:`str`, default: "START"
        Annotation text for trajectory starting point
    add_start_text : :py:class:`bool`, default: True
        Toggle for displaying start point annotation

    Returns
    -------
    :py:class:`matplotlib.axes.Axes`
        Configured axes object with MJO phase space trajectory

    Returns
    -------
    :py:class:`matplotlib.axes.Axes`
        Configured axes object with MJO phase space trajectory

    Examples
    --------
    >>> import xarray as xr
    >>> import matplotlib.pyplot as plt
    >>> # Load RMM indices dataset
    >>> mjo_ds = xr.open_dataset('http://iridl.ldeo.columbia.edu/SOURCES/.BoM/.MJO/.RMM/dods', decode_times=False)
    >>> T = mjo_ds.T.values
    >>> mjo_ds['T'] = pd.date_range("1974-06-01", periods=len(T))
    >>> mjo_ds = ecl.utility.get_compress_xarraydata(mjo_ds)
    >>> fig, ax = plt.subplots(figsize=(8,8))
    >>> ecl.field.equatorial_wave.draw_mjo_phase_space_basemap()
    >>> ecl.field.equatorial_wave.draw_mjo_phase_space(
    ...:    mjo_data = mjo_data.sel(time = slice('2024-12-01', '2024-12-31')),
    ...:    rmm1_dim = "RMM1",
    ...:    rmm2_dim = "RMM2",
    ...:    time_dim = "time"rmm_data, ax=ax, color="red")
    >>> plt.title("MJO Phase Space Trajectory")
    >>> plt.show()

    .. note::
        1. Recommended to use with :py:func:`draw_mjo_phase_space_basemap` for basic map.
        2. Temporal resolution affects trajectory smoothness (daily data recommended).

    .. seealso::

        https://www.ncl.ucar.edu/Applications/mjoclivar.shtml

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_mjo_phase.py
    """
    if ax is None:
        ax = plt.gca()

    mjo_data.plot.scatter(x=rmm1_dim, y=rmm2_dim, ax=ax, color=color)
    ax.plot(mjo_data[rmm1_dim], mjo_data[rmm2_dim], color=color)
    ax.text(
        mjo_data[rmm1_dim].isel({time_dim: 0}) + 0.1,
        mjo_data[rmm2_dim].isel({time_dim: 0}),
        start_text if add_start_text else "",
        ha="left",
        zorder=2,
    )
