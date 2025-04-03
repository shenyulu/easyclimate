"""
Madden Julian Oscillation (MJO)
"""

import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
import matplotlib.lines as lines


def draw_mjo_phase_space_basemap(
    ax=None, add_phase_text: bool = True, add_location_text: bool = True
):
    """ """
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
        ax.text(1, 3, "Phase 6", size="x-large", weight="semibold")
        ax.text(-2, 3, "Phase 7", size="x-large", weight="semibold")
        ax.text(2.8, 1, "Phase 5", size="x-large", weight="semibold", ha="center")
        ax.text(-2.8, 1, "Phase 8", size="x-large", weight="semibold", ha="center")

        ax.text(1, -3, "Phase 3", size="x-large", weight="semibold")
        ax.text(-2, -3, "Phase 2", size="x-large", weight="semibold")
        ax.text(2.8, -1, "Phase 4", size="x-large", weight="semibold", ha="center")
        ax.text(-2.8, -1, "Phase 1", size="x-large", weight="semibold", ha="center")

    if add_location_text:
        ax.text(
            0,
            3.7,
            "Pacific Ocean",
            ha="center",
            bbox=dict(facecolor="white", edgecolor="white"),
        )
        ax.text(
            0,
            -3.8,
            "Indian Ocean",
            ha="center",
            bbox=dict(facecolor="white", edgecolor="white"),
        )
        ax.text(
            -3.8,
            0,
            "West. Hem., Africa",
            va="center",
            rotation="vertical",
            bbox=dict(facecolor="white", edgecolor="white"),
        )
        ax.text(
            3.7,
            0,
            "Maritime Continent",
            va="center",
            rotation="vertical",
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
