from __future__ import annotations

from typing import Any
from ....backend import gill_simulate

import numpy as np
import xarray as xr

__all__ = ["simulate"]

CASES = {
    "gill_sym": "Gill: symmetric heating",
    "gill_asym": "Gill: asymmetric heating",
    "kelvin_wave": "Kelvin wave",
    "gauss_pert": "Gaussian-shape initial perturbation",
    "matsuno_heating": "Matsuno: forced stationary motion",
    "gill_sym_gaussian_x": "Gill: symmetric heating with Gaussian shape in zonal direction",
    "gill_sym_asym": "Gill: symmetric heating + asymmetric heating",
    "gill_delta_y1": "Gill: delta heating at y = 1",
}

_ARRAY_KEYS = {
    "x",
    "y",
    "z",
    "q",
    "eta",
    "u",
    "v",
    "v_staggered",
    "time",
    "eta_history",
    "u_history",
    "v_history",
}


def _as_numpy_dict(raw: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = dict(raw)
    for key in _ARRAY_KEYS:
        if key in out:
            out[key] = np.asarray(out[key], dtype=float)
    return out


def _as_xarray_dataset(result: dict[str, Any]) -> Any:
    """Convert a solver result to a labelled xarray dataset."""
    coords = {
        "x": result["x"],
        "y": result["y"],
        "z": result["z"],
        "y_staggered": result["y"][:-1],
    }

    time = result["time"]
    if time.size:
        coords["time"] = time

    # Solver output: (x, y) -> xarray convention: (y, x)
    data_vars = {name: (("y", "x"), result[name].T) for name in ("q", "eta", "u", "v")}

    # (x, y_staggered) -> (y_staggered, x)
    data_vars["v_staggered"] = (
        ("y_staggered", "x"),
        result["v_staggered"].T,
    )

    if time.size:
        # (time, x, y) -> (time, y, x)
        for name in ("eta_history", "u_history", "v_history"):
            data_vars[name] = (
                ("time", "y", "x"),
                result[name].transpose(0, 2, 1),
            )

    return xr.Dataset(
        data_vars=data_vars,
        coords=coords,
        attrs={key: value for key, value in result.items() if key not in _ARRAY_KEYS},
    )


def simulate(
    case: str = "gill_sym",
    *,
    t: float | None = None,
    nt: int | None = None,
    nx: int = 128,
    ny: int = 60,
    nz: int = 20,
    lx: float | None = None,
    ly: float | None = None,
    lz: float | None = None,
    damping: float | None = None,
    store_every: int = 0,
    kernel: str = "reyon",
) -> dict[str, Any] | Any:
    """Integrate the linear damped shallow-water model on an equatorial beta plane.

    The model advances the zonal and meridional velocity perturbations ``u``
    and ``v`` together with the geopotential-height perturbation ``eta``.
    Their evolution includes linear pressure-gradient, Coriolis, divergence,
    heating, and Rayleigh-damping terms. This makes the function suitable for
    examining the idealized tropical circulation excited by prescribed heating
    or by an initial height and wind perturbation; it is not a full nonlinear
    atmospheric model.

    The horizontal problem is periodic in ``x``. The solver Fourier transforms
    each field along that direction, precomputes the implicit linear transform
    for every zonal wavenumber, advances the independent spectral modes, and
    transforms the final state back to physical space. The meridional velocity
    is held on a staggered grid with ``ny`` points. ``u``, ``eta``, forcing,
    and the returned interpolated ``v`` use the ``ny + 1`` point ``y`` grid.
    The returned ``v`` is the average of adjacent staggered values, with zero
    values outside the two meridional boundaries.

    ``case`` selects the forcing or initial condition. Heating experiments
    begin from rest; free-wave experiments begin with the listed perturbation
    and no heating:

    - ``gill_sym``: compact cosine heating in ``x`` multiplied by a Gaussian
      profile symmetric about the equator.
    - ``gill_asym``: the same compact zonal heating with an odd, Gaussian-
      weighted meridional profile.
    - ``kelvin_wave``: an equatorially trapped cosine perturbation in ``eta``
      and ``u``; ``v`` and heating are initially zero.
    - ``gauss_pert``: a two-dimensional Gaussian initial ``eta`` perturbation
      with initially zero wind and heating.
    - ``matsuno_heating``: sinusoidal zonal heating with a symmetric Gaussian
      meridional envelope.
    - ``gill_sym_gaussian_x``: symmetric heating localized by a Gaussian in
      the zonal direction.
    - ``gill_sym_asym``: the compact zonal envelope multiplied by the sum of
      symmetric and asymmetric Gaussian meridional components.
    - ``gill_delta_y1``: the symmetric Gill heating plus a unit-width band of
      meridionally localized heating centred at ``y = 1``.

    Unless overridden, the integration uses ``t=50``, ``nt=500``,
    ``lx=12*pi``, ``ly=5*pi``, ``lz=pi``, and ``damping=0.1``. The Kelvin-wave
    experiment instead defaults to ``t=10``, ``nt=100``, and zero damping;
    the Gaussian-perturbation experiment defaults to zero damping; and the
    forced stationary-motion experiment defaults to ``damping=0.2``.

    Parameters
    ----------
    case : str, default: ``"gill_sym"``
        One of ``gill_sym``, ``gill_asym``, ``kelvin_wave``, ``gauss_pert``,
        ``matsuno_heating``, ``gill_sym_gaussian_x``, ``gill_sym_asym``,
        or ``gill_delta_y1``. The backend also accepts the descriptive aliases
        ``symmetric_heating``, ``asymmetric_heating``, ``gaussian_initial``,
        ``forced_stationary``, ``symmetric_gaussian_x``,
        ``symmetric_asymmetric``, and ``delta_y1``.
    t : float, optional
        Total nondimensional integration time. The time step is ``t / nt``.
        The selected case determines the default.
    nt : int, optional
        Number of equal time steps. Must be positive. The selected case
        determines the default.
    nx, ny, nz : int, default: 128, 60, 20
        Number of zonal grid points, meridional staggered-grid points, and
        vertical diagnostic levels, respectively. ``nx`` and ``ny`` control
        the dynamics; ``nz`` only controls the returned ``z`` coordinate used
        by :func:`plot_cross_sections` and must be at least 2.
    lx, ly, lz : float, optional
        Nondimensional domain lengths. They define ``dx=lx/nx``,
        ``dy=ly/(ny+1)``, and ``dz=lz/(nz-1)``. ``lz`` is used only when
        constructing vertical diagnostic modes.
    damping : float, optional
        Linear damping coefficient applied to the prognostic perturbation
        fields. The selected case determines the default.
    store_every : int, default: 0
        Save ``eta``, ``u``, and ``v`` snapshots at this interval in time
        steps. A final snapshot is included even when ``nt`` is not divisible
        by this interval. Use ``0`` to omit all time-history variables.
    kernel : {``"serial"``, ``"reyon"``}, default: ``"reyon"``
        Use ``"serial"`` for one CPU core or ``"reyon"`` to evaluate
        independent zonal Fourier modes in parallel. Both kernels perform the
        same calculation. The ``"reyon"`` kernel uses at most six CPU cores.

    Returns
    -------
    xarray.Dataset
        Labelled simulation result with the following data variables:

        - ``q(y, x)``: imposed heating, or zeros for the two free-wave cases.
        - ``eta(y, x)``: final geopotential-height perturbation.
        - ``u(y, x)``: final zonal-velocity perturbation.
        - ``v(y, x)``: final meridional velocity interpolated from the
          staggered grid to the ``eta`` grid.
        - ``v_staggered(y_staggered, x)``: final native staggered meridional
          velocity.
        - ``eta_history(time, y, x)``, ``u_history(time, y, x)``, and
          ``v_history(time, y, x)``: saved snapshots, included only when
          ``store_every`` is positive.

        Coordinates are ``x``, ``y``, ``y_staggered``, ``z``, and, when
        snapshots are saved, ``time``. Attributes retain the resolved
        ``case`` label, integration and grid settings, selected ``kernel``,
        and the derived increments ``dt``, ``dx``, ``dy``, and ``dz``.

    Raises
    ------
    ValueError
        If the case or kernel is unknown, or if ``nt``, ``nx``, or ``ny`` is
        zero, or ``nz`` is less than 2.

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

    raw = gill_simulate(
        case=case,
        t=t,
        nt=nt,
        nx=nx,
        ny=ny,
        nz=nz,
        lx=lx,
        ly=ly,
        lz=lz,
        damping=damping,
        store_every=store_every,
        kernel=kernel,
    )
    result = _as_numpy_dict(raw)
    return _as_xarray_dataset(result)
