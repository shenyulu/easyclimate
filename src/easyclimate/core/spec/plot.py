"""
Visualization for spherical harmonics analysis
"""

import numpy as np
import xarray as xr


__all__ = [
    "spectral_coefficients_to_matrix",
    "spectral_power_by_total_wavenumber",
    "plot_spectral_power",
    "plot_spectral_matrix",
    "plot_spectral_real_imag",
]


def spectral_coefficients_to_matrix(
    spec: xr.DataArray,
    value: str = "abs",
    eps: float = 1.0e-30,
    mode: str = "spec_dim",
) -> xr.DataArray:
    """
    Convert packed triangular spectral coefficients to an ``(m, n)`` matrix.

    The input coefficients are assumed to follow triangular spherical-harmonic
    ordering, where the packed dimension contains all pairs with
    ``0 <= m <= n <= ntrunc``. Cells outside the triangular domain
    (``n < m``) are filled with NaN in the returned matrix.

    Parameters
    ----------
    spec : xr.DataArray
        One-dimensional spectral coefficients. Select non-mode dimensions first,
        e.g. ``spec_data.isel(time=0, level=0)``.
    value : {"abs", "logabs", "power", "real", "imag"}
        Which quantity to place in the matrix.
    eps : float
        Small number used only for ``logabs``.
    mode : str, default: "spec_dim"
        Name of the packed spectral coefficient dimension.

    Returns
    -------
    xr.DataArray
        Two-dimensional coefficient matrix with dimensions ``("m", "n")``.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        visual_spec

    """
    if mode not in spec.dims:
        raise ValueError(f"The spec must include the 'mode' dimension ({mode}).")
    other_dims = tuple(d for d in spec.dims if d != mode)
    if other_dims:
        raise ValueError(
            f"Before drawing, please select the non-mode dimension: {other_dims}"
        )

    coef = np.asarray(spec.values)
    ntrunc = _infer_triangular_ntrunc(coef.size)
    ms, ns = _triangular_mode_indices(ntrunc)

    mat = np.full((ntrunc + 1, ntrunc + 1), np.nan, dtype=np.float64)
    if value == "abs":
        vals = np.abs(coef)
    elif value == "logabs":
        vals = np.log10(np.abs(coef) + eps)
    elif value == "power":
        vals = np.abs(coef) ** 2
    elif value == "real":
        vals = coef.real
    elif value == "imag":
        vals = coef.imag
    else:
        raise ValueError(
            "``value`` must be 'abs', 'logabs', 'power', 'real', or 'imag'"
        )

    mat[ms, ns] = vals
    return xr.DataArray(
        mat,
        dims=("m", "n"),
        coords={"m": np.arange(ntrunc + 1), "n": np.arange(ntrunc + 1)},
        name=f"{spec.name or 'spec'}_{value}_mn",
        attrs={"spectral_grid": "triangular", "ntrunc": ntrunc},
    )


def spectral_power_by_total_wavenumber(
    spec: xr.DataArray,
    mode: str = "spec_dim",
) -> xr.DataArray:
    """
    Sum spectral power over zonal wavenumber for each total wavenumber.

    The returned one-dimensional spectrum is
    :math:`E(n) = \\sum_m |\\hat{a}_{n,m}|^2` for packed triangular
    spherical-harmonic coefficients.

    Parameters
    ----------
    spec : xr.DataArray
        One-dimensional spectral coefficients. Select non-mode dimensions first,
        e.g. ``spec_data.isel(time=0, level=0)``.
    mode : str, default: "spec_dim"
        Name of the packed spectral coefficient dimension.

    Returns
    -------
    xr.DataArray
        Spectral power indexed by total wavenumber ``n``.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        visual_spec

    """
    if mode not in spec.dims:
        raise ValueError(f"The spec must include the 'mode' dimension ({mode}).")
    other_dims = tuple(d for d in spec.dims if d != mode)
    if other_dims:
        raise ValueError(
            f"Please select the non-mode dimension before calculating the 1D spectrum: {other_dims}"
        )

    coef = np.asarray(spec.values)
    ntrunc = _infer_triangular_ntrunc(coef.size)
    _, ns = _triangular_mode_indices(ntrunc)
    power = np.bincount(ns, weights=np.abs(coef) ** 2, minlength=ntrunc + 1)
    return xr.DataArray(
        power,
        dims=("n",),
        coords={"n": np.arange(ntrunc + 1)},
        name=f"{spec.name or 'spec'}_power_n",
        attrs={"spectral_grid": "triangular", "ntrunc": ntrunc},
    )


def plot_spectral_power(
    spec: xr.DataArray,
    ax=None,
    logy: bool = True,
    grid: bool = True,
    **plot_kwargs,
):
    """
    Plot spectral power as a function of total wavenumber.

    This function first sums the packed triangular spectral coefficients over
    zonal wavenumber ``m`` with :func:`spectral_power_by_total_wavenumber`, then
    draws :math:`E(n) = \\sum_m |\\hat{a}_{n,m}|^2` against total wavenumber
    ``n``.

    Parameters
    ----------
    spec : xr.DataArray
        One-dimensional packed triangular spectral coefficients. If the input
        data has time, level, or other non-spectral dimensions, select one slice
        before plotting, e.g. ``spec.isel(time=0)``.
    ax : matplotlib.axes.Axes, optional
        Axes used for drawing. If ``None``, a new figure and axes are created.
    logy : bool, default: True
        If ``True``, draw the spectrum with a logarithmic y-axis by calling
        :meth:`matplotlib.axes.Axes.semilogy`. If ``False``, call
        :meth:`matplotlib.axes.Axes.plot`.
    grid : bool, default: True
        Whether to add a light grid to the axes.
    **plot_kwargs
        Additional keyword arguments passed to ``ax.semilogy`` or ``ax.plot``.

    Returns
    -------
    tuple
        ``(fig, ax)``, where ``fig`` is the Matplotlib figure and ``ax`` is the
        axes containing the spectrum.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        visual_spec

    """
    import matplotlib.pyplot as plt

    power = spectral_power_by_total_wavenumber(spec)
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    plot_kwargs.setdefault("linewidth", 1.6)
    if logy:
        ax.semilogy(power["n"].values, power.values, **plot_kwargs)
    else:
        ax.plot(power["n"].values, power.values, **plot_kwargs)
    ax.set_xlabel("Total wavenumber n")
    ax.set_ylabel(r"$\sum_m |\hat{a}_{n,m}|^2$")
    ax.set_title(f"{spec.name or 'Spectral'} power")
    if grid:
        ax.grid(True, alpha=0.3)
    return fig, ax


def plot_spectral_matrix(
    spec: xr.DataArray,
    value: str = "logabs",
    ax=None,
    cmap: str = "viridis",
    add_colorbar: bool = True,
    **imshow_kwargs,
):
    """
    Plot a packed triangular spectral coefficient field as an m-n matrix.

    The packed spectral dimension is reshaped to a triangular matrix whose
    rows are zonal wavenumber ``m`` and columns are total wavenumber ``n``.
    Matrix cells with ``n < m`` are outside the triangular spectral domain and
    are shown as NaN.

    Parameters
    ----------
    spec : xr.DataArray
        One-dimensional packed triangular spectral coefficients. Select any
        non-spectral dimensions before plotting.
    value : {"abs", "logabs", "power", "real", "imag"}, default: "logabs"
        Quantity to visualize:

        - ``"abs"``: coefficient magnitude.
        - ``"logabs"``: base-10 logarithm of coefficient magnitude.
        - ``"power"``: squared coefficient magnitude.
        - ``"real"``: real part of the coefficient.
        - ``"imag"``: imaginary part of the coefficient.
    ax : matplotlib.axes.Axes, optional
        Axes used for drawing. If ``None``, a new axes is created.
    cmap : str, default: "viridis"
        Matplotlib colormap used by :meth:`matplotlib.axes.Axes.imshow`.
    add_colorbar : bool, default: True
        Whether to add a colorbar for the plotted matrix.
    **imshow_kwargs
        Additional keyword arguments passed to ``ax.imshow``. The defaults are
        ``origin="lower"`` and ``aspect="auto"`` unless explicitly provided.

    Returns
    -------
    matplotlib.axes.Axes
        Axes containing the coefficient matrix image.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        visual_spec

    """
    import matplotlib.pyplot as plt

    mat = spectral_coefficients_to_matrix(spec, value=value)
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5))
    imshow_kwargs.setdefault("origin", "lower")
    imshow_kwargs.setdefault("aspect", "auto")
    im = ax.imshow(mat.values, cmap=cmap, **imshow_kwargs)
    if add_colorbar:
        label = {
            "abs": r"$|\hat{a}_{n,m}|$",
            "logabs": r"$\log_{10}(|\hat{a}_{n,m}|)$",
            "power": r"$|\hat{a}_{n,m}|^2$",
            "real": "Real part",
            "imag": "Imaginary part",
        }[value]
        plt.colorbar(im, ax=ax, label=label)
    ax.set_xlabel("Total wavenumber n")
    ax.set_ylabel("Zonal wavenumber m")
    ax.set_title(f"{spec.name or 'Spectral'} coefficients: {value}")
    return ax


def plot_spectral_real_imag(
    spec: xr.DataArray,
    figsize: tuple[float, float] = (11, 4),
    cmap: str = "viridis",
    **imshow_kwargs,
):
    """
    Plot real and imaginary parts of packed triangular spectral coefficients.

    This is a convenience wrapper around :func:`plot_spectral_matrix`. It
    creates a two-panel figure and draws ``value="real"`` in the left panel and
    ``value="imag"`` in the right panel.

    Parameters
    ----------
    spec : xr.DataArray
        One-dimensional packed triangular spectral coefficients. Select any
        non-spectral dimensions before plotting.
    figsize : tuple of float, default: (11, 4)
        Matplotlib figure size passed to :func:`matplotlib.pyplot.subplots`.
    cmap : str, default: "viridis"
        Matplotlib colormap used for both matrix panels.
    **imshow_kwargs
        Additional keyword arguments passed to :func:`plot_spectral_matrix` and
        then to ``ax.imshow`` for both panels.

    Returns
    -------
    tuple
        ``(fig, ax)``, where ``fig`` is the Matplotlib figure and ``ax`` is the
        two-element array of axes containing the real and imaginary matrices.

    .. ecl-minigallery::
        :add-heading: Example(s) related to the function

        visual_spec

    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=figsize)
    plot_spectral_matrix(spec, value="real", ax=ax[0], cmap=cmap, **imshow_kwargs)
    plot_spectral_matrix(spec, value="imag", ax=ax[1], cmap=cmap, **imshow_kwargs)
    fig.tight_layout()
    return fig, ax


def _triangular_mode_indices(ntrunc: int) -> tuple[np.ndarray, np.ndarray]:
    """Return packed-mode arrays m, n for ordering: for m, for n=m..ntrunc."""
    ms: list[int] = []
    ns: list[int] = []
    for m in range(ntrunc + 1):
        for n in range(m, ntrunc + 1):
            ms.append(m)
            ns.append(n)
    return np.asarray(ms, dtype=np.int64), np.asarray(ns, dtype=np.int64)


def _infer_triangular_ntrunc(nmode: int) -> int:
    """Infer triangular truncation ntrunc from packed spectral mode length."""
    disc = 1 + 8 * int(nmode)
    root = int(round(np.sqrt(disc)))
    if root * root != disc:
        raise ValueError(
            f"The mode length {nmode} is incompatible with the triangular spectral truncation. (ntrunc+1)(ntrunc+2)/2"
        )
    ntrunc = (root - 3) // 2
    if (ntrunc + 1) * (ntrunc + 2) // 2 != nmode:
        raise ValueError(
            f"The mode length {nmode} is incompatible with the triangular spectral truncation."
        )
    return int(ntrunc)
