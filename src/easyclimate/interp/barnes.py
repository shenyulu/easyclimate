"""
Fast Barnes interpolation
"""

from __future__ import annotations

from typing import Literal, Optional
import numpy as np
import xarray as xr
import pandas as pd

from math import sqrt
from ..backend import (
    barnes_rs,
    barnesS2_rs,
    radius_mask_2d_rs,
    barnes_numba,
    barnes_S2_numba,
    kdtree_fastbarnes as kdtree,
)

__all__ = [
    "interp_spatial_barnes",
    "interp_spatial_barnes_rs",
    "interp_spatial_barnesS2",
    "interp_spatial_barnesS2_rs",
]


# -------------------------
# Helpers
# -------------------------
def _build_lonlat_grid(point, grid_x, grid_y, step_deg):
    """
    Build regular lon/lat grid.

    point: [lon0, lat0] lower-left corner (degrees)
    grid_x, grid_y: domain size (degrees)
    step_deg: grid spacing (degrees)
    """
    lon0, lat0 = float(point[0]), float(point[1])
    nx = int(round(grid_x / step_deg)) + 1
    ny = int(round(grid_y / step_deg)) + 1
    x0 = np.asarray([lon0, lat0], dtype=np.float64)
    size = np.asarray([nx, ny], dtype=np.int64)

    gridX = lon0 + np.arange(nx) * step_deg
    gridY = lat0 + np.arange(ny) * step_deg
    return x0, size, gridX, gridY


def _precheck_lambert_options_S2(
    *,
    method: str,
    x0: np.ndarray,
    step_deg: float,
    size: tuple[int, int] | np.ndarray,
    num_iter: int,
    sigma_grid: float,
    lambert_proj,
    lambert_grid,
    auto_proj: bool,
):
    """
    Front-load Lambert argument checks so errors happen in Python (readable),
    not deep inside rust/numba.

    Mirrors key checks in fastbarnes_new.interpolationS2.interpolate_opt_convol_S2_part1
    and _infer_lambert_proj (incl. equator-crossing restriction).
    """
    # Only matters for optimized_convolution_S2
    if str(method) != "optimized_convolution_S2":
        return

    # Normalize size
    if isinstance(size, np.ndarray):
        size = tuple(int(x) for x in size.tolist())
    else:
        size = (int(size[0]), int(size[1]))

    if not isinstance(step_deg, (int, float)) or float(step_deg) <= 0:
        raise ValueError("grid_res_deg/step must be > 0")

    if int(num_iter) < 1:
        raise ValueError("num_iter must be >= 1")

    # --- auto_proj=False needs both lambert_proj & lambert_grid
    if not auto_proj:
        if lambert_proj is None:
            raise ValueError("lambert_proj must be specified when auto_proj is False")
        if lambert_grid is None:
            raise ValueError("lambert_grid must be specified when auto_proj is False")

    # --- lambert_proj validation (if provided)
    if lambert_proj is not None:
        try:
            lambert_proj = tuple(lambert_proj)
        except Exception as e:
            raise ValueError(
                f"lambert_proj must be tuple-like of length 5, got {type(lambert_proj)}"
            ) from e
        if len(lambert_proj) != 5:
            raise ValueError(
                f"lambert_proj must have length 5, got {len(lambert_proj)}"
            )

    # --- lambert_grid validation (if provided)
    if lambert_grid is not None:
        if not (isinstance(lambert_grid, (tuple, list)) and len(lambert_grid) == 2):
            raise ValueError("lambert_grid must be (lam_x0, lam_size)")
        lam_x0, lam_size = lambert_grid

        lam_x0 = np.asarray(lam_x0, dtype=np.float64)
        if lam_x0.shape != (2,):
            raise ValueError("lambert_grid[0] (lam_x0) must be length-2 array-like")

        if not (isinstance(lam_size, (tuple, list)) and len(lam_size) == 2):
            raise ValueError("lambert_grid[1] (lam_size) must be length-2 tuple/list")
        lam_size = (int(lam_size[0]), int(lam_size[1]))
        if lam_size[0] < 2 or lam_size[1] < 2:
            raise ValueError("lambert_grid size must be >= (2,2)")

    # --- equator-crossing restriction (only when we need to infer projection)
    # interpolationS2._infer_lambert_proj forbids domains crossing equator
    # (optimized_convolution_S2 does not support it) :contentReference[oaicite:3]{index=3}
    if auto_proj and (lambert_proj is None):
        lon0 = float(x0[0])
        lat0 = float(x0[1])
        lat1 = float(lat0 + (size[1] - 1) * float(step_deg))
        lat_min = min(lat0, lat1)
        lat_max = max(lat0, lat1)
        if lat_min < 0.0 < lat_max:
            raise RuntimeError(
                "optimized_convolution_S2 does not support domains crossing the equator; "
                "split into hemispheres or pass a custom lambert_proj/lambert_grid"
            )


# -------------------------
# API
# -------------------------
def interp_spatial_barnes(
    data: pd.DataFrame,
    var_name: str,
    *,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    # grid (degrees)
    grid_res_deg: float = 0.25,
    auto_domain: bool = True,
    buffer_deg: float = 1.0,
    point: tuple[float, float] | None = None,
    grid_x: float | None = None,
    grid_y: float | None = None,
    # Barnes (degrees)
    sigma_deg: float | None = None,  # None -> auto
    influence_radius_deg: float | None = None,  # None -> auto (=k*sigma)
    influence_k_sigma: float = 4.0,
    num_iter: int = 2,
    method: str = "optimized_convolution",
    max_dist_sigma: float | None = None,  # optional override
    min_weight: float = 1e-6,
    # mask
    mask_radius_deg: float | None = None,  # None -> auto (=influence_radius_deg)
    # data cleaning
    missing_sentinels: tuple[float, ...] = (9999.9, 9999.0, -9999.0, -9999.9),
):
    """
    Computes Barnes interpolation for observation values ``var_name`` sampled at irregular
    locations in ``data`` (lon/lat in degrees), using Gaussian weights with width
    parameter ``sigma_deg`` and returning a regular lon/lat grid as an
    :py:class:`xarray.DataArray<xarray.DataArray>`.

    Barnes interpolation is widely used in meteorology and geosciences to remodel
    irregular point observations into a smooth gridded field. It can be written as

    .. math::
        f(\\boldsymbol{x})=\\frac{\\sum_{k=1}^N f_k\\cdot w_k(\\boldsymbol{x})}{\\sum_{k=1}^N w_k(\\boldsymbol{x})}

    with Gaussian weights

    .. math::
        w_k(\\boldsymbol{x})=\\text{e}^{-\\frac{1}{2\\sigma^2}\\left|x-\\boldsymbol{x}_k\\right|^2}

    Naive computation of Barnes interpolation leads to an algorithmic complexity of :math:`O(N \\times W \\times H)`,
    where :math:`N` is the number of sample points and :math:`W \\times H` the size of the underlying grid.

    For sufficiently large :math:`n` (in general in the range from 3 to 6) a good approximation of
    Barnes interpolation with a reduced complexity :math:`O(N + W \\times H)` can be obtained by the convolutional expression

    .. math::
        f(\\boldsymbol{x})\\approx \\frac{ (\\sum_{k=1}^{N}f_k\\cdot\\delta_{\\boldsymbol{x}_k}) *  ( r_n^{*n[x]}(x)\\cdot r_n^{*n[y]}(y) )   }{ ( \\sum_{k=1}^{N} \\delta_{\\boldsymbol{x}_k}  ) *  (  r_{n}^{*n[x]}(x)\\cdot r_{n}^{*n[y]}(y)  )   }

    where :math:`\\delta` is the Dirac impulse function and :math:`r(.)` an elementary rectangular function of a specific length that depends on :math:`\\sigma` and :math:`n`.

    This wrapper is **degree-based**:

    - ``grid_res_deg`` / ``sigma_deg`` / ``influence_radius_deg`` / ``mask_radius_deg`` are all in degrees.
    - Internally converts to the fast-barnes "grid units" used by the core implementation:

    .. math::
        \\sigma_{\\text{grid}} = \\sigma_{\deg} / \\Delta_{\\deg},
        \\qquad
        \\text{max_dist_sigma} = R_{\\deg} / \\sigma_{\\deg}.

    - If ``sigma_deg`` is not provided, it is estimated from station density using a simple spacing heuristic.
    - Optional masking keeps only grid points that have at least one station within ``mask_radius_deg`` using a radius-search KD-tree.

    Parameters
    ----------
    data : :py:class:`pandas.DataFrame<pandas.DataFrame>`
        Input table containing at least ``lon_dim``, ``lat_dim`` and ``var_name`` columns.
        Any rows with NaN/Inf in these columns are dropped after cleaning.

        There should be a similar structure as follows

        +------------+------------+-----------+
        |    lon     |    lat     |    qff    |
        +============+============+===========+
        |   -3.73    |   56.33    |   995.1   |
        +------------+------------+-----------+
        |    2.64    |   47.05    |  1012.5   |
        +------------+------------+-----------+
        |    ...     |   ...      |   ...     |
        +------------+------------+-----------+

        .. note::
            Data points should contain longitude (`lon`), latitude (`lat`) and data variables (the above data variable name is `qff`).

    var_name : :py:class:`str<str>`
        Name of the variable column to interpolate. This should match the one in the parameter `data`.
    lon_dim, lat_dim : :py:class:`str<str>`, optional
        Column names for longitude/latitude (degrees). Defaults are ``"lon"`` and ``"lat"``.
    grid_res_deg : :py:class:`float<float>`, optional
        Output grid spacing in degrees. Must be > 0. Default is 0.25.
    auto_domain : :py:class:`bool<bool>`, optional
        If True (default), the interpolation domain is set to the data bounding box
        expanded by ``buffer_deg`` on each side.
    buffer_deg : :py:class:`float<float>`, optional
        Padding (degrees) added around the data bounding box when ``auto_domain=True``.
    point : tuple(float, float), optional
        Lower-left corner (lon0, lat0) of the output grid (degrees). Required when
        ``auto_domain=False``.
    grid_x, grid_y : :py:class:`float<float>`, optional
        Domain size in degrees in x (lon) / y (lat). Required when ``auto_domain=False``.
    sigma_deg : :py:class:`float<float>`, optional
        Gaussian width in degrees. If None, an automatic estimate based on station density is used.
    influence_radius_deg : :py:class:`float<float>`, optional
        Radius of influence in degrees. If None, uses ``influence_k_sigma * sigma_deg``.
    influence_k_sigma : :py:class:`float<float>`, optional
        Multiplier used when ``influence_radius_deg`` is None. Default is 4.0.
    num_iter : :py:class:`int<int>`, optional
        Number of self-convolutions used by convolution-based methods. Must be >= 1.
        The number of performed self-convolutions of the underlying rect-kernel.
        Applies only if method is 'optimized_convolution' or 'convolution'.
        The default is 2. Applies only to Convol interpolations: one of 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50.
    method : :py:class:`str<str>`, optional
        Passed to :func:`fastbarnes_new.interpolation.barnes`. Common values include
        ``"optimized_convolution"``, ``"convolution"``, ``"radius"``, ``"naive"``.
        The possible implementations that can be chosen are 'naive' for the straightforward
        implementation (algorithm A from paper), 'radius' to consider only sample
        points within a specific radius of influence, both with an algorithmic
        complexity of :math:`O(N \\times W \\times H)`.
        The choice 'convolution' implements algorithm B specified in the paper
        and 'optimized_convolution' is its optimization by appending tail values
        to the rectangular kernel. The latter two algorithms reduce the complexity
        down to :math:`O(N + W \\times H)`.
    max_dist_sigma : :py:class:`float<float>`, optional
        Maximum distance (in units of ``sigma``) for which interpolation is computed.
        If None, uses ``influence_radius_deg / sigma_deg``.
    min_weight : :py:class:`float<float>`, optional
        Minimum Gaussian weight threshold used by radius-based methods in the backend.
    mask_radius_deg : :py:class:`float<float>`, optional
        If provided (or if None defaults to ``influence_radius_deg``), grid points farther
        than this radius (degrees) from all stations are set to NaN.
    missing_sentinels : tuple(float, ...), optional
        Values treated as missing in ``var_name`` and replaced by NaN before cleaning.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`
        Interpolated field on a regular lat/lon grid. Dimensions are ``(lat_dim, lon_dim)``.
        The returned DataArray includes useful metadata in ``attrs`` (grid spacing, sigma,
        influence radius, domain definition, etc.).

    .. seealso::
        - https://github.com/MeteoSwiss/fast-barnes-py
        - Zürcher, B. K.: Fast approximate Barnes interpolation: illustrated by Python-Numba implementation fast-barnes-py v1.0, Geosci. Model Dev., 16, 1697–1711, https://doi.org/10.5194/gmd-16-1697-2023, 2023.

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_interp.py
    """
    if grid_res_deg <= 0:
        raise ValueError("grid_res_deg must be > 0")
    if num_iter < 1:
        raise ValueError("num_iter must be >= 1")

    # -------------------------
    # Clean data
    # -------------------------
    df = data[[lon_dim, lat_dim, var_name]].copy()

    for s in missing_sentinels:
        df.loc[np.isclose(df[var_name].astype(float), s), var_name] = np.nan

    df = df.replace([np.inf, -np.inf], np.nan).dropna(
        subset=[lon_dim, lat_dim, var_name]
    )
    if len(df) < 3:
        raise ValueError(f"Not enough valid points after cleaning: {len(df)}")

    stn_lon = df[lon_dim].astype(float).to_numpy()
    stn_lat = df[lat_dim].astype(float).to_numpy()
    values = df[var_name].astype(float).to_numpy()

    # -------------------------
    # Auto domain from bbox
    # -------------------------
    if auto_domain:
        lon_min, lon_max = float(np.min(stn_lon)), float(np.max(stn_lon))
        lat_min, lat_max = float(np.min(stn_lat)), float(np.max(stn_lat))
        lon0 = lon_min - buffer_deg
        lat0 = lat_min - buffer_deg
        gx = (lon_max - lon_min) + 2 * buffer_deg
        gy = (lat_max - lat_min) + 2 * buffer_deg
    else:
        if point is None or grid_x is None or grid_y is None:
            raise ValueError(
                "When auto_domain=False, you must provide point, grid_x, grid_y"
            )
        lon0, lat0 = float(point[0]), float(point[1])
        gx, gy = float(grid_x), float(grid_y)

    step = float(grid_res_deg)
    nx = int(round(gx / step)) + 1
    ny = int(round(gy / step)) + 1
    x0 = np.asarray([lon0, lat0], dtype=np.float64)
    size = np.asarray([nx, ny], dtype=np.int64)

    gridX = lon0 + np.arange(nx) * step
    gridY = lat0 + np.arange(ny) * step

    # -------------------------
    # Auto sigma (deg) WITHOUT kNN:
    # estimate mean spacing from density in bbox.
    # spacing ~ sqrt(area / N)
    # -------------------------
    if sigma_deg is None:
        # domain area in degree^2 (rough, but sufficient for parameter scale)
        area = max(gx * gy, 1e-6)
        spacing = sqrt(area / len(df))
        # sigma as a fraction of spacing; 0.5~0.8 often reasonable
        sigma_deg = 0.6 * spacing
        # clamp to a reasonable range to avoid insane values
        sigma_deg = float(np.clip(sigma_deg, 0.2, 3.0))
    else:
        sigma_deg = float(sigma_deg)
        if sigma_deg <= 0:
            raise ValueError("sigma_deg must be > 0")

    # -------------------------
    # Influence radius (deg) and mask radius (deg)
    # -------------------------
    if influence_radius_deg is None:
        influence_radius_deg = influence_k_sigma * sigma_deg
    influence_radius_deg = float(influence_radius_deg)

    if max_dist_sigma is None:
        max_dist_sigma = influence_radius_deg / sigma_deg
    else:
        max_dist_sigma = float(max_dist_sigma)

    if mask_radius_deg is None:
        mask_radius_deg = influence_radius_deg
    mask_radius_deg = float(mask_radius_deg) if mask_radius_deg is not None else None

    # Convert to fast-barnes grid-units
    sigma_grid = sigma_deg / step  # degrees -> grid units

    # -------------------------
    # Run Barnes
    # -------------------------
    pts = np.c_[stn_lon, stn_lat].astype(np.float64)

    field = barnes_numba(
        pts,
        values.astype(np.float64),
        float(sigma_grid),
        x0,
        float(step),
        size,
        method=method,
        num_iter=int(num_iter),
        max_dist=float(max_dist_sigma),
        min_weight=float(min_weight),
    )

    da = xr.DataArray(
        field,
        dims=(lat_dim, lon_dim),
        coords={lon_dim: gridX, lat_dim: gridY},
        name=var_name,
        attrs={
            "grid_res_deg": step,
            "sigma_deg": sigma_deg,
            "sigma_grid": float(sigma_grid),
            "influence_radius_deg": influence_radius_deg,
            "max_dist_sigma": float(max_dist_sigma),
            "method": method,
            "num_iter": int(num_iter),
            "auto_domain": bool(auto_domain),
            "buffer_deg": float(buffer_deg),
            "domain_lon0": float(lon0),
            "domain_lat0": float(lat0),
            "domain_grid_x_deg": float(gx),
            "domain_grid_y_deg": float(gy),
        },
    )

    # -------------------------
    # Mask using fastbarnes kdtree radius-search
    # We *don't* compute nearest distance; we simply keep grid points
    # that have >=1 station within mask_radius_deg.
    # -------------------------
    if mask_radius_deg is not None:
        # build kd-tree once
        tree, pts_ref = kdtree.create_kdtree(pts)
        kd_search = kdtree.prepare_search(float(mask_radius_deg), tree, pts_ref)
        res_index, res_sqr_dist, _, _, _ = kd_search  # reusable arrays

        mask = np.zeros((ny, nx), dtype=bool)
        c = np.empty(2, dtype=np.float64)
        for j in range(ny):
            c[1] = lat0 + j * step
            for i in range(nx):
                c[0] = lon0 + i * step
                kdtree.radius_search(c, *kd_search)
                mask[j, i] = res_index[-1] > 0

        da = da.where(mask)
        da.attrs["mask_radius_deg"] = float(mask_radius_deg)

    return da


def interp_spatial_barnes_rs(
    data: pd.DataFrame,
    var_name: str,
    *,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    # grid (degrees)
    grid_res_deg: float = 0.25,
    auto_domain: bool = True,
    buffer_deg: float = 1.0,
    point: tuple[float, float] | None = None,
    grid_x: float | None = None,
    grid_y: float | None = None,
    # Barnes (degrees)
    sigma_deg: float | None = None,  # None -> auto
    influence_radius_deg: float | None = None,  # None -> auto (=k*sigma)
    influence_k_sigma: float = 4.0,
    num_iter: int = 2,
    method: str = "optimized_convolution",
    max_dist_sigma: float | None = None,  # optional override
    min_weight: float = 1e-6,
    # mask
    mask_radius_deg: float | None = None,  # None -> auto (=influence_radius_deg)
    # data cleaning
    missing_sentinels: tuple[float, ...] = (9999.9, 9999.0, -9999.0, -9999.9),
):
    """
    Rust-accelerated Barnes interpolation wrapper (lon/lat in degrees).

    This function is API-compatible with :func:`interp_spatial_barnes` but calls the
    Rust backend (:func:`easyclimate_rust._easyclimate_rust.barnes`) for the core
    interpolation step, and (optionally) uses the Rust radius mask
    (:func:`easyclimate_rust._easyclimate_rust.radius_mask_2d`) to avoid showing
    extrapolation far from stations.

    Mathematically, the method targets the same Barnes analysis:

    .. math::
        f(\\boldsymbol{x})=\\frac{\\sum_{k=1}^N f_k\\cdot w_k(\\boldsymbol{x})}{\\sum_{k=1}^N w_k(\\boldsymbol{x})}

    with Gaussian weights

    .. math::
        w_k(\\boldsymbol{x})=\\text{e}^{-\\frac{1}{2\\sigma^2}\\left|x-\\boldsymbol{x}_k\\right|^2}

    .. note::
        - All user-facing spatial parameters are in **degrees** (same conversion rules as :func:`interp_spatial_barnes`).
        - ``sigma_deg`` auto-estimation uses a station-density spacing heuristic.
        - If ``mask_radius_deg`` is enabled, masking is performed in Rust on the regular grid.

    Parameters
    ----------
    Identical to :func:`interp_spatial_barnes`.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`
        Interpolated field on a regular lat/lon grid. Dimensions are ``(lat_dim, lon_dim)``.
        The returned DataArray includes useful metadata in ``attrs`` (grid spacing, sigma,
        influence radius, domain definition, etc.).

    .. seealso::
        - :func:`interp_spatial_barnes` (pure Python/NumPy backend)
        - https://github.com/MeteoSwiss/fast-barnes-py
        - Zürcher, B. K.: Fast approximate Barnes interpolation: illustrated by Python-Numba implementation fast-barnes-py v1.0, Geosci. Model Dev., 16, 1697–1711, https://doi.org/10.5194/gmd-16-1697-2023, 2023.

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_interp.py
    """
    if grid_res_deg <= 0:
        raise ValueError("grid_res_deg must be > 0")
    if num_iter < 1:
        raise ValueError("num_iter must be >= 1")

    # -------------------------
    # Clean data
    # -------------------------
    df = data[[lon_dim, lat_dim, var_name]].copy()

    for s in missing_sentinels:
        df.loc[np.isclose(df[var_name].astype(float), s), var_name] = np.nan

    df = df.replace([np.inf, -np.inf], np.nan).dropna(
        subset=[lon_dim, lat_dim, var_name]
    )
    if len(df) < 3:
        raise ValueError(f"Not enough valid points after cleaning: {len(df)}")

    stn_lon = df[lon_dim].astype(float).to_numpy()
    stn_lat = df[lat_dim].astype(float).to_numpy()
    values = df[var_name].astype(float).to_numpy()

    # -------------------------
    # Auto domain from bbox
    # -------------------------
    if auto_domain:
        lon_min, lon_max = float(np.min(stn_lon)), float(np.max(stn_lon))
        lat_min, lat_max = float(np.min(stn_lat)), float(np.max(stn_lat))
        lon0 = lon_min - buffer_deg
        lat0 = lat_min - buffer_deg
        gx = (lon_max - lon_min) + 2 * buffer_deg
        gy = (lat_max - lat_min) + 2 * buffer_deg
    else:
        if point is None or grid_x is None or grid_y is None:
            raise ValueError(
                "When auto_domain=False, you must provide point, grid_x, grid_y"
            )
        lon0, lat0 = float(point[0]), float(point[1])
        gx, gy = float(grid_x), float(grid_y)

    step = float(grid_res_deg)
    nx = int(round(gx / step)) + 1
    ny = int(round(gy / step)) + 1
    x0 = np.asarray([lon0, lat0], dtype=np.float64)
    size = np.asarray([nx, ny], dtype=np.int64)

    gridX = lon0 + np.arange(nx) * step
    gridY = lat0 + np.arange(ny) * step

    # -------------------------
    # Auto sigma (deg) WITHOUT kNN:
    # estimate mean spacing from density in bbox.
    # spacing ~ sqrt(area / N)
    # -------------------------
    if sigma_deg is None:
        # domain area in degree^2 (rough, but sufficient for parameter scale)
        area = max(gx * gy, 1e-6)
        spacing = sqrt(area / len(df))
        # sigma as a fraction of spacing; 0.5~0.8 often reasonable
        sigma_deg = 0.6 * spacing
        # clamp to a reasonable range to avoid insane values
        sigma_deg = float(np.clip(sigma_deg, 0.2, 3.0))
    else:
        sigma_deg = float(sigma_deg)
        if sigma_deg <= 0:
            raise ValueError("sigma_deg must be > 0")

    # -------------------------
    # Influence radius (deg) and mask radius (deg)
    # -------------------------
    if influence_radius_deg is None:
        influence_radius_deg = influence_k_sigma * sigma_deg
    influence_radius_deg = float(influence_radius_deg)

    if max_dist_sigma is None:
        max_dist_sigma = influence_radius_deg / sigma_deg
    else:
        max_dist_sigma = float(max_dist_sigma)

    if mask_radius_deg is None:
        mask_radius_deg = influence_radius_deg
    mask_radius_deg = float(mask_radius_deg) if mask_radius_deg is not None else None

    # Convert to fast-barnes grid-units
    sigma_grid = sigma_deg / step  # degrees -> grid units

    # -------------------------
    # Run Barnes
    # -------------------------
    pts = np.c_[stn_lon, stn_lat].astype(np.float64)

    field = barnes_rs(
        pts,
        values.astype(np.float64),
        float(sigma_grid),
        x0,
        float(step),
        size,
        method=method,
        num_iter=int(num_iter),
        max_dist=float(max_dist_sigma),
        min_weight=float(min_weight),
    )

    da = xr.DataArray(
        field,
        dims=(lat_dim, lon_dim),
        coords={lon_dim: gridX, lat_dim: gridY},
        name=var_name,
        attrs={
            "grid_res_deg": step,
            "sigma_deg": sigma_deg,
            "sigma_grid": float(sigma_grid),
            "influence_radius_deg": influence_radius_deg,
            "max_dist_sigma": float(max_dist_sigma),
            "method": method,
            "num_iter": int(num_iter),
            "auto_domain": bool(auto_domain),
            "buffer_deg": float(buffer_deg),
            "domain_lon0": float(lon0),
            "domain_lat0": float(lat0),
            "domain_grid_x_deg": float(gx),
            "domain_grid_y_deg": float(gy),
        },
    )

    # -------------------------
    # Mask using fastbarnes kdtree radius-search
    # We *don't* compute nearest distance; we simply keep grid points
    # that have >=1 station within mask_radius_deg.
    # -------------------------
    if mask_radius_deg is not None:
        # pts is (N,2) float64, x0 is (2,), step is scalar, size is (2,) [nx, ny]
        mask = radius_mask_2d_rs(
            pts.astype(np.float64),
            x0.astype(np.float64),
            np.asarray([step, step], dtype=np.float64),
            np.asarray([nx, ny], dtype=np.int64),
            float(mask_radius_deg),
        )
        da = da.where(mask)
        da.attrs["mask_radius_deg"] = float(mask_radius_deg)

    return da


# -------------------------
# S2 versions
# -------------------------
def interp_spatial_barnesS2(
    data: pd.DataFrame,
    var_name: str,
    *,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    # grid (degrees)
    grid_res_deg: float = 0.25,
    auto_domain: bool = True,
    buffer_deg: float = 1.0,
    point: tuple[float, float] | None = None,
    grid_x: float | None = None,
    grid_y: float | None = None,
    # Barnes (degrees)
    sigma_deg: float | None = None,  # None -> auto
    influence_radius_deg: float | None = None,  # None -> auto (=k*sigma)
    influence_k_sigma: float = 4.0,
    num_iter: int = 4,
    method: Literal[
        "optimized_convolution_S2", "naive_S2"
    ] = "optimized_convolution_S2",
    max_dist_sigma: float | None = None,  # optional override (in sigma units)
    resample: bool = True,
    # mask
    mask_radius_deg: float | None = None,  # None -> auto (=influence_radius_deg)
    # data cleaning
    missing_sentinels: tuple[float, ...] = (9999.9, 9999.0, -9999.0, -9999.9),
    # Lambert options (passed through to barnes_S2)
    lambert_proj=None,
    lambert_grid=None,
    auto_proj: bool = True,
) -> xr.DataArray:
    """
    Computes Barnes interpolation on the unit sphere :math:`S^2` (spherical metric)
    for observations given in lon/lat degrees, returning a regular lon/lat grid as an
    :py:class:`xarray.DataArray<xarray.DataArray>`.

    Compared to planar Barnes interpolation, the :math:`S^2` variant uses spherical
    geometry internally (implemented by :func:`fastbarnes_new.interpolationS2.barnes_S2`).
    For performance, the optimized method relies on a Lambert projection and a projected
    grid; this wrapper forwards Lambert options and performs early validation so errors
    are raised in Python (more readable than deep backend errors).


    Parameters
    ----------
    data : :py:class:`pandas.DataFrame<pandas.DataFrame>`
        Input table containing at least ``lon_dim``, ``lat_dim`` and ``var_name`` columns.
        Any rows with NaN/Inf in these columns are dropped after cleaning.

        There should be a similar structure as follows

        +------------+------------+-----------+
        |    lon     |    lat     |    qff    |
        +============+============+===========+
        |   -3.73    |   56.33    |   995.1   |
        +------------+------------+-----------+
        |    2.64    |   47.05    |  1012.5   |
        +------------+------------+-----------+
        |    ...     |   ...      |   ...     |
        +------------+------------+-----------+

        .. note::
            Data points should contain longitude (`lon`), latitude (`lat`) and data variables (the above data variable name is `qff`).

    var_name : :py:class:`str<str>`
        Name of the variable column to interpolate. This should match the one in the parameter `data`.
    lon_dim, lat_dim : :py:class:`str<str>`, optional
        Column names for longitude/latitude (degrees). Defaults are ``"lon"`` and ``"lat"``.
    grid_res_deg : :py:class:`float<float>`, optional
        Output grid spacing in degrees. Must be > 0. Default is 0.25.
    auto_domain : :py:class:`bool<bool>`, optional
        If True (default), the interpolation domain is set to the data bounding box
        expanded by ``buffer_deg`` on each side.
    buffer_deg : :py:class:`float<float>`, optional
        Padding (degrees) added around the data bounding box when ``auto_domain=True``.
    point : tuple(float, float), optional
        Lower-left corner (lon0, lat0) of the output grid (degrees). Required when
        ``auto_domain=False``.
    grid_x, grid_y : :py:class:`float<float>`, optional
        Domain size in degrees in x (lon) / y (lat). Required when ``auto_domain=False``.
    sigma_deg : :py:class:`float<float>`, optional
        Gaussian width in degrees. If None, an automatic estimate based on station density is used.
    influence_radius_deg : :py:class:`float<float>`, optional
        Radius of influence in degrees. If None, uses ``influence_k_sigma * sigma_deg``.
    influence_k_sigma : :py:class:`float<float>`, optional
        Multiplier used when ``influence_radius_deg`` is None. Default is 4.0.
    num_iter : :py:class:`int<int>`, optional
        Number of self-convolutions used by convolution-based methods. Must be >= 1.
        The number of performed self-convolutions of the underlying rect-kernel.
        Applies only if method is 'optimized_convolution' or 'convolution'.
        The default is 2. Applies only to Convol interpolations: one of 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50.
    method : {'optimized_convolution_S2', 'naive_S2'}, default: 'optimized_convolution_S2'.
        Designates the Barnes interpolation method to be used. The possible
        implementations that can be chosen are 'naive_S2' for the straightforward
        implementation (algorithm A from the paper) with an algorithmic complexity
        of :math:`O(N \\times W \\times H)`.
        The choice 'optimized_convolution_S2' implements the optimized algorithm B
        specified in the paper by appending tail values to the rectangular kernel.
        The latter algorithm has a reduced complexity of :math:`O(N + W \\times H)`.
        The default is 'optimized_convolution_S2'.
    max_dist_sigma : :py:class:`float<float>`, optional
        Maximum distance (in units of ``sigma``) for which interpolation is computed.
        If None, uses ``influence_radius_deg / sigma_deg``.
    min_weight : :py:class:`float<float>`, optional
        Minimum Gaussian weight threshold used by radius-based methods in the backend.
    mask_radius_deg : :py:class:`float<float>`, optional
        If provided (or if None defaults to ``influence_radius_deg``), grid points farther
        than this radius (degrees) from all stations are set to NaN.
    missing_sentinels : tuple(float, ...), optional
        Values treated as missing in ``var_name`` and replaced by NaN before cleaning.
    resample : :py:class:`bool<bool>`, default True
        Passed through to the S2 backend (controls internal resampling behavior).
    lambert_proj, lambert_grid, auto_proj
        Lambert projection configuration forwarded to :func:`barnes_S2`.
        If ``auto_proj=True`` and no ``lambert_proj`` is provided, the backend will infer
        a projection; **domains crossing the equator are not supported** by the optimized
        S2 convolution path (split into hemispheres or pass an explicit projection/grid).

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`
        Interpolated spherical Barnes field on a regular lat/lon grid, with metadata in ``attrs``.

    .. seealso::
        - :func:`interp_spatial_barnes` (planar metric)
        - :func:`interp_spatial_barnesS2_rs` (Rust backend on :math:`S^2`)
        - https://github.com/MeteoSwiss/fast-barnes-py
        - Zürcher, B. K.: Fast approximate Barnes interpolation: illustrated by Python-Numba implementation fast-barnes-py v1.0, Geosci. Model Dev., 16, 1697–1711, https://doi.org/10.5194/gmd-16-1697-2023, 2023.
    """

    if grid_res_deg <= 0:
        raise ValueError("grid_res_deg must be > 0")
    if num_iter < 1:
        raise ValueError("num_iter must be >= 1")

    # -------------------------
    # Clean data (match interp_spatial_barnes)
    # -------------------------
    df = data[[lon_dim, lat_dim, var_name]].copy()

    for s in missing_sentinels:
        df.loc[np.isclose(df[var_name].astype(float), s), var_name] = np.nan

    df = df.replace([np.inf, -np.inf], np.nan).dropna(
        subset=[lon_dim, lat_dim, var_name]
    )
    if len(df) < 3:
        raise ValueError(f"Not enough valid points after cleaning: {len(df)}")

    stn_lon = df[lon_dim].astype(float).to_numpy()
    stn_lat = df[lat_dim].astype(float).to_numpy()
    values = df[var_name].astype(float).to_numpy()

    # -------------------------
    # Auto domain from bbox (match interp_spatial_barnes)
    # -------------------------
    if auto_domain:
        lon_min, lon_max = float(np.min(stn_lon)), float(np.max(stn_lon))
        lat_min, lat_max = float(np.min(stn_lat)), float(np.max(stn_lat))
        lon0 = lon_min - buffer_deg
        lat0 = lat_min - buffer_deg
        gx = (lon_max - lon_min) + 2 * buffer_deg
        gy = (lat_max - lat_min) + 2 * buffer_deg
        point = (lon0, lat0)
        grid_x = gx
        grid_y = gy
    else:
        if point is None or grid_x is None or grid_y is None:
            raise ValueError(
                "When auto_domain=False, you must provide point, grid_x, grid_y"
            )
        lon0, lat0 = float(point[0]), float(point[1])
        gx, gy = float(grid_x), float(grid_y)

    step = float(grid_res_deg)
    x0, size, gridX, gridY = _build_lonlat_grid(
        point, float(grid_x), float(grid_y), step
    )
    nx = int(size[0])
    ny = int(size[1])

    # -------------------------
    # Auto sigma (deg) WITHOUT kNN (same heuristic as interp_spatial_barnes)
    # -------------------------
    if sigma_deg is None:
        area = max(float(grid_x) * float(grid_y), 1e-6)  # deg^2
        spacing = sqrt(area / len(df))
        sigma_deg = float(np.clip(0.6 * spacing, 0.2, 3.0))
    else:
        sigma_deg = float(sigma_deg)
        if sigma_deg <= 0:
            raise ValueError("sigma_deg must be > 0")

    # -------------------------
    # Influence radius + max_dist in sigma units
    # -------------------------
    if influence_radius_deg is None:
        influence_radius_deg = influence_k_sigma * sigma_deg
    influence_radius_deg = float(influence_radius_deg)

    if max_dist_sigma is None:
        max_dist = influence_radius_deg / sigma_deg
    else:
        max_dist = float(max_dist_sigma)

    # mask radius default = influence radius
    if mask_radius_deg is None:
        mask_radius_deg = influence_radius_deg
    mask_radius_deg = float(mask_radius_deg) if mask_radius_deg is not None else None

    # Convert to grid units for barnes_S2
    sigma_grid = float(sigma_deg) / step

    _precheck_lambert_options_S2(
        method=str(method),
        x0=x0,
        step_deg=float(step),
        size=tuple(size.tolist()) if isinstance(size, np.ndarray) else tuple(size),
        num_iter=int(num_iter),
        sigma_grid=float(sigma_grid),
        lambert_proj=lambert_proj,
        lambert_grid=lambert_grid,
        auto_proj=bool(auto_proj),
    )

    # -------------------------
    # Run Barnes_S2
    # -------------------------
    pts = np.c_[stn_lon, stn_lat].astype(np.float64)

    field = barnes_S2_numba(
        pts,
        values.astype(np.float64),
        float(sigma_grid),
        x0.astype(np.float64),
        float(step),
        tuple(size.tolist()) if isinstance(size, np.ndarray) else tuple(size),
        method=str(method),
        num_iter=int(num_iter),
        max_dist=float(max_dist),
        resample=bool(resample),
        lambert_proj=lambert_proj,
        lambert_grid=lambert_grid,
        auto_proj=bool(auto_proj),
    )

    da = xr.DataArray(
        field,
        dims=(lat_dim, lon_dim),
        coords={lon_dim: gridX, lat_dim: gridY},
        name=var_name,
        attrs={
            "method": str(method),
            "grid_res_deg": float(step),
            "sigma_deg": float(sigma_deg),
            "sigma_grid": float(sigma_grid),
            "influence_radius_deg": float(influence_radius_deg),
            "max_dist_sigma": float(max_dist),
            "num_iter": int(num_iter),
            "resample": bool(resample),
            "auto_domain": bool(auto_domain),
            "buffer_deg": float(buffer_deg),
            "domain_lon0": float(lon0),
            "domain_lat0": float(lat0),
            "domain_grid_x_deg": float(gx),
            "domain_grid_y_deg": float(gy),
            "auto_proj": bool(auto_proj),
        },
    )

    # -------------------------
    # Optional mask (same style as interp_spatial_barnes)
    # Note: mask uses lon/lat degree distance (planar) as a pragmatic cutoff.
    # -------------------------
    if mask_radius_deg is not None:
        tree, pts_ref = kdtree.create_kdtree(pts)
        kd_search = kdtree.prepare_search(float(mask_radius_deg), tree, pts_ref)
        res_index, _, _, _, _ = kd_search

        mask = np.zeros((ny, nx), dtype=bool)
        c = np.empty(2, dtype=np.float64)
        for j in range(ny):
            c[1] = lat0 + j * step
            for i in range(nx):
                c[0] = lon0 + i * step
                kdtree.radius_search(c, *kd_search)
                mask[j, i] = res_index[-1] > 0

        da = da.where(mask)
        da.attrs["mask_radius_deg"] = float(mask_radius_deg)

    return da


def interp_spatial_barnesS2_rs(
    data: pd.DataFrame,
    var_name: str,
    *,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    # grid (degrees)
    grid_res_deg: float = 0.25,
    auto_domain: bool = True,
    buffer_deg: float = 1.0,
    point: tuple[float, float] | None = None,
    grid_x: float | None = None,
    grid_y: float | None = None,
    # Barnes (degrees)
    sigma_deg: float | None = None,  # None -> auto
    influence_radius_deg: float | None = None,  # None -> auto (=k*sigma)
    influence_k_sigma: float = 4.0,
    num_iter: int = 4,
    method: Literal[
        "optimized_convolution_S2", "naive_S2"
    ] = "optimized_convolution_S2",
    max_dist_sigma: float | None = None,  # optional override (in sigma units)
    resample: bool = True,
    # mask
    mask_radius_deg: float | None = None,  # None -> auto (=influence_radius_deg)
    # data cleaning
    missing_sentinels: tuple[float, ...] = (9999.9, 9999.0, -9999.0, -9999.9),
    # Lambert options (passed through to barnes_s2)
    lambert_proj=None,
    lambert_grid=None,
    auto_proj: bool = True,
) -> xr.DataArray:
    """
    Rust-accelerated Barnes interpolation on the unit sphere :math:`S^2` (spherical metric).

    This function mirrors :func:`interp_spatial_barnesS2` but uses the Rust backend
    (:func:`easyclimate_rust._easyclimate_rust.barnes_s2`) for the S2 interpolation.
    It keeps the same degree-based, user-friendly API and returns the same
    :py:class:`xarray.DataArray<xarray.DataArray>` layout.


    Parameters
    ----------
    Identical to :func:`interp_spatial_barnesS2`.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`
        Interpolated spherical Barnes field on a regular lat/lon grid, with metadata in ``attrs``.

    .. seealso::
        - :func:`interp_spatial_barnes` (planar metric)
        - :func:`interp_spatial_barnesS2_rs` (Rust backend on :math:`S^2`)
        - https://github.com/MeteoSwiss/fast-barnes-py
        - Zürcher, B. K.: Fast approximate Barnes interpolation: illustrated by Python-Numba implementation fast-barnes-py v1.0, Geosci. Model Dev., 16, 1697–1711, https://doi.org/10.5194/gmd-16-1697-2023, 2023.
    """

    if grid_res_deg <= 0:
        raise ValueError("grid_res_deg must be > 0")
    if num_iter < 1:
        raise ValueError("num_iter must be >= 1")

    # -------------------------
    # Clean data
    # -------------------------
    df = data[[lon_dim, lat_dim, var_name]].copy()
    for s in missing_sentinels:
        df.loc[np.isclose(df[var_name].astype(float), s), var_name] = np.nan
    df = df.replace([np.inf, -np.inf], np.nan).dropna(
        subset=[lon_dim, lat_dim, var_name]
    )
    if len(df) < 3:
        raise ValueError(f"Not enough valid points after cleaning: {len(df)}")

    stn_lon = df[lon_dim].astype(float).to_numpy()
    stn_lat = df[lat_dim].astype(float).to_numpy()
    values = df[var_name].astype(float).to_numpy()

    # -------------------------
    # Auto domain from bbox
    # -------------------------
    if auto_domain:
        lon_min, lon_max = float(np.min(stn_lon)), float(np.max(stn_lon))
        lat_min, lat_max = float(np.min(stn_lat)), float(np.max(stn_lat))
        lon0 = lon_min - buffer_deg
        lat0 = lat_min - buffer_deg
        gx = (lon_max - lon_min) + 2 * buffer_deg
        gy = (lat_max - lat_min) + 2 * buffer_deg
        point = (lon0, lat0)
        grid_x = gx
        grid_y = gy
    else:
        if point is None or grid_x is None or grid_y is None:
            raise ValueError(
                "When auto_domain=False, you must provide point, grid_x, grid_y"
            )
        lon0, lat0 = float(point[0]), float(point[1])
        gx, gy = float(grid_x), float(grid_y)

    step = float(grid_res_deg)
    x0, size, gridX, gridY = _build_lonlat_grid(
        point, float(grid_x), float(grid_y), step
    )
    nx = int(size[0])
    ny = int(size[1])

    # -------------------------
    # Auto sigma (deg) WITHOUT kNN
    # -------------------------
    if sigma_deg is None:
        area = max(float(grid_x) * float(grid_y), 1e-6)  # deg^2
        spacing = sqrt(area / len(df))
        sigma_deg = float(np.clip(0.6 * spacing, 0.2, 3.0))
    else:
        sigma_deg = float(sigma_deg)
        if sigma_deg <= 0:
            raise ValueError("sigma_deg must be > 0")

    # -------------------------
    # Influence radius + max_dist in sigma units
    # -------------------------
    if influence_radius_deg is None:
        influence_radius_deg = influence_k_sigma * sigma_deg
    influence_radius_deg = float(influence_radius_deg)

    if max_dist_sigma is None:
        max_dist = influence_radius_deg / sigma_deg
    else:
        max_dist = float(max_dist_sigma)

    if mask_radius_deg is None:
        mask_radius_deg = influence_radius_deg
    mask_radius_deg = float(mask_radius_deg) if mask_radius_deg is not None else None

    # Convert to grid units for rust barnes_s2
    sigma_grid = float(sigma_deg) / step

    pts = np.c_[stn_lon, stn_lat].astype(np.float64)

    field = barnesS2_rs(
        pts,
        values.astype(np.float64),
        float(sigma_grid),
        x0.astype(np.float64),
        float(step),
        tuple(size.tolist()) if isinstance(size, np.ndarray) else tuple(size),
        method=str(method),
        num_iter=int(num_iter),
        max_dist=float(max_dist),
        resample=bool(resample),
        lambert_proj=lambert_proj,
        lambert_grid=lambert_grid,
        auto_proj=bool(auto_proj),
    )

    da = xr.DataArray(
        field,
        dims=(lat_dim, lon_dim),
        coords={lon_dim: gridX, lat_dim: gridY},
        name=var_name,
        attrs={
            "method": str(method),
            "grid_res_deg": float(step),
            "sigma_deg": float(sigma_deg),
            "sigma_grid": float(sigma_grid),
            "influence_radius_deg": float(influence_radius_deg),
            "max_dist_sigma": float(max_dist),
            "num_iter": int(num_iter),
            "resample": bool(resample),
            "auto_domain": bool(auto_domain),
            "buffer_deg": float(buffer_deg),
            "domain_lon0": float(lon0),
            "domain_lat0": float(lat0),
            "domain_grid_x_deg": float(gx),
            "domain_grid_y_deg": float(gy),
            "auto_proj": bool(auto_proj),
        },
    )

    if mask_radius_deg is not None:
        mask = radius_mask_2d_rs(
            pts.astype(np.float64),
            x0.astype(np.float64),
            np.asarray([step, step], dtype=np.float64),
            np.asarray([nx, ny], dtype=np.int64),
            float(mask_radius_deg),
        )
        da = da.where(mask)
        da.attrs["mask_radius_deg"] = float(mask_radius_deg)

    return da
